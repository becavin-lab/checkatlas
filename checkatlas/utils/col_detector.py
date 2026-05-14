"""
Single-Cell Column Detector for CheckAtlas
===========================================

Automated detection of relevant columns and embedding keys in single-cell
AnnData objects for quality assessment across annotation, clustering, and
dimensionality reduction tasks.

Architecture
------------
The detector employs a multi-strategy scoring framework:

1. **Semantic Analysis**: Regex pattern matching against curated dictionaries
   of known single-cell tool outputs, biological conventions, and naming
   patterns with weighted confidence scores.

2. **Statistical Profiling**: Analysis of data distribution properties
   including cardinality, completeness, categorical nature, and label
   characteristics (average string length, numeric ratio).

3. **Cross-contamination Penalties**: Penalises columns whose name
   simultaneously matches multiple competing categories (e.g. a column
   named ``celltype_predicted`` cannot be a ground-truth reference).

Task Categories Detected
------------------------
* **Reference annotations** — ground-truth / author-curated cell-type labels
* **Predicted annotations** — automated / inferred cell-type predictions
  from tools such as CellTypist, SingleR, scVI, Azimuth, Garnett, etc.
* **Cluster labels** — algorithmic clustering results: Leiden, Louvain,
  k‑means, Seurat clusters, PhenoGraph, etc.
* **Embeddings** — dimensionality-reduction representations stored in
  ``.obsm`` (e.g. ``X_umap``, ``X_pca``, ``X_tsne``, ``X_scvi``)
* **Metadata** — columns explicitly excluded from analysis (batch,
  barcode, QC metrics, cell-cycle phase, etc.)
"""

import re
import logging
from typing import Dict, List, Tuple, Optional, Any

import numpy as np
import pandas as pd
from anndata import AnnData

logger = logging.getLogger("checkatlas")

# ═══════════════════════════════════════════════════════════════════════
# Pattern Registry — curated regex patterns with confidence weights
# ═══════════════════════════════════════════════════════════════════════

_REFERENCE_ANNOTATION_PATTERNS: Dict[str, float] = {
    # ── Ground truth / author-curated ──────────────────────────────
    r"\bcell[_\s]?type\b":                              1.00,
    r"\bauthor[_\s]?cell[_\s]?type\b":                  0.95,
    r"\bground[_\s]?truth\b":                           0.95,
    r"\breference\b":                                   0.85,
    r"\bcurated\b":                                     0.85,
    r"\bontology\b":                                    0.75,

    # ── Verified / known ──────────────────────────────────────────
    r"\btrue[_\s]?(label|annotation|celltype)\b":       0.90,
    r"\boriginal[_\s]?(label|annotation|celltype)\b":   0.90,
    r"\bexpert[_\s]?(label|annotation)\b":              0.90,
    r"\bmanual[_\s]?(label|annotation)\b":              0.90,
    r"\bknown[_\s]?(label|annotation)\b":               0.85,
    r"\bverified[_\s]?(label|annotation)\b":            0.85,

    # ── Annotation-specific naming ────────────────────────────────
    r"\bann(otation)?[_\s]?(finest|level|label)\b":     0.85,
    r"\bcell[_\s]?type[_\s]?annotation\b":             0.95,
    r"\bannot[_\s]?cell\b":                             0.80,

    # ── Biological groupings ──────────────────────────────────────
    r"\blineage\b":                                     0.60,
    r"\b(?<!sub)type\b":                                0.40,

    # ── Generic label terms (low weight — easy to false-positive) ─
    r"(^|[_\s])annotation($|[_\s])":                    0.70,
    r"(^|[_\s])annotated($|[_\s])":                     0.65,
    r"\blabel\b":                                       0.35,
}

_PREDICTED_ANNOTATION_PATTERNS: Dict[str, float] = {
    # ── Automated annotation tools ────────────────────────────────
    r"(?<![a-zA-Z])celltypist(?![a-zA-Z])":             1.00,
    r"(?<![a-zA-Z])single[_\s]?r(?![a-zA-Z])":          0.95,
    r"(?<![a-zA-Z])sctype(?![a-zA-Z])":                 0.95,
    r"(?<![a-zA-Z])azimuth(?![a-zA-Z])":                0.95,
    r"(?<![a-zA-Z])cellassign(?![a-zA-Z])":             0.95,
    r"(?<![a-zA-Z])garnett(?![a-zA-Z])":                0.90,
    r"(?<![a-zA-Z])scmap(?![a-zA-Z])":                  0.90,
    r"(?<![a-zA-Z])chetah(?![a-zA-Z])":                 0.90,
    r"(?<![a-zA-Z])scpred(?![a-zA-Z])":                 0.85,
    r"(?<![a-zA-Z])cello(?![a-zA-Z])":                  0.85,
    r"(?<![a-zA-Z])onclass(?![a-zA-Z])":                0.85,
    r"(?<![a-zA-Z])scarches(?![a-zA-Z])":               0.85,
    r"(?<![a-zA-Z])scclassif[yi](?![a-zA-Z])":          0.85,
    r"(?<![a-zA-Z])cellid(?![a-zA-Z])":                 0.80,
    r"(?<![a-zA-Z])scina(?![a-zA-Z])":                  0.80,
    r"(?<![a-zA-Z])itclust(?![a-zA-Z])":                0.75,
    r"(?<![a-zA-Z])cellhint(?![a-zA-Z])":               0.85,

    # ── Deep generative / probabilistic ───────────────────────────
    r"(?<![a-zA-Z])sc[_\s]?vi(?![a-zA-Z])":             0.90,
    r"(?<![a-zA-Z])scanvi(?![a-zA-Z])":                 0.90,
    r"(?<![a-zA-Z])scgen(?![a-zA-Z])":                  0.80,
    r"(?<![a-zA-Z])scvi[_\s]?tools?(?![a-zA-Z])":       0.90,

    # ── Generic predicted / inferred / automated ─────────────────
    r"(?<![a-zA-Z])predict(ed|ion)?(?![a-zA-Z])":       0.85,
    r"(?<![a-zA-Z])infer(red)?(?![a-zA-Z])":            0.85,
    r"(?<![a-zA-Z])auto(matic|mated)?(?![a-zA-Z])":     0.80,
    r"(?<![a-zA-Z])transfer(red)?[_\s]?label(?![a-zA-Z])": 0.85,
    r"(?<![a-zA-Z])mapped?[_\s]?(label|annotation)(?![a-zA-Z])": 0.80,

    # ── Ensemble / consensus ─────────────────────────────────────
    r"(?<![a-zA-Z])consensus(?![a-zA-Z])":              0.80,
    r"(?<![a-zA-Z])majority[_\s]?vote(?![a-zA-Z])":     0.80,
    r"(?<![a-zA-Z])ensemble(?![a-zA-Z])":               0.75,

    # ── Harmonised / integrated annotations ──────────────────────
    r"(?<![a-zA-Z])harmon(ized|ised|y)[_\s]?(label|annotation)?(?![a-zA-Z])": 0.60,
    r"(?<![a-zA-Z])integrat(ed|ion)[_\s]?(label|annotation)(?![a-zA-Z])": 0.60,

    # ── Cross-species / reference-based ─────────────────────────
    r"(?<![a-zA-Z])scsimilarity(?![a-zA-Z])":           0.80,
    r"(?<![a-zA-Z])scjoint(?![a-zA-Z])":                0.75,
}

_CLUSTER_LABEL_PATTERNS: Dict[str, float] = {
    # ── Graph-based community detection ──────────────────────────
    r"(?<![a-zA-Z])leiden(?![a-zA-Z])":                 1.00,
    r"(?<![a-zA-Z])louvain(?![a-zA-Z])":                1.00,
    r"(?<![a-zA-Z])walktrap(?![a-zA-Z])":               0.95,
    r"(?<![a-zA-Z])infomap(?![a-zA-Z])":                0.90,

    # ── Seurat-specific ──────────────────────────────────────────
    r"(?<![a-zA-Z])seurat[_\s]?clusters?(?![a-zA-Z])":  1.00,
    r"(?<![a-zA-Z])rna[_\s]?snn[_\s]?res":              0.95,
    r"(?<![a-zA-Z])sct[_\s]?snn[_\s]?res":              0.95,
    r"(?<![a-zA-Z])integrated[_\s]?snn[_\s]?res":      0.90,
    r"(?<![a-zA-Z])wsnn[_\s]?res":                     0.90,

    # ── Graph-based / other ─────────────────────────────────────
    r"(?<![a-zA-Z])graph[_\s]?clust":                  0.90,
    r"(?<![a-zA-Z])phenograph(?![a-zA-Z])":             0.85,

    # ── Partitional ─────────────────────────────────────────────
    r"(?<![a-zA-Z])k[_\s]?means?(?![a-zA-Z])":          0.85,
    r"(?<![a-zA-Z])kmeans(?![a-zA-Z])":                 0.85,
    r"(?<![a-zA-Z])hdbscan(?![a-zA-Z])":                0.85,
    r"(?<![a-zA-Z])dbscan(?![a-zA-Z])":                 0.80,

    # ── scVI-then-cluster ──────────────────────────────────────
    r"(?<![a-zA-Z])scvi[_\s]?clust(er)?(?![a-zA-Z])":   0.85,
    r"(?<![a-zA-Z])scanvi[_\s]?clust(er)?(?![a-zA-Z])": 0.80,

    # ── Generic clustering (lower weight, high false-positive risk)
    r"(?<![a-zA-Z])cluster(s|ing)?(?![a-zA-Z])":        0.65,
    r"(?<![a-zA-Z])partition(?![a-zA-Z])":              0.55,
}

# Embedding patterns are intentionally *ordered* — the first match wins
# (specificity-gated).  The catch‑all ``^X_`` sits at the end so that a
# key like ``X_umap`` is first matched by the exact ``^X_umap$`` pattern
# and receives full confidence instead of being slurped by the fallback.
_EMBEDDING_PATTERNS: List[Tuple[str, float]] = [
    (r"^X_pca$",                                       1.00),
    (r"^X_umap$",                                      1.00),
    (r"^X_tsne$",                                      1.00),
    (r"^X_diffmap$",                                   0.95),
    (r"^X_draw_graph",                                 0.90),
    (r"^X_scvi$",                                      0.90),
    (r"^X_scanvi$",                                    0.90),
    (r"^X_mde$",                                       0.85),
    (r"^X_phate$",                                     0.85),
    (r"^X_pacmap$",                                    0.85),
    (r"^X_trimap$",                                    0.85),
    (r"^X_.*pca",                                      0.80),
    (r"^X_.*umap",                                     0.75),
    (r"^X_.*tsne",                                     0.75),
    (r"^X_",                                           0.50),
]

_METADATA_PATTERNS: Dict[str, float] = {
    # ── Experimental metadata ────────────────────────────────────
    r"(?<![a-zA-Z])(batch|sample|donor|patient|replicate)(?![a-zA-Z])": 1.00,
    r"(?<![a-zA-Z])(library|sequencing|chemistry|plate|well)(?![a-zA-Z])": 0.90,
    r"(?<![a-zA-Z])(tissue|organ|disease|condition|treatment)(?![a-zA-Z])": 0.70,
    r"(?<![a-zA-Z])(stim|stimulated|control|case)(?![a-zA-Z])": 0.65,
    r"(?<![a-zA-Z])dataset(?![a-zA-Z])":                0.90,
    r"(?<![a-zA-Z])(assay|protocol|platform)(?![a-zA-Z])": 0.70,
    r"(?<![a-zA-Z])(sex|gender)(?![a-zA-Z])":           0.65,
    r"(?<![a-zA-Z])(suspension|capture|preparation)(?![a-zA-Z])": 0.65,
    r"(?<![a-zA-Z])(organism|species)(?![a-zA-Z])":     0.65,
    r"(?<![a-zA-Z])development[_\s]?stage(?![a-zA-Z])":  0.70,

    # ── Cell identifiers ────────────────────────────────────────
    r"(?<![a-zA-Z])(barcode|cell[_\s]?id|obs[_\s]?names)(?![a-zA-Z])": 1.00,
    r"(?<![a-zA-Z])index(?![a-zA-Z])":                  0.95,

    # ── QC metrics ──────────────────────────────────────────────
    r"(?<![a-zA-Z])n[_\s]?(genes|counts|umis?)(?![a-zA-Z])": 0.90,
    r"(?<![a-zA-Z])total[_\s]?(counts|umi)(?![a-zA-Z])": 0.90,
    r"(?<![a-zA-Z])(percent|pct)[_\s]":                0.85,
    r"(?<![a-zA-Z])(doublet|scrublet|simian)(?![a-zA-Z])": 0.85,
    r"(?<![a-zA-Z])(mito|ribo|mitochondrial|ribosomal)(?![a-zA-Z])": 0.75,
    r"(?<![a-zA-Z])n_?feature(?![a-zA-Z])":             0.85,
    r"(?<![a-zA-Z])(log1p|log|normalized)(?![a-zA-Z])": 0.60,

    # ── Cell cycle / biology ────────────────────────────────────
    r"(?<![a-zA-Z])(phase|s[_\s]?phase|g2m[_\s]?score)(?![a-zA-Z])": 0.70,

    # ── Dimensionality artifacts ────────────────────────────────
    r"(?<![a-zA-Z])n[_\s]?(pc|count)(?![a-zA-Z])":      0.70,
}

# ═══════════════════════════════════════════════════════════════════════
# Statistical thresholds — domain-guided bounds for single-cell data
# ═══════════════════════════════════════════════════════════════════════

# Reference cell-type annotations typically contain 3 – 200 unique
# categories, are categorical, and carry biologically meaningful names.
_REF_CARDINALITY_RANGE = (3, 200)
_REF_CARDINALITY_PENALTY_UPPER = 500

# Predicted annotations / automated labels usually 2 – 100 categories.
_PRED_CARDINALITY_RANGE = (2, 100)
_PRED_CARDINALITY_PENALTY_UPPER = 200

# Cluster labels (Leiden, Louvain, k‑means) are compact: 2 – 50 groups.
_CLUST_CARDINALITY_RANGE = (2, 50)
_CLUST_CARDINALITY_PENALTY_UPPER = 150


# ═══════════════════════════════════════════════════════════════════════
# CheckAtlasColumnDetector
# ═══════════════════════════════════════════════════════════════════════

class CheckAtlasColumnDetector:
    """Multi-strategy column detector for single-cell AnnData objects.

    Detects parameters for three quality-assessment tasks:

    1. **Annotation** — reference (ground-truth) *vs.* predicted
       (computational) cell-type labelling.
    2. **Clustering** — embedding representations and cluster-label columns
       for cluster-quality metrics (silhouette, Davies–Bouldin, etc.).
    3. **Dimensionality reduction** — low-dimensional embedding keys stored
       in ``.obsm``.

    The detector combines three independent analysis strands:

    * **Semantic** — regex pattern matching against curated dictionaries of
      single-cell tool outputs and naming conventions.
    * **Statistical** — distribution properties of the underlying data
      (cardinality, completeness, categorical *vs.* numeric, average label
      length).
    * **Cross-contamination** — penalises columns whose name simultaneously
      matches conflicting categories (e.g. a column containing ``predicted``
      cannot be a trustworthy reference annotation).

    Parameters
    ----------
    adata : AnnData
        Scanpy AnnData object to analyse.
    """

    def __init__(self, adata: AnnData) -> None:
        self.adata = adata
        self.obs_columns: List[str] = list(adata.obs.columns)
        self.obsm_keys: List[str] = (
            list(adata.obsm.keys()) if hasattr(adata, "obsm") else []
        )
        self._semantic_cache: Dict[str, Dict[str, float]] = {}
        self._stats_cache: Dict[str, Dict[str, float]] = {}

    # ── Semantic analysis ───────────────────────────────────────────

    def analyze_column_semantics(self, col_name: str) -> Dict[str, float]:
        """Score a column name against all semantic pattern dictionaries.

        Parameters
        ----------
        col_name : str
            Column name to analyse.

        Returns
        -------
        Dict[str, float]
            Scores for ``reference_annotation``, ``predicted_annotation``,
            ``cluster_label``, and ``metadata``.
        """
        if col_name in self._semantic_cache:
            return self._semantic_cache[col_name]

        col_lower = col_name.lower()
        scores: Dict[str, float] = {
            "reference_annotation": 0.0,
            "predicted_annotation": 0.0,
            "cluster_label": 0.0,
            "metadata": 0.0,
        }

        for pattern, weight in _REFERENCE_ANNOTATION_PATTERNS.items():
            if re.search(pattern, col_lower):
                scores["reference_annotation"] = max(
                    scores["reference_annotation"], weight
                )

        for pattern, weight in _PREDICTED_ANNOTATION_PATTERNS.items():
            if re.search(pattern, col_lower):
                scores["predicted_annotation"] = max(
                    scores["predicted_annotation"], weight
                )

        for pattern, weight in _CLUSTER_LABEL_PATTERNS.items():
            if re.search(pattern, col_lower):
                scores["cluster_label"] = max(
                    scores["cluster_label"], weight
                )

        for pattern, weight in _METADATA_PATTERNS.items():
            if re.search(pattern, col_lower):
                scores["metadata"] = max(scores["metadata"], weight)

        self._semantic_cache[col_name] = scores
        return scores

    def analyze_obsm_semantics(self) -> List[Tuple[str, float]]:
        """Score all ``.obsm`` keys against embedding patterns.

        Patterns are ordered by specificity; the first match for each key
        determines the confidence score.  A fallback ``^X_`` pattern at
        the end of the list captures any remaining ``X_``-prefixed keys
        with reduced confidence.

        Returns
        -------
        List[Tuple[str, float]]
            List of ``(key_name, confidence_score)`` sorted by score
            descending.
        """
        results: List[Tuple[str, float]] = []

        for key in self.obsm_keys:
            best_score = 0.0
            for pattern, weight in _EMBEDDING_PATTERNS:
                if re.search(pattern, key):
                    best_score = weight
                    break  # first (most-specific) match wins
            if best_score > 0:
                results.append((key, best_score))

        results.sort(key=lambda x: x[1], reverse=True)
        return results

    # ── Statistical profiling ───────────────────────────────────────

    def analyze_column_statistics(self, col_name: str) -> Dict[str, float]:
        """Compute the statistical profile of a column's data.

        Parameters
        ----------
        col_name : str
            Column name in ``adata.obs`` to profile.

        Returns
        -------
        Dict[str, float]
            Keys: ``cardinality_ratio``, ``completeness``,
            ``is_categorical``, ``label_length_avg``, ``numeric_ratio``,
            ``n_unique``.
        """
        if col_name in self._stats_cache:
            return self._stats_cache[col_name]

        if col_name not in self.adata.obs.columns:
            empty: Dict[str, float] = {
                "cardinality_ratio": 0.0,
                "completeness": 0.0,
                "is_categorical": 0.0,
                "label_length_avg": 0.0,
                "numeric_ratio": 0.0,
                "n_unique": 0,
            }
            self._stats_cache[col_name] = empty
            return empty

        data = self.adata.obs[col_name]
        total = len(data)

        stats: Dict[str, float] = {
            "cardinality_ratio": 0.0,
            "completeness": 0.0,
            "is_categorical": 0.0,
            "label_length_avg": 0.0,
            "numeric_ratio": 0.0,
            "n_unique": 0,
        }

        if total == 0:
            self._stats_cache[col_name] = stats
            return stats

        # Completeness — fraction of non-null entries
        non_null = data.notna().sum()
        stats["completeness"] = non_null / total

        # Cardinality
        n_unique = data.nunique()
        stats["n_unique"] = n_unique
        stats["cardinality_ratio"] = n_unique / total

        # Data-type heuristics
        dtype_str = str(data.dtype)
        if dtype_str == "category" or dtype_str == "object":
            stats["is_categorical"] = 1.0
            valid_data = data.dropna()
            if len(valid_data) > 0:
                try:
                    stats["label_length_avg"] = float(
                        valid_data.astype(str).str.len().mean()
                    )
                except Exception:
                    pass
        else:
            try:
                numeric_data = pd.to_numeric(data, errors="coerce")
                stats["numeric_ratio"] = numeric_data.notna().sum() / total
            except Exception:
                pass

        self._stats_cache[col_name] = stats
        return stats

    # ── Scoring functions ───────────────────────────────────────────

    def score_reference_annotation(self, col_name: str) -> float:
        """Confidence score for a ground-truth (reference) annotation column.

        Scoring weights: **semantic 40 % + statistical 60 %**.

        A strong penalty is applied when the column name also matches
        predicted-annotation patterns (cross-contamination), because a
        column named ``celltype_predicted`` or ``azimuth_celltype`` should
        never be classified as a reference (ground-truth) annotation.

        Parameters
        ----------
        col_name : str
            Column name to score.

        Returns
        -------
        float
            Confidence in [0, 1].
        """
        semantic = self.analyze_column_semantics(col_name)
        stats = self.analyze_column_statistics(col_name)

        if col_name not in self.adata.obs.columns:
            return 0.0

        score = 0.0

        # ── Semantic contribution (40 %) ──────────────────────────
        score += semantic["reference_annotation"] * 0.4
        score -= semantic["metadata"] * 0.6

        # Cross-contamination: if the column strongly matches predicted
        # annotation patterns it cannot be a trustworthy reference.
        pred_contamination = semantic["predicted_annotation"]
        if pred_contamination > 0.3:
            score -= pred_contamination * 0.4

        # ── Statistical indicators (60 %) ─────────────────────────
        # Reference annotations must be categorical in nature.
        if stats["is_categorical"] < 0.5:
            return 0.0

        n_unique = int(stats["n_unique"])

        if _REF_CARDINALITY_RANGE[0] <= n_unique <= _REF_CARDINALITY_RANGE[1]:
            score += 0.25
        elif n_unique < _REF_CARDINALITY_RANGE[0]:
            score += 0.05
        elif n_unique > _REF_CARDINALITY_PENALTY_UPPER:
            score -= 0.20

        score += stats["completeness"] * 0.20

        # Label-length heuristic: genuine cell-type names are usually
        # longer than algorithmic cluster IDs (e.g. "T cell" vs "3").
        avg_len = stats["label_length_avg"]
        if avg_len > 10:
            score += 0.25
        elif avg_len > 5:
            score += 0.15

        return float(np.clip(score, 0.0, 1.0))

    def score_predicted_annotation(self, col_name: str) -> float:
        """Confidence score for a predicted / automated annotation column.

        Scoring weights: **semantic 60 % + statistical 40 %**.

        Columns from tools like CellTypist, SingleR, scVI, Azimuth, and
        Garnett receive high semantic scores.  Generic patterns such as
        ``predicted``, ``inferred``, and ``automated`` also contribute.

        Parameters
        ----------
        col_name : str
            Column name to score.

        Returns
        -------
        float
            Confidence in [0, 1].
        """
        semantic = self.analyze_column_semantics(col_name)
        stats = self.analyze_column_statistics(col_name)

        if col_name not in self.adata.obs.columns:
            return 0.0

        score = 0.0

        # ── Semantic contribution (60 %) ──────────────────────────
        score += semantic["predicted_annotation"] * 0.6
        score -= semantic["metadata"] * 0.6

        # ── Statistical indicators (40 %) ─────────────────────────
        n_unique = int(stats["n_unique"])

        if _PRED_CARDINALITY_RANGE[0] <= n_unique <= _PRED_CARDINALITY_RANGE[1]:
            score += 0.20
        elif _PRED_CARDINALITY_RANGE[1] < n_unique <= _PRED_CARDINALITY_PENALTY_UPPER:
            score += 0.15
        elif n_unique < _PRED_CARDINALITY_RANGE[0]:
            return 0.0
        elif n_unique > _PRED_CARDINALITY_PENALTY_UPPER:
            score -= 0.15

        score += stats["completeness"] * 0.10

        if stats["numeric_ratio"] > 0.9:
            score += 0.10
        elif stats["is_categorical"] > 0.5:
            score += 0.10

        return float(np.clip(score, 0.0, 1.0))

    def score_cluster_label(self, col_name: str) -> float:
        """Confidence score for an algorithmic cluster-label column.

        Clustering outputs (Leiden, Louvain, k‑means, Seurat clusters)
        are typically categorical with 2–50 unique values and often contain
        numeric or short string labels (``"0"``, ``"1"``, …).

        Scoring weights: **semantic 70 % + statistical 30 %**.

        Parameters
        ----------
        col_name : str
            Column name to score.

        Returns
        -------
        float
            Confidence in [0, 1].
        """
        semantic = self.analyze_column_semantics(col_name)
        stats = self.analyze_column_statistics(col_name)

        if col_name not in self.adata.obs.columns:
            return 0.0

        score = 0.0

        # ── Semantic contribution (70 %) ──────────────────────────
        # Algorithm names (leiden, louvain, kmeans) are highly specific.
        score += semantic["cluster_label"] * 0.7
        score -= semantic["metadata"] * 0.6

        # ── Statistical indicators (30 %) ─────────────────────────
        n_unique = int(stats["n_unique"])

        if _CLUST_CARDINALITY_RANGE[0] <= n_unique <= _CLUST_CARDINALITY_RANGE[1]:
            score += 0.15
        elif _CLUST_CARDINALITY_RANGE[1] < n_unique <= _CLUST_CARDINALITY_PENALTY_UPPER:
            score += 0.10
        elif n_unique < _CLUST_CARDINALITY_RANGE[0]:
            return 0.0
        elif n_unique > _CLUST_CARDINALITY_PENALTY_UPPER:
            score -= 0.10

        score += stats["completeness"] * 0.05

        # Cluster labels are frequently numeric (Leiden/Louvain output
        # integer codes by default).
        if stats["numeric_ratio"] > 0.9:
            score += 0.05
        elif stats["is_categorical"] > 0.5:
            score += 0.05

        return float(np.clip(score, 0.0, 1.0))

    # ── Main detection entry-point ──────────────────────────────────

    def detect_all_parameters(
        self,
        min_reference_score: float = 0.95,
        min_predicted_score: float = 0.6,
        min_cluster_score: float = 0.5,
    ) -> Dict[str, Any]:
        """Run the full detection pipeline on all ``.obs`` columns and
        ``.obsm`` keys.

        Parameters
        ----------
        min_reference_score : float
            Minimum confidence for reference annotations (default 0.95).
        min_predicted_score : float
            Minimum confidence for predicted annotations (default 0.6).
        min_cluster_score : float
            Minimum confidence for cluster labels (default 0.5).

        Returns
        -------
        Dict
            Structured results dictionary::

                {
                    "annotation": {
                        "reference": [(col_name, score), ...],
                        "predicted": [(col_name, score), ...]
                    },
                    "clustering": {
                        "embeddings": [(key, metadata_dict), ...],
                        "cluster_labels": [(col_name, score), ...]
                    }
                }

            Each ``metadata_dict`` contains ``shape``, ``n_components``
            and ``score`` keys.
        """
        results: Dict[str, Any] = {
            "annotation": {
                "reference": [],
                "predicted": [],
            },
            "clustering": {
                "embeddings": [],
                "cluster_labels": [],
            },
        }

        # ── Scan .obs columns ────────────────────────────────────
        _skip_exact = {"index", "cell_index", "barcode", "obs_names"}

        for col in self.obs_columns:
            # Skip internal / system columns
            if col.startswith("_"):
                continue
            if col.lower() in _skip_exact:
                continue

            try:
                ref_score = self.score_reference_annotation(col)
                pred_score = self.score_predicted_annotation(col)
                clust_score = self.score_cluster_label(col)
            except Exception as exc:
                logger.debug(
                    "Skipping column %r — scoring failed: %s", col, exc
                )
                continue

            if ref_score >= min_reference_score:
                results["annotation"]["reference"].append((col, ref_score))

            if pred_score >= min_predicted_score:
                results["annotation"]["predicted"].append((col, pred_score))

            if clust_score >= min_cluster_score:
                results["clustering"]["cluster_labels"].append(
                    (col, clust_score)
                )

        # ── Scan .obsm keys ──────────────────────────────────────
        matched = self.analyze_obsm_semantics()

        for key, score in matched:
            try:
                embedding_data = self.adata.obsm[key]
                shape = embedding_data.shape
                n_components = (
                    shape[1] if len(shape) > 1 else 1
                )
            except Exception:
                shape = (-1,)
                n_components = 0

            results["clustering"]["embeddings"].append(
                (
                    key,
                    {
                        "shape": shape,
                        "n_components": n_components,
                        "score": score,
                    },
                )
            )

        # ── Sort by confidence (descending) ───────────────────────
        results["annotation"]["reference"].sort(
            key=lambda x: x[1], reverse=True
        )
        results["annotation"]["predicted"].sort(
            key=lambda x: x[1], reverse=True
        )
        results["clustering"]["cluster_labels"].sort(
            key=lambda x: x[1], reverse=True
        )

        return results
