"""Configuration for the scFM QC pipeline."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass(frozen=True)
class SCFMConfig:
    """Frozen configuration for one scfm run.

    Constructed in ``__main__.py`` from ``argparse.Namespace`` (or
    programmatically for tests). Passed to ``cal_scfm``,
    ``diagnose``, and ``composite`` in that order.
    """

    atlas_name: str = ""
    scfm_embedding: str = ""
    baseline_embeddings: tuple[str, ...] = field(default_factory=tuple)
    predicted_label: Optional[str] = None
    ref_label: Optional[str] = None
    batch_key: Optional[str] = None
    domain_key: Optional[str] = None
    patient_key: Optional[str] = None
    outcome_key: Optional[str] = None
    scaling_fractions: tuple[float, ...] = (0.01, 0.10, 0.50, 1.00)
    n_seeds: int = 5
    noise_sigma: float = 0.10
    min_domain_size: int = 50
    weights_path: Optional[str] = None
    thresholds_path: Optional[str] = None

    def all_embeddings(self) -> list[str]:
        out: list[str] = []
        if self.scfm_embedding:
            out.append(self.scfm_embedding)
        for e in self.baseline_embeddings:
            if e and e not in out:
                out.append(e)
        return out


def from_args(args) -> "SCFMConfig":
    """Build a config from an argparse Namespace (defensive on missing attrs)."""

    def _get(name: str, default=None):
        return getattr(args, name, default)

    def _tup(name: str, default=()):
        v = _get(name, default)
        if v is None or v == []:
            return tuple()
        if isinstance(v, str):
            return tuple(v.split())
        return tuple(v)

    def _ftup(name: str, default=(0.01, 0.10, 0.50, 1.00)):
        v = _get(name, default)
        if v is None or v == []:
            return tuple(default)
        if isinstance(v, str):
            return tuple(float(x) for x in v.split())
        return tuple(float(x) for x in v)

    return SCFMConfig(
        atlas_name=str(_get("atlas_name", "")),
        scfm_embedding=str(_get("scfm_embedding", "")),
        baseline_embeddings=_tup("baseline_embeddings", ()),
        predicted_label=_get("scfm_predicted_label", None),
        ref_label=_get("scfm_ref_label", None),
        batch_key=_get("scfm_batch_key", None),
        domain_key=_get("scfm_domain_key", None),
        patient_key=_get("scfm_patient_key", None),
        outcome_key=_get("scfm_outcome_key", None),
        scaling_fractions=_ftup("scfm_scaling_fractions"),
        n_seeds=int(_get("scfm_n_seeds", 5)),
        noise_sigma=float(_get("scfm_noise_sigma", 0.10)),
        min_domain_size=int(_get("scfm_min_domain_size", 50)),
        weights_path=_get("scfm_weights", None),
        thresholds_path=_get("scfm_thresholds", None),
    )
