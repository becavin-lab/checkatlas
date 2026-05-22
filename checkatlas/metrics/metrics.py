import logging
import time
import pandas as pd
import numpy as np
from tqdm import tqdm
import os

try:
    import rpy2.robjects as ro
    import rpy2.robjects as robjects
except (ImportError, OSError):
    ro = None
    robjects = None
from sklearn.preprocessing import LabelEncoder

from . import annot, cluster, dimred
# Import CheckAtlasColumnDetector inside function or at top if no circular dependency
# from ..atlas import CheckAtlasColumnDetector # This might cause circular import if atlas imports metrics

METRICS_CLUST = cluster.__all__
METRICS_ANNOT = annot.__all__
METRICS_DIMRED = dimred.__all__


def _pivot_annot_to_wide(df, atlas_name):
    """Pivot long-format annotation results to MultiQC-compatible wide format.

    Converts a DataFrame with columns [Atlas Name, Metric Name,
    Reference/Input 1, Prediction/Input 2, Value, Time (s)] into a wide
    table headed by [Annot_Sample, Reference, obs, <metric>_value,
    <metric>_running_time, ...] where each row is one unique
    (Reference/Input 1 × Prediction/Input 2) combination.
    """
    if df.empty:
        return df
    metric_names = sorted(df['Metric Name'].unique())
    wide_rows = []
    for (ref_val, pred_val), group in df.groupby(['Reference/Input 1', 'Prediction/Input 2']):
        row = {
            'Annot_Sample': f"{atlas_name}_{ref_val}_{pred_val}",
            'Reference': ref_val,
            'obs': pred_val,
        }
        for _, r in group.iterrows():
            metric = r['Metric Name']
            row[metric] = r['Value']
            row[f"{metric}_running_time"] = r['Time (s)']
        wide_rows.append(row)
    result_df = pd.DataFrame(wide_rows)
    value_cols = []
    for m in metric_names:
        value_cols.append(m)
        value_cols.append(f"{m}_running_time")
    col_order = ['Annot_Sample', 'Reference', 'obs'] + [c for c in value_cols if c in result_df.columns]
    return result_df[col_order]


def _pivot_cluster_to_wide(df, atlas_name):
    """Pivot long-format cluster results to MultiQC-compatible wide format.

    Converts a DataFrame with columns [Atlas Name, Metric Name,
    Embedding, Label Key, Value, Time (s)] into a wide table headed by
    [Clust_Sample, obs, <metric>, <metric>_running_time, ...] where each
    row is one unique (Embedding × Label Key) combination.
    """
    if df.empty:
        return df
    metric_names = sorted(df['Metric Name'].unique())
    wide_rows = []
    for (emb_val, label_val), group in df.groupby(['Embedding', 'Label Key']):
        row = {
            'Clust_Sample': f"{atlas_name}_{emb_val}_{label_val}",
            'obs': label_val,
        }
        for _, r in group.iterrows():
            metric = r['Metric Name']
            row[metric] = r['Value']
            row[f"{metric}_running_time"] = r['Time (s)']
        wide_rows.append(row)
    result_df = pd.DataFrame(wide_rows)
    value_cols = []
    for m in metric_names:
        value_cols.append(m)
        value_cols.append(f"{m}_running_time")
    col_order = ['Clust_Sample', 'obs'] + [c for c in value_cols if c in result_df.columns]
    return result_df[col_order]


def _pivot_dimred_to_wide(df, atlas_name):
    """Pivot long-format dimred results to MultiQC-compatible wide format.

    Converts a DataFrame with columns [Metric Name, Low Dim Key,
    High Dim Key, Value, Time (s)] into a wide table headed by
    [Dimred_Sample, obsm, <metric>, <metric>_running_time, …] where
    each row is one unique Low Dim Key (embedding).
    """
    if df.empty:
        return df
    metric_names = sorted(df['Metric Name'].unique())
    wide_rows = []
    for low_dim_key, group in df.groupby('Low Dim Key'):
        row = {
            'Dimred_Sample': f"{atlas_name}_{low_dim_key}",
            'obsm': low_dim_key,
        }
        for _, r in group.iterrows():
            metric = r['Metric Name']
            row[metric] = r['Value']
            row[f"{metric}_running_time"] = r['Time (s)']
        wide_rows.append(row)
    result_df = pd.DataFrame(wide_rows)
    value_cols = []
    for m in metric_names:
        value_cols.append(m)
        value_cols.append(f"{m}_running_time")
    col_order = ['Dimred_Sample', 'obsm'] + [c for c in value_cols if c in result_df.columns]
    return result_df[col_order]


logger = logging.getLogger("checkatlas")

if robjects is not None:
    R_ANNOT = robjects.r(
        "type <- function(seurat, obs_key){ "
        "return(seurat[[obs_key]][[obs_key]])}"
    )
    R_REDUCTION = robjects.r(
        "reduc <- function(seurat, obsm_key){"
        " return(Embeddings(object = seurat, "
        "reduction = obsm_key))}"
    )
else:
    R_ANNOT = None
    R_REDUCTION = None

def cal_annot(adata, atlas_name=None, metric_list=None, all=False,
              file_dir=None, n_jobs=-1, verbose=True):
    """
    Comprehensive annotation pipeline for all annotation metrics.
    
    Args:
        adata (AnnData): Annotated data matrix.
        metric_list (list): List of metric names to calculate.
                           If provided, overrides `all` parameter.
        all (bool): If True, calculate all available annotation metrics. 
                    If False, calculate a default subset.
                    Ignored if metric_list is provided.
        file_dir (str): Directory path where the results CSV will be saved.
                       If None, saves to current working directory.
        n_jobs (int): Number of parallel jobs (-1 = all cores).
        verbose (bool): Whether to print progress information.
                    
    Returns:
        pd.DataFrame: Results dataframe with columns:
                       [Atlas Name, Metric Name, Reference/Input 1, Prediction/Input 2, Value, Time (s)]
    """
    import inspect
    from ..utils.col_detector import CheckAtlasColumnDetector
    
    # Set file directory
    if file_dir is None:
        file_dir = os.getcwd()
    else:
        # Create directory if it doesn't exist
        os.makedirs(file_dir, exist_ok=True)
    
    # Detect columns
    detector = CheckAtlasColumnDetector(adata)
    params = detector.detect_all_parameters()
    
    ref_keys = [x[0] for x in params['annotation']['reference']]
    pred_keys = [x[0] for x in params['annotation']['predicted']]
    embedding_keys = [x[0] for x in params['clustering']['embeddings']]
    
    # Detect batch keys using the column detector (semantic + statistical)
    # Falls back to the simple 'batch' substring heuristic if detector returns nothing
    batch_keys = [x[0] for x in params.get("batch", [])]
    if not batch_keys:
        batch_keys = [col for col in adata.obs.columns if 'batch' in col.lower()]
    
    # Define metrics to run
    if metric_list is not None:
        metrics_list = [m for m in metric_list if m in METRICS_ANNOT]
    elif all:
        metrics_list = METRICS_ANNOT
    else:
        metrics_list = [
            'adj_rand_index', 
            'normalized_mutual_info', 
            'adj_mutual_info',
        ]
        metrics_list = [m for m in metrics_list if m in METRICS_ANNOT]

    results = []
    atlas_name = atlas_name

    # Categorize metrics based on their input requirements
    # Ref vs Pred
    ref_pred_metrics = [
        'adj_mutual_info', 'adj_rand_index', 'fowlkes_mallows', 
        'isolated_f1_score', 'mutual_info', 'normalized_mutual_info', 
        'rand_index', 'vmeasure'
    ]
    
    # Embedding + Labels
    emb_label_metrics = ['average_silhouette_width', 'dunn_index']
    
    # Batch / Integration (adata + batch/label)
    batch_metrics = ['kbet', 'pcr'] # lisi is special (iLISI vs cLISI)
    
    # Graph Connectivity (adata + neighbors)
    graph_metrics = ['graph_connectivity']
    
    # Bio Conservation (adata_before, adata_after) - Skipping for single adata pipeline
    # unless we define strategy. 
    bio_metrics = ['cell_cycle_conservation', 'highly_variable_genes']


    # Create progress bar with custom format
    pbar = tqdm(metrics_list, desc="Calculating Annotation Metrics", 
                bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]')
    
    for metric in pbar:
        # Start timing for this metric
        metric_start_time = time.time()
        
        # Update progress bar with current metric name
        pbar.set_description(f"Processing: {metric}")
        
        metric_module = getattr(annot, metric)
        
        try:
            # 1. Ref vs Pred Metrics
            if metric in ref_pred_metrics:
                if not ref_keys or not pred_keys:
                    continue
                for ref in ref_keys:
                    for pred in pred_keys:
                        # Skip if ref == pred
                        if ref == pred:
                            continue
                        try:
                            # Preprocess labels
                            labels_true = adata.obs[ref]
                            labels_pred = adata.obs[pred]
                            # Convert to numeric if needed (some metrics handle it, some don't)
                            # Most sklearn metrics handle strings, but let's be safe or rely on metric impl
                            # checkatlas metrics usually take raw inputs or handle conversion
                            # But calc_metric_annot_scanpy uses annotation_to_num.
                            # We should probably use that helper or do it here.
                            l_pred, l_true = annotation_to_num(labels_pred, labels_true)
                            
                            pair_start = time.time()
                            # Build kwargs dynamically
                            sig = inspect.signature(metric_module.run)
                            kw = {}
                            if 'n_jobs' in sig.parameters: kw['n_jobs'] = n_jobs
                            if 'verbose' in sig.parameters: kw['verbose'] = False
                            val = metric_module.run(l_pred, l_true, **kw)
                            pair_elapsed = time.time() - pair_start
                            results.append({
                                'Atlas Name': atlas_name,
                                'Metric Name': metric,
                                'Reference/Input 1': ref,
                                'Prediction/Input 2': pred,
                                'Value': val,
                                'Time (s)': round(pair_elapsed, 3)
                            })
                        except Exception as e:
                            logger.warning(f"Failed to calculate {metric} for {ref} vs {pred}: {e}")

            # 2. Embedding + Labels (ASW, Dunn)
            elif metric in emb_label_metrics:
                if not embedding_keys:
                    continue
                # Run for both ref and pred labels? Usually for predicted clusters.
                # But ASW can be run on ground truth too.
                targets = list(set(ref_keys + pred_keys))
                for emb in embedding_keys:
                    if emb not in adata.obsm:
                        continue
                    X_emb = adata.obsm[emb]
                    for label in targets:
                        try:
                            labels = adata.obs[label]
                            # Convert to numeric for ASW/Dunn?
                            # ASW handles labels.
                            pair_start = time.time()
                            sig = inspect.signature(metric_module.run)
                            kw = {}
                            if 'n_jobs' in sig.parameters: kw['n_jobs'] = n_jobs
                            if 'verbose' in sig.parameters: kw['verbose'] = False
                            val = metric_module.run(X_emb, labels, **kw)
                            pair_elapsed = time.time() - pair_start
                            results.append({
                                'Atlas Name': atlas_name,
                                'Metric Name': metric,
                                'Reference/Input 1': emb,
                                'Prediction/Input 2': label,
                                'Value': val,
                                'Time (s)': round(pair_elapsed, 3)
                            })
                        except Exception as e:
                            logger.warning(f"Failed to calculate {metric} for {emb} vs {label}: {e}")

            # 3. LISI (Special Case: iLISI and cLISI)
            elif metric == 'lisi':
                # iLISI: needs batch
                if batch_keys:
                    for batch in batch_keys:
                        try:
                            # LISI run takes X and label.
                            # For iLISI, label is batch.
                            # We should run it on embeddings usually.
                            if embedding_keys:
                                for emb in embedding_keys:
                                    X_emb = adata.obsm[emb]
                                    pair_start = time.time()
                                    sig = inspect.signature(metric_module.run)
                                    kw = {}
                                    if 'n_jobs' in sig.parameters: kw['n_jobs'] = n_jobs
                                    if 'verbose' in sig.parameters: kw['verbose'] = False
                                    val = metric_module.run(X_emb, adata.obs[batch], **kw)
                                    pair_elapsed = time.time() - pair_start
                                    results.append({
                                        'Atlas Name': atlas_name,
                                        'Metric Name': 'iLISI',
                                        'Reference/Input 1': emb,
                                        'Prediction/Input 2': batch,
                                        'Value': val,
                                        'Time (s)': round(pair_elapsed, 3)
                                    })
                            else:
                                # Run on X
                                pair_start = time.time()
                                sig = inspect.signature(metric_module.run)
                                kw = {}
                                if 'n_jobs' in sig.parameters: kw['n_jobs'] = n_jobs
                                if 'verbose' in sig.parameters: kw['verbose'] = False
                                val = metric_module.run(adata.X, adata.obs[batch], **kw)
                                pair_elapsed = time.time() - pair_start
                                results.append({
                                    'Atlas Name': atlas_name,
                                    'Metric Name': 'iLISI',
                                    'Reference/Input 1': 'X',
                                    'Prediction/Input 2': batch,
                                    'Value': val,
                                    'Time (s)': round(pair_elapsed, 3)
                                })
                        except Exception as e:
                            logger.warning(f"Failed to calculate iLISI for {batch}: {e}")
                
                # cLISI: needs cell type (ref or pred)
                targets = list(set(ref_keys + pred_keys))
                for label in targets:
                    try:
                        if embedding_keys:
                            for emb in embedding_keys:
                                X_emb = adata.obsm[emb]
                                pair_start = time.time()
                                sig = inspect.signature(metric_module.run)
                                kw = {}
                                if 'n_jobs' in sig.parameters: kw['n_jobs'] = n_jobs
                                if 'verbose' in sig.parameters: kw['verbose'] = False
                                val = metric_module.run(X_emb, adata.obs[label], **kw)
                                pair_elapsed = time.time() - pair_start
                                results.append({
                                    'Atlas Name': atlas_name,
                                    'Metric Name': 'cLISI',
                                    'Reference/Input 1': emb,
                                    'Prediction/Input 2': label,
                                    'Value': val,
                                    'Time (s)': round(pair_elapsed, 3)
                                })
                        else:
                            pair_start = time.time()
                            sig = inspect.signature(metric_module.run)
                            kw = {}
                            if 'n_jobs' in sig.parameters: kw['n_jobs'] = n_jobs
                            if 'verbose' in sig.parameters: kw['verbose'] = False
                            val = metric_module.run(adata.X, adata.obs[label], **kw)
                            pair_elapsed = time.time() - pair_start
                            results.append({
                                'Atlas Name': atlas_name,
                                'Metric Name': 'cLISI',
                                'Reference/Input 1': 'X',
                                'Prediction/Input 2': label,
                                'Value': val,
                                'Time (s)': round(pair_elapsed, 3)
                            })
                    except Exception as e:
                        logger.warning(f"Failed to calculate cLISI for {label}: {e}")

            # 4. Batch Metrics (kBET, PCR) — evaluate per (embedding × batch)
            elif metric in batch_metrics:
                if not batch_keys:
                    continue
                for batch in batch_keys:
                    if embedding_keys:
                        for emb in embedding_keys:
                            try:
                                X_emb = adata.obsm[emb]
                                pair_start = time.time()
                                sig = inspect.signature(metric_module.run)
                                kw = {}
                                if 'n_jobs' in sig.parameters: kw['n_jobs'] = n_jobs
                                if 'verbose' in sig.parameters: kw['verbose'] = False
                                val = metric_module.run(X_emb, adata.obs[batch], **kw)
                                pair_elapsed = time.time() - pair_start
                                results.append({
                                    'Atlas Name': atlas_name,
                                    'Metric Name': metric,
                                    'Reference/Input 1': emb,
                                    'Prediction/Input 2': batch,
                                    'Value': val,
                                    'Time (s)': round(pair_elapsed, 3)
                                })
                            except Exception as e:
                                logger.warning(f"Failed to calculate {metric} for {emb} vs {batch}: {e}")
                    else:
                        try:
                            pair_start = time.time()
                            sig = inspect.signature(metric_module.run)
                            kw = {}
                            if 'n_jobs' in sig.parameters: kw['n_jobs'] = n_jobs
                            if 'verbose' in sig.parameters: kw['verbose'] = False
                            val = metric_module.run(adata.X, adata.obs[batch], **kw)
                            pair_elapsed = time.time() - pair_start
                            results.append({
                                'Atlas Name': atlas_name,
                                'Metric Name': metric,
                                'Reference/Input 1': 'X',
                                'Prediction/Input 2': batch,
                                'Value': val,
                                'Time (s)': round(pair_elapsed, 3)
                            })
                        except Exception as e:
                            logger.warning(f"Failed to calculate {metric} for {batch}: {e}")

            # 5. Graph Connectivity
            elif metric in graph_metrics:
                # We need embeddings AND labels
                targets = list(set(ref_keys + pred_keys))
                
                if embedding_keys:
                    for emb in embedding_keys:
                        # Calculate neighbors for this embedding
                        key_added = f'neighbors_{emb}'
                        try:
                            # Ensure neighbors are calculated
                            import scanpy as sc
                            sc.pp.neighbors(adata, use_rep=emb, key_added=key_added)
                        except Exception as e:
                            logger.warning(f"Failed to calculate neighbors for {emb}: {e}")
                            continue

                        for label in targets:
                            try:
                                pair_start = time.time()
                                sig = inspect.signature(metric_module.run)
                                kw = {'neighbors_key': key_added, 'label_key': label}
                                if 'n_jobs' in sig.parameters: kw['n_jobs'] = n_jobs
                                if 'verbose' in sig.parameters: kw['verbose'] = False
                                val = metric_module.run(adata, **kw)
                                pair_elapsed = time.time() - pair_start
                                results.append({
                                    'Atlas Name': atlas_name,
                                    'Metric Name': metric,
                                    'Reference/Input 1': emb,
                                    'Prediction/Input 2': label,
                                    'Value': val,
                                    'Time (s)': round(pair_elapsed, 3)
                                })
                            except Exception as e:
                                logger.warning(f"Failed to calculate {metric} for {emb} vs {label}: {e}")
                else:
                    # No embeddings found, use default neighbors (X or PCA)
                    # We still need labels
                    for label in targets:
                        try:
                            # metric_module.run will calculate neighbors if 'neighbors' key missing
                            pair_start = time.time()
                            sig = inspect.signature(metric_module.run)
                            kw = {'label_key': label}
                            if 'n_jobs' in sig.parameters: kw['n_jobs'] = n_jobs
                            if 'verbose' in sig.parameters: kw['verbose'] = False
                            val = metric_module.run(adata, **kw)
                            pair_elapsed = time.time() - pair_start
                            results.append({
                                'Atlas Name': atlas_name,
                                'Metric Name': metric,
                                'Reference/Input 1': 'Default',
                                'Prediction/Input 2': label,
                                'Value': val,
                                'Time (s)': round(pair_elapsed, 3)
                            })
                        except Exception as e:
                            logger.warning(f"Failed to calculate {metric} for Default vs {label}: {e}")

            # 6. Bio Conservation
            elif metric in bio_metrics:
                # Requires adata_before. 
                # If we don't have it, we can't run it properly.
                # We'll skip for now or log warning.
                logger.info(f"Skipping {metric} as it requires 'adata_before'.")

            else:
                logger.warning(f"Metric {metric} not categorized in pipeline.")

        except Exception as e:
            logger.error(f"Error running metric {metric}: {e}")
        
        # Calculate and display metric execution time
        metric_elapsed = time.time() - metric_start_time
        pbar.set_postfix_str(f"Time: {metric_elapsed:.2f}s", refresh=True)

    df = pd.DataFrame(results)

    # Save MultiQC-compatible wide format to file_dir if provided
    if not df.empty and file_dir is not None and atlas_name is not None:
        os.makedirs(file_dir, exist_ok=True)
        wide_df = _pivot_annot_to_wide(df, atlas_name)
        wide_path = os.path.join(file_dir, f"{atlas_name}.tsv")
        wide_df.to_csv(wide_path, sep='\t', index=False)
        logger.info("MultiQC-compatible annotation table saved to %s", wide_path)
    
    # Clear LISI and kBET kNN caches to free memory
    try:
        from .annot import lisi as _lisi
        _lisi._clear_knn_cache()
    except Exception:
        pass
    try:
        from .annot import kbet as _kbet
        _kbet._clear_knn_cache()
    except Exception:
        pass
    
    return df


def cal_cluster(adata, atlas_name=None, metric_list=None, all_metrics=True,
                file_dir=None, n_jobs=-1, verbose=True, seed=42):
    """
    Comprehensive clustering assessment pipeline.
    
    Calculates all clustering metrics to evaluate the quality of cluster assignments
    against embedding representations. Uses CheckAtlasColumnDetector to auto-detect
    embedding keys and cluster label columns.
    
    Args:
        adata (AnnData): Annotated data matrix.
        atlas_name (str): Name of the atlas for labeling results.
        metric_list (list, optional): List of metric names to calculate.
                           If provided, overrides ``all_metrics``.
        all_metrics (bool): If True, calculate all available cluster metrics.
                            If False, calculate a default subset.
                            Ignored if ``metric_list`` is provided.
        file_dir (str): Directory path where the results CSV will be saved.
                       If None, saves to current working directory.
        n_jobs (int): Number of parallel jobs (-1 = all cores).
        verbose (bool): Whether to print progress.
        seed (int): Random seed for reproducibility.
    
    Returns:
        pd.DataFrame: Results dataframe with columns:
                      [Atlas Name, Metric Name, Embedding, Label Key, Value, Time (s)]
    """
    import gc
    import inspect
    from scipy.sparse import issparse
    from ..utils.col_detector import CheckAtlasColumnDetector
    
    # Set file directory
    if file_dir is None:
        file_dir = os.getcwd()
    else:
        os.makedirs(file_dir, exist_ok=True)
    
    # Detect columns using CheckAtlasColumnDetector
    if verbose: print("Detecting embeddings and cluster labels...")
    detector = CheckAtlasColumnDetector(adata)
    params = detector.detect_all_parameters()
    
    embedding_keys = [x[0] for x in params['clustering']['embeddings']]
    label_keys = [x[0] for x in params['clustering']['cluster_labels']]
    
    if verbose:
        print(f"  Detected embeddings: {embedding_keys}")
        print(f"  Detected cluster labels: {label_keys}")
    
    if not embedding_keys:
        logger.warning("No embeddings detected in adata.obsm. Trying adata.X as fallback.")
        embedding_keys = ['X']
    
    if not label_keys:
        logger.warning("No cluster labels detected in adata.obs. Cannot run cluster metrics.")
        return pd.DataFrame()
    
    # Define metrics to run
    if metric_list is not None:
        metrics_list = [m for m in metric_list if m in METRICS_CLUST]
    elif all_metrics:
        metrics_list = METRICS_CLUST
    else:
        # Default subset: fastest and most commonly used
        metrics_list = [
            'silhouette',
            'davies_bouldin',
            'calinski_harabasz',
        ]
        metrics_list = [m for m in metrics_list if m in METRICS_CLUST]
    
    results = []
    
    # Total combinations
    total_combos = len(embedding_keys) * len(label_keys) * len(metrics_list)
    if verbose:
        print(f"\nRunning {len(metrics_list)} metrics × "
              f"{len(embedding_keys)} embeddings × "
              f"{len(label_keys)} label keys = {total_combos} calculations\n")
    
    # Create progress bar
    pbar = tqdm(total=total_combos, desc="Calculating Cluster Metrics",
                bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]',
                disable=not verbose)
    
    for emb_key in embedding_keys:
        # Extract embedding data
        if emb_key == 'X':
            X_emb = adata.X
            if issparse(X_emb):
                X_emb = X_emb.toarray()
        else:
            if emb_key not in adata.obsm:
                logger.warning(f"Embedding '{emb_key}' not found in adata.obsm. Skipping.")
                pbar.update(len(label_keys) * len(metrics_list))
                continue
            X_emb = np.asarray(adata.obsm[emb_key])
        
        for label_key in label_keys:
            if label_key not in adata.obs.columns:
                logger.warning(f"Label key '{label_key}' not found in adata.obs. Skipping.")
                pbar.update(len(metrics_list))
                continue
            
            labels = np.asarray(adata.obs[label_key])
            
            # Check we have at least 2 clusters
            n_unique = len(np.unique(labels))
            if n_unique < 2:
                logger.warning(f"Label '{label_key}' has < 2 clusters. Skipping.")
                pbar.update(len(metrics_list))
                continue
            
            for metric_name in metrics_list:
                metric_start_time = time.time()
                pbar.set_description(f"{emb_key}|{label_key}: {metric_name}")
                
                try:
                    metric_module = getattr(cluster, metric_name)
                    
                    # Build kwargs dynamically based on metric's signature
                    sig = inspect.signature(metric_module.run)
                    metric_params = sig.parameters
                    
                    kwargs = {}
                    
                    if 'n_jobs' in metric_params:
                        kwargs['n_jobs'] = n_jobs
                    if 'verbose' in metric_params:
                        kwargs['verbose'] = False
                    if 'random_state' in metric_params:
                        kwargs['random_state'] = seed
                    if 'seed' in metric_params:
                        kwargs['seed'] = seed
                    if 'max_samples' in metric_params:
                        kwargs['max_samples'] = None  # disable subsampling
                    
                    # Call the metric
                    value = metric_module.run(X_emb, labels, **kwargs)
                    
                    elapsed = time.time() - metric_start_time
                    
                    results.append({
                        'Atlas Name': atlas_name,
                        'Metric Name': metric_name,
                        'Embedding': emb_key,
                        'Label Key': label_key,
                        'Value': value,
                        'Time (s)': round(elapsed, 3)
                    })
                    
                    pbar.set_postfix_str(f"={value:.4f} ({elapsed:.1f}s)", refresh=True)
                    
                except Exception as e:
                    elapsed = time.time() - metric_start_time
                    logger.warning(f"Failed {metric_name} for {emb_key}/{label_key}: {e}")
                    results.append({
                        'Atlas Name': atlas_name,
                        'Metric Name': metric_name,
                        'Embedding': emb_key,
                        'Label Key': label_key,
                        'Value': np.nan,
                        'Time (s)': round(elapsed, 3)
                    })
                
                pbar.update(1)
    
    pbar.close()
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Summary
    if verbose and not df.empty:
        total_time = df['Time (s)'].sum()
        print(f"\nTotal computation time: {total_time:.2f}s")
        print(f"Results: {len(df)} measurements across {len(df['Metric Name'].unique())} metrics")
    
    gc.collect()
    return df


def calc_metric_cluster_scanpy(
    metric, adata, obs_key, obsm_key_representation
):
    if metric in METRICS_CLUST:
        start_time = time.time()
        logger.debug(f"Start {metric} calc")
        metric_module = getattr(cluster, metric)
        annotations = adata.obs[obs_key]
        if obsm_key_representation != "":
            count_repr = adata.obsm[obsm_key_representation]
            metric_value = metric_module.run(count_repr, annotations)
            running_time = time.time() - start_time
            logger.debug(f"{metric} calc finished, duration {running_time}")
            return metric_value, running_time

        else:
            original_count = adata.X.toarray()
            metric_value = metric_module.run(original_count, annotations)
            running_time = time.time() - start_time
            logger.debug(f"{metric} calc finished, duration {running_time}")
            return metric_value, running_time
    else:
        logger.warning(
            f"{metric} is not a recognized "
            f"cluster metric.\n"
            f"List of clustering metrics: {METRICS_CLUST}"
        )
        return -1


def calc_metric_cluster_seurat(
    metric, seurat, obs_key, obsm_key_representation
):
    if metric in METRICS_CLUST:
        start_time = time.time()
        logger.debug(f"Start {metric} calc")
        metric_module = getattr(cluster, metric)
        annotations = ro.conversion.rpy2py(R_ANNOT(seurat, obs_key))
        count_repr = ro.conversion.rpy2py(
            R_REDUCTION(seurat, obsm_key_representation)
        )
        metric_value = metric_module.run(count_repr, annotations)
        running_time = time.time() - start_time
        logger.debug(f"{metric} calc finished, duration {running_time}")
        return metric_value, running_time
    else:
        logger.warning(
            f"{metric} is not a recognized "
            f"cluster metric.\n"
            f"List of clustering metrics: {METRICS_CLUST}"
        )
        return -1


def calc_metric_annot_scanpy(metric, adata, obs_key, ref_obs):
    if metric in METRICS_ANNOT:
        start_time = time.time()
        logger.debug(f"Start {metric} calc")
        metric_module = getattr(annot, metric)
        annotation = adata.obs[obs_key]
        ref_annotation = adata.obs[ref_obs]
        annotation, ref_annotation = annotation_to_num(
            annotation, ref_annotation
        )
        metric_value = metric_module.run(annotation, ref_annotation)
        running_time = time.time() - start_time
        logger.debug(f"{metric} calc finished, duration {running_time}")
        return metric_value, running_time
    else:
        logger.warning(
            f"{metric} is not a recognized annotation metric."
            f"\nList of annotation metrics: {METRICS_ANNOT}"
        )
        return -1


def calc_metric_annot_seurat(metric, seurat, obs_key, ref_obs):
    if metric in METRICS_ANNOT:
        start_time = time.time()
        logger.debug(f"Start {metric} calc")
        metric_module = getattr(annot, metric)
        annotation = ro.conversion.rpy2py(R_ANNOT(seurat, obs_key))
        ref_annotation = ro.conversion.rpy2py(R_ANNOT(seurat, ref_obs))
        # annotation, ref_annotation = annotation_to_num(
        #    annotation, ref_annotation
        # )
        metric_value = metric_module.run(annotation, ref_annotation)
        running_time = time.time() - start_time
        logger.debug(f"{metric} calc finished, duration {running_time}")
        return metric_value, running_time
    else:
        logger.warning(
            f"{metric} is not a recognized annotation metric."
            f"\nList of annotation metrics: {METRICS_ANNOT}"
        )
        return -1


def cal_dimred(adata, atlas_name=None, low_dim_keys=None, high_dim_key='X',
               metric_list=None, k_neighbors=30, n_samples=None, seed=42,
               n_jobs=-1, file_dir=None, verbose=True):
    """Calculate dimensionality reduction metrics for multiple embeddings.

    For each ``low_dim_key`` the embedding is compared against the
    reference ``high_dim_key`` (defaults to ``adata.X`` — the raw gene
    expression matrix).

    Precomputation (distance matrices + kNN graphs) is performed once for
    the high‑dimensional reference and once per low‑dimensional embedding
    key.  Memory‑mapped files are used for N×N distance matrices on
    datasets with > 10 000 cells to keep RAM usage bounded.

    Parameters
    ----------
    adata : AnnData
    atlas_name : str, optional
        Used only for logging.
    low_dim_keys : list of str, optional
        ``.obsm`` keys to evaluate.  Defaults to all keys returned by
        :meth:`CheckAtlasColumnDetector.get_dimred_embeddings`.
    high_dim_key : str
        Reference key.  ``'X'`` means ``adata.X``.  (Default ``'X'``)
    metric_list : list of str, optional
        Metrics to compute.  Defaults to :data:`METRICS_DIMRED`.
    k_neighbors : int
    n_samples : int or None
        Number of cells to subsample.  ``None`` = use all cells.
    seed : int
    n_jobs : int
    file_dir : str, optional
        Directory for temporary memmap files.  Falls back to system
        ``/tmp``.
    verbose : bool

    Returns
    -------
    pd.DataFrame
        Long‑format table with columns
        ``[Metric Name, Low Dim Key, High Dim Key, Value, Time (s)]``.
        Caller is responsible for pivoting and saving.
    """
    import gc
    import inspect
    import tempfile
    import uuid

    from sklearn.metrics import pairwise_distances
    from sklearn.neighbors import NearestNeighbors

    from ..utils.col_detector import CheckAtlasColumnDetector

    if metric_list is None:
        metric_list = METRICS_DIMRED
    metric_list = [m for m in metric_list if m in METRICS_DIMRED]
    if not metric_list:
        return pd.DataFrame()

    if low_dim_keys is None:
        detector = CheckAtlasColumnDetector(adata)
        low_dim_keys = detector.get_dimred_embeddings()

    low_dim_keys = [k for k in low_dim_keys if k != high_dim_key]
    if not low_dim_keys:
        if verbose:
            print(f"  No embeddings to compare "
                  f"(only key is {high_dim_key}).")
        return pd.DataFrame()

    n_obs = adata.n_obs
    if n_samples is not None and n_samples < n_obs:
        np.random.seed(seed)
        sample_indices = np.random.choice(n_obs, n_samples, replace=False)
        n_cells = n_samples
    else:
        sample_indices = np.arange(n_obs)
        n_cells = n_obs

    if high_dim_key == 'X':
        high_dim_data = adata.X[sample_indices]
        if hasattr(high_dim_data, "toarray"):
            high_dim_data = high_dim_data.toarray()
    else:
        if high_dim_key not in adata.obsm_keys():
            raise ValueError(
                f"High-dim key '{high_dim_key}' not found in adata.obsm.")
        high_dim_data = adata.obsm[high_dim_key][sample_indices]

    high_n_features = high_dim_data.shape[1]

    use_memmap = n_cells > 10000
    all_memmap_files = []

    if file_dir:
        temp_dir = file_dir
    else:
        temp_dir = os.path.join(os.path.expanduser("~"), ".checkatlas", "tmp")
    os.makedirs(temp_dir, exist_ok=True)

    run_id = str(uuid.uuid4())[:8]
    chunk_size = min(1000, n_cells)

    if verbose:
        print(f"  Reference: {high_dim_key}  "
              f"({n_cells:,} cells  ×  {high_n_features} features)")
        print(f"  Embeddings: {low_dim_keys}")
        print(f"  Metrics: {metric_list}")

    # ── High‑dim kNN ──
    nbrs_high = NearestNeighbors(
        n_neighbors=k_neighbors + 1, n_jobs=n_jobs
    ).fit(high_dim_data)
    high_knn_dists, high_knn_indices = nbrs_high.kneighbors(high_dim_data)

    # ── High‑dim distance matrix (only when distance‑based metrics
    #     are in the list) ──
    _DIST_METRICS = frozenset(("kruskal_stress", "spearman_rho", "dCor",
                               "trustworthiness", "continuity"))
    need_high_dists = bool(set(metric_list) & _DIST_METRICS)
    high_dim_dists = None

    if need_high_dists:
        if verbose:
            print("  Precomputing high-dim distance matrix "
                  f"({n_cells:,}×{n_cells:,})…")
        if use_memmap:
            high_dists_path = os.path.join(
                temp_dir, f"high_dists_{run_id}.dat")
            all_memmap_files.append(high_dists_path)
            high_dim_dists = np.memmap(
                high_dists_path, dtype="float32", mode="w+",
                shape=(n_cells, n_cells))
        else:
            high_dim_dists = np.zeros(
                (n_cells, n_cells), dtype=np.float32)

        for i in tqdm(range(0, n_cells, chunk_size),
                      desc="High-Dim Distances", disable=not verbose):
            end = min(i + chunk_size, n_cells)
            high_dim_dists[i:end, :] = pairwise_distances(
                high_dim_data[i:end], high_dim_data, n_jobs=n_jobs)
            if use_memmap:
                high_dim_dists.flush()
        gc.collect()

    total_combos = len(low_dim_keys) * len(metric_list)
    if verbose:
        print(
            f"\n  {len(metric_list)} metrics × {len(low_dim_keys)} embeddings"
            f" = {total_combos} calculations\n")

    results = []
    pbar = tqdm(
        total=total_combos,
        desc="Calculating Dimred Metrics",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} "
                   "[{elapsed}<{remaining}]",
        disable=not verbose)

    for low_dim_key in low_dim_keys:
        low_dim_data = adata.obsm[low_dim_key][sample_indices]

        nbrs_low = NearestNeighbors(
            n_neighbors=k_neighbors + 1, n_jobs=n_jobs
        ).fit(low_dim_data)
        low_knn_dists, low_knn_indices = nbrs_low.kneighbors(low_dim_data)

        low_dim_dists = None
        low_dists_path = None
        need_low_dists = bool(set(metric_list) & _DIST_METRICS)

        if need_low_dists:
            if use_memmap:
                low_dists_path = os.path.join(
                    temp_dir,
                    f"low_dists_{low_dim_key.replace('/', '_')}_{run_id}.dat")
                all_memmap_files.append(low_dists_path)
                low_dim_dists = np.memmap(
                    low_dists_path, dtype="float32", mode="w+",
                    shape=(n_cells, n_cells))
            else:
                low_dim_dists = np.zeros(
                    (n_cells, n_cells), dtype=np.float32)

            for i in tqdm(range(0, n_cells, chunk_size),
                          desc=f"Low-Dim Distances ({low_dim_key})",
                          disable=not verbose):
                end = min(i + chunk_size, n_cells)
                low_dim_dists[i:end, :] = pairwise_distances(
                    low_dim_data[i:end], low_dim_data, n_jobs=n_jobs)
                if use_memmap:
                    low_dim_dists.flush()
            gc.collect()

        for metric_name in metric_list:
            pbar.set_description(f"{low_dim_key}: {metric_name}")
            t_start = time.time()

            try:
                metric_module = getattr(dimred, metric_name)
                sig = inspect.signature(metric_module.run)
                params = sig.parameters

                kwargs = {}
                if "low_dim_key" in params:
                    kwargs["low_dim_key"] = low_dim_key
                if "high_dim_key" in params:
                    kwargs["high_dim_key"] = high_dim_key
                if "k_neighbors" in params:
                    kwargs["k_neighbors"] = k_neighbors
                if "n_samples" in params:
                    kwargs["n_samples"] = n_samples
                if "n_jobs" in params:
                    kwargs["n_jobs"] = n_jobs
                if "seed" in params:
                    kwargs["seed"] = seed
                if "verbose" in params:
                    kwargs["verbose"] = False

                if "precomputed_high_knn" in params:
                    kwargs["precomputed_high_knn"] = high_knn_indices
                if "precomputed_low_knn" in params:
                    kwargs["precomputed_low_knn"] = low_knn_indices
                if "precomputed_high_knn_dists" in params:
                    kwargs["precomputed_high_knn_dists"] = high_knn_dists
                if "precomputed_low_knn_dists" in params:
                    kwargs["precomputed_low_knn_dists"] = low_knn_dists
                if "precomputed_high_dists" in params:
                    kwargs["precomputed_high_dists"] = high_dim_dists
                if "precomputed_low_dists" in params:
                    kwargs["precomputed_low_dists"] = low_dim_dists

                value = metric_module.run(adata, **kwargs)
                elapsed = round(time.time() - t_start, 3)

                results.append({
                    "Metric Name": metric_name,
                    "Low Dim Key": low_dim_key,
                    "High Dim Key": high_dim_key,
                    "Value": value,
                    "Time (s)": elapsed,
                })
            except Exception as e:
                elapsed = round(time.time() - t_start, 3)
                logger.warning(
                    "Failed %s on %s: %s", metric_name, low_dim_key, e)
                results.append({
                    "Metric Name": metric_name,
                    "Low Dim Key": low_dim_key,
                    "High Dim Key": high_dim_key,
                    "Value": np.nan,
                    "Time (s)": elapsed,
                })

            pbar.update(1)

        if low_dim_dists is not None:
            try:
                del low_dim_dists
            except Exception:
                pass
            if low_dists_path and os.path.exists(low_dists_path):
                try:
                    os.remove(low_dists_path)
                except Exception:
                    pass
                if low_dists_path in all_memmap_files:
                    all_memmap_files.remove(low_dists_path)
        gc.collect()

    pbar.close()

    if high_dim_dists is not None:
        try:
            del high_dim_dists
        except Exception:
            pass
        gc.collect()

    for fpath in all_memmap_files:
        try:
            if os.path.exists(fpath):
                os.remove(fpath)
        except Exception as exc:
            logger.debug(
                "Failed to clean up temp file %s: %s", fpath, exc)

    if "_dimred_cache" in adata.uns:
        del adata.uns["_dimred_cache"]

    gc.collect()

    df = pd.DataFrame(results)

    if verbose and not df.empty:
        total_time = df["Time (s)"].sum()
        n_results = len(df)
        print(f"\nTotal computation time: {total_time:.2f}s")
        print(f"Results: {n_results} measurements across "
              f"{len(metric_list)} metrics")

    return df


def calc_metric_dimred(metric, adata, obsm_key):
    """
    Calculate a single dimensionality reduction metric (legacy function).
    
    For comprehensive assessment, use cal_dimred() instead.
    """
    if metric in METRICS_DIMRED:
        start_time = time.time()
        logger.debug(f"Start {metric} calc")
        metric_module = getattr(dimred, metric)
        high_dim_counts = adata.X
        low_dim_counts = adata.obsm[obsm_key]
        metric_value = metric_module.run(high_dim_counts, low_dim_counts)
        running_time = time.time() - start_time
        logger.debug(f"{metric} calc finished, duration {running_time}")
        return metric_value, running_time
    else:
        logger.warning(
            f"{metric} is not a recognized "
            f"dimensionality reduction metric."
            f"\nList of dim. red. metrics: {METRICS_DIMRED}"
        )
        return -1


def run_all_metrics(adata=None, atlas_path=None, atlas_name=None,
                    file_dir=None, n_jobs=-1, verbose=True, seed=42,
                    # Annotation params
                    run_annotation=True, all_annot_metrics=True,
                    # Clustering params
                    run_clustering=True, all_cluster_metrics=True,
                    # Dimred params
                    run_dimred=True, all_dimred_metrics=True,
                    low_dim_key='X_umap', high_dim_key='X',
                    k_neighbors=30, n_samples=None,
                    use_memmap=True, temp_dir=None, chunk_size=1000):
    """
    Unified metrics pipeline for CheckAtlas.
    
    The standalone entry point that runs ALL metric tasks (annotation, clustering,
    dimensionality reduction) on a single-cell atlas and produces a consolidated
    results CSV. Designed for biologists and bioinformaticians — just provide
    an atlas path and get comprehensive quality metrics.
    
    Usage:
        >>> from checkatlas.metrics.metrics import run_all_metrics
        >>> results = run_all_metrics(atlas_path="my_atlas.h5ad")
        
        >>> # Or with a preloaded AnnData:
        >>> results = run_all_metrics(adata=my_adata, atlas_name="my_atlas")
    
    Args:
        adata (AnnData, optional): Preloaded AnnData object. If None, loads from atlas_path.
        atlas_path (str, optional): Path to .h5ad atlas file. Used if adata is None.
        atlas_name (str, optional): Name for the atlas (used in output). If None,
                                   derived from atlas_path filename.
        file_dir (str, optional): Directory for output CSV. Defaults to current directory.
        n_jobs (int): Number of parallel jobs (-1 = all cores).
        verbose (bool): Whether to print progress and summary.
        seed (int): Random seed for reproducibility.
        
        run_annotation (bool): Whether to run annotation metrics.
        all_annot_metrics (bool): If True, run all annotation metrics.
        
        run_clustering (bool): Whether to run clustering metrics.
        all_cluster_metrics (bool): If True, run all clustering metrics.
        
        run_dimred (bool): Whether to run dimensionality reduction metrics.
        all_dimred_metrics (bool): If True, run all dimred metrics.
        low_dim_key (str): Key in adata.obsm for low-dimensional embedding.
        high_dim_key (str): Key for high-dimensional data ('X' = adata.X).
        k_neighbors (int): Number of neighbors for kNN-based metrics.
        n_samples (int, optional): Number of samples for dimred metrics.
        use_memmap (bool): Whether to use memory-mapped files for large distance matrices.
        temp_dir (str, optional): Directory for temporary memmap files.
        chunk_size (int): Chunk size for batched distance computation.
    
    Returns:
        pd.DataFrame: Consolidated results with columns:
                      [Atlas Name, Task, Metric Name, Input 1, Input 2, Value, Time (s)]
    """
    import gc
    import scanpy as sc
    
    overall_start = time.time()
    
    # ────────────────────────────────────────────────────────────────────
    # 1. Load atlas
    # ────────────────────────────────────────────────────────────────────
    if adata is None and atlas_path is not None:
        if verbose: print(f"Loading atlas from: {atlas_path}")
        adata = sc.read_h5ad(atlas_path)
        if atlas_name is None:
            atlas_name = os.path.splitext(os.path.basename(atlas_path))[0]
    elif adata is not None:
        if atlas_name is None:
            atlas_name = "atlas"
    else:
        raise ValueError("Provide either `adata` (AnnData) or `atlas_path` (str).")
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"  CheckAtlas — Unified Metrics Pipeline")
        print(f"  Atlas: {atlas_name}")
        print(f"  Cells: {adata.n_obs:,}  |  Genes: {adata.n_vars:,}")
        print(f"  Tasks: ", end="")
        tasks = []
        if run_annotation: tasks.append("Annotation")
        if run_clustering: tasks.append("Clustering")
        if run_dimred: tasks.append("Dimred")
        print(" → ".join(tasks))
        print(f"{'='*70}\n")
    
    # Set output directory
    if file_dir is None:
        file_dir = os.getcwd()
    else:
        os.makedirs(file_dir, exist_ok=True)
    
    # Consolidated results
    all_results = []
    task_summary = {}
    
    # ────────────────────────────────────────────────────────────────────
    # 2. Run Annotation Metrics
    # ────────────────────────────────────────────────────────────────────
    if run_annotation:
        task_start = time.time()
        if verbose:
            print(f"\n{'─'*50}")
            print(f"  TASK 1/3: Annotation Metrics")
            print(f"{'─'*50}")
        
        try:
            df_annot = cal_annot(
                adata, atlas_name=atlas_name,
                all=all_annot_metrics,
                file_dir=file_dir
            )
            
            if not df_annot.empty:
                # Normalize columns to unified schema
                unified = pd.DataFrame({
                    'Atlas Name': df_annot['Atlas Name'],
                    'Task': 'annotation',
                    'Metric Name': df_annot['Metric Name'],
                    'Input 1': df_annot['Reference/Input 1'],
                    'Input 2': df_annot['Prediction/Input 2'],
                    'Value': df_annot['Value'],
                    'Time (s)': df_annot.get('Time (s)', np.nan)
                })
                all_results.append(unified)
            
            task_elapsed = time.time() - task_start
            n_metrics = len(df_annot) if not df_annot.empty else 0
            task_summary['Annotation'] = {
                'count': n_metrics,
                'time': round(task_elapsed, 2)
            }
            if verbose:
                print(f"  ✓ Annotation: {n_metrics} results in {task_elapsed:.1f}s")
        
        except Exception as e:
            task_elapsed = time.time() - task_start
            logger.error(f"Annotation pipeline failed: {e}")
            task_summary['Annotation'] = {'count': 0, 'time': round(task_elapsed, 2), 'error': str(e)}
            if verbose:
                print(f"  ✗ Annotation failed: {e}")
    
    # ────────────────────────────────────────────────────────────────────
    # 3. Run Clustering Metrics
    # ────────────────────────────────────────────────────────────────────
    if run_clustering:
        task_start = time.time()
        if verbose:
            print(f"\n{'─'*50}")
            print(f"  TASK 2/3: Clustering Metrics")
            print(f"{'─'*50}")
        
        try:
            df_cluster = cal_cluster(
                adata, atlas_name=atlas_name,
                all_metrics=all_cluster_metrics,
                file_dir=file_dir,
                n_jobs=n_jobs,
                verbose=verbose,
                seed=seed
            )
            
            if not df_cluster.empty:
                unified = pd.DataFrame({
                    'Atlas Name': df_cluster['Atlas Name'],
                    'Task': 'clustering',
                    'Metric Name': df_cluster['Metric Name'],
                    'Input 1': df_cluster['Embedding'],
                    'Input 2': df_cluster['Label Key'],
                    'Value': df_cluster['Value'],
                    'Time (s)': df_cluster['Time (s)']
                })
                all_results.append(unified)
            
            task_elapsed = time.time() - task_start
            n_metrics = len(df_cluster) if not df_cluster.empty else 0
            task_summary['Clustering'] = {
                'count': n_metrics,
                'time': round(task_elapsed, 2)
            }
            if verbose:
                print(f"  ✓ Clustering: {n_metrics} results in {task_elapsed:.1f}s")
        
        except Exception as e:
            task_elapsed = time.time() - task_start
            logger.error(f"Clustering pipeline failed: {e}")
            task_summary['Clustering'] = {'count': 0, 'time': round(task_elapsed, 2), 'error': str(e)}
            if verbose:
                print(f"  ✗ Clustering failed: {e}")
    
    # ────────────────────────────────────────────────────────────────────
    # 4. Run Dimensionality Reduction Metrics
    # ────────────────────────────────────────────────────────────────────
    if run_dimred:
        task_start = time.time()
        if verbose:
            print(f"\n{'─'*50}")
            print(f"  TASK 3/3: Dimensionality Reduction Metrics")
            print(f"{'─'*50}")
        
        try:
            metric_list = METRICS_DIMRED if all_dimred_metrics else None
            df_dimred = cal_dimred(
                adata, atlas_name=atlas_name,
                high_dim_key=high_dim_key,
                metric_list=metric_list,
                k_neighbors=k_neighbors,
                n_samples=n_samples,
                seed=seed,
                n_jobs=n_jobs,
                file_dir=file_dir,
                verbose=verbose,
            )
            
            if not df_dimred.empty:
                unified = pd.DataFrame({
                    'Atlas Name': atlas_name,
                    'Task': 'dimred',
                    'Metric Name': df_dimred['Metric Name'],
                    'Input 1': df_dimred['Low Dim Key'],
                    'Input 2': df_dimred['High Dim Key'],
                    'Value': df_dimred['Value'],
                    'Time (s)': df_dimred['Time (s)']
                })
                all_results.append(unified)
            
            task_elapsed = time.time() - task_start
            n_metrics = len(df_dimred) if not df_dimred.empty else 0
            task_summary['Dimred'] = {
                'count': n_metrics,
                'time': round(task_elapsed, 2)
            }
            if verbose:
                print(f"  ✓ Dimred: {n_metrics} results in {task_elapsed:.1f}s")
        
        except Exception as e:
            task_elapsed = time.time() - task_start
            logger.error(f"Dimred pipeline failed: {e}")
            task_summary['Dimred'] = {'count': 0, 'time': round(task_elapsed, 2), 'error': str(e)}
            if verbose:
                print(f"  ✗ Dimred failed: {e}")
    
    # ────────────────────────────────────────────────────────────────────
    # 5. Consolidate & Save
    # ────────────────────────────────────────────────────────────────────
    if all_results:
        df_all = pd.concat(all_results, ignore_index=True)
    else:
        df_all = pd.DataFrame(columns=[
            'Atlas Name', 'Task', 'Metric Name', 'Input 1', 'Input 2', 'Value', 'Time (s)'
        ])
    
    # Save unified CSV
    if not df_all.empty:
        filename = os.path.join(file_dir, f"checkatlas_all_metrics_{atlas_name}.csv")
        df_all.to_csv(filename, index=False)
        if verbose:
            print(f"\nSaved unified results to: {filename}")
        logger.info(f"Saved unified metrics to {filename}")
    
    # ────────────────────────────────────────────────────────────────────
    # 6. Print Summary
    # ────────────────────────────────────────────────────────────────────
    overall_elapsed = time.time() - overall_start
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"  CheckAtlas — Results Summary")
        print(f"{'='*70}")
        print(f"  {'Task':<20} {'Metrics':>10} {'Time (s)':>12} {'Status':>10}")
        print(f"  {'─'*52}")
        
        for task_name, info in task_summary.items():
            status = "✗ Error" if 'error' in info else "✓ Done"
            print(f"  {task_name:<20} {info['count']:>10} {info['time']:>12.2f} {status:>10}")
        
        print(f"  {'─'*52}")
        total_count = sum(info['count'] for info in task_summary.values())
        print(f"  {'TOTAL':<20} {total_count:>10} {overall_elapsed:>12.2f}")
        print(f"{'='*70}")
        
        if not df_all.empty:
            print(f"\n  Unique metrics computed: {df_all['Metric Name'].nunique()}")
            print(f"  Total measurements: {len(df_all)}")
            
            # Top 5 slowest metrics
            if 'Time (s)' in df_all.columns:
                slowest = df_all.nlargest(5, 'Time (s)')[['Task', 'Metric Name', 'Time (s)']]
                if not slowest.empty:
                    print(f"\n  Slowest metrics:")
                    for _, row in slowest.iterrows():
                        print(f"    {row['Task']}/{row['Metric Name']}: {row['Time (s)']:.3f}s")
        
        print()
    
    gc.collect()
    return df_all


def annotation_to_num(annotation, ref_annotation):

    """
    Transforms the annotations from categorical to numerical

    Parameters
    ----------
    adata
    partition_key
    reference

    Returns
    -------

    """
    annotation = annotation.to_numpy()
    ref_annotation = ref_annotation.to_numpy()
    le = LabelEncoder()
    le.fit(annotation)
    annotation = le.transform(annotation)
    le2 = LabelEncoder()
    le2.fit(ref_annotation)
    ref_annotation = le2.transform(ref_annotation)
    return annotation, ref_annotation
