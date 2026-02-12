import logging
import time
import pandas as pd
import numpy as np
from tqdm import tqdm
import os

import rpy2.robjects as ro
import rpy2.robjects as robjects
from sklearn.preprocessing import LabelEncoder

from . import annot, cluster, dimred
# Import CheckAtlasColumnDetector inside function or at top if no circular dependency
# from ..atlas import CheckAtlasColumnDetector # This might cause circular import if atlas imports metrics

METRICS_CLUST = cluster.__all__
METRICS_ANNOT = annot.__all__
METRICS_DIMRED = dimred.__all__

logger = logging.getLogger("checkatlas")

R_ANNOT = robjects.r(
    "type <- function(seurat, obs_key){ "
    "return(seurat[[obs_key]][[obs_key]])}"
)
R_REDUCTION = robjects.r(
    "reduc <- function(seurat, obsm_key){"
    " return(Embeddings(object = seurat, "
    "reduction = obsm_key))}"
)

def cal_annot(adata, atlas_name=None, all=False, file_dir=None):
    """
    Comprehensive annotation pipeline for all annotation metrics.
    
    Args:
        adata (AnnData): Annotated data matrix.
        all (bool): If True, calculate all available annotation metrics. 
                    If False, calculate a default subset.
        file_dir (str): Directory path where the results CSV will be saved.
                       If None, saves to current working directory.
                    
    Returns:
        pd.DataFrame: Results dataframe with columns:
                      [Atlas Name, Metric Name, Reference/Input 1, Prediction/Input 2, Value]
    """
    # Import here to avoid circular dependency
    from ..atlas import CheckAtlasColumnDetector
    
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
    
    # Detect batch keys (heuristic: look for 'batch' in metadata or column name)
    # CheckAtlasColumnDetector identifies metadata, but doesn't explicitly label 'batch'.
    # We'll look for columns containing 'batch' or use a default list if provided in adata.uns
    batch_keys = [col for col in adata.obs.columns if 'batch' in col.lower()]
    
    # Define metrics to run
    if all:
        metrics_list = METRICS_ANNOT
    else:
        # Default subset of metrics
        ## default is only ARI
        metrics_list = [
            'adj_rand_index', 
            'normalized_mutual_info', 
            'adj_mutual_info',
            # 'lisi',
            # 'kbet'
        ]
        # Filter to ensure they exist in METRICS_ANNOT
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
                            
                            val = metric_module.run(l_pred, l_true)
                            results.append({
                                'Atlas Name': atlas_name,
                                'Metric Name': metric,
                                'Reference/Input 1': ref,
                                'Prediction/Input 2': pred,
                                'Value': val
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
                            val = metric_module.run(X_emb, labels)
                            results.append({
                                'Atlas Name': atlas_name,
                                'Metric Name': metric,
                                'Reference/Input 1': emb,
                                'Prediction/Input 2': label,
                                'Value': val
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
                                    val = metric_module.run(X_emb, adata.obs[batch])
                                    results.append({
                                        'Atlas Name': atlas_name,
                                        'Metric Name': 'iLISI',
                                        'Reference/Input 1': emb,
                                        'Prediction/Input 2': batch,
                                        'Value': val
                                    })
                            else:
                                # Run on X
                                val = metric_module.run(adata.X, adata.obs[batch])
                                results.append({
                                    'Atlas Name': atlas_name,
                                    'Metric Name': 'iLISI',
                                    'Reference/Input 1': 'X',
                                    'Prediction/Input 2': batch,
                                    'Value': val
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
                                val = metric_module.run(X_emb, adata.obs[label])
                                results.append({
                                    'Atlas Name': atlas_name,
                                    'Metric Name': 'cLISI',
                                    'Reference/Input 1': emb,
                                    'Prediction/Input 2': label,
                                    'Value': val
                                })
                        else:
                            val = metric_module.run(adata.X, adata.obs[label])
                            results.append({
                                'Atlas Name': atlas_name,
                                'Metric Name': 'cLISI',
                                'Reference/Input 1': 'X',
                                'Prediction/Input 2': label,
                                'Value': val
                            })
                    except Exception as e:
                        logger.warning(f"Failed to calculate cLISI for {label}: {e}")

            # 4. Batch Metrics (kBET, PCR)
            elif metric in batch_metrics:
                if not batch_keys:
                    continue
                for batch in batch_keys:
                    try:
                        # kBET and PCR take adata and batch_label
                        # They might use X or embedding internally.
                        # kBET usually uses X or embedding.
                        # My updated kBET uses adata.X.
                        # PCR uses adata.X.
                        # If they support embeddings, we should loop embeddings?
                        # The signatures I saw: run(adata, batch_label, ...)
                        # They use adata.X by default.
                        # If we want to test embeddings, we might need to swap adata.X or pass embedding?
                        # For now, run on default adata.X (or whatever run uses).
                        val = metric_module.run(adata, batch_label=batch)
                        results.append({
                            'Atlas Name': atlas_name,
                            'Metric Name': metric,
                            'Reference/Input 1': 'adata',
                            'Prediction/Input 2': batch,
                            'Value': val
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
                                val = metric_module.run(adata, neighbors_key=key_added, label_key=label)
                                results.append({
                                    'Atlas Name': atlas_name,
                                    'Metric Name': metric,
                                    'Reference/Input 1': emb,
                                    'Prediction/Input 2': label,
                                    'Value': val
                                })
                            except Exception as e:
                                logger.warning(f"Failed to calculate {metric} for {emb} vs {label}: {e}")
                else:
                    # No embeddings found, use default neighbors (X or PCA)
                    # We still need labels
                    for label in targets:
                        try:
                            # metric_module.run will calculate neighbors if 'neighbors' key missing
                            val = metric_module.run(adata, label_key=label) 
                            results.append({
                                'Atlas Name': atlas_name,
                                'Metric Name': metric,
                                'Reference/Input 1': 'Default',
                                'Prediction/Input 2': label,
                                'Value': val
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
    
    # Save to CSV if results exist
    if not df.empty:
        # Define a filename with full path
        filename = os.path.join(file_dir, f"checkatlas_annotation_metrics_{atlas_name}.csv")
        df.to_csv(filename, index=False)
        logger.info(f"Saved annotation metrics to {filename}")
        
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


def cal_dimred(adata, atlas_name=None, low_dim_key='X_umap', high_dim_key='X',
               k_neighbors=30, n_samples=None, seed=42, n_jobs=-1, 
               all_metrics=True, file_dir=None, verbose=True, use_memmap=True,
               temp_dir=None, chunk_size=1000):
    """
    Comprehensive dimensionality reduction assessment pipeline.
    
    Calculates all dimred metrics to assess the quality of embeddings (UMAP, t-SNE, etc.)
    compared to the high-dimensional reference (adata.X).
    
    Optimizations:
    1. Precomputes pairwise distances and kNN graphs for both high and low dim spaces.
    2. Uses memory-mapped files for distance matrices to avoid RAM overload.
    3. Caches kNN results in adata.uns['_dimred_cache'] for efficiency.
    
    Args:
        adata (AnnData): Annotated data matrix.
        atlas_name (str): Name of the atlas for labeling results.
        low_dim_key (str): Key in adata.obsm for embedding (default: 'X_umap').
        high_dim_key (str): Key for reference data. 'X' means adata.X. (default: 'X').
        k_neighbors (int): Number of neighbors for neighborhood-based metrics.
        n_samples (int or None): Number of cells to subsample. None = use all cells.
        seed (int): Random seed for reproducibility.
        n_jobs (int): Number of parallel jobs (-1 = all cores).
        all_metrics (bool): If True, calculate all metrics. If False, use a default subset.
        file_dir (str): Directory to save results CSV. If None, saves to cwd.
        verbose (bool): Show progress bars and status messages.
        use_memmap (bool): If True, use memory-mapped arrays for distance matrices (recommended for large datasets).
        temp_dir (str): Directory for temporary memmap files. If None, uses system temp.
        chunk_size (int): Chunk size for distance computation (default: 1000).
                     
    Returns:
        pd.DataFrame: Results dataframe with columns:
                      [Atlas Name, Metric Name, Low Dim Key, High Dim Key, Value, Time (s)]
    """
    # Import here to avoid circular dependency
    from ..atlas import CheckAtlasColumnDetector
    import scanpy as sc
    from sklearn.metrics import pairwise_distances
    from sklearn.neighbors import NearestNeighbors
    from joblib import Parallel, delayed

    ## if the temp_dir is not provided, use the 'file_dir'
    if temp_dir is None:
        temp_dir = file_dir
    
    # Set file directory
    if file_dir is None:
        file_dir = os.getcwd()
    else:
        os.makedirs(file_dir, exist_ok=True)
    
    # Check if low dim key exists, compute if missing
    if low_dim_key not in adata.obsm.keys():
        if verbose: 
            print(f"Key '{low_dim_key}' not found. Calculating UMAP...")
        sc.tl.umap(adata, n_components=2, random_state=seed)

    # Prepare data for high dimension
    if high_dim_key == 'X':
        if verbose: print("Using adata.X as high-dimensional reference.")
        # Check if X is sparse, convert to dense if needed for distance calc
        # Note: pairwise_distances handles sparse, but dense might be faster for smaller subsamples
        pass 
    elif high_dim_key in adata.obsm.keys():
        pass # Use as is
    else:
        # Fallback or error? User requested adata.X mostly.
        if verbose: print(f"Key '{high_dim_key}' not found. Using adata.X.")
        high_dim_key = 'X'

    # --- Precomputation Step ---
    if verbose: print("Precomputing distances and neighbors to optimize metric calculation...")
    
    # 1. Subsampling (if specified)
    n_obs = adata.n_obs
    if n_samples is not None and n_samples < n_obs:
        if verbose: print(f"Subsampling {n_samples} cells out of {n_obs}...")
        np.random.seed(seed)
        sample_indices = np.random.choice(n_obs, n_samples, replace=False)
    else:
        if verbose: print(f"Using all {n_obs} cells...")
        sample_indices = np.arange(n_obs)

    # Get data subsets
    if high_dim_key == 'X':
        high_dim_data = adata.X[sample_indices]
        # Convert sparse to dense if needed
        if hasattr(high_dim_data, "toarray"):
            high_dim_data = high_dim_data.toarray()
    else:
        high_dim_data = adata.obsm[high_dim_key][sample_indices]
        
    low_dim_data = adata.obsm[low_dim_key][sample_indices]

    n_cells = high_dim_data.shape[0]
    
    # Setup temp directory for memmap files
    import tempfile
    import gc
    
    if temp_dir is None:
        temp_dir = tempfile.gettempdir()
    os.makedirs(temp_dir, exist_ok=True)
    
    # Generate unique filenames for this run
    import uuid
    run_id = str(uuid.uuid4())[:8]
    high_dists_path = os.path.join(temp_dir, f"high_dim_dists_{run_id}.dat")
    low_dists_path = os.path.join(temp_dir, f"low_dim_dists_{run_id}.dat")
    
    # Track memmap files for cleanup
    memmap_files = [high_dists_path, low_dists_path]
    
    # Initialize variables to None to prevent UnboundLocalError
    high_dim_dists = None
    low_dim_dists = None
    high_knn_indices = None
    low_knn_indices = None
    high_knn_dists = None
    low_knn_dists = None

    try:
        # 2. Compute Pairwise Distances using Memory-Mapped Arrays
        if use_memmap:
            if verbose: 
                print(f"Computing pairwise distances using memory-mapped arrays...")
                print(f"  Memory-efficient mode: storing {n_cells}x{n_cells} matrices on disk")
            
            # Create memory-mapped arrays on disk
            high_dim_dists = np.memmap(high_dists_path, dtype='float32', mode='w+', shape=(n_cells, n_cells))
            low_dim_dists = np.memmap(low_dists_path, dtype='float32', mode='w+', shape=(n_cells, n_cells))
            
            # Fill in chunks to minimize RAM usage
            for i in tqdm(range(0, n_cells, chunk_size), desc="High-Dim Distances", disable=not verbose):
                end = min(i + chunk_size, n_cells)
                high_dim_dists[i:end, :] = pairwise_distances(
                    high_dim_data[i:end], high_dim_data, metric='euclidean', n_jobs=n_jobs
                )
                high_dim_dists.flush()  # Write to disk immediately
            
            for i in tqdm(range(0, n_cells, chunk_size), desc="Low-Dim Distances", disable=not verbose):
                end = min(i + chunk_size, n_cells)
                low_dim_dists[i:end, :] = pairwise_distances(
                    low_dim_data[i:end], low_dim_data, metric='euclidean', n_jobs=n_jobs
                )
                low_dim_dists.flush()  # Write to disk immediately
            
            # Force garbage collection to free any temporary memory
            gc.collect()
            
        else:
            # Standard in-memory computation (for smaller datasets)
            if verbose: print("Computing pairwise distances in memory...")
            
            high_dim_dists = np.zeros((n_cells, n_cells), dtype=np.float32)
            for i in tqdm(range(0, n_cells, chunk_size), desc="High-Dim Distances", disable=not verbose):
                end = min(i + chunk_size, n_cells)
                high_dim_dists[i:end, :] = pairwise_distances(
                    high_dim_data[i:end], high_dim_data, metric='euclidean', n_jobs=n_jobs
                )
            
            low_dim_dists = np.zeros((n_cells, n_cells), dtype=np.float32)
            for i in tqdm(range(0, n_cells, chunk_size), desc="Low-Dim Distances", disable=not verbose):
                end = min(i + chunk_size, n_cells)
                low_dim_dists[i:end, :] = pairwise_distances(
                    low_dim_data[i:end], low_dim_data, metric='euclidean', n_jobs=n_jobs
                )
        
        # 3. Compute k-NN Graphs (always in memory - small footprint)
        if verbose: print(f"Computing k-NN graphs (k={k_neighbors}) for High-Dim and Low-Dim...")
        
        # Query k+1 to exclude self
        nbrs_high = NearestNeighbors(n_neighbors=k_neighbors+1, n_jobs=n_jobs).fit(high_dim_data)
        high_knn_dists_list, high_knn_indices_list = [], []
        for i in tqdm(range(0, n_cells, chunk_size), desc="High-Dim KNN", disable=not verbose):
            end = min(i + chunk_size, n_cells)
            d, ind = nbrs_high.kneighbors(high_dim_data[i:end])
            high_knn_dists_list.append(d)
            high_knn_indices_list.append(ind)
        high_knn_dists = np.vstack(high_knn_dists_list)
        high_knn_indices = np.vstack(high_knn_indices_list)
        
        nbrs_low = NearestNeighbors(n_neighbors=k_neighbors+1, n_jobs=n_jobs).fit(low_dim_data)
        low_knn_dists_list, low_knn_indices_list = [], []
        for i in tqdm(range(0, n_cells, chunk_size), desc="Low-Dim KNN", disable=not verbose):
            end = min(i + chunk_size, n_cells)
            d, ind = nbrs_low.kneighbors(low_dim_data[i:end])
            low_knn_dists_list.append(d)
            low_knn_indices_list.append(ind)
        low_knn_dists = np.vstack(low_knn_dists_list)
        low_knn_indices = np.vstack(low_knn_indices_list)
        
        # Free intermediate lists
        del high_knn_dists_list, high_knn_indices_list, low_knn_dists_list, low_knn_indices_list
        gc.collect()
        
        # Cache kNN data in adata.uns (distance matrices are in memmap)
        adata.uns['_dimred_cache'] = {
            'high_dim_dists': high_dim_dists,  # memmap reference or numpy array
            'low_dim_dists': low_dim_dists,    # memmap reference or numpy array
            'high_knn_indices': high_knn_indices,
            'low_knn_indices': low_knn_indices,
            'high_knn_dists': high_knn_dists,
            'low_knn_dists': low_knn_dists,
            'sample_indices': sample_indices,
            'use_memmap': use_memmap,
            'memmap_files': memmap_files if use_memmap else []
        }

    except Exception as e:
        logger.warning(f"Precomputation failed: {e}. Metrics will attempt local calculation.")
        if verbose: print(f"Warning: Precomputation failed ({e}). Proceeding with local calculation (slower).")
        pass

    # Define metrics to run
    if all_metrics:
        metrics_list = METRICS_DIMRED
    else:
        # Default subset of metrics
        metrics_list = [
            'trustworthiness',
            'continuity',
            'kruskal_stress',
            # 'spearman_rho',
            # 'dCor',
            # 'lcmc',
            # 'entourage',
            # 'coknn',
            # 'ged',
            'avg_jaccard_dis',
            # 'den_pre'
        ]
        # Filter to ensure they exist in METRICS_DIMRED
        metrics_list = [m for m in metrics_list if m in METRICS_DIMRED]

    results = []
    
    # Create progress bar
    pbar = tqdm(metrics_list, desc="Calculating Dimred Metrics", 
                bar_format='{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]',
                disable=not verbose)
    
    for metric in pbar:
        # Update progress bar with current metric name
        pbar.set_description(f"Processing: {metric}")
        
        # Start timing
        metric_start_time = time.time()
        
        try:
            metric_module = getattr(dimred, metric)
            
            # Call the metric's run function with appropriate parameters
            # All dimred metrics now have a consistent signature:
            # run(adata, low_dim_key, high_dim_key, k_neighbors, n_samples, seed, n_jobs, verbose, **kwargs)
            
            import inspect
            sig = inspect.signature(metric_module.run)
            params = sig.parameters
            
            # Build kwargs based on what the metric accepts
            kwargs = {
                'low_dim_key': low_dim_key,
                'high_dim_key': high_dim_key,
                'seed': seed,
            }
            
            if 'k_neighbors' in params:
                kwargs['k_neighbors'] = k_neighbors
            if 'n_samples' in params:
                kwargs['n_samples'] = n_samples
            if 'n_jobs' in params:
                kwargs['n_jobs'] = n_jobs
            if 'verbose' in params:
                kwargs['verbose'] = verbose 
            
            # Pass precomputed data if accepted
            if 'precomputed_high_dists' in params:
                kwargs['precomputed_high_dists'] = high_dim_dists
            if 'precomputed_low_dists' in params:
                kwargs['precomputed_low_dists'] = low_dim_dists
            if 'precomputed_high_knn' in params:
                kwargs['precomputed_high_knn'] = high_knn_indices
            if 'precomputed_low_knn' in params:
                kwargs['precomputed_low_knn'] = low_knn_indices
            if 'precomputed_high_knn_dists' in params:
                kwargs['precomputed_high_knn_dists'] = high_knn_dists
            if 'precomputed_low_knn_dists' in params:
                kwargs['precomputed_low_knn_dists'] = low_knn_dists
            
            # Run the metric
            value = metric_module.run(adata, **kwargs)
            
            # Calculate elapsed time
            elapsed_time = time.time() - metric_start_time
            
            results.append({
                'Atlas Name': atlas_name,
                'Metric Name': metric,
                'Low Dim Key': low_dim_key,
                'High Dim Key': high_dim_key,
                'Value': value,
                'Time (s)': round(elapsed_time, 3)
            })
            
            # Update progress bar with time
            pbar.set_postfix_str(f"Time: {elapsed_time:.2f}s", refresh=True)
            
        except Exception as e:
            elapsed_time = time.time() - metric_start_time
            logger.warning(f"Failed to calculate {metric}: {e}")
            results.append({
                'Atlas Name': atlas_name,
                'Metric Name': metric,
                'Low Dim Key': low_dim_key,
                'High Dim Key': high_dim_key,
                'Value': np.nan,
                'Time (s)': round(elapsed_time, 3)
            })
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Calculate total time
    total_time = df['Time (s)'].sum()
    if verbose:
        print(f"\nTotal computation time: {total_time:.2f}s")
    
    # Save to CSV if results exist
    if not df.empty:
        filename = os.path.join(file_dir, f"checkatlas_dimred_metrics_{atlas_name}.csv")
        df.to_csv(filename, index=False)
        if verbose:
            print(f"Saved dimred metrics to {filename}")
        logger.info(f"Saved dimred metrics to {filename}")
        
    # Cleanup memory-mapped files before returning
    if use_memmap:
        if verbose: print("Cleaning up temporary memmap files...")
        # Close memmap references
        try:
            if 'high_dim_dists' in dir() and hasattr(high_dim_dists, '_mmap'):
                del high_dim_dists
            if 'low_dim_dists' in dir() and hasattr(low_dim_dists, '_mmap'):
                del low_dim_dists
        except:
            pass
        gc.collect()
        
        # Remove files from disk
        for fpath in memmap_files:
            try:
                if os.path.exists(fpath):
                    os.remove(fpath)
            except Exception as e:
                logger.warning(f"Failed to cleanup temp file {fpath}: {e}")
                
        # Clear cache reference
        if '_dimred_cache' in adata.uns:
            adata.uns['_dimred_cache']['memmap_files'] = []
            
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
