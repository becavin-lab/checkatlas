import numpy as np
import scanpy as sc
from sklearn.metrics import pairwise_distances
from sklearn.neighbors import NearestNeighbors
from joblib import Parallel, delayed
from tqdm import tqdm
import gc

def run(adata, low_dim_key='X_umap', high_dim_key='X', k_neighbors=30, 
        n_samples=None, seed=42, verbose=True, n_jobs=-1,
        precomputed_high_dists=None, precomputed_low_dists=None,
        use_memmap=True):
    """
    Trustworthiness (Memory-Optimized)
    Measures the preservation of local neighborhood by penalizing false neighbors (points that are neighbors in low-dim but not in high-dim).
    
    Optimized for memory usage by processing rows in parallel without storing full N*N rank matrices.

    Args:
        adata (AnnData): Annotated data matrix.
        low_dim_key (str): Key for low-dimensional embedding in adata.obsm.
        high_dim_key (str): Key for high-dimensional data (default: 'X').
        k_neighbors (int): Number of neighbors to consider.
        n_samples (int): Number of samples to subsample for calculation. None = all.
        seed (int): Random seed for reproducibility.
        verbose (bool): Whether to print progress.
        n_jobs (int): Number of parallel jobs for computation.
        precomputed_high_dists (np.ndarray or memmap): Precomputed high-dim distance matrix.
        precomputed_low_dists (np.ndarray or memmap): Precomputed low-dim distance matrix.
        use_memmap (bool): Hint that we are using memory mapped files (not strictly used logic-wise but good for API).

    Returns:
        float: The Trustworthiness score.

    Interpretation:
        Range 0 to 1.
        Higher is better (1 means perfect trustworthiness).
    """
    
    # Determine number of workers
    if n_jobs == -1:
        import os
        n_workers = os.cpu_count() or 4
    else:
        n_workers = max(1, n_jobs)

    # 1. Use Precomputed Distances if available
    if precomputed_high_dists is not None and precomputed_low_dists is not None:
        if verbose: print("Using precomputed distance matrices...")
        high_dists = precomputed_high_dists
        low_dists = precomputed_low_dists
        n_cells = high_dists.shape[0]
        
    else:
        # Fallback to local computation
        if verbose: print("Precomputed distances not provided. Calculating locally...")
        
        # Check keys
        if low_dim_key not in adata.obsm.keys():
            if verbose: print(f"Calculating {low_dim_key}...")
            sc.tl.umap(adata, n_components=2, random_state=seed)

        # Prepare Data
        n_obs = adata.n_obs
        if n_samples is not None and n_samples < n_obs:
            if verbose: print(f"Subsampling {n_samples} cells...")
            np.random.seed(seed)
            indices = np.random.choice(n_obs, n_samples, replace=False)
        else:
            indices = np.arange(n_obs)

        # High Dim Data
        if high_dim_key == 'X':
            high_dim_data = adata.X[indices]
            if hasattr(high_dim_data, "toarray"): high_dim_data = high_dim_data.toarray()
        else:
            if high_dim_key not in adata.obsm.keys():
                if verbose: print(f"Calculating {high_dim_key}...")
                sc.tl.pca(adata, random_state=seed)
            high_dim_data = adata.obsm[high_dim_key][indices]

        low_dim_data = adata.obsm[low_dim_key][indices]
        n_cells = high_dim_data.shape[0]

        # Compute Distances (Chunked)
        # If not provided, we compute them on the fly? 
        # For Trustworthiness, we need random access to rows. 
        # If we don't have memmap, we must compute full matrix in RAM or recompute rows.
        # Computing full matrix in RAM defeats the purpose.
        # Recomputing rows on demand is slow but memory safe.
        
        # NOTE: cal_dimred usually provides precomputed memmaps now.
        # If run standalone, we warn.
        if verbose and n_cells > 10000:
            print("Warning: Running on large dataset without precomputed memmap distances. usage may be high.")
            
        high_dists = pairwise_distances(high_dim_data, n_jobs=n_workers)
        low_dists = pairwise_distances(low_dim_data, n_jobs=n_workers)

    # 2. Calculate Trustworthiness Row-by-Row
    # Metric: T(k) = 1 - (2 / (N*k*(2N-3k-1))) * sum_{i} sum_{j in U_k(i)} max(0, r(i,j) - k)
    
    if verbose: print(f"Computing Trustworthiness (k={k_neighbors}) row-by-row...")

    # We process in batches to allow joblib parallelism
    batch_size = 100 # Process 100 rows per job
    n_batches = (n_cells + batch_size - 1) // batch_size
    batches = [range(i, min(i + batch_size, n_cells)) for i in range(0, n_cells, batch_size)]

    def _process_batch(row_indices):
        batch_penalty = 0.0
        for i in row_indices:
            # 1. Get neighbors in LOW dimension
            # We need the k neighbors of i in low_dists
            # Sorting the whole row takes O(N log N)
            # argpartition takes O(N)
            l_row = low_dists[i]
            # indices of k+1 smallest distances (including self at 0)
            knn_indices = np.argpartition(l_row, k_neighbors+1)[:k_neighbors+1]
            
            # Sort these specifically to find exact order? No, just need set of neighbors.
            # But we need to exclude self (d=0). argpartition doesn't guarantee order.
            # Let's verify self is in there. Usually yes.
            
            # Refine: get exact k neighbors excluding self
            # We can just sort the small set of k+1 candidates
            candidate_dists = l_row[knn_indices]
            sorted_args = np.argsort(candidate_dists)
            knn_indices_sorted = knn_indices[sorted_args]
            
            # Exclude self (closest, index i)
            # If i is in list (dist ~ 0), it should be first.
            if knn_indices_sorted[0] == i:
                neighbors = knn_indices_sorted[1:k_neighbors+1]
            else:
                # Fallback if i not found (numerical noise?)
                # Just take first k that are not i
                neighbors = knn_indices[knn_indices != i][:k_neighbors]

            # 2. Calculate rank of these neighbors in HIGH dimension
            h_row = high_dists[i]
            
            # We need rank of each neighbor j in h_row
            # Rank = count of elements smaller than h_row[j]
            # Since we need ranks for strictly k items, we can iterate.
            # But calculating rank for one item is O(N). Doing it k times is O(k*N).
            # Sorting h_row is O(N log N).
            # If K is small (30) and N is large (50k), K*N ~ 1.5M ops. N log N ~ 800k ops.
            # Sorting is faster.
            
            # However, we only need to know if rank > k.
            # Actually we need the value of (rank - k).
            
            # FULL SORT APPROACH (Memory Efficient than storing rank matrix)
            # Just calculating argsort of h_row gives us the order.
            # Ranks are positions in argsort.
            # rank[j] is where j appears in argsort.
            
            # argsort of h_row: indices ordered by distance
            # e.g. [i, closest1, closest2, ...]
            # We can invert this permutation to get rank of each index.
            
            # This is O(N) memory per thread. Safe.
            
            # Optimization: We only care about rank of 'neighbors'.
            # Extract their distances: d_j for j in neighbors.
            # Rank of j = count(h_row < d_j). 
            # This handles ties by using 'min' rank or similar?
            # Sklearn implementation defines rank r(i, j) as "the number of sample points 
            # closer to i than j in the high-dimensional space".
            
            neighbor_dists_high = h_row[neighbors]
            
            for dist_j in neighbor_dists_high:
                # Count points closer than j
                # Subtract 1 because 'i' itself has dist 0? 
                # Definition usually excludes self.
                # If we count all points < dist_j, that includes self (dist=0).
                # So rank = count - 1?
                # Actually, if we just count valid neighbors (excluding self), 
                # then rank 1 is closest neighbor.
                
                # count_nonzero is fast C-level loop
                r = np.count_nonzero(h_row < dist_j)
                
                # h_row contains self at 0. So r includes self. 
                # So rank amongst others is r. (e.g. 1 smaller means only self is smaller -> rank 1).
                # If rank > k, penalty.
                
                if r > k_neighbors:
                    batch_penalty += (r - k_neighbors)
                    
        return batch_penalty

    # Execute parallel
    penalties = Parallel(n_jobs=n_workers)(
        delayed(_process_batch)(batch) for batch in tqdm(batches, desc="Computing penalties", disable=not verbose)
    )
    
    total_penalty = sum(penalties)
    
    # 3. Normalize Score
    if n_cells <= k_neighbors + 1:
        return 1.0

    max_penalty = n_cells * k_neighbors * (2 * n_cells - 3 * k_neighbors - 1) / 2.0
    
    if max_penalty == 0:
        return 1.0
        
    score = 1.0 - (2.0 / max_penalty) * total_penalty
    return max(0.0, min(1.0, score))
