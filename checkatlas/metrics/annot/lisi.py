import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse, csr_matrix
from scipy.sparse.csgraph import dijkstra
from sklearn.neighbors import NearestNeighbors
from joblib import Parallel, delayed
from numba import njit
import warnings


def run(X, label, perplexity=30, n_jobs=-1, verbose=True):
    """
    Calculate Local Inverse Simpson's Index (LISI) for batch mixing evaluation.
    
    This implementation uses a graph-based approach similar to scib (Single-Cell
    Integration Benchmark). It constructs a kNN graph and computes the Inverse
    Simpson's Index of the batch labels in the local neighborhood of each cell,
    defined by shortest paths on the graph.
    
    The output is scaled to [0, 1] where 1 indicates good mixing (iLISI).
    
    :param X: array-like of shape (n_samples, n_features)
        Feature matrix (e.g., PCA embedding, integrated space)
    :param label: batch column name in ``adata.obs`` or array-like
        Batch labels to evaluate mixing
    :param perplexity: int, default=30
        Perplexity parameter for local neighborhood size
    :param n_jobs: int, default=-1
        Number of parallel jobs. -1 uses all cores.
    :param verbose: bool, default=True
        Whether to print progress information.
    :return: float
        Mean scaled LISI score across all cells.
        Range: [0, 1], where higher values indicate better mixing.
    """
    # Handle sparse matrices
    if issparse(X):
        X = X.toarray()
    
    # Ensure label is a numpy array
    if isinstance(label, (pd.Series, pd.Index)):
        label = label.values
    else:
        label = np.asarray(label)
    
    if verbose:
        print(f"Computing LISI ({len(X):,} samples, n_jobs={n_jobs})...")
    
    # Create AnnData object for scanpy compatibility
    adata = sc.AnnData(X=X)
    adata.obs['batch'] = label
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    
    # Compute kNN graph with parallel neighbor search
    n_neighbors_graph = 15
    if verbose:
        print(f"  Computing kNN graph (k={n_neighbors_graph})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors_graph, use_rep='X')
    
    # Compute LISI with parallelized chunk processing
    lisi_scores = lisi_graph_py(
        adata=adata,
        obs_key='batch',
        n_neighbors=90,
        perplexity=perplexity,
        n_jobs=n_jobs,
        verbose=verbose
    )
    
    # Scale iLISI
    n_batches = adata.obs['batch'].nunique()
    if n_batches > 1:
        mean_lisi = np.nanmean(lisi_scores)
        scaled_lisi = (mean_lisi - 1) / (n_batches - 1)
    else:
        scaled_lisi = 0.0

    if verbose:
        print(f"  Scaled LISI = {scaled_lisi:.6f}")

    return scaled_lisi


def ilisi_graph(adata, batch_key, k0=90, perplexity=None, scale=True,
                n_jobs=-1, verbose=False):
    """
    Integration LISI (iLISI) score.
    Wrapper function to calculate iLISI on an AnnData object.
    """
    scores = lisi_graph_py(adata, batch_key, n_neighbors=k0, perplexity=perplexity,
                           n_jobs=n_jobs, verbose=verbose)
    ilisi = np.nanmedian(scores)
    
    if scale:
        nbatches = adata.obs[batch_key].nunique()
        if nbatches > 1:
            ilisi = (ilisi - 1) / (nbatches - 1)
        else:
            ilisi = 0.0
            
    return ilisi


def clisi_graph(adata, label_key, k0=90, perplexity=None, scale=True,
                n_jobs=-1, verbose=False):
    """
    Cell-type LISI (cLISI) score.
    Wrapper function to calculate cLISI on an AnnData object.
    """
    scores = lisi_graph_py(adata, label_key, n_neighbors=k0, perplexity=perplexity,
                           n_jobs=n_jobs, verbose=verbose)
    clisi = np.nanmedian(scores)
    
    if scale:
        nlabs = adata.obs[label_key].nunique()
        if nlabs > 1:
            clisi = (nlabs - clisi) / (nlabs - 1)
        else:
            clisi = 0.0
            
    return clisi


def _process_chunk(graph, chunk_indices, batch_labels, n_batches,
                   n_neighbors, perplexity):
    """
    Process a chunk of cells: compute Dijkstra shortest paths and
    Simpson index for each cell. Designed to run in parallel.
    """
    n_cells_total = graph.shape[0]
    chunk_scores = np.full(len(chunk_indices), np.nan)
    
    try:
        dists = dijkstra(graph, indices=chunk_indices, return_predecessors=False)
    except TypeError:
        dists = dijkstra(graph, indices=chunk_indices, return_predecessors=False)
    
    for j, cell_idx in enumerate(chunk_indices):
        row_dists = dists[j]
        # argpartition is faster than full sort for finding k smallest
        if n_neighbors < len(row_dists):
            knn_indices = np.argpartition(row_dists, n_neighbors)[:n_neighbors]
        else:
            knn_indices = np.arange(len(row_dists))
        knn_dists = row_dists[knn_indices]
        
        simpson = compute_simpson_index_cell(
            knn_dists, knn_indices, batch_labels, n_batches, perplexity
        )
        chunk_scores[j] = 1.0 / simpson
    
    return chunk_indices, chunk_scores


def lisi_graph_py(
    adata,
    obs_key,
    n_neighbors=90,
    perplexity=None,
    subsample=None,
    n_jobs=-1,
    verbose=False,
):
    """
    Compute LISI score on shortest path based on kNN graph provided in the
    adata object. Parallelized using joblib for chunk processing.
    """
    if "neighbors" not in adata.uns:
        raise AttributeError(
            "Key 'neighbors' not found. Please make sure that a kNN graph has been computed"
        )
    
    # Get distances (CSR matrix)
    if "distances" in adata.obsp:
        graph = adata.obsp["distances"]
    else:
        graph = adata.obsp["connectivities"]
        
    # Ensure labels are codes
    adata.obs[obs_key] = adata.obs[obs_key].astype("category")
    batch_labels = adata.obs[obs_key].cat.codes.values
    n_batches = adata.obs[obs_key].nunique()

    if perplexity is None or perplexity >= n_neighbors:
        perplexity = np.floor(n_neighbors / 3)

    n_cells = adata.shape[0]
    
    # Subsampling
    subset_indices = np.arange(n_cells)
    if subsample is not None:
        if subsample < 100:
            size = int(n_cells * subsample / 100)
            subset_indices = np.random.choice(n_cells, size, replace=False)
    
    # Split into chunks for parallel processing
    chunk_size = 1000
    chunks = [subset_indices[i:i + chunk_size]
              for i in range(0, len(subset_indices), chunk_size)]
    
    if verbose:
        print(f"  Computing LISI ({len(chunks)} chunks, n_jobs={n_jobs})...")
    
    # Resolve n_jobs for joblib
    import os
    effective_jobs = n_jobs if n_jobs != -1 else os.cpu_count()
    
    if effective_jobs == 1 or len(chunks) == 1:
        # Serial fallback
        lisi_scores = np.full(n_cells, np.nan)
        for chunk in chunks:
            idxs, scores = _process_chunk(
                graph, chunk, batch_labels, n_batches, n_neighbors, perplexity
            )
            lisi_scores[idxs] = scores
    else:
        # Parallel execution
        results = Parallel(n_jobs=n_jobs, prefer="threads")(
            delayed(_process_chunk)(
                graph, chunk, batch_labels, n_batches, n_neighbors, perplexity
            )
            for chunk in chunks
        )
        
        lisi_scores = np.full(n_cells, np.nan)
        for idxs, scores in results:
            lisi_scores[idxs] = scores

    return lisi_scores


@njit
def compute_simpson_index_cell(distances, indices, batch_labels, n_batches,
                               perplexity, tol=1e-5):
    """
    Compute Simpson index for a single cell.
    """
    # Get batch labels of neighbors
    batches = batch_labels[indices]
    
    # Calculate probabilities P using Hbeta (perplexity based)
    beta = 1.0
    betamin = -np.inf
    betamax = np.inf
    logU = np.log(perplexity)
    
    H, P = Hbeta(distances, beta)
    Hdiff = H - logU
    tries = 0
    
    while np.abs(Hdiff) > tol and tries < 50:
        if Hdiff > 0:
            betamin = beta
            if betamax == np.inf:
                beta *= 2.0
            else:
                beta = (beta + betamax) / 2.0
        else:
            betamax = beta
            if betamin == -np.inf:
                beta /= 2.0
            else:
                beta = (beta + betamin) / 2.0
        
        H, P = Hbeta(distances, beta)
        Hdiff = H - logU
        tries += 1
        
    if H == 0:
        return 1.0
        
    # Sum P per batch
    sumP = np.zeros(n_batches)
    for k in range(len(batches)):
        b = batches[k]
        sumP[b] += P[k]
        
    # Simpson index = sum(sumP^2)
    simpson = np.dot(sumP, sumP)
    return simpson


@njit
def Hbeta(D, beta):
    """
    Helper function for simpson index computation.
    Computes entropy and probabilities for a given beta.
    """
    P = np.exp(-D * beta)
    sumP = np.sum(P)
    if sumP == 0:
        H = 0.0
        P = np.zeros(len(D))
    else:
        H = np.log(sumP) + beta * np.sum(D * P) / sumP
        P /= sumP
    return H, P
