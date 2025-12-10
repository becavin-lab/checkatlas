import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse, csr_matrix
from scipy.sparse.csgraph import dijkstra
from sklearn.neighbors import NearestNeighbors
import multiprocessing as mp
from numba import njit
import warnings

def run(X, label, perplexity=30):
    """
    Calculate Local Inverse Simpson's Index (LISI) for batch mixing evaluation.
    
    This implementation uses a graph-based approach similar to scib (Single-Cell Integration Benchmark).
    It constructs a kNN graph and computes the Inverse Simpson's Index of the batch labels
    in the local neighborhood of each cell, defined by shortest paths on the graph.
    
    The output is scaled to [0, 1] where 1 indicates good mixing (iLISI).
    
    :param X: array-like of shape (n_samples, n_features)
        Feature matrix (e.g., PCA embedding, integrated space)
    :param label: batch column name in ``adata.obs`` or array-like
        Batch labels to evaluate mixing
    :param perplexity: int, default=30
        Perplexity parameter for local neighborhood size
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
    
    # Create AnnData object for scanpy compatibility
    adata = sc.AnnData(X=X)
    adata.obs['batch'] = label
    adata.obs['batch'] = adata.obs['batch'].astype('category')
    
    # Compute kNN graph if not present (using Euclidean distance on X)
    # scib uses k=15 for the graph construction usually
    n_neighbors_graph = 15
    sc.pp.neighbors(adata, n_neighbors=n_neighbors_graph, use_rep='X')
    
    # Compute LISI
    # We use the graph-based LISI implementation
    # For iLISI (batch mixing), we want high diversity.
    
    lisi_scores = lisi_graph_py(
        adata=adata,
        obs_key='batch',
        n_neighbors=90, # scib default for LISI calculation
        perplexity=perplexity,
        verbose=False
    )
    
    # Scale iLISI
    n_batches = adata.obs['batch'].nunique()
    if n_batches > 1:
        mean_lisi = np.nanmean(lisi_scores)
        scaled_lisi = (mean_lisi - 1) / (n_batches - 1)
    else:
        scaled_lisi = 0.0

    return scaled_lisi

def ilisi_graph(adata, batch_key, k0=90, perplexity=None, scale=True, verbose=False):
    """
    Integration LISI (iLISI) score.
    Wrapper function to calculate iLISI on an AnnData object.
    """
    scores = lisi_graph_py(adata, batch_key, n_neighbors=k0, perplexity=perplexity, verbose=verbose)
    ilisi = np.nanmedian(scores)
    
    if scale:
        nbatches = adata.obs[batch_key].nunique()
        if nbatches > 1:
            ilisi = (ilisi - 1) / (nbatches - 1)
        else:
            ilisi = 0.0
            
    return ilisi

def clisi_graph(adata, label_key, k0=90, perplexity=None, scale=True, verbose=False):
    """
    Cell-type LISI (cLISI) score.
    Wrapper function to calculate cLISI on an AnnData object.
    """
    scores = lisi_graph_py(adata, label_key, n_neighbors=k0, perplexity=perplexity, verbose=verbose)
    clisi = np.nanmedian(scores)
    
    if scale:
        nlabs = adata.obs[label_key].nunique()
        if nlabs > 1:
            clisi = (nlabs - clisi) / (nlabs - 1)
        else:
            clisi = 0.0
            
    return clisi

def lisi_graph_py(
    adata,
    obs_key,
    n_neighbors=90,
    perplexity=None,
    subsample=None,
    n_cores=1,
    verbose=False,
):
    """
    Compute LISI score on shortest path based on kNN graph provided in the adata object.
    Re-implementation of scib's logic using scipy.sparse.csgraph.dijkstra.
    """
    if "neighbors" not in adata.uns:
        raise AttributeError(
            "Key 'neighbors' not found. Please make sure that a kNN graph has been computed"
        )
    
    # Get distances (CSR matrix)
    # scib uses connectivities, but for shortest path we usually want distances.
    # If 'distances' is available, use it.
    if "distances" in adata.obsp:
        graph = adata.obsp["distances"]
    else:
        # Fallback to connectivities (similarity), inverted? 
        # Or just treat as distance (if unweighted).
        # For now, let's assume distances are present if neighbors was run.
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
    
    # Compute shortest paths using Dijkstra
    # We need to find k nearest neighbors for each cell in the subset
    # limit=np.inf is default, but we can set a limit if we knew the radius.
    # Since we don't, we might have to compute more.
    # However, dijkstra with return_predecessors=False returns a matrix.
    # For large N, this is N*N.
    # We can process in chunks to avoid OOM.
    
    chunk_size = 1000 # Adjust based on memory
    lisi_scores = np.zeros(n_cells)
    lisi_scores[:] = np.nan # Initialize with NaN
    
    # Process subset
    for i in range(0, len(subset_indices), chunk_size):
        chunk = subset_indices[i : i + chunk_size]
        
        # Compute distances from chunk nodes to all other nodes
        # This is still expensive (chunk_size * N).
        # If N=100k, 1000 * 100k * 8 bytes = 800MB. Manageable.
        
        try:
            # limit parameter was added in scipy 1.10.0
            # We can use it to speed up if we had a radius, but we need exactly k neighbors.
            # So we compute full row and sort.
            dists = dijkstra(graph, indices=chunk, return_predecessors=False)
        except TypeError:
             # Fallback for older scipy if needed, but we checked 1.10.1
             dists = dijkstra(graph, indices=chunk, return_predecessors=False)

        # For each cell in chunk, find k nearest
        for j, cell_idx in enumerate(chunk):
            row_dists = dists[j]
            # Get indices of k nearest neighbors
            # argpartition is faster than sort
            knn_indices = np.argpartition(row_dists, n_neighbors)[:n_neighbors]
            # We also need the actual distances for these neighbors for LISI calculation
            knn_dists = row_dists[knn_indices]
            
            # Compute Simpson Index for this cell
            simpson = compute_simpson_index_cell(
                knn_dists, knn_indices, batch_labels, n_batches, perplexity
            )
            lisi_scores[cell_idx] = 1.0 / simpson

    return lisi_scores

@njit
def compute_simpson_index_cell(distances, indices, batch_labels, n_batches, perplexity, tol=1e-5):
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
    # We can use a simple loop or array accumulation
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
