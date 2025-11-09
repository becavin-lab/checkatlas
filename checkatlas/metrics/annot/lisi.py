import numpy as np
from sklearn.neighbors import NearestNeighbors


def run(X, labels, perplexity=30):
    """
    Calculate Local Inverse Simpson's Index (LISI) for batch mixing evaluation.
    
    LISI quantifies the diversity of batches/labels within local neighborhoods.
    It measures how well different batches are mixed in the embedding space.
    
    `LISI readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/integration/lisi/>`__

    :param X: array-like of shape (n_samples, n_features)
        Feature matrix (e.g., PCA embedding, integrated space)
    :param labels: array-like of shape (n_samples,)
        Batch labels or categorical labels to evaluate mixing
    :param perplexity: int, default=30
        Perplexity parameter for local neighborhood size
    :return: float
        Mean LISI score across all cells.
        Range: [1, n_batches], where higher values indicate better mixing.
        
    """
    X = np.asarray(X)
    labels = np.asarray(labels)
    
    n_samples = X.shape[0]
    
    # Get unique labels
    unique_labels = np.unique(labels)
    n_labels = len(unique_labels)
    
    if n_labels == 1:
        # Only one batch, LISI is 1
        return 1.0
    
    # Determine number of neighbors based on perplexity
    # Typical relationship: n_neighbors ≈ 3 * perplexity
    n_neighbors = min(int(3 * perplexity), n_samples - 1)
    
    # Find k-nearest neighbors
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, algorithm='auto')
    nbrs.fit(X)
    distances, indices = nbrs.kneighbors(X)
    
    # Calculate LISI for each cell
    lisi_scores = []
    
    for i in range(n_samples):
        neighbor_indices = indices[i]
        neighbor_labels = labels[neighbor_indices]
        
        # Count occurrences of each label in the neighborhood
        label_counts = {}
        for label in unique_labels:
            label_counts[label] = np.sum(neighbor_labels == label)
        
        # Calculate Simpson's Index
        total_neighbors = len(neighbor_labels)
        simpson_index = 0.0
        
        for count in label_counts.values():
            if count > 0:
                proportion = count / total_neighbors
                simpson_index += proportion ** 2
        
        # LISI is the inverse of Simpson's Index
        if simpson_index > 0:
            lisi = 1.0 / simpson_index
        else:
            lisi = 1.0
        
        lisi_scores.append(lisi)
    
    # Return mean LISI across all cells
    mean_lisi = np.mean(lisi_scores)
    
    return mean_lisi
