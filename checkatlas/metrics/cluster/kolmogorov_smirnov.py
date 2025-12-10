import numpy as np
from scipy.stats import ks_2samp
from scipy.sparse import issparse
import warnings

def run(X, labels):
    """
    Calculate the Kolmogorov-Smirnov (KS) Statistic for clustering separation.
    
    This metric evaluates the degree of separation between clusters by comparing
    the distributions of features across different clusters. For each pair of clusters,
    it computes the KS statistic for each feature and averages them. A higher KS score
    indicates that the clusters have distinct feature distributions.
    
    Range: [0, 1], where 1 indicates distinct distributions (good separation).
    
    :param X: array-like of shape (n_samples, n_features)
        Feature matrix.
    :param labels: array-like of shape (n_samples,)
        Cluster labels.
    :return: float
        Average KS statistic over all cluster pairs.
    """
    if issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X)
        
    labels = np.asarray(labels)
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)
    
    if n_clusters < 2:
        return 0.0
        
    ks_scores = []
    
    # Iterate over all pairs of clusters
    for i in range(n_clusters):
        for j in range(i + 1, n_clusters):
            label_i = unique_labels[i]
            label_j = unique_labels[j]
            
            # Get data for each cluster
            X_i = X[labels == label_i]
            X_j = X[labels == label_j]
            
            # Compute KS statistic for each feature
            # We can average over features
            n_features = X.shape[1]
            feature_ks_scores = []
            
            # Optimization: If many features, maybe select top variable ones?
            # Or just average all. For scRNA-seq PCA (50 dims), it's fast.
            # For full gene expression (20k dims), it might be slow.
            # Assuming X is usually PCA/embedding (low dim).
            
            for k in range(n_features):
                # ks_2samp returns (statistic, pvalue)
                stat, _ = ks_2samp(X_i[:, k], X_j[:, k])
                feature_ks_scores.append(stat)
            
            # Average KS for this pair of clusters
            avg_ks_pair = np.mean(feature_ks_scores)
            ks_scores.append(avg_ks_pair)
            
    # Average over all pairs
    if not ks_scores:
        return 0.0
        
    final_score = np.mean(ks_scores)
    
    return final_score
