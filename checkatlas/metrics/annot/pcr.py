import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from scipy.sparse import issparse


def run(adata, batch_label = "batch", n_components=50, cv=5):
    """
    Calculate Principal Component Regression (PCR) batch effect score.
    
    PCR measures how much variance in the principal components can be explained
    by batch labels. Lower scores indicate better batch correction.
    
    `PCR readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/integration/pcr/>`__

    :param adata: AnnData object
        AnnData object containing the feature matrix (raw or integrated data)
        Can be sparse or dense matrix
    :param label: batch column name in ``adata.obs`` or array-like
        Batch labels to evaluate mixing
    :param n_components: int, default=50
        Number of principal components to use
    :param cv: int, default=5
        Number of cross-validation folds
    :return: float
        PCR score (mean cross-validated accuracy).
        Range: [0, 1], where lower values indicate better batch correction.
    """
    # Handle sparse matrices
    if issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = np.asarray(adata.X)
    
    ## 'batch' col check in adata.obs
    if batch_label not in adata.obs:
        print(f"Batch label '{batch_label}' not found in adata.obs")
        return None
    
    batch_labels = np.asarray(adata.obs[batch_label])
    
    n_samples = X.shape[0]
    
    # Check if we have multiple batches
    unique_batches = np.unique(batch_labels)
    n_batches = len(unique_batches)
    
    if n_batches == 1:
        # Only one batch, perfect batch correction
        return 1.0 / n_batches  # Equivalent to random guessing
    
    # Adjust n_components if needed
    max_components = min(n_samples, X.shape[1])
    n_components = min(n_components, max_components)
    
    # Perform PCA
    pca = PCA(n_components=n_components, random_state=0)
    X_pca = pca.fit_transform(X)
    
    # Encode batch labels as integers
    label_encoder = {label: i for i, label in enumerate(unique_batches)}
    y = np.array([label_encoder[label] for label in batch_labels])
    
    # Train logistic regression classifier with cross-validation
    clf = LogisticRegression(max_iter=1000, random_state=0, multi_class='auto')
    
    # Cross-validated accuracy
    try:
        scores = cross_val_score(clf, X_pca, y, cv=cv, scoring='accuracy')
        pcr_score = np.mean(scores)
    except Exception as e:
        # Fallback: simple train/test split if CV fails
        from sklearn.model_selection import train_test_split
        X_train, X_test, y_train, y_test = train_test_split(
            X_pca, y, test_size=0.2, random_state=0
        )
        clf.fit(X_train, y_train)
        pcr_score = clf.score(X_test, y_test)
    
    return pcr_score
