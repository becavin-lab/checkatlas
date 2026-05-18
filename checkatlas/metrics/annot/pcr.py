import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from scipy.sparse import issparse
from anndata import AnnData


def run(X, labels, n_components=50, cv=5, n_jobs=-1, verbose=True):
    """
    Calculate Principal Component Regression (PCR) batch effect score.

    PCR measures how much variance in the principal components can be
    explained by batch labels.  Lower scores indicate better batch
    correction.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Feature matrix (e.g. PCA embedding, integrated space).
    labels : array-like of shape (n_samples,)
        Batch labels.
    n_components : int, default=50
        Number of principal components to use.
    cv : int, default=5
        Number of cross-validation folds.
    n_jobs : int, default=-1
        Number of parallel jobs.  -1 uses all cores.
    verbose : bool, default=True
        Print progress information.

    Returns
    -------
    float
        PCR score (mean cross-validated accuracy).  Range [0, 1];
        lower values indicate better batch correction.
    """
    # ── backward compat: if first arg is AnnData ──────────────────
    if isinstance(X, AnnData):
        adata = X
        batch_label = labels if isinstance(labels, str) else "batch"
        if issparse(adata.X):
            X_arr = adata.X.toarray()
        else:
            X_arr = np.asarray(adata.X)
        labels_arr = np.asarray(adata.obs[batch_label])
        return run(X_arr, labels_arr, n_components=n_components, cv=cv,
                   n_jobs=n_jobs, verbose=verbose)

    if issparse(X):
        X = X.toarray()
    else:
        X = np.asarray(X, dtype=np.float64)

    labels = np.asarray(labels)
    n_samples = X.shape[0]

    unique_batches = np.unique(labels)
    n_batches = len(unique_batches)
    if n_batches <= 1:
        return 1.0 / max(n_batches, 1)

    if verbose:
        print(f"Computing PCR ({n_samples:,} samples, "
              f"{n_batches} batches)...")

    # Adjust n_components
    max_components = min(n_samples, X.shape[1])
    n_components = min(n_components, max_components)

    # Run PCA (skip if data is already low-dim)
    if X.shape[1] > n_components:
        if verbose:
            print(f"  Running PCA (n_components={n_components})...")
        pca = PCA(n_components=n_components, random_state=0)
        X_pca = pca.fit_transform(X)
    else:
        X_pca = X.copy()

    # Encode labels
    label_encoder = {lbl: i for i, lbl in enumerate(unique_batches)}
    y = np.array([label_encoder[lbl] for lbl in labels])

    # LogisticRegression CV
    if verbose:
        print(f"  Running {cv}-fold cross-validation (n_jobs={n_jobs})...")

    clf = LogisticRegression(max_iter=1000, random_state=0,
                             multi_class='auto', n_jobs=n_jobs)
    try:
        scores = cross_val_score(clf, X_pca, y, cv=cv, scoring='accuracy',
                                 n_jobs=n_jobs)
        pcr_score = float(np.mean(scores))
    except Exception:
        if verbose:
            print(f"  CV failed, falling back to train/test split...")
        from sklearn.model_selection import train_test_split
        X_train, X_test, y_train, y_test = train_test_split(
            X_pca, y, test_size=0.2, random_state=0
        )
        clf.fit(X_train, y_train)
        pcr_score = float(clf.score(X_test, y_test))

    if verbose:
        print(f"  PCR score = {pcr_score:.4f}")

    return pcr_score
