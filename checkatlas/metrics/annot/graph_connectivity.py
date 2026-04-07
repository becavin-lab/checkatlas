import numpy as np
from scipy.sparse import issparse
from scipy.sparse.csgraph import connected_components
import scanpy as sc
import re



def _detect_embedding_keys(adata):
    """
    Detect embedding keys in adata.obsm using regex patterns.
    Adapted from checkatlas.atlas.CheckAtlasColumnDetector.analyze_obsm_semantics
    """
    obsm_keys = list(adata.obsm.keys()) if hasattr(adata, 'obsm') else []
    embedding_patterns = [
        r'^(X_pca|X_umap|X_tsne|X_scvi|X_emb)$',
        r'^X_diffmap$',
        r'^X_draw_graph',
        r'^X_.*pca',
        r'^X_.*umap',
        r'^X_.*tsne',
        r'^X_',  # Match any X_ prefixed key
    ]
    
    matched_keys = []
    
    for key in obsm_keys:
        for pattern in embedding_patterns:
            if re.search(pattern, key):
                matched_keys.append(key)
                break
    
    return matched_keys


def run(adata, neighbors_key='neighbors', label_key=None, n_jobs=-1, verbose=True):
    """
    Calculate Graph Connectivity metric for batch integration quality.
    
    This metric evaluates whether cells from the same biological group remain
    connected in the nearest neighbor graph after batch correction. It measures
    the number of connected components per biological group.
    
    `Graph Connectivity readthedocs
    <https://checkatlas.readthedocs.io/en/latest/metrics/integration/graph_connectivity/>`__

    :param adata: AnnData object
        Annotated data matrix with computed neighbor graph
    :param neighbors_key: str, default='neighbors'
        Key in adata.uns containing the neighbor graph information
    :param label_key: str, optional
        Key in adata.obs containing biological labels (e.g., cell types).
        If provided, connectivity is calculated per group and averaged.
        If None, calculates connectivity for the entire graph.
    :return: float
        Graph connectivity score. Range: [0, 1], where 1 indicates perfect connectivity.
        
    """
    # Check if neighbors graph exists
    if neighbors_key not in adata.uns:
        print(f"Neighbors key '{neighbors_key}' not found in adata.uns. Calculating neighbors...")
        # Detect embedding key
        embedding_keys = _detect_embedding_keys(adata)
        if len(embedding_keys) > 0:
            use_rep = embedding_keys[0]
            print(f"Using embedding key: {use_rep}")
            # Calculate neighbors using the detected embedding
            sc.pp.neighbors(adata, use_rep=use_rep, key_added=neighbors_key if neighbors_key != 'neighbors' else None)
        else:
            print("No embedding key found. Using default behavior.")
            # Fallback to default behavior (likely uses X or X_pca if available)
            sc.pp.neighbors(adata, key_added=neighbors_key if neighbors_key != 'neighbors' else None)
            
        # Re-check if neighbors key exists after calculation
        if neighbors_key not in adata.uns:
            raise ValueError(f"Neighbors key '{neighbors_key}' not found in adata.uns. "
                            f"Run sc.pp.neighbors() first.")
    
    # Get connectivity matrix
    connectivities = None
    if neighbors_key in adata.uns and 'connectivities' in adata.uns[neighbors_key]:
        connectivities = adata.uns[neighbors_key]['connectivities']
    elif f'{neighbors_key}_connectivities' in adata.obsp:
        connectivities = adata.obsp[f'{neighbors_key}_connectivities']
    elif neighbors_key == 'neighbors' and 'connectivities' in adata.obsp:
        connectivities = adata.obsp['connectivities']
    
    if connectivities is None:
        raise ValueError(f"Connectivity matrix not found for key '{neighbors_key}' in adata.uns or adata.obsp")
    
    if label_key:
        if label_key not in adata.obs:
             raise ValueError(f"Label key '{label_key}' not found in adata.obs")
        
        labels = adata.obs[label_key]
        # Ensure categorical
        if not hasattr(labels, 'cat'):
             labels = labels.astype('category')
             
        scores = []
        for cat in labels.cat.categories:
            mask = (labels == cat).values
            if np.sum(mask) == 0:
                continue
                
            # Subset matrix
            # connectivities is usually CSR or CSC. Slicing [mask, :][:, mask] works.
            sub_conn = connectivities[mask, :][:, mask]
            
            n_components, component_labels = connected_components(
                csgraph=sub_conn,
                directed=False,
                connection="strong"
            )
            
            largest_c = np.max(np.bincount(component_labels))
            scores.append(largest_c / np.sum(mask))
            
        return np.mean(scores)
    else:
        # Calculate connected components (Global)
        n_components, labels = connected_components(
            csgraph=connectivities,
            directed=False,
            return_labels=True,
            connection="strong"
        )
        
        # Calculate connectivity score
        # Perfect connectivity = 1 (all cells in one component)
        # Poor connectivity = close to 0 (many small components)
        n_cells = adata.n_obs
        largest_component_size = np.max(np.bincount(labels))
        
        # Score based on fraction of cells in largest component
        connectivity_score = largest_component_size / n_cells
        
        return connectivity_score
