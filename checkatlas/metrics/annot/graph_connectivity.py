import numpy as np
from scipy.sparse import issparse
from scipy.sparse.csgraph import connected_components


def run(adata, neighbors_key='neighbors'):
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
    :return: float
        Graph connectivity score. Range: [0, 1], where 1 indicates perfect connectivity.
        
    """
    # Check if neighbors graph exists
    if neighbors_key not in adata.uns:
        raise ValueError(f"Neighbors key '{neighbors_key}' not found in adata.uns. "
                        f"Run sc.pp.neighbors() first.")
    
    # Get connectivity matrix
    if 'connectivities' in adata.obsp:
        connectivities = adata.obsp['connectivities']
    elif 'connectivities' in adata.uns[neighbors_key]:
        connectivities = adata.uns[neighbors_key]['connectivities']
    else:
        raise ValueError("Connectivity matrix not found in adata")
    
    # Calculate connected components
    n_components, labels = connected_components(
        csgraph=connectivities,
        directed=False,
        return_labels=True
    )
    
    # Calculate connectivity score
    # Perfect connectivity = 1 (all cells in one component)
    # Poor connectivity = close to 0 (many small components)
    n_cells = adata.n_obs
    largest_component_size = np.max(np.bincount(labels))
    
    # Score based on fraction of cells in largest component
    connectivity_score = largest_component_size / n_cells
    
    return connectivity_score
