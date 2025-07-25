from sklearn.metrics import roc_auc_score
import numpy as np
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import NearestNeighbors

def run(high_dim_counts, low_dim_counts):
    # for i in range(4):
    k_neighbors = 4

    # high_dim_counts = high_dim_counts[0:100,]
    # low_dim_conts = low_dim_counts[0:100,]

    X = high_dim_counts
    neigh = NearestNeighbors(n_neighbors=k_neighbors)
    neigh.fit(X)
    A = neigh.kneighbors_graph(X, mode="distance")

    X_low = low_dim_conts
    neigh_low = NearestNeighbors(n_neighbors=k_neighbors)
    print(f"number k : {k_neighbors}")
    neigh.fit(X_low)
    B = neigh.kneighbors_graph(X_low, mode="distance")

    n = X.shape[0]
    k=k_neighbors

    total_inter = 0
    for i in range(n):
        if i%100 == 0:
            print(f"i{i} n{n}") 
        indice_high=A[i].indices.tolist()
        indice_low=B[i].indices.tolist()
        
        new_ind_high = np.array(indice_high[1:])
        new_ind_low = np.array(indice_low[1:])

        inter=np.intersect1d(new_ind_high, new_ind_low)
        card=len(inter)
        total_inter= total_inter + card
    

    entourage_score = total_inter / (n*k)
    print(f"Entourage : {entourage_score}")

    return entourage_score
    
    #l_high =[]
    #for row in A : 
     #   print(f"indice row de A : {row}")
      #  print(f"indice row de A : {row.indices}")
       # indices_high = row.indices.tolist()
        #new_ind_high = np.array(indices_high[1:])
        #print(f"new indice est : {new_ind_high}")
#        l_high.append(new_ind_high)
 #       print(f"La liste est : {l_high}")
 #       lh = np.array(l_high)
  #      print(f"l'array est : {lh}")
    

 #   l_low = []
#    for row in B : 
  #      print(f"indice row de B : {row}")
  #      print(f"indice row de B : {row.indices}")
  #      indices = row.indices.tolist()
   #     new_indices = np.array(indices[1:])
   #     print(f"new indice est : {new_indices}")
   #     l_low.append(new_indices)
  #      print(f"La liste est : {l_low}")
  #      l = np.array(l_low)
  #  print(lh.size)

   # total_inter=0
  #  for i in range(len(lh)):
  #      inter = np.intersect1d(lh[i],l[i])
  #      print(f"l'intersection est : {inter}")
  #      card=len(inter)
  #      print(f"le cardinal est : {card}")
  #      total_inter=total_inter + card
  #      print(f"L'intersection totale est : {total_inter}")
        
  #  n=len(X)
  #  k=k_neighbors
  #  entourage_score = total_inter / (n*k)
  #  print(f"Entourage : {entourage_score}")



    

    # ligne par ligne
        # la liste ordonnée des plus proche voisins
    # cree une matrice avec listes A'

    # A - B

    # Faut que je regarde dans la citation 6 ce qu'il faut faire 
    # Regarder entourage aussi 

    # metric_value = roc_auc_score(annotation,ref_annotation)
    #metric_value = entourage_score
    
    #return metric_value