## Mutual Information 

## Description ##

La Mutual Information (MI) est une mesure externe fondée sur la théorie de l’information qui quantifie la dépendance entre deux partitions ou variables aléatoires.
Elle mesure la quantité d’information mutuellement partagée : en pratique, combien la connaissance de la partition U réduit l’incertitude sur la partition V, et vice versa.
Autrement dit, MI évalue la similarité entre un jeu de labels de référence (par exemple des types cellulaires annotés) et une partition calculée par un algorithme de clustering. Cette mesure est symétrique, non négative et s’appuie sur les entropies marginales et conjointe des partitions.
Une MI élevée signifie que les deux partitions partagent beaucoup d’information (les clusters coïncident bien avec les vrais types cellulaires), alors qu’une MI nulle indique qu’elles sont complètement indépendantes.

## Formules ##

Pour deux partitions $U$ et $V$ de $n$ objets (avec $R$ clusters dans $U$ et $C$ dans $V$), on définit les probabilités marginales 

$$p(i)=\frac{|U_i|}{n} \text{ et } p(j)\frac{|V_j|}{n}$$

,ainsi que la probabilité conjointe 

$$p(i,j)=\frac{|U_i \cup V_j|}{n}$$.

