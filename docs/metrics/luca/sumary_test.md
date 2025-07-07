<table>
    <tr>
        <td>Name</td>
        <td>Source</td>
        <td>Upload</td>
        <td>Date</td>
        <td>Name on github</td>
        <td>Verification</td>
        <td>Date</td>
        <td>Type</td>
        <td>Summary</td>
        <td>Interval</td>
    </tr>
    <tr>
        <td>Adjusted Rand Index (ARI)</td>
        <td>Luecken  et. al // OP</td>
        <td>20/06/25</td>
        <td>adjusted_rand_index.md</td>
        <td>Batch integration</td>
        <td>Adjusted Rand Index compares clustering overlap, correcting for random labels and considering correct overlaps and disagreements.</td>
        <td>[-1,1] the higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Average Silhoeutte Width (ASW)</td>
        <td>Luecken  et. al // OP</td>
        <td>23/06/25</td>
        <td>asw.md</td>
        <td>Batch integration</td>
        <td>ASW Batch : Evaluates batch mixing quality after integration. // ASW Label : Evaluates cell type clustering preservation.</td>
        <td>[0,1] the higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>iLISI &amp; cLISI</td>
        <td>Luecken  et. al // OP</td>
        <td>23/06/25</td>
        <td>lisi.md</td>
        <td>Batch integration</td>
        <td>iLISI &amp; cLISI measure batch mixing and cell-type preservation quality in single-cell integration.</td>
        <td>The higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Cell Cycle Conservation</td>
        <td>Luecken  et. al // OP</td>
        <td>23/06/25</td>
        <td>cc_conservation.md</td>
        <td>Batch integration</td>
        <td>Measures preservation of cell-cycle biological signals after batch effect correction in integration</td>
        <td>[0,1] the higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Graph Connectivity</td>
        <td>Luecken  et. al // OP</td>
        <td>23/06/25</td>
        <td>graph_connectivity.md</td>
        <td>Batch integration</td>
        <td>Measures preservation of biological cell-type neighborhoods in integrated kNN graphs.</td>
        <td>[0,1] the higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>HVG overlap</td>
        <td>Luecken  et. al // OP</td>
        <td>24/06/25</td>
        <td>hvg.md</td>
        <td>Batch integration</td>
        <td>Measures how well biological gene variability is preserved after batch integration.</td>
        <td>[0,1] the higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Isolated Label F1 score</td>
        <td>Luecken  et. al // OP</td>
        <td>24/06/25</td>
        <td>f1_score.md</td>
        <td>Batch integration</td>
        <td>Measures how well rare or isolated cell types are preserved after batch integration.</td>
        <td>[0,1] the higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>kBET</td>
        <td>Luecken  et. al // OP</td>
        <td>24/06/25</td>
        <td>kbet.md</td>
        <td>Batch integration</td>
        <td>Quantifies batch mixing by testing local neighborhood batch label distributions against global proportions.</td>
        <td>[0,1] the higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>NMI</td>
        <td>Luecken  et. al // OP</td>
        <td>24/06/25</td>
        <td>nmi.md</td>
        <td>Batch integration</td>
        <td>NMI quantifies clustering accuracy by comparing predicted and true labels in integrated datasets.</td>
        <td>[0,1] the higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>PCR</td>
        <td>Luecken  et. al // OP</td>
        <td>24/06/25</td>
        <td>pcr.md</td>
        <td>Batch integration</td>
        <td>Quantifies how well continuous biological variation is preserved after data integration.</td>
        <td>[0,1] the higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Silhouette</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>silhouette.md</td>
        <td>Clustering</td>
        <td>Internal metric evaluating a trade-off between average inter and intra-cluster distances.</td>
        <td>[-1;1]</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Davies_Bouldin</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>dbi.md</td>
        <td>Clustering</td>
        <td>Internal Metric Based on a Similarity Notion</td>
        <td>The lower the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Dunn Index</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>dunn.md</td>
        <td>Clustering</td>
        <td>Internal metric evaluating a trade-off between characteristic inter-and intra-cluster distances.</td>
        <td>The higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>DBCV</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>dbcv.md</td>
        <td>Clustering</td>
        <td>Internal metric for evaluating density-based clustering methods.</td>
        <td>[-1;1] The higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Rand_index</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>rand.md</td>
        <td>Clustering</td>
        <td>External metric based on consistency between two clusterings</td>
        <td>[0,1] The higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Fowlkes_Mallows</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>fowlkes_mallows.md</td>
        <td>Clustering</td>
        <td>External metric comparing agreement and disagreement rates between two clusterings</td>
        <td>[0,1] The higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Mutual_Information</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>mutual_information.md</td>
        <td>Clustering</td>
        <td>External metric based on joint entropy</td>
        <td>The higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>V-Measure</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>v_measure.md</td>
        <td>Clustering</td>
        <td>External metric based on the concepts of completeness and homogeneity</td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Mean Sqaured Error</td>
        <td>Batson et. al // OP</td>
        <td>26/06/25</td>
        <td>mse.md</td>
        <td>Denoising</td>
        <td>« The mean squared error between the denoised counts of the training dataset and the true counts of the test dataset after reweighting by the train/test ratio. »</td>
        <td>The lower the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Poisson Loss</td>
        <td>Batson et. al // OP</td>
        <td>26/06/25</td>
        <td>poisson.md</td>
        <td>Denoising</td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Kruskal&#39;s stress</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>kruskal_stress.md</td>
        <td>Dimension reduction</td>
        <td>Metric based on the difference between distances in high and low dimensions.</td>
        <td>[0,1] The lower the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Spearman&#39;s rho</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>spearman_rho.md</td>
        <td>Dimension reduction</td>
        <td>Metric based on the difference between the orders of distances in high and low dimensions</td>
        <td>The lower the better [0, inf]</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>LCMC</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>lcmc.md</td>
        <td>Dimension reduction</td>
        <td>Metric evaluating changes in the nearest neighbours matrix after dimension reduction</td>
        <td>The higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Trustworthiness and continuity</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>trust_continuity.md</td>
        <td>Dimension reduction</td>
        <td>Metric based on notions of reliability and continuity linked to nearest neighbours</td>
        <td>[0,1] The higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Average Jaccard Distance (AJD)</td>
        <td>Cooley et. al</td>
        <td>02/07/25</td>
        <td>ajd.md</td>
        <td>Dimension reduction</td>
        <td>AJD quantifies how much local neighborhood structure is preserved after dimensionality reduction in scRNA-seq.</td>
        <td>[0,1] The lower the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Graph Edit Distance</td>
        <td>Cooley et. al</td>
        <td>02/07/25</td>
        <td>ged.md</td>
        <td>Dimension reduction</td>
        <td>GED quantifies graph dissimilarity by counting minimum edits needed to transform one graph into another.</td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Kolmogorov–Smirnov</td>
        <td>Pachter et al.</td>
        <td>03/07/25</td>
        <td>kolmogorov_smirnov.md</td>
        <td>Clustering</td>
        <td>The KS statistic is the maximum difference between two empirical CDFs; values near 1 signal highly dissimilar distributions.</td>
        <td>[0,1] The higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011288</td>
        <td>Pachter et al.</td>
        <td>Dimension reduction</td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Entourage</td>
        <td>Becavin et al.</td>
        <td>02/07/25</td>
        <td>entourage.md</td>
        <td>Dimension reduction</td>
        <td>Local conservation of entourage // Quantifies local neighborhood preservation quality during dimensionality reduction for structural fidelity assessment.</td>
        <td>[0,1] The higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Trustworthiness</td>
        <td>Venna &amp; Kasi // OP</td>
        <td>26/06/25</td>
        <td>trust_continuity.md</td>
        <td>Dimension reduction</td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Continuity</td>
        <td>Zhang et. al. // OP</td>
        <td>26/06/25</td>
        <td>trust_continuity.md</td>
        <td>Dimension reduction</td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>LCMC</td>
        <td>Zhang et. al. // OP</td>
        <td>26/06/25</td>
        <td>lcmc.md</td>
        <td>Dimension reduction</td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>co-KNN (AUC &amp; size)</td>
        <td>Zhang et. al. // OP</td>
        <td>26/06/25</td>
        <td>co_knn.md</td>
        <td>Dimension reduction</td>
        <td>These metrics evaluate how well local neighborhood structures are preserved after dimensionality reduction</td>
        <td>The higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Density Preservation</td>
        <td>Narayan et al. // OP</td>
        <td>26/06/25</td>
        <td>density_preservation.md</td>
        <td>Dimension reduction</td>
        <td>Density Preservation measures how well local point densities are maintained after dimensionality reduction for accurate visualization.</td>
        <td>[-1,1] the higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Distance correlation</td>
        <td>Schober et. al. // OP</td>
        <td>26/06/25</td>
        <td>distance_correlation.md</td>
        <td>Dimension reduction</td>
        <td>Distance Correlation detects both linear and nonlinear dependencies, equaling zero only when variables are truly independent.</td>
        <td>[0;1]</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Local &amp; Global co-KNN</td>
        <td>Zhang et. al. // OP</td>
        <td>27/06/25</td>
        <td>local_global_coknn.md</td>
        <td>Dimension reduction</td>
        <td>Global &amp; Local co-KNN measures how well both local and global neighborhood structures are preserved after dimensionality reduction.</td>
        <td>[0,1] the higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Coefficient of determination</td>
        <td>Miles et al. // OP</td>
        <td>26/06/25</td>
        <td>r2.md</td>
        <td>Spatial decomposition</td>
        <td>Evaluates how well gene expression variability is explained by models in single-cell data.</td>
        <td>[0,1] the higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Kendall&#39;s correlation</td>
        <td>Kendall et. al. // OP</td>
        <td>26/06/25</td>
        <td>tau.md</td>
        <td>Spatially variable genes</td>
        <td>[0;1] or [-1;1]</td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>One-vs-All</td>
        <td>Report AC</td>
        <td>17/06/25</td>
        <td>ova_ovm.md</td>
        <td>Specificity</td>
        <td>Gene expression in one cell type versus the rest.</td>
        <td>[0,1] The higher the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>One-vs-Max</td>
        <td>Report AC</td>
        <td>17/06/25</td>
        <td>ova_ovm.md</td>
        <td>Specificity</td>
        <td>Gene expression in one cell type vs the second more highly expressed cell type.</td>
        <td>[0,+∞[</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Shannon Entropy</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>shannon.md</td>
        <td>Specificity</td>
        <td>Entropy of gene distribution across cell types.</td>
        <td>[0;1]</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Gini coefficient</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>gini.md</td>
        <td>Specificity</td>
        <td>Evaluates the inequality of the gene&#39;s distribution.</td>
        <td>The lower the better</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>Kendall&#39;s tau</td>
        <td>Report AC</td>
        <td>16/06/25</td>
        <td>tau.md</td>
        <td>Specificity</td>
        <td>Similar to Gini and Shannon, but less restrictive.</td>
        <td>[0;1] or [-1;1]</td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
    <tr>
        <td>ZADU</td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
        <td></td>
    </tr>
</table>
