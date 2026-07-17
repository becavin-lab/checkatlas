# Figure Captions

## Main Figures

### Figure 1 — Overview heatmap of scaled metric scores per atlas

Heatmap of min-max scaled metric scores (0 = worst, 1 = best) for each atlas across four tasks: Cluster (purple), Annotation (red), Batch Correction (blue), and Dimensionality Reduction (green). Hatched cells indicate metrics not computed for that atlas. Scores are aggregated as the mean across all embeddings and label keys within each atlas.

### Figure 2 — Batch mixing vs. cell-type purity across atlases and embeddings

Scatter plot of iLISI (batch mixing, higher = better) versus 1 − cLISI (cell-type purity, higher = better) for each embedding in each atlas. Each point represents one embedding in one atlas, colored by embedding type. The top-right region represents ideal integration: batches are well-mixed while cell-type neighborhoods remain pure.

### Figure 3 — Scalability: runtime vs. atlas size

Log-log plot of total pipeline runtime (seconds) against number of cells. Circles indicate atlases evaluated on full gene space; squares indicate atlases where PCA subsampling was used (≥300k cells). Blue = full pipeline (with cluster metrics); red = partial pipeline (no cluster labels detected). Dashed and dotted lines show O(N) and O(N log N) reference slopes.

### Figure 4 — Dimensionality reduction quality across atlases and embeddings

Heatmap of scaled dimred metric scores per atlas × embedding. Columns are split into global structure metrics (Spearman ρ, distance correlation) and local neighborhood metrics (coKNN, entourage, LCMC, density preservation, Kruskal stress, Jaccard distance, graph edit distance). Higher scores indicate better preservation of the high-dimensional structure in the low-dimensional embedding.

### Figure 5 — Per-task and overall atlas quality scores

Bar chart showing scaled task-level scores (Cluster, Annotation, Batch Correction, Dimred) for each atlas. Diamonds show the unweighted mean overall score. Scores are min-max scaled per metric within each task. Atlases missing cluster or annotation metrics receive 0 for those tasks.

## Supplementary Figures

### Figure S1 — Cluster quality: mean silhouette width per embedding per atlas

Horizontal bar charts of mean silhouette width across all label keys for each embedding, per atlas. Positive values indicate compact, well-separated clusters; negative values indicate overlapping clusters. Embeddings are colored consistently with Figure 2.

### Figure S2 — Runtime breakdown per atlas per task

Stacked bar chart showing the time spent in each pipeline phase: dimred precomputation (distance matrices + kNN), other precomputation, cluster metrics, annotation metrics, batch correction metrics, and dimred metrics. Total runtime in seconds is annotated above each bar.

### Figure S3 — Inter-metric correlations within each task

Pearson correlation matrices of raw metric values within each task (Cluster, Dimred, Batch Correction). Diagonal entries are omitted. Values near +1 indicate redundant metrics; values near −1 indicate metrics capturing opposing properties.

### Figure S4 — Rare cell-type preservation: Isolated F1 scores per atlas

Distribution of isolated F1 scores across embeddings and label keys per atlas. The isolated F1 measures how well rare cell types (present in fewest batches) are captured as distinct clusters. The dashed line marks the scIB threshold of 0.50.

### Figure S5 — Principal Component Regression: batch predictability before and after scVI integration

Per-atlas comparison of PCR scores (lower = less batch signal in embedding) for PCA vs. scVI embeddings across the top batch keys. A reduction in PCR from PCA to scVI indicates successful batch effect removal by the integration method.
