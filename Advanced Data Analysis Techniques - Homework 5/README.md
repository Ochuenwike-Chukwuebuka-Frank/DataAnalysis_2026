## Project Description
 
This project applies multivariate ordination and clustering techniques to a community ecology dataset. The work included NMDS ordination across three distance metrics, UPGMA hierarchical clustering, species vector fitting using envfit, and PERMANOVA to test cluster differences.
 
---
 
## Data Description
 
The dataset `data.txt` contains species abundance records from 6 sites and 34 species. Abundance values are semi-quantitative counts per site.
 
**Sites:** Ledyanaya Gora, Karaul village, Ladyginskie Yary, Sopochnaya Karga, Sibiryakov Island Srednee Lake, Chernyi Bay
 
---
 
## R Environment
 
- **IDE:** RStudio
- **R version:** 4.4.0
- **Key package:** `vegan`
---


## Task 1 — Comparison of Distance Metrics (NMDS and UPGMA)
 
NMDS ordination and UPGMA clustering were performed using three distance metrics:
 
- **Euclidean** — quantitative distances, sensitive to double zeros
- **Bray–Curtis** — standard ecological dissimilarity, robust to double zeros
- **Jaccard** — presence/absence-weighted dissimilarity
For each metric, `metaMDS()` produced ordination plots and `vegdist()` + `hclust(method = "average")` produced dendrograms in two styles: labels aligned (`hang = -1`) and default R style.
 
All three metrics produced identical clustering topology — a consistent 2-group split across all methods.
 
**Output files:**
- `Task1_Euclidean_NMDS_plot.png`, `Task1_UPGMA_Euclidean_labels_aligned.png`, `Task1_UPGMA_Euclidean_standard.png`
- `Task1_Bray_Curtis_NMDS_plot.png`, `Task1_UPGMA_BrayCurtis_labels_aligned.png`, `Task1_UPGMA_BrayCurtis_standard.png`
- `Task1_Jaccard_NMDS_plot.png`, `Task1_UPGMA_Jaccard_labels_aligned.png`, `Task1_UPGMA_Jaccard_standard.png`
---


## Task 2 — Detailed Bray–Curtis Analysis
 
### 2.1 — NMDS
NMDS computed using `metaMDS(distance = "bray", trymax = 50)`. Stress approached zero due to small sample size (n = 6); results interpreted with caution.
 
### 2.2 — Significant Species Vectors (envfit)
`envfit()` applied with 1,000 permutations. No significant species vectors were found (p ≤ 0.05), likely due to the small sample size.
 
### 2.3 — UPGMA Clustering
Bray–Curtis dissimilarity matrix computed via `vegdist()` and clustered with `hclust(method = "average")`.
 
### 2.4 — Cluster Extraction (cutree)
Dendrogram cut into k = 2 clusters:
- **Cluster 1:** Ledyanaya Gora, Ladyginskie Yary, Sopochnaya Karga, Chernyi Bay
- **Cluster 2:** Karaul village, Sibiryakov Island Srednee Lake
### 2.5 — Full NMDS Plot
Final plot produced with cluster-coloured points, 95% confidence ellipses (`ordiellipse`), non-overlapping site labels (`orditorp`), and significant species arrows (none plotted).
 
**Output file:** `Task2.5_NMDS_BrayCurtis_full.png`
 
---


## Task 3 — PERMANOVA
 
PERMANOVA run using `adonis2(data ~ clusters, method = "bray", permutations = 999)`.
 
| | Df | SumOfSqs | R² | F | Pr(>F) |
|---|---|---|---|---|---|
| Model | 1 | 0.49182 | 0.4088 | 2.7661 | 0.0667 |
| Residual | 4 | 0.71121 | 0.5912 | | |
| Total | 5 | 1.20303 | 1.0000 | | |
 
- **R² = 0.4088** — cluster membership explains ~41% of species composition variation
- **p-value = 0.0667** — not significant at the 0.05 threshold
**Conclusion:** Differences between clusters are not statistically significant (p ≥ 0.05), likely due to limited statistical power with only 6 sites.
 
**Output files:** `Task3.2_PERMANOVA_adonis_results_console_output.png`, `Task3.3_PERMANOVA_summary_console_output.png`
 
