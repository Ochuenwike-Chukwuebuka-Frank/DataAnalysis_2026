## Project Description
Morphometric analysis of plant specimens from 6 collection sites, covering data standardization, correlation analysis, and dimensionality reduction via PCA.

## Data Description
Dataset (data_morphometry.txt) was provided by lecturer. It Contains 100 observations across 12 variables (11 morphometric traits and 1 grouping variable)  and collection sites: Г, И, Л, О, Х, Б. 
Traits cover shoot height, leaf dimensions (first and second leaf), perianth dimensions (outer and inner), and reproductive organ heights (stamen and pistil).

## R Environment
- **IDE:** RStudio
- **R version:** 4.6.0

## Key Packages
- `vegan` — range normalization
- `factoextra` — PCA biplots
- `plotly` — interactive 3D visualization

## Procedures

**Data Standardization**
All 11 traits scaled to [0, 1] using min-max normalization via `decostand()`.

**Spearman Correlation**
Pairwise Spearman rank correlations computed between all traits. Non-significant results (p > 0.05) set to zero.
Output: `correlation_table.csv`

**PCA**
Run on standardized data using `prcomp()`.
- Simple biplot with points colored by collection site → `simple_pca_biplot_no_grouping.png`
- Biplot with 95% confidence ellipses per group → `PCA_biplot_colored_pnts_95_confi_ellipses.png`
- Interactive 3D biplot (PC1, PC2, PC3) with loading vectors → `interactive_3D_scatter_plot_of_observations.png`

## Key Results
- PC1 (51%) captures overall plant size — driven by perianth length and shoot height
- PC2 (17.2%) contrasts leaf width against perianth length
- PC3 (10.3%) captures additional residual variation
- PC1 + PC2 + PC3 explain 78.5% of total variance (PC1 + PC2 alone: 68.2%)
- Groups overlap substantially with no clean morphological separation between sites
