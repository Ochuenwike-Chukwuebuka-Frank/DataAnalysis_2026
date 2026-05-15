Project Description:

This project covers statistical analysis techniques applied to a clinical dataset. The work included 
descriptive statistics, normality and variance testing, hypothesis testing using multiple statistical 
tests, and exploratory correlation analysis.

Data Description:

The dataset used is data_for_analysis, derived from the cleaned and imputed dataset produced in 
Practice 1. It contains 1,148 records and 31 variables, covering hormone variables (hormone1 to 
hormone10_generated), lipid variables, antioxidant indices, carbohydrate metabolism, and 
outcome/factor variables. The outcome variable indicates presence (1) or absence (0) of a tumour, 
with 160 positive cases and 987 negative cases.

R Environment:
IDE: RStudio (Mac)
R version: 4.5.3

Key Packages Used:
- gtsummary  — for publication-ready descriptive statistics tables
- car        — for Levene's Test for homogeneity of variance
- lawstat    — for the Brunner-Munzel test
- ggplot2    — for correlation heatmap visualisation
- reshape2   — for data reshaping for heatmaps
- DataExplorer — for automated EDA report generation


Task 1 — Create a Table with Descriptive Statistics by Group (All Hormones):

A descriptive statistics table was created for all hormone variables grouped by outcome using 
tbl_summary() from the gtsummary package. Since all hormones were found to be non-normally 
distributed (confirmed by Shapiro-Wilk test), parameters are reported as median and IQR. 
A Wilcoxon-based p-value was added to the table to indicate whether differences between groups 
were statistically significant.

Output file: task1_hormone_descriptive_stats_table.csv

________________

Task 2 — Levene's Test for Homogeneity of Variance and Shapiro-Wilk Test:

The Shapiro-Wilk test was performed for all hormone variables in both outcome groups. All hormones 
were found to be non-normally distributed (p < 0.05) in both groups. Levene's Test for Homogeneity 
of Variance was then performed, and all hormones showed equal variances between groups (p > 0.05).

Output files: shapiro_wilk_results.csv and levene_test_results.csv

_______________

Task 3 — Q-Q Plots and Histograms (All Hormones):

Histograms and Q-Q plots were generated for all hormone variables separately for each outcome group. 
All distributions showed visible deviation from normality, consistent with the Shapiro-Wilk results. 
Most hormones showed right-skewed distributions. hormone10_generated showed the most extreme 
deviation, with a near-zero spike in outcome group 0.

Output files: PNG files for histograms and Q-Q plots per hormone per group.
__________________

Task 4 — Brunner-Munzel Test, t.test, and Wilcox.test for 2 Independent Groups (All Hormones):

Three statistical tests were performed for each hormone variable comparing outcome group 0 
(no tumour) vs outcome group 1 (tumour): Brunner-Munzel test, Welch t-test, and Wilcoxon rank-sum 
test. P-values from all three tests were compiled into a single table alongside median (IQR) 
for each group.

Conclusion: The Wilcoxon rank-sum test is the most applicable test for this dataset, as all hormone 
variables are non-normally distributed (Shapiro-Wilk, p < 0.05) and have equal variances between 
groups (Levene's test, p > 0.05).

Output file: task4_descriptive_with_pvalues.csv

____________________

Task 5 — Correlation Heatmap for All Hormones by Group:

Spearman correlation was selected as the method for all heatmaps, since all hormone variables were 
found to be non-normally distributed. Correlation heatmaps were constructed separately for outcome 
group 0 (no tumour) and outcome group 1 (tumour) to compare hormone correlation patterns between 
the two groups. Correlation matrices were also exported as CSV files.

Output files: task5_correlation_heatmap_outcome0.png, task5_correlation_heatmap_outcome1.png,
task5_correlation_matrix_outcome0.csv, task5_correlation_matrix_outcome1.csv

