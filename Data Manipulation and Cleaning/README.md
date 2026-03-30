1. Project Description:
   
	This project covers data manipulation and cleaning techniques applied to a clinical dataset. The work included exploratory data analysis, missing value identification and visualisation, removal of high missingness variables, imputation using multiple methods, and outlier detection.

Three homework tasks were completed:

  Little's MCAR Test to determine the missing data mechanism
  MICE imputation using Random Forest and PMM methods
  Outlier detection using boxplots and the LOF algorithm

_______________________________________________
_______________________________________________

2. Data Description:

	The dataset was provided by my lecturer and contains clinical and medical measurements across 1,148 records and 41 variables.
	It covers several variable groups including hormone variables (hormone1 to hormone14), lipid variables (lipids1 to lipids5, lipid_pero1 to lipid_pero5), antioxidant indices (antioxidant1 to antioxidant5), carbohydrate metabolism, and outcome/factor variables (outcome, factor_eth, factor_h, factor_pcos, factor_prl).
The dataset contained missing values across multiple variables, which were addressed through imputation as part of this project.

_______________________________________________
_______________________________________________

3. R Environment:
	IDE: RStudio (Mac)
	R version: 4.5.3

	Key Packages Used:

	dplyr – for data manipulation
	tidyr – for data reshaping
	ggplot2 – for visualization
	skimr – for detailed data summaries
	visdat – to visualize missing data
	naniar – for missing data analysis and MCAR testing
	mice – for multiple imputation
	dbscan – for Local Outlier Factor (LOF) calculations

_______________________________________________
_______________________________________________

4. Procedures:

- Data Exploration:
  
The dataset was first explored using str() to check variable types and skim() from the skimr package to get a summary of distributions across all variables.

- Dataset Preparation:
  
The dataset was split into two parts. MD_df kept most of the variables for analysis, while factor_df held the outcome and categorical variables separately.

- Missing Value Analysis:
  
Missingness was calculated per column as a percentage. Variables with more than 35% missing data were flagged for removal, while those at or below 35% were kept for imputation.

- Visualising Missingness:
  
Two plots were produced using the visdat and naniar packages. vis_miss() gave a visual map of where NAs appear across the dataset, and gg_miss_var() showed missingness per variable as a barplot.

- Removing High Missingness Columns:

Variables hormone9 to hormone14 exceeded the 35% threshold and were removed, producing the cleaned dataset handle_MD_df used in all three tasks.

_______________________________________________

Task 1 — Little's MCAR Test:

MCAR (Missing Completely At Random) means missing values occur purely by chance, with no relationship to any variable in the dataset. Little's MCAR test was run to determine the missing data mechanism, as this decides which imputation method is appropriate.

Hypotheses:
	•	H₀: Data is MCAR
	•	H₁: Data is not MCAR (MAR or MNAR)
	
Result:

  Statistic: 1809
  df: 1102
  p-value: 0
  Missing patterns: 54

Conclusion: 
	The p-value is effectively zero, so we reject H₀. The data is not MCAR, indicating missingness is systematic and likely MAR. 
	This justifies the use of MICE imputation in Task 2.
	
Output file: mcar_test_results.txt

_______________________________________________

Task 2 — Imputation with MICE:

  Two imputation methods were tested using the mice package:
  •	Random Forest (rf): uses decision trees to predict missing values based on relationships between variables.
  •	Predictive Mean Matching (pmm): fills missing values by matching them to observed values with similar predicted characteristics, keeping the original distribution intact.
	
Looking at the density plots, both methods did a reasonable job of preserving the distribution of hormone10_generated. 
PMM was the closer match because the imputed curve sits almost on top of the original. 
Random Forest drifted slightly in the 0.5–1.0 range. 

Conclusion: PMM is the better fit for this dataset.

Output Files:

•	Density_Plot_PMM.png and Density_Plot_Random_Forest.png (density plot comparisons)
•	imputed_dataset_rf.csv (final imputed dataset using Random Forest)
•	imputed_dataset_pmm.csv (final imputed dataset using PMM)

_______________________________________________

Task 3 — Outlier Detection:

Boxplots were used to identify outliers across numeric variables in imputed_handle_MD_df_final. 
Two plots were produced: one for lipid variables and one covering the full dataset.

Among the lipid variables, lipids2 and lipids4 showed the most notable outliers. 
Across the full dataset, hormone variables (hormone2, hormone3, hormone4, hormone5) showed the most extreme values, with several observations well above the interquartile range. 

Isolated outliers were also found in antioxidant4, hormone10_generated, and carb_metabolism.

Output Files:
Outlier_Detection_lipid_variables.png and Outlier_Detection_full_datasets.png

_______________________________________________

Task 3b — LOF Outlier Detection

The Local Outlier Factor (LOF) algorithm from the dbscan package was applied to detect multivariate outliers in imputed_handle_MD_df_final. 
Observations with LOF scores above 2 were flagged as outliers.

The histogram shows most observations scoring close to 1, with only a few exceeding 2 and one reaching nearly 5. 
The scatterplot confirms only 3 outliers were detected across lipids1 and lipids2, sitting within the main cluster rather than at the edges, showing LOF identifies outliers based on local density patterns across all variables, not just individual ones.

Output Files:
Histogram_of_LOF_scores.png and Bivariate_scatterplot_with_outliers.png



