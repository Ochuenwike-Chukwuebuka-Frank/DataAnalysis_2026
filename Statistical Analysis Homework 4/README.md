Project Description:

This project covers correlation analysis and regression analysis techniques applied to a clinical 
dataset. The work included normality testing, Spearman correlation with permutation-based 
significance assessment, multiple regression model fitting and comparison using BIC, and logistic 
regression to predict a binary outcome using hormone variables.

Data Description:

The dataset used is data_for_analysis, derived from the cleaned and imputed dataset produced in 
Practice 1. It contains 1,148 records and 31 variables, covering hormone variables (hormone1 to 
hormone10_generated), lipid variables, antioxidant indices, carbohydrate metabolism, and 
outcome/factor variables. The outcome variable indicates presence (1) or absence (0) of a tumour, 
with 160 positive cases and 987 negative cases.

R Environment:
IDE: RStudio (Mac)
R version: 4.5.3 (updated during practice)

Key Packages Used:
- wPerm       — for permutation-based correlation significance testing
- pROC        — for ROC curve and AUC computation
- ggplot2     — for visualization


Task 1 — Perform Correlation Analysis Between Hormone Variables:

Since all hormone variables were confirmed as non-normally distributed in Practice 3 
(Shapiro-Wilk, p < 0.05), Spearman rank correlation was selected as the appropriate method. 
A full Spearman correlation matrix was computed for all 9 hormone variables 
(hormone1 to hormone10_generated). The strongest correlation was found between hormone3 and 
hormone4 (rho = 0.584), indicating a moderate positive relationship. Most other hormone pairs 
showed weak or negligible correlations.

_____________

Task 2 — Obtain a Table with Correlation Coefficients with Significance Assessment 
(Permutation Method):

A permutation test (1,000 permutations) was applied to all 36 hormone pairs using the 
perm.relation() function from the wPerm package. Results were stored in a table containing 
the Spearman rho and permutation-based p-value for each pair.

Significant pairs (p < 0.05) included:
- hormone3 & hormone4: rho = 0.584 (strongest)
- hormone4 & hormone7: rho = -0.325 (negative)
- hormone5 & hormone8: rho = 0.259
- hormone3 & hormone5: rho = 0.241
- hormone7 & hormone10_generated: rho = 0.200

Output file: task1_2_hormone_correlation_permutation.csv

____________________

Task 3 — Perform Regression Analysis Between Hormone Variables:

The hormone pair with the strongest Spearman correlation was selected for regression analysis: 
hormone3 (predictor) and hormone4 (outcome), with rho = 0.584 (p = 0.002).

Five regression models were fitted:
- Linear regression
- 2nd-degree polynomial regression
- 3rd-degree polynomial regression
- Exponential regression (log-transformed response)
- Logarithmic regression (exp-transformed response)

The exponential model achieved the highest R-squared (0.297), while the logarithmic model 
performed very poorly due to extreme values produced by exp(hormone4).

________________________

Task 4 — Select the Best Model (BIC):

BIC was used to compare all five regression models. Lower BIC indicates a better model.

Results:
- Exponential model: BIC = 2049.49 (BEST)
- 2nd-degree polynomial: BIC = 8476.35
- 3rd-degree polynomial: BIC = 8477.46
- Linear model: BIC = 8525.91
- Logarithmic model: BIC = 296869.92 (WORST)

Conclusion: The exponential model (log(hormone4) ~ hormone3) is the best model, confirmed 
by both the lowest BIC and the highest R-squared (0.297).

Output file: task4_hormone_BIC_comparison.csv

______________

Task 5 — Fit Logistic Regression Models Using Hormone Variables to Predict Binary Outcome, 
Compare Model Performance via AIC/BIC, Compute Odds Ratios, and Summarize the Results:

Three binomial logistic regression models were fitted to predict the binary outcome 
(0 = no tumour, 1 = tumour) using hormone variables:
- model_logit_h1: single predictor (hormone1)
- model_logit_h2: two predictors (hormone3 and hormone4)
- model_logit_hall: all nine hormone variables

AIC/BIC Comparison:
- model_logit_h1:   AIC = 928.83, BIC = 938.92
- model_logit_h2:   AIC = 929.69, BIC = 944.82
- model_logit_hall: AIC = 919.31, BIC = 969.76 (lowest AIC - best by AIC)

Stepwise variable selection (AIC, direction = "both") reduced the full model to four 
predictors: hormone1, hormone2, hormone5, and hormone8, with AIC = 913.12.

hormone8 was identified as the strongest and only statistically significant predictor 
(p = 0.00122), with an odds ratio of 0.996 (95% CI: 0.994–0.998), indicating that higher 
hormone8 levels are associated with lower odds of tumour presence.

ROC Curve and AUC:
The full hormone model achieved an AUC of 0.6262, indicating moderate discriminative ability 
above random chance (AUC = 0.5).

Confusion Matrix:
The model predicted all cases as non-tumour (class 0), reflecting the class imbalance in the 
dataset (987 non-tumour vs 160 tumour cases).

Output files:
- task5_logistic_AIC_BIC.csv
- task5_confusion_matrix.csv
- task5_odds_ratios.csv
- ROC_Curve_Hormone_logistic_regression.png
