# <ins>Statistics_Biological_data</ins>
## This respiratory contains scripts generated for statistical testing of biological data.
### 2 scripts are available - Statistics_analysis runs on a csv file with an index column and feature columns, it then performs multiple statistical testing and visualzation of the data. Batch_analysis takes as input a meta data file and a feature data file ( see file examples ) and then performs batch effect detection analysis to find batch paramters in the data.
### All output files will be saved in the main folder fro both scripts.
### Analysis steps for Statistics_analysis include : 
* Processing data- removing unwanted column prior to analysis and choosing index column for the analysis.
* Multi-group analysis, this step includes -
  - Outlier detection uisng mad/iqr.
  - Multiple regression anakysis.
  - Box plots for each feature.
  - Multi-group statistical testing using Anova with Tukey or Kruskal-Wallis with Dunns.
  - Manova analysis.
  - Bootstrap confidence interval analysis and visualization.
* Two-group analysis ( by user choice ), this step includes -
  - Bootstrap confidence interval analysis and visualization.
  - Variance testing ( Levene's test ) & Normality testing ( Shapiro-Wilk test ).
  - Two-group statistical testing using T-test, Wilcoxon or Mann-Whitney U.
  - Optional FDR correction & p value histograms.
  - Statistical summary with many paramters including cohens effect size calculation and plotting/
  - Power analysis.
  - Bayes factor analysis for each feature.
  - Permutation testing for each feature.
  - Logisitics regression for each feature and combination of features.
  - Linear regression for each feature and combination of features.
### Analysis steps for Batch_analysis include : 
* Processing data- choosing the row index for the feature data, choosing the treatment and batch parameters to test and data normalization. 
* Batch-effect analysis, this step includes -
  - PCA for batch detection.
  - BOX & Density plots for each feature.
  - RLE plots for each group in the treatment parameter groups.
  - Linear Regression analysis for batch effect detection.
  - Heatmap on all features and on selected features.
  - Variance calculation.
  - Batch effect type analysis. 
## References : 
- Cumming et al. The New Statistics: Why and How.
- Halsey et al. The fickle P value generates irreproducible results.
- Ho et al. Moving beyond P values: data analysis with estimation graphics.
- Krzywinski et al. Power and sample size.
- Cumming et al. Replication and p Intervals.
- Simpson et al. Package ‘permute’.
- Morey et al. Package ‘BayesFactor’.
- Wang et al. Managing batch effects in microbiome data.
- Leek et al. Tackling the widespread and critical impact of batch effects in high-throughput data.
