### **Data-driven prioritisation of high-risk individuals for weight loss interventions**

This repository contains the code used for the analyses in the paper (doi: XXX, (unpublished)). Functions and code build up on Carrasco-Zanini, J. et al. Proteomic signatures improve risk prediction for common and rare diseases. Nat Med. 30, 2489–2498 (2024).


### Project Structure

| **Script**           | **Description** |
|----------------------|-----------------|
| [`01_imputation.R`](01_prep_imputation/01_imputation.R)                        | Handles missing data imputation. |
| [`01_submit_imputation.sh`](01_prep_imputation/01_submit_imputation.sh)        | SLURM submission script for imputation. |
| [`01.2_post_imputation_assembly.R`](01_prep_imputation/01.2_post_imputation_assembly.R)        | Assemble imputed files. |
| [`coxnet_ridge_optim_boot.R`](functions/coxnet_ridge_optim_boot.R) | Function to run model optimization, and testing. |
| [`01_feature_selection_all_stepwise.R`](02_feature_selection/01_feature_selection_all_stepwise.R) | Performs stepwise multidomain feature selection. |
| [`01_feature_selection_all_submit.sh`](02_feature_selection/01_feature_selection_all_submit.sh) | SLURM submission script for stepwise multidomain feature selection. |
| [`01.2_feature_selection_all_individual.R`](02_feature_selection/01.2_feature_selection_all_individual.R) | Performs top 20 individual domain feature selection. |
| [`01.2_feature_selection_all_individual_submit.sh`](02_feature_selection/01.2_feature_selection_all_individual_submit.sh) | SLURM submission script top 20 individual domain feature selection. |
| [`01.3_feature_selection_all_individual_nonrandom.R`](02_feature_selection/01.3_feature_selection_all_individual_nonrandom.R) | Performs non-random individual domain feature selection. |
| [`01.3_feature_selection_all_individual_nonrandom_submit.sh`](02_feature_selection/01.3_feature_selection_all_individual_nonrandom_submit.sh) | SLURM submission script non-random individual domain feature selection. |
| [`02_LOFO.R`](02_feature_selection/02_LOFO.R)                  | Performs Leave-One-Feature-Out importance analysis. |
| [`02_LOFO_submit.sh`](02_feature_selection/02_LOFO_submit.sh)  | SLURM submission script for LOFO analysis. |
| [`03_shared_model.R`](02_feature_selection/03_shared_model.R)  | Trains and tests shared model across outcomes. |
| [`01_cumulative_incidence_predictors.R`](03_analysis/01_cumulative_incidence_predictors.R) | Calculates and visualizes cumulative incidences and predictors considered. |
| [`02_unimodal_performance.R`](03_analysis/02_unimodal_performance.R) | Visualizes unimodal model performances. |
| [`03_multimodal_performance.R`](03_analysis/03_multimodal_performance.R) | Visualizes multimodal model performances. |
| [`04_shared_model_performance.R`](03_analysis/04_shared_model_performance.R) | Computes and visualizes performance metrics for the shared model. |
| [`05_risk_stratification.R`](03_analysis/05_risk_stratification.R) | Computes and visualizes AR, Risk distribution, FDR/DR, NNT, pre-/post-test probabilities. |
| [`01_AdditionalParameters.R`](04_additional_analysis/01_AdditionalParameters.R) | Calculates number of deaths occuring during follow-up |
| [`02_Calibration.R`](04_additional_analysis/02_Calibration.R) | Generates calibration plots and metrics. |
| [`03_subgroup_analyses.R`](04_additional_analysis/03_subgroup_analyses.R) | Plots performances of OBSCORE in subgroup analyses. |
| [`04_time_dep_ROC.R`](04_additional_analysis/04_time_dep_ROC.R) | Time-dependent ROC analyses and performances of OBSCORE. |
| [`05_earlyexclusions.R`](04_additional_analysis/05_earlyexclusions.R) | Effect of extended exclusions of early incident cases on performance. |
| [`06_whr_whtr_performance.R`](04_additional_analysis/06_whr_whtr_performance.R) | Comparison of performance of OBSCORE with WHR vs. with WHtR |
| [`07_performance_mini.R`](04_additional_analysis/07_performance_mini.R) | Performance of OBSCORE with down to 10 features. |
| [`08_correlation_matrix_features.R`](04_additional_analysis/08_correlation_matrix_features.R) | Correlation matrix of OBSCORE features. |
| [`01_standardize_compute_score.R`](05_external_data/01_standardize_compute_score.R) | Standardizes predictors in SMT-1. |
| [`02_delta_endpoints.R`](05_external_data/02_delta_endpoints.R) | Calculates changes in endpoints. |
| [`03_filter_stratify.R`](05_external_data/03_filter_stratify.R) | Creates complete data and stratifies baseline risks. |
| [`04_linear_regressions.R`](05_external_data/04_linear_regressions.R) | Calculates adiposity measure change across risk groups. |
| [`05_delta_scores.R`](05_external_data/05_delta_scores.R) | Calculates changes in predicted risk. |
| [`06_plot_baseline_scores.R`](05_external_data/06_plot_baseline_scores.R) | Plots baseline risk across treatment arms. |
| [`07_plot_change_in_score.R`](05_external_data/07_plot_change_in_score.R) | Plots changes in score. |
| [`08_SMT1_Change_AdiposityMeasures.R`](05_external_data/08_SMT1_Change_AdiposityMeasures.R) | Publication ready plot adiposity measure changes across risk quartiles. |
| [`09_SMT1_ChangesRiskScore.R`](05_external_data/09_SMT1_ChangesRiskScore.R) | Publication ready plot changes in risk after treatment. |
| [`10_GH_val.R`](05_external_data/10_GH_val.R) | External validation in G&H. |
| [`11_EPICN_val.R`](05_external_data/11_EPICN_val.R) | External validation in EPIC-Norfolk. |

#### Note
*The provided scripts are not designed to work out of the box, but rather illustrate the main steps done for each of the analysis performed.*

#### Citation
*XXX, (unpublished)*
