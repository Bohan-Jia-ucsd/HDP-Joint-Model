# HDP-Joint-Model
This repository contains the R code used to implement the proposed **Hierarchical Dirichlet Process Joint Model** on CCTG 594 and ACTG 320 data and to reproduce the simulations presented in our study. The repository also includes implementations of several alternative models used for comparison.

## Real-Data Analysis

The `Real_Data_runs` directory contains the R files and datasets used for the real-data application.

* `data_cleaning.R`: prepares and combines the original CCTG 594 and ACTG 320 data used in the analysis.
* `HDP.R`: fits the proposed hierarchical Dirichlet process joint model.
* `HDP_POST.R`: fits the proposed hierarchical Dirichlet process joint model considering event-induced, structural change in the biomarker trajectory.
* `DP_ALL.R`: fits a Dirichlet process joint model to the combined HIV cohorts without modeling the subpopulation structure.
* `DP_EACH1.R` and `DP_EACH2.R`: fit separate Dirichlet process joint models to CCTG 594 and ACTG 320 cohorts.
* `dat_CD4_combine.rdata`: processed longitudinal data.
* `dat_surv_combine.rdata`: processed event data.

## Simulation Study

The `Simulation` directory contains the R files used to generate the simulated datasets, fit the competing models, and evaluate their performance.

### Data generation and model fitting

* `sim_from_dp_improved_newton_tarin.R` generates the training datasets under the simulation design.
* `sim_from_dp_improved_newton_test.R` generates the test datasets under the simulation design.

The following R files fit the candidate models across the simulated datasets:

* `run_100traindata_HDP.R`: proposed HDP joint model.
* `run_100traindata_DP_All.R`: DP joint model fitted to the combined subpopulation.
* `run_100traindata_DP_Each1.R` and `run_100traindata_DP_Each2.R`: separate DP joint models fitted to each subpopulation.
* `run_100traindata_JM.R`: conventional joint model.
* `run_100traindata_Surv.R`: survival model in the Hierarchical Dirichlet mixture model setting.

### Model evaluation

**Clustering performance**

* `ARI_train.R`: calculates the Adjusted Rand Index and Fuzzy Adjusted Rand Index for evaluating recovery of the latent class structure based on the training datasets.

**Dynamic prediction discrimination**

* `AUC_HDP_test.R`
* `AUC_DP_All_test.R`
* `AUC_DP_Each_test.R`
* `AUC_JM_test.R`
* `AUC_Surv_test.R`

These files calculate time-dependent AUC for the corresponding models based on the test datasets. `AUC_results.R` summarizes the resulting AUC estimates across simulation replicates.

**Dynamic prediction error**

* `Brier_HDP_test.R`
* `Brier_DP_All_test.R`
* `Brier_DP_Each_test.R`
* `Brier_JM_test.R`
* `Brier_Surv_test.R`

These files calculate the Brier score for the corresponding models based on the test datasets. `Brier_results.R` summarizes the results across simulation replicates.

**Goodness of fit**

* `GOF_train.R`: evaluates goodness of fit for the fitted models based on the training dataset.
