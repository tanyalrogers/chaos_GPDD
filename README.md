# Chaos in the Global Population Dynamics Database
This repository contains the data and code associated with the paper "Chaos is not rare in natural ecosystems" by Tanya L. Rogers, Bethany J. Johnson, and Stephan B. Munch.  

We employed 6 chaos detection methods:  
1. Direct LE estimaton (DLE)
2. Jacobian LE estimaton (JLE)
3. Recurrence Quantification Analysis (RQA)
4. Permutation Entropy (PE)
5. Horizontal Visibility Graphs (HVG)
6. Chaos Decision Tree (CDT)

All methods were applied to simulated data. The most reliable methods (JLE, RQA, PE) were applied to empirical data in the Global Population Dynamics Database (GPDD).

DLE and JLE are implemented in R. For JLE, the code uses rEDM version 0.7.4. Newer versions of the package are not backwards compatible.  

The remaining methods are implemented in MATLAB. The RQA code requires access to the CRP toolbox which must be requested [here](https://tocsy.pik-potsdam.de/CRPtoolbox/).


## Code files
Directory | Name | Description
-------- | ------- | -------
Methods | `LE_ChaosDetectionMethods.R` | Functions for E and tau selection, JLE, and DLE.
Methods | `ChaosClassification_MethodsB3toB6.m` | Function applying RQA, PE, HVG, and CDT methods. Code for individual methods are included in subdirectories.
Simulations | `Simulated_Models/` | Contains code for individual models used in simulations. See Supplementary Materials in the paper for more detail and the parameter values we used. Note that exact outputs will vary depending on the random seed. The simulated datasets we used in the paper are provided in the `data` folder. 
Simulations | `ChaosClassification_Example.m` | Simple example of how to generate a simulation and apply RQA, PE, HVG, and CDT methods.
Simulations | `Simulation_run.R` | Applies JLE and DLE to simulated datasets. Also obtains E and tau for use in other analyses. Note: This file takes >24 hours to run in full.
Simulations | `Simulation_run_MethodsB3toB6.m` | Applies RQA, RE, HVG, and CDT methods to simulated data. Note: This file takes >10 hours to run in full.
Simulations | `Simulation_plot.R` | Loads simulation results and generates plots.
Simulations | `ggplot_themes_rogers.R` | Custom ggplot themes for plotting.
GPDD | `GPDD_stability_dataprocessing.R` | Obtains and filters GPDD time series from package rgpdd and merges with information in `gpdd_lifehistory.csv`.
GPDD | `GPDD_stability_run.R` | Applies JLE to GPDD data. Also obtains E and tau for use in other analyses. 
GPDD | `GPDD_stability_run_RQA_PE.m` | Applies RQA and PE to GPDD data.
GPDD | `GPDD_stability_plot.R` | Loads GPDD results and generates plots.

## Data files
Name | Description
------- | -------
`simulation_dataset_embedding.csv` | Simulated *embedding* dataset.
`simulation_dataset_test.csv` | Simulated *test* dataset.
`simulation_dataset_validation.csv` | Simulated *validation* dataset.
`sims_test_E.csv` `sims_test_tau.csv` | E and tau values for test dataset.
`sims_validation_E.csv` `sims_validation_tau.csv` | E and tau values for validation dataset.
`sims_test_results.csv` | JLE and DLE results, test dataset.
`sims_test_results_othermethods.csv` | RQA, PE, HVG, and CDT results, test dataset.
`sims_validation_results.csv` | JLE and DLE results, validation dataset.
`sims_validation_results_othermethods.csv` | RQA, PE, HVG, and CDT results, validation dataset.
`gpdd_lifehistory.csv` | Life history and other additional metadata for GPDD species used.
`gpdd_timeseries.csv` | Cleaned GPDD time series in long format.
`gpdd_ts_metadata.csv` | All metadata for cleaned GPDD time series.
`gpdd_Etau_smap.csv` | E and tau values for GPDD time series.
`gpdd_results_smap.csv` | JLE results for GPDD time series.
`gpdd_results_truncation_smap.csv` | JLE results for GPDD time series (truncated chaotic series).
`gpdd_results_othermethods.csv` | RQA and PE results for GPDD time series.
`AndersonGilloolyLEdata.csv` | Empirical LE data from Anderson and Gillooly 2020.
