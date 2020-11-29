# Chaos in the Global Population Dynamics Database
This repository contains the data and code associated with the paper "Chaos is not rare in natural ecosystems" by Tanya L. Rogers, Bethany J. Johnson, and Stephan B. Munch.

## Code files
Name | Description
------- | -------
`GPDD_stability_dataprocessing.R` | Obtains and filters GPDD time series from package `rgpdd` and merges with information in `gpdd_lifehistory.csv`.
`GPDD_stability_functions.R` | Contains code for E and tau selection using s-map regression, Jacobian LE estimation, and direct LE estimation.
`GPDD_stability_run.R` | Applies analyses in `GPDD_stability_functions.R` to GPDD data.
`.mat` | Applies RQA and PE to GPDD data.
`GPDD_stability_plot.R` | Loads GPDD results and generates plots.
`.mat` | Generates simulated time series.
`Simulation_run.R` | Applies analyses in `GPDD_stability_functions.R` to simulated datasets. Note: This file takes >24 hours to run in full.
`.mat` | Applies RQA, RE, HVG, and CDT methods to simulated data. 
`Simulation_plot.R` | Loads simulation results and generates plots.
`ggplot themes rogers.R` | Custom ggplot themes.

## Data files
Name | Description
------- | -------
`gpdd_lifehistory.csv` | Life history and other additional metadata for GPDD species used.
`gpdd_timeseries.csv` | Cleaned GPDD time series in long format.
`gpdd_ts_metadata.csv` | All metadata for cleaned GPDD time series.
`gpdd_results_smap.csv` | Jacobian LE results for GPDD time series.
`gpdd_results_truncation_smap.csv` | Jacobian LE results for GPDD time series (truncated chaotic series).
`gpdd_othermethods.csv` | RQA and PE results for GPDD time series.
`AndersonGilloolyLEdata.csv` | Empirical LE data from Anderson and Gillooly 2020.
`KnownEandTau2.csv` | Simulated *embedding* dataset.
`ChaosMetaAnalysisSimulatedDataCORRECTED3.csv` | Simulated *test* dataset.
`ChaosMetaAnalysisSimulatedDataVALIDATION.csv` | Simulated *validation* dataset.
`sims_test_E.csv` `sims_test_tau.csv` | E and tau values for test dataset.
`sims_validation_E.csv` `sims_validation_tau.csv` | E and tau values for validation dataset.
`sims_test_results_allmethods.csv` | Results from all chaos detection methods, test dataset.
`sims_validation_results_allmethods.csv` | Results from all chaos detection methods, validation dataset.
