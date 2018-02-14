# Structural Estimation of Behavioral Heterogeneity

by Zhentao Shi and Huanhuan Zheng, forthcoming at *Journal of Applied Econometrics*

We provide code of XMM and ELXM that replicates the empirical results.

## Environment

The numerical optimization depends on `NLPOT` and `RMOSEK`. The software
must be installed to run the code.

## Files

* `raw_data.csv`: the raw data set.
* `master1_data_cleaning.R`: generate the detrended time series for each period. After running this script, three Rdata files will be saved in the folder.
* `master2_estimation.R`: estimation by XMM and ELXM using the Rdata files generated above.

The following are the functions used in `master2_estimation.R`
* `CUE.R`: XMM
* `EL.R`: EXLM
* `func5.R`: collection of all the functions for the estimation procedures.
