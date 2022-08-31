# Code for model estimation and sample size determination in 
"Multiple assessments of non-inferiority trials with ordinal endpoints" by Xu et al. (2022)

# Content:
* demo.R: load core functions to analysis Example 1 in the aforementioned paper. 

* FS.R:
R code for fitting latent variable model via Fisher scoring algorithm.

* delta.R: 
R code for computing delta0, delta1, delta2, delta and v defined in the aforementioned paper. 

* CritVal.R
R code for computing critical value.

* Power.R
R code for computing testing power.

* OptSize.R
R code for determing sample size via L-BFGS-B algorithm and bisection algorithm, see Aloorithm 1 in the aforementioned paper.

* common.R
R code for common functions frequently used in the above R script files.

* cd.cpp
Rcpp functions for implementing coordinate descent and block coordinate descent algorithms.

* simdata.R
R code for generating simulated dataset.

