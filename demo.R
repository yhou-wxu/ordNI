##########################################################################################
# This script use dataset of Example 1 in 
# "Multiple assessments of non-inferiority trials with ordinal endpoints" by Xu et al
# to demonstrate how to use our code to implement the method proposed in the above paper.
##########################################################################################
rm(list=ls())
source('FS.R') # fisher scoring algorithm for model estimation.
source('Delta.R') # compute delta and v ( = sigma^2 x delta).
source('CritVal.R') # compute critical value.
source('Power.R') # compute theoretical power.
source('OptSize.R') # compute optimal design confugaration.

##########################################################################################
# Argument initialization. 
# Notice that `qmvnorm` function in R package `mvtnorm` uses a random algorithm to compute
# quantiles of a multivariate normal distribution. To reduce the randomness, we calculate
# critical value `nrepeat` times and use their average as the final critical value.
# For illustration purpose, we set nrepeat = 1. To get the result of our paper, set
# nrepeat = 100.
##########################################################################################
alpha <- 0.025 # type 1 error, significant level.
m <- 2 # number of experimental treatment.
nrepeat <- 1 # number of compute critical values.
beta <- 0.2 # type 2 error.
Dlt1 <- 0.3 # assay sensitivity margin
Dlt2 <- 0.15 # NI margin.

##########################################################################################
# Dataset for Example 1, see Table 3 in Xu et al.
##########################################################################################
y.E1 <- c(rep(1, 2), rep(2, 7), rep(3, 21), rep(4, 20)) # response of treatment E1
y.E2 <- c(rep(1, 0), rep(2, 2), rep(3, 17), rep(4, 31)) # response of treatment E2
y.P <- c(rep(1, 20), rep(2, 22), rep(3, 6), rep(4, 2)) # response of treatment P
y.R <- c(rep(1, 7), rep(2, 13), rep(3, 18), rep(4, 12)) # response of treatment R
y <- c(y.E1, y.E2, y.P, y.R)
tr <- rep(1:(m+2), each=50) # treatment indicator vector.

##########################################################################################
# Model estimation and Critical value.
##########################################################################################
fit <- FS(tr, y) # model estimation
cat('tau is fixed at =', fit$tau, '\n')
cat('Estimated (mu_E1, mu_E2, mu_P, mu_R) =', fit$mu, '\n')
cat('Estimated (sigma_E1, sigma_E2, sigma_P, sigma_R) =', fit$sigma, '\n')

v <- Delta(fit$mu, fit$sigma, fit$tau)$v # compute v for obtaining critical value
Sig <- Sigma(ratio=rep(1, m+2), v=v) # correlation matrix of test statistics
set.seed(123)
cv <- CritVal(Sig, alpha, nrepeat)
cat('critical value =', cv$cv_mean, '\n')

##########################################################################################
# Compute optimal allocation ratio such that the power is maximized for a given total 
# sample size N = length(tr) = 150.  
##########################################################################################
Lambda1 <- fit$mu[1:m] + Dlt2 - fit$mu[m+2] 
Lambda2 <- fit$mu[m+2] - Dlt1 - fit$mu[m+1]
set.seed(123)
allo <- OptAllo(N=length(tr), v=v, Lambda1=Lambda1, Lambda2=Lambda2, alpha=alpha, nrepeat=nrepeat)
cat('The optimal allocation proportion nR : nE1: nE2 : nP = 1:', 
    allo$cE[1], ':', allo$cE[2], ':', allo$cP, '\n')
nR <- length(tr) / (1 + sum(allo$cE) + allo$cP)
nE <- length(tr) * allo$cE / (1 + sum(allo$cE) + allo$cP)
nP <- length(tr) * allo$cP / (1 + sum(allo$cE) + allo$cP)
cat('The corresponding (nR, nE1, nE2, nP) = (', nR, ',', nE[1], ',', nE[2], ',', nP, ')\n')
cat('The corresponding power =', allo$P, '\n')

# power of complete balance design
Sig <- Sigma(rep(1, m+2), v) # correlation matrix of test statistics
set.seed(123)
cv <- CritVal(Sig, alpha, nrepeat)$cv_mean
P.bal <- Power(length(tr), rep(1, m+2), Lambda1, Lambda2, cv, v, Sig)
cat('Power of complete balance allocation scheme =', P.bal$P, '\n')

##########################################################################################
# Compute optimal design configuration for a given level of power = 1 - beta
##########################################################################################
set.seed(123)
out <- Search(v=v, Lambda1=Lambda1, Lambda2=Lambda2, alpha=alpha, beta=beta, nrepeat=nrepeat)
cat('The minimu total sample size =', out$N, 'for a given level of power =', 1 - beta, '\n')
cat('The optimal allocation proportion nR : nE1: nE2 : nP = 1:', 
    out$cE[1], ':', out$cE[2], ':', out$cP, '\n')
nR <- length(tr) / (1 + sum(out$cE) + out$cP)
nE <- length(tr) * out$cE / (1 + sum(out$cE) + out$cP)
nP <- length(tr) * out$cP / (1 + sum(out$cE) + out$cP)
cat('The corresponding (nR, nE1, nE2, nP) = (', nR, ',', nE[1], ',', nE[2], ',', nP, ')\n')
cat('The corresponding power =', out$P, '\n')







