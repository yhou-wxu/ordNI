library(mvtnorm) # for using `qmvnorm` function.
##########################################################################################
# Compute correlation matrix of the test statistics.
##########################################################################################
Sigma <- function(ratio, v){
    # Input:
    # ratio: An (m+2) x 1 vector, where m is number of experimental treatments. 
    #   ratio[i] = n_{E_i} : n_R for i = 1, ..., m, ratio[m+1] = n_P : n_R, and 
    #   ratio[m+2] = n_R : n_R = 1, where n_i denotes sample size of treatment i, for 
    #   i = E_1, ..., E_m, P and R.
    # v: An (m+2) x 1 vector in which the element is sample size times the corresponding
    #   variance of the estimates of the treatment mean.
    #
    # Output:
    # Sig: An (m+1) x (m+1) matrix of correlation matrix.
    m <- length(ratio) - 2 # number of experimental treatments
    rho <- ratio[-(m+2)] * v[m+2] 
    rho <- rho / (rho + ratio[m+2] * v[-(m+2)])
    rho <- sqrt(rho)
    Sig <- outer(rho, rho)
    Sig[, m+1] <- -Sig[, m+1]
    Sig[m+1, ] <- -Sig[m+1, ]
    diag(Sig) <- 1
    return(Sig)
}

##########################################################################################
# Compute critical Value. 
# Notice that `qmvnorm` function in R package `mvtnorm` uses a random algorithm to compute
# quantiles of a multivariate normal distribution. To reduce the randomness, we calculate
# critical value `nrepeat` times and use their average as the final critical value.
##########################################################################################
CritVal <- function(Sig, alpha, nrepeat=1, ratio=NULL, v=NULL){
    # Input:
    # Sig: An (m+1) x (m+1) matrix of correlation matrix of the test statistics.
    # alpha: A scalar of significant level.
    # nrepeat: An integer of how many times to run `qmvnorm` function to compute quantile
    #   of the multivarite normal distribution. 
    # ratio: An (m+2) x 1 vector, where m is number of experimental treatments. 
    #   ratio[i] = n_{E_i} : n_R for i = 1, ..., m, ratio[m+1] = n_P : n_R, and 
    #   ratio[m+2] = n_R : n_R = 1, where n_i denotes sample size of treatment i, for 
    #   i = E_1, ..., E_m, P and R.
    # v: An (m+2) x 1 vector in which the element is sample size times the corresponding
    #   variance of the estimates of the treatment mean.
    #
    # Notice: (ratio, v) is only used to compute Sig if Sig is null.
    #
    # Output:
    # cv_mean: A scalar of average of the critical value for the `nrepeat` number of calculations.
    # cv_sd: A scalar of standard deviation of the critical value for the `nrepeat` number of calculations`.
    if(is.null(Sig)) Sig <- Sigma(ratio, v) 
    cv <- rep(0, nrepeat)
    for(i in 1:nrepeat){
        cv[i] <- qmvnorm(1-alpha, tail='lower.tail', corr=Sig)$quantile
    }
    list(cv_mean=mean(cv), cv_sd=sd(cv))
}



