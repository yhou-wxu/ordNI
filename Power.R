library(mvtnorm) # for using pmvnorm function.
##########################################################################################
# Translate sample ratio to sample size. 
##########################################################################################
ratio2n <- function(N, ratio){ 
    # N: total sample size.
    # ratio: An (m+2) x 1 vector, where m is number of experimental treatments. 
    #   ratio[i] = n_{E_i} : n_R for i = 1, ..., m, ratio[m+1] = n_P : n_R, and 
    #   ratio[m+2] = n_R : n_R = 1, where n_i denotes sample size of treatment i, for 
    #   i = E_1, ..., E_m, P and R.
    m <- length(ratio) - 2
    ratio.sum <- sum(ratio)
    n <- N * ratio[1:m] / ratio.sum
    n <- c(n, N * ratio[m+1] / ratio.sum, N * ratio[m+2] / ratio.sum)
    return(n)
}

##########################################################################################
# Compute testing power.
##########################################################################################
Power <- function(N, ratio, Lambda1, Lambda2, cv, v, Sig=NULL){
    # Input:
    # N: total sample size.
    # ratio: An (m+2) x 1 vector, where m is number of experimental treatments. 
    #   ratio[i] = n_{E_i} : n_R for i = 1, ..., m, ratio[m+1] = n_P : n_R, and 
    #   ratio[m+2] = n_R : n_R = 1, where n_i denotes sample size of treatment i, for 
    #   i = E_1, ..., E_m, P and R.
    # Lambda1: A (m+1) x 1 vector. Lambda1[i] = mu_{E_i} + Delta2 - mu_R for i = 1, ...., m, 
    #   where mu_j is treatment mean of treatment j for j = E_1, ..., E_m and R, and Delta2
    #   is NI margin.
    # Lambda2: = mu_R - mu_P - Delta1, where Delta1 is assay sensitivity margin.
    # cv: critical value.
    # v: An (m+2) x 1 vector in which the element is sample size times the corresponding
    #   variance of the estimates of the treatment mean.
    # Sig: An (m+1) x (m+1) matrix of correlation matrix of the test statistics, where m
    #   is number of experimental treatments. If Sig is not null, argument v need to be
    #   provided. Otherwise, argument v is useless.
    # Notice that (ratio, v) is only used to compute Sig if Sig is null.
    # Output:
    # Ps: Power for (Z_{E_l}, Z_P: l = 1, ..., m), whose minimum is Power. See Eq (16) in 
    #   Multiple assessments of non-inferiority trials with ordinal endpoints" by Xu et al.
    # P: power.
    m <- length(v) - 2
    n <- ratio2n(N, ratio)
    # means of Z_{E_l}(l = 1, ..., m) and Z_P,  where Z_i is the test statistics
    mu_Z <- c(Lambda1, Lambda2) / sqrt(v[-(m+2)] / n[-(m+2)] + v[m+2] / n[m+2]) 
    # correlations between Z_E and Z_P
    m <- length(n) - 2 # number of experimental treatments
    if(is.null(Sig)) Sig <- Sigma(ratio, v)
    rho <- Sig[1:m, m+1]
    # compute P(Z_E >  cv and Z_P > cv), which equals
    # P{(Z_E - mu_{Z_E}) > (cv - mu_{Z_E}) and (Z_P - mu_{Z_P}) > (cv - mu_{Z_P})} = 
    # P{-(Z_E - mu_{Z_E}) < (-cv + mu_{Z_E}) and -(Z_P - mu_{Z_P}) < (-cv + mu_{Z_P})}.
    u <- -cv + mu_Z
    Ps <- rep(0, m)
    for(l in 1:m){
        corr <- rbind(c(1, rho[l]), c(rho[l], 1))
        Ps[l] <- pmvnorm(lower=rep(-Inf, 2), upper=c(u[l], u[m+1]), corr=corr)
    }
    # compute Power
    P <- min(Ps)
    list(Ps=Ps, P=P)
}


