library(mvtnorm)
##########################################################################################
# Transform `Power` function in R script Power.R into a new form to use `optim` function.
# Notice that `qmvnorm` function in R package `mvtnorm` uses a random algorithm to compute
# quantiles of a multivariate normal distribution. To reduce the randomness, we calculate
# critical value `nrepeat` times and use their average as the final critical value.
##########################################################################################
Power4optim <- function(cEprimeP, N, v, alpha, Lambda1, Lambda2, nrepeat){
    # cEprimeP: A 2 x 1 vector. The first element in cEprimeP is cEprime and the second
    #   element is cP.
    m <- length(v) - 2
    cEprime <- cEprimeP[1]
    cP <- cEprimeP[2]
    cE <- cEprime * sqrt(v[1:m]) / min(sqrt(v[1:m]))
    ratio <- c(cE, cP, 1)
    Sig <- Sigma(ratio, v)
    cv <- CritVal(Sig, alpha, nrepeat)$cv_mean
    P <- Power(N, ratio, Lambda1, Lambda2, cv, v, Sig)$P
    return(-P)
}

##########################################################################################
# Compute the optimal allocation proportion cEprime and cP given total sample size by 
# L-BFGS-B algorithm implemented in R function optim, where cEprime and cP are allocation
# proportions such that n_R : n_{E_1} : ... : n_{E_m} : n_P = 1 : cE_1: ... : cE_m : cP with
# cE_i = cEprime x sqrt(v[i]) / min(sqrt(v[i])) for i = 1, ..., m, where m is number of 
# experimental treatments, and v[i] = n_{E_i} x var( \hat{ mu_{E_i} } ) for i = 1, ..., m. See 
# "Multiple assessments of non-inferiority trials with ordinal endpoints" by Xu et al for 
# detial information of v.
# Notice that `qmvnorm` function in R package `mvtnorm` uses a random algorithm to compute
# quantiles of a multivariate normal distribution. To reduce the randomness, we calculate
# critical value `nrepeat` times and use their average as the final critical value.
##########################################################################################
OptAllo <- function(N, v, Lambda1, Lambda2, init=NULL, lower=c(0.1, 0.1), upper=c(1, 0.8), alpha=0.025, nrepeat=1){
    # Input:
    # N: Total sample size.
    # v: An (m+2) x 1 vector in which the element is sample size times the corresponding
    #   variance of the estimates of the treatment mean.
    # Lambda1: A (m+1) x 1 vector. Lambda1[i] = mu_{E_i} + Delta2 - mu_R for i = 1, ...., m, 
    #   where mu_j is treatment mean of treatment j for j = E_1, ..., E_m and R, and Delta2
    #   is NI margin.
    # Lambda2: = mu_R - mu_P - Delta1, where Delta1 is assay sensitivity margin.
    # init: A 2 x 1 vector of initial values of cEprime and cP.
    # lower: A 2 x 1 vector. The first element is the lower bound of min( n_{E_l} ):n_R for 
    #   l = 1, ..., m. The second element is the lower bound of n_P : n_R.
    # upper: A 2 x 1 vector. The first element is the upper bound of max( n_{E_l} ):n_R for 
    #   l = 1, ..., m. The second element is the upper bound of n_P : n_R.
    # alpha: The significant level.
    # nrepeat: An integer of how many times to run `qmvnorm` function to compute quantile
    #   of the multivarite normal distribution. 
    #
    # Output:
    # P: the maximum of power given total sample size N.
    # cE: A vector of m x 1 vector of optimal allocation proportion of experimental treatment
    #   compared to treatment R.
    # cP: optimal allocation proportion compared to treatment R.
    # cEprime: optimal allocation proportion of experimental treatment such that 
    #   cE[i] = cEprime[i] * sqrt(v[i]) / min(sqrt(v[i])).
    # conv: = 0 if convergence is obtained.
    m <- length(v) - 2
    cEprime_max <- min(sqrt(v[1:m])) / max(sqrt(v[1:m])) 
    upper[1] <- upper[1] * cEprime_max
    if(is.null(init)) init <- c(0.8 * cEprime_max, 0.2) 
    fit <- optim(init, Power4optim, N=N, v=v, alpha=alpha, Lambda1=Lambda1, 
                 Lambda2=Lambda2, nrepeat=nrepeat, method='L-BFGS-B', 
                 lower=lower, upper=upper, control=list(trace=0))
    cEprime <- fit$par[1]
    cP <- fit$par[2]
    P <- -fit$value
    cE <- cEprime * sqrt(v[1:m]) / min(sqrt(v[1:m])) 
    # conv = 0 convergent.
    # conv = 1 iteration limit maxit had been reached.
    list(P=P, cE=cE, cP=cP, cEprime=cEprime, conv=fit$convergence)
}

##########################################################################################
# Backward and forward search to find initial smaple sizes Nl and Nu such that
# power at Nl is less than 1 - beta and power at Nu is greater than 1 - beta, where 
# 1 - beta is the specified power level.
##########################################################################################
BF <- function(beta, Nl, Nu, fl, fu, func, ...){
    while(fl$P > (1 - beta)){ # make sure Pl - (1 - beta) < 0
        Nl <- Nl - 10
        fl <- func(Nl, ...)
    }
    while(fu$P < (1 - beta)){ # make sure Pu - (1 - beta) > 0
        Nu <- Nu + 10
        fu <- func(Nu, ...)
    }
    list(Nl=Nl, Nu=Nu, fl=fl, fu=fu)
}

##########################################################################################
# Compute root of equation func = 1 - beta, where func is a function, using bisection algorithm
##########################################################################################
Bisec <- function(eps, beta, Nl, Nu, Pl, Pu, func, ..., trace=FALSE){
    # Input:
    # eps: convergence tolerance.
    # beta: the value that we want find the root of equation `func` = 1 - beta.
    # Nl: lower value in bisection algorithm.
    # Nu: upper value in bisection algorithm.
    # Pl: function value at Nl.
    # Pu: function value at Nu.
    # func: the function we need to evaluate.
    # ...: other arguments for `func`.
    # trace: Logical. True means print out the iteration process of the bisection algorithm.
    #
    # Output:
    # Nm: the root of equation `func` = 1 - beta.
    # fm: the return of the function `func` at Nm.
    while((Nu - Nl) > 1){
        Nm <- (Nl + Nu) / 2
        fm <- func(Nm, ...)
        if(trace){
            cat('(Nl, Nm, Nu) =', '(', Nl, ',', Nm, ',', Nu, ')\n')
            cat('(Pl, Pm, Pu) =', '(', Pl, ',', fm$P, ',', Pu, ')\n\n')   
        }
        if(abs(fm$P - 1 + beta) < eps) break
        if(((fm$P - 1 + beta) * (Pl - 1 + beta))< 0){
            Nu <- Nm
            Pu <- fm$P
        }else{
            Nl <- Nm
            Pl <- fm$P
        }
    }
    list(Nm=Nm, fm=fm)
}

##########################################################################################
# L-BFGS-B bisection algorithm to find optimal design configuration.
# This function implements Algorithm 1 in 
# "Multiple assessments of non-inferiority trials with ordinal endpoints" by Xu et al.
# Notice that `qmvnorm` function in R package `mvtnorm` uses a random algorithm to compute
# quantiles of a multivariate normal distribution. To reduce the randomness, we calculate
# critical value `nrepeat` times and use their average as the final critical value.
##########################################################################################
Search <- function(v, Lambda1, Lambda2, lower=c(0.1, 0.1), upper=c(1, 0.8), mu=NULL, sigma=NULL, 
                   tau=NULL, Delta1=NULL, Delta2=NULL, alpha=0.025, beta=0.2, nrepeat=1){
    # Input:
    # v: An (m+2) x 1 vector in which the element is sample size times the corresponding
    #   variance of the estimates of the treatment mean.
    # Lambda1: A (m+1) x 1 vector. Lambda1[i] = mu_{E_i} + Delta2 - mu_R for i = 1, ...., m, 
    #   where mu_j is treatment mean of treatment j for j = E_1, ..., E_m and R, and Delta2
    #   is NI margin.
    # Lambda2: = mu_R - mu_P - Delta1, where Delta1 is assay sensitivity margin.
    # lower: A 2 x 1 vector. The first element is the lower bound of min( n_{E_l} ):n_R for 
    #   l = 1, ..., m. The second element is the lower bound of n_P : n_R.
    # upper: A 2 x 1 vector. The first element is the upper bound of max( n_{E_l} ):n_R for 
    #   l = 1, ..., m. The second element is the upper bound of n_P : n_R.
    # mu: An (m+2) x 1 vector of means of treatment E_1, ..., E_m, P and R, where m is number
    #   of experimental treatments.
    # sigma: An (m+2) x 1 vector of standard deviations of treatment E_1, ..., E_m, P and R,
    #   where m is number of experimental treatments.
    # tau: An (K+1) x 1 vector of thresholds in which the first and last elements are 
    #   infinity and -infinity, where K is number of category.
    # Delta1: Assay sensitivity margin.
    # Delta2: NI margin.
    # alpha: significant level.
    # beta: Type 2 error.
    # nrepeat: An integer of how many times to run `qmvnorm` function to compute quantile
    #   of the multivarite normal distribution. 
    #
    # Notice:
    # (mu, Delta1, Delta2) are only used to compute Lambda1 and Lambda2 if they are nulls.
    # (mu, sigma, tau) are only used to compute v if v is null.
    #
    # Output: A list contains the following elements
    # N: the minimum total sample size.
    # P: the power at N, which should close to 1 - beta.
    # cE: A vector of m x 1 vector of optimal allocation proportion of experimental treatment
    #   compared to treatment R.
    # cP: optimal allocation proportion compared to treatment R.
    # cEprime: optimal allocation proportion of experimental treatment such that 
    #   cE[i] = cEprime[i] * sqrt(v[i]) / min(sqrt(v[i])).
    m <- ifelse(is.null(mu), length(v) - 2, length(mu) - 2)
    # compute Lambda1 and Lambda2 for computing power
    if(is.null(Lambda1)) Lambda1 <- mu[1:m] - mu[m+2] + Delta2
    if(is.null(Lambda2)) Lambda2 <- mu[m+2] - mu[m+1] - Delta1
    # find maximum of c_E^prime
    if(is.null(v)) v <- Delta(mu, sigma, tau)$v
    cEprime_max <- min(sqrt(v[1:m])) / max(sqrt(v[1:m]))
    # Obtain initial value of N with max(cE_l) = 0.8 and cP = 0.2.
    cEprime0 <- 0.8 * cEprime_max; cP0 <- 0.2 
    cE0 <- cEprime0 * sqrt(v[1:m]) / min(sqrt(v[1:m]))
    #   compute critical value given ratio
    ratio <- c(cE0, cP0, 1)
    Sig <- Sigma(ratio, v)
    cv <- CritVal(Sig, alpha, nrepeat)$cv_mean
    #   use bisection method to find N s.t. P = 1 - beta. 
    #       compute lower and upper bound of N for bisection algorithm
    Nl <- 100; Nu <- 1000 # initial lower and upper of N
    fl <- Power(Nl, ratio, Lambda1, Lambda2, cv, v, Sig)
    fu <- Power(Nu, ratio, Lambda1, Lambda2, cv, v, Sig)
    # backward and forward search make sure Pl - (1 - beta) < 0 and Pu - (1 - beta) > 0
    Nf <- BF(beta, Nl, Nu, fl, fu, Power, ratio, Lambda1, Lambda2, cv, v, Sig)
    Nl <- Nf$Nl; Nu <- Nf$Nu; Pl <- Nf$fl$P; Pu <- Nf$fu$P
    if(Pu - 1 + beta < 1e-2){ # Pu close to beta (diff is less than 0.01)
        N0 <- Nu; P0 <- Pu
    }else{ #  bisection algorithm
        root <- Bisec(1e-2, beta, Nl, Nu, Pl, Pu, Power, ratio, Lambda1, Lambda2, cv, v, Sig)
        N0 <- root$Nm; P0 <- root$fm$P
    }
    #cat('Initial sample size and corr power =', N0, ',', P0, '\n')
    Nu <- N0
    Nl <- Nu / 2
    # theoretically, Pu should not less than 1 - beta because our initial N
    init <- c(cEprime0, cP0)
    fu <- OptAllo(Nu, v, Lambda1, Lambda2, init, lower, upper, alpha, nrepeat)
    #cat('With the initial sample size and optimal allocation, power =', fu$P, '\n')
    fl <- OptAllo(Nl, v, Lambda1, Lambda2, init, lower, upper, alpha, nrepeat)
    Nf <- BF(beta, Nl, Nu, fl, fu, OptAllo, v, Lambda1, Lambda2, init, lower, upper, alpha, nrepeat)
    Nl <- Nf$Nl; Nu <- Nf$Nu; Pl <- Nf$fl$P; Pu <- Nf$fu$P
    #   bisection algorithm
    if(Pu - 1 + beta < 1e-3){ # Pu close to beta (diff is less than 0.01)
        Opt <- Nf$fu
        Opt$N <- Nf$Nu
    }else{ #  bisection algorithm
        root <- Bisec(1e-3, beta, Nl, Nu, Pl, Pu, OptAllo, v, Lambda1, Lambda2, init, lower, upper, alpha, nrepeat)
        Opt <- root$fm
        Opt$N <- root$Nm
    }
    Opt$N0 <- N0
    Opt$P0 <- P0
    return(Opt)
}

