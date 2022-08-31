##########################################################################################
# Loss function for estimating parameters of one treatment.
# This function gives: (1) val, loss value; (2) grad, gradient; (3) hessian, expectationof
# the Hessian matrix.
##########################################################################################
loss <- function(n, tau, mu, sigma){
    K <- length(n) # number of classes
    r.u <- tau[-1] - mu
    r.l <- tau[-(K+1)] - mu
    cdf <- pnorm(r.u/sigma) - pnorm(r.l/sigma)
    # loss value, i.e., negative log likelihood function
    val <- -sum(n * log(cdf))
    # compute gradient w.r.t. mu
    pdf.u <- dnorm(r.u/sigma)
    pdf.l <- dnorm(r.l/sigma)
    pdf <- pdf.u - pdf.l
    grad_mu <- sum(n * pdf / cdf) / sigma
    # compute gradient w.r.t. sigma
    r.u_x_pdf.u <- r.u * pdf.u # r.u x pdf.u
    r.u_x_pdf.u[K] <- 0 # inf x phi(inf) = 0
    r.l_x_pdf.l <- r.l * pdf.l # r.l x pdf.l
    r.l_x_pdf.l[1] <- 0 # -inf x phi(-inf) = 0
    r_x_pdf <- r.u_x_pdf.u - r.l_x_pdf.l
    grad_sigma <- sum(n * r_x_pdf / cdf) / sigma^2
    grad <- c(grad_mu, grad_sigma)
    # compute expected Hessian w.r.t. mu
    hessian_mu <- sum(pdf^2 / cdf) * sum(n) / sigma^2
    # compute expected Hessian w.r.t. sigma
    hessian_sigma <- sum(r_x_pdf^2 / cdf) * sum(n) / sigma^4
    # compute expected Hessian w.r.t. mu sigma.
    hessian_mu_sigma <- sum(pdf * r_x_pdf / cdf) * sum(n) / sigma^3
    hessian <- matrix(0, 2, 2)
    hessian[1, 1] <- hessian_mu
    hessian[1, 2] <- hessian_mu_sigma
    hessian[2, 1] <- hessian_mu_sigma
    hessian[2, 2] <- hessian_sigma
    list(val=val, grad=grad, hessian=hessian)
}

##########################################################################################
# Fisher scoring algorithm for estimating parameters of one treatment.
##########################################################################################
FS4one <- function(n, tau, mu, sigma, itmax=1000, eps=1e-6){
    # Input:
    # n: A vector of K x 1, where K is number of category. The i-th element in n is sample
    #   size of category i.
    # tau: An (K+1) x 1 vector of thresholds in which the first and last elements are 
    #   infinity and -infinity, where K is number of category.
    # mu: mean of one treatment.
    # sigma: standard deviation of one treatment.
    # itmax: The maximum number of iterations allowed.
    # eps: The tolerance at which the difference in the objective values is considered close 
    #   enough to 0 to declare convergence.
    #
    # Output:
    # mu: The estimated treatment mean.
    # sigma: The estimated treatment standard deviation.
    # conv: Logical; TRUE if convergence is obtained.
    # objs: The vector of loss at each iteration.
    it <- 1; convergence <- FALSE
    info <- loss(n, tau, mu, sigma)
    obj <- info$val
    objs <- info$val
    while(it <= itmax && !convergence){
        hessian.inverse <- solve(info$hessian)
        d <- drop(hessian.inverse %*% info$grad)
        mu <- mu - d[1]
        sigma <- sigma - d[2]
        if(sigma <= 0) break
        info <- loss(n, tau, mu, sigma)
        objs <- c(objs, info$val)
        #if(is.infinite(info$val)) break
        if(abs(info$val - obj) <= eps){
            convergence <- TRUE
        }else{
            obj <- info$val
            it <- it + 1
        }
    }
    list(mu=mu, sigma=sigma, conv=convergence, objs=objs)
}

##########################################################################################
# Fisher scoring algorithm for estimating all the parameters.
##########################################################################################
FS <- function(tr, y, mu0=NULL, sigma0=NULL, itmax=1000, eps=1e-6){
    # Input:
    # tr: An n x 1 vector of treatment label, where n is total samole size. The i-th element 
    #   in tr equals j if the i-th sample belongs to treatment E_j for j = 1, ..., m, and 
    #   equals m+1 if the i-th sample belongs to treatment P; Otherwise, equals m+2, where m 
    #   is number of experimental treatments.
    # y: An n x 1 vector of categories. The i-th element in y equals k if its category is k
    # for k = 1, ..., K with K being total number of categories.
    # mu0: An (m+2) x 1 vector of initial means of treatment E_1, ..., E_m, P and R, where m 
    #   is number of experimental treatments.
    # sigma0: An (m+2) x 1 vector of initial standard deviations of treatment E_1, ..., E_m, P and R,
    #   where m is number of experimental treatments.
    # itmax: The maximum number of iterations allowed.
    # eps: The tolerance at which the difference in the objective values is considered close 
    #   enough to 0 to declare convergence.
    # Output:
    # tau: The (K+1) x 1 vector of thresholds including -infinity and infinity.
    # mu: The (m+2) x 1 vector of estimated treatment means.
    # sigma: The (m+2) x 1 vector of estimated treatment standard deviation.
    # conv: A vector indicating if algorithm is convergent, along the treatment.
    # objs: A list of loss at each iteration, along the treatment.
    
    # compute tau based on observations in treatment R
    freq <- table(y[tr==(max(tr))])
    cdf <- cumsum(freq) / sum(freq)
    cdf <- c(0, cdf) # include -inf
    tau <- qnorm(cdf)
    # compute sample size of treatment i with label k, i = E_1, ..., E_m, P, R.
    n <- table(tr, y)
    # compute mu and sigma using FS algorithm
    nt <- dim(n)[1] # number of treatments
    if(is.null(mu0)) mu0 <- rep(0, nt)
    if(is.null(sigma0)) sigma0 <- rep(1, nt)
    mu <- rep(NA, nt)
    sigma <- rep(NA, nt)
    conv <- rep(NA, nt)
    objs <- list()
    for(i in 1:nt){
        fit <- FS4one(n[i, ], tau, mu0[i], sigma0[i], itmax, eps)
        mu[i] <- fit$mu
        sigma[i] <- fit$sigma
        conv[i] <- fit$conv
        objs[[i]] <- fit$objs
    }
    list(tau=tau, mu=mu, sigma=sigma, conv=conv, objs=objs)
}