##########################################################################################
# Compute delta_0, delta_1, delta_2, delta and v defined in the paper 
# "Multiple assessments of non-inferiority trials with ordinal endpoints" by Xu et al.
##########################################################################################
Delta <- function(mu, sigma, tau){
    # Input:
    # mu: An (m+2) x 1 vector of means of treatment E_1, ..., E_m, P and R, where m is number
    #   of experimental treatments.
    # sigma: An (m+2) x 1 vector of standard deviations of treatment E_1, ..., E_m, P and R.
    # tau: An (K+1) x 1 vector of thresholds including -1 x infinity and infinity, where 
    #   where K is number of category.
    #
    # Output:
    # delta0: An (m+2) x 1 vector of delta_0.
    # delta1: An (m+2) x 1 vector of delta_1.
    # delta2: An (m+2) x 1 vector of delta_2.
    # delta: An (m+2) x 1 vector of delta. 
    # v: An (m+2) x 1 vector in which the element is sample size times the corresponding
    #   variance of the estimates of the treatment mean.
    r.u <- sapply(mu, function(x) tau[-1] - x)
    r.u <- t(r.u) # (m + 2) x K, row i = tau[2:(K+1)] - mu[i]
    r.l <- sapply(mu, function(x) tau[-length(tau)] - x)
    r.l <- t(r.l) # (m + 2) x K, row i = tau[1:K] - mu[i]
    cdf <- pnorm(r.u/sigma) - pnorm(r.l/sigma)
    pdf.u <- dnorm(r.u/sigma)
    pdf.l <- dnorm(r.l/sigma)
    pdf <- pdf.u - pdf.l
    # delta0
    delta0 <- 1.0 / rowSums(pdf^2 / cdf)
    # delta1
    r.u_x_pdf.u <- r.u * pdf.u # r.u x pdf.u
    r.u_x_pdf.u[, dim(r.u_x_pdf.u)[2]] <- 0 # inf x phi(inf) = 0
    r.l_x_pdf.l <- r.l * pdf.l # r.l x pdf.l
    r.l_x_pdf.l[, 1] <- 0 # -inf x phi(-inf) = 0
    r_x_pdf <- r.u_x_pdf.u - r.l_x_pdf.l
    delta1 <- rowSums((pdf / cdf) * r_x_pdf)
    # delta2
    delta2 <- rowSums(r_x_pdf^2 / cdf)
    # delta
    delta <- delta0 + (delta0^2 * delta1^2) / (delta2 - delta0 * delta1^2)
    # v, which is sigma^2 x delta
    v = sigma^2 * delta
    list(delta0=delta0, delta1=delta1, delta2=delta2, delta=delta, v=v)
}

