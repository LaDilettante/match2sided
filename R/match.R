#' Two sided matching model
#'
#' @param opp A matrix of opportunity
#' @param choice A vector of observed choice
#' @param alpha Multinomial preference
#' @param beta Binomial preference
#' @param mu_beta Hyperparameter of beta
#' @param sigma_beta Hyperparameter of beta
#' @param prior A list of prior parameters
#' @name match2sided
NULL

# ---- MCMC ----

#' Calculate the log joint pdf, useful for testing
#'
#' @rdname match2sided
#' @importFrom mvtnorm dmvnorm
joint_lpdf <- function(opp, choice, alpha, beta,
                       mu_beta, Tau_beta,
                       prior, ww, xx) {
  n_i <- length(choice)
  w <- as.matrix(ww[choice, ])
  logp_A <- sum(w %*% alpha) - sum(log(opp %*% exp(ww %*% alpha)))

  XB <- xx %*% beta
  logp_O <- sum(opp * XB) - sum(log1p(exp(XB)))

  return(logp_A + logp_O +
    dmvnorm(alpha, prior$alpha$mu, solve(prior$alpha$Tau), log = TRUE) +
    sum(dmvnorm(t(beta), mu_beta, solve(Tau_beta), log = TRUE)) +
    dmvnorm(mu_beta, prior$mu_beta$mu, solve(prior$mu_beta$Tau), log = TRUE) +
    log(MCMCpack::dwish(Tau_beta, prior$Tau_beta$nu, prior$Tau_beta$S)))
}

f_pA_den <- function(opp, ww, alpha) {
  exp_WA <- exp(ww %*% alpha)
  return(c(opp %*% exp_WA))
}

#' log_mh opp_set
#' @param indnew Indexes of jobs to be flipped
logmh_opp <- function(opp, new, alpha, beta, ww, xx) {
  # browser()
  n_i <- nrow(xx)
  num_new_offers <- length(new) / n_i
  n_j <- nrow(ww)

  indnew <- cbind(rep(1:n_i, each = num_new_offers), new)
  oo <- opp[indnew] # current offer status of the jobs under consideration
  plusminus <- ifelse(oo, -1, 1)

  pA_den <- f_pA_den(opp, ww, alpha)

  # Part of MH ratio from P(A|O,alpha)
  pA_denstar <- pA_den +
    colSums(matrix(exp(ww %*% alpha)[new] * plusminus, ncol = n_i))
  # Part of MH ratio from P(O|beta)  (logistic model)
  xb <- (xx %*% beta)[indnew]
  exp_xb <- colSums(matrix(exp(plusminus * xb), ncol = n_i))

  return(log(pA_den) - log(pA_denstar) + log(exp_xb))
}

#' log_mh alpha
logmh_alpha <- function(alpha, alphastar, ww, opp, wa, prior) {

  # Calculate likelihood ratio
  exp_WA <- exp(ww %*% alpha)
  pA_den <- opp %*% exp_WA
  exp_WA_star <-  exp(ww %*% alphastar)
  pA_denstar <- opp %*% exp_WA_star

  logmh_alpha <- sum(wa * (alphastar - alpha)) + sum(log(pA_den) - log(pA_denstar)) +
    dmvnorm(alphastar, prior$alpha$mu, solve(prior$alpha$Tau), log = TRUE) -
    dmvnorm(alpha, prior$alpha$mu, solve(prior$alpha$Tau), log = TRUE)
  return(logmh_alpha)
}

#' log_mh beta
logmh_beta <- function(beta, betastar, xx, opp, mu_beta, Tau_beta) {
  XB <- xx %*% beta
  XB_star <- xx %*% betastar
  # Calculate likelihood ratio (using logistic structure)(don't count unemp.
  lrat <- sum((opp * (XB_star - XB)))    # the `canonical' part (unemp. cancels)
  logmh_beta <- lrat + sum(log(1 + exp(XB)) - log(1 + exp(XB_star))) +
    sum(dmvnorm(t(betastar), mu_beta, solve(Tau_beta), log = TRUE)) -
    sum(dmvnorm(t(beta), mu_beta, solve(Tau_beta), log = TRUE))
  return(logmh_beta)
}

cond_mu_beta <- function(Tau_beta, prior, beta) {
  # browser()
  n_j <- ncol(beta)
  Vinv <- prior$mu_beta$Tau + n_j * Tau_beta
  V <- solve(Vinv)
  m <- c(V %*% (prior$mu_beta$Tau %*% prior$mu_beta$mu +
                n_j * Tau_beta %*% apply(beta, 1, mean)))
  return(list(m=m, V=V))
}

cond_Tau_beta <- function(mu_beta, prior, beta) {
  n_j <- ncol(beta)
  nu <- prior$Tau_beta$nu + n_j
  S <- solve(prior$Tau_beta$Sinv) + (beta - mu_beta) %*% t(beta - mu_beta)
  return(list(nu = nu, Sinv = solve(S)))
}
