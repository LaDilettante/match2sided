#' Run the MCMC estimation
#'
#' @export
match2sided <- function(iter, eps_alpha, eps_beta, frac_beta, frac_opp,
                 ww, xx, choice, opp) {
  n_i <- nrow(xx)
  n_j <- nrow(ww)
  p_i <- ncol(xx)
  p_j <- ncol(ww)

  prior <- list(alpha = list(mu = rnorm(p_j),
                             Tau = solve(diag(abs(rnorm(p_j))))),
                mu_beta = list(mu = rnorm(p_i),
                               Tau = solve(diag(abs(rnorm(p_i))))),
                Tau_beta = list(nu = p_i + 2,
                                Sinv = solve(diag(abs(rnorm(p_i))))))

  # ---- Starting values ----
  alpha <- rep(0, p_j)            # worker preferences
  exp_WA <- exp(ww %*% alpha) # linear predictors of
  pA_den <- opp %*% exp_WA # vector of denominators in p(A | O, alpha)

  # beta starting values (from 1-sided logit estimates)
  beta  <- matrix(0, p_i, n_j)   # employer preferences
  for(j in 1:n_j) {
    y <- as.numeric(opp[, j])
    mod <- glm(y ~ . - 1, family=binomial,
               data=as.data.frame(xx) )
    beta[, j] <- mod$coef
  }
  XB <- xx %*% beta # worker side linear predictors (big matrix), n_i x n_j, same as opp

  # mu, Sigma starting values
  mu_beta <- prior$mu_beta$mu
  Tau_beta <- prior$Tau_beta$Sinv

  # ---- Pre compute ----
  tmp <- as.matrix(ww[choice, ])
  wa <- apply(tmp, 2, sum) # sum of characteristics of accepted jobs; used in alpha update

  bmat <- eps_beta * matrix(1, p_i, n_j)

  # ---- Initialize storage ----
  acceptance_rate <- rep(0, 3)                  # Metropolis acceptance rates
  asave <- matrix(NA, iter, p_j)   # saved alphas
  bsave <- array(NA, dim = c(iter, p_i, n_j)) # saved betas
  mu_betasave <- matrix(NA, iter, p_i)
  Tau_betasave <- array(NA, dim = c(iter, p_i, p_i))
  logpost <- numeric(iter) # log posterior density (unnormalized)

  # ---- Loop ----
  for (i in 1:iter) {
    # Update opp
    num_new_offers <- floor(frac_opp * n_j)
    new <- replicate(n_i,
                     sample(2:n_j, size=floor(frac_opp * n_j), replace = FALSE))
    new <- c(new) # Flatten 1-column matrix into a vector
    ind <- cbind(rep(1:n_i, each = num_new_offers), new)

    my_logmh_opp <- logmh_opp(opp, new, alpha, beta, ww, xx)
    ok_opp <- log(runif(n_i)) <= my_logmh_opp

    # Don't change an offer for an accepted job
    accepted_jobs_sampled <- colSums(matrix(new == rep(choice, each = num_new_offers),
                                            nrow = num_new_offers)) > 0
    ok_opp[accepted_jobs_sampled] <- F  # don't change an offer for an accepted job
    if (any(ok_opp)) {
      opp[ind][rep(ok_opp, each = num_new_offers)] <- !(opp[ind][rep(ok_opp, each=num_new_offers)]) # Update the opportunity set
      acceptance_rate[1] <- acceptance_rate[1] + 1
    }

    # Update alpha
    deviation <- eps_alpha * runif(p_j, min=-1, max=1) # Symmetric proposal
    alphastar <- alpha + deviation

    my_logmh_alpha <- logmh_alpha(alpha, alphastar, ww, opp, wa, prior)
    ok_alpha <- ifelse(log(runif(1)) <= my_logmh_alpha, T, F)
    if (ok_alpha) {
      alpha <- alphastar
      acceptance_rate[2] <- acceptance_rate[2] + 1
    }

    # Update beta
    whichones <- sample(c(0, 1), size = p_i * n_j, replace = TRUE,
                        prob = c(1 - frac_beta, frac_beta))
    # Sample betastar from a [-eps_beta, eps_beta] box around beta
    rmat <- matrix(runif(p_i * n_j, min=-1, max=1) * whichones,
                   nrow = p_i, ncol = n_j)
    deviation <- bmat * rmat
    betastar <- beta + deviation
    my_logmh_beta <- logmh_beta(beta, betastar, xx, opp, mu_beta, Tau_beta)
    ok_beta <- ifelse(log(runif(1)) <= my_logmh_beta, T, F)
    if (ok_beta) {
      beta <- betastar
      acceptance_rate[3] <- acceptance_rate[3] + 1
    }

    # Update mu (multivariate normal)
    mu_beta_posterior <- cond_mu_beta(Tau_beta, prior, beta)
    mu_beta <- mvtnorm::rmvnorm(1, mu_beta_posterior$m,
                                   mu_beta_posterior$V)
    mu_beta <- c(mu_beta) # Flatten into a vector

    # Update Tau (Wishart)
    Tau_beta_posterior <- cond_Tau_beta(mu_beta, prior, beta)
    Tau_beta <- MCMCpack::rwish(Tau_beta_posterior$nu,
                                Tau_beta_posterior$Sinv)

    # Store results
    logpost[i] <- joint_lpdf(opp, choice,
                             alpha, beta, mu_beta, Tau_beta, prior,
                             ww, xx)
    asave[i, ] <- alpha
    bsave[i, ,] <- beta # vectorize
    mu_betasave[i, ] <- mu_beta
    Tau_betasave[i, ,] <- Tau_beta

    # Increment
    if (i %% 100 == 0) cat("Iteration", i, "done", "\n")
  }

  if (!is.null(colnames(ww))) colnames(asave) <- colnames(ww)
  if (!is.null(colnames(xx))) {
    dimnames(bsave)[[2]] <- colnames(xx)
    colnames(mu_betasave) <- colnames(xx)
    dimnames(Tau_betasave)[[2]] <- colnames(xx)
    dimnames(Tau_betasave)[[3]] <- rownames(xx)
  }
  if (!is.null(rownames(ww))) {
    dimnames(bsave)[[3]] <- rownames(ww)
  }

  return(list(alpha = asave, beta = bsave,
              mu_beta = mu_betasave, Tau_beta = Tau_betasave,
              lp = logpost,
              acceptance_rate = acceptance_rate / iter,
              mcmc_settings = list(iter = iter,
                                  eps_alpha = eps_alpha,
                                  eps_beta = eps_beta,
                                  frac_beta = frac_beta,
                                  frac_opp = frac_opp)))
}
