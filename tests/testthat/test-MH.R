context("Metropolis-Hasting acceptance ratio")

n_i <- sample(seq(10, 20), 1) # multinomial logit side
p_i <- sample(seq(2, 3), 1)
n_j <- sample(seq(5, 10), 1) # logit side
p_j <- sample(seq(2, 3), 1)

choice <- sample(seq(1, n_j), size = n_i, replace = TRUE) # Each i picks a j
opp <- matrix(FALSE, n_i, n_j)
opp[cbind(1:n_i, choice)] <- TRUE # the choice is always in the opportunity set
# randomly flip 1/10 to TRUE
opp[sample(1:(n_i * n_j), size = n_i * n_j / 10, replace = TRUE)] <- TRUE

alpha <- rnorm(p_j)
beta <- matrix(rnorm(p_i * n_j), p_i, n_j)
mu_beta <- rnorm(p_i)
Tau_beta <- solve(diag(abs(rnorm(p_i))))
prior <- list(alpha = list(mu = rnorm(p_j), Tau = solve(diag(abs(rnorm(p_j))))),
              mu_beta = list(mu = rnorm(p_i), Tau = solve(diag(abs(rnorm(p_i))))),
              Tau_beta = list(nu = p_i + 2, Sinv = solve(diag(abs(rnorm(p_i))))))

ww <- matrix(rnorm(n_j * p_j), n_j, p_j)
tmp <- as.matrix(ww[choice, ])
wa <- apply(tmp, 2, sum)

xx <- matrix(rnorm(n_i * p_i), n_i, p_i) # Characteristics of the multinomial side
xx[, 1] <- 1 # The first column is the intercept

test_that("logmh_opp is correct", {
  # Randomly choose a proportion to be updated of the opp set
  frac_opp <- runif(1, min = 0.5, max = 1)
  num_new_offers <- floor(frac_opp * n_j)
  new <- c(replicate(n_i,
                     sample(2:n_j, size=floor(frac_opp * n_j), replace = FALSE)))
  # new <- sample(2:n_j, size=n_i, replace=T) # don't sample unemp. (1)
  indnew <- cbind(rep(1:n_i, each = num_new_offers), new)
  oppstar <- opp
  oppstar[indnew] <- !(opp[indnew])

  expected <- joint_lpdf(oppstar, choice, alpha, beta, mu_beta, Tau_beta, prior, ww, xx) -
              joint_lpdf(opp, choice, alpha, beta, mu_beta, Tau_beta, prior, ww, xx)
  # We need to sum because this is across workers
  observed <- sum(logmh_opp(opp, new, alpha, beta, ww, xx))
  expect_equal(expected, observed)
})

test_that("logmh_alpha is correct", {
  eps <- rnorm(1)
  deviation <- eps * runif(length(alpha), min=-1, max=1) # Symmetric proposal
  alphastar <- alpha + deviation

  expected <- joint_lpdf(opp, choice, alphastar, beta, mu_beta, Tau_beta, prior, ww, xx) -
              joint_lpdf(opp, choice, alpha, beta, mu_beta, Tau_beta, prior, ww, xx)
  observed <- logmh_alpha(alpha, alphastar, ww, opp, wa, prior)
  expect_equal(expected, observed)
})

test_that("logmh_beta is correct", {
  c <- ncol(beta)
  r <- nrow(beta)
  # Choosing which beta to update
  frac <- runif(1)
  whichones <- sample(c(0, 1), size = c * r, replace = TRUE,
                      prob = c(1 - frac, frac))
  # Sample betastar from a [-eps_beta, eps_beta] box around beta
  rmat <- matrix(runif(c * r, min=-1, max=1) * whichones,
                 nrow = r, ncol = c)

  eps <- rnorm(1)
  bmat <- eps * matrix(1, p_i, n_j)
  deviation <- bmat * rmat
  betastar <- beta + deviation

  expected <- joint_lpdf(opp, choice, alpha, betastar, mu_beta, Tau_beta, prior, ww, xx) -
              joint_lpdf(opp, choice, alpha, beta, mu_beta, Tau_beta, prior, ww, xx)
  observed <- logmh_beta(beta, betastar, xx, opp, mu_beta, Tau_beta)
  expect_equal(expected, observed)
})

test_that("cond_mu_beta is correct", {
  mu_beta2 <- rnorm(p_i)

  expected <- joint_lpdf(opp, choice, alpha, beta, mu_beta, Tau_beta, prior, ww, xx) -
              joint_lpdf(opp, choice, alpha, beta, mu_beta2, Tau_beta, prior, ww, xx)
  posterior <- cond_mu_beta(Tau_beta, prior, beta)
  observed <- mvtnorm::dmvnorm(mu_beta, posterior$m, posterior$V, log = TRUE) -
              mvtnorm::dmvnorm(mu_beta2, posterior$m, posterior$V, log = TRUE)
  expect_equal(expected, observed)
})

test_that("cond_Tau_beta is correct", {
  Tau_beta2 <- solve(diag(abs(rnorm(p_i))))

  expected <- joint_lpdf(opp, choice, alpha, beta, mu_beta, Tau_beta, prior, ww, xx) -
              joint_lpdf(opp, choice, alpha, beta, mu_beta, Tau_beta2, prior, ww, xx)
  posterior <- cond_Tau_beta(mu_beta, prior, beta)
  observed <- log(MCMCpack::dwish(Tau_beta, posterior$nu, posterior$Sinv)) -
              log(MCMCpack::dwish(Tau_beta2, posterior$nu, posterior$Sinv))
  expect_equal(expected, observed)
})
