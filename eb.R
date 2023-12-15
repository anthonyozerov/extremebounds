# generate covariates
X <- rnorm(1000 * 10) |> matrix(ncol = 10)
X_df <- as.data.frame(X)
colnames(X_df) <- paste0("X", 1:10)

# generate treatment indicator
alpha <- X %*% seq(2.5, -2.5, -5/9) + rnorm(1000)
Z <- rbinom(1000, 1, 1/(1+exp(-alpha)))

tau <- 10 # the true ATE
Y <- X %*% seq(-2.5, 2.5, 5 / 9) + tau * Z + rnorm(1000)

# define several linear outcome and propensity score models.
f_mu0s <- c()
f_mu1s <- c()
f_es <- c()
# for every combination of variables in X containing at least 6 variables
# create a function that contains a linear term for each variable
for (i in 0:(2**10 - 1)) {
  if (sum(as.integer(intToBits(i))) >= 7) {
    X_terms <- paste0(colnames(X_df)[as.logical(intToBits(i))], collapse = "+")
    f_mu <- paste0("Y ~ Z+", X_terms)
    f_mu0s <- c(f_mu0s, f_mu)
    f_mu1s <- c(f_mu1s, f_mu)
    f_e <- paste0("Z ~ ", X_terms)
    f_es <- c(f_es, f_e)
  }
}

# E is a vector of estimators
E <- c()

# add IPW to E
truncs <- c(0.01, 0.05, 0.1, 0.15)
ipw_ht <- function(Z, Y, pscores, trunc) {
  # truncate the propensity scores
  pscores <- pmin(pscores, 1 - trunc)
  pscores <- pmax(pscores, trunc)
  # calculate the weights
  weights <- Z / pscores + (1 - Z) / (1 - pscores)
  # calculate the weighted average
  return(sum(weights * Y) / sum(weights))
}
ipw_ht_t <- function(trunc) {
  f <- function(Z, Y, pscores) {
    ipw_ht(Z, Y, pscores, trunc)
  }
  return(f)
}

ipw_hajek <- function(Z, Y, pscores, trunc) {
  # truncate the propensity scores
  pscores <- pmin(pscores, 1 - trunc)
  pscores <- pmax(pscores, trunc)

  t1 <- Z / pscores
  t2 <- (1 - Z) / (1 - pscores)

  return(sum(t1 * Y) / sum(t1) - sum(t2 * Y) / sum(t2))
}
ipw_hajek_t <- function(trunc) {
  f <- function(Z, Y, pscores) {
    ipw_hajek(Z, Y, pscores, trunc)
  }
  return(f)
}

for (t in truncs) {
  E <- c(E, ipw_ht_t(t), ipw_hajek_t(t))
}

# PS stratification estimator

# number of strata
Ks <- c(2, 5, 10, 20)
pscore_strata <- function(pscore, K) {
  strata <- cut(pscore, breaks = K, labels = 1:K)
  return(strata)
}

# strata is stratum assignments
strat <- function(Z, Y, s) {
  K <- nlevels(s)

  # get the number of treated and untreated in each stratum
  n <- length(Z)
  n_1 <- sapply(1:K, function(k) sum(Z[s == k]))
  n_0 <- sapply(1:K, function(k) sum(1 - Z[s == k]))


  # get the mean of Y in each stratum
  mu_1 <- sapply(1:K, function(k) mean(Y[s == k & Z == 1]))
  mu_0 <- sapply(1:K, function(k) mean(Y[s == k & Z == 0]))

  # get the weighted average of the difference in means
  return(sum((mu_1 - mu_0) * (n_1 + n_0) / n))
}

# measure covariate balance in each stratification using absolute difference
covariate_balance <- function(X, Z, strata) {
  # within each stratum
  # get the mean of each covariate for treated and untreated
  # get the absolute difference between the means
  # return the average absolute difference
  K <- length(unique(strata))
  covs <- colnames(X)

  # get the number of treated and untreated in each stratum
  n <- length(Z)
  n_1 <- sapply(1:K, function(k) sum(Z[strata == k]))
  n_0 <- sapply(1:K, function(k) sum(1 - Z[strata == k]))

  # get the mean of each covariate in each stratum
  mu_1 <- sapply(1:K, function(k) colMeans(X[strata == k & Z == 1, ]))
  mu_0 <- sapply(1:K, function(k) colMeans(X[strata == k & Z == 0, ]))

  # get the sd of each covariate in each stratum
  sd_1 <- sapply(1:K, function(k) apply(X[strata == k & Z == 1, ], 2, sd))
  sd_0 <- sapply(1:K, function(k) apply(X[strata == k & Z == 0, ], 2, sd))

  # get the absolute difference between the means
  abs_diff <- abs(mu_1 - mu_0) / sqrt((sd_1^2 + sd_0^2) / 2)

  # weighted average across strata
  abs_diff <- rowSums(abs_diff * c(n_1, n_0)) / n

  return(mean(abs_diff))
}

# define strata-making functions
S <- lapply(Ks, function(K) function(pscore) pscore_strata(pscore, K))

# function to check if a stratification contains one unit from each treatment group
# in each stratum
check_stratification <- function(Z, strata) {
  K <- length(unique(strata))
  # get the number of treated and untreated in each stratum
  n <- length(Z)
  n_1 <- sapply(1:K, function(k) sum(Z[strata == k]))
  n_0 <- sapply(1:K, function(k) sum(1 - Z[strata == k]))

  # return true if there is at least one treated and one untreated in each stratum
  return(all(n_1 > 0) && all(n_0 > 0))
}

produce_estimates <- function(X_df, Z, Y) {
  print("fitting pscore models")
  # get vector of fitted pscore models
  pscore_models <- lapply(f_es, function(f) glm(f, data = X_df, family = binomial))

  print("getting fitted pscores and screening")
  # get vector of fitted pscores
  pscores <- lapply(pscore_models, function(m) predict(m, X_df, type = "response"))
  # get vector of aics
  aics <- sapply(pscore_models, function(m) AIC(m))
  # remove pscores in bottom half of aics
  good_pscores <- aics > median(aics)

  print("computing strata")
  # get vector of stratum assignments using pscores
  strata <- lapply(S, function(s) lapply(pscores, function(p) s(p)))
  good_strata <- unlist(lapply(strata, function(s) good_pscores))
  strata <- unlist(strata, recursive = FALSE)

  print("measuring covariate balance")
  # get vector of covariate balance
  # balance checking is a fairly expensive procedure, so we will skip it for now
  # mean_diff <- sapply(strata, function(s) covariate_balance(X_df, Z, s))
  # discard strata with poor covariate balance or NA
  # good_strata2 <- (mean_diff > quantile(mean_diff, 0.25, na.rm = TRUE)) & !is.na(mean_diff)
  # good_strata <- good_strata & good_strata2
  # check that the strata each have a treatment and control unit
  good_strata3 <- sapply(strata, function(s) check_stratification(Z, s))
  good_strata <- good_strata & good_strata3

  print("computing ipw estimates")
  # get vector of estimates using pscores in ipw_hajek and ipw_hajek_t
  tauhats_ipw <- lapply(E, function(e) lapply(pscores, function(p) e(Z, Y, p)))
  good_tauhats_ipw <- unlist(lapply(E, function(e) good_pscores))
  # flatten nested list into a vector
  tauhats_ipw <- unlist(tauhats_ipw)

  print("computing stratified estimates")
  tauhats_strat <- sapply(strata, function(s) strat(Z, Y, s))
  good_tauhats_strat <- good_strata
  tauhats <- c(tauhats_ipw, tauhats_strat)
  good_tauhats <- c(good_tauhats_ipw, good_tauhats_strat)

  return(list(tauhats = tauhats, good_tauhats = good_tauhats))
}

result <- produce_estimates(X_df, Z, Y)
good_tauhats <- result$good_tauhats
tauhats <- result$tauhats
tau_min <- min(tauhats[good_tauhats])
tau_max <- max(tauhats[good_tauhats])

I1 <- c(tau_min, tau_max)

B <- 100
tauhat_mat <- matrix(tauhats, nrow = 1, ncol = length(tauhats))
good_tauhat_mat <- matrix(good_tauhats, nrow = 1, ncol = length(good_tauhats))
for (b in 1:B) {
  print(b)
  n <- length(Z)
  samples <- sample(1:n, n, replace = TRUE)
  # sample with replacement from the data
  X_dfb <- X_df[samples, ]
  Zb <- Z[samples]
  Yb <- Y[samples]

  result_b <- produce_estimates(X_dfb, Zb, Yb)
  tauhat_mat <- rbind(tauhat_mat, result_b$tauhats)
  good_tauhat_mat <- rbind(good_tauhat_mat, result_b$good_tauhats)
}

bootstrap_param <- function(est, ests_b, alpha = 0.05) {
  sigmahat <- sd(ests_b)
  return(c(est + sigmahat * qnorm(alpha / 2), est + sigmahat * qnorm(1 - alpha / 2)))
}

bootstrap_basic <- function(est, ests_b, alpha = 0.05) {
  est_lower <- unname(quantile(ests_b, alpha / 2))
  est_upper <- unname(quantile(ests_b, 1 - alpha / 2))
  return(c(2 * est - est_upper, 2 * est - est_lower))
  return(quantile(ests_b, c(alpha / 2, 1 - alpha / 2)))
}

bootstrap_percentile <- function(est, ests_b, alpha = 0.05) {
  return(unname(quantile(ests_b, c(alpha / 2, 1 - alpha / 2))))
}

tau_min_b <- sapply(2:nrow(tauhat_mat), function(i) min(tauhat_mat[i, good_tauhat_mat[i, ]]))
tau_max_b <- sapply(2:nrow(tauhat_mat), function(i) max(tauhat_mat[i, good_tauhat_mat[i, ]]))

# compute different versions of the outer interval depending on the bootstrap method
I3_1 <- c(bootstrap_param(tau_min, tau_min_b)[1], bootstrap_param(tau_max, tau_max_b)[2])
I3_2 <- c(bootstrap_basic(tau_min, tau_min_b)[1], bootstrap_basic(tau_max, tau_max_b)[2])
I3_3 <- c(bootstrap_percentile(tau_min, tau_min_b)[1], bootstrap_percentile(tau_max, tau_max_b)[2])

# now we need to compute the middle interval.

# find which estimators were "good" at least 20 times.
# This is a heuristic strategy.
middle_tauhat_mat <- good_tauhat_mat[, (colSums(good_tauhat_mat) >= 20)]
# for each estimator, get a confidence interval using each bootstrap method
lower_bounds <- c()
upper_bounds <- c()
for (i in 1:ncol(middle_tauhat_mat)) {
  tauhat <- tauhat_mat[1, i]
  tauhat_b <- tauhat_mat[-1, i][good_tauhat_mat[-1, i]]

  int1 <- bootstrap_param(tauhat, tauhat_b)
  int2 <- bootstrap_basic(tauhat, tauhat_b)
  int3 <- bootstrap_percentile(tauhat, tauhat_b)
  # if tauhat is not good, the parametric and basic intervals don't make sense
  if (!good_tauhat_mat[1, i]) {
    lower_bounds <- c(lower_bounds, int3[1])
    upper_bounds <- c(upper_bounds, int3[2])
  } else {
    lower_bounds <- c(lower_bounds, int1[1], int2[1], int3[1])
    upper_bounds <- c(upper_bounds, int1[2], int2[2], int3[2])
  }
}
# together, lower_bounds and upper_bounds are our bag of confidence intervals
# we can now make the middle interval.
I2 <- c(min(lower_bounds), max(upper_bounds))

# print the intervals nicely
print("model interval:")
print(round(I1, 2))
print("CI cover:")
print(round(I2, 2))
print("Extreme CI:")
print(round(I3_1, 2))
print(round(I3_2, 2))
print(round(I3_3, 2))