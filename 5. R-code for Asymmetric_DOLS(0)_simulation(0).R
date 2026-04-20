#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MASTER SCRIPT: Asymmetric-DOLS Monte Carlo (1000 Reps)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++
setwd("C:/R/fmols")
library(cointReg)

# --- Simulation Parameters ---
T_val <- 200     # Sample size
reps <- 1000     # Number of repetitions
true_beta <- 2.5
bandwidth <- 8   # Fixed bandwidth for comparison

# Matrices to store t-statistics and slopes for 4 scenarios
sim_plus          <- matrix(NA, nrow = reps, ncol = 4)
sim_minus         <- matrix(NA, nrow = reps, ncol = 4)
sim_tstats_plus   <- matrix(NA, nrow = reps, ncol = 4)
sim_tstats_minus  <- matrix(NA, nrow = reps, ncol = 4)
colnames(sim_tstats_plus)  <- c("Strong", "Weak", "Very_Weak", "None")
colnames(sim_tstats_minus) <- c("Strong", "Weak", "Very_Weak", "None")

n_leads <- 2
n_lags  <- 2

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: Refined Asymmetric DLOS Function
#++++++++++++++++++++++++++++++++++++++++++++++++
asymmetric_dols <- function(y, x, n_leads, n_lags) {
  # Ensure inputs are matrices
  y <- as.matrix(y)
  x <- as.matrix(x)
  T <- nrow(x)

  # 1. Asymmetric Decomposition (Shin 2014)
  dx <- diff(x)
  dx_pos <- ifelse(dx > 0, dx, 0)
  dx_neg <- ifelse(dx < 0, dx, 0)

  # Calculate partial sums (x+ and x-)
  x_pos <- as.matrix(cumsum(c(0, dx_pos)))
  x_neg <- as.matrix(cumsum(c(0, dx_neg)))

  # 2. Alignment and Trimming
  # Valid indices accounting for differencing and leads/lags
  valid_idx <- (n_lags + 2):(T - n_leads)

  y_trim <- y[valid_idx, , drop = FALSE]
  # x_trim <- x[valid_idx, , drop = FALSE]
  x_pos_trim <- x_pos[valid_idx, , drop = FALSE]
  x_neg_trim <- x_neg[valid_idx, , drop = FALSE]

  # Initialize augmented matrix with the two long-run components
  # X_aug <- cbind(x_trim, x_pos_trim, x_neg_trim)
  X_aug <- cbind(x_pos_trim, x_neg_trim)

  # 3. Add Leads and Lags for BOTH positive and negative changes
  for (j in (-n_leads):n_lags) {
    # Lead/Lag for positive differences
    dx_pos_shifted <- dx_pos[(valid_idx - 1 - j), , drop = FALSE]
    # Lead/Lag for negative differences
    dx_neg_shifted <- dx_neg[(valid_idx - 1 - j), , drop = FALSE]

    X_aug <- cbind(X_aug, dx_pos_shifted, dx_neg_shifted)
  }

  # 4. Add Intercept
  X_aug <- cbind(1, X_aug)
  XtX_inv <- tryCatch(solve(t(X_aug) %*% X_aug), error = function(e) NULL)

  if (is.null(XtX_inv)) {
    return(list(beta_plus = NA, beta_minus = NA, se_plus = NA, se_minus = NA))
  }
  # 5. OLS Estimation: (X'X)^(-1) X'Y
  beta_hat <- solve(t(X_aug) %*% X_aug) %*% t(X_aug) %*% y_trim

  # Calculate Residuals and Conventional Variance
  u_hat <- y_trim - X_aug %*% beta_hat
  df <- nrow(X_aug) - ncol(X_aug)
  sigma_sq <- sum(u_hat^2) / df

  cov_matrix <- sigma_sq * XtX_inv
  se_hat <- sqrt(diag(cov_matrix))

  # Results: 1st=Intercept, 2nd=Beta+, 3rd=Beta-
  results <- list(
    beta_plus  = beta_hat[2, 1],
    beta_minus = beta_hat[3, 1],
    se_plus    = se_hat[2],
    se_minus   = se_hat[3]
  )
  return(results)
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: The Simulation Loop
#++++++++++++++++++++++++++++++++++++++++++++++++
set.seed(42)
cat("Starting 1000 Repetitions...\n")

for (i in 1:reps) {
  # Generate X (Random Walk)
  x <- cumsum(rnorm(T_val))

  # Generate 4 Scenarios
  y1 <- 10 + true_beta * x + rnorm(T_val)                            # Strong
  y2 <- 10 + true_beta * x + arima.sim(list(ar=0.8), T_val)          # Weak
  y3 <- 10 + true_beta * x + arima.sim(list(ar=0.95), T_val)         # Very Weak
  y4 <- 10 + true_beta * x + cumsum(rnorm(T_val))                    # None

  for (j in 1:4) {
    y_current <- get(paste0("y", j))
    res <- asymmetric_dols(y_current, x, n_leads, n_lags)
    sim_plus[i, j]  <- res$beta_plus
    sim_minus[i, j] <- res$beta_minus
    # Calculate and Store T-statistics: t = (beta_hat - true_beta) / SE
    sim_tstats_plus[i, j]  <- (res$beta_plus  - true_beta) / res$se_plus
    sim_tstats_minus[i, j] <- (res$beta_minus - true_beta) / res$se_minus
  }
  if (i %% 100 == 0) cat("Repetition:", i, "\n")
}

#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: Analysis of Results
#++++++++++++++++++++++++++++++++++++++++++++++++
# Rejection Rate at 5% (Testing H0: Beta = 2.5)
# If H0 is TRUE, this should be 0.05.
power_size1 <- apply(sim_tstats_plus, 2, function(x) mean(abs(x) > 1.96))
power_size2 <- apply(sim_tstats_minus, 2, function(x) mean(abs(x) > 1.96))

cat("\nEmpirical Rejection Rates (H0: Beta = 2.5):\n")
print(power_size1)
print(power_size2)

# Summary of Slopes
cat("\nSummary of Slope Estimates:\n")
print(round(apply(sim_plus,  2, summary), 4))
print(round(apply(sim_minus, 2, summary), 4))

