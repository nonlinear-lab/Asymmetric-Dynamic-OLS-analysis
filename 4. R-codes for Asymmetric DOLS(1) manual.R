#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: Libraries and Setup
#++++++++++++++++++++++++++++++++++++++++++++++++
setwd("C:/R/fmols")
library(cointReg)
mout<- matrix(data=NA,nrow=4,ncol=2)
#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: Custom Asymmetric DOLS Function
#++++++++++++++++++++++++++++++++++++++++++++++++
asymmetric_dols <- function(y, x, n_leads, n_lags) {
  y <- as.matrix(y)
  x <- as.matrix(x)
  n <- nrow(x)
  # 1. Asymmetric Decomposition (Shin 2014)
  dx <- diff(x)
  dx_pos <- ifelse(dx > 0, dx, 0)
  dx_neg <- ifelse(dx < 0, dx, 0)

  # Calculate partial sums (x+ and x-)
  x_pos <- as.matrix(cumsum(c(0, dx_pos)))
  x_neg <- as.matrix(cumsum(c(0, dx_neg)))

  # 2. Alignment and Trimming
  # Valid indices accounting for differencing and leads/lags
  valid_idx <- (n_lags + 2):(n - n_leads)

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

  # 5. OLS Estimation: (X'X)^(-1) X'Y
  beta_hat <- solve(t(X_aug) %*% X_aug) %*% t(X_aug) %*% y_trim

  # Results: 1st is Intercept, 2nd is Beta+, 3rd is Beta-
  results <- list(
    intercept  = beta_hat[1, 1],
    beta_plus  = beta_hat[2, 1],
    beta_minus = beta_hat[3, 1],
    all_coefs  = beta_hat
  )

  return(results)
}
#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 3: load data and run
#++++++++++++++++++++++++++++++++++++++++++++++++
my_data <- read.csv("data(4).csv")
my_data <- as.matrix(my_data)
n <- nrow(my_data)

for (j in 1:4) {
  if (j == 1) {
    y <- as.matrix(my_data[2:n, 3])
    x <- as.matrix(my_data[2:n, 2])
  } else if (j == 2) {
    y <- as.matrix(my_data[2:n, 4])
    x <- as.matrix(my_data[2:n, 2])
  } else if (j == 3) {
    y <- as.matrix(my_data[2:n, 5])
    x <- as.matrix(my_data[2:n, 2])
  } else {
    y <- as.matrix(my_data[2:n, 6])
    x <- as.matrix(my_data[2:n, 2])
  }
  # Call the function using your 'y' and 'x' data, specifying leads and lags
  my_results <- asymmetric_dols(y = y, x = x, n_leads = 2, n_lags = 2)
  dols1   <- my_results$beta_plus
  dols2   <- my_results$beta_minus
  mout[j,1]   <- dols1
  mout[j,2]   <- dols2
  print(my_results)
}
print(mout)
