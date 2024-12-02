# Required Libraries
library(mvtnorm)   # For multivariate normal simulations
library(glasso)    # For graphical lasso
library(xdcc)      # For xdcc package methods (COV, LS, NLS)
library(quantmod)  # For time series data manipulation

# Step 1: Adapted DCC-GARCH Method (As discussed earlier)
garch11_fit <- function(y, omega, alpha, beta) {
  T <- length(y)
  sigma2 <- rep(0, T)
  sigma2[1] <- var(y)
  
  for (t in 2:T) {
    sigma2[t] <- omega + alpha * y[t-1]^2 + beta * sigma2[t-1]
  }
  
  return(sigma2)
}

estimate_garch11 <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, garch11_loglik, y = y, method = "L-BFGS-B", 
                        lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  
  return(optim_result$par)
}

dcc_update <- function(e, qbar_inv, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  for (t in 2:T) {
    e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)
    diff_cov <- t(e_t_minus_1) %*% e_t_minus_1
    Qt_inv[,,t] <- (1 - a - b) * qbar_inv + a * diff_cov + b * Qt_inv[,,t-1]
    
    Rt[,,t] <- diag(1 / sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1 / sqrt(diag(solve(Qt_inv[,,t]))))
  }
  
  return(Rt)
}

dcc_garch <- function(data, weights) {
  T <- nrow(data)
  train_size <- floor(0.8 * T)
  train_data <- data[1:train_size, ]
  test_data <- data[(train_size + 1):T, ]
  
  garch_params <- t(apply(train_data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(train_data), ncol(train_data))
  e <- matrix(0, nrow(train_data), ncol(train_data))
  
  for (i in 1:ncol(train_data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(train_data[,i], params[1], params[2], params[3])
    e[,i] <- train_data[,i] / sqrt(sigma2[,i])
  }
  
  qbar_inv <- estimate_inverse_covariance(e)
  
  start_params <- c(a = 0.01, b = 0.95)
  dcc_params <- optim(start_params, dcc_garch_loglik, data = train_data, method = "L-BFGS-B", 
                      lower = c(0, 0), upper = c(1, 1))$par
  a <- dcc_params[1]
  b <- dcc_params[2]
  
  forecast_sigma2 <- matrix(0, nrow(test_data), ncol(test_data))
  forecast_correlation <- array(0, dim = c(ncol(test_data), ncol(test_data), nrow(test_data)))
  
  for (i in 1:ncol(test_data)) {
    params <- garch_params[i, ]
    forecast_sigma2[,i] <- garch11_fit(test_data[,i], params[1], params[2], params[3])
  }
  
  Rt <- dcc_update(e, qbar_inv, a, b)
  
  for (t in 1:nrow(test_data)) {
    forecast_correlation[,,t] <- Rt[,,t]
  }
  
  return(list(sigma2 = forecast_sigma2, correlation = forecast_correlation, dcc_params = dcc_params))
}

portfolio_metrics <- function(forecast_sigma2, forecast_correlation, weights, data) {
  T <- nrow(data)
  portfolio_returns <- rowSums(data %*% weights)
  
  portfolio_variance <- numeric(T)
  for (t in 1:T) {
    cov_matrix_t <- diag(sqrt(forecast_sigma2[t,])) %*% forecast_correlation[,,t] %*% diag(sqrt(forecast_sigma2[t,]))
    portfolio_variance[t] <- t(weights) %*% cov_matrix_t %*% weights
  }
  
  portfolio_std_dev <- sqrt(portfolio_variance)
  annualized_return <- mean(portfolio_returns) * 252
  annualized_std_dev <- sqrt(mean(portfolio_std_dev^2)) * sqrt(252)
  
  info_ratio <- annualized_return / annualized_std_dev
  
  return(list(annualized_return = annualized_return, annualized_std_dev = annualized_std_dev, info_ratio = info_ratio))
}

# Step 2: XDCC Methods (COV, LS, NLS)
xdcc_methods <- function(data, weights) {
  # Fit the xdcc models: COV, LS, NLS
  cov_est <- xdcc(data, method = "COV")
  ls_est <- xdcc(data, method = "LS")
  nls_est <- xdcc(data, method = "NLS")
  
  # Extract forecasted correlation matrices
  cov_corr <- cov_est$corr
  ls_corr <- ls_est$corr
  nls_corr <- nls_est$corr
  
  # For simplicity, assume constant variances (identity matrix) since xdcc doesn't estimate variances
  n <- ncol(data)
  forecast_sigma2 <- matrix(1, nrow = nrow(data), ncol = n)
  
  # Compute portfolio metrics for each method
  cov_metrics <- portfolio_metrics(forecast_sigma2, cov_corr, weights, data)
  ls_metrics <- portfolio_metrics(forecast_sigma2, ls_corr, weights, data)
  nls_metrics <- portfolio_metrics(forecast_sigma2, nls_corr, weights, data)
  
  return(list(COV = cov_metrics, LS = ls_metrics, NLS = nls_metrics))
}

# Step 3: Compare Adapted DCC with XDCC Methods
set.seed(123)
n <- 2
T <- 500
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.8, 0.8, 1), n, n)

# Simulate correlated returns
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# Portfolio weights (for example, equal weights)
weights <- c(0.5, 0.5)

# Estimate the DCC-GARCH model
dcc_result <- dcc_garch(returns, weights)

# Compute portfolio metrics for adapted DCC
dcc_metrics <- portfolio_metrics(dcc_result$sigma2, dcc_result$correlation, weights, returns[(nrow(returns) - 100):nrow(returns), ])

# Compute portfolio metrics for xdcc methods (COV, LS, NLS)
xdcc_metrics <- xdcc_methods(returns[(nrow(returns) - 100):nrow(returns), ], weights)

# Step 4: Print Comparison Results
print("Adapted DCC Metrics:")
print(dcc_metrics)

print("XDCC COV Metrics:")
print(xdcc_metrics$COV)

print("XDCC LS Metrics:")
print(xdcc_metrics$LS)

print("XDCC NLS Metrics:")
print(xdcc_metrics$NLS)
