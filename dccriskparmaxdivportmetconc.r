# Required Libraries
library(mvtnorm)   # For multivariate normal simulations
library(glasso)    # For graphical lasso
library(xdcc)      # For xdcc package methods (COV, LS, NLS)

# Step 1: Adapted DCC-GARCH Method (Same as before)
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

# Step 2: Maximum Diversification Portfolio
max_div_portfolio <- function(cov_matrix) {
  inv_vol <- 1 / sqrt(diag(cov_matrix))
  weights <- inv_vol / sum(inv_vol)
  return(weights)
}

# Step 3: Risk Parity Portfolio
risk_parity_portfolio <- function(cov_matrix) {
  inv_var <- 1 / diag(cov_matrix)
  weights <- inv_var / sum(inv_var)
  return(weights)
}

# Step 4: XDCC Methods (COV, LS, NLS)
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

# Step 5: Run Simulation for Different Sample Sizes
set.seed(123)
n <- 500  # 500 stocks
Sigma <- matrix(0.5, n, n) + diag(0.5, n)  # Create a covariance matrix with moderate correlations

sample_sizes <- seq(50, 1000, by = 50)  # Sample sizes from 50 to 1000
results_table <- data.frame(SampleSize = sample_sizes, Strategy = NA, 
                            DCC_Return = NA, DCC_StdDev = NA, DCC_InfoRatio = NA, 
                            COV_Return = NA, COV_StdDev = NA, COV_InfoRatio = NA, 
                            LS_Return = NA, LS_StdDev = NA, LS_InfoRatio = NA,
                            NLS_Return = NA, NLS_StdDev = NA, NLS_InfoRatio = NA)

strategies <- c("Max Diversification", "Risk Parity")

for (s in sample_sizes) {
  cat("Running for sample size:", s, "\n")
  returns <- rmvnorm(s, mean = rep(0, n), sigma = Sigma)  # Simulate stock returns
  
  for (strategy in strategies) {
    cat("Running strategy:", strategy, "\n")
    
    if (strategy == "Max Diversification") {
      cov_matrix <- cov(returns)
      weights <- max_div_portfolio(cov_matrix)
    } else if (strategy == "Risk Parity") {
      cov_matrix <- cov(returns)
      weights <- risk_parity_portfolio(cov_matrix)
    }
    
    # Adapted DCC
    dcc_result <- dcc_garch(returns, weights)
    dcc_metrics <- portfolio_metrics(dcc_result$sigma2, dcc_result$correlation, weights, returns)
    
    # XDCC Methods
    xdcc_metrics <- xdcc_methods(returns, weights)
    
    # Store results in the table
    results_table[results_table$SampleSize == s & results_table$Strategy == strategy, "DCC_Return"] <- dcc_metrics$annualized_return
    results_table[results_table$SampleSize == s & results_table$Strategy == strategy, "DCC_StdDev"] <- dcc_metrics$annualized_std_dev
    results_table[results_table$SampleSize == s & results_table$Strategy == strategy, "DCC_InfoRatio"] <- dcc_metrics$info_ratio
    
    results_table[results_table$SampleSize == s & results_table$Strategy == strategy, "COV_Return"] <- xdcc_metrics$COV$annualized_return
    results_table[results_table$SampleSize == s & results_table$Strategy == strategy, "COV_StdDev"] <- xdcc_metrics$COV$annualized_std_dev
    results_table[results_table$SampleSize == s & results_table$Strategy == strategy, "COV_InfoRatio"] <- xdcc_metrics$COV$info_ratio
    
    results_table[results_table$SampleSize == s & results_table$Strategy == strategy, "LS_Return"] <- xdcc_metrics$LS$annualized_return
    results_table[results_table$SampleSize == s & results_table$Strategy == strategy, "LS_StdDev"] <- xdcc_metrics$LS$annualized_std_dev
    results_table[results_table$SampleSize == s & results_table$Strategy == strategy, "LS_InfoRatio"] <- xdcc_metrics$LS$info_ratio
    
    results_table[results_table$SampleSize == s & results_table$Strategy == strategy, "NLS_Return"] <- xdcc_metrics$NLS$annualized_return
    results_table[results_table$SampleSize == s & results_table$Strategy == strategy, "NLS_StdDev"] <- xdcc_metrics$NLS$annualized_std_dev
    results_table[results_table$SampleSize == s & results_table$Strategy == strategy, "NLS_InfoRatio"] <- xdcc_metrics$NLS$info_ratio
  }
}

# Display the final results table
print(results_table)
