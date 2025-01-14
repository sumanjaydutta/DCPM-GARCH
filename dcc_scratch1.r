# Load necessary libraries
library(mvtnorm)   # For multivariate normal simulations
library(xts)       # For time series data manipulation
library(quantmod)  # For downloading stock data

# Step 1: Univariate GARCH(1,1) estimation
garch11_fit <- function(y, omega, alpha, beta) {
  T <- length(y)
  sigma2 <- rep(0, T)
  sigma2[1] <- var(y)
  
  for (t in 2:T) {
    sigma2[t] <- omega + alpha * y[t-1]^2 + beta * sigma2[t-1]
  }
  
  return(sigma2)
}

# Step 2: Log-likelihood function for GARCH(1,1)
garch11_loglik <- function(params, y) {
  omega <- params[1]
  alpha <- params[2]
  beta  <- params[3]
  
  sigma2 <- garch11_fit(y, omega, alpha, beta)
  
  # Return negative log-likelihood
  return(-sum(dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE)))
}

# Step 3: Estimate parameters for each asset using numerical optimization
estimate_garch11 <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, garch11_loglik, y = y, method = "L-BFGS-B", 
                        lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  
  return(optim_result$par)
}

# Step 4: Dynamic Conditional Correlation (DCC) function
dcc_update <- function(e, qbar, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))
  Qt[,,1] <- qbar
  
  for (t in 2:T) {
    Qt[,,t] <- (1 - a - b) * qbar + a * (e[t-1,] %*% t(e[t-1,])) + b * Qt[,,t-1]
    Rt[,,t] <- diag(1/sqrt(diag(Qt[,,t]))) %*% Qt[,,t] %*% diag(1/sqrt(diag(Qt[,,t])))
  }
  
  return(Rt)
}

# Step 5: DCC-GARCH estimation function
dcc_garch <- function(data, a, b) {
  N <- ncol(data)
  T <- nrow(data)
  
  # Fit GARCH(1,1) to each time series
  garch_params <- apply(data, 2, estimate_garch11)
  sigma2 <- matrix(0, T, N)
  e <- matrix(0, T, N)
  
  for (i in 1:N) {
    params <- garch_params[,i]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  # Step 6: Estimate the unconditional correlation matrix
  qbar <- cov(e)
  
  # Step 7: Update correlation matrix using DCC method
  Rt <- dcc_update(e, qbar, a, b)
  
  return(list(sigma2 = sigma2, correlation = Rt))
}

# Step 8: Simulating data to test DCC-GARCH
set.seed(123)
n <- 2
T <- 500
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.8, 0.8, 1), n, n)

# Simulate correlated returns
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# Step 9: Run the DCC-GARCH estimation
result <- dcc_garch(returns, a = 0.01, b = 0.95)

# Print results
result$sigma2  # Conditional variances
result$correlation[,,500]  # Correlation matrix at time T

#for params

# Step 5: Log-likelihood function for DCC
dcc_loglik <- function(params, e, qbar) {
  a <- params[1]
  b <- params[2]
  
  T <- nrow(e)
  N <- ncol(e)
  
  # Initialize Qt and Rt matrices
  Qt <- array(0, dim = c(N, N, T))
  Rt <- array(0, dim = c(N, N, T))
  Qt[,,1] <- qbar
  
  loglik <- 0  # To store log-likelihood
  
  # Calculate Qt and Rt matrices and compute log-likelihood
  for (t in 2:T) {
    Qt[,,t] <- (1 - a - b) * qbar + a * (e[t-1,] %*% t(e[t-1,])) + b * Qt[,,t-1]
    Rt[,,t] <- diag(1/sqrt(diag(Qt[,,t]))) %*% Qt[,,t] %*% diag(1/sqrt(diag(Qt[,,t])))
    
    # Log-likelihood based on standardized residuals and Rt
    loglik <- loglik + log(det(Rt[,,t])) + t(e[t,]) %*% solve(Rt[,,t]) %*% e[t,]
  }
  
  return(0.5 * loglik)  # Return half of the log-likelihood for symmetry
}

#
# Step 6: DCC parameter estimation function
estimate_dcc_params <- function(e, qbar) {
  # Initial guesses for a and b
  start_params <- c(a = 0.01, b = 0.95)
  
  # Optimization to estimate a and b
  optim_result <- optim(start_params, dcc_loglik, e = e, qbar = qbar, 
                        method = "L-BFGS-B", lower = c(1e-6, 1e-6), upper = c(0.999, 0.999))
  
  return(optim_result$par)
}

# Step 7: DCC-GARCH estimation function with parameter estimation
dcc_garch <- function(data) {
  N <- ncol(data)
  T <- nrow(data)
  
  # Fit GARCH(1,1) to each time series
  garch_params <- apply(data, 2, estimate_garch11)
  sigma2 <- matrix(0, T, N)
  e <- matrix(0, T, N)
  
  for (i in 1:N) {
    params <- garch_params[,i]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  # Estimate the unconditional correlation matrix
  qbar <- cov(e)
  
  # Estimate DCC parameters
  dcc_params <- estimate_dcc_params(e, qbar)
  a <- dcc_params[1]
  b <- dcc_params[2]
  
  # Update correlation matrix using DCC method with estimated parameters
  Rt <- dcc_update(e, qbar, a, b)
  
  return(list(sigma2 = sigma2, correlation = Rt, dcc_params = dcc_params))
}

# Simulating data to test DCC-GARCH
set.seed(123)
n <- 2
T <- 500
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.8, 0.8, 1), n, n)

# Simulate correlated returns
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# Run the DCC-GARCH estimation
result <- dcc_garch(returns)

# Print results
print(result$dcc_params)  # Estimated DCC parameters (a, b)
result$sigma2  # Conditional variances
result$correlation[,,500]  # Correlation matrix at time T


#experiment for comparison with dcc garch function

# Load necessary libraries
library(mvtnorm)
library(rmgarch)

# Step 1: Simulate correlated returns
set.seed(123)
n <- 2  # number of assets
T <- 500  # number of time points
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.8, 0.8, 1), n, n)
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# Step 2: Estimate DCC parameters with custom DCC-GARCH
result_custom <- dcc_garch(returns)
dcc_params_custom <- result_custom$dcc_params
print(paste("Custom DCC parameters: a =", dcc_params_custom[1], "b =", dcc_params_custom[2]))

# Step 3: Estimate DCC parameters with rmgarch package
uspec <- ugarchspec(mean.model = list(armaOrder = c(0, 0)), 
                    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                    distribution.model = "norm")
spec <- dccspec(uspec = multispec(replicate(n, uspec)), 
                dccOrder = c(1, 1), 
                distribution = "mvnorm")
fit <- dccfit(spec, data = returns)
dcc_params_package <- fit@mfit$matcoef[1:2, 1]
print(paste("Package DCC parameters: a =", dcc_params_package[1], "b =", dcc_params_package[2]))

# Step 4: Compare results
comparison <- data.frame(
  Method = c("Custom", "Package"),
  a = c(dcc_params_custom[1], dcc_params_package[1]),
  b = c(dcc_params_custom[2], dcc_params_package[2])
)
print(comparison)

# Step 5: Evaluate differences
diff_a <- abs(dcc_params_custom[1] - dcc_params_package[1])
diff_b <- abs(dcc_params_custom[2] - dcc_params_package[2])
print(paste("Difference in 'a':", diff_a))
print(paste("Difference in 'b':", diff_b))

#dcc garch using inverse covariance estimation 

# Load required libraries
library(glasso)      # For graphical lasso (Glasso)
library(MASS)        # For data simulation
library(Matrix)      # For sparse matrices
library(stats)       # For optimization

# Step 1: Simulate synthetic data for asset returns
set.seed(123)
n_assets <- 10
n_obs <- 200

# Simulate factor returns and idiosyncratic returns
factor_returns <- matrix(rnorm(n_obs * 5), ncol = 5)
factor_loadings <- matrix(runif(n_assets * 5, -1, 1), ncol = 5)
idiosyncratic_returns <- matrix(rnorm(n_obs * n_assets), ncol = n_assets)
asset_returns <- factor_returns %*% t(factor_loadings) + idiosyncratic_returns

# Step 2: Define functions for Graphical Lasso
estimate_precision_matrix <- function(data, lambda) {
  glasso_fit <- glasso(cov(data), rho = lambda)
  return(glasso_fit$wi)  # Return precision matrix
}

# Step 3: Define functions for GARCH(1,1)
garch11_fit <- function(y, omega, alpha, beta) {
  T <- length(y)
  sigma2 <- rep(0, T)
  sigma2[1] <- var(y)
  
  for (t in 2:T) {
    sigma2[t] <- omega + alpha * y[t-1]^2 + beta * sigma2[t-1]
  }
  
  return(sigma2)
}

# Step 4: Define functions for DCC GARCH
dcc_update <- function(e, qbar_inv, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  for (t in 2:T) {
    if (t > 2) {
      Qt_inv[,,t] <- (1 - b) * qbar_inv + a * (crossprod(e[t-1, , drop = FALSE]) - crossprod(e[t-2, , drop = FALSE])) + b * Qt_inv[,,t-1]
    } else {
      Qt_inv[,,t] <- (1 - b) * qbar_inv + a * (crossprod(e[t-1, , drop = FALSE])) + b * Qt_inv[,,t-1]
    }
    Rt[,,t] <- diag(1/sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1/sqrt(diag(solve(Qt_inv[,,t]))))
  }
  
  return(Rt)
}

# Step 5: Log-likelihood function for DCC-GARCH
dcc_loglik <- function(params, data) {
  a <- params[1]
  b <- params[2]
  lambda <- params[3]
  
  N <- ncol(data)
  T <- nrow(data)
  
  # Fit GARCH(1,1) to each time series
  garch_params <- t(apply(data, 2, function(y) {
    start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
    optim(start_params, function(params) {
      omega <- params[1]
      alpha <- params[2]
      beta <- params[3]
      sigma2 <- garch11_fit(y, omega, alpha, beta)
      -sum(dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE))
    }, method = "L-BFGS-B", lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))$par
  }))
  
  sigma2 <- matrix(0, T, N)
  e <- matrix(0, T, N)
  
  for (i in 1:N) {
    params <- garch_params[i,]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  # Estimate the unconditional precision matrix
  qbar_inv <- estimate_precision_matrix(e, lambda)
  
  # Update correlation matrix using DCC method
  Rt <- dcc_update(e, qbar_inv, a, b)
  
  # Compute log-likelihood for DCC-GARCH
  loglik <- 0
  for (t in 2:T) {
    loglik <- loglik - sum(dnorm(e[t,], mean = 0, sd = sqrt(diag(Rt[,,t])), log = TRUE))
  }
  
  return(-loglik)  # Return negative log-likelihood
}

# Step 6: Optimize DCC-GARCH parameters
estimate_dcc_parameters <- function(data) {
  start_params <- c(a = 0.01, b = 0.95, lambda = 0.1)
  optim_result <- optim(start_params, dcc_loglik, data = data, method = "L-BFGS-B",
                        lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  return(optim_result$par)
}

# Step 7: Run the DCC-GARCH estimation with parameter optimization
set.seed(123)
dcc_params <- estimate_dcc_parameters(asset_returns)
result <- estimate_dcc_garch(asset_returns, a = dcc_params[1], b = dcc_params[2], lambda = dcc_params[3])

# Print results
print("DCC-GARCH model estimation complete.")
print("Estimated parameters:")
print(dcc_params)
print("Correlation matrix at the last time point:")
print(result$correlation[,,n_obs])

#comparison with standard dcc garch parameters

# Install and load the rmgarch package if not already installed

library(rmgarch)

# Step 1: Define custom DCC-GARCH functions (as in previous code)
# Ensure to use the functions provided above

# Step 2: Fit the custom DCC-GARCH model and extract parameters
set.seed(123)
custom_dcc_params <- estimate_dcc_parameters(asset_returns)
custom_result <- estimate_dcc_garch(asset_returns, a = custom_dcc_params[1], b = custom_dcc_params[2], lambda = custom_dcc_params[3])

# Step 3: Fit the DCC-GARCH model using rmgarch
# Specify GARCH(1,1) specification for univariate GARCH
garch_spec <- ugarchspec(variance.model = list(garchOrder = c(1, 1)),
                         mean.model = list(armaOrder = c(0, 0)))

# Specify DCC specification
dcc_spec <- dccspec(uspec = multispec(replicate(n_assets, garch_spec)),
                    dccOrder = c(1, 1), distribution = "mvnorm")

# Fit the DCC-GARCH model
dcc_fit <- dccfit(dcc_spec, data = asset_returns)

# Extract DCC parameters
dcc_params_rmgarch <- coef(dcc_fit)
print("DCC-GARCH parameters from rmgarch:")
print(dcc_params_rmgarch)

# Compare parameters
print("Custom DCC-GARCH parameters:")
print(custom_dcc_params)

# Step 4: Compare custom parameters with rmgarch parameters
# Extract relevant parts for comparison
custom_a <- custom_dcc_params[1]
custom_b <- custom_dcc_params[2]

rmgarch_a <- dcc_params_rmgarch["a1"]
rmgarch_b <- dcc_params_rmgarch["b1"]

print("Comparison of parameters:")
cat(sprintf("Custom a: %f\n", custom_a))
cat(sprintf("Rmgarch a: %f\n", rmgarch_a))
cat(sprintf("Custom b: %f\n", custom_b))
cat(sprintf("Rmgarch b: %f\n", rmgarch_b))


#GGM adapted DCC

# Load necessary libraries
library(mvtnorm)   # For multivariate normal simulations
library(xts)       # For time series data manipulation
library(quantmod)  # For downloading stock data
library(glasso)    # For graphical lasso

# Step 1: Univariate GARCH(1,1) estimation
garch11_fit <- function(y, omega, alpha, beta) {
  T <- length(y)
  sigma2 <- rep(0, T)
  sigma2[1] <- var(y)
  
  for (t in 2:T) {
    sigma2[t] <- omega + alpha * y[t-1]^2 + beta * sigma2[t-1]
  }
  
  return(sigma2)
}

# Step 2: Log-likelihood function for GARCH(1,1)
garch11_loglik <- function(params, y) {
  omega <- params[1]
  alpha <- params[2]
  beta  <- params[3]
  
  sigma2 <- garch11_fit(y, omega, alpha, beta)
  
  # Return negative log-likelihood
  return(-sum(dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE)))
}

# Step 3: Estimate parameters for each asset using numerical optimization
estimate_garch11 <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, garch11_loglik, y = y, method = "L-BFGS-B", 
                        lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  
  return(optim_result$par)
}

# Step 4: Estimate the inverse covariance matrix using Graphical Lasso (Glasso)
estimate_inverse_covariance <- function(e) {
  # Apply Glasso to the idiosyncratic errors
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 5: Dynamic Conditional Correlation (DCC) function
dcc_update <- function(e, qbar_inv, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  for (t in 2:T) {
    if (t > 2) {
      e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)  # Make sure it's a row vector
      e_t_minus_2 <- matrix(e[t-2, ], nrow = 1)  # Make sure it's a row vector
      
      diff_cov <- t(e_t_minus_1) %*% e_t_minus_1 - t(e_t_minus_2) %*% e_t_minus_2
    } else {
      e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)  # Make sure it's a row vector
      diff_cov <- t(e_t_minus_1) %*% e_t_minus_1
    }
    
    # Ensure dimensions match
    if (dim(diff_cov)[1] != dim(qbar_inv)[1] || dim(diff_cov)[2] != dim(qbar_inv)[2]) {
      stop("Dimension mismatch between diff_cov and qbar_inv")
    }
    
    Qt_inv[,,t] <- (1 - a - b) * qbar_inv + a * diff_cov + b * Qt_inv[,,t-1]
    
    # Check if Qt_inv is positive definite
    if (any(eigen(Qt_inv[,,t])$values <= 0)) {
      Qt_inv[,,t] <- qbar_inv  # Revert to qbar_inv if not positive definite
      warning("Qt_inv is not positive definite. Using qbar_inv instead.")
    }
    
    Rt[,,t] <- diag(1 / sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1 / sqrt(diag(solve(Qt_inv[,,t]))))
  }
  
  return(Rt)
}

# Step 6: DCC-GARCH estimation function
dcc_garch <- function(data, a, b) {
  N <- ncol(data)
  T <- nrow(data)
  
  # Fit GARCH(1,1) to each time series
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, T, N)
  e <- matrix(0, T, N)
  
  for (i in 1:N) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  # Step 7: Estimate the inverse covariance matrix
  qbar_inv <- estimate_inverse_covariance(e)
  
  # Step 8: Update correlation matrix using DCC method
  Rt <- dcc_update(e, qbar_inv, a, b)
  
  return(list(sigma2 = sigma2, correlation = Rt))
}

# Step 9: Simulating data to test DCC-GARCH
set.seed(123)
n <- 2
T <- 500
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.8, 0.8, 1), n, n)

# Simulate correlated returns
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# Step 10: Run the DCC-GARCH estimation
result <- dcc_garch(returns, a = 0.01, b = 0.95)

# Print results
print("Conditional variances:")
print(result$sigma2)

print("Correlation matrix at time T:")
print(result$correlation[,,T])
#adapted DCC with param estimation
library(mvtnorm)   # For multivariate normal simulations
library(xts)       # For time series data manipulation
library(quantmod)  # For downloading stock data
library(glasso)    # For graphical lasso

# Step 1: Univariate GARCH(1,1) estimation
garch11_fit <- function(y, omega, alpha, beta) {
  T <- length(y)
  sigma2 <- rep(0, T)
  sigma2[1] <- var(y)
  
  for (t in 2:T) {
    sigma2[t] <- omega + alpha * y[t-1]^2 + beta * sigma2[t-1]
  }
  
  return(sigma2)
}

# Step 2: Log-likelihood function for GARCH(1,1)
garch11_loglik <- function(params, y) {
  omega <- params[1]
  alpha <- params[2]
  beta  <- params[3]
  
  sigma2 <- garch11_fit(y, omega, alpha, beta)
  
  # Return negative log-likelihood
  return(-sum(dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE)))
}

# Step 3: Estimate parameters for each asset using numerical optimization
estimate_garch11 <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, garch11_loglik, y = y, method = "L-BFGS-B", 
                        lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  
  return(optim_result$par)
}

# Step 4: Estimate the inverse covariance matrix using Graphical Lasso (Glasso)
estimate_inverse_covariance <- function(e) {
  # Apply Glasso to the idiosyncratic errors
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 5: Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, qbar_inv, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  for (t in 2:T) {
    if (t > 2) {
      e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)  # Make sure it's a row vector
      e_t_minus_2 <- matrix(e[t-2, ], nrow = 1)  # Make sure it's a row vector
      
      diff_cov <- t(e_t_minus_1) %*% e_t_minus_1 - t(e_t_minus_2) %*% e_t_minus_2
    } else {
      e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)  # Make sure it's a row vector
      diff_cov <- t(e_t_minus_1) %*% e_t_minus_1
    }
    
    # Ensure dimensions match
    if (dim(diff_cov)[1] != dim(qbar_inv)[1] || dim(diff_cov)[2] != dim(qbar_inv)[2]) {
      stop("Dimension mismatch between diff_cov and qbar_inv")
    }
    
    Qt_inv[,,t] <- (1 - a- b) * qbar_inv + a* diff_cov + b * Qt_inv[,,t-1]
    
    # Check if Qt_inv is positive definite
    if (any(eigen(Qt_inv[,,t])$values <= 0)) {
      Qt_inv[,,t] <- qbar_inv  # Revert to qbar_inv if not positive definite
      warning("Qt_inv is not positive definite. Using qbar_inv instead.")
    }
    
    Rt[,,t] <- diag(1 / sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1 / sqrt(diag(solve(Qt_inv[,,t]))))
  }
  
  return(Rt)
}

# Step 6: Define a Log-Likelihood Function for the DCC-GARCH Model
dcc_garch_loglik <- function(params, data) {
  a <- params[1]
  b <- params[2]
  
  # Fit GARCH(1,1) to each time series
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  # Estimate the inverse covariance matrix
  qbar_inv <- estimate_inverse_covariance(e)
  
  # Update correlation matrix using DCC method
  Rt <- dcc_update(e, qbar_inv, a, b)
  
  # Compute log-likelihood (negative to minimize)
  loglik <- -sum(dnorm(e, mean = 0, sd = sqrt(sigma2), log = TRUE))
  
  return(loglik)
}

# Step 7: Optimize DCC parameters
estimate_dcc_params <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1))
  
  return(optim_result$par)
}

# Step 8: DCC-GARCH estimation function
dcc_garch <- function(data) {
  # Estimate DCC parameters
  dcc_params <- estimate_dcc_params(data)
  a <- dcc_params[1]
  b <- dcc_params[2]
  
  # Fit GARCH(1,1) to each time series
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  # Estimate the inverse covariance matrix
  qbar_inv <- estimate_inverse_covariance(e)
  
  # Update correlation matrix using DCC method
  Rt <- dcc_update(e, qbar_inv, a, b)
  
  return(list(sigma2 = sigma2, correlation = Rt, dcc_params = dcc_params))
}

# Step 9: Simulate data to test DCC-GARCH
set.seed(123)
n <- 2
T <- 500
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.8, 0.8, 1), n, n)

# Simulate correlated returns
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# Step 10: Run the DCC-GARCH estimation
result <- dcc_garch(returns)

# Print results
print("Conditional variances:")
print(result$sigma2)

print("Correlation matrix at time T:")
print(result$correlation[,,T])

print("Estimated DCC parameters:")
print(result$dcc_params)

#adapted dcc with inverse cov est for idiosyncratic errors

library(mvtnorm)   # For multivariate normal simulations
library(xts)       # For time series data manipulation
library(quantmod)  # For downloading stock data
library(glasso)    # For graphical lasso
library(MASS)      # For pseudo-inverse (ginv)

# Step 1: Univariate GARCH(1,1) estimation
garch11_fit <- function(y, omega, alpha, beta) {
  T <- length(y)
  sigma2 <- rep(0, T)
  sigma2[1] <- var(y)  # Initial variance
  
  for (t in 2:T) {
    sigma2[t] <- omega + alpha * y[t-1]^2 + beta * sigma2[t-1]
  }
  
  return(sigma2)
}

# Step 2: Log-likelihood function for GARCH(1,1)
garch11_loglik <- function(params, y) {
  omega <- params[1]
  alpha <- params[2]
  beta  <- params[3]
  
  sigma2 <- garch11_fit(y, omega, alpha, beta)
  
  # Return negative log-likelihood
  return(-sum(dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE)))
}

# Step 3: Estimate parameters for each asset using numerical optimization
estimate_garch11 <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, garch11_loglik, y = y, method = "L-BFGS-B", 
                        lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  
  return(optim_result$par)
}

# Step 4: Estimate the inverse covariance matrix using Graphical Lasso (Glasso)
estimate_inverse_covariance <- function(e) {
  # Apply Glasso to the idiosyncratic errors
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 5: Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, qbar_inv, a, b, reg_lambda = 1e-5) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  for (t in 2:T) {
    if (t > 2) {
      e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)  # Row vector for t-1
      e_t_minus_2 <- matrix(e[t-2, ], nrow = 1)  # Row vector for t-2
      
      inv_cov_t_minus_1 <- tryCatch({
        solve(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k))
      }, error = function(err) {
        return(ginv(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k)))  # Use pseudo-inverse
      })
      
      inv_cov_t_minus_2 <- tryCatch({
        solve(t(e_t_minus_2) %*% e_t_minus_2 + diag(reg_lambda, k))
      }, error = function(err) {
        return(ginv(t(e_t_minus_2) %*% e_t_minus_2 + diag(reg_lambda, k)))  # Use pseudo-inverse
      })
      
      diff_inv_cov <- inv_cov_t_minus_1 - inv_cov_t_minus_2
    } else {
      e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)  # Row vector for t-1
      inv_cov_t_minus_1 <- solve(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k))
      diff_inv_cov <- inv_cov_t_minus_1
    }
    
    Qt_inv[,,t] <- (1 - b) * qbar_inv + b * diff_inv_cov + a * Qt_inv[,,t-1]
    
    # Check if Qt_inv is positive definite
    if (any(eigen(Qt_inv[,,t])$values <= 0)) {
      Qt_inv[,,t] <- qbar_inv  # Revert to qbar_inv if not positive definite
      warning("Qt_inv is not positive definite. Using qbar_inv instead.")
    }
    
    Rt[,,t] <- diag(1 / sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1 / sqrt(diag(solve(Qt_inv[,,t]))))
  }
  
  return(Rt)
}

# Step 6: Define a Log-Likelihood Function for the DCC-GARCH Model
dcc_garch_loglik <- function(params, data) {
  a <- params[1]
  b <- params[2]
  
  # Fit GARCH(1,1) to each time series
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  # Estimate the inverse covariance matrix
  qbar_inv <- estimate_inverse_covariance(e)
  
  # Update correlation matrix using DCC method
  Rt <- dcc_update(e, qbar_inv, a, b)
  
  # Compute log-likelihood (negative to minimize)
  loglik <- -sum(dnorm(e, mean = 0, sd = sqrt(sigma2), log = TRUE))
  
  return(loglik)
}

# Step 7: Optimize DCC parameters
estimate_dcc_params <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1))
  
  return(optim_result$par)
}

# Step 8: DCC-GARCH estimation function
dcc_garch <- function(data) {
  # Estimate DCC parameters
  dcc_params <- estimate_dcc_params(data)
  a <- dcc_params[1]
  b <- dcc_params[2]
  
  # Fit GARCH(1,1) to each time series
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  # Estimate the inverse covariance matrix
  qbar_inv <- estimate_inverse_covariance(e)
  
  # Update correlation matrix using DCC method
  Rt <- dcc_update(e, qbar_inv, a, b)
  
  return(list(sigma2 = sigma2, correlation = Rt, dcc_params = dcc_params, garch_params = garch_params))
}

# Step 9: Simulate data to test DCC-GARCH
set.seed(123)
n <- 2  # Number of assets
T <- 500  # Time points
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.8, 0.8, 1), n, n)

# Simulate correlated returns
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# Step 10: Run the DCC-GARCH estimation
result <- dcc_garch(returns)

# Print the estimated DCC parameters (a, b)
print(result$dcc_params)

# Print the estimated GARCH(1,1) parameters for each asset (omega, alpha, beta)
print(result$garch_params)

# Print the dynamic correlation matrices Rt
print(result$correlation)

#adapted dcc latest version with example for 30 stocks
# Load Required Libraries
library(mvtnorm)   # For multivariate normal simulations
library(glasso)    # For graphical lasso
library(xts)       # For time series data manipulation

# Step 1: Univariate GARCH(1,1) estimation
garch11_fit <- function(y, omega, alpha, beta) {
  T <- length(y)
  sigma2 <- rep(0, T)
  sigma2[1] <- var(y)
  
  for (t in 2:T) {
    sigma2[t] <- omega + alpha * y[t-1]^2 + beta * sigma2[t-1]
  }
  
  return(sigma2)
}

# Step 2: Log-likelihood function for GARCH(1,1)
garch11_loglik <- function(params, y) {
  omega <- params[1]
  alpha <- params[2]
  beta  <- params[3]
  
  sigma2 <- garch11_fit(y, omega, alpha, beta)
  
  # Return negative log-likelihood
  return(-sum(dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE)))
}

# Step 3: Estimate parameters for each asset using numerical optimization
estimate_garch11 <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, garch11_loglik, y = y, method = "L-BFGS-B", 
                        lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  
  return(optim_result$par)
}

# Step 4: Estimate the inverse covariance matrix using Graphical Lasso (Glasso)
estimate_inverse_covariance <- function(e) {
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 5: Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, qbar_inv, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  for (t in 2:T) {
    e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)
    
    diff_inv_cov <- t(e_t_minus_1) %*% e_t_minus_1
    
    Qt_inv[,,t] <- (1 - b) * qbar_inv + b * diff_inv_cov + a * Qt_inv[,,t-1]
    
    # Ensure Qt_inv is positive definite
    if (any(eigen(Qt_inv[,,t])$values <= 0)) {
      Qt_inv[,,t] <- qbar_inv  # Revert to qbar_inv if not positive definite
      warning("Qt_inv is not positive definite. Using qbar_inv instead.")
    }
    
    # Calculate the correlation matrix
    Rt[,,t] <- diag(1 / sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1 / sqrt(diag(solve(Qt_inv[,,t]))))
  }
  
  return(Rt)
}

# Step 6: Define a Log-Likelihood Function for the DCC-GARCH Model
dcc_garch_loglik <- function(params, data) {
  a <- params[1]
  b <- params[2]
  
  # Fit GARCH(1,1) to each time series
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  # Estimate the inverse covariance matrix
  qbar_inv <- estimate_inverse_covariance(e)
  
  # Update correlation matrix using DCC method
  Rt <- dcc_update(e, qbar_inv, a, b)
  
  # Compute log-likelihood (negative to minimize)
  loglik <- -sum(dnorm(e, mean = 0, sd = sqrt(sigma2), log = TRUE))
  
  return(loglik)
}

# Step 7: Optimize DCC parameters
estimate_dcc_params <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1))
  
  return(optim_result$par)
}

# Step 8: DCC-GARCH estimation function
dcc_garch <- function(data) {
  # Estimate DCC parameters
  dcc_params <- estimate_dcc_params(data)
  a <- dcc_params[1]
  b <- dcc_params[2]
  
  # Fit GARCH(1,1) to each time series
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  # Estimate the inverse covariance matrix
  qbar_inv <- estimate_inverse_covariance(e)
  
  # Update correlation matrix using DCC method
  Rt <- dcc_update(e, qbar_inv, a, b)
  
  # Calculate the conditional covariance matrix at time T
  conditional_covariance <- Rt[,,T] * diag(sqrt(sigma2[T, ])) %*% diag(sqrt(sigma2[T, ]))
  
  return(list(sigma2 = sigma2, correlation = Rt, conditional_covariance = conditional_covariance, dcc_params = dcc_params))
}

# Step 9: Simulate data for 30 stocks with stronger correlations
set.seed(123)
n <- 30
T <- 500
mu <- rep(0, n)
# Adjust the covariance matrix for stronger off-diagonal correlations
Sigma <- matrix(0.5, n, n)  # Base correlation
diag(Sigma) <- 1

# Simulate correlated returns
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# Step 10: Run the DCC-GARCH estimation
result <- dcc_garch(returns)

# Output the results
print("Conditional variances for each stock at time T:")
print(result$sigma2[T, ])

print("Correlation matrix at time T:")
print(result$correlation[,,T])

print("Conditional covariance matrix at time T:")
print(result$conditional_covariance)

print("Estimated DCC parameters:")
print(result$dcc_params)


#corrected adapted dcc code

# Load necessary libraries
library(mvtnorm)   # For multivariate normal simulations
library(xts)       # For time series data manipulation
library(quantmod)  # For downloading stock data
library(glasso)    # For graphical lasso
library(MASS)      # For pseudo-inverse (ginv)

# Step 1: Univariate GARCH(1,1) estimation
garch11_fit <- function(y, omega, alpha, beta) {
  T <- length(y)
  sigma2 <- rep(0, T)
  sigma2[1] <- var(y)  # Initial variance
  
  for (t in 2:T) {
    sigma2[t] <- omega + alpha * y[t-1]^2 + beta * sigma2[t-1]
  }
  
  return(sigma2)
}

# Step 2: Log-likelihood function for GARCH(1,1)
garch11_loglik <- function(params, y) {
  omega <- params[1]
  alpha <- params[2]
  beta  <- params[3]
  
  sigma2 <- garch11_fit(y, omega, alpha, beta)
  
  # Return negative log-likelihood
  return(-sum(dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE)))
}

# Step 3: Estimate parameters for each asset using numerical optimization
estimate_garch11 <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, garch11_loglik, y = y, method = "L-BFGS-B", 
                        lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  
  return(optim_result$par)
}

# Step 4: Estimate the inverse covariance matrix using Graphical Lasso (Glasso)
estimate_inverse_covariance <- function(e) {
  # Apply Glasso to the idiosyncratic errors
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 5: Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, qbar_inv, a, b, reg_lambda = 1e-5) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  # Store conditional covariance matrices
  conditional_cov <- array(0, dim = c(k, k, T))
  
  for (t in 2:T) {
    if (t > 2) {
      e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)  # Row vector for t-1
      e_t_minus_2 <- matrix(e[t-2, ], nrow = 1)  # Row vector for t-2
      
      inv_cov_t_minus_1 <- tryCatch({
        solve(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k))
      }, error = function(err) {
        return(ginv(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k)))  # Use pseudo-inverse
      })
      
      inv_cov_t_minus_2 <- tryCatch({
        solve(t(e_t_minus_2) %*% e_t_minus_2 + diag(reg_lambda, k))
      }, error = function(err) {
        return(ginv(t(e_t_minus_2) %*% e_t_minus_2 + diag(reg_lambda, k)))  # Use pseudo-inverse
      })
      
      diff_inv_cov <- inv_cov_t_minus_1 - inv_cov_t_minus_2
    } else {
      e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)  # Row vector for t-1
      inv_cov_t_minus_1 <- solve(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k))
      diff_inv_cov <- inv_cov_t_minus_1
    }
    
    Qt_inv[,,t] <- (1 - b) * qbar_inv + b * diff_inv_cov + a * Qt_inv[,,t-1]
    
    # Check if Qt_inv is positive definite
    if (any(eigen(Qt_inv[,,t])$values <= 0)) {
      Qt_inv[,,t] <- qbar_inv  # Revert to qbar_inv if not positive definite
      warning("Qt_inv is not positive definite. Using qbar_inv instead.")
    }
    
    Rt[,,t] <- diag(1 / sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1 / sqrt(diag(solve(Qt_inv[,,t]))))
    
    # Store the conditional covariance matrix
    conditional_cov[,,t] <- diag(sqrt(sigma2[t, ])) %*% Rt[,,t] %*% diag(sqrt(sigma2[t, ]))
  }
  
  return(list(Rt = Rt, conditional_cov = conditional_cov))
}

# Step 6: Define a Log-Likelihood Function for the DCC-GARCH Model
dcc_garch_loglik <- function(params, data) {
  a <- params[1]
  b <- params[2]
  
  # Fit GARCH(1,1) to each time series
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  # Estimate the inverse covariance matrix
  qbar_inv <- estimate_inverse_covariance(e)
  
  # Update correlation matrix using DCC method
  dcc_result <- dcc_update(e, qbar_inv, a, b)
  
  # Compute log-likelihood (negative to minimize)
  loglik <- -sum(dnorm(e, mean = 0, sd = sqrt(sigma2), log = TRUE))
  
  return(loglik)
}

# Step 7: Optimize DCC parameters
estimate_dcc_params <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1))
  
  return(optim_result$par)
}

# Step 8: DCC-GARCH estimation function
dcc_garch <- function(data) {
  # Estimate DCC parameters
  dcc_params <- estimate_dcc_params(data)
  a <- dcc_params[1]
  b <- dcc_params[2]
  
  # Fit GARCH(1,1) to each time series
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  # Estimate the inverse covariance matrix
  qbar_inv <- estimate_inverse_covariance(e)
  
  # Update correlation matrix using DCC method
  dcc_result <- dcc_update(e, qbar_inv, a, b)
  
  return(list(sigma2 = sigma2, correlation = dcc_result$Rt, conditional_cov = dcc_result$conditional_cov, dcc_params = dcc_params, garch_params = garch_params))
}

# Step 9: Simulate data to test DCC-GARCH with five stocks
set.seed(123)
n <- 5  # Number of assets
T <- 500  # Time points
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.8, 0.8, 0.8, 0.8,
                  0.8, 1, 0.8, 0.8, 0.8,
                  0.8, 0.8, 1, 0.8, 0.8,
                  0.8, 0.8, 0.8, 1, 0.8,
                  0.8, 0.8, 0.8, 0.8, 1), n, n)

# Simulate correlated returns
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# Apply the DCC-GARCH model
result <- dcc_garch(returns)

# Print the results
print(result$sigma2)               # Conditional variances
print(result$correlation)           # Conditional correlations
print(result$conditional_cov)       # Conditional covariances
print(result$dcc_params)            # DCC parameters
print(result$garch_params)          # GARCH parameters

# Step 10: Extract conditional covariance matrix at time T
conditional_cov_at_T <- result$conditional_cov[,,T]
print("Conditional Covariance Matrix at Time T:")
print(conditional_cov_at_T)

#adapted dcc - L-BFGS-B vs nelder mead

# Load necessary libraries
library(mvtnorm)

# Step 1: GARCH(1,1) log-likelihood function
garch11_loglik <- function(params, y) {
  omega <- params[1]
  alpha <- params[2]
  beta <- params[3]
  
  n <- length(y)
  sigma2 <- numeric(n)
  sigma2[1] <- var(y)
  for (t in 2:n) {
    sigma2[t] <- omega + alpha * y[t-1]^2 + beta * sigma2[t-1]
  }
  
  loglik <- -sum(dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE))
  return(loglik)
}

# Step 2: Fit GARCH(1,1) for each asset using Nelder-Mead
estimate_garch11_nm <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, garch11_loglik, y = y, method = "Nelder-Mead")
  return(list(par = optim_result$par, value = optim_result$value, convergence = optim_result$convergence))
}

# Step 2: Fit GARCH(1,1) for each asset using L-BFGS-B
estimate_garch11_lbfgsb <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, garch11_loglik, y = y, method = "L-BFGS-B", 
                        lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  return(list(par = optim_result$par, value = optim_result$value, convergence = optim_result$convergence))
}

# Step 4: DCC-GARCH log-likelihood function
dcc_garch_loglik <- function(params, data) {
  a <- params[1]
  b <- params[2]
  
  n <- ncol(data)
  T <- nrow(data)
  S <- cov(data)
  Qbar <- cov(data)
  Qt <- Qbar
  
  loglik <- 0
  for (t in 2:T) {
    Qt <- (1 - a - b) * Qbar + a * (data[t-1, ] %*% t(data[t-1, ])) + b * Qt
    Rt <- cov2cor(Qt)
    loglik <- loglik - 0.5 * (log(det(Rt)) + t(data[t, ]) %*% solve(Rt) %*% data[t, ])
  }
  
  return(-loglik)
}

# Step 5: Optimize DCC parameters using Nelder-Mead
estimate_dcc_params_nm <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "Nelder-Mead")
  return(list(par = optim_result$par, value = optim_result$value, convergence = optim_result$convergence))
}

# Step 5: Optimize DCC parameters using L-BFGS-B
estimate_dcc_params_lbfgsb <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1))
  return(list(par = optim_result$par, value = optim_result$value, convergence = optim_result$convergence))
}

# Step 6: DCC-GARCH estimation using Nelder-Mead
dcc_garch_nm <- function(data) {
  dcc_params <- estimate_dcc_params_nm(data)
  
  garch_params <- t(apply(data, 2, function(col) estimate_garch11_nm(col)$par))
  return(list(dcc_params = dcc_params, garch_params = garch_params))
}

# Step 6: DCC-GARCH estimation using L-BFGS-B
dcc_garch_lbfgsb <- function(data) {
  dcc_params <- estimate_dcc_params_lbfgsb(data)
  
  garch_params <- t(apply(data, 2, function(col) estimate_garch11_lbfgsb(col)$par))
  return(list(dcc_params = dcc_params, garch_params = garch_params))
}

# Step 7: Simulate data for testing
set.seed(123)
n <- 5  # Number of assets
T <- 500  # Time points
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.8, 0.8, 0.8, 0.8,
                  0.8, 1, 0.8, 0.8, 0.8,
                  0.8, 0.8, 1, 0.8, 0.8,
                  0.8, 0.8, 0.8, 1, 0.8,
                  0.8, 0.8, 0.8, 0.8, 1), n, n)

# Simulate multivariate normal returns
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# Step 8: Run Nelder-Mead DCC-GARCH
system.time(result_nm <- dcc_garch_nm(returns))  # Measure execution time
cat("Nelder-Mead Optimization:\n")
cat("DCC Parameters (a, b):", result_nm$dcc_params$par, "\n")
cat("Log-likelihood value:", -result_nm$dcc_params$value, "\n")
cat("Convergence status:", result_nm$dcc_params$convergence, "\n\n")

# Step 8: Run L-BFGS-B DCC-GARCH
system.time(result_lbfgsb <- dcc_garch_lbfgsb(returns))  # Measure execution time
cat("L-BFGS-B Optimization:\n")
cat("DCC Parameters (a, b):", result_lbfgsb$dcc_params$par, "\n")
cat("Log-likelihood value:", -result_lbfgsb$dcc_params$value, "\n")
cat("Convergence status:", result_lbfgsb$dcc_params$convergence, "\n\n")

#comparison using different dcc methods


