# Load necessary packages
library(rugarch)

# Step 1: Fit GARCH(1,1) model to each time series
fit_garch <- function(returns) {
  spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                     mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                     distribution.model = "norm")
  garch_fit <- ugarchfit(spec = spec, data = returns)
  
  sigma2 <- sigma(garch_fit)^2  # Conditional variance
  residuals <- residuals(garch_fit, standardize = TRUE)  # Standardized residuals
  
  return(list(sigma2 = sigma2, residuals = residuals))
}

# Step 2: DCC-GARCH Model

# Function to compute DCC log-likelihood and return time-varying Qt matrices
dcc_garch_log_likelihood <- function(params, residuals) {
  a <- params[1]
  b <- params[2]
  
  T <- nrow(residuals)
  N <- ncol(residuals)
  
  Q_bar <- cov(residuals)  # Unconditional covariance matrix of residuals
  Qt <- Q_bar
  Qt_list <- list()  # List to store Qt matrices over time
  
  log_likelihood <- 0
  
  for (t in 2:T) {
    Qt <- (1 - a - b) * Q_bar + a * (residuals[t-1, ] %*% t(residuals[t-1, ])) + b * Qt
    Rt <- diag(1 / sqrt(diag(Qt))) %*% Qt %*% diag(1 / sqrt(diag(Qt)))  # Dynamic correlation matrix
    
    log_likelihood <- log_likelihood + log(det(Rt)) + t(residuals[t, ]) %*% solve(Rt) %*% residuals[t, ]
    
    Qt_list[[t]] <- Qt  # Store Qt for averaging later
  }
  
  return(list(log_likelihood = log_likelihood, Qt_list = Qt_list))
}

# Fit DCC model to residuals
fit_dcc_garch <- function(residuals) {
  initial_params <- c(0.05, 0.9)  # Initial guesses for DCC parameters
  result <- optim(initial_params, function(params) dcc_garch_log_likelihood(params, residuals)$log_likelihood,
                  method = "L-BFGS-B", lower = c(1e-6, 1e-6), upper = c(1, 1))
  
  # Return both the parameters and the list of Qt matrices
  Qt_result <- dcc_garch_log_likelihood(result$par, residuals)
  
  return(list(params = result$par, Qt_list = Qt_result$Qt_list))
}

# Function to compute the time-average of the Qt matrices
average_Qt <- function(Qt_list) {
  T <- length(Qt_list)
  Qt_sum <- Qt_list[[2]]  # Initialize with the first Qt
  
  for (t in 3:T) {
    Qt_sum <- Qt_sum + Qt_list[[t]]
  }
  
  Qt_avg <- Qt_sum / (T - 1)  # Time-average (excluding the first value)
  return(Qt_avg)
}

# Function to compute Frobenius norm of the difference between matrices
frobenius_norm <- function(A, B) {
  return(norm(A - B, type = "F"))
}

# Full DCC-GARCH model
dcc_garch <- function(returns_matrix) {
  N <- ncol(returns_matrix)  # Number of time series
  T <- nrow(returns_matrix)
  
  residuals_matrix <- matrix(0, T, N)
  
  # Step 1: Fit univariate GARCH models and get standardized residuals
  for (i in 1:N) {
    garch_result <- fit_garch(returns_matrix[, i])
    residuals_matrix[, i] <- garch_result$residuals
  }
  
  # Step 2: Fit DCC model to standardized residuals
  dcc_result <- fit_dcc_garch(residuals_matrix)
  
  # Calculate time-averaged Qt
  Qt_avg <- average_Qt(dcc_result$Qt_list)
  
  # Compare Qt_avg with Q_bar (unconditional covariance matrix)
  Q_bar <- cov(residuals_matrix)
  
  # Compute Frobenius norm of the difference
  frob_norm_diff <- frobenius_norm(Qt_avg, Q_bar)
  
  return(list(dcc_params = dcc_result$params, Qt_avg = Qt_avg, Q_bar = Q_bar, frob_norm_diff = frob_norm_diff))
}

# Example usage
set.seed(42)
T <- 1000  # Number of observations
N <- 50    # Number of time series

# Generate synthetic data
returns_matrix <- matrix(rnorm(T * N), ncol = N)

# Fit the DCC-GARCH model
dcc_result <- dcc_garch(returns_matrix)

# Print the results
print("Estimated DCC parameters:")
print(dcc_result$dcc_params)

print("Frobenius norm of the difference between time-averaged Qt and Q_bar:")
print(dcc_result$frob_norm_diff)

#test for inverse Qbar

# Load necessary packages
library(rugarch)

# Step 1: Fit GARCH(1,1) model to each time series
fit_garch <- function(returns) {
  spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                     mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                     distribution.model = "norm")
  garch_fit <- ugarchfit(spec = spec, data = returns)
  
  sigma2 <- sigma(garch_fit)^2  # Conditional variance
  residuals <- residuals(garch_fit, standardize = TRUE)  # Standardized residuals
  
  return(list(sigma2 = sigma2, residuals = residuals))
}

# Step 2: DCC-GARCH Model

# Function to compute DCC log-likelihood and return time-varying Qt matrices
dcc_garch_log_likelihood <- function(params, residuals) {
  a <- params[1]
  b <- params[2]
  
  T <- nrow(residuals)
  N <- ncol(residuals)
  
  Q_bar <- cov(residuals)  # Unconditional covariance matrix of residuals
  Qt <- Q_bar
  Qt_list <- list()  # List to store Qt matrices over time
  
  log_likelihood <- 0
  
  for (t in 2:T) {
    Qt <- (1 - a - b) * Q_bar + a * (residuals[t-1, ] %*% t(residuals[t-1, ])) + b * Qt
    Rt <- diag(1 / sqrt(diag(Qt))) %*% Qt %*% diag(1 / sqrt(diag(Qt)))  # Dynamic correlation matrix
    
    log_likelihood <- log_likelihood + log(det(Rt)) + t(residuals[t, ]) %*% solve(Rt) %*% residuals[t, ]
    
    Qt_list[[t]] <- Qt  # Store Qt for averaging later
  }
  
  return(list(log_likelihood = log_likelihood, Qt_list = Qt_list))
}

# Fit DCC model to residuals
fit_dcc_garch <- function(residuals) {
  initial_params <- c(0.05, 0.9)  # Initial guesses for DCC parameters
  result <- optim(initial_params, function(params) dcc_garch_log_likelihood(params, residuals)$log_likelihood,
                  method = "L-BFGS-B", lower = c(1e-6, 1e-6), upper = c(1, 1))
  
  # Return both the parameters and the list of Qt matrices
  Qt_result <- dcc_garch_log_likelihood(result$par, residuals)
  
  return(list(params = result$par, Qt_list = Qt_result$Qt_list))
}

# Function to compute the time-average of the inverse of Qt matrices
average_inverse_Qt <- function(Qt_list) {
  T <- length(Qt_list)
  inv_Qt_sum <- solve(Qt_list[[2]])  # Initialize with the inverse of the first Qt
  
  for (t in 3:T) {
    inv_Qt_sum <- inv_Qt_sum + solve(Qt_list[[t]])
  }
  
  inv_Qt_avg <- inv_Qt_sum / (T - 1)  # Time-average of the inverse (excluding the first value)
  return(inv_Qt_avg)
}

# Function to compute Frobenius norm of the difference between matrices
frobenius_norm <- function(A, B) {
  return(norm(A - B, type = "F"))
}

# Full DCC-GARCH model
dcc_garch <- function(returns_matrix) {
  N <- ncol(returns_matrix)  # Number of time series
  T <- nrow(returns_matrix)
  
  residuals_matrix <- matrix(0, T, N)
  
  # Step 1: Fit univariate GARCH models and get standardized residuals
  for (i in 1:N) {
    garch_result <- fit_garch(returns_matrix[, i])
    residuals_matrix[, i] <- garch_result$residuals
  }
  
  # Step 2: Fit DCC model to standardized residuals
  dcc_result <- fit_dcc_garch(residuals_matrix)
  
  # Calculate time-averaged inverse Qt
  inv_Qt_avg <- average_inverse_Qt(dcc_result$Qt_list)
  
  # Compare the inverse of Q_bar (unconditional covariance matrix)
  Q_bar <- cov(residuals_matrix)
  inv_Q_bar <- solve(Q_bar)
  
  # Compute Frobenius norm of the difference between inverse averages
  frob_norm_diff <- frobenius_norm(inv_Qt_avg, inv_Q_bar)
  
  return(list(dcc_params = dcc_result$params, inv_Qt_avg = inv_Qt_avg, inv_Q_bar = inv_Q_bar, frob_norm_diff = frob_norm_diff))
}

# Example usage
set.seed(42)
T <- 1000  # Number of observations
N <- 50    # Number of time series

# Generate synthetic data
returns_matrix <- matrix(rnorm(T * N), ncol = N)

# Fit the DCC-GARCH model
dcc_result <- dcc_garch(returns_matrix)

# Print the results
print("Estimated DCC parameters:")
print(dcc_result$dcc_params)

print("Frobenius norm of the difference between the average of inverse Qt and inverse Q_bar:")
print(dcc_result$frob_norm_diff)

#comparison for inverse Qbar with rmgarch DCC function

# Load necessary packages
library(rugarch)
library(rmgarch)

# Function to fit a univariate GARCH(1,1) model and extract the residuals
fit_garch <- function(returns) {
  spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                     mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                     distribution.model = "norm")
  garch_fit <- ugarchfit(spec = spec, data = returns)
  
  residuals <- residuals(garch_fit, standardize = TRUE)  # Standardized residuals
  return(residuals)
}

# Function to compute the Frobenius norm of the difference between matrices
frobenius_norm <- function(A, B) {
  return(norm(A - B, type = "F"))
}

# Step 1: Fit DCC-GARCH using rmgarch
dcc_garch_rmgarch <- function(returns_matrix) {
  N <- ncol(returns_matrix)  # Number of time series
  
  # Step 2: Define and fit DCC-GARCH model
  uspec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                      mean.model = list(armaOrder = c(0, 0)),
                      distribution.model = "norm")
  
  # DCC-GARCH specification
  spec <- dccspec(uspec = multispec(replicate(N, uspec)), dccOrder = c(1, 1), distribution = "mvnorm")
  
  # Fit the DCC model to the returns data
  fit <- dccfit(spec, data = returns_matrix)
  
  # Step 3: Extract the Q_t matrices
  Qt_list <- rcor(fit, type = "Q")
  
  # Calculate the unconditional covariance matrix (Q_bar)
  Q_bar <- cov(returns_matrix)
  inv_Q_bar <- solve(Q_bar)
  
  # Step 4: Compute the average of the inverses of Q_t matrices
  T <- dim(Qt_list)[3]
  inv_Qt_sum <- solve(Qt_list[, , 1])  # Initialize with the inverse of the first Qt matrix
  
  for (t in 2:T) {
    inv_Qt_sum <- inv_Qt_sum + solve(Qt_list[, , t])
  }
  
  inv_Qt_avg <- inv_Qt_sum / T  # Time-average of the inverses
  
  # Step 5: Compute the Frobenius norm of the difference between inv_Qt_avg and inv_Q_bar
  frob_norm_diff <- frobenius_norm(inv_Qt_avg, inv_Q_bar)
  
  return(list(dcc_fit = fit, inv_Qt_avg = inv_Qt_avg, inv_Q_bar = inv_Q_bar, frob_norm_diff = frob_norm_diff))
}

# Example usage
set.seed(42)
T <- 1000  # Number of observations
N <- 50    # Number of time series

# Generate synthetic data for 50 variables
returns_matrix <- matrix(rnorm(T * N), ncol = N)

# Fit the DCC-GARCH model and perform the comparison
dcc_result <- dcc_garch_rmgarch(returns_matrix)

# Print the results
print("Frobenius norm of the difference between the average of inverse Qt and inverse Q_bar:")
print(dcc_result$frob_norm_diff)

#test for inverse Qbar using Adapted DCC (using cumulative information)
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
  # Check for NA, NaN, or Inf in the input
  if (any(is.na(e)) || any(is.nan(e)) || any(is.infinite(e))) {
    stop("Input contains NA, NaN, or Inf values.")
  }
  
  # Apply Glasso to the idiosyncratic errors
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 5: Dynamic Conditional Correlation (DCC) update function with Graphical Lasso
dcc_update <- function(e, qbar_inv, a, b, sigma2) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  # Store conditional covariance matrices
  conditional_cov <- array(0, dim = c(k, k, T))
  
  for (t in 2:T) {
    if (t > 2) {
      # Use cumulative information up to t-1 and t-2
      e_t_minus_1 <- matrix(e[1:(t-1), ], nrow = t-1)  # Cumulative up to t-1
      e_t_minus_2 <- matrix(e[1:(t-2), ], nrow = t-2)  # Cumulative up to t-2
      
      # Estimate inverse covariance matrix using Glasso for t-1
      if (nrow(e_t_minus_1) >= 2) {
        glasso_t_minus_1 <- glasso(cov(e_t_minus_1), rho = 0.1)
        inv_cov_t_minus_1 <- glasso_t_minus_1$wi  # Inverse covariance matrix for t-1
      } else {
        inv_cov_t_minus_1 <- qbar_inv  # Fallback to qbar_inv if not enough data
      }

      # Estimate inverse covariance matrix using Glasso for t-2
      if (nrow(e_t_minus_2) >= 2) {
        glasso_t_minus_2 <- glasso(cov(e_t_minus_2), rho = 0.1)
        inv_cov_t_minus_2 <- glasso_t_minus_2$wi  # Inverse covariance matrix for t-2
      } else {
        inv_cov_t_minus_2 <- qbar_inv  # Fallback to qbar_inv if not enough data
      }
      
      diff_inv_cov <- inv_cov_t_minus_1 - inv_cov_t_minus_2
    } else {
      # Use cumulative information up to t-1
      e_t_minus_1 <- matrix(e[1:(t-1), ], nrow = t-1)  # Cumulative up to t-1
      
      # Estimate inverse covariance matrix using Glasso for t-1
      if (nrow(e_t_minus_1) >= 2) {
        glasso_t_minus_1 <- glasso(cov(e_t_minus_1), rho = 0.1)
        inv_cov_t_minus_1 <- glasso_t_minus_1$wi  # Inverse covariance matrix for t-1
      } else {
        inv_cov_t_minus_1 <- qbar_inv  # Fallback to qbar_inv if not enough data
      }
      
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
  
  return(list(Rt = Rt, conditional_cov = conditional_cov, Qt_inv = Qt_inv))
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
  dcc_result <- dcc_update(e, qbar_inv, a, b, sigma2)
  
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
  dcc_result <- dcc_update(e, qbar_inv, a, b, sigma2)
  
  # Calculate average of inverse Qt matrices
  avg_inverse_Qt <- apply(dcc_result$Qt_inv, c(1, 2), mean)
  
  # Calculate Frobenius norm of the difference between avg_inverse_Qt and qbar_inv
  frobenius_norm <- norm(avg_inverse_Qt - qbar_inv, type = "F")
  
  return(list(dcc_result = dcc_result, frobenius_norm = frobenius_norm))
}

# Step 9: Synthetic Data Generation
set.seed(123)
n <- 1000  # Number of observations
k <- 5     # Number of variables

# Generate synthetic data from a multivariate normal distribution
mu <- rep(0, k)
sigma <- diag(rep(1, k))  # Identity covariance matrix
data <- rmvnorm(n, mean = mu, sigma = sigma)

# Step 10: Run DCC-GARCH estimation
dcc_results <- dcc_garch(data)

# Step 11: Display Frobenius norm result
print(paste("Frobenius Norm of the Difference:", dcc_results$frobenius_norm))

#test for inverse qbar using adapted dcc (single period info)


