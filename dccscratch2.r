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
    e[,i] <- data[,i] / sqrt(sigma2[,i])  # Standardized residuals
  }
  
  # Estimate the inverse covariance matrix
  qbar_inv <- estimate_inverse_covariance(e)
  
  # Update correlation matrix using DCC method
  dcc_result <- dcc_update(e, qbar_inv, a, b)
  
  # Compute log-likelihood using the inverse of the conditional covariance
  T <- nrow(data)
  k <- ncol(data)
  loglik <- 0
  
  for (t in 1:T) {
    Qt_inv <- dcc_result$Rt[,,t]  # Use the inverse of the conditional covariance matrix
    e_t <- matrix(e[t, ], ncol = 1)  # Convert to column vector
    
    # Calculate the log-likelihood for the t-th observation
    loglik_t <- -0.5 * (k * log(2 * pi) - log(det(Qt_inv)) + t(e_t) %*% Qt_inv %*% e_t)
    loglik <- loglik + loglik_t
  }
  
  return(-loglik)  # Return negative log-likelihood
}

# Step 7: Optimize DCC parameters
estimate_dcc_params <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1), 
                        constraints = list(a + b < 1))
  
  return(optim_result$par)
}

# Example usage
# Assuming 'data' is your multivariate time series matrix
# data <- your_data_matrix_here

# Estimate DCC parameters
dcc_params <- estimate_dcc_params(data)
cat("Optimal DCC parameters (a, b):", dcc_params, "\n")

# Get Q_t and H_t estimates using the last DCC update
final_dcc_result <- dcc_update(e, qbar_inv, dcc_params[1], dcc_params[2])
Q_ts <- final_dcc_result$conditional_cov
H_ts <- apply(Q_ts, c(1, 2), mean)  # Mean of conditional covariances across time

# Print estimates
cat("Q_t estimates:\n", Q_ts, "\n")
cat("H_t estimates:\n", H_ts, "\n")

#dcc garch basic implementation

# Load necessary libraries
library(mvtnorm)  # For multivariate normal simulations

# Univariate GARCH(1,1) fitting function
garch11_fit <- function(y, omega, alpha, beta) {
  T <- length(y)
  sigma2 <- rep(0, T)
  sigma2[1] <- var(y)  # Initial variance
  
  for (t in 2:T) {
    sigma2[t] <- omega + alpha * y[t-1]^2 + beta * sigma2[t-1]
  }
  
  return(sigma2)
}

# Log-likelihood function for GARCH(1,1)
garch11_loglik <- function(params, y) {
  omega <- params[1]
  alpha <- params[2]
  beta  <- params[3]
  
  sigma2 <- garch11_fit(y, omega, alpha, beta)
  
  # Return negative log-likelihood
  return(-sum(dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE)))
}

# Function to estimate GARCH(1,1) parameters
estimate_garch11 <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, garch11_loglik, y = y, method = "L-BFGS-B", 
                        lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  
  return(optim_result$par)
}

# DCC update function
dcc_update <- function(e, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  # Initialize matrices for conditional covariance
  Q <- array(0, dim = c(k, k, T))
  R <- array(0, dim = c(k, k, T))
  
  # Initial Q matrix based on sample covariance
  Q[,,1] <- cov(e)
  
  for (t in 2:T) {
    # Calculate the correlation matrix at time t-1
    R[,,t-1] <- diag(1 / sqrt(diag(Q[,,t-1]))) %*% Q[,,t-1] %*% diag(1 / sqrt(diag(Q[,,t-1])))
    
    # Update Q using DCC equations
    Q[,,t] <- (1 - a - b) * cov(e) + a * (e[t-1, ] %*% t(e[t-1, ])) + b * Q[,,t-1]
    
    # Ensure Q is positive definite
    if (any(eigen(Q[,,t])$values <= 0)) {
      Q[,,t] <- cov(e)  # Revert to sample covariance if not positive definite
    }
  }
  
  return(list(Q = Q, R = R))
}

# Log-likelihood function for DCC-GARCH
dcc_garch_loglik <- function(params, data) {
  a <- params[1]
  b <- params[2]
  
  # Fit GARCH(1,1) to each time series
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[, i] <- garch11_fit(data[, i], params[1], params[2], params[3])
    e[, i] <- data[, i] / sqrt(sigma2[, i])  # Standardized residuals
  }
  
  # Update DCC parameters
  dcc_result <- dcc_update(e, a, b)
  
  # Compute log-likelihood using the inverse of the conditional covariance
  T <- nrow(data)
  k <- ncol(data)
  loglik <- 0
  
  for (t in 1:T) {
    Qt <- dcc_result$Q[,,t]  # Use the conditional covariance matrix
    e_t <- matrix(e[t, ], ncol = 1)  # Convert to column vector
    
    # Calculate the log-likelihood for the t-th observation
    loglik_t <- -0.5 * (k * log(2 * pi) - log(det(Qt)) + t(e_t) %*% solve(Qt) %*% e_t)
    loglik <- loglik + loglik_t
  }
  
  return(-loglik)  # Return negative log-likelihood
}

# Function to estimate DCC parameters
estimate_dcc_params <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1), 
                        constraints = list(a + b < 1))
  
  return(optim_result$par)
}

# Example usage with simulated data
set.seed(123)
n <- 500  # Number of observations
k <- 3    # Number of assets

# Simulate some multivariate normal data
data <- matrix(rnorm(n * k), ncol = k)

# Estimate DCC parameters
dcc_params <- estimate_dcc_params(data)
cat("Optimal DCC parameters (a, b):", dcc_params, "\n")

# Get Q_t estimates using the last DCC update
final_dcc_result <- dcc_update(data, dcc_params[1], dcc_params[2])
Q_ts <- final_dcc_result$Q
cat("Q_t estimates:\n", Q_ts[,,n], "\n")  # Q_t for the last time step

#dcc basic without dnorm

# Load necessary libraries
library(mvtnorm)  # For multivariate normal distributions

# Step 1: GARCH(1,1) model estimation for each time series
garch11_fit <- function(y, omega, alpha, beta) {
  T <- length(y)
  sigma2 <- numeric(T)
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
  beta <- params[3]
  
  sigma2 <- garch11_fit(y, omega, alpha, beta)
  
  # Return negative log-likelihood
  return(-sum(dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE)))
}

# Step 3: Estimate GARCH(1,1) parameters for each asset
estimate_garch11 <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, garch11_loglik, y = y, method = "L-BFGS-B", 
                        lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  
  return(optim_result$par)
}

# Step 4: DCC update function
dcc_update <- function(e, qbar, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))  # Correlation matrices
  
  Qt[,,1] <- qbar
  
  for (t in 2:T) {
    Qt[,,t] <- (1 - a - b) * qbar + a * (t(e[t-1, ]) %*% e[t-1, ]) + b * Qt[,,t-1]
    diag_inv <- diag(1 / sqrt(diag(Qt[,,t])))
    Rt[,,t] <- diag_inv %*% Qt[,,t] %*% diag_inv  # Normalize to get correlation matrix
  }
  
  return(Rt)
}

# Step 5: Log-likelihood function for the DCC-GARCH model
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
    e[,i] <- data[,i] / sqrt(sigma2[,i])  # Standardized residuals
  }
  
  # Calculate Q-bar (average covariance matrix of residuals)
  qbar <- cov(e)
  
  # Update correlation matrices using DCC method
  Rt <- dcc_update(e, qbar, a, b)
  
  # Compute log-likelihood using time-varying R_t
  T <- nrow(data)
  k <- ncol(data)
  loglik <- 0
  
  for (t in 1:T) {
    e_t <- matrix(e[t, ], ncol = 1)  # Convert to column vector
    R_t <- Rt[,,t]
    # Calculate the log-likelihood for the t-th observation using multivariate normal distribution
    loglik_t <- -0.5 * (k * log(2 * pi) + log(det(R_t)) + t(e_t) %*% solve(R_t) %*% e_t)
    loglik <- loglik + loglik_t
  }
  
  return(-loglik)  # Return negative log-likelihood
}

# Step 6: Optimize DCC parameters
estimate_dcc_params <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1), 
                        constraints = list(a + b < 1))
  
  # Final estimates for a and b
  return(optim_result$par)
}

# Example usage:
# Generate synthetic data for two assets
set.seed(123)
n <- 1000
data <- cbind(rnorm(n), rnorm(n))  # Two independent standard normal variables

# Estimate DCC parameters
dcc_params <- estimate_dcc_params(data)
print(dcc_params)

#dcpm 
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

# Step 4: Modified Dynamic Conditional Correlation (DCC) update function to return R_t inverse
dcc_update <- function(e, qbar, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt <- array(0, dim = c(k, k, T))
  Rt_inv <- array(0, dim = c(k, k, T))  # Array to store R_t inverses
  Qt[,,1] <- qbar
  
  for (t in 2:T) {
    Qt[,,t] <- (1 - a - b) * qbar + a * (t(e[t-1, ]) %*% e[t-1, ]) + b * Qt[,,t-1]
    diag_inv <- diag(1 / sqrt(diag(Qt[,,t])))
    Rt <- diag_inv %*% Qt[,,t] %*% diag_inv  # Normalized correlation matrix R_t
    Rt_inv[,,t] <- diag(sqrt(diag(Qt[,,t]))) %*% solve(Qt[,,t]) %*% diag(sqrt(diag(Qt[,,t])))  # R_t inverse
  }
  
  return(list(Qt = Qt, Rt_inv = Rt_inv))
}

# Step 5: Modified Log-Likelihood Function for the DCC-GARCH Model
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
    e[,i] <- data[,i] / sqrt(sigma2[,i])  # Standardized residuals
  }
  
  # Calculate Q-bar (average covariance matrix of residuals)
  qbar <- cov(e)
  
  # Update correlation matrix using DCC method with R_t inverse
  dcc_result <- dcc_update(e, qbar, a, b)
  Qt <- dcc_result$Qt
  Rt_inv <- dcc_result$Rt_inv
  
  # Compute log-likelihood using precomputed R_t inverse
  T <- nrow(data)
  k <- ncol(data)
  loglik <- 0
  
  for (t in 1:T) {
    Q_t <- Qt[,,t]
    R_t_inv <- Rt_inv[,,t]
    e_t <- matrix(e[t, ], ncol = 1)  # Convert to column vector
    
    # Calculate the log-likelihood for the t-th observation using R_t inverse
    loglik_t <- -0.5 * (k * log(2 * pi) + log(det(Q_t)) + t(e_t) %*% R_t_inv %*% e_t)
    loglik <- loglik + loglik_t
  }
  
  return(-loglik)  # Return negative log-likelihood
}

# Step 6: Optimize DCC parameters
estimate_dcc_params <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1), 
                        constraints = list(a + b < 1))
  
  return(optim_result$par)
}

# Example usage with simulated data
# Simulate some multivariate data or load real data as needed
# For example:
# data <- rmvnorm(100, mean = rep(0, 5), sigma = diag(5))

# Estimate DCC parameters for the data
# dcc_params <- estimate_dcc_params(data)

#dcc params with Ht, Qt and Rts estimates
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

# Step 4: Modified Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, qbar, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))  # Array to store R_t matrices
  Rt_inv <- array(0, dim = c(k, k, T))  # Array to store R_t inverses
  Ht <- array(0, dim = c(k, k, T))  # Array to store H_t matrices
  
  Qt[,,1] <- qbar
  
  for (t in 2:T) {
    Qt[,,t] <- (1 - a - b) * qbar + a * (t(e[t-1, ]) %*% e[t-1, ]) + b * Qt[,,t-1]
    diag_inv <- diag(1 / sqrt(diag(Qt[,,t])))
    Rt[,,t] <- diag_inv %*% Qt[,,t] %*% diag_inv  # Normalized correlation matrix R_t
    Rt_inv[,,t] <- diag(sqrt(diag(Qt[,,t]))) %*% solve(Qt[,,t]) %*% diag(sqrt(diag(Qt[,,t])))  # R_t inverse
    Ht[,,t] <- diag(sqrt(diag(Qt[,,t]))) %*% Rt[,,t] %*% diag(sqrt(diag(Qt[,,t])))  # Conditional covariance matrix H_t
  }
  
  return(list(Qt = Qt, Rt = Rt, Rt_inv = Rt_inv, Ht = Ht))
}

# Step 5: Modified Log-Likelihood Function for the DCC-GARCH Model
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
    e[,i] <- data[,i] / sqrt(sigma2[,i])  # Standardized residuals
  }
  
  # Calculate Q-bar (average covariance matrix of residuals)
  qbar <- cov(e)
  
  # Update correlation and covariance matrices using DCC method
  dcc_result <- dcc_update(e, qbar, a, b)
  Qt <- dcc_result$Qt
  Rt <- dcc_result$Rt
  Rt_inv <- dcc_result$Rt_inv
  Ht <- dcc_result$Ht
  
  # Compute log-likelihood using precomputed R_t inverse
  T <- nrow(data)
  k <- ncol(data)
  loglik <- 0
  
  for (t in 1:T) {
    Q_t <- Qt[,,t]
    R_t_inv <- Rt_inv[,,t]
    e_t <- matrix(e[t, ], ncol = 1)  # Convert to column vector
    
    # Calculate the log-likelihood for the t-th observation using R_t inverse
    loglik_t <- -0.5 * (k * log(2 * pi) + log(det(Q_t)) + t(e_t) %*% R_t_inv %*% e_t)
    loglik <- loglik + loglik_t
  }
  
  return(-loglik)  # Return negative log-likelihood
}

# Step 6: Optimize DCC parameters and output R_t, H_t, and Q_t at the end
estimate_dcc_params <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1), 
                        constraints = list(a + b < 1))
  
  # Final estimates for a and b
  a <- optim_result$par[1]
  b <- optim_result$par[2]
  
  # Recalculate DCC matrices using final parameters
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])  # Standardized residuals
  }
  
  # Calculate Q-bar and final DCC matrices
  qbar <- cov(e)
  dcc_result <- dcc_update(e, qbar, a, b)
  
  # Return the final values for Rt, Ht, and Qt
  list(a = a, b = b, Rt = dcc_result$Rt, Ht = dcc_result$Ht, Qt = dcc_result$Qt)
}

#dcpm final code

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

# Step 4: Modified Dynamic Conditional Correlation (DCC) update function to return time-varying R_t
dcc_update <- function(e, qbar, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))  # Array to store R_t matrices
  Rt_inv <- array(0, dim = c(k, k, T))  # Array to store R_t inverses
  Ht <- array(0, dim = c(k, k, T))  # Array to store H_t matrices
  
  Qt[,,1] <- qbar
  
  for (t in 2:T) {
    Qt[,,t] <- (1 - a - b) * qbar + a * (t(e[t-1, ]) %*% e[t-1, ]) + b * Qt[,,t-1]
    diag_inv <- diag(1 / sqrt(diag(Qt[,,t])))
    Rt[,,t] <- diag_inv %*% Qt[,,t] %*% diag_inv  # Normalized correlation matrix R_t
    Rt_inv[,,t] <- diag(sqrt(diag(Qt[,,t]))) %*% solve(Qt[,,t]) %*% diag(sqrt(diag(Qt[,,t])))  # R_t inverse
    Ht[,,t] <- diag(sqrt(diag(Qt[,,t]))) %*% Rt[,,t] %*% diag(sqrt(diag(Qt[,,t])))  # Conditional covariance matrix H_t
  }
  
  return(list(Qt = Qt, Rt = Rt, Rt_inv = Rt_inv, Ht = Ht))
}

# Step 5: Modified Log-Likelihood Function for the DCC-GARCH Model using time-varying R_t
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
    e[,i] <- data[,i] / sqrt(sigma2[,i])  # Standardized residuals
  }
  
  # Calculate Q-bar (average covariance matrix of residuals)
  qbar <- cov(e)
  
  # Update correlation and covariance matrices using DCC method
  dcc_result <- dcc_update(e, qbar, a, b)
  Rt <- dcc_result$Rt
  Rt_inv <- dcc_result$Rt_inv
  Qt <- dcc_result$Qt
  
  # Compute log-likelihood using time-varying R_t
  T <- nrow(data)
  k <- ncol(data)
  loglik <- 0
  
  for (t in 1:T) {
    R_t <- Rt[,,t]
    R_t_inv <- Rt_inv[,,t]
    e_t <- matrix(e[t, ], ncol = 1)  # Convert to column vector
    
    # Calculate the log-likelihood for the t-th observation using R_t and its inverse
    loglik_t <- -0.5 * (k * log(2 * pi) + log(det(R_t)) + t(e_t) %*% R_t_inv %*% e_t)
    loglik <- loglik + loglik_t
  }
  
  return(-loglik)  # Return negative log-likelihood
}

# Step 6: Optimize DCC parameters and output R_t, H_t, and Q_t at the end
estimate_dcc_params <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1), 
                        constraints = list(a + b < 1))
  
  # Final estimates for a and b
  a <- optim_result$par[1]
  b <- optim_result$par[2]
  
  # Recalculate DCC matrices using final parameters
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])  # Standardized residuals
  }
  
  # Calculate Q-bar and final DCC matrices
  qbar <- cov(e)
  dcc_result <- dcc_update(e, qbar, a, b)
  
  # Return the final values for Rt, Ht, and Qt
  list(a = a, b = b, Rt = dcc_result$Rt, Ht = dcc_result$Ht, Qt = dcc_result$Qt)
}

#dcpm final corrected code

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

# Step 2: Estimate parameters for each asset using numerical optimization
estimate_garch11 <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, function(params) {
    omega <- params[1]
    alpha <- params[2]
    beta  <- params[3]
    sigma2 <- garch11_fit(y, omega, alpha, beta)
    # Return negative log-likelihood
    return(-sum(log(sigma2) + (y^2) / sigma2) / 2)
  }, method = "L-BFGS-B", lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  
  return(optim_result$par)
}

# Step 3: Estimate the inverse covariance matrix using Graphical Lasso (Glasso)
estimate_inverse_covariance <- function(e) {
  # Apply Glasso to the idiosyncratic errors
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 4: Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, qbar_inv, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  # Store conditional covariance matrices
  conditional_cov <- array(0, dim = c(k, k, T))
  
  for (t in 2:T) {
    e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)  # Row vector for t-1
    inv_cov_t_minus_1 <- tryCatch({
      solve(t(e_t_minus_1) %*% e_t_minus_1)
    }, error = function(err) {
      return(ginv(t(e_t_minus_1) %*% e_t_minus_1))  # Use pseudo-inverse
    })
    
    Qt_inv[,,t] <- (1 - b) * qbar_inv + b * inv_cov_t_minus_1 + a * Qt_inv[,,t-1]
    
    # Check if Qt_inv is positive definite
    if (any(eigen(Qt_inv[,,t])$values <= 0)) {
      Qt_inv[,,t] <- qbar_inv  # Revert to qbar_inv if not positive definite
      warning("Qt_inv is not positive definite. Using qbar_inv instead.")
    }
    
    # Calculate Rt from Qt_inv
    Rt[,,t] <- diag(1 / sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1 / sqrt(diag(solve(Qt_inv[,,t]))))
    
    # Store the conditional covariance matrix
    conditional_cov[,,t] <- diag(sqrt(sigma2[t, ])) %*% Rt[,,t] %*% diag(sqrt(sigma2[t, ]))
  }
  
  return(list(Rt = Rt, conditional_cov = conditional_cov))
}

# Step 5: Define a Log-Likelihood Function for the DCC-GARCH Model
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
  Rt <- dcc_result$Rt
  
  # Compute log-likelihood using the full expression
  k <- ncol(data)  # Number of assets
  loglik <- -0.5 * sum(log(det(Rt[,,2:nrow(Rt)])) + rowSums(e[2:nrow(e),] %*% (solve(Rt[,,2:nrow(Rt)]) * e[2:nrow(e),])) + k * log(2 * pi))
  
  return(-loglik)  # Return negative log-likelihood for optimization
}

# Step 6: Optimize DCC parameters
estimate_dcc_params <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1), 
                        method = "L-BFGS-B", 
                        control = list(fnscale = -1), # Maximize log-likelihood
                        gr = NULL)
  
  return(optim_result$par)
}

# Step 7: DCC-GARCH estimation function
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

#dcpm final corrected code:

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

# Step 2: Estimate parameters for each asset using numerical optimization
estimate_garch11 <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, function(params) {
    omega <- params[1]
    alpha <- params[2]
    beta  <- params[3]
    sigma2 <- garch11_fit(y, omega, alpha, beta)
    # Return negative log-likelihood
    return(-sum(log(sigma2) + (y^2) / sigma2) / 2)
  }, method = "L-BFGS-B", lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  
  return(optim_result$par)
}

# Step 3: Estimate the inverse covariance matrix using Graphical Lasso (Glasso)
estimate_inverse_covariance <- function(e) {
  # Apply Glasso to the idiosyncratic errors
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 4: Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, qbar_inv, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Qt <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  Qt[,,1] <- solve(qbar_inv)  # Q_t for the first time period
  
  # Store conditional covariance matrices
  conditional_cov <- array(0, dim = c(k, k, T))
  
  for (t in 2:T) {
    e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)  # Row vector for t-1
    inv_cov_t_minus_1 <- tryCatch({
      solve(t(e_t_minus_1) %*% e_t_minus_1)
    }, error = function(err) {
      return(ginv(t(e_t_minus_1) %*% e_t_minus_1))  # Use pseudo-inverse
    })
    
    Qt_inv[,,t] <- (1 - b) * qbar_inv + b * inv_cov_t_minus_1 + a * Qt_inv[,,t-1]
    
    # Calculate D_t
    D_t <- diag(sqrt(diag(Qt[,,t])))
    
    # Calculate R_t and R_t_inv
    R_t <- solve(D_t) %*% Qt[,,t] %*% solve(D_t)
    R_t_inv <- D_t %*% Qt_inv[,,t] %*% D_t
    
    # Store the conditional covariance matrix
    conditional_cov[,,t] <- R_t
    
    # Check if Qt_inv is positive definite
    if (any(eigen(Qt_inv[,,t])$values <= 0)) {
      Qt_inv[,,t] <- qbar_inv  # Revert to qbar_inv if not positive definite
      warning("Qt_inv is not positive definite. Using qbar_inv instead.")
    }
  }
  
  return(list(Rt = R_t, Rt_inv = R_t_inv, conditional_cov = conditional_cov))
}

# Step 5: Define a Log-Likelihood Function for the DCC-GARCH Model
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
  Rt <- dcc_result$Rt
  Rt_inv <- dcc_result$Rt_inv
  
  # Compute log-likelihood using the full expression with Rt and Rt_inv
  k <- ncol(data)  # Number of assets
  loglik <- -0.5 * sum(log(det(Rt[,,2:nrow(Rt)])) + 
                       rowSums(e[2:nrow(e),] %*% (Rt_inv[,,2:nrow(Rt_inv)] * e[2:nrow(e),])) + 
                       k * log(2 * pi))
  
  return(-loglik)  # Return negative log-likelihood for optimization
}

# Step 6: Optimize DCC parameters
estimate_dcc_params <- function(data) {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1), 
                        control = list(fnscale = -1), # Maximize log-likelihood
                        gr = NULL)
  
  return(optim_result$par)
}

# Step 7: DCC-GARCH estimation function
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
  
  return(list(sigma2 = sigma2, correlation = dcc_result$Rt_inv, conditional_cov = dcc_result$conditional_cov, dcc_params = dcc_params, garch_params = garch_params))
}

# Example usage with simulated data
set.seed(123)
n <- 1000  # Number of observations
p <- 3     # Number of assets

# Simulate correlated asset returns
data <- matrix(rnorm(n * p), nrow = n, ncol = p)
dcc_results <- dcc_garch(data)

# View results
print(dcc_results$sigma2)         # GARCH variances
print(dcc_results$correlation)     # Conditional correlations (R_t^{-1})
print(dcc_results$conditional_cov) # Conditional covariance matrices
