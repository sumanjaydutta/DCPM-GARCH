# Load required libraries
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
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 5: BEKK parameter estimation and covariance update function
bekk_update <- function(e, omega, A, B) {
  T <- nrow(e)
  k <- ncol(e)
  
  H <- array(0, dim = c(k, k, T))
  
  # Initial covariance matrix
  H[,,1] <- diag(apply(e, 2, var))
  
  for (t in 2:T) {
    H[,,t] <- omega + A %*% (e[t-1, ] %*% t(e[t-1, ])) + B %*% H[,,t-1]
    
    # Ensure positive definiteness
    if (any(eigen(H[,,t])$values <= 0)) {
      H[,,t] <- H[,,1]  # Revert to initial covariance if not positive definite
      warning("Covariance matrix is not positive definite. Using initial covariance instead.")
    }
  }
  
  return(H)
}

# Step 6: Define a Log-Likelihood Function for the BEKK-GARCH Model
bekk_garch_loglik <- function(params, data) {
  omega <- params[1]
  A <- diag(params[2:3])  # Assuming A is a diagonal matrix
  B <- diag(params[4:5])  # Assuming B is a diagonal matrix
  
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
  
  # Update covariance matrix using BEKK method
  H <- bekk_update(e, omega, A, B)
  
  # Compute log-likelihood (negative to minimize)
  loglik <- -sum(dnorm(e, mean = 0, sd = sqrt(sigma2), log = TRUE))
  
  return(loglik)
}

# Step 7: Optimize BEKK parameters
estimate_bekk_params <- function(data) {
  start_params <- c(omega = 0.1, A1 = 0.1, A2 = 0.1, B1 = 0.1, B2 = 0.1)
  optim_result <- optim(start_params, bekk_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0, 0, 0, 0), upper = c(1, 1, 1, 1, 1))
  
  return(optim_result$par)
}

# Step 8: BEKK-GARCH estimation function
bekk_garch <- function(data) {
  # Estimate BEKK parameters
  bekk_params <- estimate_bekk_params(data)
  omega <- bekk_params[1]
  A <- diag(bekk_params[2:3])
  B <- diag(bekk_params[4:5])
  
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
  
  # Update covariance matrix using BEKK method
  H <- bekk_update(e, omega, A, B)
  
  return(list(sigma2 = sigma2, covariance = H, bekk_params = bekk_params))
}

# Step 9: Simulate data to test BEKK-GARCH
set.seed(123)
n <- 2
T <- 500
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.8, 0.8, 1), n, n)

# Simulate correlated returns
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# Step 10: Run the BEKK-GARCH estimation
result <- bekk_garch(returns)

# Print results
print("Conditional variances:")
print(result$sigma2)

print("Covariance matrix at time T:")
print(result$covariance[,,T])

print("Estimated BEKK parameters:")
print(result$bekk_params)
