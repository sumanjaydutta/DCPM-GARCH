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

