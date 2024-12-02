# Load necessary libraries
library(mvtnorm)   # For multivariate normal simulations
library(xts)       # For time series data manipulation
library(quantmod)  # For downloading stock data
library(glasso)    # For graphical lasso
library(MASS)      # For pseudo-inverse (ginv)
library(ggplot2)   # For plotting
library(gridExtra) # For arranging multiple plots

# Function to simulate data using DCC-GARCH
simulate_dcc_garch_data <- function(num_assets, num_points, garch_params, dcc_params) {
  # Unpack GARCH parameters
  omega <- garch_params[1]
  alpha <- garch_params[2]
  beta <- garch_params[3]
  
  # Unpack DCC parameters
  a <- dcc_params[1]
  b <- dcc_params[2]

  # Initialize matrices
  returns <- matrix(0, nrow = num_points, ncol = num_assets)
  sigma2 <- matrix(0, nrow = num_points, ncol = num_assets)
  e <- matrix(0, nrow = num_points, ncol = num_assets)
  qbar_inv <- diag(num_assets)  # Initial inverse covariance matrix
  
  # Simulate GARCH(1,1) process for each asset
  for (i in 1:num_assets) {
    # Initialize variance
    sigma2[1, i] <- 0.1  # Initial variance
    for (t in 2:num_points) {
      sigma2[t, i] <- omega + alpha * returns[t - 1, i]^2 + beta * sigma2[t - 1, i]
      e[t, i] <- rnorm(1) * sqrt(sigma2[t, i])
      returns[t, i] <- e[t, i]
    }
  }
  
  # Dynamic Conditional Correlation
  Rt <- array(0, dim = c(num_assets, num_assets, num_points))
  Rt[,,1] <- diag(num_assets)  # Start with identity matrix

  for (t in 2:num_points) {
    e_t <- e[t, ]
    
    # Update inverse covariance matrix using DCC
    qbar_inv <- (1 - b) * qbar_inv + b * (1 / (e_t %*% t(e_t)))
    Rt[,,t] <- diag(1 / sqrt(diag(solve(qbar_inv)))) %*% solve(qbar_inv) %*% diag(1 / sqrt(diag(solve(qbar_inv))))
    
    # Ensure positive definiteness
    if (any(eigen(Rt[,,t])$values <= 0)) {
      Rt[,,t] <- diag(num_assets)
      warning("Rt is not positive definite. Using identity matrix instead.")
    }
  }
  
  # Return both the simulated returns and the time-varying covariances
  list(returns = returns, Rt = Rt)
}

# Example of usage
set.seed(123)  # For reproducibility
num_assets <- 5
num_points <- 500

# Specify GARCH parameters
garch_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)

# Specify DCC parameters
dcc_params <- c(a = 0.01, b = 0.95)

# Simulate DCC-GARCH data
simulation_results <- simulate_dcc_garch_data(num_assets, num_points, garch_params, dcc_params)
simulated_data <- simulation_results$returns
dcc_cov_matrices <- simulation_results$Rt

# 1. Compute the unconditional covariance matrix of the returns
unconditional_cov <- cov(simulated_data)

# 2. Compute the average of the DCC covariance matrices across time
average_dcc_cov <- apply(dcc_cov_matrices, c(1, 2), mean)

# 3. Compare the two matrices
difference <- unconditional_cov - average_dcc_cov
print("Difference between unconditional covariance and average DCC covariance:")
print(difference)

# Check if the matrices are approximately equal (small differences due to numerical errors)
if (all(abs(difference) < 1)) {
  print("The average DCC covariance is approximately equal to the unconditional covariance.")
} else {
  print("The average DCC covariance is not equal to the unconditional covariance.")
}

# Plotting simulated data
volatility_plots <- lapply(1:num_assets, function(i) {
  data.frame(Time = 1:num_points, Returns = simulated_data[, i]) %>%
    ggplot(aes(x = Time, y = Returns)) +
    geom_line() +
    ggtitle(paste("Simulated Returns for Asset", i)) +
    theme_minimal()
})

# Combine plots
grid.arrange(grobs = volatility_plots, ncol = 2)
