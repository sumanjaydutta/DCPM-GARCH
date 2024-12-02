#testing dcc function from rmgarch for Qt and Qbar

# Required Libraries

library(fGarch)
library(rmgarch)  # For DCC-GARCH modeling
library(MASS)     # For pseudo-inverse (ginv)

# Step 1: GARCH(1,1) Estimation Function
estimate_garch11 <- function(data) {
  garch_fit <- garchFit(~ garch(1, 1), data = data, trace = FALSE)
  params <- c(mu = coef(garch_fit)["mu"], 
              alpha = coef(garch_fit)["alpha1"], 
              beta = coef(garch_fit)["beta1"])
  return(params)
}

# Step 2: GARCH(1,1) Conditional Variance Calculation
garch11_fit <- function(data, mu, alpha, beta) {
  T <- length(data)
  sigma2 <- numeric(T)
  sigma2[1] <- var(data)
  
  for (t in 2:T) {
    sigma2[t] <- alpha * (data[t-1] - mu)^2 + beta * sigma2[t-1]
  }
  
  return(sigma2)
}

# Step 3: DCC-GARCH estimation function
dcc_garch <- function(data) {
  # Define the GARCH model for each series
  spec <- dccspec(uspec = multispec(replicate(ncol(data), 
    garchSpec(model = list(variance.model = list(model = "sGARCH", 
    garchOrder = c(1, 1)), mean.model = list(armaOrder = c(0, 0)))))), 
    dccOrder = c(1, 1), distribution = "mvnorm")
  
  fit <- dccfit(spec, data)
  return(fit)
}

# Step 4: Calculate Frobenius norm of differences between Q_t and Q_bar
calculate_frobenius_norm <- function(fit, qbar) {
  # Extract the conditional covariance matrices Q_t
  Q_t <- fit@mfit$Q
  
  # Calculate the Frobenius norm of the difference between each Q_t and Q_bar
  frobenius_norms <- numeric(dim(Q_t)[3])  # Create a vector to store norms
  for (t in 1:dim(Q_t)[3]) {
    frobenius_norms[t] <- norm(Q_t[,,t] - qbar, type = "F")
  }
  
  return(frobenius_norms)
}

# Step 5: Main function to simulate data and perform DCC-GARCH estimation
simulate_and_analyze <- function(n_assets, n_obs, qbar) {
  # Simulate multivariate normal returns
  set.seed(123)  # For reproducibility
  data <- matrix(rnorm(n_assets * n_obs), nrow = n_obs, ncol = n_assets)
  
  # Fit the DCC-GARCH model
  fit <- dcc_garch(data)
  
  # Calculate Frobenius norms
  frobenius_norms <- calculate_frobenius_norm(fit, qbar)
  
  return(frobenius_norms)
}

plot_frobenius_norms <- function(frobenius_norms) {
  # Create a time index for the x-axis
  time_index <- 1:length(frobenius_norms)
  
  # Create a data frame for plotting
  df <- data.frame(Time_Index = time_index, Frobenius_Norm = frobenius_norms)
  
  # Create the plot
  ggplot(df, aes(x = Time_Index, y = Frobenius_Norm)) +
    geom_line(color = "blue") +
    geom_point(color = "red") +
    labs(title = "Frobenius Norms vs. Time Index", x = "Time Index", y = "Frobenius Norm") +
    theme_minimal()
}

# Usage: Set parameters for simulation
n_assets <- 3  # Number of assets
n_obs <- 500   # Number of observations
qbar <- diag(1, n_assets)  # Example: Identity matrix for Q_bar

# Run simulation and analyze results
frobenius_norms <- simulate_and_analyze(n_assets, n_obs, qbar)

# Print Frobenius norms
print(frobenius_norms)

# Plot Frobenius norms vs. the defined time index
plot_frobenius_norms(frobenius_norms)

#testing adapted dcc or dcp method

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
  
  return(list(sigma2 = sigma2, correlation = dcc_result$Rt, conditional_cov = dcc_result$conditional_cov, 
              dcc_params = dcc_params, garch_params = garch_params, Qt_inv = dcc_result$Qt_inv, qbar_inv = qbar_inv))
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
                  0.8, 0.8, 0.8, 0.8, 1), nrow = n)

# Simulate multivariate normal returns
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# Step 10: Apply the DCC-GARCH model
result <- dcc_garch(returns)

# Print results
print("Estimated Conditional Correlation Matrices:")
print(result$correlation)

print("Estimated Conditional Covariance Matrices:")
print(result$conditional_cov)

print("DCC Parameters:")
print(result$dcc_params)

print("GARCH Parameters:")
print(result$garch_params)

# Step 11: Frobenius Norm Calculation
frobenius_norms <- numeric(nrow(result$Qt_inv) - 1)  # Only from t=2 to T
for (t in 2:nrow(result$Qt_inv)) {
  frobenius_norms[t - 1] <- norm(result$Qt_inv[,,t] - solve(result$qbar_inv), type = "F")
}

print("Frobenius Norms of the Differences:")
print(frobenius_norms)


