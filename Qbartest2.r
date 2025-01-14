#test Qbar and average of Qts over multiple simulations
# Load necessary libraries
library(mvtnorm)   # For multivariate normal simulations
library(glasso)    # For graphical lasso

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
  glasso_fit <- glasso(cov(e), rho = 0.1)  # Apply Glasso with tuned rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 5: Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, qbar_inv, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  for (t in 2:T) {
    Qt_inv[,,t] <- (1 - b) * qbar_inv + b * (t(e[t-1, , drop=FALSE]) %*% e[t-1, , drop=FALSE]) + a * Qt_inv[,,t-1]
  }
  
  return(Qt_inv)
}

# Step 6: DCC-GARCH estimation function
dcc_garch <- function(data, a, b) {
  # Estimate GARCH parameters for each asset
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
  
  # Update correlation matrices using DCC method
  Qt_inv <- dcc_update(e, qbar_inv, a, b)
  
  return(list(Qt_inv = Qt_inv, garch_params = garch_params))
}

# Step 7: Function to simulate and average Qts over multiple runs
simulate_dcc_garch <- function(n_simulations, T, N, a, b, Sigma) {
  Qt_sum <- 0
  qbar_list <- list()
  
  for (sim in 1:n_simulations) {
    # Simulate correlated returns
    returns <- rmvnorm(T, sigma = Sigma)
    
    # Apply DCC-GARCH model
    dcc_result <- dcc_garch(returns, a, b)
    Qt_sum <- Qt_sum + apply(dcc_result$Qt_inv, c(1, 2), mean)  # Average across time for this run
    
    # Store Qbar (unconditional covariance) for comparison
    if (sim == 1) {
      qbar_list <- solve(dcc_result$Qt_inv[,,1])
    }
  }
  
  # Compute the average of Qt across all simulations
  Qt_avg <- Qt_sum / n_simulations
  
  return(list(Qt_avg = Qt_avg, Qbar_inv = qbar_list))
}

# Step 8: Frobenius norm of the difference
frobenius_norm <- function(A, B) {
  return(norm(A - B, type = "F"))
}

# Step 9: Run the simulation and compute results
set.seed(123)
n_simulations <- 100  # Number of simulations
T <- 500  # Number of time points
N <- 5  # Number of assets

# Correlation matrix (Sigma) for the variables
Sigma <- matrix(c(1, 0.8, 0.8, 0.8, 0.8,
                  0.8, 1, 0.8, 0.8, 0.8,
                  0.8, 0.8, 1, 0.8, 0.8,
                  0.8, 0.8, 0.8, 1, 0.8,
                  0.8, 0.8, 0.8, 0.8, 1), N, N)

# Parameters for DCC
a <- 0.05
b <- 0.9

# Simulate and calculate average Qt
result <- simulate_dcc_garch(n_simulations, T, N, a, b, Sigma)

# Step 10: Compare Qt_avg with Qbar_inv using Frobenius norm
frobenius_diff <- frobenius_norm(result$Qt_avg, result$Qbar_inv)

# Print results
print(paste("Frobenius norm of the difference:", frobenius_diff))

#plot for varying samples
# Load necessary libraries
library(mvtnorm)   # For multivariate normal simulations
library(glasso)    # For graphical lasso

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
  glasso_fit <- glasso(cov(e), rho = 0.1)  # Apply Glasso with tuned rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 5: Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, qbar_inv, a, b) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  for (t in 2:T) {
    Qt_inv[,,t] <- (1 - b) * qbar_inv + b * (t(e[t-1, , drop=FALSE]) %*% e[t-1, , drop=FALSE]) + a * Qt_inv[,,t-1]
  }
  
  return(Qt_inv)
}

# Step 6: DCC-GARCH estimation function
dcc_garch <- function(data, a, b) {
  # Estimate GARCH parameters for each asset
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
  
  # Update correlation matrices using DCC method
  Qt_inv <- dcc_update(e, qbar_inv, a, b)
  
  return(list(Qt_inv = Qt_inv, garch_params = garch_params))
}

# Step 7: Function to simulate and average Qts over multiple runs
simulate_dcc_garch <- function(n_simulations, T, N, a, b, Sigma) {
  Qt_sum <- 0
  qbar_list <- list()
  
  for (sim in 1:n_simulations) {
    # Simulate correlated returns
    returns <- rmvnorm(T, sigma = Sigma)
    
    # Apply DCC-GARCH model
    dcc_result <- dcc_garch(returns, a, b)
    Qt_sum <- Qt_sum + apply(dcc_result$Qt_inv, c(1, 2), mean)  # Average across time for this run
    
    # Store Qbar (unconditional covariance) for comparison
    if (sim == 1) {
      qbar_list <- solve(dcc_result$Qt_inv[,,1])
    }
  }
  
  # Compute the average of Qt across all simulations
  Qt_avg <- Qt_sum / n_simulations
  
  return(list(Qt_avg = Qt_avg, Qbar_inv = qbar_list))
}

# Step 8: Frobenius norm of the difference
frobenius_norm <- function(A, B) {
  return(norm(A - B, type = "F"))
}

# Step 9: Vary sample sizes and compute Frobenius norm differences
sample_sizes <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200)  # Different sample sizes
n_simulations <- 100  # Number of simulations
N <- 5  # Number of assets

# Correlation matrix (Sigma) for the variables
Sigma <- matrix(c(1, 0.8, 0.8, 0.8, 0.8,
                  0.8, 1, 0.8, 0.8, 0.8,
                  0.8, 0.8, 1, 0.8, 0.8,
                  0.8, 0.8, 0.8, 1, 0.8,
                  0.8, 0.8, 0.8, 0.8, 1), N, N)

# Parameters for DCC
a <- 0.05
b <- 0.9

# Store Frobenius norms for each sample size
frobenius_diffs <- numeric(length(sample_sizes))

for (i in seq_along(sample_sizes)) {
  T <- sample_sizes[i]
  
  # Simulate and calculate average Qt for the given sample size
  result <- simulate_dcc_garch(n_simulations, T, N, a, b, Sigma)
  
  # Calculate Frobenius norm difference
  frobenius_diffs[i] <- frobenius_norm(result$Qt_avg, result$Qbar_inv)
}

# Step 10: Plot Frobenius norm vs sample size
plot(sample_sizes, frobenius_diffs, type = "b", col = "blue", pch = 19,
     xlab = "Sample Size", ylab = "Frobenius Norm of Difference",
     main = "Frobenius Norm vs Sample Size")


#test for inverse Qts and Qbar (varying samples)
# Load necessary libraries
library(mvtnorm)   # For multivariate normal simulations
library(xts)       # For time series data manipulation
library(quantmod)  # For downloading stock data
library(glasso)    # For graphical lasso
library(MASS)      # For pseudo-inverse (ginv)
library(ggplot2)   # For plotting

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
dcc_update <- function(e, sigma2, qbar_inv, a, b, reg_lambda = 1e-5) {
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
  dcc_result <- dcc_update(e, sigma2, qbar_inv, a, b)
  
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
  dcc_result <- dcc_update(e, sigma2, qbar_inv, a, b)
  
  return(list(sigma2 = sigma2, correlation = dcc_result$Rt, conditional_cov = dcc_result$conditional_cov, dcc_params = dcc_params, garch_params = garch_params, qbar_inv = qbar_inv))
}

# Step 9: Calculate Frobenius norm and plot results
run_frobenius_experiment <- function(sample_sizes, n_simulations, N) {
  frob_norms <- numeric(length(sample_sizes))
  
  for (i in seq_along(sample_sizes)) {
    T <- sample_sizes[i]
    Qt_sum <- 0
    
    for (sim in 1:n_simulations) {
      mu <- rep(0, N)
      Sigma <- diag(1, N)  # Identity matrix as covariance
      diag(Sigma) <- 1  # Diagonal elements are variances
      returns <- mvtnorm::rmvnorm(T, mean = mu, sigma = Sigma)
      
      # Fit DCC-GARCH model
      dcc_fit <- dcc_garch(returns)
      Qt_sum <- Qt_sum + dcc_fit$qbar_inv
    }
    
    # Calculate average inverse covariance matrix
    Q_avg_inv <- Qt_sum / n_simulations
    
    # Calculate the difference from Q_bar and Frobenius norm
    frob_norms[i] <- norm(Q_avg_inv - dcc_fit$qbar_inv, type = "F")
  }
  
  # Create a data frame for plotting
  results <- data.frame(sample_size = sample_sizes, frobenius_norm = frob_norms)
  
  # Plot Frobenius norm vs sample size
  ggplot(results, aes(x = sample_size, y = frobenius_norm)) +
    geom_line() +
    labs(title = "Frobenius Norm of Average Inverse Covariance Matrix vs Sample Size",
         x = "Sample Size",
         y = "Frobenius Norm") +
    theme_minimal()
}

# Run the experiment with specified parameters
sample_sizes <- seq(50, 500, by = 50)  # Sample sizes from 50 to 500
n_simulations <- 100  # Number of simulations
N <- 5  # Number of assets

run_frobenius_experiment(sample_sizes, n_simulations, N)

#inverse comparison with cumulative info adapted dcc garch
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
    Qt_inv[,,t] <- (1 - b) * qbar_inv + b * (t(Qt_inv[,,t-1])) + a * (t(Qt_inv[,,t-1]))
    Rt[,,t] <- diag(1 / sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1 / sqrt(diag(solve(Qt_inv[,,t]))))
    
    # Store the conditional covariance matrix
    conditional_cov[,,t] <- diag(sqrt(sigma2[t, ])) %*% Rt[,,t] %*% diag(sqrt(sigma2[t, ]))
  }
  
  return(list(Rt = Rt, Qt_inv = Qt_inv))
}

# Step 6: Simulate data, run DCC-GARCH, and calculate Frobenius norm vs samples
simulate_dcc_garch <- function(num_simulations, num_samples_list, num_vars) {
  frobenius_norms <- numeric(length(num_samples_list))
  
  for (s in 1:length(num_samples_list)) {
    num_samples <- num_samples_list[s]
    Qt_inv_list <- list()
    
    for (sim in 1:num_simulations) {
      # Simulate multivariate normal data
      data <- rmvnorm(num_samples, mean = rep(0, num_vars), sigma = diag(num_vars))
      
      # Fit GARCH(1,1) to each time series and estimate Qt_inv
      garch_params <- t(apply(data, 2, estimate_garch11))
      sigma2 <- matrix(0, nrow = num_samples, ncol = num_vars)
      e <- matrix(0, nrow = num_samples, ncol = num_vars)
      
      for (i in 1:num_vars) {
        params <- garch_params[i, ]
        sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
        e[,i] <- data[,i] / sqrt(sigma2[,i])
      }
      
      # Estimate qbar_inv using the Glasso method
      qbar_inv <- estimate_inverse_covariance(e)
      
      # Run DCC-GARCH
      dcc_result <- dcc_update(e, qbar_inv, a = 0.01, b = 0.95, sigma2)
      
      # Store the inverse Qt for the final time point
      Qt_inv_list[[sim]] <- dcc_result$Qt_inv[,,num_samples]
    }
    
    # Calculate average inverse Qt
    avg_Qt_inv <- Reduce("+", Qt_inv_list) / num_simulations
    
    # Compute Frobenius norm between average inverse Qt and inverse qbar
    frobenius_diff <- norm(avg_Qt_inv - qbar_inv, type = "F")
    frobenius_norms[s] <- frobenius_diff
  }
  
  # Plot Frobenius norm vs samples
  plot(num_samples_list, frobenius_norms, type = "b", col = "blue", pch = 19, 
       xlab = "Number of Samples", ylab = "Frobenius Norm of Difference",
       main = "Frobenius Norm vs Number of Samples")
}

# Step 7: Run the simulation
set.seed(123)  # For reproducibility
num_simulations <- 50
num_samples_list <- seq(50, 500, 50)  # Varying sample sizes
num_vars <- 100  # Number of variables

simulate_dcc_garch(num_simulations, num_samples_list, num_vars)

#correct version

# Required libraries
library(xts)       # For time series data manipulation
library(quantmod)  # For downloading stock data
library(glasso)    # For graphical lasso
library(ggplot2)   # For plotting
library(rmgarch)   # For DCC-GARCH simulation
library(MASS)      # For pseudo-inverse (ginv)

# Step 1: Set up a DCC-GARCH specification using the 'rmgarch' package
dcc_garch_spec <- function(N) {
  # Univariate GARCH specification (GARCH(1,1) for each asset)
  univariate_spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                                variance.model = list(garchOrder = c(1, 1)),
                                distribution.model = "norm")
  
  # DCC specification with the same univariate GARCH for each asset
  dcc_spec <- dccspec(uspec = multispec(replicate(N, univariate_spec)), 
                      dccOrder = c(1, 1), 
                      distribution = "mvnorm")
  return(dcc_spec)
}

# Step 2: Simulate returns using DCC-GARCH
dcc_simulate_returns <- function(dcc_spec, T, N) {
  # Simulate from DCC-GARCH model
  dcc_fit <- dccfit(dcc_spec, data = NULL)  # No input data needed for simulation
  sim <- dccsim(dcc_fit, n.sim = T, m.sim = 1)
  returns <- fitted(sim)  # Get the simulated returns
  return(returns)
}

# Step 3: Estimate the inverse covariance matrix using Graphical Lasso (Glasso)
estimate_inverse_covariance <- function(e) {
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 4: Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, sigma2, qbar_inv, a, b, reg_lambda = 1e-5) {
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

# Step 5: DCC-GARCH estimation function
dcc_garch <- function(data) {
  # Fit GARCH(1,1) to each time series
  garch_fits <- lapply(1:ncol(data), function(i) estimate_garch11_rugarch(data[, i]))
  sigma2 <- do.call(cbind, lapply(garch_fits, function(fit) fit$sigma2))
  e <- do.call(cbind, lapply(garch_fits, function(fit) fit$residuals / sqrt(fit$sigma2)))
  
  # Estimate the inverse covariance matrix
  qbar_inv <- estimate_inverse_covariance(e)
  
  # DCC parameters
  a <- 0.01
  b <- 0.95
  
  # Update correlation matrix using DCC method
  dcc_result <- dcc_update(e, sigma2, qbar_inv, a, b)
  
  return(list(sigma2 = sigma2, correlation = dcc_result$Rt, conditional_cov = dcc_result$conditional_cov, qbar_inv = qbar_inv))
}

# Step 6: Frobenius norm calculation and experiment with DCC-GARCH simulation
run_frobenius_experiment_dcc <- function(sample_size, n_simulation_range, N) {
  frob_norms <- numeric(length(n_simulation_range))
  
  # Set up DCC-GARCH specification
  dcc_spec <- dcc_garch_spec(N)
  
  for (i in seq_along(n_simulation_range)) {
    n_simulations <- n_simulation_range[i]
    Qt_sum <- 0
    
    for (sim in 1:n_simulations) {
      # Simulate returns using DCC-GARCH process
      returns <- dcc_simulate_returns(dcc_spec, sample_size, N)
      
      # Fit DCC-GARCH model on the simulated returns
      dcc_fit <- dcc_garch(returns)
      Qt_sum <- Qt_sum + dcc_fit$qbar_inv
    }
    
    # Calculate average inverse covariance matrix
    Q_avg_inv <- Qt_sum / n_simulations
    
    # Calculate the difference from the initial qbar_inv and Frobenius norm
    frob_norms[i] <- norm(Q_avg_inv - dcc_fit$qbar_inv, type = "F")
  }
  
  # Create a data frame for plotting
  results <- data.frame(n_simulations = n_simulation_range, frobenius_norm = frob_norms)
  
  # Plot Frobenius norm vs number of simulations
  ggplot(results, aes(x = n_simulations, y = frobenius_norm)) +
    geom_line(color = "blue") +
    geom_point(color = "red") +
    labs(title = "Frobenius Norm vs Number of Simulations", 
         x = "Number of Simulations", 
         y = "Frobenius Norm") +
    theme_minimal()
}

# Step 7: Run the experiment with fixed sample size and varying number of simulations
sample_size <- 100  # Fixed sample size
n_simulation_range <- seq(10, 100, by = 10)  # Varying number of simulations
N <- 3  # Number of assets

run_frobenius_experiment_dcc(sample_size, n_simulation_range, N)

#correct version
# Required libraries
library(xts)       # For time series data manipulation
library(quantmod)  # For downloading stock data
library(glasso)    # For graphical lasso
library(ggplot2)   # For plotting
library(rmgarch)   # For DCC-GARCH simulation
library(MASS)      # For pseudo-inverse (ginv)

# Step 1: Set up a DCC-GARCH specification using the 'rmgarch' package
dcc_garch_spec <- function(N) {
  # Univariate GARCH specification (GARCH(1,1) for each asset)
  univariate_spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                                variance.model = list(garchOrder = c(1, 1)),
                                distribution.model = "norm")
  
  # DCC specification with the same univariate GARCH for each asset
  dcc_spec <- dccspec(uspec = multispec(replicate(N, univariate_spec)), 
                      dccOrder = c(1, 1), 
                      distribution = "mvnorm")
  return(dcc_spec)
}

# Step 2: Simulate returns using DCC-GARCH
dcc_simulate_returns <- function(dcc_spec, T, N) {
  # Simulate from DCC-GARCH model
  dcc_fit <- dccfit(dcc_spec, data = NULL)  # No input data needed for simulation
  sim <- dccsim(dcc_fit, n.sim = T, m.sim = 1)
  returns <- fitted(sim)  # Get the simulated returns
  return(returns)
}

# Step 3: Estimate the inverse covariance matrix using Graphical Lasso (Glasso)
estimate_inverse_covariance <- function(e) {
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 4: Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, sigma2, qbar_inv, a, b, reg_lambda = 1e-5) {
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

# Step 5: DCC-GARCH estimation function
dcc_garch <- function(data) {
  # Fit GARCH(1,1) to each time series
  garch_fits <- lapply(1:ncol(data), function(i) estimate_garch11_rugarch(data[, i]))
  sigma2 <- do.call(cbind, lapply(garch_fits, function(fit) fit$sigma2))
  e <- do.call(cbind, lapply(garch_fits, function(fit) fit$residuals / sqrt(fit$sigma2)))
  
  # Estimate the inverse covariance matrix
  qbar_inv <- estimate_inverse_covariance(e)
  
  # DCC parameters
  a <- 0.01
  b <- 0.95
  
  # Update correlation matrix using DCC method
  dcc_result <- dcc_update(e, sigma2, qbar_inv, a, b)
  
  return(list(sigma2 = sigma2, correlation = dcc_result$Rt, conditional_cov = dcc_result$conditional_cov, qbar_inv = qbar_inv))
}

# Step 6: Frobenius norm calculation and experiment with DCC-GARCH simulation
run_frobenius_experiment_dcc <- function(sample_size, n_simulation_range, N) {
  frob_norms <- numeric(length(n_simulation_range))
  
  # Set up DCC-GARCH specification
  dcc_spec <- dcc_garch_spec(N)
  
  for (i in seq_along(n_simulation_range)) {
    n_simulations <- n_simulation_range[i]
    Qt_sum <- 0
    
    for (sim in 1:n_simulations) {
      # Simulate returns using DCC-GARCH process
      returns <- dcc_simulate_returns(dcc_spec, sample_size, N)
      
      # Fit DCC-GARCH model on the simulated returns
      dcc_fit <- dcc_garch(returns)
      Qt_sum <- Qt_sum + dcc_fit$qbar_inv
    }
    
    # Calculate average inverse covariance matrix
    Q_avg_inv <- Qt_sum / n_simulations
    
    # Calculate the difference from the initial qbar_inv and Frobenius norm
    frob_norms[i] <- norm(Q_avg_inv - dcc_fit$qbar_inv, type = "F")
  }
  
  # Create a data frame for plotting
  results <- data.frame(n_simulations = n_simulation_range, frobenius_norm = frob_norms)
  
  # Plot Frobenius norm vs number of simulations
  ggplot(results, aes(x = n_simulations, y = frobenius_norm)) +
    geom_line(color = "blue") +
    geom_point(color = "red") +
    labs(title = "Frobenius Norm vs Number of Simulations", 
         x = "Number of Simulations", 
         y = "Frobenius Norm") +
    theme_minimal()
}

# Step 7: Run the experiment with fixed sample size and varying number of simulations
sample_size <- 100  # Fixed sample size
n_simulation_range <- seq(10, 100, by = 10)  # Varying number of simulations
N <- 3  # Number of assets

run_frobenius_experiment_dcc(sample_size, n_simulation_range, N)


#comparing glasso and greedyprune
# Required libraries
library(xts)       # For time series data manipulation
library(quantmod)  # For downloading stock data
library(glasso)    # For graphical lasso
 # For GreedyPrune inverse covariance estimation
library(ggplot2)   # For plotting
library(rmgarch)   # For DCC-GARCH simulation
library(MASS)      # For pseudo-inverse (ginv)

# Step 1: Set up a DCC-GARCH specification using the 'rmgarch' package
dcc_garch_spec <- function(N) {
  # Univariate GARCH specification (GARCH(1,1) for each asset)
  univariate_spec <- ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                                variance.model = list(garchOrder = c(1, 1)),
                                distribution.model = "norm")
  
  # DCC specification with the same univariate GARCH for each asset
  dcc_spec <- dccspec(uspec = multispec(replicate(N, univariate_spec)), 
                      dccOrder = c(1, 1), 
                      distribution = "mvnorm")
  return(dcc_spec)
}

# Step 2: Simulate returns using DCC-GARCH
dcc_simulate_returns <- function(dcc_spec, T, N) {
  # Simulate from DCC-GARCH model
  dcc_fit <- dccfit(dcc_spec, data = NULL)  # No input data needed for simulation
  sim <- dccsim(dcc_fit, n.sim = T, m.sim = 1)
  returns <- fitted(sim)  # Get the simulated returns
  return(returns)
}

# Step 3: Estimate the inverse covariance matrix using Graphical Lasso (Glasso)
estimate_inverse_covariance_glasso <- function(e) {
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  return(glasso_fit$wi)  # Inverse covariance matrix
}

# Step 4: Estimate the inverse covariance matrix using GreedyPrune
estimate_inverse_covariance_greedyprune <- function(e) {
  greedyprune_fit <- GreedyPrune(t(e), lambda = 0.1, standardize = TRUE)  # You may need to tune lambda
  return(greedyprune_fit$Omega)  # Inverse covariance matrix
}

# Step 5: Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, sigma2, qbar_inv, a, b, reg_lambda = 1e-5) {
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

# Step 6: DCC-GARCH estimation function
dcc_garch <- function(data, method = "glasso") {
  # Fit GARCH(1,1) to each time series
  garch_fits <- lapply(1:ncol(data), function(i) estimate_garch11_rugarch(data[, i]))
  sigma2 <- do.call(cbind, lapply(garch_fits, function(fit) fit$sigma2))
  e <- do.call(cbind, lapply(garch_fits, function(fit) fit$residuals / sqrt(fit$sigma2)))
  
  # Estimate the inverse covariance matrix based on method
  if (method == "glasso") {
    qbar_inv <- estimate_inverse_covariance_glasso(e)
  } else if (method == "greedyprune") {
    qbar_inv <- estimate_inverse_covariance_greedyprune(e)
  }
  
  # DCC parameters
  a <- 0.01
  b <- 0.95
  
  # Update correlation matrix using DCC method
  dcc_result <- dcc_update(e, sigma2, qbar_inv, a, b)
  
  return(list(sigma2 = sigma2, correlation = dcc_result$Rt, conditional_cov = dcc_result$conditional_cov, qbar_inv = qbar_inv))
}

# Step 7: Frobenius norm calculation and experiment with DCC-GARCH simulation for both methods
run_frobenius_experiment_dcc <- function(sample_size, n_simulation_range, N) {
  frob_norms_glasso <- numeric(length(n_simulation_range))
  frob_norms_greedyprune <- numeric(length(n_simulation_range))
  
  # Set up DCC-GARCH specification
  dcc_spec <- dcc_garch_spec(N)
  
  for (i in seq_along(n_simulation_range)) {
    n_simulations <- n_simulation_range[i]
    Qt_sum_glasso <- 0
    Qt_sum_greedyprune <- 0
    
    for (sim in 1:n_simulations) {
      # Simulate returns using DCC-GARCH process
      returns <- dcc_simulate_returns(dcc_spec, sample_size, N)
      
      # Fit DCC-GARCH model using Glasso
      dcc_fit_glasso <- dcc_garch(returns, method = "glasso")
      Qt_sum_glasso <- Qt_sum_glasso + dcc_fit_glasso$qbar_inv
      
      # Fit DCC-GARCH model using GreedyPrune
      dcc_fit_greedyprune <- dcc_garch(returns, method = "greedyprune")
      Qt_sum_greedyprune <- Qt_sum_greedyprune + dcc_fit_greedyprune$qbar_inv
    }
    
    # Calculate average inverse covariance matrix for both methods
    Q_avg_inv_glasso <- Qt_sum_glasso / n_simulations
    Q_avg_inv_greedyprune <- Qt_sum_greedyprune / n_simulations
    
    # Calculate the difference from the initial qbar_inv and Frobenius norm
    frob_norms_glasso[i] <- norm(Q_avg_inv_glasso - dcc_fit_glasso$qbar_inv, type = "F")
    frob_norms_greedyprune[i] <- norm(Q_avg_inv_greedyprune - dcc_fit_greedyprune$qbar_inv, type = "F")
  }
  
  # Create a data frame for plotting
  results <- data.frame(n_simulations = n_simulation_range, 
                        frobenius_norm_glasso = frob_norms_glasso,
                        frobenius_norm_greedyprune = frob_norms_greedyprune)
  
  # Plot the results
  plot <- ggplot(results, aes(x = n_simulations)) + 
    geom_line(aes(y = frobenius_norm_glasso, color = "Glasso")) + 
    geom_line(aes(y = frobenius_norm_greedyprune, color = "GreedyPrune")) + 
    labs(x = "Number of Simulations", y = "Frobenius Norm Error", 
         title = "Frobenius Norm Error vs. Number of Simulations") +
    scale_color_manual(values = c("Glasso" = "blue", "GreedyPrune" = "red"), name = "Method") +
    theme_minimal()
  
  return(plot)
}

# Example usage: running the experiment with sample size = 100, N = 50 (number of assets), 
# and n_simulation_range from 10 to 100
sample_size <- 100
n_simulation_range <- seq(10, 100, by = 10)
N <- 50

# Run the experiment and plot the results
plot_result <- run_frobenius_experiment_dcc(sample_size, n_simulation_range, N)
print(plot_result)
