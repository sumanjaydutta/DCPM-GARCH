#test for qbar using average of Qt over different times

#using multi period info
# Load necessary libraries 
library(mvtnorm)   # For multivariate normal simulations
library(xts)       # For time series data manipulation
library(quantmod)  # For downloading stock data
library(glasso)    # For graphical lasso
library(GreedyPrune) # For GreedyPrune (replacing Clime)
library(MASS)      # For pseudo-inverse (ginv)
library(Matrix)    # For Frobenius norm
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

# Step 4: Estimate the inverse covariance matrix using GreedyPrune and Glasso
estimate_inverse_covariance <- function(e) {
  # Estimate with graphical lasso
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  inv_glasso <- glasso_fit$wi  # Inverse covariance matrix from Glasso
  
  # Estimate with GreedyPrune
  greedyprune_fit <- GreedyPrune(e)  # Assuming default parameters work well
  inv_greedyprune <- greedyprune_fit$Omega  # Extract inverse covariance matrix from GreedyPrune
  
  return(list(inv_glasso = inv_glasso, inv_greedyprune = inv_greedyprune))  # Return both
}

# Step 5: Dynamic Conditional Correlation (DCC) update function using cumulative information
dcc_update <- function(e, qbar_inv, a, b, reg_lambda = 1e-5) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  # Store conditional covariance matrices
  conditional_cov <- array(0, dim = c(k, k, T))
  
  for (t in 2:T) {
    # Use cumulative information from time 0 to t-1
    e_t_minus_1 <- matrix(apply(e[1:(t-1), ], 2, sum), nrow = 1)
    
    # Calculate inverse covariance for t-1 period using cumulative info
    inv_cov_t_minus_1 <- tryCatch({
      solve(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k))
    }, error = function(err) {
      return(ginv(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k)))  # Use pseudo-inverse
    })
    
    # Update Qt_inv with cumulative info
    Qt_inv[,,t] <- (1 - b) * qbar_inv + b * inv_cov_t_minus_1 + a * Qt_inv[,,t-1]
    
    # Calculate Rt from Qt_inv
    Rt[,,t] <- diag(1 / sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1 / sqrt(diag(solve(Qt_inv[,,t]))))
    
    # Store the conditional covariance matrix
    conditional_cov[,,t] <- Rt[,,t]
  }
  
  return(Qt_inv)
}

# Step 6: Frobenius norm calculation
frobenius_norm <- function(matrix1, matrix2) {
  return(norm(as.matrix(matrix1 - matrix2), type = "F"))
}

# Step 7: Run experiment for varying sample sizes and plot Frobenius norm difference
run_experiment <- function(sample_sizes, k, a, b, n_trials = 10) {
  frobenius_diff_glasso <- numeric(length(sample_sizes))
  frobenius_diff_greedyprune <- numeric(length(sample_sizes))
  
  for (i in seq_along(sample_sizes)) {
    sample_size <- sample_sizes[i]
    avg_qt_inv_glasso <- 0
    avg_qt_inv_greedyprune <- 0
    
    for (trial in 1:n_trials) {
      # Generate synthetic data
      e <- matrix(rnorm(sample_size * k), ncol = k)
      
      # Estimate the inverse covariance matrices using both methods
      inv_covs <- estimate_inverse_covariance(e)
      qbar_inv_glasso <- inv_covs$inv_glasso
      qbar_inv_greedyprune <- inv_covs$inv_greedyprune
      
      # Update Qt matrices using DCC
      qt_inv_glasso <- dcc_update(e, qbar_inv_glasso, a, b)
      qt_inv_greedyprune <- dcc_update(e, qbar_inv_greedyprune, a, b)
      
      # Average Qt_inv across time
      avg_qt_inv_glasso <- avg_qt_inv_glasso + apply(qt_inv_glasso, c(1, 2), mean)
      avg_qt_inv_greedyprune <- avg_qt_inv_greedyprune + apply(qt_inv_greedyprune, c(1, 2), mean)
    }
    
    # Take average across trials
    avg_qt_inv_glasso <- avg_qt_inv_glasso / n_trials
    avg_qt_inv_greedyprune <- avg_qt_inv_greedyprune / n_trials
    
    # Compute Frobenius norm of the difference with Qbar_inv
    frobenius_diff_glasso[i] <- frobenius_norm(avg_qt_inv_glasso, inv_covs$inv_glasso)
    frobenius_diff_greedyprune[i] <- frobenius_norm(avg_qt_inv_greedyprune, inv_covs$inv_greedyprune)
  }
  
  # Plot the results
  plot_data <- data.frame(
    Sample_Size = rep(sample_sizes, 2),
    Frobenius_Norm = c(frobenius_diff_glasso, frobenius_diff_greedyprune),
    Method = factor(rep(c("Glasso", "GreedyPrune"), each = length(sample_sizes)))
  )
  
  ggplot(plot_data, aes(x = Sample_Size, y = Frobenius_Norm, color = Method)) +
    geom_line(size = 1.2) +
    labs(title = "Frobenius Norm of Difference between Avg(Qt_inv) and Qbar_inv",
         x = "Sample Size", y = "Frobenius Norm") +
    theme_minimal()
}

# Step 8: Run the experiment with varying sample sizes
sample_sizes <- seq(50, 500, by = 50)  # Example sample sizes
k <- 10  # Number of variables
a <- 0.01  # DCC parameter
b <- 0.95  # DCC parameter

run_experiment(sample_sizes, k, a, b)

	   
#qbartest for single period info:

# Load necessary libraries
library(mvtnorm)   # For multivariate normal simulations
library(xts)       # For time series data manipulation
library(quantmod)  # For downloading stock data
library(glasso)    # For graphical lasso
library(GreedyPrune) # For GreedyPrune (replacing Clime)
library(MASS)      # For pseudo-inverse (ginv)
library(Matrix)    # For Frobenius norm
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

# Step 4: Estimate the inverse covariance matrix using GreedyPrune and Glasso
estimate_inverse_covariance <- function(e) {
  # Estimate with graphical lasso
  glasso_fit <- glasso(cor(e), rho = 0.1)  # You may need to tune rho
  inv_glasso <- glasso_fit$wi  # Inverse covariance matrix from Glasso
  
  # Estimate with GreedyPrune
  greedyprune_fit <- GreedyPrune(e)  # Assuming default parameters work well
  inv_greedyprune <- greedyprune_fit$Omega  # Extract inverse covariance matrix from GreedyPrune
  
  return(list(inv_glasso = inv_glasso, inv_greedyprune = inv_greedyprune))  # Return both
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
    e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)
    inv_cov_t_minus_1 <- tryCatch({
      solve(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k))
    }, error = function(err) {
      return(ginv(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k)))  # Use pseudo-inverse
    })
    
    Qt_inv[,,t] <- (1 - b) * qbar_inv + b * inv_cov_t_minus_1 + a * Qt_inv[,,t-1]
    
    Rt[,,t] <- diag(1 / sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1 / sqrt(diag(solve(Qt_inv[,,t]))))
    
    # Store the conditional covariance matrix
    conditional_cov[,,t] <- Rt[,,t]
  }
  
  return(Qt_inv)
}

# Step 6: Frobenius norm calculation
frobenius_norm <- function(matrix1, matrix2) {
  return(norm(as.matrix(matrix1 - matrix2), type = "F"))
}

# Step 7: Run experiment for varying sample sizes and plot Frobenius norm difference
run_experiment <- function(sample_sizes, k, a, b, n_trials = 10) {
  frobenius_diff_glasso <- numeric(length(sample_sizes))
  frobenius_diff_greedyprune <- numeric(length(sample_sizes))
  
  for (i in seq_along(sample_sizes)) {
    sample_size <- sample_sizes[i]
    avg_qt_inv_glasso <- 0
    avg_qt_inv_greedyprune <- 0
    
    for (trial in 1:n_trials) {
      # Generate synthetic data
      e <- matrix(rnorm(sample_size * k), ncol = k)
      
      # Estimate the inverse covariance matrices using both methods
      inv_covs <- estimate_inverse_covariance(e)
      qbar_inv_glasso <- inv_covs$inv_glasso
      qbar_inv_greedyprune <- inv_covs$inv_greedyprune
      
      # Update Qt matrices using DCC
      qt_inv_glasso <- dcc_update(e, qbar_inv_glasso, a, b)
      qt_inv_greedyprune <- dcc_update(e, qbar_inv_greedyprune, a, b)
      
      # Average Qt_inv across time
      avg_qt_inv_glasso <- avg_qt_inv_glasso + apply(qt_inv_glasso, c(1, 2), mean)
      avg_qt_inv_greedyprune <- avg_qt_inv_greedyprune + apply(qt_inv_greedyprune, c(1, 2), mean)
    }
    
    # Take average across trials
    avg_qt_inv_glasso <- avg_qt_inv_glasso / n_trials
    avg_qt_inv_greedyprune <- avg_qt_inv_greedyprune / n_trials
    
    # Compute Frobenius norm of the difference with Qbar_inv
    frobenius_diff_glasso[i] <- frobenius_norm(avg_qt_inv_glasso, inv_covs$inv_glasso)
    frobenius_diff_greedyprune[i] <- frobenius_norm(avg_qt_inv_greedyprune, inv_covs$inv_greedyprune)
  }
  
  # Plot the results
  plot_data <- data.frame(
    Sample_Size = rep(sample_sizes, 2),
    Frobenius_Norm = c(frobenius_diff_glasso, frobenius_diff_greedyprune),
    Method = factor(rep(c("Glasso", "GreedyPrune"), each = length(sample_sizes)))
  )
  
  ggplot(plot_data, aes(x = Sample_Size, y = Frobenius_Norm, color = Method)) +
    geom_line(size = 1.2) +
    labs(title = "Frobenius Norm of Difference between Avg(Qt_inv) and Qbar_inv",
         x = "Sample Size", y = "Frobenius Norm") +
    theme_minimal()
}

# Step 8: Run the experiment with varying sample sizes
sample_sizes <- seq(50, 500, by = 50)  # Example sample sizes
k <- 10  # Number of variables
a <- 0.01  # DCC parameter
b <- 0.95  # DCC parameter

run_experiment(sample_sizes, k, a, b)


	   
	   
	   
	   
