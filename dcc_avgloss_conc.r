# Load Required Libraries
library(xts)
library(mvtnorm)
library(quantmod)
library(glasso)
library(MASS)
library(xdccgarch)  # For CDCC-GARCH and NL-DCC

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

# Step 9: CDCC-GARCH estimation function
cdcc_garch <- function(data) {
  # Use xdccgarch's CDCC-GARCH function
  result <- xdccgarch::cdcc_garch(data, order = c(1, 1), dcc_order = c(1, 1))
  return(result)
}

# Step 10: Nonlinear DCC estimation function
nl_dcc <- function(data) {
  # Use xdccgarch's nonlinear DCC function
  result <- xdccgarch::nl_dcc(data)
  return(result)
}

# Step 11: Compare all methods for different sample sizes
compare_methods <- function(sample_size, n_assets) {
  set.seed(123)
  mu <- rep(0, n_assets)
  Sigma <- matrix(0.8, n_assets, n_assets) + diag(n_assets) * 0.2  # Example covariance matrix
  returns <- rmvnorm(sample_size, mean = mu, sigma = Sigma)

  # Estimate conditional covariance with Adapted DCC-GARCH
  dcc_result <- dcc_garch(returns)
  
  # Estimate conditional covariance with CDCC-GARCH
  cdcc_result <- cdcc_garch(returns)
  
  # Estimate conditional covariance with Nonlinear DCC
  nl_result <- nl_dcc(returns)
  
  # Collect results
  comparison_results <- list(
    adapted_dcc = dcc_result$conditional_cov[,,sample_size],
    cdcc = cdcc_result$conditional_cov,
    nl_dcc = nl_result$conditional_cov
  )
  
  return(comparison_results)
}

# Step 12: Loop over different sample sizes and collect results
sample_sizes <- c(100, 200, 300, 400, 500)  # Adjust sample sizes as needed
n_assets <- 5  # Keep fixed dimension

comparison_all <- lapply(sample_sizes, function(size) compare_methods(size, n_assets))

# Print out the conditional covariance matrices for each method
for (i in seq_along(sample_sizes)) {
  cat("Sample Size:", sample_sizes[i], "\n")
  cat("Adapted DCC Conditional Covariance:\n")
  print(comparison_all[[i]]$adapted_dcc)
  cat("CDCC-GARCH Conditional Covariance:\n")
  print(comparison_all[[i]]$cdcc)
  cat("NL-DCC Conditional Covariance:\n")
  print(comparison_all[[i]]$nl_dcc)
  cat("\n----------------------------------\n")
}
