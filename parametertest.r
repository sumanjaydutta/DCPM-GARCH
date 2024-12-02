#using multvariate norm sim
# Required Libraries
library(glasso)        # For Graphical Lasso
library(clime)         # For GreedyPrune (Clime)
library(mvtnorm)       # For simulating multivariate normal data
library(MASS)          # For pseudo-inverse (ginv)

# Step 1: GARCH(1,1) Estimation Function
estimate_garch11 <- function(data) {
  library(fGarch)
  garch_fit <- garchFit(~ garch(1,1), data = data, trace = FALSE)
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

# Step 3: Estimate the inverse covariance matrix using Graphical Lasso or GreedyPrune
estimate_inverse_covariance <- function(e, method = "glasso") {
  corr_matrix <- cor(e)
  
  if (method == "glasso") {
    glasso_fit <- glasso(corr_matrix, rho = 0.1)  # Adjust rho as needed
    return(glasso_fit$wi)  # Inverse covariance matrix
  } else if (method == "greedyprune") {
    clime_fit <- greedyprune(corr_matrix, lambda = 0.1)  # Adjust lambda as needed
    return(clime_fit$icov)  # Inverse covariance matrix
  } else {
    stop("Invalid method. Use 'glasso' or 'greedyprune'.")
  }
}

# Step 4: Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, qbar_inv, a, b, reg_lambda = 1e-5) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  for (t in 2:T) {
    e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)
    inv_cov_t_minus_1 <- solve(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k))
    
    Qt_inv[,,t] <- (1 - b) * qbar_inv + b * inv_cov_t_minus_1 + a * Qt_inv[,,t-1]
    
    if (any(eigen(Qt_inv[,,t])$values <= 0)) {
      Qt_inv[,,t] <- qbar_inv  # Use qbar_inv if not positive definite
    }
    
    Rt[,,t] <- diag(1 / sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1 / sqrt(diag(solve(Qt_inv[,,t]))))
  }
  
  return(list(Rt = Rt, Qt_inv = Qt_inv))
}

# Step 5: Log-likelihood function for DCC-GARCH
dcc_garch_loglik <- function(params, data, method = "glasso") {
  a <- params[1]
  b <- params[2]
  
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  qbar_inv <- estimate_inverse_covariance(e, method)
  
  dcc_result <- dcc_update(e, qbar_inv, a, b)
  
  loglik <- -sum(dnorm(e, mean = 0, sd = sqrt(sigma2), log = TRUE))
  
  return(loglik)
}

# Step 6: Optimize DCC parameters
estimate_dcc_params <- function(data, method = "glasso") {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1), method = method)
  
  return(optim_result$par)
}

# Step 7: DCC-GARCH estimation function
dcc_garch <- function(data, method = "glasso") {
  dcc_params <- estimate_dcc_params(data, method)
  a <- dcc_params[1]
  b <- dcc_params[2]
  
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  qbar_inv <- estimate_inverse_covariance(e, method)
  
  dcc_result <- dcc_update(e, qbar_inv, a, b)
  
  return(list(sigma2 = sigma2, correlation = dcc_result$Rt, Qt_inv = dcc_result$Qt_inv, dcc_params = dcc_params, garch_params = garch_params))
}

# Step 8: Simulate and estimate DCC parameters multiple times
simulate_dcc_params <- function(n_sim, true_cov_matrix, method = "glasso") {
  n <- nrow(true_cov_matrix)  # Number of assets
  T <- 500  # Time points
  mu <- rep(0, n)
  
  dcc_params_list <- matrix(0, nrow = n_sim, ncol = 2)  # Store a and b for each simulation
  
  for (sim in 1:n_sim) {
    # Simulate correlated returns
    returns <- rmvnorm(T, mean = mu, sigma = true_cov_matrix)
    
    # Apply the DCC-GARCH model
    result <- dcc_garch(returns, method = method)
    
    # Store the estimated DCC parameters (a and b)
    dcc_params_list[sim, ] <- result$dcc_params
  }
  
  # Compute the average and standard deviation of DCC parameters
  avg_dcc_params <- colMeans(dcc_params_list)
  sd_dcc_params <- apply(dcc_params_list, 2, sd)
  
  return(list(avg_dcc_params = avg_dcc_params, sd_dcc_params = sd_dcc_params))
}

# Step 9: Set up the simulation
set.seed(123)
true_cov_matrix <- matrix(c(1, 0.6, 0.6, 0.6, 1, 0.6, 0.6, 0.6, 1), nrow = 3)
n_sim <- 100

# Run simulations and get the average and standard deviation of DCC parameters for both methods
results_glasso <- simulate_dcc_params(n_sim, true_cov_matrix, method = "glasso")
results_greedyprune <- simulate_dcc_params(n_sim, true_cov_matrix, method = "greedyprune")

# Print results
cat("Average DCC parameters (Graphical Lasso):", results_glasso$avg_dcc_params, "\n")
cat("Standard deviation of DCC parameters (Graphical Lasso):", results_glasso$sd_dcc_params, "\n")
cat("Average DCC parameters (GreedyPrune):", results_greedyprune$avg_dcc_params, "\n")
cat("Standard deviation of DCC parameters (GreedyPrune):", results_greedyprune$sd_dcc_params, "\n")


#using dcc sim

# Required Libraries 
library(glasso)        # For Graphical Lasso
library(clime)         # For GreedyPrune (Clime)
library(fGarch)        # For GARCH modeling
library(MASS)          # For pseudo-inverse (ginv)

# Step 1: GARCH(1,1) Estimation Function
estimate_garch11 <- function(data) {
  # Estimating GARCH(1,1) parameters (mu, alpha, beta)
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

# Step 3: Estimate the inverse covariance matrix using Graphical Lasso or GreedyPrune
estimate_inverse_covariance <- function(e, method = "glasso") {
  # Correlation matrix of the residuals
  corr_matrix <- cor(e)
  
  if (method == "glasso") {
    # Apply Graphical Lasso
    glasso_fit <- glasso(corr_matrix, rho = 0.1)  # Adjust rho as needed
    return(glasso_fit$wi)  # Inverse covariance matrix
  } else if (method == "greedyprune") {
    # Apply GreedyPrune (Clime)
    clime_fit <- greedyprune(corr_matrix, lambda = 0.1)  # Adjust lambda as needed
    return(clime_fit$icov)  # Inverse covariance matrix
  } else {
    stop("Invalid method. Use 'glasso' or 'greedyprune'.")
  }
}

# Step 4: Dynamic Conditional Correlation (DCC) update function
dcc_update <- function(e, qbar_inv, a, b, reg_lambda = 1e-5) {
  T <- nrow(e)
  k <- ncol(e)
  
  Qt_inv <- array(0, dim = c(k, k, T))
  Rt <- array(0, dim = c(k, k, T))
  Qt_inv[,,1] <- qbar_inv
  
  conditional_cov <- array(0, dim = c(k, k, T))
  
  for (t in 2:T) {
    if (t > 2) {
      e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)  
      e_t_minus_2 <- matrix(e[t-2, ], nrow = 1)  
      
      inv_cov_t_minus_1 <- tryCatch({
        solve(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k))
      }, error = function(err) {
        return(ginv(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k)))  
      })
      
      inv_cov_t_minus_2 <- tryCatch({
        solve(t(e_t_minus_2) %*% e_t_minus_2 + diag(reg_lambda, k))
      }, error = function(err) {
        return(ginv(t(e_t_minus_2) %*% e_t_minus_2 + diag(reg_lambda, k)))  
      })
      
      diff_inv_cov <- inv_cov_t_minus_1 - inv_cov_t_minus_2
    } else {
      e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)  
      inv_cov_t_minus_1 <- solve(t(e_t_minus_1) %*% e_t_minus_1 + diag(reg_lambda, k))
      diff_inv_cov <- inv_cov_t_minus_1
    }
    
    Qt_inv[,,t] <- (1 - b) * qbar_inv + b * diff_inv_cov + a * Qt_inv[,,t-1]
    
    if (any(eigen(Qt_inv[,,t])$values <= 0)) {
      Qt_inv[,,t] <- qbar_inv  
      warning("Qt_inv is not positive definite. Using qbar_inv instead.")
    }
    
    Rt[,,t] <- diag(1 / sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1 / sqrt(diag(solve(Qt_inv[,,t]))))
    
    conditional_cov[,,t] <- diag(sqrt(sigma2[t, ])) %*% Rt[,,t] %*% diag(sqrt(sigma2[t, ]))
  }
  
  return(list(Rt = Rt, conditional_cov = conditional_cov))
}

# Step 5: Log-likelihood function for DCC-GARCH
dcc_garch_loglik <- function(params, data, method = "glasso") {
  a <- params[1]
  b <- params[2]
  
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  qbar_inv <- estimate_inverse_covariance(e, method)
  
  dcc_result <- dcc_update(e, qbar_inv, a, b)
  
  loglik <- -sum(dnorm(e, mean = 0, sd = sqrt(sigma2), log = TRUE))
  
  return(loglik)
}

# Step 6: Optimize DCC parameters
estimate_dcc_params <- function(data, method = "glasso") {
  start_params <- c(a = 0.01, b = 0.95)
  optim_result <- optim(start_params, dcc_garch_loglik, data = data, method = "L-BFGS-B", 
                        lower = c(0, 0), upper = c(1, 1), method = method)
  
  return(optim_result$par)
}

# Step 7: DCC-GARCH estimation function
dcc_garch <- function(data, method = "glasso") {
  dcc_params <- estimate_dcc_params(data, method)
  a <- dcc_params[1]
  b <- dcc_params[2]
  
  garch_params <- t(apply(data, 2, estimate_garch11))
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  qbar_inv <- estimate_inverse_covariance(e, method)
  
  dcc_result <- dcc_update(e, qbar_inv, a, b)
  
  return(list(sigma2 = sigma2, correlation = dcc_result$Rt, conditional_cov = dcc_result$conditional_cov, dcc_params = dcc_params, garch_params = garch_params))
}

# Step 8: Simulate DCC-GARCH and collect DCC parameters
simulate_dcc_parameters <- function(n_sim, true_cov_matrix, method = "glasso") {
  n <- nrow(true_cov_matrix)  # Number of assets
  T <- 500  # Time points
  mu <- rep(0, n)
  
  dcc_params_list <- matrix(0, n_sim, 2)  # To store DCC parameters (a and b) from each simulation
  
  for (sim in 1:n_sim) {
    # Simulate returns using DCC-GARCH
    returns <- matrix(0, T, n)
    for (i in 1:n) {
      # Generate GARCH(1,1) returns for each asset
      garch_params <- estimate_garch11(rnorm(T))
      sigma2 <- garch11_fit(rnorm(T), garch_params[1], garch_params[2], garch_params[3])
      returns[,i] <- rnorm(T, mean = 0, sd = sqrt(sigma2))
    }
    
    # Apply the DCC-GARCH model
    result <- dcc_garch(returns, method = method)
    
    # Store the DCC parameters (a, b)
    dcc_params_list[sim, ] <- result$dcc_params
  }
  
  # Calculate mean and standard deviation for DCC parameters (a and b)
  dcc_params_mean <- colMeans(dcc_params_list)
  dcc_params_sd <- apply(dcc_params_list, 2, sd)
  
  return(list(mean = dcc_params_mean, sd = dcc_params_sd))
}

# Step 9: Set parameters for simulation
true_cov_matrix <- diag(1, 3)  # Example: Identity matrix for three assets
n_sim <- 100  # Number of simulations

# Step 10: Run simulation and collect DCC parameters for 'glasso' and 'greedyprune'
glasso_results <- simulate_dcc_parameters(n_sim, true_cov_matrix, method = "glasso")
greedyprune_results <- simulate_dcc_parameters(n_sim, true_cov_matrix, method = "greedyprune")

# Step 11: Print results
cat("Glasso Method:\n")
cat("Mean DCC parameters (a, b):", glasso_results$mean, "\n")
cat("Standard deviation of DCC parameters (a, b):", glasso_results$sd, "\n")

cat("\nGreedyPrune Method:\n")
cat("Mean DCC parameters (a, b):", greedyprune_results$mean, "\n")
cat("Standard deviation of DCC parameters (a, b):", greedyprune_results$sd, "\n")
