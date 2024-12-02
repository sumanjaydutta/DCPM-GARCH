library(mvtnorm)   # For multivariate normal simulations
library(glasso)    # For graphical lasso
library(xts)       # For time series handling

# Regularize the covariance matrix to avoid singularity
regularize_cov <- function(cov_matrix, epsilon = 1e-6) {
  diag(cov_matrix) <- diag(cov_matrix) + epsilon
  return(cov_matrix)
}

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

# Step 2: Estimate parameters for each asset using numerical optimization
estimate_garch11 <- function(y) {
  start_params <- c(omega = 0.1, alpha = 0.1, beta = 0.8)
  optim_result <- optim(start_params, function(params) {
    omega <- params[1]
    alpha <- params[2]
    beta  <- params[3]
    
    sigma2 <- garch11_fit(y, omega, alpha, beta)
    
    # Return negative log-likelihood
    -sum(dnorm(y, mean = 0, sd = sqrt(sigma2), log = TRUE))
  }, y = y, method = "L-BFGS-B", lower = c(1e-6, 1e-6, 1e-6), upper = c(1, 1, 1))
  
  return(optim_result$par)
}

# Step 3: Estimate the inverse covariance matrix using Graphical Lasso (Glasso)
estimate_inverse_covariance <- function(e) {
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
  
  for (t in 2:T) {
    e_t_minus_1 <- matrix(e[t-1, ], nrow = 1)
    diff_cov <- t(e_t_minus_1) %*% e_t_minus_1
    
    Qt_inv[,,t] <- (1 - a - b) * qbar_inv + a * diff_cov + b * Qt_inv[,,t-1]
    
    # Check if Qt_inv is positive definite
    if (any(eigen(Qt_inv[,,t])$values <= 0)) {
      Qt_inv[,,t] <- qbar_inv  # Revert to qbar_inv if not positive definite
      warning("Qt_inv is not positive definite. Using qbar_inv instead.")
    }
    
    Rt[,,t] <- diag(1 / sqrt(diag(solve(Qt_inv[,,t])))) %*% solve(Qt_inv[,,t]) %*% diag(1 / sqrt(diag(solve(Qt_inv[,,t]))))
  }
  
  return(Rt)
}

# Adapted DCC-GARCH estimation function with loss function

# Step 5: DCC-GARCH estimation function (updated with conditional covariance)
dcc_garch <- function(data) {
  # Estimate DCC parameters
  start_params <- c(a = 0.01, b = 0.95)
  
  # Define the objective function to minimize (negative log-likelihood)
  objective_fn <- function(params) {
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
    
    qbar_inv <- estimate_inverse_covariance(e)
    Rt <- dcc_update(e, qbar_inv, a, b)
    
    # Compute log-likelihood (negative to minimize)
    return(-sum(dnorm(e, mean = 0, sd = sqrt(sigma2), log = TRUE)))
  }
  
  optim_result <- optim(start_params, objective_fn, method = "L-BFGS-B", lower = c(0, 0), upper = c(1, 1))
  
  dcc_params <- optim_result$par
  garch_params <- t(apply(data, 2, estimate_garch11))
  
  sigma2 <- matrix(0, nrow(data), ncol(data))
  e <- matrix(0, nrow(data), ncol(data))
  
  for (i in 1:ncol(data)) {
    params <- garch_params[i, ]
    sigma2[,i] <- garch11_fit(data[,i], params[1], params[2], params[3])
    e[,i] <- data[,i] / sqrt(sigma2[,i])
  }
  
  qbar_inv <- estimate_inverse_covariance(e)
  Rt <- dcc_update(e, qbar_inv, dcc_params[1], dcc_params[2])
  
  # Calculate the conditional covariance matrix
  conditional_cov <- array(0, dim = c(ncol(data), ncol(data), nrow(data)))
  for (t in 1:nrow(data)) {
    conditional_cov[,,t] <- diag(sqrt(sigma2[t,])) %*% Rt[,,t] %*% diag(sqrt(sigma2[t,]))
  }
  
  return(list(sigma2 = sigma2, correlation = Rt, conditional_cov = conditional_cov, dcc_params = dcc_params))
}

# Loss Function for comparison
average_loss <- function(estimated_cov, true_cov, m) {
  estimated_cov_inv <- solve(estimated_cov)
  true_cov_inv <- solve(true_cov)
  
  num <- t(m) %*% estimated_cov_inv %*% true_cov %*% estimated_cov_inv %*% m
  denom_est <- (t(m) %*% estimated_cov_inv %*% m)^2
  denom_true <- t(m) %*% true_cov_inv %*% m
  
  loss <- (num / denom_est) - (1 / denom_true)
  return(as.numeric(loss))
}

# Simulation data
set.seed(123)
n <- 2  # Number of variables
T <- 500  # Number of observations
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.8, 0.8, 1), n, n)

# Simulate correlated returns
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# True Covariance Matrix
true_cov <- Sigma

# Mean vector for the loss function
mean_returns <- colMeans(returns)

# Step 6: Run the adapted DCC-GARCH estimation
result <- dcc_garch(returns)

# Calculate the conditional covariance at the final time step
estimated_cov <- result$conditional_cov[,,T]

# Step 7: Calculate the average loss
loss <- average_loss(estimated_cov, true_cov, mean_returns)
 # Load Required Libraries
library(mvtnorm)  # For multivariate normal simulations
library(rugarch)  # For GARCH model
library(xdcclarge) # For CDCC estimation

# Step 1: Simulate Synthetic Data
set.seed(123)
n <- 3  # Number of assets
T <- 500  # Number of observations

# Generate correlated returns
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.5, 0.5,
                  0.5, 1, 0.3,
                  0.5, 0.3, 1), n, n)  # Example covariance matrix

# Simulate correlated returns
returns <- rmvnorm(T, mean = mu, sigma = Sigma)

# Calculate log returns
Rtn <- log(returns[-1,] / returns[-T,])

# Step 2: GARCH Parameter Estimation with rugarch
spec = ugarchspec()
mspec = multispec(replicate(spec, n = n))
fitlist = multifit(multispec = mspec, data = Rtn)

# Check for residuals and conditional variances
ht <- sigma(fitlist)^2  # Conditional variances
residuals <- residuals(fitlist)  # Standardized residuals

# Ensure that ht and residuals are not NULL
if (is.null(ht) || is.null(residuals)) {
  stop("ht or residuals are NULL. Check the GARCH fitting.")
}

# Step 3: Check Residuals for CDCC Estimation
if (nrow(residuals) < 1 || ncol(residuals) < 1) {
  stop("Residuals matrix is empty. Check GARCH fitting.")
}

# Step 4: Define the Average Loss Function
average_loss <- function(estimated_cov, true_cov, m) {
  estimated_cov_inv <- solve(estimated_cov)
  true_cov_inv <- solve(true_cov)
  
  num <- t(m) %*% estimated_cov_inv %*% true_cov %*% estimated_cov_inv %*% m
  denom_est <- (t(m) %*% estimated_cov_inv %*% m)^2
  denom_true <- t(m) %*% true_cov_inv %*% m
  
  loss <- (num / denom_est) - (1 / denom_true)
  return(as.numeric(loss))
}

# True Covariance Matrix
true_cov <- Sigma

# Mean vector for the loss function
mean_returns <- colMeans(Rtn)

# Step 5: Calculate Average Loss for Different Methods
methods <- c("COV", "LS", "NLS")
loss_results <- list()

for (method in methods) {
  cat("\nRunning CDCC estimation for method:", method, "\n")
  
  # CDCC Estimation using the xdcc package
  cDCC <- cdcc_estimation(ini.para = c(0.05, 0.93), ht = ht, residuals = residuals, method = method)
  
  # Check if the conditional covariance is returned correctly
  if (is.null(cDCC$cdcc_Cov)) {
    cat("Warning: cdcc_Cov is NULL for method", method, ". Continuing to next method.\n")
    next
  }
  
  # Extract the conditional covariance matrix at the final time step
  estimated_cov <- cDCC$cdcc_Cov[, , T]
  
  # Calculate the average loss
  loss <- average_loss(estimated_cov, true_cov, mean_returns)
  
  # Store the loss for this method
  loss_results[[method]] <- loss
  
  # Print the average loss for this method
  cat("Average loss for method", method, ":", loss, "\n")
}

# Output the loss results
cat("\nFinal Loss Results:\n")
print(loss_results)

# Output the results
cat("Conditional variances:\n")
print(result$sigma2)

cat("Conditional covariance matrix at time T:\n")
print(estimated_cov)

cat("Estimated DCC parameters:\n")
print(result$dcc_params)

cat("Average loss:\n")
print(loss)
