# Load required libraries
library(xts)
library(quantmod)
library(xdcclarge)
library(quadprog)

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


# Step 1: Define Portfolio Strategies
# Minimum Variance Portfolio
min_variance_weights <- function(H) {
  k <- ncol(H)
  dvec <- rep(0, k)
  Amat <- cbind(rep(1, k), diag(k))
  bvec <- c(1, rep(0, k))
  result <- solve.QP(Dmat = H, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
  return(result$solution)
}

# Risk Parity Portfolio
risk_parity_weights <- function(H) {
  inv_vol <- 1 / sqrt(diag(H))
  return(inv_vol / sum(inv_vol))
}

# Maximum Diversification Portfolio
max_diversification_weights <- function(H) {
  avg_correlation <- colSums(H) / sqrt(diag(H) %*% t(diag(H)))
  return(1 / avg_correlation)
}

# Step 2: Conditional Covariance Estimation Using DCC
dcc_cov_estimation <- function(returns, method) {
  model <- cdcc(returns, method = method)
  return(model$Ht)
}

# Step 3: Implement Rolling Window & Portfolio Rebalancing
rolling_portfolio_variance <- function(returns, window_size, rebalance_freq, strategy_func, method) {
  n <- ncol(returns)
  T <- nrow(returns)
  
  # Define matrices to store results
  out_of_sample_variances <- numeric(T - window_size)
  
  for (t in seq(window_size, T - rebalance_freq, by = rebalance_freq)) {
    in_sample <- returns[(t - window_size + 1):t, ]
    out_of_sample <- returns[(t + 1):(t + rebalance_freq), ]
    
    # Step 4: Calculate Conditional Covariance Matrix
    Ht <- dcc_cov_estimation(in_sample, method = method)
    
    # Step 5: Portfolio Weights
    weights <- strategy_func(Ht[,,dim(Ht)[3]])
    
    # Step 6: Out-of-Sample Portfolio Variance
    out_of_sample_variances[t - window_size + 1] <- mean(rowSums((out_of_sample %*% weights)^2))
  }
  return(out_of_sample_variances)
}

# Step 7: Experiment Over Multiple Sample Sizes
portfolio_experiment <- function(returns, strategies, methods, sample_sizes, rebalance_freq) {
  results <- list()
  
  for (method in methods) {
    method_results <- list()
    
    for (sample_size in sample_sizes) {
      strategy_variances <- sapply(strategies, function(strategy) {
        rolling_variance <- rolling_portfolio_variance(returns, sample_size, rebalance_freq, strategy, method)
        return(c(mean(rolling_variance), sd(rolling_variance)))
      })
      
      method_results[[as.character(sample_size)]] <- strategy_variances
    }
    results[[method]] <- method_results
  }
  
  return(results)
}

# Step 8: Run Experiment and Plot Results
run_experiment <- function(returns, strategies, methods, sample_sizes, rebalance_freq) {
  results <- portfolio_experiment(returns, strategies, methods, sample_sizes, rebalance_freq)
  
  concentration_ratios <- ncol(returns) / sample_sizes
  avg_variances <- sapply(results, function(method_results) {
    sapply(method_results, function(strategy_variances) strategy_variances[1, ])
  })
  
  plot(concentration_ratios, avg_variances, type = "o", col = 1:length(methods),
       xlab = "Concentration Ratio (Dimension/Sample Size)", ylab = "Average Out-of-Sample Variance")
  legend("topright", legend = methods, col = 1:length(methods), lty = 1)
}

# Parameters
strategies <- list("Minimum Variance" = min_variance_weights, 
                   "Risk Parity" = risk_parity_weights, 
                   "Max Diversification" = max_diversification_weights)

methods <- c("cov", "ls", "nls")
sample_sizes <- c(50, 100, 150, 200)
rebalance_freq <- 20  # Monthly rebalancing for daily returns (approx. 20 days)

# Sample Data (Example: SP500 returns for demonstration)
data <- getSymbols("^GSPC", from = "2010-01-01", to = "2020-01-01", auto.assign = FALSE)
returns <- diff(log(Cl(data))) * 100
returns <- na.omit(returns)

# Run the Experiment
run_experiment(returns, strategies, methods, sample_sizes, rebalance_freq)

#incorporating DCCGL

# Load necessary libraries
library(mvtnorm)   
library(xts)       
library(quantmod)  
library(glasso)    
library(MASS)      

# Define all methods, including DCCGL, to estimate H_t on a rolling window basis
estimate_Ht_methods <- function(data) {
  # Run rolling estimation for each method and store results
  window_size <- 252  # Adjust based on data sample size and rebalancing frequency
  
  results <- list(
    DCCGL = list(),
    DCCNL = list(),
    # Add other methods if needed: DCCGP, DCCL1, etc.
  )
  
  for (t in (window_size + 1):nrow(data)) {
    window_data <- data[(t - window_size):(t - 1), ]
    
    # Estimate H_t for DCCGL
    dccgl_result <- dcc_garch(window_data)  # DCCGL estimation function from your code
    results$DCCGL[[t]] <- dccgl_result$conditional_cov[,,t]
    
    # Include other methods similarly (e.g., DCCNL using `cdcc` function if previously set up)
    # Example placeholder: results$DCCNL[[t]] <- dcc_nl_estimation(window_data)
  }
  
  return(results)
}

# Function to calculate out-of-sample variance for each portfolio strategy
calc_portfolio_variance <- function(Ht, weights) {
  apply(Ht, 3, function(H) {
    t(weights) %*% H %*% weights
  })
}

# Compare strategies across methods
compare_methods <- function(data, methods_results) {
  # Define portfolio strategies here: minimum variance, risk parity, maximum diversification
  portfolios <- list(
    min_variance = function(H) solve(H) %*% rep(1, ncol(H)) / sum(solve(H) %*% rep(1, ncol(H))),
    # Define risk parity and maximum diversification similarly
  )
  
  # Calculate out-of-sample variances for each strategy and method
  variances <- list()
  
  for (method in names(methods_results)) {
    Ht_method <- methods_results[[method]]
    
    variances[[method]] <- sapply(names(portfolios), function(port) {
      weights <- portfolios[[port]](Ht_method)
      calc_portfolio_variance(Ht_method, weights)
    })
  }
  
  return(variances)
}

# Apply the comparison
rolling_Ht_results <- estimate_Ht_methods(data)  # Use your dataset
comparison_results <- compare_methods(data, rolling_Ht_results)

# Plot comparison
plot_comparison <- function(results) {
  # Code to plot results; e.g., boxplot for variance comparisons across methods
  boxplot(results, main = "Out-of-Sample Portfolio Variance Comparison")
}

plot_comparison(comparison_results)

