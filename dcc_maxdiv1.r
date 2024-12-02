# Load necessary libraries
library(rmgarch)
library(PerformanceAnalytics)
library(xts)
library(quadprog)  # For quadratic optimization

# Step 1: Load or simulate the data (we'll use random data for illustration)
set.seed(123)
n <- 1000  # number of observations
assets <- 5  # number of assets

# Simulate returns for 5 assets
returns <- matrix(rnorm(n * assets, mean = 0.0005, sd = 0.01), ncol = assets)
colnames(returns) <- c("Asset1", "Asset2", "Asset3", "Asset4", "Asset5")
returns <- xts(returns, order.by = as.Date('2000-01-01') + 1:n)

# Step 2: Specify and fit the DCC-GARCH model
# Univariate GARCH specification
uspec <- ugarchspec(mean.model = list(armaOrder = c(0,0)),
                    variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
                    distribution.model = "norm")

# Multivariate DCC specification
multispec <- multispec(replicate(assets, uspec))
dccspec <- dccspec(uspec = multispec, dccOrder = c(1,1), distribution = "mvnorm")

# Fit the DCC-GARCH model
dccfit <- dccfit(dccspec, data = returns)

# Step 3: Extract the conditional covariance matrix from the DCC model
cov_matrices <- rcov(dccfit)  # Extract time-varying covariance matrices

# Step 4: Function to solve for maximum diversification portfolio weights
max_div_weights <- function(cov_matrix) {
  # Calculate volatilities (sqrt of diagonal elements of the covariance matrix)
  volatilities <- sqrt(diag(cov_matrix))
  
  # Objective: maximize diversification ratio, equivalent to minimizing portfolio variance
  # Subject to sum of weights = 1, w_i >= 0
  # Use quadratic programming to minimize the quadratic form w' * cov_matrix * w
  
  Dmat <- 2 * cov_matrix  # Quadratic term (2 times the covariance matrix)
  dvec <- rep(0, assets)  # Linear term (no linear component)
  Amat <- cbind(rep(1, assets), diag(assets))  # Constraints matrix
  bvec <- c(1, rep(0, assets))  # Constraints: sum of weights = 1 and w_i >= 0
  
  # Solve quadratic programming problem
  qp_result <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  weights <- qp_result$solution  # Extract the portfolio weights
  
  return(weights)
}

# Step 5: Calculate out-of-sample variance and returns for Maximum Diversification Portfolio
n_out_of_sample <- 250  # Number of out-of-sample observations

# Initialize vectors for out-of-sample returns and variances
out_sample_returns <- numeric(n_out_of_sample)
out_sample_variances <- numeric(n_out_of_sample)

for (i in 1:n_out_of_sample) {
  # Get the covariance matrix at time t
  cov_matrix_t <- cov_matrices[,,n - n_out_of_sample + i]
  
  # Get maximum diversification weights at time t
  weights_t <- max_div_weights(cov_matrix_t)
  
  # Calculate portfolio returns
  out_sample_returns[i] <- sum(weights_t * returns[n - n_out_of_sample + i, ])
  
  # Calculate portfolio variance
  out_sample_variances[i] <- t(weights_t) %*% cov_matrix_t %*% weights_t
}

# Convert to time series
out_sample_returns <- xts(out_sample_returns, order.by = index(returns[(n-n_out_of_sample+1):n]))
out_sample_variances <- xts(out_sample_variances, order.by = index(returns[(n-n_out_of_sample+1):n]))

# Step 6: Annualize return and variance, and calculate Information Ratio
annualized_return <- mean(out_sample_returns) * 252  # Assuming 252 trading days
annualized_variance <- mean(out_sample_variances) * 252  # Assuming 252 trading days
portfolio_sd <- sqrt(annualized_variance)
information_ratio <- annualized_return / portfolio_sd

# Print results
print(paste("Annualized Return:", round(annualized_return, 4)))
print(paste("Annualized Portfolio Standard Deviation:", round(portfolio_sd, 4)))
print(paste("Information Ratio:", round(information_ratio, 4)))
