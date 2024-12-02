# Load necessary libraries
library(rmgarch)
library(PerformanceAnalytics)
library(xts)

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

# Step 4: Estimate the Risk Parity Portfolio (equal risk contribution portfolio)
# Define a function to compute portfolio variance given weights and covariance matrix
portfolio_variance <- function(weights, cov_matrix) {
  return(t(weights) %*% cov_matrix %*% weights)
}

# Function to solve for risk parity weights
risk_parity_weights <- function(cov_matrix) {
  inv_vols <- 1 / sqrt(diag(cov_matrix))  # Inverse volatility for each asset
  weights <- inv_vols / sum(inv_vols)  # Normalize to sum to 1
  return(weights)
}

# Step 5: Calculate out-of-sample variance and returns
n_out_of_sample <- 250  # Number of out-of-sample observations

# Initialize vectors for out-of-sample returns and variances
out_sample_returns <- numeric(n_out_of_sample)
out_sample_variances <- numeric(n_out_of_sample)

for (i in 1:n_out_of_sample) {
  # Get the covariance matrix at time t
  cov_matrix_t <- cov_matrices[,,n - n_out_of_sample + i]
  
  # Get risk parity weights at time t
  weights_t <- risk_parity_weights(cov_matrix_t)
  
  # Calculate portfolio returns
  out_sample_returns[i] <- sum(weights_t * returns[n - n_out_of_sample + i, ])
  
  # Calculate portfolio variance
  out_sample_variances[i] <- portfolio_variance(weights_t, cov_matrix_t)
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
