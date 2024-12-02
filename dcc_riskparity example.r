#use conditional covariance matrix estimate from DCC-GARCH funcion from rmgrach package for the risk-prity portfolio
library(rmgarch)
library(quadprog)
library(mvtnorm)

# Simulate correlated returns for five assets
set.seed(123)
n <- 5  # Number of assets
T <- 500  # Time points
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.8, 0.8, 0.8, 0.8,
                  0.8, 1, 0.8, 0.8, 0.8,
                  0.8, 0.8, 1, 0.8, 0.8,
                  0.8, 0.8, 0.8, 1, 0.8,
                  0.8, 0.8, 0.8, 0.8, 1), n, n)

# Simulate correlated returns
returns <- mvtnorm::rmvnorm(T, mean = mu, sigma = Sigma)

# Fit DCC-GARCH Model
spec <- dccspec(uspec = multispec(replicate(n, ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                                            mean.model = list(armaOrder = c(0, 0))))),
                 dccOrder = c(1, 1), distribution = "mvnorm")

dcc_fit <- dccfit(spec, data = returns)

# Get the conditional covariance matrix
conditional_cov_matrix <- rcov(dcc_fit)

# Function to calculate risk contributions
risk_contribution <- function(weights, cov_matrix) {
  portfolio_variance <- as.numeric(t(weights) %*% cov_matrix %*% weights)
  
  # Marginal Risk Contribution
  mrc <- cov_matrix %*% weights
  risk_contributions <- (weights * mrc) / sqrt(portfolio_variance)
  return(risk_contributions)
}

# Function to optimize for risk parity portfolio
risk_parity_optimization <- function(cov_matrix) {
  n <- ncol(cov_matrix)
  initial_weights <- rep(1/n, n)
  
  # Objective function to minimize
  objective <- function(weights) {
    risk_contributions <- risk_contribution(weights, cov_matrix)
    diff <- risk_contributions - mean(risk_contributions)
    return(sum(diff^2))
  }
  
  # Optimize weights
  result <- optim(initial_weights, objective, method = "L-BFGS-B",
                  lower = rep(0, n), upper = rep(1, n))
  
  # Normalize weights to sum to 1
  optimal_weights <- result$par / sum(result$par)
  return(optimal_weights)
}

# Optimize for risk parity weights using the covariance matrix at the last time point
risk_parity_weights <- risk_parity_optimization(conditional_cov_matrix[,,T])
print("Risk Parity Weights:")
print(risk_parity_weights)

# Print conditional covariance matrix at the last time point
print("Conditional Covariance Matrix at Time T:")
print(conditional_cov_matrix[,,T])
