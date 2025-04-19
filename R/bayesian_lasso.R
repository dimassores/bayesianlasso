#' Bayesian LASSO Implementation
#' 
#' Functions for implementing the Bayesian LASSO method through Gibbs sampling
#' as presented by Park and Casella (2008).
#' 
#' @name bayesian_lasso
#' @docType package
NULL

#' Calculate the Mode of a Vector
#'
#' @param x A vector (numeric or character)
#' @return The mode (most frequent value) of the input vector
#' @export
#' @examples
#' Mode(c(1, 2, 2, 3, 4))
#' Mode(c("a", "b", "b", "c"))
Mode <- function(x) {
  if (!is.vector(x)) {
    stop("Input must be a vector")
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Bayesian LASSO using Gibbs Sampling
#'
#' Implements the Bayesian LASSO method through Gibbs sampling as presented by
#' Park and Casella (2008). The function performs MCMC sampling to estimate
#' regression coefficients and hyperparameters.
#'
#' @param y Response vector
#' @param X Design matrix
#' @param n_iterations Total number of MCMC iterations (default: 3000)
#' @param n_samples Number of MCMC samples to keep (default: 2000)
#' @param thin Thinning interval for MCMC samples (default: 2)
#' @param a Shape parameter for sigma2 prior (default: 1)
#' @param b Scale parameter for sigma2 prior (default: 1)
#' @param r Shape parameter for lambda2 prior (default: 1)
#' @param delta Scale parameter for lambda2 prior (default: 1)
#'
#' @return A list containing:
#' \itemize{
#'   \item beta_chain - Matrix of MCMC samples for regression coefficients
#'   \item sigma2_chain - Vector of MCMC samples for error variance
#'   \item lambda2_chain - Vector of MCMC samples for lambda2 parameter
#'   \item tau2_chain - Matrix of MCMC samples for tau2 parameters
#' }
#' @export
#' @importFrom stats rgamma rnorm
#'
#' @examples
#' X <- matrix(rnorm(100 * 5), 100, 5)
#' y <- X %*% c(1, 0.5, 0, -0.8, 0.3) + rnorm(100, 0, 0.5)
#' result <- blasso(y, X)
blasso <- function(y, X, n_iterations = 3000, n_samples = 2000, thin = 2,
                  a = 1, b = 1, r = 1, delta = 1) {
  # Input validation
  if (!is.matrix(X)) stop("Design matrix X must be a matrix")
  if (!is.vector(y)) stop("Response vector y must be a vector")
  if (nrow(X) != length(y)) stop("Number of rows in X must match length of y")
  if (n_iterations <= 0) stop("Number of iterations must be positive")
  if (n_samples <= 0) stop("Number of samples must be positive")
  if (thin <= 0) stop("Thinning interval must be positive")
  
  # Get dimensions
  n <- length(y)
  p <- ncol(X)
  
  # Initialize parameters
  beta <- matrix(0, nrow = n_samples, ncol = p)
  sigma2 <- rep(1, n_samples)
  lambda2 <- rep(1, n_samples)
  tau2 <- matrix(1, nrow = n_samples, ncol = p)
  
  # Initialize current values
  beta_current <- rep(0, p)
  sigma2_current <- 1
  lambda2_current <- 1
  tau2_current <- rep(1, p)
  
  # Center response for numerical stability
  y_centered <- y - mean(y)
  
  # MCMC iterations
  burn_in <- n_iterations - n_samples * thin
  if (burn_in < 0) stop("Number of iterations must be greater than n_samples * thin")
  
  sample_idx <- 1
  
  for(iter in 1:n_iterations) {
    # Update beta
    for(j in 1:p) {
      Vj <- 1 / (t(X[,j]) %*% X[,j] / sigma2_current + 1/tau2_current[j])
      mj <- Vj * t(X[,j]) %*% (y_centered - X[,-j] %*% beta_current[-j]) / sigma2_current
      beta_current[j] <- rnorm(1, mj, sqrt(Vj))
    }
    
    # Update sigma2
    RSS <- sum((y_centered - X %*% beta_current)^2)
    sigma2_current <- 1/rgamma(1, shape = (n+p)/2 + a, 
                              rate = (RSS + sum(beta_current^2/tau2_current))/2 + b)
    
    # Update tau2
    for(j in 1:p) {
      tau2_current[j] <- 1/rgamma(1, shape = 1, 
                                 rate = 1/lambda2_current + beta_current[j]^2/(2*sigma2_current))
    }
    
    # Update lambda2
    lambda2_current <- rgamma(1, shape = p+r, 
                            rate = sum(1/tau2_current) + delta)
    
    # Store samples
    if(iter > burn_in && (iter - burn_in) %% thin == 0) {
      beta[sample_idx,] <- beta_current
      sigma2[sample_idx] <- sigma2_current
      lambda2[sample_idx] <- lambda2_current
      tau2[sample_idx,] <- tau2_current
      sample_idx <- sample_idx + 1
    }
  }
  
  # Return results
  return(list(
    beta_chain = beta,
    sigma2_chain = sigma2,
    lambda2_chain = lambda2,
    tau2_chain = tau2
  ))
}

#' Process MCMC Chain Output
#'
#' Process MCMC chain output by applying burn-in and thinning.
#'
#' @param chain Matrix where each column is a parameter and each row is a sampled value
#' @param burn_in Number of initial samples to discard (optional)
#' @param thin Keep every nth sample (optional)
#'
#' @return Matrix of processed MCMC samples with dimensions (n_samples, n_parameters)
#' @export
#'
#' @examples
#' X <- matrix(rnorm(100 * 5), 100, 5)
#' y <- X %*% c(1, 0.5, 0, -0.8, 0.3) + rnorm(100, 0, 0.5)
#' result <- blasso(y, X)
#' processed_chain <- process_chain(result$beta_chain, burn_in = 200, thin = 5)
process_chain <- function(chain, burn_in = NULL, thin = NULL) {
  if (!is.matrix(chain)) {
    stop("Input chain must be a matrix")
  }
  
  n_samples <- nrow(chain)
  
  # Apply burn-in if specified
  if (!is.null(burn_in)) {
    if (burn_in < 0 || burn_in >= n_samples) {
      stop("Burn-in period must be between 0 and number of samples")
    }
    chain <- chain[(burn_in + 1):n_samples, , drop = FALSE]
  }
  
  # Apply thinning if specified
  if (!is.null(thin)) {
    if (thin < 1) {
      stop("Thinning interval must be positive")
    }
    chain <- chain[seq(1, nrow(chain), by = thin), , drop = FALSE]
  }
  
  return(chain)
} 