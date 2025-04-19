test_that("Mode function works correctly", {
  # Test numeric vectors
  expect_equal(Mode(c(1, 2, 2, 3)), 2)
  expect_equal(Mode(c(1, 1, 2, 2, 2)), 2)
  
  # Test character vectors
  expect_equal(Mode(c("a", "b", "b", "c")), "b")
  
  # Test error handling
  expect_error(Mode(matrix(1:4, 2, 2)), "Input must be a vector")
})

test_that("process_chain works correctly", {
  # Create test matrix
  test_matrix <- matrix(1:100, ncol = 4)
  
  # Test basic functionality
  processed <- process_chain(test_matrix, burn_in = 5, thin = 2)
  expect_equal(nrow(processed), floor((nrow(test_matrix) - 5) / 2) + 1)
  expect_equal(ncol(processed), ncol(test_matrix))
  
  # Test error handling
  expect_error(process_chain("not a matrix"), "Input chain must be a matrix")
  expect_error(process_chain(test_matrix, burn_in = -1), "Burn-in period must be between 0 and number of samples")
  expect_error(process_chain(test_matrix, thin = 0), "Thinning interval must be positive")
})

test_that("blasso input validation works", {
  # Create test data
  X <- matrix(rnorm(100 * 5), 100, 5)
  y <- X %*% c(1, 0.5, 0, -0.8, 0.3) + rnorm(100, 0, 0.5)
  
  # Test input validation
  expect_error(blasso(y, as.data.frame(X)), "Design matrix X must be a matrix")
  expect_error(blasso(as.matrix(y), X), "Response vector y must be a vector")
  expect_error(blasso(y, X[1:50,]), "Number of rows in X must match length of y")
  expect_error(blasso(y, X, n_iterations = 0), "Number of iterations must be positive")
  expect_error(blasso(y, X, n_samples = 0), "Number of samples must be positive")
  expect_error(blasso(y, X, thin = 0), "Thinning interval must be positive")
  expect_error(blasso(y, X, n_iterations = 100, n_samples = 100, thin = 2), 
               "Number of iterations must be greater than n_samples * thin")
})

test_that("blasso produces correct output structure", {
  # Create test data
  X <- matrix(rnorm(100 * 5), 100, 5)
  y <- X %*% c(1, 0.5, 0, -0.8, 0.3) + rnorm(100, 0, 0.5)
  
  # Run with minimal iterations for testing
  result <- blasso(y, X, n_iterations = 100, n_samples = 50, thin = 2)
  
  # Check output structure
  expect_type(result, "list")
  expect_equal(length(result), 4)
  expect_true(all(c("beta_chain", "sigma2_chain", "lambda2_chain", "tau2_chain") %in% names(result)))
  
  # Check dimensions
  expect_equal(nrow(result$beta_chain), 50)
  expect_equal(ncol(result$beta_chain), 5)
  expect_equal(length(result$sigma2_chain), 50)
  expect_equal(length(result$lambda2_chain), 50)
  expect_equal(nrow(result$tau2_chain), 50)
  expect_equal(ncol(result$tau2_chain), 5)
})

test_that("blasso produces reasonable estimates", {
  # Create test data with known coefficients
  set.seed(123)
  n <- 100
  p <- 5
  true_beta <- c(1, 0.5, 0, -0.8, 0.3)
  X <- matrix(rnorm(n * p), n, p)
  y <- X %*% true_beta + rnorm(n, 0, 0.5)
  
  # Run MCMC
  result <- blasso(y, X, n_iterations = 1000, n_samples = 500, thin = 2)
  
  # Check if estimates are close to true values
  beta_means <- colMeans(result$beta_chain)
  expect_true(all(abs(beta_means - true_beta) < 0.5))
  
  # Check if variance estimates are positive
  expect_true(all(result$sigma2_chain > 0))
  expect_true(all(result$lambda2_chain > 0))
  expect_true(all(result$tau2_chain > 0))
})

test_that("blasso handles different prior parameters", {
  # Create test data
  X <- matrix(rnorm(100 * 5), 100, 5)
  y <- X %*% c(1, 0.5, 0, -0.8, 0.3) + rnorm(100, 0, 0.5)
  
  # Test with different prior parameters
  result1 <- blasso(y, X, a = 2, b = 2, r = 2, delta = 2)
  result2 <- blasso(y, X, a = 0.5, b = 0.5, r = 0.5, delta = 0.5)
  
  # Check that results are different but still valid
  expect_false(identical(result1$beta_chain, result2$beta_chain))
  expect_true(all(result1$sigma2_chain > 0))
  expect_true(all(result2$sigma2_chain > 0))
}) 