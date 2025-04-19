# The Bayesian LASSO

This package implements the Bayesian LASSO method through Gibbs sampling as presented by Park and Casella (2008). It provides tools for performing Bayesian LASSO regression, including MCMC sampling, diagnostics, and visualization tools.

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("dimassores/bayesianlasso")
```

## Usage

Here's a basic example of how to use the Bayesian LASSO:

```r
library(bayesianlasso)

# Generate example data
X <- matrix(rnorm(100 * 5), 100, 5)
y <- X %*% c(1, 0.5, 0, -0.8, 0.3) + rnorm(100, 0, 0.5)

# Fit the Bayesian LASSO model
result <- bayesian_lasso_MCMC(
  maxit = 10000,    # Number of MCMC iterations
  X = X,            # Design matrix
  Y = y             # Response vector
)

# Process the MCMC chain
cleaned_chain <- chain_cleaner(
  result$beta,
  burn_in = 2000,   # Number of initial samples to discard
  thinning = 5      # Keep every 5th sample
)

# Calculate posterior means
colMeans(cleaned_chain)
```

## Features

- Implementation of the Bayesian LASSO using Gibbs sampling
- MCMC chain processing with burn-in and thinning
- Input validation and error checking
- Comprehensive documentation and examples

## Documentation

For more detailed information about the functions and their usage, see the package documentation:

```r
?bayesian_lasso_MCMC
?chain_cleaner
```

## References

Park, Trevor, and George Casella. "The bayesian lasso." Journal of the American Statistical Association 103.482 (2008): 681-686.

Tibshirani, Robert. "Regression shrinkage and selection via the lasso." Journal of the Royal Statistical Society: Series B (Methodological) 58.1 (1996): 267-288.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
