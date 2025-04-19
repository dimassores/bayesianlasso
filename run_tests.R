# Load required packages
library(testthat)
library(devtools)

# Load the package
message("Loading package...")
devtools::load_all()

# Run tests
message("\nRunning tests...")
test_results <- test_dir("tests/testthat", reporter = "progress")

# Print summary
message("\nTest Summary:")
print(test_results) 