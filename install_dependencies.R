# Function to install packages with error handling
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing %s...", pkg))
    install.packages(pkg)
  } else {
    message(sprintf("%s is already installed.", pkg))
  }
}

# List of required packages
required_packages <- c(
  "mvtnorm",
  "invgamma",
  "SuppDists",
  "glmnet",
  "testthat",
  "knitr",
  "rmarkdown",
  "devtools"
)

# Install required packages
message("Installing required packages...")
for (pkg in required_packages) {
  install_if_missing(pkg)
}

# Install package dependencies
message("\nInstalling package dependencies...")
devtools::install_deps()

message("\nAll dependencies installed successfully!") 