# Clear environment variables and free memory
rm(list = ls(globalenv()))
gc()
# Two package required
library(httr)
library(RcppProgress)

# Define file URL and download location
url <- "https://drive.google.com/uc?authuser=0&id=1ERM3z5G24DTKebOzyhGIFxfekjdx3vBA&export=download"
temp <- file.path(getwd(), "downloaded_package.tar.gz") # Save in the working directory (C:/Users/Bastien/OneDrive/Documents)

# Download the file
cat("Attempting to download the file...\n")
tryCatch({
  response <- GET(url, write_disk(temp, overwrite = TRUE))
  if (file.exists(temp)) {
    cat("File downloaded successfully to:", temp, "\n")
  } else {
    stop("Failed to download the file. Check the URL or internet connection.")
  }
}, error = function(e) {
  stop("Error during file download:", e$message)
})

# Inspect the file
cat("Inspecting the file...\n")
tryCatch({
  file_contents <- untar(temp, list = TRUE)
  cat("File contains:\n")
  print(file_contents)
}, error = function(e) {
  stop("Error inspecting the file. It may be corrupted.")
})

# Install missing dependencies (if any)
cat("Installing missing dependencies...\n")
required_dependencies <- c("RcppProgress") # Add more dependencies here if needed
installed_dependencies <- installed.packages()[, "Package"]

for (dep in required_dependencies) {
  if (!(dep %in% installed_dependencies)) {
    cat("Installing dependency:", dep, "\n")
    install.packages(dep)
  }
}

# Attempt to install the package
cat("Installing the main package...\n")
tryCatch({
  install.packages(temp, repos = NULL, type = "source")
  cat("Package installed successfully!\n")
}, error = function(e) {
  stop("Installation failed. Error:", e$message)
})

# Cleanup
unlink(temp)
cat("Script completed successfully.\n")

# Load library
library(lidR)
library(lidRtls)



