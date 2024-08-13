# Set working directory to the directory where this script is
if (Sys.getenv("RSTUDIO") == "1")
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Function to find all required packages in an R or Rmd script
findReqdLibraries <- function(dir, r_file){
  x <- grep('library\\(', readLines(file.path(dir, r_file)), value=TRUE)
  y <- grep('^\\s*#', x, value=TRUE, invert=TRUE) # remove commented out libraries
  z <- sub("^\\s*library\\((.+)\\)\\s*$", "\\1", y) # get library names
  return(z)
}

# Function to uninstall all required packages so script can be run from scratch
uninstallReqdPkgs <- function(pkgs, directories){
  message("Uninstalling all required packages")
  remove.packages(pkgs)
  # Loop over all directories
  for (dir in directories){
    message(paste("Changing directory to:", dir))
    r_files <- list.files(dir, "[:.:](r|R|Rmd)$")
    # Loop over R and Rmd scripts
    for (r_file in r_files){
      message(paste("Uninstalling packages from:", r_file))
      # Uninstall all required packages
      z <- findReqdLibraries(dir, r_file)
      remove.packages(z)
    }
  }
}

# Function to install required packages
getReqdPkgs <- function(req_pkgs) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    inst_pkgs <- installed.packages()[,"Package"]
    needed_pkgs <- setdiff(req_pkgs,inst_pkgs)
    std_pkgs <- needed_pkgs[
      needed_pkgs %in% tools::CRAN_package_db()[,"Package"]]
    bioc_pkgs <- needed_pkgs[needed_pkgs %in% BiocManager::available()]
    other_pkgs <- setdiff(needed_pkgs,c(std_pkgs,bioc_pkgs))

    if(!length(other_pkgs)==0) {
      message(cat("Packages not found in CRAN or BioC still needed:",
                    paste(other_pkgs, collapse=", ")))
      require(devtools)
      # Install VISION package
      if ('VISION' %in% other_pkgs) 
        install_github("YosefLab/VISION", quiet=TRUE)
      # Install presto package
      if ('presto' %in% other_pkgs)
        install_github("immunogenomics/presto", quiet=TRUE)
      # Install diprate package
      if ('diprate' %in% other_pkgs)
        install_github("QuLab-VU/dipDRC", subdir="diprate", dependencies=TRUE, 
                       quiet=TRUE)
    }

    if(length(std_pkgs)==0 & length(bioc_pkgs)==0) {
      message('All required CRAN and BioC packages have previously been installed')
    } 
    else {
      if(!is.null(std_pkgs)) install.packages(std_pkgs)
      if(!is.null(bioc_pkgs)) BiocManager::install(bioc_pkgs, ask=FALSE)
    }
}

##### MAIN SCRIPT #####

# Set to TRUE if you want to remove all required packages and run the script 
# from scratch, reinstalling everything along the way
START_FROM_SCRATCH <- FALSE

# Paths to all directories
directories <- c(file.path("barcoding", "code"), 
                 file.path("dose_response", "code"),
                 file.path("GO_correlations", "code"), 
                 file.path("ion_flux", "code"), 
                 file.path("RNA", "code"),
                 file.path("scRNA", "code"), 
                 file.path("ATAC", "code"))

# Some required packages
otherpkgs <- c("glmGamPoi", "presto")

# Uninstall all required packages if starting from scratch
if (START_FROM_SCRATCH){
  try(uninstallReqdPkgs(otherpkgs, directories))
  # reinstall the following packages, which are needed to run this script
  install.packages("purrr")
  install.packages("stringr")
}

# Install all required packages
getReqdPkgs(otherpkgs)
# Loop over all directories
for (dir in directories){
  message(paste("Changing directory to:", dir))
  r_files <- list.files(dir, "[:.:](r|R|Rmd)$")
  # Loop over R and Rmd scripts
  for (r_file in r_files){
    message(paste("Inspecting script:", r_file))
    # Get all required libraries
    getReqdPkgs(findReqdLibraries(dir, r_file))
  }
}
