library(tools)
library(knitr)
library(hash)

# Set working directory to directory where this script is
if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Set this variable to TRUE if you want to remove all required packages and run 
# the script from scratch, reinstalling everything along the way
START_FROM_SCRATCH <- TRUE

# Function to find all required packages in an R or Rmd script
findReqdLibraries <- function(dir, r_file){
  x <- grep('library\\(', readLines(file.path(dir, r_file)), value=TRUE)
  y <- grep('#', x, value=TRUE, invert=TRUE) # remove commented out libraries
  z <- sub("^\\s*library\\((.+)\\)\\s*$", "\\1", y) # get library names
  return(z)
}

# Function to uninstall all required packages so we can run this script from scratch
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

# Some required packages
otherpkgs <- c("glmGamPoi", "presto")

# Paths to all directories
directories <- c("ATAC/code", "barcoding/code", "dose_response/code",
                 "GO_correlations/code", "ion_flux/code", "RNA/code",
                 "scRNA/code")

# Create the hash for copying plots created by the R scripts to the figure numbers for the paper
figure_names <- hash()
# ATAC/code
# ...
"barcoding/code"
# ...
"dose_response/code"
# ...
"GO_correlations/code"
# ...
"ion_flux/code"
# ...
"RNA/code"
# ...
# scRNA/code
figure_names['UMAP_combined_SKMEL5_hg38_qcCCReg.svg'] <- 'Fig_1A.svg'
figure_names['SKMEL5_allClusters_CellCycleState_proportion.pdf'] <- 'Fig_1E.pdf'
# ...

if (START_FROM_SCRATCH == TRUE)
  try(uninstallReqdPkgs(otherpkgs, directories))

source('getReqdPkgs.r')
getReqdPkgs(otherpkgs)

# Create the '_HAYFORD_2024_PAPER_PLOTS' directory if it doesn't already exist
if (!dir.exists("_HAYFORD_2024_PAPER_PLOTS")){
  dir.create("_HAYFORD_2024_PAPER_PLOTS")
}

########## RUN ALL R AND Rmd SCRIPTS ##########

# Loop over all directories
for (dir in directories){
  message(paste("Changing directory to:", dir))
  r_files <- list.files(dir, "[:.:](r|R|Rmd)$")
  # Loop over R and Rmd scripts
  for (r_file in r_files){
    message(paste("Running script:", r_file))
    # Get all required libraries
    z <- findReqdLibraries(dir, r_file)
    getReqdPkgs(z)
    # Run the R or Rmd script
    if (file_ext(r_file) == 'Rmd'){
      outfile <- sub(".Rmd", ".R", r_file)
      try(source(knitr::purl(file.path(dir, r_file), 
                             output=file.path(dir, outfile), 
                             quiet=TRUE), chdir=TRUE))
      # remove temporarily created .R file
      if (file.exists(file.path(dir, outfile))){
        file.remove(file.path(dir, outfile))
      }
    }
    else{
      try(source(file.path(dir, r_file), chdir=TRUE))
    }
  }
  # Copy the files to the '_HAYFORD_2024_PAPER_PLOTS' directory
  for (old_name in keys(figure_names)){
    new_name <- values(figure_names, keys=old_name)
    if (file.exists(file.path(dir, old_name))){
      print(file.path(dir, old_name))
      file.copy(file.path(dir, old_name), "_HAYFORD_2024_PAPER_PLOTS")
      file.rename(file.path("_HAYFORD_2024_PAPER_PLOTS", old_name), 
                  file.path("_HAYFORD_2024_PAPER_PLOTS", new_name))
    }
  }
}


