library(tools)
library(knitr)

# Set working directory to directory where this script is
if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

source('getReqdPkgs.r')

# Some required packages
getReqdPkgs(c("glmGamPoi", "presto"))

# Paths to all directories
directories <- c("ATAC/code", "barcoding/code", "dose_response/code",
                 "GO_correlations/code", "ion_flux/code", "RNA/code",
                 "scRNA/code")

# Create the 'HAYFORD_2024_PAPER_PLOTS' directory if it doesn't already exist
if (!dir.exists("HAYFORD_2024_PAPER_PLOTS")){
  dir.create("HAYFORD_2024_PAPER_PLOTS")
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
    x <- grep('library\\(', readLines(file.path(dir, r_file)), value=TRUE)
    y <- grep('#', x, value=TRUE, invert=TRUE) # remove commented out libraries
    z <- sub("^\\s*library\\((.+)\\)\\s*$", "\\1", y) # get library names
    getReqdPkgs(z)
    # Run the R or Rmd script
    if (file_ext(r_file) == 'Rmd'){
      outfile <- sub(".Rmd", ".R", r_file)
      try(source(knitr::purl(file.path(dir, r_file), 
                             output=file.path(dir, outfile), 
                             quiet=TRUE), chdir=TRUE))
      if (file.exists(file.path(dir, outfile))){
        file.remove(file.path(dir, outfile))
      }
    }
    else{
      try(source(file.path(dir, r_file), chdir=TRUE))
    }
  }
  # Copy the files to the 'HAYFORD_2024_PAPER_PLOTS' directory
  # ...
}
