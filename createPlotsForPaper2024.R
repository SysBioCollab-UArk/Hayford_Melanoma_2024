library(tools)
library(knitr)
library(hash)
library(devtools)
library(purrr)
library(stringr)

# Clear the session memory 
# NOTE: Need to do this to prevent memory exhaustion while running the script
rm(list=ls())

# Set the following variable to TRUE if you want to remove all required packages
# and run the script from scratch, reinstalling everything along the way
# *****
# NOTE: This basically works but often requires restarting the R session and
# then rerunning the script with START_FROM_SCRATCH == FALSE. We may eventually
# want to remove this in favor of a virtual environment --LAH 
START_FROM_SCRATCH <- FALSE

# Set the following variable to TRUE to run the script in 'debug' mode 
# (i.e., press [enter] before each script in the repo runs)
DEBUG <- TRUE

# Create the hash for copying plots created by the R scripts to the figure 
# numbers for the paper
figure_names <- hash()

# ATAC/code
figure_names['GOenrichment_sub25_I_MF.pdf'] <- 'Fig_3C.pdf'
figure_names['cCurve_colors_SKMEL5.pdf'] <- 'Fig_S3A.pdf'
figure_names['ISM_colors_SKMEL5.pdf'] <- 'Fig_S3B.pdf'
figure_names['venn_sharedPeaks.pdf'] <- 'Fig_S3C.pdf'
figure_names['ATAC_sub25_distanceToTSS_UniqueShared.pdf'] <- 'Fig_S3D.pdf'
figure_names['ATAC_sub25_annotationDistribution_UniqueShared.pdf'] <- 'Fig_S3E.pdf'

# barcoding/code
figure_names['SKMEL5_barcode_propRankAbundance_comparison.pdf'] <- 'Fig_2A.pdf'
figure_names['SKMEL5_barcode_FCdensity_bcOverlay.pdf'] <- 'Fig_2B.pdf'
figure_names['SKMEL5_barcode_RPM_rank.pdf'] <- 'Fig_S2A.pdf'
figure_names['SKMEL5_numUniqueBC_byCondition.pdf'] <- 'Fig_S2B.pdf'
figure_names['SKMEL5_propShared_byReplicate.pdf'] <- 'Fig_S2C.pdf'

# dose_response/code
figure_names['RSL3-all_DRC_2024.pdf'] <- 'Fig_4B.pdf'

# GO_correlations/code
figure_names['RNA-ATACsub25_MF.svg'] <- 'Fig_3D.svg'
figure_names['RNA-ATACsub25_BP.svg'] <- 'Fig_S3F.svg'
figure_names['RNA-ATACsub25_CC.svg'] <- 'Fig_S3G.svg'

# ion_flux/code
figure_names['SOCE_untreated_idling.pdf'] <- 'Fig_3E.pdf'

# RNA/code
figure_names['SKMEL5_sublines_timeSeriesRNA_rld_Fig.pdf'] <- 'Fig_1B.pdf'
figure_names['GOenrichment_genesUp_MF.pdf'] <- 'Fig_3A.pdf'
figure_names['ca-genes_fold-change_ComplexHeatmap.pdf'] <- 'Fig_3B.pdf'
figure_names['Ferroptosis_heatmap_annot.pdf'] <- 'Fig_4A_main.pdf'
figure_names['Ferroptosis_heatmap_annot_legend.pdf'] <- 'Fig_4A_leg.pdf'

# scRNA/code
figure_names['UMAP_combined_SKMEL5_hg38_qcCCReg_treatmentPoint.svg'] <- 'Fig_1A.svg'
figure_names['subsetHallmarks_acrossSKMEL5clusters.svg'] <- 'Fig_1C.svg'
figure_names['UMAP_combined_SKMEL5_hg38_qcCCReg_treatmentDensity_CCStatePoint.svg'] <- 'Fig_1D.svg'
figure_names['SKMEL5_allClusters_CellCycleState_proportion.pdf'] <- 'Fig_1E.pdf'
figure_names['UMAP_combined_lineageID_tinted_BCs_2_5_13_9.svg'] <- 'Fig_2C.svg'
figure_names['SKMEL5_I_CellCycleState_onlyDividing_proportionWAverage.pdf'] <- 'Fig_2D.pdf'
figure_names['Idling_bcFC_correlation.pdf'] <- 'Fig_2E.pdf'
figure_names['UMAP_combined_lineageID_tinted_BCs_others.svg'] <- 'Fig_S2D.svg'
figure_names['Untreated_smallClusterGO_BP.pdf'] <- 'Fig_S1A.pdf'
figure_names['Idling_smallClusterGO_BP.pdf'] <- 'Fig_S1B.pdf'

# Set working directory to directory where this script is
if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Function to find all required packages in an R or Rmd script
findReqdLibraries <- function(dir, r_file){
  x <- grep('library\\(', readLines(file.path(dir, r_file)), value=TRUE)
  y <- grep('#', x, value=TRUE, invert=TRUE) # remove commented out libraries
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

# Some required packages
otherpkgs <- c("glmGamPoi", "presto")

# Paths to all directories
directories <- c("barcoding/code", "dose_response/code",
                 "GO_correlations/code", "ion_flux/code", "RNA/code",
                 "scRNA/code", "ATAC/code")

# directories <- c("scRNA/code", "ATAC/code")

# Uninstall all required packages if starting from scratch
if (START_FROM_SCRATCH == TRUE)
  try(uninstallReqdPkgs(otherpkgs, directories))

# Install all required packages
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
    
    # To prevent memory exhaustion, delete all objects created in the previous 
    # R script before running the next one
    if (dir == directories[1]){
      # Don't delete objects from the main script
      save_list <- c(rep(ls()), "save_list") 
    }
    else{
      # Delete all objects in the workspace except the ones in the main script
      rm(list=ls()[!(ls() %in% save_list)])
    }
    
    # This pauses the code before each script in the repo runs, so the user
    # can keep track of what's going on
    if (DEBUG == TRUE) readline(prompt="Press [enter] to continue")
    
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
    else{ # r or R script
      try(source(file.path(dir, r_file), chdir=TRUE))
    }
  }
  
  # COPY FIGURE FILES TO '_HAYFORD_2024_PAPER_PLOTS' DIRECTORY
  for (old_name in hash::keys(figure_names)){
    if (file.exists(file.path(dir, old_name))){
      # Copy the file to the new directory
      file.copy(file.path(dir, old_name), "_HAYFORD_2024_PAPER_PLOTS")
      # Now change its name
      new_name <- hash::values(figure_names, keys=old_name)
      file.rename(file.path("_HAYFORD_2024_PAPER_PLOTS", old_name),
                  file.path("_HAYFORD_2024_PAPER_PLOTS", new_name))
      # Tell the user about it
      message(paste("Copied", old_name, 
                    "to '_HAYFORD_2024_PAPER_PLOTS' and renamed it", new_name))
    }
  }
}
