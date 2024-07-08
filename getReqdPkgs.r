
getReqdPkgs <- function(req_pkgs) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    
    inst_pkgs <- installed.packages()[,"Package"]
    
    needed_pkgs <- setdiff(req_pkgs,inst_pkgs)
    std_pkgs <- needed_pkgs[needed_pkgs %in% 
                                tools::CRAN_package_db()[,"Package"]]
    bioc_pkgs <- needed_pkgs[needed_pkgs %in% BiocManager::available()]
    other_pkgs <- setdiff(needed_pkgs,c(std_pkgs,bioc_pkgs))
    
    if(!length(other_pkgs)==0) {
      message(cat("Packages not found in CRAN or BioC still needed:",
                    paste(other_pkgs, collapse=", ")))
      require(devtools)
      # Install VISION package
      if ('VISION' %in% other_pkgs){
        install_github("YosefLab/VISION")
      }
      # Install presto package
      if ('presto' %in% other_pkgs){
        install_github("immunogenomics/presto")
      }
      # Install diprate package
      if ('diprate' %in% other_pkgs){
        install_github("QuLab-VU/dipDRC", subdir="diprate", dependencies=TRUE)
      }
    }
    
    if(length(std_pkgs)==0 & length(bioc_pkgs)==0) {
        message('All required CRAN and BioC packages have previously been installed')
    } else {
        if(!is.null(std_pkgs)) install.packages(std_pkgs)
        if(!is.null(bioc_pkgs)) BiocManager::install(bioc_pkgs)
    }
}
