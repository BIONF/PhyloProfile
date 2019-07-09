#' Startup script for PhyloProfile
#' 1) install and load packages
#' 2) start the PhyloProfile app

source("R/functions.R")

# List of dependent packages --------------------------------------------------
packages <- c("ape", "colourpicker", "data.table",
              "plyr", "DT", "energy", "ggplot2", "gplots",
              "grid", "gridExtra", "gtable", "RCurl", "RColorBrewer",
              "reshape2", "scales", "shiny", "shinyBS", "shinyjs",
              "tidyr", "zoo")

# Find & install missing packages ---------------------------------------------
installPackages(packages)

# Load packages
lapply(packages, library, character.only = TRUE)

# Install packages from bioconductor ------------------------------------------
bioconductor_pkgs <- c("Biostrings", "bioDist")
installPackagesBioconductor(bioconductor_pkgs)
lapply(bioconductor_pkgs, library, character.only = TRUE)

# Install OmaDB's dependencies "GO.db", "GenomeInfoDbData"
# and ExperimentHub
oma_pkgs <- c("GO.db", "GenomeInfoDbData", "ExperimentHub")
installPackagesBioconductor(oma_pkgs)
if (packageVersion("ExperimentHub") < "1.11.1")
    BiocManager::install(pkgs = "ExperimentHub", version = "devel")
lapply(oma_pkgs, library, character.only = TRUE)

# Load demo data from PhyloProfileData package
eh = ExperimentHub()
myData <- query(eh, "PhyloProfileData")
