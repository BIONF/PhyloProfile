#' Startup script for PhyloProfile
#' 1) install and load packages
#' 2) start the PhyloProfile app

source("R/functions.R")

# List of dependent packages --------------------------------------------------
packages <- c(
    "ape", "bioDist", "Biostrings", "colourpicker", "data.table", "energy", 
    "GenomeInfoDbData", "ggplot2", "GO.db", "grid", "gridExtra", "RColorBrewer",
    "shiny", "shinyBS", "shinyjs", "OmaDB", "zoo"
)

# Find & install missing packages ---------------------------------------------
# installPackages(packages)

# Load packages
lapply(packages, library, character.only = TRUE)

# Install packages from bioconductor ------------------------------------------
# bioconductor_pkgs <- c("Biostrings", "bioDist")
# installPackagesBioconductor(bioconductor_pkgs)
# lapply(bioconductor_pkgs, library, character.only = TRUE)

# Install OmaDB's dependencies "GO.db", "GenomeInfoDbData"
# and ExperimentHub
# oma_pkgs <- c("GO.db", "GenomeInfoDbData", "ExperimentHub")
# installPackagesBioconductor(oma_pkgs)
if (!requireNamespace("ExperimentHub"))
    BiocManager::install("ExperimentHub")
if (packageVersion("ExperimentHub") < "1.11.1")
    BiocManager::install(pkgs = "ExperimentHub", version = "devel")
library(ExperimentHub)

# Load demo data from PhyloProfileData package
eh = ExperimentHub()
myData <- query(eh, "PhyloProfileData")
