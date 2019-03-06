#' Startup script for PhyloProfile
#' 1) install and load packages
#' 2) start the PhyloProfile app

source("R/functions.R")

# List of dependent packages --------------------------------------------------
packages <- c("shiny", "shinyBS", "shinyjs", "colourpicker", "DT",
              "devtools", "ggplot2", "reshape2",
              "plyr", "dplyr", "tidyr", "scales", "grid",
              "gridExtra", "ape", "stringr", "gtable",
              "dendextend", "ggdendro", "gplots", "data.table",
              "taxize", "zoo", "RCurl", "energy",
              "RColorBrewer")

# Set path for install packages while deploy into shiny server ----------------
# (from https://gist.github.com/wch/c3653fb39a00c63b33cf)

# Find & install missing packages ---------------------------------------------
installPackages(packages)

# Load packages
lapply(packages, library, character.only = TRUE)

# Check version and install ggplot2 (require v >= 2.2.0) ----------------------
version_above <- function(pkg, than) {
    compareVersion(as.character(packageVersion(pkg)), than)
}

if ("ggplot2" %in% rownames(installed.packages())) {
    installPackages("ggplot2")
    library(ggplot2)
}

# Install packages from bioconductor ------------------------------------------
bioconductor_pkgs <- c("Biostrings", "bioDist")
installPackagesBioconductor(bioconductor_pkgs)
lapply(bioconductor_pkgs, library, character.only = TRUE)

# Install OmaDB and its dependencies
oma_pkgs <- c("GO.db", "GenomeInfoDbData")
installPackagesBioconductor(oma_pkgs)
lapply(oma_pkgs, library, character.only = TRUE)

if (!("OmaDB" %in% rownames(installed.packages()))) {
    devtools::install_github("trvinh/OmaDB", force = TRUE)
}
library(OmaDB)

# Install shinycssloaders from github -----------------------------------------
if (!("shinycssloaders" %in% rownames(installed.packages()))) {
    devtools::install_github("andrewsali/shinycssloaders", force = TRUE)
    library(shinycssloaders)
}
