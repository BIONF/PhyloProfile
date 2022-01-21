#' Startup script for PhyloProfile
#' 1) install and load packages
#' 2) start the PhyloProfile app

source("R/functions.R")

# List of dependent packages --------------------------------------------------
packages <- c(
    "ape", "bioDist", "Biostrings", "colourpicker", "data.table", "energy",
    "GenomeInfoDbData", "ggplot2", "GO.db", "grid", "gridExtra", "RColorBrewer",
    "shiny", "shinyBS", "shinyFiles", "shinyjs", "OmaDB", "zoo"
)

# Load packages
lapply(packages, library, character.only = TRUE)

# Install ExperimentHub to load demo data sets
if (hasInternet() == TRUE) {
    if (!requireNamespace("ExperimentHub"))
        BiocManager::install("ExperimentHub")
    if (packageVersion("ExperimentHub") < "1.11.1")
        BiocManager::install(pkgs = "ExperimentHub", version = "devel")
    library(ExperimentHub)
    eh = ExperimentHub(localHub = TRUE)
    if ("EH2549" %in% eh$ah_id) {
        myData <- query(eh, "PhyloProfileData")
    } else {
        eh = ExperimentHub()
        myData <- query(eh, "PhyloProfileData")
    }
}
