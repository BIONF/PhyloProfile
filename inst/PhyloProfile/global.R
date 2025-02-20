#' Startup script for PhyloProfile
#' 1) install and load packages
#' 2) start the PhyloProfile app
library(PhyloProfile)
source("R/functions.R")

# List of dependent packages --------------------------------------------------
packages <- c("ape", "BiocStyle","data.table", "dplyr", "ggplot2", "shiny", "shinyBS",
    "shinyFiles", "shinyjs", "scattermore", "svglite", "plotly"
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
    eh = ExperimentHub(localHub = FALSE)
    if ("EH2549" %in% eh$ah_id) {
        myData <- query(eh, "PhyloProfileData")
    } else {
        eh = ExperimentHub()
        myData <- query(eh, "PhyloProfileData")
    }
}
