#' Startup script for PhyloProfile
#' 1) install and load packages
#' 2) start the PhyloProfile app
library(PhyloProfile)
source("R/functions.R")

# List of dependent packages --------------------------------------------------
packages <- c(
    "ape", "BiocStyle", "bsplus", "data.table", "dplyr", "ggplot2", "gridExtra",
    "htmlwidgets", "shiny", "shinyFiles", "shinyjs", "scattermore", 
    "svglite", "plotly")

# Load packages
lapply(packages, library, character.only = TRUE)
