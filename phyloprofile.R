#' Startup script for PhyloProfile
#' 1) install and load packages
#' 2) start the PhyloProfile app

# List of dependent packages --------------------------------------------------
packages <- c("shiny", "shinyBS", "shinyjs", "colourpicker", "DT",
              "devtools", "ggplot2", "reshape2", 
              "plyr", "dplyr", "tidyr", "scales", "grid", 
              "gridExtra", "ape", "stringr", "gtable", 
              "dendextend", "ggdendro", "gplots", "data.table", 
              "taxize", "zoo", "RCurl")

# Find & install missing packages ---------------------------------------------
missingPkg <- packages[!packages %in% rownames(installed.packages())]
if (length(missingPkg)) {
	install.packages(
	  missingPkg,
	  dependencies = TRUE,
	  repos = "http://cran.us.r-project.org"
	)
}

# Check version and install ggplot2 (require v >= 2.2.0) ----------------------
version_above <- function(pkg, than) {
  compareVersion(as.character(packageVersion(pkg)), than)
}

if ("ggplot2" %in% rownames(installed.packages())) {
  if (version_above("ggplot2","2.2.0") == -1) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("ggplot2")
  }
} else {
  source("https://bioconductor.org/biocLite.R")
  biocLite("ggplot2")
}

# Install biostrings from bioconductor ----------------------------------------
if (!("Biostrings" %in% rownames(installed.packages()))) {
	source("https://bioconductor.org/biocLite.R")
	biocLite("Biostrings")
	library(Biostrings)
}

# Install shinycssloaders from github -----------------------------------------
if (!("shinycssloaders" %in% rownames(installed.packages()))) {
	devtools::install_github('andrewsali/shinycssloaders', force = TRUE)
  library(shinycssloaders)
}

### load require packages
lapply(packages, library, character.only = TRUE)

### run phyloprofile shiny app
shiny::runApp(appDir = getwd(), launch.browser = TRUE)
