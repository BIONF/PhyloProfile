#####################################
### 1) install and load packages  ###
### 2) start the PhyloProfile app ###
#####################################

### list of dependent packages
packages <- c("shiny","shinyBS","shinyjs","colourpicker","ggplot2","reshape2","DT","plyr","dplyr","tidyr","scales","grid","gridExtra","ape","stringr","gtable","dendextend","ggdendro","gplots","data.table","taxize","rdrop2")

### find missing packages and install them
missingPkg <- packages[!packages %in% rownames(installed.packages())]
if(length(missingPkg)){
	install.packages(missingPkg, dependencies = TRUE, repos="http://cran.us.r-project.org")
}

### check version and install ggplot2 (gplot2 v2.1 and below has some issues with the ordering)
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

### install biostrings from bioconductor
if (!("Biostrings" %in% rownames(installed.packages()))) {
	source("https://bioconductor.org/biocLite.R")
	biocLite("Biostrings")
}

### load require packages
sapply(packages, require, character.only = TRUE)
require(Biostrings)

### run phyloprofile shiny app
shiny::runApp(,launch.browser=TRUE)
