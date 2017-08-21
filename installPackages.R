# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.
# downloaded from https://gist.github.com/stevenworthington/3178163

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg))
        install.packages(new.pkg, dependencies = TRUE, repos="http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE)
}

# usage
packages <- c("shiny","shinyBS","shinyjs","colourpicker","ggplot2","reshape2","DT","plyr","dplyr","tidyr","scales","grid","gridExtra","ape","stringr","gtable","dendextend","ggdendro","gplots","data.table","taxize","rdrop2")
ipak(packages)

### install biostrings from bioconductor
if ("ggplot2" %in% rownames(installed.packages())) {
	source("https://bioconductor.org/biocLite.R")
	biocLite("Biostrings")
}

### check version and install ggplot2
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

### load require packages
lapply(packages, require, character.only = TRUE)

### run phyloprofile shiny app
shiny::runApp()
