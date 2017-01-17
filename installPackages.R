libs=c("shiny","shinyBS","ggplot2","reshape","plyr","dplyr","scales‌​","grid","gridExtra","ape","colourpicker","shinyjs","stringr","svglite")
type=getOption("pkgType")

    CheckInstallPackage <- function(packages, repos="http://cran.r-project.org",
       depend="Depends", ...) {
         installed=as.data.frame(installed.packages())
    for(p in packages) {
        if(is.na(charmatch(p, installed[,1]))) { 
          install.packages(p, repos=repos, dependencies=depend, ...) 
        }
      }
    } 
    CheckInstallPackage(packages=libs)

source("https://bioconductor.org/biocLite.R")
biocLite("Biostrings")