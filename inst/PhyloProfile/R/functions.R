#' Function to keep user defined geneID order
#' @param data data frame contains gene ID column
#' @param order TRUE or FALSE (from input$ordering)
#' @return data either sorted or non-sorted
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
unsortID <- function(data, order){
    data$geneID <- as.factor(data$geneID)
    if (order == FALSE) {
        # keep user defined geneID order
        data$geneID <- factor(data$geneID, levels = unique(data$geneID))
    }
    return(data)
}

#' Check installed packages
#' and install missing packages automatically
#' @param packages list of packages need to be checked
#' @return none
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

installPackages <- function(packages){
    missingPackages <-
        packages[!(packages %in% installed.packages()[, "Package"])]
    if (length(missingPackages))
        install.packages(
            missingPackages,
            dependencies = TRUE,
            repos = "http://cran.us.r-project.org"
        )
}

installPackagesBioconductor <- function(packages){
    missingPackages <-
        packages[!(packages %in% installed.packages()[, "Package"])]
    if (length(missingPackages)) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install(missingPackages, ask = FALSE)
    }
}

#' Check internet connection
#' @return status of internet connection
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

hasInternet <- function(){
    !is.null(curl::nslookup("r-project.org", error = FALSE))
}

# FUNCTIONS FOR RENDER UI ELEMENTS ============================================
createSliderCutoff <- function(id, title, start, stop, varID){
    if (is.null(varID)) return()
    if (varID == "") {
        sliderInput(id, title,
                    min = 1,
                    max = 1,
                    step = 0.025,
                    value = 1,
                    width = 200)
    } else {
        sliderInput(id, title,
                    min = 0,
                    max = 1,
                    step = 0.025,
                    value = c(start, stop),
                    width = 200)
    }
}

updateSliderCutoff <- function(session, id, title, newVar, varID){
    if (is.null(varID) || varID == "") return()
    updateSliderInput(session, id, title,
                      value = newVar,
                      min = 0,
                      max = 1,
                      step = 0.025)
}

createPlotSize <- function(id, title, value) {
    numericInput(id,
                 title,
                 min = 100,
                 max = 3200,
                 step = 50,
                 value = value,
                 width = 100)
}

createTextSize <- function(id, title, value, width) {
    numericInput(id,
                 title,
                 min = 3,
                 max = 99,
                 step = 1,
                 value = value,
                 width = width)
}

createSelectGene <- function(id, list, selected) {
    selectInput(id,
                "",
                list,
                selected = selected,
                multiple = TRUE,
                selectize = FALSE)
}
