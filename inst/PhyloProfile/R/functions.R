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

#' Create link to public database
#' @return link to public database
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

createDBlink <- function(id, source, type = ""){
    linkText <- ""
    url <- ""
    if (source == "NCBI") {
        url <- paste0("https://www.ncbi.nlm.nih.gov/protein/", id)
    } else if (source == "UniProt") {
        url <- paste0("https://www.uniprot.org/uniprot/", id)
    } else if (source == "OrthoDB") {
        url <- paste0("https://www.orthodb.org/?query=", id)
        if (type == "gene") {
            idMod <- gsub(":", "%3A", id)
            url <- paste0("https://www.orthodb.org/?gene=", idMod) 
        }
    } else if (source == "OMA") {
        url <- paste0("https://omabrowser.org/oma/omagroup/", id, "/members/")
        if (type == "gene") {
            url <- paste0("https://omabrowser.org/oma/info/", id)
        }
    }
    
    if (length(url > 0)) {
        if (RCurl::url.exists(url)) {
            linkText <- paste0(
                "<p><a href='", url, "' target='_blank'>",
                source, " entry for <strong>", id, "</strong></a></p>"
            )
        }
    }
    return(linkText)
}

#' Get the most specific taxonomy rank from input taxa
#' @return the most specific taxonomy rank
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
getLowestRank <- function(longInput, taxDB) {
    inputTaxa <- gsub("ncbi","",levels(as.factor(longInput$ncbiID)))
    allDB <- getNameList(taxDB)
    allRanks <- getTaxonomyRanks()
    inputRanks <- levels(as.factor(allDB[allDB$ncbiID %in% inputTaxa,]$rank))
    if ("subspecies" %in% inputRanks) inputRanks <- c(inputRanks, "strain")
    inputRanks <- inputRanks[inputRanks %in% allRanks]
    return(allRanks[min(match(inputRanks, allRanks))][1])
}

#' Get colors for gene categoies
#' @param geneCategoryFile gene category file
#' @param type Either from manually uploaded file or predefined in config file
#' @return named vector of gene categories and their colors
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

getCatColors <- function(geneCategoryFile, type = "file"){
    catColors <- NULL
    if (type == "file") {
        inputCatDt <- read.table(
            file = geneCategoryFile$datapath,
            sep = "\t",
            header = FALSE,
            check.names = FALSE,
            comment.char = "",
            fill = TRUE
        )
        
    } else if (type == "config"){
        inputCatDt <- read.table(
            file = geneCategoryFile,
            sep = "\t",
            header = FALSE,
            check.names = FALSE,
            comment.char = "",
            fill = TRUE
        )
    }
    if (ncol(inputCatDt) == 3) {
        catColorDf <- unique(inputCatDt[,c(2,3)])
        catColors <- as.character(catColorDf$V3)
        names(catColors) <- catColorDf$V2
    }
    return(catColors)
}

#' Get if ortho ID is in BIONF format (e.g. Q6PCB6|SACCE@4932@qfo|PROTID|1)
#' @param orthoID one ortholog ID
#' @param seedID seed ID of that ortholog
#' @param ncbiID ncbi ID of that ortholog
#' @return TRUE (if it is in BIONF format) or FALSE 
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

checkBionfFormat <- function(orthoID, seedID, ncbiID) {
    ortho <- strsplit(as.character(orthoID),'|',fixed = TRUE)[[1]]
    if (length(ortho) >= 3 && ortho[1] == seedID && grepl(ncbiID, ortho[2])) 
        return(TRUE)
    return(FALSE)
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

#' Replace ~ symbol by the name of home folder
#' @export
#' @param fullPath complete path to input file or folder
#' @return complete path of input file or folder without ~ symbol
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @import stringr
#' @examples 
#' replaceHomeCharacter("~/path/to/something")
replaceHomeCharacter <- function (fullPath = NULL) {
    if (!(Sys.info()['sysname'] == "Windows")){
        homeName <- system("echo $HOME", intern = TRUE)
        return(stringr::str_replace(fullPath, "~", homeName))
    } else return(fullPath)
}

substrRight <- function(x, n) {
    substr(x, nchar(x)-n+1, nchar(x))
}

substrLeft <- function(x, n) {
    substr(x, 1, n)
}
