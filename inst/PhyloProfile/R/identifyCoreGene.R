#' Core gene identification
#' @param filteredData full processed main data
#' (from reactive fn "getDataFiltered")
#' @param rankSelect selected taxonomy rank (input$rankSelect)
#' @param taxaCore selected list of taxa (input$taxaCore)
#' @param percentCore cutoff of percentage taxa present in a supertaxon
#' (input$percentCore)
#' @return list of core genes
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

# source("R/functions.R")

identifyCoreGeneUI <- function(id){
    ns <- NS(id)
    DT::dataTableOutput(ns("coreGene.table"))
}

identifyCoreGene <- function(
    input, output, session,
    filteredData, rankSelect, taxaCore, percentCore,
    var1Cutoff, var2Cutoff, coreCoverage
){

    output$coreGene.table <- DT::renderDataTable({
        data <- coreGeneDf()
        if (is.null(data)) stop("Input data is NULL!")
        else {
            data <- as.data.frame(data)
            data
        }
    })

    coreGeneDf <- reactive({
        coreGeneDf <- getCoreGene(
            rankSelect(), taxaCore(), filteredData(),
            var1Cutoff(), var2Cutoff(), percentCore(), coreCoverage()
        )
        return(coreGeneDf)
    })

    if (!is.null(coreGeneDf)) return(coreGeneDf)
}
