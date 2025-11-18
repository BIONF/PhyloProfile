#' Core gene identification
#' @param filteredData full processed main data
#' (from reactive fn "getDataFiltered")
#' @param rankSelect selected taxonomy rank (input$rankSelect)
#' @param taxaCore selected list of taxa (input$taxaCore)
#' @param percentCore cutoff of percentage taxa present in a supertaxon
#' @param var1Cutoff variable 1 cutoff
#' @param var2Cutoff variable 2 cutoff
#' @param coreCoverage the least percentage of selected taxa should be
#' considered
#' @param taxDB path to the taxonomy DB files
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
    filteredData, taxaCount, rankSelect, taxaCore, percentCore,
    var1Cutoff, var2Cutoff, coreCoverage, taxDB
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
        req(taxaCore())
        req(filteredData())
        coreGeneDf <- getCoreGene(
            rankSelect(), taxaCore(), filteredData(), taxaCount(),
            var1Cutoff(), var2Cutoff(), percentCore(), coreCoverage(), taxDB()
        )
        return(coreGeneDf)
    })

    if (!is.null(coreGeneDf)) return(coreGeneDf)
}
