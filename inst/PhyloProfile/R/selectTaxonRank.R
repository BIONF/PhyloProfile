#' Popup windows for (sub-)selecting rank & (super)taxon of interest
#' @param rankSelect initial selected taxonomy rank (input$rankSelect)
#' @param inputTaxonID list of all input taxon IDs
#' @taxDB Path to the taxonomy DB files
#' (from reactive fn inputTaxonID)
#' @return list of all taxa in the initial taxonomy rank based on selected
#' supertaxon from input$taxaSelectCus
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

selectTaxonRankUI <- function(id) {
    ns <- NS(id)
    tagList(
        uiOutput(ns("rankSelectCus.ui")),
        uiOutput(ns("taxaSelectCus.ui"))
    )
}

selectTaxonRank <- function(
    input, output, session, rankSelect, inputTaxonID, taxDB
){

    # render list of available taxonomy ranks
    # (where the lowest rank is the same as the chosen main rank)
    output$rankSelectCus.ui <- renderUI({
        ns <- session$ns
        mainRank <- rankSelect()
        mainChoices <- list("Strain " = "strain",
                            "Species" = "species",
                            "Genus" = "genus",
                            "Family" = "family",
                            "Order" = "order",
                            "Class" = "class",
                            "Phylum" = "phylum",
                            "Kingdom" = "kingdom",
                            "Superkingdom" = "superkingdom",
                            "unselected" = "")
        cusChoices <-
            mainChoices[match(mainRank, mainChoices):(length(mainChoices)-1)]
        selectInput(
            ns("rankSelectCus"),
            label = h5("Select taxonomy rank:"),
            choices = as.list(cusChoices),
            selected = mainRank
        )
    })

    # get list of possible taxon names for each selected rank from
    # rankSelectCus
    getTaxaCus <- reactive({
        rankSelectCus <- input$rankSelectCus

        if (length(rankSelectCus) == 0) return()
        else {
            # load list of unsorted taxa
            Dt <- getTaxonomyMatrix(taxDB(), TRUE, inputTaxonID())
            # load list of taxon name
            nameList <- getNameList(taxDB())
            # get rank name from rankSelect
            rankName <- rankSelectCus

            choice <- data.frame(
                ncbiID = unlist(Dt[rankName]), stringsAsFactors = FALSE
            )
            choice <- merge(choice, nameList, by = "ncbiID", all = FALSE)
            return(choice)
        }
    })

    # render list of possible taxon names from getTaxaCus()
    output$taxaSelectCus.ui <- renderUI({
        ns <- session$ns
        choice <- getTaxaCus()
        choice$fullName <- as.factor(choice$fullName)
        selectInput(
            ns("taxaSelectCus"),
            h5("Choose (super)taxon of interest:"),
            as.list(levels(choice$fullName)),
            levels(choice$fullName)[1]
        )
    })

    # get list of all taxa of the current taxonomy rank (from input$rankSelect)
    # based on selected supertaxon from input$taxaSelectCus
    cusTaxaName <- reactive({
        taxaSelectCus <- input$taxaSelectCus
        rankName <- input$rankSelectCus
        if (taxaSelectCus == "") return()

        # load list of unsorted taxa
        Dt <- getTaxonomyMatrix(taxDB(), TRUE, inputTaxonID())

        # get ID of selected (super)taxon from input$taxaSelectCus
        taxaList <- getNameList(taxDB())
        superID <- taxaList$ncbiID[taxaList$fullName == taxaSelectCus
                                    & taxaList$rank %in% c(rankName, "norank")]

        # from that ID, get list of all taxa for main selected taxon
        mainRankName <- rankSelect()
        customizedtaxaID <-
            levels(as.factor(Dt[mainRankName][Dt[[rankName]] == superID, ]))

        cusTaxaName <-
            taxaList$fullName[taxaList$ncbiID %in% customizedtaxaID]
        return(cusTaxaName)
    })
    return(cusTaxaName)
}
