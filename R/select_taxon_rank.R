#' Popup windows for (sub-)selecting rank & (super)taxon of interest
#'
#' @export
#' @param rank_select initial selected taxonomy rank (input$rank_select)
#' @param subset_taxa list of all input taxon IDs (from reactive fn subset_taxa)
#' @return list of all taxa in the initial taxonomy rank based on selected
#' supertaxon from input$taxa_select_cus
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

select_taxon_rank_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("rank_select_cus.ui")),
    uiOutput(ns("taxa_select_cus.ui"))
  )
}

select_taxon_rank <- function(input, output, session,
                              rank_select, subset_taxa) {

  # render list of available taxonomy ranks
  # (where the lowest rank is the same as the chosen main rank)
  output$rank_select_cus.ui <- renderUI({
    ns <- session$ns
    mainRank <- rank_select()
    mainChoices <- get_taxonomy_ranks()
    cusChoices <- mainChoices[mainChoices >= mainRank]
    selectInput(
      ns("rank_select_cus"),
      label = h5("Select taxonomy rank:"),
      choices = as.list(cusChoices),
      selected = mainRank
    )
  })

  # get list of possible taxon names for each selected rank from rank_select_cus
  get_taxa_cus <- reactive({
    rank_select_cus <- input$rank_select_cus

    if (length(rank_select_cus) == 0) return()
    else {
      # load list of unsorted taxa
      Dt <- get_taxa_list(TRUE, subset_taxa())
      # load list of taxon name
      nameList <- get_name_list(TRUE, FALSE)
      # get rank name from rank_select
      rankName <- substr(rank_select_cus, 4, nchar(rank_select_cus))
      
      choice <- as.data.frame
      choice <- rbind(Dt[rankName])
      colnames(choice) <- "ncbiID"
      choice <- merge(choice, nameList, by = "ncbiID", all = FALSE)
      return(choice)
    }
  })
  
  # render list of possible taxon names from get_taxa_cus()
  output$taxa_select_cus.ui <- renderUI({
    ns <- session$ns
    choice <- get_taxa_cus()
    choice$fullName <- as.factor(choice$fullName)
    selectInput(
      ns("taxa_select_cus"),
      h5("Choose (super)taxon of interest:"),
      as.list(levels(choice$fullName)),
      levels(choice$fullName)[1]
    )
  })

  # get list of all taxa of the current taxonomy rank (from input$rank_select)
  # based on selected supertaxon from input$taxa_select_cus
  cus_taxaName <- reactive({

    taxa_select_cus <- input$taxa_select_cus
    rankName <- substr(input$rank_select_cus, 4, nchar(input$rank_select_cus))
    if (taxa_select_cus == "") return()

    # load list of unsorted taxa
    Dt <- get_taxa_list(TRUE, subset_taxa())

    # get ID of selected (super)taxon from input$taxa_select_cus
    taxa_list <- get_name_list(FALSE, FALSE)
    superID <- taxa_list$ncbiID[taxa_list$fullName == taxa_select_cus
                                & taxa_list$rank %in% c(rankName, "norank")]

    # from that ID, get list of all taxa for main selected taxon
    mainRankName <- substr(rank_select(), 4, nchar(rank_select()))
    customizedtaxa_id <-
      levels(as.factor(Dt[mainRankName][Dt[rankName] == superID, ]))

    cus_taxaName <-
      taxa_list$fullName[taxa_list$rank %in% c(mainRankName, "norank")
                         & taxa_list$ncbiID %in% customizedtaxa_id]
    
    return(cus_taxaName)
  })

  return(cus_taxaName)
}
