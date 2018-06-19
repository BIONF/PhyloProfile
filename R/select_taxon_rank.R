#' Windows for selecting rank & (super)taxon of interest
#'

source("R/functions.R")

select_taxon_rank_ui <- function(id){
  ns <- NS(id)
  tagList(
    uiOutput(ns("rank_select_cus")),
    uiOutput(ns("taxa_select_cus"))
  )
}

select_taxon_rank <- function(input, output, session,
                              rank_select, subset_taxa){

  # print list of available customized taxonomy ranks -------------------------
  # (the lowest rank is the same as the chosen main rank)
  output$rank_select_cus <- renderUI({
    ns <- session$ns
    mainRank <- rank_select()
    mainChoices <- get_taxonomy_ranks()
    cusChoices <- mainChoices[mainChoices >= mainRank]
    selectInput(ns("rank_select_cus"),
                label = h5("Select taxonomy rank:"),
                choices = as.list(cusChoices),
                selected = mainRank)
  })

  output$taxa_select_cus <- renderUI({
    ns <- session$ns
    choice <- taxa_select_cus()
    choice$fullName <- as.factor(choice$fullName)
    selectInput(ns("taxa_select_cus"),
                h5("Choose (super)taxon of interest:"),
                as.list(levels(choice$fullName)),
                levels(choice$fullName)[1])
  })

  # print list of available taxa for customized plot --------------------------
  # (based on rank from rank_select_cus)
  taxa_select_cus <- reactive({
    rank_select_cus <- input$rank_select_cus

    if (length(rank_select_cus) == 0) return()
    else{
      ### load list of unsorted taxa
      Dt <- get_taxa_list(TRUE, subset_taxa())

      ### load list of taxon name
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

  # get list of taxa based on selected taxa_select_cus ------------------------
  cus_taxaName <- reactive({

    taxa_select_cus <- input$taxa_select_cus
    rankName <- substr(input$rank_select_cus, 4, nchar(input$rank_select_cus))

    if (taxa_select_cus == "") return()

    # load list of unsorted taxa
    Dt <- get_taxa_list(TRUE, subset_taxa())

    # get ID of customized (super)taxon
    taxa_list <- get_name_list(FALSE, FALSE)
    superID <- taxa_list$ncbiID[taxa_list$fullName == taxa_select_cus
                                & taxa_list$rank %in% c(rankName, "norank")]

    # from that ID, get list of all taxa for main selected taxon
    mainRankName <- substr(rank_select(), 4, nchar(rank_select()))
    customizedtaxa_id <- {
      levels(as.factor(Dt[mainRankName][Dt[rankName] == superID, ]))
    }

    cus_taxaName <- {
      taxa_list$fullName[taxa_list$rank %in% c(mainRankName, "norank")
                         & taxa_list$ncbiID %in% customizedtaxa_id]
    }
    return(cus_taxaName)
  })

  return(cus_taxaName)
}
