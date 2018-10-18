#' Core gene identification
#'
#' @export
#' @param filtered_data full processed main data
#' (from reactive fn "get_data_filtered")
#' @param rank_select selected taxonomy rank (input$rank_select)
#' @param taxa_core selected list of taxa (input$taxa_core)
#' @param percent_core cutoff of percentage taxa present in a supertaxon
#' (input$percent_core)
#' @return list of core genes
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

identify_core_gene_ui <- function(id){
  ns <- NS(id)
  dataTableOutput(ns("core_gene.table"))
}

identify_core_gene <- function(input, output, session,
                               filtered_data,
                               rank_select, taxa_core, percent_core,
                               var1_cutoff, var2_cutoff,
                               core_coverage){

  output$core_gene.table <- renderDataTable({
    data <- core_geneDf()
    if (is.null(data)) return()
    else {
      data <- as.data.frame(data)
      data
    }
  })

  core_geneDf <- reactive({
    rankName <- substr(rank_select(), 4, nchar(rank_select()))

    # get ID list of chosen taxa
    taxa_list <- get_name_list(FALSE, FALSE)

    if ("none" %in% taxa_core()) superID <- NA
    else{
      superID <-
        taxa_list$ncbiID[taxa_list$fullName %in% taxa_core()
                         & taxa_list$rank %in% c(rankName, "norank")]
    }

    # get main input data
    mdData <- filtered_data()
    mdData <- mdData[, c("geneID",
                         "ncbiID",
                         "fullName",
                         "supertaxon",
                         "supertaxonID",
                         "rank",
                         "presSpec",
                         "mVar1",
                         "mVar2")]

    # filter by var1 and var2 cutoffs
    var1_cutoff_min <- var1_cutoff()[1]
    var1_cutoff_max <- var1_cutoff()[2]
    var2_cutoff_min <- var2_cutoff()[1]
    var2_cutoff_max <- var2_cutoff()[2]
    
    if (!is.null(var1_cutoff_max)) {
      if (!is.na(var1_cutoff_max)) {
        mdData <- subset(mdData, supertaxonID %in% superID
                       & mVar1 >= var1_cutoff_min)
        mdData <- subset(mdData, supertaxonID %in% superID
                       & mVar1 <= var1_cutoff_max)
      }
    }
    
    if (!is.null(var2_cutoff_max)) {
      if (!is.na(var2_cutoff_max)) {
        mdData <- subset(mdData, supertaxonID %in% superID
                       & mVar2 >= var2_cutoff_min)
        mdData <- subset(mdData, supertaxonID %in% superID
                       & mVar2 <= var2_cutoff_max)
      }
    }
    
    # filter by selecting taxa
    if (is.na(superID[1])) mdData <- NULL
    else {
      data <- subset(mdData, supertaxonID %in% superID
                     & presSpec >= percent_core()[1])
      data <- subset(data, supertaxonID %in% superID
                     & presSpec <= percent_core()[2])

      # get supertaxa present in each geneID
      supertaxonCount <- {
        as.data.frame(plyr::count(data,
                                  c("geneID", "supertaxonID")))
      }
      
      # count number of supertaxa present in each geneID
      # and get min number of supertaxa muss be taken into account
      count <- as.data.frame(table(supertaxonCount$geneID))
      require_coverage <- length(superID) * (core_coverage() / 100)
      
      # get only gene that contains orthologs in that coverage # of taxa
      core_gene <- subset(count, Freq >= require_coverage)
      core_gene$Var1 <- factor(core_gene$Var1)

      return(levels(core_gene$Var1))
    }
  })

  if (!is.null(core_geneDf)) return(core_geneDf)
}
