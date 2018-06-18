#' Core gene identification
#'

source("R/functions.R")

identify_core_gene_ui <- function(id){
  ns <- NS(id)

dataTableOutput(ns("core_gene.table"))
}

identify_core_gene <- function(input, output, session,
                               filtered_data,
                               rank_select, taxa_core, percent_core){

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
      superID <- {
        taxa_list$ncbiID[taxa_list$fullName %in% taxa_core()
                         & taxa_list$rank %in% c(rankName, "norank")]
      }
    }

    # get main input data
    mdData <- filtered_data() # <--- get_data_filtered()
    mdData <- mdData[, c("geneID",
                         "ncbiID",
                         "fullName",
                         "supertaxon",
                         "supertaxonID",
                         "rank",
                         "presSpec",
                         "mVar1",
                         "mVar2")]

    # filter by selecting taxa
    if (is.na(superID[1])) data <- NULL
    else{
      data <- subset(mdData, supertaxonID %in% superID
                     & presSpec >= percent_core())
      # get supertaxa present in each geneID
      supertaxonCount <- {
        as.data.frame(plyr::count(data,
                                  c("geneID", "supertaxonID")))
      }
      # count number of supertaxa present in each geneID
      # and get only gene that contains all chosen taxa
      count <- as.data.frame(table(supertaxonCount$geneID))
      core_gene <- subset(count, Freq == length(superID))
      core_gene$Var1 <- factor(core_gene$Var1)

      return(levels(core_gene$Var1))
    }
  })

  return(core_geneDf)
}
