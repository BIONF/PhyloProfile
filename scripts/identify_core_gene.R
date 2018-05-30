#' Core gene identification
#' 

source("scripts/functions.R")

identify_core_gene_ui <- function(id){
  ns <- NS(id)

  dataTableOutput(ns("cons_gene.table"))
}

identify_core_gene <- function(input, output, session, 
                               filtered_data,
                               rank_select, taxaCons, percent_cons){
  
  output$cons_gene.table <- renderDataTable({
    data <- cons_geneDf()
    if (is.null(data)) return()
    else {
      data <- as.data.frame(data)
      data
    }
  })
  
  cons_geneDf <- reactive({
    rankName <- substr(rank_select(), 4, nchar(rank_select()))
    
    # get ID list of chosen taxa
    taxa_list <- get_name_list(FALSE, FALSE)
    
    if ("none" %in% taxaCons()) superID <- NA
    else{
      superID <- {
        taxa_list$ncbiID[taxa_list$fullName %in% taxaCons()
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
                     & presSpec >= percent_cons())
      # get supertaxa present in each geneID
      supertaxonCount <- {
        as.data.frame(plyr::count(data,
                                  c("geneID", "supertaxonID")))
      }
      # count number of supertaxa present in each geneID
      # and get only gene that contains all chosen taxa
      count <- as.data.frame(table(supertaxonCount$geneID))
      cons_gene <- subset(count, Freq == length(superID))
      cons_gene$Var1 <- factor(cons_gene$Var1)
      
      return(levels(cons_gene$Var1))
    }
  })
  
  return(cons_geneDf)
}


