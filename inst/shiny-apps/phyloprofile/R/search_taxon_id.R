#' Search NCBI taxonomy IDs for a list of taxon names
#'
#' @export
#' @param
#' @return
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

search_taxon_id_ui <- function(id){
  ns <- NS(id)

  tabPanel(
    "Search for NCBI taxonomy IDs",
    column(
      3,
      fileInput(ns("taxa_list"), h4("Upload taxa list")),
      checkboxInput(
        ns("from_online"),
        strong("DO NOT search taxonomy IDs using online NCBI database",
               style = "color:red"),
        value = FALSE
      ),
      shinyBS::bsButton(ns("id_search"), "Search")
    ),
    column(
      9,
      h4("Mismatch(es):"),
      dataTableOutput(ns("not_found_taxa")),
      downloadButton(ns("download_not_found_taxa"), "Download"),

      hr(),
      h4("Retrieved taxonomy ID(s):"),
      dataTableOutput(ns("taxa_id")),
      downloadButton(ns("download_taxa_id"), "Download")
    )
  )
}

search_taxon_id <- function(input, output, session){
  # retrieve ID for list of taxa names
  taxa_id <- reactive({
    if (input$id_search > 0) {
      taxain <- input$taxa_list
      if (is.null(taxain)) return()

      taxa_name_df <- as.data.frame(read.table(file = taxain$datapath,
                                               sep = "\t",
                                               header = FALSE,
                                               check.names = FALSE,
                                               comment.char = ""))
      colnames(taxa_name_df) <- c("name")
      
      tax_df <- as.data.frame(read.table("data/taxonomyMatrix.txt",
                                         sep = "\t",
                                         header = T,
                                         stringsAsFactors = T))
      id_df <- tax_df[tax_df$fullName %in% taxa_name_df$name,
                      c("abbrName","fullName")]
      colnames(id_df) <- c("id","name")
      
      id_df <- merge(taxa_name_df, id_df, all.x = TRUE)
      
      id_df$type <- "retrieved"
      id_df$type[is.na(id_df$id)] <- "notfound"
      id_df$new_name <- id_df$name
      
      notfound_df <- id_df[is.na(id_df$id),]
      if (nrow(notfound_df) > 0) {
        if (input$from_online == FALSE) {
          withProgress(
            message = "Retrieving IDs for unknown taxa...", 
            value = 0, {
              for (i in 1:nrow(notfound_df)) {
                id_df_tmp <- search_ncbi_online(
                  as.character(notfound_df[i,]$name)
                )
                id_df <- rbind(id_df, id_df_tmp)
              }
              # Increment the progress bar, and update the detail text.
              incProgress(1 / nrow(notfound_df),
                          detail = paste(i, "/", nrow(notfound_df)))
            }
          )
          id_df <- id_df[!is.na(id_df$id),]
        }
      }
    
      # return
      return(id_df)
    }
  })

  # output retrieved taxa IDs
  output$taxa_id <- renderDataTable(option = list(searching = FALSE), {
    if (input$id_search > 0) {
      if (length(taxa_id()) > 0) {
        tb <- as.data.frame(taxa_id())
        tb_filtered <- tb[tb$type == "retrieved", ]
        retrieved_dt <- tb_filtered[, c("name", "id")]
        colnames(retrieved_dt) <- c("Taxon_name", "Taxon_ID")
        retrieved_dt
      }
    }
  }, rownames = FALSE)

  # download retrieved taxa IDs
  output$download_taxa_id <- downloadHandler(
    filename = function(){
      c("retrievedtaxa_id.txt")
    },
    content = function(file){
      tb <- as.data.frame(taxa_id())
      tb_filtered <- tb[tb$type == "retrieved", ]
      retrieved_dt <- tb_filtered[, c("name", "id")]
      colnames(retrieved_dt) <- c("Taxon name", "Taxon ID")

      write.table(retrieved_dt, file,
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)
    }
  )

  # output mismatched taxa
  output$not_found_taxa <- renderDataTable(option = list(searching = FALSE), {
    if (input$id_search > 0) {
      if (length(taxa_id()) > 0) {
        tb <- as.data.frame(taxa_id())
        tb_filtered <- tb[tb$type == "notfound", ]
        not_found_dt <- tb_filtered[, c("name", "new_name", "id")]
        colnames(not_found_dt) <- c("Summitted name",
                                    "Alternative name",
                                    "Alternative ID")
        not_found_dt
      }
    }
  }, rownames = FALSE)

  # download mismatched taxa
  output$download_not_found_taxa <- downloadHandler(
    filename = function(){
      c("mismatchedTaxa.txt")
    },
    content = function(file){
      tb <- as.data.frame(taxa_id())
      tb_filtered <- tb[tb$type == "notfound", ]
      not_found_dt <- tb_filtered[, c("name", "new_name", "id")]
      colnames(not_found_dt) <- c("Summitted name",
                                  "Alternative name",
                                  "Alternative ID")

      write.table(not_found_dt, file,
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)
    }
  )
}

search_ncbi_online <- function(tax_name){
  id <- get_uid(sciname = tax_name)[1]
  
  id_df <- data.frame("name" = character(),
                      "new_name" = character(),
                      "id" = character(),
                      "type" = character(),
                      stringsAsFactors = FALSE)
  
  if (is.na(id)) {
    temp <- gnr_resolve(names = as.character(tax_name))
    if (nrow(temp) > 0) {
      new_id <- get_uid(sciname = temp[1, 3])[1]
      if (is.na(new_id)) {
        id_df[1, ] <- c(as.character(tax_name),
                        as.character(temp[1, 3]),
                        paste0("NA"), "notfound")
      } else {
        id_df[1, ] <- c(as.character(tax_name),
                        as.character(temp[1, 3]),
                        paste0("ncbi", new_id),
                        "notfound")
      }
    } else {
      id_df[1, ] <- c(as.character(tax_name),
                      paste0("no alternative"),
                      paste0("NA"),
                      "notfound")
    }
  } else {
    id_df[1, ] <- c(as.character(tax_name),
                    "NA",
                    paste0("ncbi", id),
                    "retrieved")
  }
  return(id_df)
}
