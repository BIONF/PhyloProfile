#' search for NCBI taxonomy IDs for a list of taxon names
#' 10.05.2018
#' 

search_taxon_id_ui <- function(id){
  ns <- NS(id)
  
  tabPanel("Search for NCBI taxonomy IDs",
           column(3,
                  fileInput(ns("taxa_list"), h4("Upload taxa list")),
                  shinyBS::bsButton(ns("id_search"), "Search")
           ),
           column(9,
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
  # get ncbi taxa IDs ---------------------------------------------------------
  # retrieve ID for list of taxa names
  taxa_id <- reactive({
    if (input$id_search > 0){
      taxain <- input$taxa_list
      if (is.null(taxain)) return()
      
      taxa_name_df <- as.data.frame(read.table(file = taxain$datapath,
                                               sep = "\t",
                                               header = F,
                                               check.names = FALSE,
                                               comment.char = ""))
      
      id_df <- data.frame("name" = character(),
                          "new_name" = character(),
                          "id" = character(),
                          "type" = character(),
                          stringsAsFactors = FALSE)
      
      withProgress(message = "Retrieving IDs...", value = 0, {
        for (i in 1:nrow(taxa_name_df)){
          id <- get_uid(sciname = taxa_name_df[i, ])[1]
          if (is.na(id)){
            temp <- gnr_resolve(names = as.character(taxa_name_df[i, ]))
            if (nrow(temp) > 0){
              new_id <- get_uid(sciname = temp[1, 3])[1]
              if (is.na(new_id)){
                id_df[i, ] <- c(as.character(taxa_name_df[i, ]),
                                as.character(temp[1, 3]),
                                paste0("NA"), "notfound")
              } else {
                id_df[i, ] <- c(as.character(taxa_name_df[i, ]),
                                as.character(temp[1, 3]),
                                paste0("ncbi", new_id),
                                "notfound")
              }
            } else {
              id_df[i, ] <- c(as.character(taxa_name_df[i, ]),
                              paste0("no alternative"),
                              paste0("NA"),
                              "notfound")
            }
          } else {
            id_df[i, ] <- c(as.character(taxa_name_df[i, ]),
                            "NA",
                            paste0("ncbi", id),
                            "retrieved")
          }
          # Increment the progress bar, and update the detail text.
          incProgress(1 / nrow(taxa_name_df),
                      detail = paste(i, "/", nrow(taxa_name_df)))
        }
      })
      # return
      id_df
    }
  })
  
  # output retrieved taxa IDs -------------------------------------------------
  output$taxa_id <- renderDataTable(option = list(searching = FALSE), {
    if (input$id_search > 0){
      if (length(taxa_id()) > 0){
        tb <- as.data.frame(taxa_id())
        tb_filtered <- tb[tb$type == "retrieved", ]
        retrieved_dt <- tb_filtered[, c("name", "id")]
        colnames(retrieved_dt) <- c("Taxon_name", "Taxon_ID")
        retrieved_dt
      }
    }
  })
  
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
  
  # output mismatched taxa ----------------------------------------------------
  output$not_found_taxa <- renderDataTable(option = list(searching = FALSE), {
    if (input$id_search > 0){
      if (length(taxa_id()) > 0){
        tb <- as.data.frame(taxa_id())
        tb_filtered <- tb[tb$type == "notfound", ]
        not_found_dt <- tb_filtered[, c("name", "new_name", "id")]
        colnames(not_found_dt) <- c("Summitted name",
                                    "Alternative name",
                                    "Alternative ID")
        not_found_dt
      }
    }
  })
  
  # MISSING COMMENT -----------------------------------------------------------
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

