#' Functions for downloading filtered data from customized profile
#' 

download_filtered_customized_ui <- function(id){
  ns <- NS(id)
  
  tabPanel(
    "Customized data",
    column(
      12,
      strong(em("NOTE: Depending on your choice in [Download filtered data -> Main data], 
                  either all or only representative sequences will be downloaded!"), 
             style = "color:red"),
      hr()
    ),
    column(
      12,
      dataTableOutput(ns("filtered_custom_data"))
    ),
    column(
      3,
      downloadButton(ns("download_custom_data"),
                          "Download customized data")
    ),
    column(
      9,
      downloadButton(ns("download_custom_fasta"),
                     "Download FASTA sequences"),
      uiOutput(ns("download_custom_fasta.ui"))
    )
  )
}


# render variable used for identifying representative genes -----------------
download_filtered_customized <- function(input, output, session, data, fasta,
                                         in_seq, in_taxa){

  # filtered data for downloading (Customized Profile) ------------------------
  download_custom_data <- reactive({
    # check input
    # if (v$doPlot == FALSE) return()
    data <- as.data.frame(data()) # <--- download_data()
    
    # get subset of data according to selected genes/taxa
    if (!is.null(in_seq()) | !is.null(in_taxa())){
      if (in_seq()[1] != "all" & in_taxa()[1] == "all"){
        # select data for selected sequences only
        custom_data <- subset(data, geneID %in% in_seq())
      } else if (in_seq()[1] == "all" & in_taxa()[1] != "all"){
        # select data for selected taxa only
        custom_data <- subset(data, supertaxon %in% in_taxa())
      } else if (in_seq()[1] != "all" & in_taxa()[1] != "all") {
        # select data for selected sequences and taxa
        custom_data <- subset(data, geneID %in% in_seq()
                              & supertaxon %in% in_taxa())
      } else {
        custom_data <- data
      }
    } else {
      custom_data <- data
    }
    # return data
    custom_data <- as.matrix(custom_data)
    custom_data
  })
  
  # download data -------------------------------------------------------------
  output$download_custom_data <- downloadHandler(
    filename = function(){
      c("customFilteredData.out")
    },
    content = function(file){
      data_out <- download_custom_data()
      write.table(data_out, file,
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)
    }
  )
  
  # data table ui tab ---------------------------------------------------------
  output$filtered_custom_data <- renderDataTable(rownames = FALSE, {
    # if (v$doPlot == FALSE) return()
    data <- download_custom_data()
    data
  })
  
  # download FASTA ------------------------------------------------------------
  output$download_custom_fasta <- downloadHandler(
    filename = function(){
      c("customFilteredSeq.fa")
    },
    content = function(file){
      fasta_out_df <- fasta()
      write.table(fasta_out_df,
                  file,
                  sep = "\t",
                  col.names = FALSE,
                  row.names = FALSE,
                  quote = FALSE)
    }
  )
  
  return(download_custom_data)
}