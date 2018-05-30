#' Functions for downloading filtered data from main profile
#' 

download_filtered_main_ui <- function(id){
  ns <- NS(id)
  
  tabPanel(
    "Main data",
    column(
      4,
      checkboxInput(
        ns("get_representative_main"),
        strong(em("Download representative sequences")),
        value = FALSE,
        width = NULL
      )
    ),
    column(
      4,
      conditionalPanel(
        condition = {
          sprintf("input['%s'] == true", ns("get_representative_main"))
        },
        uiOutput(ns("ref_var_main.ui"))
      )
    ),
    column(
      4,
      conditionalPanel(
        condition = {
          sprintf("input['%s'] == true", ns("get_representative_main"))
        },
        radioButtons(
          inputId = ns("ref_type_main"),
          label = {
            "Select representative by"
          },
          choices = list("max", "min"),
          selected = "max",
          inline = T
        )
      )
    ),
    column(
      12,
      dataTableOutput(ns("filtered_main_data"))
    ),
    column(
      3,
      downloadButton(ns("download_data"),
                     "Download filtered data")
    ),
    column(
      9,
      downloadButton(ns("download_fasta"),
                     "Download FASTA sequences")
      # ,
      # uiOutput(ns("download_fasta.ui"))
    )
  )
}


# render variable used for identifying representative genes -----------------
download_filtered_main <- function(input, output, session, data, fasta,
                                   var1_id, var2_id,
                                   var1, var2, percent){
  
  output$ref_var_main.ui <- renderUI({
    ns <- session$ns
    if (nchar(var2_id()) < 1 & nchar(var1_id()) < 1){
      radioButtons(inputId = ns("ref_var_main"), label = "Reference variable",
                   choices = list(var1_id(), var2_id()),
                   selected = var1_id())
    } else if (nchar(var2_id()) < 1){
      radioButtons(inputId = ns("ref_var_main"),
                   label = "Reference variable",
                   choices = list(var1_id()),
                   selected = var1_id())
    } else {
      radioButtons(inputId = ns("ref_var_main"),
                   label = "Reference variable",
                   choices = list(var1_id(), var2_id()),
                   selected = var1_id())
    }
  })
  
  # filtered data for downloading (Main Profile ) -----------------------------
  download_data <- reactive({
    ### filtered data
    data_out <- data() # <---- get_data_filtered()
    data_out <- as.data.frame(data_out[data_out$presSpec > 0, ])
    data_out <- data_out[!is.na(data_out$geneID), ]
    
    data_out <- as.data.frame(data_out[data_out$presSpec >= percent()[1], ])
    data_out <- as.data.frame(data_out[data_out$var1 >= var1()[1]
                                       & data_out$var1 <= var1()[2], ])
    
    if (!all(is.na(data_out$var2))){
      data_out <- as.data.frame(data_out[data_out$var2 >= var2()[1]
                                         & data_out$var2 <= var2()[2], ])
    } else {
      data_out$var2 <- 0
    }

    ### select only representative genes if chosen
    if (input$get_representative_main == TRUE){
      if (is.null(input$ref_var_main)) return()
      else{
        if (input$ref_var_main == var1_id()){
          data_out_agg <- aggregate(as.numeric(data_out$var1),
                                    by = list(data_out$geneID, data_out$ncbiID),
                                    FUN = input$ref_type_main)
        } else if (input$ref_var_main == var2_id()){
          data_out_agg <- aggregate(as.numeric(data_out$var2),
                                    by = list(data_out$geneID, data_out$ncbiID),
                                    FUN = input$ref_type_main)
        } else {
          data_out_agg <- data_out[data_out, c("geneID", "ncbiID", "var1")]
        }
        colnames(data_out_agg) <- c("geneID", "ncbiID", "var_best")
        
        data_out_representative <- merge(data_out, data_out_agg,
                                         by = c("geneID", "ncbiID"),
                                         all.x = TRUE)
        
        if (input$ref_var_main == var1_id()){
          data_out <- {
            data_out_representative[data_out_representative$var1 == data_out_representative$var_best, ]
          }
        } else if (input$ref_var_main == var2_id()){
          data_out <- {
            data_out_representative[data_out_representative$var2 == data_out_representative$var_best, ]
          }
        } else {
          data_out <- data_out
        }
        # used to select only one ortholog,
        # if there exist more than one "representative"
        data_out$dup <- paste0(data_out$geneID, "#", data_out$ncbiID)
        data_out <- data_out[!duplicated(c(data_out$dup)), ]
      }
    }
    
    # sub select columns of dataout
    data_out <- data_out[, c("geneID",
                             "orthoID",
                             "fullName",
                             "ncbiID",
                             "supertaxon",
                             "var1",
                             "var2",
                             "presSpec")] #,"numberSpec"
    data_out <- data_out[order(data_out$geneID, data_out$supertaxon), ]
    data_out <- data_out[complete.cases(data_out), ]
    
    data_out$geneID <- as.character(data_out$geneID)
    data_out$fullName <- as.character(data_out$fullName)
    data_out$ncbiID <- substr(data_out$ncbiID,
                              5,
                              nchar(as.character(data_out$ncbiID)))
    data_out$supertaxon <- substr(data_out$supertaxon,
                                  6,
                                  nchar(as.character(data_out$supertaxon)))
    data_out$var1 <- as.character(data_out$var1)
    data_out$var2 <- as.character(data_out$var2)
    # data_out$numberSpec <- as.numeric(data_out$numberSpec)
    data_out$presSpec <- as.numeric(data_out$presSpec)
    
    # rename columns
    names(data_out)[names(data_out) == "presSpec"] <- "%Spec"
    # names(data_out)[names(data_out)=="numberSpec"] <- "totalSpec"
    if (nchar(var1_id()) > 0){
      names(data_out)[names(data_out) == "var1"] <- var1_id()
    } else {
      data_out <- subset(data_out, select = -c(var1) )
    }
    if (nchar(var2_id()) > 0){
      names(data_out)[names(data_out) == "var2"] <- var2_id()
    } else {
      data_out <- subset(data_out, select = -c(var2) )
    }
    
    # return data for downloading
    data_out <- as.matrix(data_out)
    return(data_out)
  })
  
  # download data -------------------------------------------------------------
  output$download_data <- downloadHandler(
    filename = function(){
      c("filteredData.out")
    },
    content = function(file){
      data_out <- download_data()
      write.table(data_out, file,
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)
    }
  )

  # data table ui tab ---------------------------------------------------------
  output$filtered_main_data <- renderDataTable(rownames = FALSE, {
    
    # if (is.null(data())) return()
    data <- download_data()
    data
  })
  
  # download FASTA ------------------------------------------------------------
  output$download_fasta <- downloadHandler(
    filename = function(){
      c("filteredSeq.fa")
    },
    content = function(file){
      fasta_out_df <- fasta()
      write.table(fasta_out_df, file,
                  sep = "\t",
                  col.names = FALSE,
                  row.names = FALSE,
                  quote = FALSE)
    }
  )
  
  return(download_data)
}