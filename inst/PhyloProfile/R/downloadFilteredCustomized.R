#' Download filtered data from customized profile
#' @param data MAIN data for downloading (from module "downloadFilteredMain.R"
#' @param fasta fasta sequences (from reactive fn "customizedFastaDownload")
#' @param inSeq selected sequences in customized profile (from input$inSeq)
#' @param inTaxa selected taxa in customized profile (from input$inTaxa)
#' @return data of customized profile for downloading
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

downloadFilteredCustomizedUI <- function(id) {
    ns <- NS(id)

    tabPanel(
        "Customized data",
        column(
            12,
            strong(
                em(
                    "NOTE: Depending on your choice in Download filtered data
                    -> Main data, either all or only representative sequences
                    will be downloaded!"
                ),
                style = "color:red"
            ),
            hr()
        ),
        column(
            12,
            dataTableOutput(ns("filteredCustomData"))
        ),
        column(
            4,
            downloadButton(ns("downloadCustomData"),
                           "Download customized data")
        ),
        column(
            4,
            downloadButton(ns("downloadCustomFasta"),
                           "Download FASTA sequences"),
            uiOutput(ns("downloadCustomFasta.ui"))
        ),
        column(
            4,
            downloadButton(ns("downloadCustomLong"),
                           "Download data as PhyloProfile input format")
        )
    )
}

downloadFilteredCustomized <- function(input, output, session,
                                        data, fasta, inSeq, inTaxa){

    # filtered data for downloading (Customized Profile) -----------------------
    downloadCustomData <- reactive({
        data <- as.data.frame(data())

        # get subset of data according to selected genes/taxa
        if (!is.null(inSeq()) | !is.null(inTaxa())) {
            if (inSeq()[1] != "all" & inTaxa()[1] == "all") {
                # select data for selected sequences only
                customData <- subset(data, geneID %in% inSeq())
            } else if (inSeq()[1] == "all" & inTaxa()[1] != "all") {
                # select data for selected taxa only
                customData <- subset(data, supertaxon %in% inTaxa())
            } else if (inSeq()[1] != "all" & inTaxa()[1] != "all") {
                # select data for selected sequences and taxa
                customData <- subset(data, geneID %in% inSeq()
                                      & supertaxon %in% inTaxa())
            } else {
                customData <- data
            }
        } else {
            customData <- data
        }
        # return data
        customData <- as.matrix(customData)
        return(customData)
    })

    # download data ------------------------------------------------------------
    output$downloadCustomData <- downloadHandler(
        filename = function(){
            c("customFilteredData.out")
        },
        content = function(file){
            dataOut <- downloadCustomData()
            write.table(dataOut, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    # render download data table -----------------------------------------------
    output$filteredCustomData <- renderDataTable({
        data <- downloadCustomData()
        data
    })

    # download FASTA -----------------------------------------------------------
    output$downloadCustomFasta <- downloadHandler(
        filename = function(){
            c("customFilteredSeq.fa")
        },
        content = function(file){
            fastaOutDf <- fasta()
            write.table(fastaOutDf,
                        file,
                        sep = "\t",
                        col.names = FALSE,
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    # download data as long format ---------------------------------------------
    downloadCustomDataLong <- reactive({
        downloadCustomData <- downloadCustomData()

        if (ncol(downloadCustomData) == 6) {
            downloadCustomDataLong <- downloadCustomData[,c(1,4,2)]
        } else if (ncol(downloadCustomData) == 7) {
            downloadCustomDataLong <- downloadCustomData[,c(1,4,2,6)]
        } else if (ncol(downloadCustomData) == 8) {
            downloadCustomDataLong <- downloadCustomData[,c(1,4,2,6,7)]
        }

        return(downloadCustomDataLong)
    })

    output$downloadCustomLong <- downloadHandler(
        filename = function(){
            c("customFilteredData.phyloprofile")
        },
        content = function(file){
            dataOut <- downloadCustomDataLong()
            write.table(dataOut, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    return(downloadCustomData)
}
