#' Download filtered data from main profile
#' @param full processed main data (from reactive fn "getDataFiltered")
#' @param fasta fasta sequences (from reactive fn "mainFastaDownload")
#' @param var1ID name of 1st variable (from input$var1ID)
#' @param var2ID name of 2nd variable (from input$var2ID)
#' @param var1 cutoff value of 1st variable (from input$var1)
#' @param var2 cutoff value of 2nd variable (from input$var2)
#' @param percent cutoff value of percentage species in each supertaxon
#' (from input$percent)
#' @return data of main profile for downloading
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

downloadFilteredMainUI <- function(id) {
    ns <- NS(id)

    tabPanel(
        "Main data",
        column(
            4,
            checkboxInput(
                ns("getRepresentativeMain"),
                strong(em("Download representative sequences")),
                value = FALSE,
                width = NULL
            )
        ),
        column(
            4,
            conditionalPanel(
                condition = {
                    sprintf("input['%s'] == true",
                            ns("getRepresentativeMain"))
                },
                uiOutput(ns("refVarMain.ui"))
            )
        ),
        column(
            4,
            conditionalPanel(
                condition = {
                    sprintf("input['%s'] == true",
                            ns("getRepresentativeMain"))
                },
                radioButtons(
                    inputId = ns("refTypeMain"),
                    label = {
                        "Select representative by"
                    },
                    choices = list("max", "min"),
                    selected = "max",
                    inline = TRUE
                )
            )
        ),
        column(
            12,
            DT::dataTableOutput(ns("filteredMainData"))
        ),
        column(
            4,
            downloadButton(ns("downloadData"),
                           "Download filtered data")
        ),
        column(
            4,
            downloadButton(ns("downloadFasta"),
                           "Download FASTA sequences")
        ),
        column(
            4,
            downloadButton(ns("downloadLong"),
                           "Download data as PhyloProfile input format")
        )
    )
}

downloadFilteredMain <- function(
    input, output, session, data, taxaCount, fasta, var1ID, var2ID, var1, var2, 
    percent
){

    # render options for downloading -------------------------------------------
    output$refVarMain.ui <- renderUI({
        ns <- session$ns
        if (nchar(var2ID()) < 1 & nchar(var1ID()) < 1) {
            radioButtons(
                inputId = ns("refVarMain"), label = "Reference variable",
                choices = list(var1ID(), var2ID()),
                selected = var1ID()
            )
        } else if (nchar(var2ID()) < 1) {
            radioButtons(
                inputId = ns("refVarMain"),
                label = "Reference variable",
                choices = list(var1ID()),
                selected = var1ID()
            )
        } else {
            radioButtons(
                inputId = ns("refVarMain"),
                label = "Reference variable",
                choices = list(var1ID(), var2ID()),
                selected = var1ID()
            )
        }
    })

    # filtered data for downloading (Main Profile ) ----------------------------
    downloadData <- reactive({
        if (is.null(data())) stop("Data for downloading is NULL!")
        ### filtered data
        dataOut <- data()
        if (length(var1()) == 1) {
            dataOut <- dataOut[dataOut$var1 >= var1()[1], ]
        } else {
            dataOut <- dataOut[
                dataOut$var1 >= var1()[1] & dataOut$var1 <= var1()[2], 
            ]
        }
        if (!all(is.na(dataOut$var2))) {
            if (length(var2()) == 1) {
                dataOut <- dataOut[dataOut$var2 >= var2()[1], ]
            } else {
                dataOut <- dataOut[
                    dataOut$var2 >= var2()[1] & dataOut$var2 <= var2()[2], 
                ]
            }
        } else {
            dataOut$var2 <- 0
        }
        # calculate presSpec
        finalPresSpecDt <- calcPresSpec(dataOut, taxaCount())
        dataOut <- Reduce(
            function(x, y) merge(x, y, by = c("geneID", "supertaxon"), all.x=TRUE),
            list(dataOut, finalPresSpecDt))
        # dataOut <- as.data.frame(dataOut[dataOut$presSpec > 0, ])
        dataOut <- dataOut[
            dataOut$presSpec >= percent()[1] & dataOut$presSpec <= percent()[2], 
        ]
        dataOut <- dataOut[!is.na(dataOut$geneID), ]
        ### select only representative genes if chosen
        if (input$getRepresentativeMain == TRUE) {
            if (is.null(input$refVarMain)) return()
            else {
                if (input$refVarMain == var1ID()) {
                    dataOutAgg <- aggregate(
                        dataOut$var1,
                        by = list(dataOut$geneID, dataOut$ncbiID),
                        FUN = input$refTypeMain
                    )
                } else if (input$refVarMain == var2ID()) {
                    dataOutAgg <- aggregate(
                        dataOut$var2,
                        by = list(dataOut$geneID, dataOut$ncbiID),
                        FUN = input$refTypeMain
                    )
                } else {
                    dataOutAgg <-
                        dataOut[dataOut, c("geneID", "ncbiID", "var1")]
                }
                colnames(dataOutAgg) <- c("geneID", "ncbiID", "varBest")

                dataOutRepresentative <- merge(dataOut, dataOutAgg,
                                                 by = c("geneID", "ncbiID"),
                                                 all.x = TRUE)

                if (input$refVarMain == var1ID()) {
                    dataOut <-
                        dataOutRepresentative[
                            dataOutRepresentative$var1 ==
                                dataOutRepresentative$varBest,
                            ]
                } else if (input$refVarMain == var2ID()) {
                    dataOut <-
                        dataOutRepresentative[
                            dataOutRepresentative$var2 ==
                                dataOutRepresentative$varBest,
                            ]
                } else {
                    dataOut <- dataOut
                }
                # used to select only one ortholog,
                # if there exist more than one "representative"
                dataOut$dup <- paste0(dataOut$geneID, "#", dataOut$ncbiID)
                dataOut <- dataOut[!duplicated(c(dataOut$dup)), ]
            }
        }

        # sub select columns of dataout
        dataOut <- dataOut[, c("geneID",
                                 "orthoID",
                                 "fullName",
                                 "ncbiID",
                                 "supertaxon",
                                 "var1",
                                 "var2",
                                 "presSpec")]
        dataOut <- dataOut[order(dataOut$geneID, dataOut$supertaxon), ]
        dataOut <- dataOut[complete.cases(dataOut), ]

        dataOut$geneID <- as.character(dataOut$geneID)
        dataOut$fullName <- as.character(dataOut$fullName)
        dataOut$ncbiID <- as.character(dataOut$ncbiID)
        dataOut$supertaxon <- substr(dataOut$supertaxon,
                                      6,
                                      nchar(as.character(dataOut$supertaxon)))
        dataOut$var1 <- as.character(dataOut$var1)
        dataOut$var2 <- as.character(dataOut$var2)
        dataOut$presSpec <- dataOut$presSpec

        # rename columns
        names(dataOut)[names(dataOut) == "presSpec"] <- "%Spec"
        if (nchar(var1ID()) > 0) {
            names(dataOut)[names(dataOut) == "var1"] <- var1ID()
        } else {
            dataOut <- subset(dataOut, select = -c(var1) )
        }
        if (nchar(var2ID()) > 0) {
            names(dataOut)[names(dataOut) == "var2"] <- var2ID()
        } else {
            dataOut <- subset(dataOut, select = -c(var2) )
        }

        # return data for downloading
        dataOut <- as.matrix(dataOut)
        return(dataOut)
    })

    # download data ------------------------------------------------------------
    output$downloadData <- downloadHandler(
        filename = function(){
            c("filteredData.out")
        },
        content = function(file){
            dataOut <- downloadData()
            write.table(dataOut, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    # render download data table -----------------------------------------------
    output$filteredMainData <- DT::renderDataTable({
        data <- downloadData()
        data
    })

    # download FASTA -----------------------------------------------------------
    output$downloadFasta <- downloadHandler(
        filename = function(){
            c("filteredSeq.fa")
        },
        content = function(file){
            fastaOutDf <- fasta()
            write.table(fastaOutDf, file,
                        sep = "\t",
                        col.names = FALSE,
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    # download data as long format ---------------------------------------------
    downloadDataLong <- reactive({
        downloadData <- downloadData()

        if (ncol(downloadData) == 6) {
            downloadDataLong <- downloadData[,c(1,4,2)]
        } else if (ncol(downloadData) == 7) {
            downloadDataLong <- downloadData[,c(1,4,2,6)]
        } else if (ncol(downloadData) == 8) {
            downloadDataLong <- downloadData[,c(1,4,2,6,7)]
        }

        return(downloadDataLong)
    })

    output$downloadLong <- downloadHandler(
        filename = function(){
            c("filteredData.phyloprofile")
        },
        content = function(file){
            dataOut <- downloadDataLong()
            write.table(dataOut, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    return(downloadData)
}
