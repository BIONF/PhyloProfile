#' Search NCBI taxonomy IDs for a list of taxon names
#'
#' @export
#' @param
#' @return
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

searchTaxonIDUI <- function(id){
    ns <- NS(id)

    tabPanel(
        "Search for NCBI taxonomy IDs",
        column(
            3,
            fileInput(ns("taxaList"), h4("Upload taxa list")),
            checkboxInput(
                ns("fromOnline"),
                strong("DO NOT search taxonomy IDs using online NCBI database",
                       style = "color:red"),
                value = FALSE
            ),
            shinyBS::bsButton(ns("idSearch"), "Search")
        ),
        column(
            9,
            h4("Mismatch(es):"),
            dataTableOutput(ns("notFoundTaxa")),
            downloadButton(ns("downloadNotFoundTaxa"), "Download"),

            hr(),
            h4("Retrieved taxonomy ID(s):"),
            dataTableOutput(ns("taxaID")),
            downloadButton(ns("downloadTaxaID"), "Download")
        )
    )
}

searchTaxonID <- function(input, output, session){
    # retrieve ID for list of taxa names
    taxaID <- reactive({
        if (input$idSearch > 0) {
            taxain <- input$taxaList
            if (is.null(taxain)) return()

            taxaNameDf <- as.data.frame(read.table(file = taxain$datapath,
                                                    sep = "\t",
                                                    header = FALSE,
                                                    check.names = FALSE,
                                                    comment.char = ""))
            colnames(taxaNameDf) <- c("name")

            taxDf <- as.data.frame(read.table("data/taxonomyMatrix.txt",
                                               sep = "\t",
                                               header = TRUE,
                                               stringsAsFactors = TRUE))
            idDf <- taxDf[taxDf$fullName %in% taxaNameDf$name,
                            c("abbrName","fullName")]
            colnames(idDf) <- c("id","name")

            idDf <- merge(taxaNameDf, idDf, all.x = TRUE)

            idDf$type <- "retrieved"
            idDf$type[is.na(idDf$id)] <- "notfound"
            idDf$newName <- idDf$name

            notfoundDf <- idDf[is.na(idDf$id),]
            if (nrow(notfoundDf) > 0) {
                if (input$fromOnline == FALSE) {
                    withProgress(
                        message = "Retrieving IDs for unknown taxa...",
                        value = 0, {
                            for (i in 1:nrow(notfoundDf)) {
                                idDfTmp <-
                                    PhyloProfile::searchTaxonIDOnline(
                                        as.character(notfoundDf[i,]$name)
                                    )
                                idDf <- rbind(idDf, idDfTmp)
                            }
                            # Increment the progress bar
                            incProgress(1 / nrow(notfoundDf),
                                        detail = paste(
                                            i, "/", nrow(notfoundDf))
                                        )
                        }
                    )
                    idDf <- idDf[!is.na(idDf$id),]
                }
            }

            # return
            return(idDf)
        }
    })

    # output retrieved taxa IDs
    output$taxaID <- renderDataTable(option = list(searching = FALSE), {
        if (input$idSearch > 0) {
            if (length(taxaID()) > 0) {
                tb <- as.data.frame(taxaID())
                tbFiltered <- tb[tb$type == "retrieved", ]
                retrievedDt <- tbFiltered[, c("name", "id")]
                colnames(retrievedDt) <- c("TaxonName", "TaxonID")
                retrievedDt
            }
        }
    }, rownames = FALSE)

    # download retrieved taxa IDs
    output$downloadTaxaID <- downloadHandler(
        filename = function(){
            c("retrievedtaxaID.txt")
        },
        content = function(file){
            tb <- as.data.frame(taxaID())
            tbFiltered <- tb[tb$type == "retrieved", ]
            retrievedDt <- tbFiltered[, c("name", "id")]
            colnames(retrievedDt) <- c("Taxon name", "Taxon ID")

            write.table(retrievedDt, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    # output mismatched taxa
    output$notFoundTaxa <- renderDataTable(option = list(searching = FALSE), {
        if (input$idSearch > 0) {
            if (length(taxaID()) > 0) {
                tb <- as.data.frame(taxaID())
                tbFiltered <- tb[tb$type == "notfound", ]
                notFoundDt <- tbFiltered[, c("name", "newName", "id")]
                colnames(notFoundDt) <- c("Summitted name",
                                            "Alternative name",
                                            "Alternative ID")
                notFoundDt
            }
        }
    }, rownames = FALSE)

    # download mismatched taxa
    output$downloadNotFoundTaxa <- downloadHandler(
        filename = function(){
            c("mismatchedTaxa.txt")
        },
        content = function(file){
            tb <- as.data.frame(taxaID())
            tbFiltered <- tb[tb$type == "notfound", ]
            notFoundDt <- tbFiltered[, c("name", "newName", "id")]
            colnames(notFoundDt) <- c("Summitted name",
                                        "Alternative name",
                                        "Alternative ID")

            write.table(notFoundDt, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )
}
