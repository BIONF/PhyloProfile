#' Gene age estimation plot
#' @param data for gene age plot (from reactive fn geneAgeDf)
#' @param geneAgeWidth plot width (from input$geneAgeWidth)
#' @param geneAgeHeight plot width (from input$geneAgeHeight)
#' @param geneAgeText text size (from input$geneAgeText)
#' @return list of genes of a selected age
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

plotGeneAgeUI <- function(id) {
    ns <- NS(id)
    tagList(
        column(
            2,
            downloadButton(ns("geneAgePlotDownload"),
                           "Download plot", class = "butDL"),
            tags$head(
                tags$style(HTML(
                    ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
            )
        ),
        column(
            10,
            uiOutput(ns("geneAge.ui")),
            br(),
            em(
                h6(
                    "01-Species; 02-Family; 03-Class; 04-Phylum; 05-Kingdom;
                    06-Superkingdom; 07-Last universal common ancestor;
                    Undef-Genes have been filtered out"
                )
            )
        ),
        hr(),
        column(
            4,
            downloadButton(ns("geneAgeTableDownload"), "Download gene list")
        ),
        column(
            8,
            tableOutput(ns("geneAge.table"))
        )
    )
}

plotGeneAge <- function(input, output, session,
                        data, geneAgeWidth, geneAgeHeight, geneAgeText) {
    # render gene age plot -----------------------------------------------------
    output$geneAgePlot <- renderPlot({
        if (is.null(data())) stop("No input data available!")
        createGeneAgePlot(geneAgePlotDf(data()), geneAgeText())
    })

    output$geneAge.ui <- renderUI({
        ns <- session$ns
        shinycssloaders::withSpinner(
            plotOutput(
                ns("geneAgePlot"),
                width = 600 * geneAgeWidth(),
                height = 150 * geneAgeHeight(),
                click = ns("plotClickGeneAge")
            )
        )
    })

    output$geneAgePlotDownload <- downloadHandler(
        filename = function() {
            "geneAgePlot.pdf"
        },
        content = function(file) {
            ggsave(
                file,
                plot = createGeneAgePlot(geneAgePlotDf(data()), geneAgeText()),
                width = 600 * geneAgeWidth() * 0.056458333,
                height = 150 * geneAgeHeight() * 0.056458333,
                units = "cm", dpi = 300, device = "pdf"
            )
        }
    )

    # render genAge.table based on clicked point on geneAgePlot ---------------
    selectedgeneAge <- reactive({
        if (is.null(data())) stop("No input data available!")
        selectedGene <- getSelectedGeneAge(data(), input$plotClickGeneAge$x)
        return(selectedGene)
    })

    output$geneAge.table <- renderTable({
        if (is.null(data())) stop("No input data available!")
        if (is.null(input$plotClickGeneAge$x)) return()

        data <- as.data.frame(selectedgeneAge())
        data$number <- rownames(data)
        colnames(data) <- c("geneID", "No.")
        data <- data[, c("No.", "geneID")]
        data
    })

    output$geneAgeTableDownload <- downloadHandler(
        filename = function(){
            c("selectedGeneList.out")
        },
        content = function(file){
            dataOut <- selectedgeneAge()
            write.table(dataOut, file,
                        sep = "\t",
                        row.names = FALSE,
                        quote = FALSE)
        }
    )

    return(selectedgeneAge)
}

#' get list of gene for a selected age
#' @param geneAgeDf data of estimated gene age (from fn "estimateGeneAge")
#' @param clickedX x coordinate of selected age

getSelectedGeneAge <- function(geneAgeDf, clickedX){
    data <- geneAgeDf

    # calculate the coordinate range for each age group
    rangeDf <- plyr::count(data, c("age"))

    rangeDf$percentage <- round(rangeDf$freq / sum(rangeDf$freq) * 100)
    rangeDf$rangeStart[1] <- 0
    rangeDf$rangeEnd[1] <- rangeDf$percentage[1]
    if (nrow(rangeDf) > 1) {
        for (i in 2:nrow(rangeDf)) {
            rangeDf$rangeStart[i] <- rangeDf$rangeEnd[i - 1] + 1
            rangeDf$rangeEnd[i] <-
                rangeDf$percentage[i] + rangeDf$rangeEnd[i - 1]
        }
    }

    # get list of gene for selected age
    if (is.null(clickedX)) return()
    else{
        corX <- 100 - round(clickedX)
        selectAge <- {
            as.character(rangeDf[rangeDf$rangeStart <= corX
                                 & rangeDf$rangeEnd >= corX, ]$age)
        }
        subData <- subset(data, age == selectAge)
        data <- data[data$age == selectAge, ]
    }

    # return list of genes
    geneList <- levels(as.factor(subData$geneID))
    return(geneList)
}
