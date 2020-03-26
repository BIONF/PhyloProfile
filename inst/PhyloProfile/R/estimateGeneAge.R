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
            8,
            downloadButton(ns("geneAgePlotDownload"),
                           "Download plot", class = "butDL"),
            tags$head(
                tags$style(HTML(
                    ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
            ),
            uiOutput(ns("geneAge.ui"))
        ),
        column(1),
        column(
            3,
            downloadButton(ns("geneAgeTableDownload"), "Download gene list"),
            br(), br(),
            DT::dataTableOutput(ns("geneAge.table"))
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
        # shinycssloaders::withSpinner(
            plotOutput(
                ns("geneAgePlot"),
                width = 800 * geneAgeWidth(),
                height = 300 * geneAgeHeight(),
                click = ns("plotClickGeneAge")
            )
        # )
    })

    output$geneAgePlotDownload <- downloadHandler(
        filename = function() {
            "geneAgePlot.pdf"
        },
        content = function(file) {
            ggsave(
                file,
                plot = createGeneAgePlot(geneAgePlotDf(data()), geneAgeText()),
                width = 800 * geneAgeWidth() * 0.056458333,
                height = 300 * geneAgeHeight() * 0.056458333,
                units = "cm", dpi = 300, device = "pdf"
            )
        }
    )

    # render genAge.table based on clicked point on geneAgePlot ---------------
    selectedgeneAge <- reactive({
        if (is.null(data())) stop("No input data available!")
        selectedGene <- getSelectedGeneAge(data(), input$plotClickGeneAge$y)
        return(selectedGene)
    })
    
    output$geneAge.table <- DT::renderDataTable(
        options = list(searching = FALSE, pageLength = 10
    ), {
        if (is.null(data())) stop("No input data available!")
        if (is.null(input$plotClickGeneAge$x)) return()
        data <- as.data.frame(selectedgeneAge())
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

getSelectedGeneAge <- function(geneAgeDf, clickedY){
    data <- geneAgeDf
    # calculate the coordinate range for each age group
    rangeDf <- plyr::count(data, c("age"))
    rangeDf <- rangeDf[seq(dim(rangeDf)[1],1),]
    
    if (is.null(clickedY)) return()
    else {
        corY <- round(clickedY)
        selectAge <- as.character(rangeDf$age[corY])
        subData <- subset(data, age == selectAge)
        data <- data[data$age == selectAge, ]
        geneList <- list(levels(as.factor(subData$geneID)))
        names(geneList) <- substring(selectAge, 4)
        return(geneList)
    }
}
