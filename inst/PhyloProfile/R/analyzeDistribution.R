#' Distribution plots
#' @param data data for plotting (from reactive fn "presSpecAllDt")
#' @param varID name of variable (either input$var1ID, input$var2ID or
#' "% present taxa"; from input$selectedDist)
#' @param varType type of variable (either var1, var2 or presSpec)
#' @param percent percentage cutoff (from input$percent)
#' @param distTextSize text size of distribution plot
#' @param distWidth width of distribution plot
#' (from input$distTextSize)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

analyzeDistributionUI <- function(id) {
    ns <- NS(id)
    tagList(
        column(
            2,
            downloadButton(ns("plotDownloadDist"), "Download plot")
        ),
        column(
            10,
            uiOutput(ns("distPlot.ui"))
        )
    )
}

analyzeDistribution <- function(input, output, session,
                                 data,
                                 varID, varType,
                                 percent,
                                 distTextSize, distWidth){

    # render distPlot.ui ------------------------------------------------------
    output$distPlot.ui <- renderUI({
        ns <- session$ns
        # shinycssloaders::withSpinner(
            plotOutput(ns("distributionPlot"),  width = distWidth())
        # )
    })

    output$distributionPlot <- renderPlot(width = distWidth(), height = 356, {
        createVarDistPlot(data(), varID(), varType(), percent(), distTextSize())
    })

    output$plotDownloadDist <- downloadHandler(
        filename = function() {
            paste0("distributionPlot.pdf")
        },
        content = function(file) {
            ggsave(
                file,
                plot = createVarDistPlot(
                    data(), varID(), varType(), percent(), distTextSize()
                ),
                width = distWidth() * 0.056458333,
                height = 356 * 0.056458333,
                unit = "cm",
                dpi = 300, device = "pdf", limitsize = FALSE)
        }
    )
}
