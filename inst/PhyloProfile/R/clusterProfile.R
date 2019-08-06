#' Profile clustering
#' @param distanceMatrix
#' @param clusterMethod Method to cluster the distances (input$clusterMethod)
#' @param plotWidth Width of the generated plot (input$clusterPlot.width)
#' @param plotHeight Height of the generated plot (input$clusterPlot.height)

source("R/functions.R")

clusterProfileUI <- function(id) {
    ns <- NS(id)
    tagList(
        column(
            8,
            downloadButton("downloadCluster",
                           "Download plot", class = "butDL"),
            tags$head(tags$style(
                HTML(".butDL{background-color:#476ba3;} .butDL{color: white;}")
            )),
            uiOutput(ns("cluster.ui"))
        ),
        column(
            4,
            downloadButton(ns("downloadDistanceMatrix"),
                           "Download distance matrix"),
            downloadButton(ns("downloadClusterGenes"),
                           "Download gene list"),
            tableOutput(ns("brushedCluster.table"))
        )
    )
}

clusterProfile <- function(input, output, session,
                            distanceMatrix,
                            clusterMethod,
                            plotWidth, plotHeight) {
    # Reactive function holding data for clustering ============================
    clusterData <- reactive({
        if (is.null(distanceMatrix)) stop("Distance matrix is NULL!")
        df <- clusterDataDend(distanceMatrix(), clusterMethod())
        return(df)
    })

    # Dendrogram ===============================================================
    output$dendrogram <- renderPlot({
        getDendrogram(clusterData())
    })

    output$cluster.ui <- renderUI({
        ns <- session$ns
        shinycssloaders::withSpinner(plotOutput(
            ns("dendrogram"),
            width = plotWidth(),
            height = plotHeight(),
            brush = brushOpts(
                id = ns("plotBrush"),
                delay = input$brushDelay,
                delayType = input$brushPolicy,
                direction = input$brushDir,
                resetOnNew = input$brushReset
            )
        ))
    })

    # download clustered plot ==================================================
    output$downloadCluster <- downloadHandler(
        filename = function() {
            "clusteredPlot.pdf"
        },
        content = function(file) {
            ggsave(
                file,
                plot = getDendrogram(clusterData()),
                dpi = 300,
                device = "pdf",
                limitsize = FALSE
            )
        }
    )

    # Brushed cluster table ====================================================
    #' render brushedCluster.table based on clicked point on dendrogram plot
    brushedClusterGene <- reactive({
		clusteredTree <- as.phylo(clusterData())
		labels <- rev(sortTaxaFromTree(clusteredTree))

        # get list of selected gene(s)
        if (is.null(input$plotBrush))
            return()
        else {
            top <- round(input$plotBrush$ymin)
            bottom <- round(input$plotBrush$ymax)
            df <- labels[bottom:top]
        }

        # return list of genes
        df <- df[!is.na(df)]
        return(df)
    })

    output$brushedCluster.table <- renderTable({
        if (is.null(input$plotBrush$ymin))
            return()

        data <- as.data.frame(brushedClusterGene())
        data$number <- rownames(data)
        colnames(data) <- c("geneID", "No.")
        data <- data[, c("No.", "geneID")]
        data
    })

    #' download gene list from brushedCluster.table
    output$downloadClusterGenes <- downloadHandler(
        filename = function() {
            c("selectedClusteredGeneList.out")
        },
        content = function(file) {
            dataOut <- brushedClusterGene()
            write.table(
                dataOut,
                file,
                sep = "\t",
                row.names = FALSE,
                quote = FALSE
            )
        }
    )

    #' download distance matrix
    output$downloadDistanceMatrix <- downloadHandler(
        filename = function() {
            c("distanceMatrixClustering.out")
        },
        content = function(file) {
            dataOut <- distanceMatrix()
            dataOut <- as.matrix(dataOut)
            write.table(
                dataOut,
                file,
                col.names = TRUE,
                row.names = TRUE,
                quote = FALSE,
                sep = " \t"
            )
        }
    )

    #' Return the brushed genes
    return(brushedClusterGene)
}
