#' Profile clustering
#' 
#' @param distance_matrix
#' @param plot_width
#' @param plot_height

source("R/functions.R")

cluster_profile_ui <- function(id){
  ns <- NS(id)
  tagList(
    column(8,
           downloadButton("download_cluster",
                          "Download plot", class = "butDL"),
           tags$head(
             tags$style(HTML(
               ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
           ),
           uiOutput(ns("cluster.ui"))
    ),
    column(4,
           downloadButton(ns("download_cluster_genes"), "Download gene list"),
           tableOutput(ns("brushed_cluster.table"))
    )
  )
}

cluster_profile <- function(input, output, session,
                            distance_matrix,
                            cluster_method,
                            plot_width, plot_height
                            ){
  # Reactive function holding data for clustering =========================
  cluster_data <- reactive({
    df <- clusterDataDend(distance_matrix(), cluster_method())
    return(df)
  })

  # Dendrogram =========================
  output$dendrogram <- renderPlot({
    if (is.null(data())) return()
    dendrogram(cluster_data())
  })

  output$cluster.ui <- renderUI({
    ns <- session$ns
    withSpinner(
      plotOutput(ns("dendrogram"),
                 width = plot_width(),
                 height = plot_height(),
                 brush = brushOpts(
                   id = ns("plot_brush"),
                   delay = input$brush_delay,
                   delayType = input$brush_policy,
                   direction = input$brush_dir,
                   resetOnNew = input$brush_reset)
      )
    )
  })

  # download clustered plot =========================
  output$download_cluster <- downloadHandler(
    filename = function() {
      "clustered_plot.pdf"
    },
    content = function(file) {
      ggsave(file, plot = dendrogram(cluster_data()),
             dpi = 300, device = "pdf",
             limitsize = FALSE)
    }
  )

  # Brushed cluster table =========================
  #' render brushed_cluster.table based on clicked point on dendrogram plot 
  brushed_clusterGene <- reactive({
    # if (v$doPlot == FALSE) return()

    dd.col <- cluster_data()
    dt <- dendro_data(dd.col)
    dt$labels$label <- levels(dt$labels$label)

    # get list of selected gene(s)
    if (is.null(input$plot_brush)) return()
    else{
      top <- as.numeric(-round(input$plot_brush$ymin))
      bottom <- as.numeric(-round(input$plot_brush$ymax))

      df <- dt$labels[bottom:top, ]
    }

    # return list of genes
    df <- df[complete.cases(df), 3]
  })

  output$brushed_cluster.table <- renderTable({
    if (is.null(input$plot_brush$ymin)) return()

    data <- as.data.frame(brushed_clusterGene())
    data$number <- rownames(data)
    colnames(data) <- c("geneID", "No.")
    data <- data[, c("No.", "geneID")]
    data
  })

  #' download gene list from brushed_cluster.table 
  output$download_cluster_genes <- downloadHandler(
    filename = function(){
      c("selectedClusteredGeneList.out")
    },
    content = function(file){
      data_out <- brushed_clusterGene()
      write.table(data_out, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  #' Return the brushed genes 
  return(brushed_clusterGene)
}


#' cluster data ---------------------------------------------------------------
#' @export
#' @param distance_matrix 
#' @return new data frame with % of present species
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
clusterDataDend <- function(distance_matrix, cluster_method){
  # if (v$doPlot == FALSE) return()
  if (is.null(distance_matrix)) return() 
  dd.col <- as.dendrogram(hclust(distance_matrix,
                                 method = cluster_method))
  return(dd.col)
}

#' plot clustered profiles ----------------------------------------------------
#' @export
#' @param dd.col
#' @return plot clustered profiles
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

dendrogram <- function(dd.col){
  if (is.null(dd.col)) return()
  py <- as.ggdend(dd.col)
  p <- ggplot(py, horiz = TRUE)
  #, theme = theme_minimal()) +
  #  theme(axis.title = element_blank(), axis.text.y = element_blank())
  return(p) 
}
