#' Profile clustering
#' 

cluster_profile_ui <- function(id){
  ns <- NS(id)
  tagList(
    column(8,
           tags$head(
             tags$style(HTML(
               "#download_cluster{background-color:#A9E2F3}"))
           ),
           downloadButton("download_cluster", "Download plot"),
           uiOutput(ns("cluster.ui"))
    ),
    column(4,
           downloadButton(ns("download_cluster_genes"), "Download gene list"),
           tableOutput(ns("brushed_cluster.table"))
    )
  )
}

cluster_profile <- function(input, output, session,
                            data,
                            dist_method, cluster_method,
                            cluster_plot.width, cluster_plot.height){
  
  cluster_data <- reactive({
    df <- clusterDataDend(data(), dist_method(), cluster_method())
    return(df)
  }) 
  
  output$dendrogram <- renderPlot({
    # if (v$doPlot == FALSE) return()
    dendrogram(cluster_data())
  })
  
  output$cluster.ui <- renderUI({
    ns <- session$ns
    withSpinner(
      plotOutput(ns("dendrogram"),
                 width = cluster_plot.width(),
                 height = cluster_plot.height(),
                 brush = brushOpts(
                   id = ns("plot_brush"),
                   delay = input$brush_delay,
                   delayType = input$brush_policy,
                   direction = input$brush_dir,
                   resetOnNew = input$brush_reset)
      )
    )
  })
  
  # download clustered plot ---------------------------------------------------
  output$download_cluster <- downloadHandler(
    filename = function() {
      "clustered_plot.pdf"
    },
    content = function(file) {
      ggsave(file, plot = dendrogram(),
             dpi = 300, device = "pdf",
             limitsize = FALSE)
    }
  )
  
  # render brushed_cluster.table based on clicked point on dendrogram plot ----
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
  
  # download gene list from brushed_cluster.table -----------------------------
  output$download_cluster_genes <- downloadHandler(
    filename = function(){
      c("selectedClusteredGeneList.out")
    },
    content = function(file){
      data_out <- brushed_clusterGene()
      write.table(data_out, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  return(brushed_clusterGene)
}

# cluster data --------------------------------------------------------------
clusterDataDend <- function(data, dist_method, cluster_method){
  # if (v$doPlot == FALSE) return()
  # dataframe for calculate distance matrix
  dataHeat <- data
  
  sub_data_heat <- subset(dataHeat, dataHeat$presSpec > 0)
  sub_data_heat <- sub_data_heat[, c("geneID", "supertaxon", "presSpec")]
  sub_data_heat <- sub_data_heat[!duplicated(sub_data_heat), ]
  
  wide_data <- spread(sub_data_heat, supertaxon, presSpec)
  dat <- wide_data[, 2:ncol(wide_data)]  # numerical columns
  rownames(dat) <- wide_data[, 1]
  dat[is.na(dat)] <- 0
  
  dd.col <- as.dendrogram(hclust(dist(dat, method = dist_method),
                                 method = cluster_method))
  
  return(dd.col)
}

# plot clustered profiles -----------------------------------------------------
dendrogram <- function(dd.col){
  py <- as.ggdend(dd.col)
  p <- ggplot(py, horiz = TRUE, theme = theme_minimal()) +
    theme(axis.title = element_blank(), axis.text.y = element_blank())
  p
}
