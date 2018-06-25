#' Profile clustering
#'
#' @param data
#' @param dist_method
#' @param cluster_method
#' @param cluster_plot.width
#' @param cluster_plot.height
#' @param var1_aggregate_by

if (!require("bioDist")) install.packages("bioDist") # for the mutual information
if (!require("energy")) install.packages("energy") # for the mutual information, pearson

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
                            data,
                            dist_method, cluster_method,
                            cluster_plot.width, cluster_plot.height,
                            var1_aggregate_by){
  
  cluster_data <- reactive({
    df <- clusterDataDend(data(),
                          dist_method(),
                          cluster_method(),
                          var1_aggregate_by())

    return(df)
  })

  output$dendrogram <- renderPlot({
    if (is.null(data())) return()
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


#' cluster data ----------------------------------------------------------------
#' @export
#' @param data
#' @param dist_method
#' @param cluster_method
#' @param var1_aggregate_by
#' @return new data frame with % of present species
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

clusterDataDend <- function(data,
                            dist_method,
                            cluster_method,
                            var1_aggregate_by){
  # if (v$doPlot == FALSE) return()
  # dataframe for calculate distance matrix
  dat <- get_data_clustering(data, dist_method, var1_aggregate_by)
  dd.col <- as.dendrogram(hclust(get_distance_matrix(dat, dist_method),
                                 method = cluster_method))
  return(dd.col)
}

#' plot clustered profiles -----------------------------------------------------
#' @export
#' @param dd.col
#' @return plot clustered profiles
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

dendrogram <- function(dd.col){
  py <- as.ggdend(dd.col)
  p <- ggplot(py, horiz = TRUE, theme = theme_minimal()) +
    theme(axis.title = element_blank(), axis.text.y = element_blank())
  return(p) 
}


#' Get the distance matrix depending on the distance method---------------------
#' @export
#' @param profiles datafram containing phylogenetic profiles
#' @param dist_method distance method
#' @return distance matrix
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_distance_matrix <- function(profiles, method){
  dist_methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary")
  if (method %in% dist_methods) {
    distance_matrix <- dist(profiles, method = method)
  } else if (method %in% c("fisher", "distance_correlation")) {
    matrix <- data.frame()
    for (i in 1:nrow(profiles)) { # rows
      for (j in 1:nrow(profiles)) { # columns
        if (i == j) {
          matrix[i,i] = 1 # if this cell is NA as.dist does not work probably 
          break
        }
        if (method == "fisher") {
          contigency_table <- get_contengency_table(profiles[i,], profiles[j,])
          dist <- fisher.test(contigency_table)
        } else if (method == "distance_correlation") {
          dist <- dcor(unlist(profiles[i,]), unlist(profiles[j,]))
        }
        matrix[i,j] <- dist 
      }
    }
    profile_names <- rownames(profiles)
    colnames(matrix) <- profile_names[1:length(profile_names) - 1]
    rownames(matrix) <- profile_names
    distance_matrix <- as.dist(matrix)
  } else if (method == "mutual_information") {
    distance_matrix <- mutualInfo(as.matrix(profiles))
  } else if (method == "pearson") {
    distance_matrix <-  cor.dist(as.matrix(profiles))
  }
  return(distance_matrix)
}
