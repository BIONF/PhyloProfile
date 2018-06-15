#' Profile clustering
#' 
if (!require("bioDist")) install.packages("bioDist") # for the mutual information
if (!require("energy")) install.packages("energy") # for the mutual information, pearson

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
  dat <- get_data_clustering(data, dist_method)
  dd.col <- as.dendrogram(hclust(get_distance_matrix(dat, dist_method),
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

# get the phylogenetic profiles -----------------------------------------------
get_data_clustering <- function(data, dist_method){
  sub_data_heat <- subset(data, data$presSpec > 0)
  if (dist_method %in% c("mutual_information", "distance_correlation")){
    # Profiles with FAS scores
    sub_data_heat <- sub_data_heat[, c("geneID", "supertaxon", "var1")]
    sub_data_heat <- sub_data_heat[!duplicated(sub_data_heat), ]
    sub_data_heat <- {
      sub_data_heat[order(sub_data_heat$geneID, -abs(sub_data_heat$var1)),]}
    sub_data_heat <- {
      sub_data_heat[!duplicated(sub_data_heat[c("geneID", "supertaxon")]),]}
    
    wide_data <- spread(sub_data_heat, supertaxon, var1)
  }else {
    # Binary Profiles 
    sub_data_heat <- sub_data_heat[, c("geneID", "supertaxon", "presSpec")]
    sub_data_heat <- sub_data_heat[!duplicated(sub_data_heat), ]
    wide_data <- spread(sub_data_heat, supertaxon, presSpec)
  }
  dat <- wide_data[, 2:ncol(wide_data)]  # numerical columns
  rownames(dat) <- wide_data[, 1]
  dat[is.na(dat)] <- 0
  return(dat)
}

# Get the distance matrix depending on the distance method --------------------
get_distance_matrix <- function(profiles, method){
  dist_methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary")
  if(method %in% dist_methods){
    distance_matrix <- dist(profiles, method = method)
  } else if (method %in% c("fisher", "distance_correlation")){
    matrix <- data.frame()
    for(i in 1:nrow(profiles)){ # rows
      for(j in 1:nrow(profiles)){ # columns
        if (i == j){
          matrix[i,i] = 1 # if this cell is NA as.dist does not work probably 
          break
        }
        if(method == "fisher") {
          contigency_table <- get_table(profiles[i,], profiles[j,])
          dist <- fisher.test(contigency_table)
        } else if(method == "distance_correlation"){
          dist <- dcor(unlist(profiles[i,]), unlist(profiles[j,]))
        }
        matrix[i,j] <- dist 
      }
    }
    profile_names <- rownames(profiles)
    colnames(matrix) <- profile_names[1:length(profile_names)-1]
    rownames(matrix) <- profile_names
    distance_matrix <- as.dist(matrix)
  } else if (method == "mutual_information"){
    distance_matrix <- mutualInfo(as.matrix(profiles))
  } else if (method == "pearson"){
    distance_matrix <-  cor.dist(as.matrix(profiles))
  }
  return(distance_matrix)
}

# Calculate the contigency table for the fisher exact test --------------------
get_table <- function(profile_1, profile_2){
  contigency_table <- data.frame(c(0,0), c(0,0))
  for(i in 1:length(profile_1)){
    if(profile_1[i] == 1){
      if(profile_2[i] == 1) {
        contigency_table[1,1] <- contigency_table[1,1] + 1
      } else {
        contigency_table[2,1] <- contigency_table[2,1] + 1
      }
    } else{
      if(profile_2[i] == 1) {
        contigency_table[1,2] <- contigency_table[1,2] + 1
      } else {
        contigency_table[2,2] <- contigency_table[2,2] + 1
      }
    }
  }
  contigency_table
}

