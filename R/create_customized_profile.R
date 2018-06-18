#' Functions for creating heatmap & profile plots
#' 
source("R/create_profile_heatmap.R")
source("R/functions.R")

create_customized_profile_ui <- function(id) {
  ns <- NS(id)
  tagList(
    downloadButton(ns("selected_download"), "Download profile"),
    br(),
    uiOutput(ns("selected_plot.ui"))
  )
}

create_customized_profile <- function(input, output, session,
                                      data,
                                      clusteredDataHeat,
                                      apply_cluster,
                                      parameters,
                                      in_seq, in_taxa,
                                      rank_select, in_select,
                                      selected_width, selected_height,
                                      x_axis_selected) {
  # * data for heatmap
  dataHeat <- reactive({
    if (is.null(data()) | is.null(in_taxa()) | is.null(in_seq())) return()
    data_heat <- data_customized_plot(data(), in_taxa(), in_seq())
    # cluster dataHeat (if selected)
    if (apply_cluster() == TRUE) {
      data_heat <- data_customized_plot(clusteredDataHeat(), in_taxa(), in_seq())
    }
    return(data_heat)
  })
  
  # * plot customized profile -------------------------------------------
  output$selected_plot <- renderPlot({
    if (in_seq()[1] == "all" & in_taxa()[1] == "all") return()
    profile_plot(dataHeat(), parameters(), "none", rank_select(), "none")
  })
  
  output$selected_plot.ui <- renderUI({
    ns <- session$ns
    
    if (is.null(in_seq()[1]) | is.null(in_taxa()[1]))  return()
    else if (in_seq()[1] == "all" & in_taxa()[1] == "all") return()
    else{
      # if (input$auto_update_selected == FALSE){
      #   input$plot_custom
      #   isolate({
      #     withSpinner(
      #       plotOutput(ns("selected_plot"),
      #                  width = selected_width(),
      #                  height = selected_height(),
      #                  click = ns("plot_click_selected")
      #       )
      #     )
      #   })
      # } else {
        withSpinner(
          plotOutput(ns("selected_plot"),
                     width = selected_width(),
                     height = selected_height(),
                     click = ns("plot_click_selected")
          )
        )
      # }
    }
  })
  
  # * download customized plot ----------------------------------------------------
  output$selected_download <- downloadHandler(
    filename = function() {
      c("selected_plot.pdf")
    },
    content = function(file) {
      ggsave(file, plot = profile_plot(dataHeat(), parameters(), "none", rank_select(), "none"),#selected_plot(),
             width = selected_width() * 0.056458333,
             height = selected_height() * 0.056458333,
             units = "cm", dpi = 300, device = "pdf", limitsize = FALSE)
    }
  )
  
  # * get info of clicked point on customized profile ------------------------------
  selectedpoint_info <- reactive({
    # check input
    # if (vCt$doPlotCustom == FALSE) return()
    
    # get selected supertaxon name
    taxa_list <- get_name_list(FALSE, FALSE)
    rank_select <- rank_select()
    rankName <- substr(rank_select, 4, nchar(rank_select))
    in_select <- {
      as.numeric(taxa_list$ncbiID[taxa_list$fullName == in_select()])
    }
    
    dataHeat <- dataHeat()
    if (is.null(dataHeat)) return()
    
    # if (apply_cluster() == TRUE){
    #   dataHeat <- clusteredDataHeat()
    # }
    
    # get sub-dataframe of selected taxa and sequences
    dataHeat$supertaxonMod <- substr(dataHeat$supertaxon,
                                     6,
                                     nchar(as.character(dataHeat$supertaxon)))
    
    if (is.null(in_seq()[1]) | is.null(in_taxa()[1]))  return()
    if (in_taxa()[1] == "all" & in_seq()[1] != "all"){
      # select data from dataHeat for selected sequences only
      dataHeat <- subset(dataHeat, geneID %in% in_seq())
    } else if (in_seq()[1] == "all" & in_taxa()[1] != "all"){
      # select data from dataHeat for selected taxa only
      dataHeat <- subset(dataHeat, supertaxonMod %in% in_taxa())
    } else {
      # select data from dataHeat for selected sequences and taxa
      dataHeat <- subset(dataHeat, geneID %in% in_seq()
                         & supertaxonMod %in% in_taxa())
    }
    
    # drop all other supertaxon that are not in sub-dataframe
    dataHeat$supertaxon <- factor(dataHeat$supertaxon)
    dataHeat$geneID <- factor(dataHeat$geneID)
    
    # get values
    if (is.null(input$plot_click_selected$x)) return()
    else{
      # get cooridiate point
      if (x_axis_selected() == "genes"){
        corX <- round(input$plot_click_selected$y);
        corY <- round(input$plot_click_selected$x)
      } else {
        corX <- round(input$plot_click_selected$x);
        corY <- round(input$plot_click_selected$y)
      }
      
      # get geneID
      genes <- levels(dataHeat$geneID)
      geneID <- toString(genes[corY])
      # get supertaxon (spec)
      supertaxa <- levels(dataHeat$supertaxon)
      spec <- toString(supertaxa[corX])
      # get var1, percentage of present species and var2 score
      var1 <- NA
      if (!is.na(dataHeat$var1[dataHeat$geneID == geneID
                               & dataHeat$supertaxon == spec][1])){
        var1 <- max(na.omit(dataHeat$var1[dataHeat$geneID == geneID
                                          & dataHeat$supertaxon == spec]))
      }
      Percent <- NA
      if (!is.na(dataHeat$presSpec[dataHeat$geneID == geneID
                                   & dataHeat$supertaxon == spec][1])){
        Percent <- {
          max(na.omit(dataHeat$presSpec[dataHeat$geneID == geneID
                                        & dataHeat$supertaxon == spec]))
        }
      }
      var2 <- NA
      if (!is.na(dataHeat$var2[dataHeat$geneID == geneID
                               & dataHeat$supertaxon == spec][1])){
        var2 <- {
          max(na.omit(dataHeat$var2[dataHeat$geneID == geneID
                                    & dataHeat$supertaxon == spec]))
        }
      }
      
      # get ortholog ID
      orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID
                                  & dataHeat$supertaxon == spec]
      if (length(orthoID) > 1){
        orthoID <- paste0(orthoID[1], ",...")
      }
      
      if (is.na(as.numeric(Percent))) return()
      else{
        info <- c(geneID,
                  as.character(orthoID),
                  as.character(spec),
                  round(as.numeric(var1), 2),
                  round(as.numeric(Percent), 2),
                  round(as.numeric(var2), 2))
        return(info)
      }
    }
  })
  
  return(selectedpoint_info)
}

# create data for customized profile -------------------------------------------
data_customized_plot <- function(data_heat, in_taxa, in_seq){
  # process data
  data_heat$supertaxonMod <- {
    substr(data_heat$supertaxon,
           6,
           nchar(as.character(data_heat$supertaxon)))
  }
  
  if (in_taxa[1] == "all" & in_seq[1] != "all") {
    # select data from dataHeat for selected sequences only
    data_heat <- subset(data_heat, geneID %in% in_seq)
  } else if (in_seq[1] == "all" & in_taxa[1] != "all") {
    # select data from dataHeat for selected taxa only
    data_heat <- subset(data_heat, supertaxonMod %in% in_taxa)
  } else {
    # select data from dataHeat for selected sequences and taxa
    data_heat <- subset(data_heat,
                        geneID %in% in_seq
                        & supertaxonMod %in% in_taxa)
  }
  
  # remove unneeded dots
  data_heat$presSpec[data_heat$presSpec == 0] <- NA
  data_heat$paralog[data_heat$presSpec < 1] <- NA
  data_heat$paralog[data_heat$paralog == 1] <- NA
  
  return(data_heat)
}