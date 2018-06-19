#' Profile plot
#' 
#' @export
#' @param data data for heatmap plot (from reactive fn "dataHeat")
#' @param clusteredDataHeat clustered data (from reactive fn "clusteredDataHeat"
#' @param apply_cluster choose clustered data or not (from input$apply_cluster)
#' @param parameters plot parameters (colors, size, variable names, ...)
#' @param in_seq subset sequences for customized profile (input$in_seq)
#' @param in_taxa subset taxa for customized profile (input$in_taxa)
#' @param rank_select selected taxonomy rank (input$rank_select)
#' @param in_select selected taxon name (input$in_select)
#' @param taxon_highlight highlighted taxon (input$taxon_highlight)
#' @param gene_highlight highlighted gene (input$gene_highlight)
#' @param width plot width
#' @param height plot height
#' @param x_axis type of x_axis (either "genes" or "taxa", from input$x_axis)
#' @param type_profile either "main_profile" or "customized_profile"
#' @return info for selected point on the profile
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

create_profile_plot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    downloadButton(ns("profile_download"), "Download profile"),
    br(),
    uiOutput(ns("plot.ui"))
  )
}

create_profile_plot <- function(input, output, session,
                                data, clusteredDataHeat,
                                apply_cluster,
                                parameters,
                                in_seq, in_taxa,
                                rank_select, in_select,
                                taxon_highlight, gene_highlight,
                                width, height,
                                x_axis,
                                type_profile) {
  # data for heatmap -----------------------------------------------------------
  dataHeat <- reactive({
    if (is.null(data())) return()
    
    if (type_profile() == "customized_profile") {
      if (is.null(in_taxa()) | is.null(in_seq())) return()
      
      data_heat <- data_customized_plot(data(), in_taxa(), in_seq())
      if (apply_cluster() == TRUE) {
        data_heat <- data_customized_plot(clusteredDataHeat(),
                                          in_taxa(), in_seq())
      }
    } else {
      data_heat <- data_main_plot(data())
      if (apply_cluster() == TRUE) {
        data_heat <- data_main_plot(clusteredDataHeat())
      }
    }
    
    return(data_heat)
  })
  
  # render heatmap profile -----------------------------------------------------
  output$plot <- renderPlot({
    if (is.null(data())) return()
    if (type_profile() == "customized_profile") {
      if (in_seq()[1] == "all" & in_taxa()[1] == "all") return()
    }
    
    profile_plot(
      dataHeat(),
      parameters(),
      taxon_highlight(),
      rank_select(),
      gene_highlight()
    )
  })
  
  output$plot.ui <- renderUI({
    ns <- session$ns
    
    if (type_profile() == "customized_profile") {
      if (is.null(in_seq()[1]) | is.null(in_taxa()[1]))  return()
      else if (in_seq()[1] == "all" & in_taxa()[1] == "all") return()
    }
    
    withSpinner(
      plotOutput(
        ns("plot"),
        width = width(),
        height = height(),
        click = ns("plot_click")
      )
    )
  })
  
  output$profile_download <- downloadHandler(
    filename = function() {
      c("profile.pdf")
    },
    content = function(file) {
      ggsave(
        file,
        plot = profile_plot(dataHeat(),
                            parameters(),
                            "none",
                            rank_select(),
                            "none"),
        width = width() * 0.056458333,
        height = height() * 0.056458333,
        units = "cm", dpi = 300, device = "pdf", limitsize = FALSE
      )
    }
  )
  
  # get info of clicked point on heatmap plot ----------------------------------
  selectedpoint_info <- reactive({
    
    # get selected supertaxon name
    taxa_list <- get_name_list(FALSE, FALSE)
    rank_select <- rank_select()
    rankName <- substr(rank_select, 4, nchar(rank_select))
    in_select <- {
      as.numeric(taxa_list$ncbiID[taxa_list$fullName == in_select()])
    }
    
    dataHeat <- dataHeat()
    if (is.null(dataHeat)) return()
    
    if (type_profile() == "customized_profile") {
      # get sub-dataframe of selected taxa and sequences
      dataHeat$supertaxonMod <- substr(dataHeat$supertaxon,
                                       6,
                                       nchar(as.character(dataHeat$supertaxon)))
      
      if (is.null(in_seq()[1]) | is.null(in_taxa()[1]))  return()
      if (in_taxa()[1] == "all" & in_seq()[1] != "all") {
        # select data from dataHeat for selected sequences only
        dataHeat <- subset(dataHeat, geneID %in% in_seq())
      } else if (in_seq()[1] == "all" & in_taxa()[1] != "all") {
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
    }
    
    # get values
    if (is.null(input$plot_click$x)) return()
    else {
      # get cooridiate point
      if (x_axis() == "genes") {
        corX <- round(input$plot_click$y);
        corY <- round(input$plot_click$x)
      } else {
        corX <- round(input$plot_click$x);
        corY <- round(input$plot_click$y)
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
                               & dataHeat$supertaxon == spec][1])) {
        var1 <- max(na.omit(dataHeat$var1[dataHeat$geneID == geneID
                                          & dataHeat$supertaxon == spec]))
      }
      Percent <- NA
      if (!is.na(dataHeat$presSpec[dataHeat$geneID == geneID
                                   & dataHeat$supertaxon == spec][1])) {
        Percent <- {
          max(na.omit(dataHeat$presSpec[dataHeat$geneID == geneID
                                        & dataHeat$supertaxon == spec]))
        }
      }
      var2 <- NA
      if (!is.na(dataHeat$var2[dataHeat$geneID == geneID
                               & dataHeat$supertaxon == spec][1])) {
        var2 <- {
          max(na.omit(dataHeat$var2[dataHeat$geneID == geneID
                                    & dataHeat$supertaxon == spec]))
        }
      }
      
      # get ortholog ID
      orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID
                                  & dataHeat$supertaxon == spec]
      if (length(orthoID) > 1) {
        orthoID <- paste0(orthoID[1], ",...")
      }
      
      if (is.na(as.numeric(Percent))) return()
      else {
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


# create data for main profile -------------------------------------------------
data_main_plot <- function(data_heat){
  # reduce number of inparalogs based on filtered dataHeat
  data_heat_tb <- data.table(na.omit(data_heat))
  data_heat_tb[, paralogNew := .N, by = c("geneID", "supertaxon")]
  data_heat_tb <- data.frame(data_heat_tb[, c("geneID",
                                              "supertaxon",
                                              "paralogNew")])
  
  data_heat <- merge(data_heat, data_heat_tb,
                     by = c("geneID", "supertaxon"),
                     all.x = TRUE)
  data_heat$paralog <- data_heat$paralogNew
  data_heat <- data_heat[!duplicated(data_heat), ]
  
  # remove unneeded dots
  data_heat$presSpec[data_heat$presSpec == 0] <- NA
  data_heat$paralog[data_heat$presSpec < 1] <- NA
  data_heat$paralog[data_heat$paralog == 1] <- NA
  
  return(data_heat)
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


# create heatmap plot ----------------------------------------------------------
heatmap_plotting <- function(data,
                             x_axis,
                             var1_id, var2_id,
                             low_color_var1, high_color_var1,
                             low_color_var2, high_color_var2,
                             para_color,
                             x_size, y_size, legend_size,
                             main_legend,
                             dot_zoom,
                             x_angle,
                             guideline){
  data_heat <- data
  
  # rescale numbers of paralogs
  data_heat$paralog <- as.numeric(data_heat$paralog)
  if (length(unique(na.omit(data_heat$paralog))) > 0) {
    max_paralog <- max(na.omit(data_heat$paralog))
    data_heat$paralogSize <- (data_heat$paralog / max_paralog) * 3
  }
  
  # remove prefix number of taxa names but keep the order
  data_heat$supertaxon <- {
    mapvalues(
      warn_missing = F,
      data_heat$supertaxon,
      from = as.character(data_heat$supertaxon),
      to = substr(as.character(data_heat$supertaxon),
                  6,
                  nchar(as.character(data_heat$supertaxon)))
    )
  }
  
  # format plot
  if (x_axis == "genes") {
    p <- ggplot(data_heat, aes(x = geneID, y = supertaxon))
  } else{
    p <- ggplot(data_heat, aes(y = geneID, x = supertaxon))
  }
  if (length(unique(na.omit(data_heat$var2))) != 1) {
    p <- p + scale_fill_gradient(low = low_color_var2,
                                 high = high_color_var2,
                                 na.value = "gray95",
                                 limits = c(0, 1)) +  # fill color (var2)
      geom_tile(aes(fill = var2))    # filled rect (var2 score)
  }
  if (length(unique(na.omit(data_heat$presSpec))) < 3) {
    if (length(unique(na.omit(data_heat$var1))) == 1) {
      # geom_point for circle illusion (var1 and presence/absence)
      p <- p + geom_point(aes(colour = var1),
                          size = data_heat$presSpec * 5 * (1 + dot_zoom),
                          na.rm = TRUE, show.legend = F)
    } else {
      # geom_point for circle illusion (var1 and presence/absence)
      p <- p + geom_point(aes(colour = var1),
                          size = data_heat$presSpec * 5 * (1 + dot_zoom),
                          na.rm = TRUE)
      # color of the corresponding aes (var1)
      p <- p + scale_color_gradient(low = low_color_var1,
                                    high = high_color_var1,
                                    limits = c(0, 1))
    }
  } else {
    if (length(unique(na.omit(data_heat$var1))) == 1) {
      # geom_point for circle illusion (var1 and presence/absence)
      p <- p + geom_point(aes(size = presSpec),
                          color = "#336a98",
                          na.rm = TRUE)
    } else {
      # geom_point for circle illusion (var1 and presence/absence)
      p <- p + geom_point(aes(colour = var1, size = presSpec),
                          na.rm = TRUE)
      # color of the corresponding aes (var1)
      p <- p + 
        scale_color_gradient(
          low = low_color_var1, high = high_color_var1,
          limits = c(0, 1)
        )
    }
  }
  
  # plot inparalogs (if available)
  if (length(unique(na.omit(data_heat$paralog))) > 0) {
    p <- p + geom_point(data = data_heat,
                        aes(size = paralog),
                        color = para_color,
                        na.rm = TRUE,
                        show.legend = TRUE)
    p <- p + guides(size = guide_legend(title = "# of co-orthologs"))
    
    # to tune the size of circles
    p <- p + 
      scale_size_continuous(
        range = c(min(na.omit(data_heat$paralogSize)) * (1 + dot_zoom),
                  max(na.omit(data_heat$paralogSize)) * (1 + dot_zoom))
      )
  } else {
    # remain the scale of point while filtering
    present_vl <- data_heat$presSpec[!is.na(data_heat$presSpec)]
    
    # to tune the size of circles;
    # use "floor(value*10)/10" to round "down" the value with one decimal number
    p <- p +
      scale_size_continuous(
        range = c((floor(min(present_vl) * 10) / 10 * 5) * (1 + dot_zoom),
                  (floor(max(present_vl) * 10) / 10 * 5) * (1 + dot_zoom))
      )
  }
  p <- p + guides(fill = guide_colourbar(title = var2_id),
                  color = guide_colourbar(title = var1_id))
  base_size <- 9
  
  # guideline for separating ref species
  if (guideline == 1) {
    if (x_axis == "genes") {
      p <- p + labs(y = "Taxon")
      p <- p + geom_hline(yintercept = 0.5, colour = "dodgerblue4")
      p <- p + geom_hline(yintercept = 1.5, colour = "dodgerblue4")
    } else{
      p <- p + labs(x = "Taxon")
      p <- p + geom_vline(xintercept = 0.5, colour = "dodgerblue4")
      p <- p + geom_vline(xintercept = 1.5, colour = "dodgerblue4")
    }
  }
  
  # format theme
  p <- p + theme_minimal()
  p <- p + theme(axis.text.x = element_text(angle = x_angle,
                                            hjust = 1,
                                            size = x_size),
                 axis.text.y = element_text(size = y_size),
                 axis.title.x = element_text(size = x_size),
                 axis.title.y = element_text(size = y_size),
                 legend.title = element_text(size = legend_size),
                 legend.text = element_text(size = legend_size),
                 legend.position = main_legend)
  
  # return plot
  return(p)
}

# highlight gene and/or taxon of interest --------------------------------------
profile_plot <- function(data_heat,
                         plot_parameter,
                         taxon_name,
                         rank_select,
                         gene_highlight){
  # get heatmap
  p <- heatmap_plotting(data_heat,
                        plot_parameter$x_axis,
                        plot_parameter$var1_id,
                        plot_parameter$var2_id,
                        plot_parameter$low_color_var1,
                        plot_parameter$high_color_var1,
                        plot_parameter$low_color_var2,
                        plot_parameter$high_color_var2,
                        plot_parameter$para_color,
                        plot_parameter$x_size,
                        plot_parameter$y_size,
                        plot_parameter$legend_size,
                        plot_parameter$main_legend,
                        plot_parameter$dot_zoom,
                        plot_parameter$x_angle,
                        plot_parameter$guideline)
  
  # highlight taxon
  if (taxon_name != "none") {
    # get selected highlight taxon ID
    rank_select <- rank_select
    # get rank name from rank_select
    rank_rame <- substr(rank_select,
                        4,
                        nchar(rank_select))
    taxa_list <- as.data.frame(read.table("data/taxonNamesReduced.txt",
                                          sep = "\t",
                                          header = T))
    taxon_highlight_id <- {
      taxa_list$ncbiID[taxa_list$fullName == taxon_name
                       & taxa_list$rank == rank_rame]
    }
    
    if (length(taxon_highlight_id) == 0L) {
      taxon_highlight_id <- {
        taxa_list$ncbiID[taxa_list$fullName == taxon_name]
      }
    }
    
    # get taxonID together with it sorted index
    highlight_taxon <- {
      toString(data_heat[data_heat$supertaxonID == taxon_highlight_id, 2][1])
    }
    
    # get index
    selected_index <- as.numeric(as.character(substr(highlight_taxon, 2, 4)))
    
    # draw a rect to highlight this taxon's column
    if (plot_parameter$x_axis == "taxa") {
      rect <- data.frame(xmin = selected_index - 0.5,
                         xmax = selected_index + 0.5,
                         ymin = -Inf,
                         ymax = Inf)
    } else {
      rect <- data.frame(ymin = selected_index - 0.5,
                         ymax = selected_index + 0.5,
                         xmin = -Inf,
                         xmax = Inf)
    }
    
    p <- p + geom_rect(data = rect,
                       aes(xmin = xmin, xmax = xmax,
                           ymin = ymin, ymax = ymax),
                       color = "yellow",
                       alpha = 0.3,
                       inherit.aes = FALSE)
  }
  
  # highlight gene
  if (gene_highlight != "none") {
    # get selected highlight gene ID
    gene_highlight <- gene_highlight
    
    # get index
    all_genes <- levels(data_heat$geneID)
    selected_index <- match(gene_highlight, all_genes)
    
    # draw a rect to highlight this taxon's column
    if (plot_parameter$x_axis == "taxa") {
      rect <- data.frame(ymin = selected_index - 0.5,
                         ymax = selected_index + 0.5,
                         xmin = -Inf,
                         xmax = Inf)
    } else {
      rect <- data.frame(xmin = selected_index - 0.5,
                         xmax = selected_index + 0.5,
                         ymin = -Inf,
                         ymax = Inf)
    }
    
    p <- p + geom_rect(data = rect,
                       aes(xmin = xmin, xmax = xmax,
                           ymin = ymin, ymax = ymax),
                       color = "yellow",
                       alpha = 0.3,
                       inherit.aes = FALSE)
  }
  
  return(p)
}