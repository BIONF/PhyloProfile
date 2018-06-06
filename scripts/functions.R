# =============================================================================
# FUNCTIONS ===================================================================
# =============================================================================

# unsorting function to keep user defined geneID order ------------------------
unsort_id <- function(data, order){
  data$geneID <- as.factor(data$geneID)
  if (order == FALSE){
    # keep user defined geneID order
    data$geneID <- factor(data$geneID, levels = unique(data$geneID))
  }
  return(data)
}



# HEATMAP (main and selected plot) ============================================

# plot profile heatmap --------------------------------------------------------
heatmap_plotting <- function(data,
                             x_axis,
                             var1_id,
                             var2_id,
                             low_color_var1,
                             high_color_var1,
                             low_color_var2,
                             high_color_var2,
                             para_color,
                             x_size,
                             y_size,
                             legend_size,
                             main_legend,
                             dot_zoom,
                             x_angle,
                             guideline){
  data_heat <- data

  # rescale numbers of paralogs
  data_heat$paralog <- as.numeric(data_heat$paralog)
  if (length(unique(na.omit(data_heat$paralog))) > 0){
    max_paralog <- max(na.omit(data_heat$paralog))
    data_heat$paralogSize <- (data_heat$paralog / max_paralog) * 3
  }

  # remove prefix number of taxa names but keep the order
  data_heat$supertaxon <- {
    mapvalues(warn_missing = F,
              data_heat$supertaxon,
              from = as.character(data_heat$supertaxon),
              to = substr(as.character(data_heat$supertaxon),
                          6,
                          nchar(as.character(data_heat$supertaxon))))
  }

  # format plot
  if (x_axis == "genes"){
    p <- ggplot(data_heat, aes(x = geneID, y = supertaxon)) # global aes
  } else{
    p <- ggplot(data_heat, aes(y = geneID, x = supertaxon)) # global aes
  }
  if (length(unique(na.omit(data_heat$var2))) != 1){
    p <- p + scale_fill_gradient(low = low_color_var2,
                                high = high_color_var2,
                                na.value = "gray95",
                                limits = c(0, 1)) +  # fill color (var2)
      geom_tile(aes(fill = var2))    # filled rect (var2 score)
  }
  if (length(unique(na.omit(data_heat$presSpec))) < 3){
    if (length(unique(na.omit(data_heat$var1))) == 1){
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
    if (length(unique(na.omit(data_heat$var1))) == 1){
      # geom_point for circle illusion (var1 and presence/absence)
      p <- p + geom_point(aes(size = presSpec),
                         color = "#336a98",
                         na.rm = TRUE)
    } else {
      # geom_point for circle illusion (var1 and presence/absence)
      p <- p + geom_point(aes(colour = var1, size = presSpec),
                         na.rm = TRUE)
      # color of the corresponding aes (var1)
      p <- p + scale_color_gradient(low = low_color_var1, high = high_color_var1,
                                   limits = c(0, 1))
    }
  }

  # plot inparalogs (if available)
  if (length(unique(na.omit(data_heat$paralog))) > 0){
    p <- p + geom_point(data = data_heat,
                        aes(size = paralog),
                        color = para_color,
                        na.rm = TRUE,
                        show.legend = TRUE)
    p <- p + guides(size = guide_legend(title = "# of co-orthologs"))

    # to tune the size of circles;
    # "floor(value*10)/10" is used to round "down" the value with one decimal number
    p <- p + scale_size_continuous(range = c(min(na.omit(data_heat$paralogSize)) * (1 + dot_zoom),
                                             max(na.omit(data_heat$paralogSize)) * (1 + dot_zoom)))
  } else {
    # remain the scale of point while filtering
    present_vl <- data_heat$presSpec[!is.na(data_heat$presSpec)]

    # to tune the size of circles;
    # "floor(value*10)/10" is used to round "down" the value with one decimal number
    p <- p + scale_size_continuous(range = c( (floor(min(present_vl) * 10) / 10 * 5) * (1 + dot_zoom),
                                             (floor(max(present_vl) * 10) / 10 * 5) * (1 + dot_zoom)))
  }
  p <- p + guides(fill = guide_colourbar(title = var2_id),
                 color = guide_colourbar(title = var1_id))
  base_size <- 9

  # guideline for separating ref species
  if (guideline == 1){
    if (x_axis == "genes"){
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

# create profile heatmap ------------------------------------------------------
# v, dataHeat(), clusteredDataHeat(), get_input_main ()
main_plot <- function(v, data_heat, clustered_data_heat, input){
  if (v$doPlot == FALSE) return()

  # cluster dataHeat (if selected)
  if (input$apply_cluster == TRUE){
    data_heat <- clustered_data_heat #()
  }

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

  p <- heatmap_plotting(data_heat,
                        input$x_axis,
                        input$var1_id,
                        input$var2_id,
                        input$low_color_var1,
                        input$high_color_var1,
                        input$low_color_var2,
                        input$high_color_var2,
                        input$para_color,
                        input$x_size,
                        input$y_size,
                        input$legend_size,
                        input$main_legend,
                        input$dot_zoom,
                        input$x_angle,
                        1)

  # highlight taxon
  if (input$taxon_highlight != "none"){
    # get selected highlight taxon ID
    rank_select <- input$rank_select
    # get rank name from rank_select
    rank_rame <- substr(rank_select,
                      4,
                      nchar(rank_select))
    taxa_list <- as.data.frame(read.table("data/taxonNamesReduced.txt",
                                          sep = "\t",
                                          header = T))
    taxon_highlight <- {
      taxa_list$ncbiID[taxa_list$fullName == input$taxon_highlight
                       & taxa_list$rank == rank_rame]
    }
    if (length(taxon_highlight) == 0L){
      taxon_highlight <- {
        taxa_list$ncbiID[taxa_list$fullName == input$taxon_highlight]
      }
    }

    # get taxonID together with it sorted index
    highlight_taxon <- {
      toString(data_heat[data_heat$supertaxonID == taxon_highlight, 2][1])
    }

    # get index
    selected_index <- as.numeric(as.character(substr(highlight_taxon, 2, 4)))

    # draw a rect to highlight this taxon's column
    if (input$x_axis == "taxa"){
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
  if (input$gene_highlight != "none"){
    # get selected highlight gene ID
    gene_highlight <- input$gene_highlight

    # get index
    all_genes <- levels(data_heat$geneID)
    selected_index <- match(gene_highlight, all_genes)

    # draw a rect to highlight this taxon's column
    if (input$x_axis == "taxa"){
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

  # do plotting
  if (input$auto_update == FALSE){
    # Add dependency on the update button
    # (only update when button is clicked)
    input$update_btn

    # Add all the filters to the data based on the user inputs
    # wrap in an isolate() so that the data won't update every time an input
    # is changed
    isolate({
      p
    })
  } else {
    p
  }
}

# create plot (same as main plot) -------------------------------------------
# vCt, dataHeat(), clusteredDataHeat(), get_input_selected()
selected_plot <- function(v_ct, data_heat, clustered_data_heat, input){
  if (v_ct$doPlotCustom == FALSE) return()
  if (input$in_seq[1] == "all" & input$in_taxa[1] == "all") return()
  else{
    # cluster dataHeat (if selected)
    if (input$apply_cluster == TRUE){
      data_heat <- clustered_data_heat # ()
    }

    # process data
    data_heat$supertaxonMod <- {
      substr(data_heat$supertaxon,
             6,
             nchar(as.character(data_heat$supertaxon)))
    }

    if (input$in_taxa[1] == "all" & input$in_seq[1] != "all"){
      # select data from dataHeat for selected sequences only
      data_heat <- subset(data_heat, geneID %in% input$in_seq)
    } else if (input$in_seq[1] == "all" & input$in_taxa[1] != "all"){
      # select data from dataHeat for selected taxa only
      data_heat <- subset(data_heat, supertaxonMod %in% input$in_taxa)
    } else {
      # select data from dataHeat for selected sequences and taxa
      data_heat <- subset(data_heat,
                         geneID %in% input$in_seq
                         & supertaxonMod %in% input$in_taxa)
    }

    # remove unneeded dots
    data_heat$presSpec[data_heat$presSpec == 0] <- NA
    data_heat$paralog[data_heat$presSpec < 1] <- NA
    data_heat$paralog[data_heat$paralog == 1] <- NA

    # create plot
    p <- heatmap_plotting(data_heat,
                          input$x_axis_selected,
                          input$var1_id,
                          input$var2_id,
                          input$low_color_var1,
                          input$high_color_var1,
                          input$low_color_var2,
                          input$high_color_var2,
                          input$para_color,
                          input$x_size_select,
                          input$y_size_select,
                          input$legend_size_select,
                          input$selected_legend,
                          input$dot_zoom_select,
                          input$x_angle_select, 0)

    ### do plotting
    if (input$auto_update_selected == FALSE){
      # Add dependency on the update button (only update when button is clicked)
      input$plot_custom

      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        p
      })
    } else {
      p
    }
  }
}

# render detailed plot ------------------------------------------------------
# v, detai_plotDt(), input$detailed_text
detail_plot <- function(v, sel_df, detailed_text, var1_id, var2_id){
  if (v$doPlot == FALSE) return()

  sel_df$x_label <- paste(sel_df$orthoID,
                         " (",
                         sel_df$fullName,
                         ")",
                         sep = "")

  # if(input$detailed_remove_na == TRUE){
  #   sel_df <- sel_df[!is.na(sel_df$orthoID),]
  # }

  # create joined DF for plotting var1 next to var2
  var1Df <- subset(sel_df, select = c("x_label", "var1"))
  var1Df$type <- var1_id
  colnames(var1Df) <- c("id", "score", "var")

  var2Df <- subset(sel_df, select = c("x_label", "var2"))
  var2Df$type <- var2_id
  colnames(var2Df) <- c("id", "score", "var")

  detailed_df <- rbind(var1Df, var2Df)

  # remove ONE missing variable
  if (nlevels(as.factor(detailed_df$var)) > 1){
    detailed_df <- detailed_df[nchar(detailed_df$var) > 0, ]
  }

  # keep order of ID (x_label)
  detailed_df$id <- factor(detailed_df$id, levels = unique(detailed_df$id))

  # create plot
  gp <- ggplot(detailed_df, aes(y = score, x = id, fill = var)) +
    geom_bar(stat = "identity", position = position_dodge(), na.rm = TRUE) +
    coord_flip() +
    labs(x = "") +
    labs(fill = "") +
    theme_minimal()
  #geom_text(aes(label=var1), vjust=3)
  gp <- gp + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                axis.text = element_text(size = detailed_text),
                axis.title = element_text(size = detailed_text),
                legend.text = element_text(size = detailed_text)
  )
  gp
}

# ARCHITECTURE PLOT ===========================================================

# create domain plot ----------------------------------------------------------
# v3, point_infoDetail(), getDomainFile(), input$demo_data,
# input$one_seq_fasta, input$label_archi_size, input$title_archi_size
archi_plot <- function(v3, 
                       info, domain_df, 
                       one_seq_fasta, label_archi_size, title_archi_size){
  if (v3$doPlot3 == FALSE) return()
  # info
  group <- as.character(info[1])
  ortho <- as.character(info[2])
  var1 <- as.character(info[3])

  # get sub dataframe based on selected group_id and orthoID
  ortho <- gsub("\\|", ":", ortho)
  grepID <- paste(group, "#", ortho, sep = "")
  
  subdomain_df <- domain_df[grep(grepID, domain_df$seedID), ]
  subdomain_df$feature <- as.character(subdomain_df$feature)

  if (nrow(subdomain_df) < 1){
    v3$doPlot3 <- FALSE
    return()
  } else {

    # ortho domains df
    ortho_df <- filter(subdomain_df, orthoID == ortho)

    # seed domains df
    seed_df <- filter(subdomain_df, orthoID != ortho)

    if (nrow(seed_df) == 0) seed_df <- ortho_df

    seed <- as.character(seed_df$orthoID[1])

    # change order of one dataframe's features
    # based on order of other df's features
    if (length(ortho_df$feature) < length(seed_df$feature)){
      ordered_ortho_df <- ortho_df[order(ortho_df$feature), ]
      ordered_seed_df <- sort_domains(ordered_ortho_df, seed_df)
    } else {
      ordered_seed_df <- seed_df[order(seed_df$feature), ]
      ordered_ortho_df <- sort_domains(ordered_seed_df, ortho_df)
    }

    # join weight values and feature names
    if ("weight" %in% colnames(ordered_ortho_df)){
      ordered_ortho_df$yLabel <- paste0(ordered_ortho_df$feature,
                                      " (",
                                      round(ordered_ortho_df$weight, 2),
                                      ")")
      ordered_ortho_df$feature <- ordered_ortho_df$yLabel
    }
    if ("weight" %in% colnames(ordered_seed_df)){
      ordered_seed_df$yLabel <- paste0(ordered_seed_df$feature,
                                     " (",
                                     round(ordered_seed_df$weight, 2),
                                     ")")
      ordered_seed_df$feature <- ordered_seed_df$yLabel
    }

    # plotting
    sep <- ":"
    if (!is.null(one_seq_fasta)) sep <- "|"
    if ("length" %in% colnames(subdomain_df)){
      plot_ortho <- domain_plotting(ordered_ortho_df,
                                    ortho,
                                    sep,
                                    label_archi_size,
                                    title_archi_size,
                                    min(subdomain_df$start),
                                    max(c(subdomain_df$end,
                                          subdomain_df$length)))
      plot_seed <- domain_plotting(ordered_seed_df,
                                   seed,
                                   sep,
                                   label_archi_size,
                                   title_archi_size,
                                   min(subdomain_df$start),
                                   max(c(subdomain_df$end,
                                         subdomain_df$length)))

    } else{
      plot_ortho <- domain_plotting(ordered_ortho_df,
                                    ortho,
                                    sep,
                                    label_archi_size,
                                    title_archi_size,
                                    min(subdomain_df$start),
                                    max(subdomain_df$end))
      plot_seed <- domain_plotting(ordered_seed_df,
                                   seed,
                                   sep,
                                   label_archi_size,
                                   title_archi_size,
                                   min(subdomain_df$start),
                                   max(subdomain_df$end))
    }

    # grid.arrange(plot_seed,plot_ortho,ncol=1)

    if (ortho == seed){
      arrangeGrob(plot_seed, ncol = 1)
    } else {
      seed_height <- length(levels(as.factor(ordered_seed_df$feature)))
      ortho_height <- length(levels(as.factor(ordered_ortho_df$feature)))

      arrangeGrob(plot_seed, plot_ortho, ncol = 1,
                  heights = c(seed_height, ortho_height))
    }
  }
}

# sort one domain dataframe (ortho) based on the other domain Df (seed) -------
sort_domains <- function(seed_df, ortho_df){
  # get list of features in seed_df
  feature_list <- as.data.frame(levels(as.factor(seed_df$feature)))
  colnames(feature_list) <- c("feature")
  # and add order number to each feature
  feature_list$orderNo <- seq(length(feature_list$feature))

  # merge those info to ortho_df
  ordered_ortho_df <- merge(ortho_df, feature_list, all.x = TRUE)

  # sort ortho_df
  index <- with(ordered_ortho_df, order(orderNo))
  ordered_ortho_df <- ordered_ortho_df[index, ]

  #turn feature column into a character vector
  ordered_ortho_df$feature <- as.character(ordered_ortho_df$feature)
  #then turn it back into an ordered factor (to keep this order while plotting)
  ordered_ortho_df$feature <- factor(ordered_ortho_df$feature,
                                   levels = unique(ordered_ortho_df$feature))
  #return sorted df
  ordered_ortho_df
}

# plot domain architecture ----------------------------------------------------
domain_plotting <- function(df,
                            geneID,
                            sep,
                            label_size,
                            title_size,
                            min_start,
                            max_end){
  gg <- ggplot(df, aes(y = feature, x = end, color = feature)) +
    geom_segment(data = df,
                 aes(y = feature, yend = feature,
                     x = min_start, xend = max_end),
                 color = "white",
                 size = 0)

  # draw lines for representing sequence length
  gg <- gg + geom_segment(data = df,
                          aes(x = 0, xend = length,
                              y = feature, yend = feature),
                          size = 1,
                          color = "#b2b2b2")

  # draw line and points
  gg <- gg + geom_segment(data = df,
                          aes(x = start, xend = end,
                              y = feature, yend = feature),
                          size = 1.5)
  gg <- gg + geom_point(data = df,
                        aes(y = feature, x = start),
                        color = "#b2b2b2",
                        size = 3,
                        shape = 3)
  gg <- gg + geom_point(data = df,
                        aes(y = feature, x = end),
                        color = "#edae52",
                        size = 3,
                        shape = 5)

  # draw dashed line for domain path
  gg <- gg + geom_segment(data = df[df$path == "Y", ],
                          aes(x = start, xend = end,
                              y = feature, yend = feature),
                          size = 3,
                          linetype = "dashed")

  # # add text above
  # gg <- gg + geom_text(data = df,
  #                      aes(x = (start + end) / 2,
  #                          y = feature, label = round(weight,2)),
  #                        color = "#9fb059",
  #                        size = descSize,
  #                        vjust = -0.75,
  #                        fontface = "bold",
  #                        family = "serif")

  # theme format
  title_mod <- gsub(":", sep, geneID)
  gg <- gg + scale_y_discrete(expand = c(0.075, 0))
  gg <- gg + labs(title = paste0(title_mod), y = "Feature")
  gg <- gg + theme_minimal()
  gg <- gg + theme(panel.border = element_blank())
  gg <- gg + theme(axis.ticks = element_blank())
  gg <- gg + theme(plot.title = element_text(face = "bold", size = title_size))
  gg <- gg + theme(plot.title = element_text(hjust = 0.5))
  gg <- gg + theme(legend.position = "none", axis.title.x = element_blank(),
                   axis.text.y = element_text(size = label_size),
                   axis.title.y = element_text(size = label_size),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.x = element_blank())
  # return plot
  return(gg)
}

# gene age plot ---------------------------------------------------------------
# gene_ageDfMod(), input$gene_age_text
gene_age_plot <- function(count_df, gene_age_text){
  p <- ggplot(count_df, aes(fill = age, y = percentage, x = 1)) +
    geom_bar(stat = "identity") +
    scale_y_reverse() +
    coord_flip() +
    theme_minimal()
  p <- p + geom_text(data = count_df,
                     aes(x = 1, y = 100 - pos,
                         label = paste0(freq, "\n", percentage, "%")),
                     size = 4 * gene_age_text)
  p <- p + theme(legend.position = "bottom",
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12 * gene_age_text),
                 axis.title = element_blank(), axis.text = element_blank()) +
    scale_fill_brewer(palette = "Spectral") +
    guides(fill = guide_legend(nrow = round(nrow(count_df) / 3, 0),
                               byrow = TRUE))
  p
}

# GROUP COMPARISON ============================================================

# print list of available taxa ------------------------------------------------
taxa_select_gc <- function(rank_select_gc, subset_taxa){

  # if there is no rank set, there can not be any available taxa
  if (length(rank_select_gc) == 0) return()
  else{

    # load list of unsorted taxa
    dt <- get_taxa_list(TRUE, subset_taxa)

    # load list of taxon name
    name_list <- get_name_list(TRUE, FALSE)

    # get rank name from rank_select
    rank_name <- substr(rank_select_gc, 4, nchar(rank_select_gc))
    choice <- as.data.frame
    choice <- rbind(dt[rank_name])
    colnames(choice) <- "ncbiID"
    choice <- merge(choice, name_list, by = "ncbiID", all = FALSE)
    return(choice)
  }
}

# Generate the list with all plots --------------------------------------------
get_plot_output_list <- function(genes, input, interesting_features) {
  # if we dont have parameters we can not generate plots
  if (is.null(genes)) return()
  # Insert plot output objects the list
    plot_output_list <- lapply(1:nrow(genes), function(i) {
      plotname <- paste(genes[i, 1])
      plot_output_object <- renderPlot(get_multiplot(genes[i, ], input,
                                                     interesting_features),
                                       height = 650, width = 700)
    })
    do.call(tagList, plot_output_list) # needed to display properly.
    return(plot_output_list)

}

# Put the plots for one spicific gene in one multiplot ------------------------
# get_input_gc
get_multiplot <- function(gene_info, input, interesting_features){

  # Sorting the information to the selected gene
  gene <- as.character(gene_info$geneID)
  in_group <- as.data.frame(gene_info$in_group)
  out_group <- as.data.frame(gene_info$out_group)
  features <- as.data.frame(gene_info$features)
  pvalues <- gene_info$pvalues
  var <- gene_info$var

  # the barplot does not depent on the selected variable(s)
  barplot <-  get_barplot_gc(gene,
                             in_group,
                             out_group,
                             features,
                             input,
                             interesting_features)

  if (is.null(barplot)){
    barplot <- textGrob("The selected domains are not found in the gene")
  }
  # if both variables are selected there are going to be 2 boxplots
  if (var == "Both"){
    pvalues <- unlist(pvalues, recursive = FALSE)

    p1 <- unlist(pvalues[1])
    p2 <- unlist(pvalues[2])
    # Check if the p_values should be printed
    if (input$show_p_value == TRUE){
      info_p1 <- get_info_p_values(p1)
      info_p2 <- get_info_p_values(p2)
    }
    else{
      info_p1 <- " "
      info_p2 <- " "
    }
    # check if the significant plot should be highlighted
    if (input$highlight_significant == TRUE){
      if (is.null (p1[1])) c1 <- "grey"
      else if (p1[length(p1)] < input$significance) c1 <- "indianred2"
      else c1 <- "grey"

      if (is.null (p2[1])) c2 <- "grey"
      else if (p2[length(p2)] < input$significance) c2 <- "indianred2"
      else c2 <- "grey"
    }
    else{
      c1 <- "grey"
      c2 <- "grey"
    }

    boxplot1 <- get_boxplot_gc(in_group,
                               out_group,
                               input$var1_id,
                               gene,
                               c1,
                               info_p1, input)

    boxplot2 <- get_boxplot_gc(in_group,
                               out_group,
                               input$var2_id,
                               gene,
                               c2,
                               info_p2, input)

    m <- grid.arrange(textGrob(gene),
                      arrangeGrob(boxplot1, boxplot2, ncol = 2),
                      barplot,
                      heights = c(0.02, 0.45, 0.458), ncol = 1)
  }else {
    p <- unlist(pvalues)
    if (input$show_p_value == TRUE){
      info <- get_info_p_values(p)
    }else{
      info <- " "
    }

    boxplot <- get_boxplot_gc(in_group,
                              out_group,
                              var,
                              gene,
                              "grey",
                              info,
                              input)

    m <- grid.arrange(textGrob(gene),
                      boxplot,
                      barplot,
                      heights = c(0.02, 0.45, 0.458), ncol = 1)
  }
  return(m)
}

# Create a Boxplot ------------------------------------------------------------
#
get_boxplot_gc <- function (in_group_df,
                            out_group_df,
                            var,
                            gene,
                            colour,
                            info,
                            input){

  if (var == input$var1_id){
    in_g <- in_group_df$var1
    out_g <- out_group_df$var1
  }
  else if (var == input$var2_id) {
    in_g <- in_group_df$var2
    out_g <- out_group_df$var2
  }

  a <- length(in_g)
  b <- length(out_g)

  in_g <- as.data.frame(in_g)
  names(in_g)[1] <- paste("values")
  in_g$group <- "in_group"

  out_g <- as.data.frame(out_g)
  names(out_g)[1] <- paste("values")
  out_g$group <- "Out-Group"

  data_boxplot <- rbind(in_g, out_g)
  data_boxplot <- data_boxplot[complete.cases(data_boxplot), ]

  names <- c(paste("In-Group \n n=", a, sep = ""),
             paste("Out-Group \n n=", b, sep = ""))

  boxplot_gc <- ggplot(data_boxplot, aes(group, values)) +
    geom_boxplot (stat = "boxplot",
                  position = position_dodge(),
                  width = 0.5,
                  fill = colour) +
    labs(x = "", y = var, caption = paste(info)) +
    scale_x_discrete(labels = names) +
    theme_minimal()

  boxplot_gc <- boxplot_gc +
    theme(axis.text.x = element_text(size = input$x_size_gc, hjust = 1),
          axis.text.y = element_text(size = input$y_size_gc),
          axis.title.y = element_text(size = input$y_size_gc))
  boxplot_gc
}

# Create Barplot ------------------------------------------------------------
get_barplot_gc <- function(selected_gene,
                           in_group, out_group,
                           features, input, interesting_features){
  subdomain_df <- features
  subdomain_df$feature <- as.character(subdomain_df$feature)

  # only show features that interest the user
  if (!("all" %in% input$interesting_features)){
    ifeatures <- NULL
    for (x in input$interesting_features){
      f <- subset(subdomain_df$feature,
                  startsWith(subdomain_df$feature, x))
      ifeatures <- append(ifeatures, f)
    }

    if (is.null(ifeatures))return()
    # only keep rows in which the feature begins with a element out of the
    # interesing Features
    subdomain_df <- subset(subdomain_df, subdomain_df$feature %in% ifeatures)
  }

  # part in in_group and Out-Group
  in_group_domain_df  <-  {
    subset(subdomain_df, subdomain_df$orthoID %in% in_group$orthoID)
  }
  out_group_domain_df <- {
    subset(subdomain_df, subdomain_df$orthoID %in% out_group$orthoID)
  }

  # Get list of all seeds
  seeds <- unique (subdomain_df$seedID)
  in_not_empty <- 0
  out_not_empty <- 0

  if (nrow(in_group_domain_df) == 0) {
    data_in <- NULL
    } else{
    feature <- unique(in_group_domain_df$feature)
    data_in <- as.data.frame(feature)
    data_in$amount <- 0
  }

  if (nrow(out_group_domain_df) == 0) {
    data_out <- NULL
    } else{
    feature <- unique(out_group_domain_df$feature)
    data_out <- as.data.frame(feature)
    data_out$amount <- 0
  }

  for (seed in seeds){
    if (!is.null(data_in)){
      in_g <- subset(in_group_domain_df, in_group_domain_df$seedID == seed)
      if (!empty(in_g)){
        in_not_empty <- in_not_empty + 1
        in_group_features <-  plyr::count(in_g, "feature")
        for (i in 1:nrow(in_group_features)){
          for (j in 1:nrow(data_in)){
            if (data_in[j, 1] == in_group_features[i, 1]){
              data_in[j, 2] <- data_in[j, 2] + in_group_features[i, 2]
            }
          }
        }
      }
    }

    if (!is.null(data_out)){
      out_g <- {
        subset(out_group_domain_df, out_group_domain_df$seedID == seed)
      }

      if (!empty(out_g)){
        out_not_empty <- out_not_empty + 1
        out_group_features <-  plyr::count(out_g, "feature")
        for (i in 1:nrow(out_group_features)){
          for (j in 1:nrow(data_out)){
            if (data_out[j, 1] == out_group_features[i, 1]){
              data_out[j, 2] <- data_out[j, 2] + out_group_features[i, 2]
            }
          }
        }
      }
    }
  }

  if (!is.null(data_in)){
    data_in$amount <- data_in$amount / in_not_empty
    data_in$type <- "in_group"
  }

  if (!is.null(data_out)){
    data_out$amount <- data_out$amount / out_not_empty
    data_out$type <- "Out-Group"
  }

  if (is.null(data_in) & !is.null(data_out)) {
    data_barplot <- data_out
    }  else if (is.null(data_out) & !is.null(data_in)){
      data_barplot <- data_in
      } else if (!is.null(data_in) & !is.null(data_out)){
    data_barplot <- rbind(data_in, data_out)
  } else{
    data_barplot <- NULL
  }

  if (!is.null(data_barplot)){
    # generate Barplot
    barplot_gc <- ggplot(data_barplot,
                         aes(x = feature, y = amount, fill = type ),
                         main  = " ") +
      geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
      scale_fill_grey() +
      labs(x = " ", y = "Average instances per protein", fill = "Group") +
      theme_minimal()

    barplot_gc <- barplot_gc +
      theme(axis.text.x = element_text(size = input$x_size_gc,
                                       angle = input$angle_gc, hjust = 1),
            axis.text.y = element_text(size = input$y_size_gc),
            axis.title.y = element_text(size = input$y_size_gc),
            legend.position = input$legend_gc,
            legend.text = element_text(size = input$legend_size_gc ),
            legend.title = element_text(size = input$legend_size_gc))
    barplot_gc
  } else (return(NULL))
}

# Get the p_values to print under the plot ----------------------------------
get_info_p_values <- function (p){

  if (is.na(p[1]))info_p_values <- "not enough information"
  else if (length(p) == 1){
    info_p_values <- paste("Kolmogorov-Smirnov-Test:",
                           as.character(p[1]), sep = " " )
  }
  else{
    info_p_values1 <- paste("Kolmogorov-Smirnov-Test:",
                            as.character(p[1]), sep = " " )
    info_p_values2 <- paste("Wilcoxon-Mann-Whitney-Test: ",
                            as.character(round(p[2], 6)), sep = " ")
    info_p_values <- paste(info_p_values1, info_p_values2, sep = "\n")
  }

  info_p_values <- paste("p_values:", info_p_values, sep = "\n")
  info_p_values
}

# Get the list with all the significant genes and the dataset ---------------
# input$selected_in_group_gc, input$list_selected_genes_gc, 
# input$rank_select, input$var_name_gc, input$use_common_anchestor,
# input$inSelect
get_significant_genes <- function (in_group,
                                   list_selected_genes_gc,
                                   rank,
                                   var,
                                   use_common_anchestor,
                                   in_select,
                                   input,
                                   demo_data, anno_location, file_domain,
                                   subset_taxa,
                                   data_full,
                                   session,
                                   right_format_features,
                                   domains){
  if (is.null(in_group)
      | length(list_selected_genes_gc) == 0)return()

  # load name List
  name_list <- get_name_list(TRUE, TRUE)

  # load list of unsorted taxa
  dt <- get_taxa_list(FALSE, subset_taxa)

  # Updateing of the Input ==================================================
  # if there is more than one element in the in_group
  # we look at the next common anchstor

  if (use_common_anchestor == TRUE){
    ancestor <- get_common_ancestor(in_group, rank, name_list, dt)
    if (is.null(ancestor))return()
    in_group <- ancestor[1]
    rank <- ancestor[2]
  } else{
    rank <- substring(rank, 4)
  }

  # List with all significant genes
  significantGenes <- data.frame(
    geneID = character(),
    in_group = I(list()),
    out_group = I(list()),
    pvalues = I(list()),
    features = I(list()),
    databases = I(list()))

  # List of genes to look at
  if (is.element("all", list_selected_genes_gc)){
    genes <- data_full$geneID
    genes <- genes[!duplicated(genes)]
  } else {
    genes <- list_selected_genes_gc
  }
  genes <- sort(genes)

  # Subset depending on the rank and the in_group
  selected_subset <- get_selected_subset(rank, in_group, name_list, dt)
  selected_subset <- subset(selected_subset,
                            !selected_subset$fullName == in_select)

  # Check for each of the selected genes if they are significant
  # If a gene is significant generate the plots
  # and save them in significantGenesGroupComparison
  for (gene in genes){

    # creates Substet only for current Gene
    selected_gene_df <- subset(data_full, data_full$geneID == gene)
    # In- and Out-Group depending on the Gene
    in_group_df <- {
      subset(selected_gene_df,
             selected_gene_df$abbrName %in% selected_subset$abbrName)
    }
    out_group_df <- {
      subset (selected_gene_df,
              !(selected_gene_df$abbrName %in% selected_subset$abbrName))
    }
    out_group_df <- {
      subset(out_group_df, !out_group_df$fullName == in_select)
    }

    # Generate the p_values for the gene
    pvalues <- get_significant(in_group_df, out_group_df, var, gene, input)

    if (!is.null(pvalues)){
      new_row <- data.frame(geneID = gene,
                            in_group = NA,
                            out_group = NA,
                            pvalues = NA,
                            features = NA)
      new_row$in_group <- list(in_group_df)
      new_row$out_group <- list(out_group_df)
      new_row$pvalues <- list(pvalues)
      features  <- get_features(gene,
                                domains,
                                session)
      new_row$features <- list(features)
      if (right_format_features){
        new_row$databases <- list(get_prefix_features(features))
      }
      significantGenes <- rbind(significantGenes, new_row)
    }
  }
  if (nrow(significantGenes) != 0){
    significantGenes$var <- var
    significantGenes$rank <- rank
    return(significantGenes)
  }
}

# Get the database for each feature in a specific gene ----------------------
# f is dataframe in $features
get_prefix_features <- function(f){
  features <- f$feature
  choices <- gsub("_.*", "", features)
  choices <- choices[!duplicated(choices)]
  choices
}

# Get the Subset depending on the choosen rank ------------------------------
get_selected_subset <- function (rank, in_group, name_list, dt){
  # Look if the fullName is in the in_group
  name_list$fullName <- as.character(name_list$fullName)

  #rank <- substring(rank, 4)

  name_list_rank <- subset(name_list, name_list$rank == rank)
  x <- subset(name_list, name_list$fullName %in% in_group)
  # Look if it has the right rank
  x1 <- dt[dt[, rank] %in% x$ncbiID, ]
  x1
}

# Generate the in_group -----------------------------------------------------
get_common_ancestor <- function(in_group,
                                rank,
                                name_list,
                                all_selected_in_group_gc){
  # Get the list with all taxonomy ranks
  all_ranks <- get_taxonomy_ranks()
  all_selected_in_group_gc <- {
    all_selected_in_group_gc[!duplicated(all_selected_in_group_gc), ]
  }

  # all ranks higher than the selected rank
  # ranks were all elements of the in_group might be in the same taxon
  possible_ranks <- all_ranks [all_ranks >= rank]
  i <-  1
  if (length(in_group) == 1) rank <- substring(rank, 4)

  # find the common ancestor of all taxa in the in_group
  # and use it as in_group
  while (length(in_group) > 1 & i < length(possible_ranks)){
    current_rank <- substring(possible_ranks[i], 4)
    next_rank <- substring(possible_ranks[i + 1], 4)

    # dataframe with all elements with fitting rank
    in_g <- subset(name_list, name_list$rank == current_rank)

    # subset of in_g with elements that belong to the in_group
    in_g <- subset(in_g, in_g$fullName %in% in_group)

    # dataset with all elements that could belong to the in_group
    possible_in_group <- subset(all_selected_in_group_gc,
                                select = c(current_rank, next_rank))

    possible_in_group <- {
      possible_in_group[possible_in_group[, current_rank] %in% in_g$ncbiID, ]
    }

    possible_in_group <- possible_in_group[!duplicated(possible_in_group), ]
    # only consider elements of the list that have the next higher rank
    x <- subset(name_list, name_list$rank == next_rank)
    x <- subset(x, x$ncbiID %in% possible_in_group[, next_rank] )

    in_group <- x$fullName
    i <- i + 1
    rank <- next_rank
  }
  # If there is no common ancestor
  if (i > length(possible_ranks)) return()
  c(in_group, rank)
}

# Decide if the gene is significant -----------------------------------------
# if the gene is significant return the pvalues
get_significant <- function(in_g, out_g, var, gene, input){
   significance_level <- input$significance

  if (var == "Both"){
    var1 <- input$var1
    var2 <- input$var2

    significant <-  FALSE
    # get the pValues for both variables
    pvalues1 <- get_p_values(in_g$var1, out_g$var1, significance_level)
    pvalues2 <- get_p_values(in_g$var2, out_g$var2, significance_level)

    # check if the gene has a significant difference
    # in the distribioution of In- and Out-Group
    # in at least one of the variables
    # If there is not enough data to calculate a p-value
    # consider the gene as not intereisting
    if (is.null(pvalues1)){

    }
    # check if the at last calculated p-value
    #is smaller than the significane level
    else if (pvalues1[length(pvalues1)] < significance_level){
      significant <- TRUE
    }

    # analog to pvalues in the first variable
    if (is.null(pvalues2)){

    } else if (pvalues2[length(pvalues2)] < significance_level){
      significant <-  TRUE
    }

    # if the gene is interisting return the p_values
    if (significant){
      pvalues <- list(pvalues1, pvalues2)
      return(pvalues)
    }
    # if the gene is not interisting return NULL
    else return(NULL)

  } else{
    # Check which variable is selected and get the p_values
    if (var == input$var1_id){
      pvalues <- get_p_values(in_g$var1, out_g$var1, significance_level)
    }
    else {
      pvalues <- get_p_values(in_g$var2, out_g$var2, significance_level)
    }

    # Analog to getting the significance with both variables
    if (is.null(pvalues)) return(NULL)
    else if (is.nan(pvalues[length(pvalues)])){
      return(NULL)
    }else if (pvalues[length(pvalues)] < significance_level){
      pvalues <-  list(pvalues)
      return(pvalues)
    } else{
      return(NULL)
    }
  }
}

# calculate the p_values ----------------------------------------------------
get_p_values <- function(var_in,
                         var_out,
                         significance_level){

  # delete all entrys that are NA
  var_in <- var_in[!is.na(var_in)]
  var_out <- var_out[!is.na(var_out)]

  # if there is no data in one of the groups the p-value is NA
  if (length(var_in) == 0)return(NULL)
  else if (length(var_out) == 0)return(NULL)

  # if there is data the p-value is calculatet
  #with a combination of two non-parametric test
  else{
    # Kolmogorov-Smirnov Test
    # H0 : The two samples have the same distribution
    ks <- ks.boot(var_in, var_out, alternative = "two.sided")

    # p-value = Probability(test statistic >= value for the samples)
    # probabilitiy to recet H0 if it is correct
    p.value <- ks$ks.boot.pvalue
    # You reject the null hypothesis that the two samples were drawn
    # from the same distribution
    # if the p-value is less than your significance level
    if (p.value < significance_level) pvalue <- c(p.value)

    # if we accept H0 it does not mean that the samples are uninteresting
    # they could still have different location values
    else{
      # Wilcoxon-Mann-Whitney Test
      # H0: the samples have the same location parameters
      wilcox <- wilcox.test(var_in,
                            var_out,
                            alternative = "two.sided",
                            exact = FALSE)
      p_value_wilcox <- wilcox$p.value
      pvalue <- c(p.value, p_value_wilcox)
    }
    pvalue
  }
}

# get the list with all the features in the gene ----------------------------
# input$demo_data
get_features <- function(selected_gene,
                         domains,
                         session){

  subdomain_df <- {
    subset(domain_df,
           substr(domain_df$seedID,
                  1,
                  nchar(as.character(selected_gene))) == selected_gene)
  }
  subdomain_df <- subdomain_df[!duplicated(subdomain_df), ]
  subdomain_df
}

# get the data where to find the features -----------------------------------
# input$demo_data, input$anno_location, input$file_domain_input
# get_domain_file_gc <- function(group,
#                                demo_data,
#                                anno_location,
#                                file_domain,
#                                session,
#                                domain_path){
#   # domain file
#   if (demo_data == "lca-micros" | demo_data == "ampk-tor"){
#     updateButton(session, "do_domain_plot", disabled = FALSE)
#     if (demo_data == "lca-micros"){
#       file_domain <- {
#         suppressWarnings(paste0("https://github.com/BIONF/phyloprofile-data/blob/master/demo/domain_files/",
#                                 group,
#                                 ".domains?raw=true"))
#       }
#     } else {
#       file_domain <- {
#         suppressWarnings(paste0("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.domains_F"))
#       }
#     }
# 
#   }else {
#     if (anno_location == "from file"){
#       # file_domain <- input$fileDomain_input
#       if (is.null(file_domain)){
#         file_domain <- "noFileInput"
#       } else {
#         if (is.null(info)){
#           file_domain <- "noSelectHit"
#           updateButton(session, "do_domain_plot", disabled = TRUE)
#         } else {
#           updateButton(session, "do_domain_plot", disabled = FALSE)
#           file_domain <- file_domain$datapath
#         }
#       }
#     } else {
#       if (is.null(info)){
#         file_domain <- "noSelectHit"
#         updateButton(session, "do_domain_plot", disabled = TRUE)
#       } else {
#         ### check file extension
#         all_extension <- c("txt", "csv", "list", "domains", "architecture")
#         flag <- 0
#         for (i in 1:length(all_extension)){
# 
#           file_domain <- paste0(domain_path,
#                                 "/",
#                                 group,
#                                 ".",
#                                 all_extension[i])
#           if (file.exists(file_domain) == TRUE){
#             updateButton(session,
#                          "do_domain_plot",
#                          disabled = FALSE)
#             flag <- 1
#             break ()
#           }
#         }
# 
#         if (flag == 0){
#           file_domain <- "noFileInFolder"
#           updateButton(session,
#                        "do_domain_plot",
#                        disabled = TRUE)
#         }
#       }
#     }
#   }
# 
#   return (file_domain)
# }

get_multiplot_download_gc <- function(gene, input, interesting_features){
  x <- subset(significant_genes_gc,
              significant_genes_gc$geneID == gene)
  get_multiplot(x, input, interesting_features)
}

# Essential functions =========================================================

# get list with all taxanomy ranks --------------------------------------------
get_taxonomy_ranks <- function(){
  all_ranks <- list("Strain " = "05_strain",
                    "Species" = "06_species",
                    "Genus" = "10_genus",
                    "Family" = "14_family",
                    "Order" = "19_order",
                    "Class" = "23_class",
                    "Phylum" = "26_phylum",
                    "Kingdom" = "28_kingdom",
                    "Superkingdom" = "29_superkingdom",
                    "unselected" = "")
  return(all_ranks)
}

# Get name list ("data/taxonNamesReduced.txt") --------------------------------
get_name_list <- function (as_character, delete_duplicated){
  name_list <- as.data.frame(read.table("data/taxonNamesReduced.txt",
                                        sep = "\t",
                                        header = T,
                                        fill = TRUE))
  if (as_character) {
    name_list$fullName <- as.character(name_list$fullName)
    name_list$rank <- as.character(name_list$rank)
  }
  if (delete_duplicated){
    name_list <- name_list[!duplicated(name_list), ]
  }
  return(name_list)
}

# Get taxa list ("data/taxonomyMatrix.txt") -----------------------------------
get_taxa_list <- function(subset_taxa_check, subset_taxa){
  dt <- as.data.frame(read.table("data/taxonomyMatrix.txt",
                                 sep = "\t",
                                 header = T,
                                 stringsAsFactors = T))
  if (subset_taxa_check){
    dt <- dt[dt$abbrName  %in% subset_taxa, ]
  }
  return(dt)
}

# reverse string --------------------------------------------------------------
str_reverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")

# get last n characters from string x  ----------------------------------------
substr_right <- function(x, n){
  substr(x, nchar(x) - n + 1, nchar(x))
}

# check internet connection  --------------------------------------------------
has_internet <- function(){
  !is.null(curl::nslookup("r-project.org", error = FALSE))
}

# FUNCTION FOR CLUSTERING PROFILES --------------------------------------------
clustered_gene_list <- function(data, dist_method, cluster_method){

  # do clustering
  row.order <- hclust(dist(data, method = dist_method),
                      method = cluster_method)$order
  col.order <- hclust(dist(t(data), method = dist_method),
                      method = cluster_method)$order

  # re-order data accoring to clustering
  dat_new <- data[row.order, col.order]

  # return clustered gene ID list
  clustered_gene_ids <- as.factor(row.names(dat_new))
  clustered_gene_ids
}

# calculate percentage of present species -------------------------------------
calc_pres_spec <- function(taxa_md_data, taxa_count){
   # taxa_md_data = df("geneID",
   #                 "ncbiID",
   #                 "orthoID",
   #                 "var1",
   #                 "var2",
   #                 "paralog",
   #                 ....,
   #                 "supertaxon")
  taxa_md_data <- taxa_md_data[taxa_md_data$orthoID != "NA", ]

  # get geneID and supertaxon
  gene_id_supertaxon <- subset(taxa_md_data,
                             select = c("geneID",
                                        "supertaxon",
                                        "paralog",
                                        "abbrName"))
  # remove duplicated rows
  gene_id_supertaxon <- gene_id_supertaxon[!duplicated(gene_id_supertaxon), ]

  # remove NA rows from taxa_md_data
  taxa_md_data_no_na <- taxa_md_data[taxa_md_data$orthoID != "NA", ]

  # count present frequency of supertaxon for each gene
  gene_supertaxon_count <- plyr::count(taxa_md_data_no_na,
                                       c("geneID", "supertaxon"))

  # merge with taxa_count to get total number of species of each supertaxon
  # and calculate presSpec
  pres_spec_dt <- merge(gene_supertaxon_count,
                      taxa_count,
                      by = "supertaxon",
                      all.x = TRUE)

  spec_count <- plyr::count(gene_id_supertaxon, c("geneID",
                                              "supertaxon"))
  pres_spec_dt <- merge(pres_spec_dt,
                      spec_count, by = c("geneID",
                                       "supertaxon"))

  pres_spec_dt$presSpec <- pres_spec_dt$freq / pres_spec_dt$freq.y

  pres_spec_dt <- pres_spec_dt[pres_spec_dt$presSpec <= 1, ]
  pres_spec_dt <- pres_spec_dt[order(pres_spec_dt$geneID), ]
  pres_spec_dt <- pres_spec_dt[, c("geneID", "supertaxon", "presSpec")]

  # add absent supertaxon into pres_spec_dt
  gene_id_supertaxon <- subset(gene_id_supertaxon,
                             select = -c(paralog, abbrName))
  final_pres_spec_dt <- merge(pres_spec_dt,
                           gene_id_supertaxon,
                           by = c("geneID", "supertaxon"),
                           all.y = TRUE)
  final_pres_spec_dt$presSpec[is.na(final_pres_spec_dt$presSpec)] <- 0

  # remove duplicated rows
  final_pres_spec_dt <- final_pres_spec_dt[!duplicated(final_pres_spec_dt), ]

  # return final_pres_spec_dt
  return(final_pres_spec_dt)
}

###################### FUNCTIONS FOR RENDER UI ELEMENTS #######################

create_slider_cutoff <- function(id, title, start, stop, var_id){
  if(is.null(var_id)) return()
  if(var_id == ""){
    sliderInput(id, title,
                min = 1,
                max = 1,
                step = 0.025,
                value = 1,
                width = 200)
  } else {
    sliderInput(id, title,
                min = 0,
                max = 1,
                step = 0.025,
                value = c(start, stop),
                width = 200)
  }
}

update_slider_cutoff <- function(session, id, title, new_var, var_id){
  if(is.null(var_id) || var_id == "") return()
  
  updateSliderInput(session, id, title,
                    value = new_var,
                    min = 0,
                    max = 1,
                    step = 0.025)
}
