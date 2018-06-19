# unsorting function to keep user defined geneID order ------------------------
unsort_id <- function(data, order){
  data$geneID <- as.factor(data$geneID)
  if (order == FALSE) {
    # keep user defined geneID order
    data$geneID <- factor(data$geneID, levels = unique(data$geneID))
  }
  return(data)
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
get_name_list <- function (as_character, delete_duplicated) {
  name_list <- as.data.frame(read.table("data/taxonNamesReduced.txt",
                                        sep = "\t",
                                        header = T,
                                        fill = TRUE))
  if (as_character) {
    name_list$fullName <- as.character(name_list$fullName)
    name_list$rank <- as.character(name_list$rank)
  }
  if (delete_duplicated) {
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
  if (subset_taxa_check) {
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
  if (is.null(var_id)) return()
  if (var_id == "") {
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
  if (is.null(var_id) || var_id == "") return()
  
  updateSliderInput(session, id, title,
                    value = new_var,
                    min = 0,
                    max = 1,
                    step = 0.025)
}

create_plot_size <- function(id, title, value) {
  numericInput(id,
               title,
               min = 100,
               max = 3200,
               step = 50,
               value = value,
               width = 100)
}

create_text_size <- function(id, title, value, width) {
  numericInput(id,
               title,
               min = 3,
               max = 99,
               step = 1,
               value = value,
               width = width)
}

create_select_gene <- function(id, list, selected) {
  selectInput(id,
              "",
              list,
              selected = selected,
              multiple = TRUE,
              selectize = FALSE)
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
  if (is.null(var_id)) return()
  if (var_id == "") {
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
  if (is.null(var_id) || var_id == "") return()
  
  updateSliderInput(session, id, title,
                    value = new_var,
                    min = 0,
                    max = 1,
                    step = 0.025)
}

create_plot_size <- function(id, title, value) {
  numericInput(id,
               title,
               min = 100,
               max = 3200,
               step = 50,
               value = value,
               width = 100)
}

create_text_size <- function(id, title, value, width) {
  numericInput(id,
               title,
               min = 3,
               max = 99,
               step = 1,
               value = value,
               width = width)
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
  
  if (is.null(barplot)) {
    barplot <- textGrob("The selected domains are not found in the gene")
  }
  # if both variables are selected there are going to be 2 boxplots
  if (var == "Both") {
    pvalues <- unlist(pvalues, recursive = FALSE)
    
    p1 <- unlist(pvalues[1])
    p2 <- unlist(pvalues[2])
    # Check if the p_values should be printed
    if (input$show_p_value == TRUE) {
      info_p1 <- get_info_p_values(p1)
      info_p2 <- get_info_p_values(p2)
    }
    else{
      info_p1 <- " "
      info_p2 <- " "
    }
    # check if the significant plot should be highlighted
    if (input$highlight_significant == TRUE) {
      if (is.null(p1[1])) c1 <- "grey"
      else if (p1[length(p1)] < input$significance) c1 <- "indianred2"
      else c1 <- "grey"
      
      if (is.null(p2[1])) c2 <- "grey"
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
    if (input$show_p_value == TRUE) {
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
get_boxplot_gc <- function(in_group_df,
                           out_group_df,
                           var,
                           gene,
                           colour,
                           info,
                           input){
  
  if (var == input$var1_id) {
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
    geom_boxplot(stat = "boxplot",
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
  if (!("all" %in% input$interesting_features)) {
    ifeatures <- NULL
    for (x in input$interesting_features) {
      f <- subset(subdomain_df$feature,
                  startsWith(subdomain_df$feature, x))
      ifeatures <- append(ifeatures, f)
    }
    
    if (is.null(ifeatures)) return()
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
  seeds <- unique(subdomain_df$seedID)
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
  
  for (seed in seeds) {
    if (!is.null(data_in)) {
      in_g <- subset(in_group_domain_df, in_group_domain_df$seedID == seed)
      if (!empty(in_g)) {
        in_not_empty <- in_not_empty + 1
        in_group_features <-  plyr::count(in_g, "feature")
        for (i in 1:nrow(in_group_features)) {
          for (j in 1:nrow(data_in)) {
            if (data_in[j, 1] == in_group_features[i, 1]) {
              data_in[j, 2] <- data_in[j, 2] + in_group_features[i, 2]
            }
          }
        }
      }
    }
    
    if (!is.null(data_out)) {
      out_g <- {
        subset(out_group_domain_df, out_group_domain_df$seedID == seed)
      }
      
      if (!empty(out_g)) {
        out_not_empty <- out_not_empty + 1
        out_group_features <-  plyr::count(out_g, "feature")
        for (i in 1:nrow(out_group_features)) {
          for (j in 1:nrow(data_out)) {
            if (data_out[j, 1] == out_group_features[i, 1]) {
              data_out[j, 2] <- data_out[j, 2] + out_group_features[i, 2]
            }
          }
        }
      }
    }
  }
  
  if (!is.null(data_in)) {
    data_in$amount <- data_in$amount / in_not_empty
    data_in$type <- "in_group"
  }
  
  if (!is.null(data_out)) {
    data_out$amount <- data_out$amount / out_not_empty
    data_out$type <- "Out-Group"
  }
  
  if (is.null(data_in) & !is.null(data_out)) {
    data_barplot <- data_out
  }  else if (is.null(data_out) & !is.null(data_in)) {
    data_barplot <- data_in
  } else if (!is.null(data_in) & !is.null(data_out)) {
    data_barplot <- rbind(data_in, data_out)
  } else{
    data_barplot <- NULL
  }
  
  if (!is.null(data_barplot)) {
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
get_info_p_values <- function(p) {
  
  if (is.na(p[1])) info_p_values <- "not enough information"
  else if (length(p) == 1) {
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
get_significant_genes <- function(in_group,
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
      | length(list_selected_genes_gc) == 0) return()
  
  # load name List
  name_list <- get_name_list(TRUE, TRUE)
  
  # load list of unsorted taxa
  dt <- get_taxa_list(FALSE, subset_taxa)
  
  # Updateing of the Input ==================================================
  # if there is more than one element in the in_group
  # we look at the next common anchstor
  
  if (use_common_anchestor == TRUE) {
    ancestor <- get_common_ancestor(in_group, rank, name_list, dt)
    if (is.null(ancestor)) return()
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
  if (is.element("all", list_selected_genes_gc)) {
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
  for (gene in genes) {
    
    # creates Substet only for current Gene
    selected_gene_df <- subset(data_full, data_full$geneID == gene)
    # In- and Out-Group depending on the Gene
    in_group_df <- {
      subset(selected_gene_df,
             selected_gene_df$abbrName %in% selected_subset$abbrName)
    }
    out_group_df <- {
      subset(selected_gene_df,
             !(selected_gene_df$abbrName %in% selected_subset$abbrName))
    }
    out_group_df <- {
      subset(out_group_df, !out_group_df$fullName == in_select)
    }
    
    # Generate the p_values for the gene
    pvalues <- get_significant(in_group_df, out_group_df, var, gene, input)
    
    if (!is.null(pvalues)) {
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
      if (right_format_features) {
        new_row$databases <- list(get_prefix_features(features))
      }
      significantGenes <- rbind(significantGenes, new_row)
    }
  }
  if (nrow(significantGenes) != 0) {
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
get_selected_subset <- function(rank, in_group, name_list, dt){
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
  possible_ranks <- all_ranks[all_ranks >= rank]
  i <-  1
  if (length(in_group) == 1) rank <- substring(rank, 4)
  
  # find the common ancestor of all taxa in the in_group
  # and use it as in_group
  while (length(in_group) > 1 & i < length(possible_ranks)) {
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
  
  if (var == "Both") {
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
    if (is.null(pvalues1)) {
      
    }
    # check if the at last calculated p-value
    #is smaller than the significane level
    else if (pvalues1[length(pvalues1)] < significance_level) {
      significant <- TRUE
    }
    
    # analog to pvalues in the first variable
    if (is.null(pvalues2)) {
      
    } else if (pvalues2[length(pvalues2)] < significance_level) {
      significant <-  TRUE
    }
    
    # if the gene is interisting return the p_values
    if (significant) {
      pvalues <- list(pvalues1, pvalues2)
      return(pvalues)
    }
    # if the gene is not interisting return NULL
    else return(NULL)
    
  } else{
    # Check which variable is selected and get the p_values
    if (var == input$var1_id) {
      pvalues <- get_p_values(in_g$var1, out_g$var1, significance_level)
    }
    else {
      pvalues <- get_p_values(in_g$var2, out_g$var2, significance_level)
    }
    
    # Analog to getting the significance with both variables
    if (is.null(pvalues)) return(NULL)
    else if (is.nan(pvalues[length(pvalues)])) {
      return(NULL)
    } else if (pvalues[length(pvalues)] < significance_level) {
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
  if (length(var_in) == 0) return(NULL)
  else if (length(var_out) == 0) return(NULL)
  
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
