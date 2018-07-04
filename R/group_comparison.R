#' Group Comparison 
#'
#' @export
#' @param selected_in_group selected In-group (input$selected_in_group_gc)
#' @param selected_genes_list list of genes to calculate the plots for
#'  (input$list_selected_genes_gc)
#' @param main_rank in input and settings selected taxonomy rank
#'  (input$rank_select)
#' @param selected_variable variable(s) to claculate the plots for
#'  (input$var_name_gc)
#' @param use_common_ancestor boolean if the next common anchestor should be
#'  used (input$use_common_ancestor)
#' @param reference_taxon selected taxon name (input$in_select)
#' @param ncbi_id_list list of ncbi ids (from reactive fn "subset_taxa")
#' @param filtered_data full processed main data
#' (from reactive fn "get_data_filtered")
#' @param right_format_features boolean if the features have the right format
#' (input$right_format_features)
#' @param domain_information dataframe holding the domain input 
#' (from reactive fn "get_domain_information)
#' @param plot information if the plots should be generated (input$plot_gc)
#' @param parameter list of parameters needed to generate the plots 
#' (from reactive fn "get_parameter_input_gc")
#' @param changed_rank rank changed for the group comparison function
#'  (input$rank_select_gc)
#' @return list of candidate genes
#' @author Carla Mölbert {carla.moelbert@gmx.de}

if (!require("Matching")) install.packages("Matching") 
source("R/functions.R")

group_comparison_ui <- function(id){
  ns <- NS(id)
  fluidPage(
    sidebarPanel(
      withSpinner(uiOutput(ns("candidate_genes"))),
      bsPopover(
        "candidate_genes",
        "",
        "Select gene to show the plots",
        "right"
      ),
      withSpinner(uiOutput(ns("features_of_interest_ui"))),
      bsPopover(
        "features_of_interest_ui",
        "",
        "This function is only use full if the features are
              saved in the right format: featuretype_featurename"
      ),
      
      downloadButton(ns("download_plots"), "Download plots"),
      width = 3
    ),
    mainPanel(
      tags$style(
        HTML("#plots_ui { height:650px; overflow-y:scroll}")
      ),
      withSpinner(uiOutput(ns("plots_ui"))),
      width = 9
    )
  )
  
  
  
}

group_comparison <- function(input, output, session,
                             selected_in_group,
                             selected_genes_list,
                             main_rank,
                             selected_variable,
                             use_common_ancestor,
                             reference_taxon,
                             ncbi_id_list,
                             filtered_data,
                             right_format_features,
                             domain_information,
                             plot,
                             parameter,
                             changed_rank){
  # Dataframe for the significant Genes =========================
  #' contains geneID, in_group, out_group, pvalues, features, databases,
  #' rank, var
  candidate_genes <- reactiveValues(plots = NULL)
  
  # List with all candidate genes ========================= 
  output$candidate_genes <- renderUI({
    ns <- session$ns
    plot()
    
    isolate({
      candidate_genes$plots <- {
        get_significant_genes(selected_in_group(),
                              selected_genes_list(),
                              main_rank(),
                              selected_variable(),
                              use_common_ancestor(),
                              reference_taxon(),
                              parameter(),
                              ncbi_id_list(),
                              filtered_data(),
                              right_format_features(),
                              domain_information(),
                              changed_rank())
      }
      
      if (is.data.frame(candidate_genes$plots)) {
        significant_genes <- candidate_genes$plots
        x <- as.vector(significant_genes$geneID)
        choices <- c("all", x)
        
        selectInput(ns("selected_gene"), "Candidate gene(s):",
                    choices,
                    selected = choices[2],
                    multiple = FALSE)
        
      } else{
        selectInput(ns("selected_gene"), "Candidate gene(s):",
                    NULL,
                    selected = NULL,
                    multiple = FALSE)
      }
    })
  })
  
  # Output of the plots for the selected gene(s) =========================
  output$plots_ui <- renderUI({
    if (is.character(candidate_genes$plots)) return(candidate_genes$plots)
    get_plots()
  })
  
  # List with possible features for the selected gene =========================
  output$features_of_interest_ui <- renderUI({
    ns <- session$ns
    input$selected_gene
    isolate({
      gene <- input$selected_gene
      if (!right_format_features()) {
        selectInput(ns("interesting_features"), "Feature type(s) of interest:",
                    NULL,
                    selected = NULL,
                    multiple = TRUE,
                    selectize = FALSE)
      } else if (is.null(gene)) {
        selectInput(ns("interesting_features"), "Feature type(s) of interest:",
                    NULL,
                    selected = NULL,
                    multiple = TRUE,
                    selectize = FALSE)
      } else if (gene == "") {
        selectInput(ns("interesting_features"), "Feature type(s) of interest:",
                    NULL,
                    selected = NULL,
                    multiple = TRUE,
                    selectize = FALSE)
      }
      else{
        significant_genes <- candidate_genes$plots
        choices <- c("all")
        if (gene == "all") {
          for (current_gene in significant_genes$geneID) {
            subset_current_gene <- subset(significant_genes,
                                          significant_genes$geneID == current_gene)
            choices <- append(choices, unlist(subset_current_gene$databases))
          }
          # show each database only once
          choices <- choices[!duplicated(choices)]
        }
        else {
          
          subset_gene <- subset(significant_genes,
                                significant_genes$geneID == gene)
          
          choices <- append(choices, unlist(subset_gene$databases))
          
        }
        selectInput(ns("interesting_features"), "Feature type(s) of interest:",
                    choices,
                    selected = choices[1],
                    multiple = TRUE,
                    selectize = FALSE)
      }
    })
  })
  
  # download file with the shown plots =========================
  output$download_plots <- downloadHandler(
    filename = "plotSignificantGenes.zip",
    content = function(file){
      genes <- input$selected_gene
      significant_genes <- candidate_genes$plots
      if ("all" %in% genes) {
        genes <- significant_genes$geneID
      }
      
      fs <- c()
      #tmpdir <- tempdir()
      setwd(tempdir())
      
      for (gene in genes) {
        path <- paste(gene, ".pdf", sep = "")
        fs <- c(fs, path)
        pdf(path)
        get_plots_to_download(gene,
                              parameter(),
                              input$interesting_features,
                              significant_genes)
        dev.off()
      }
      zip(zipfile = file, files = fs)
    },
    contentType = "application/zip"
  )
  
  #' observer for the download functions
  observe({
    if (is.null(selected_in_group())
        | length(selected_genes_list()) == 0) {
      shinyjs::disable("download_plots")
    }else if (plot() == FALSE) {
      shinyjs::disable("download_plots")
    }else if (input$selected_gene == "") {
      shinyjs::disable("download_plots")
    }else{
      shinyjs::enable("download_plots")
    }
  })
  
  # Deciding which plots will be shown =========================
  get_plots <- reactive({
    input$interesting_features
    gene <- as.character(input$selected_gene)
    plot()
    if (is.null(candidate_genes$plots)) return()
    significant_genes <- candidate_genes$plots
    if (gene == "all") {
      plot_output_list <- get_plot_output_list(significant_genes,
                                               parameter(),
                                               input$interesting_features)
    }else{
      gene_info <- {
        significant_genes[significant_genes$geneID == gene, ]
      }
      if (nrow(gene_info) == 0) return()
      plot_output_list <- get_plot_output_list(gene_info,
                                               parameter(),
                                               input$interesting_features)
    }
    #' List with all plots that will be shown
    return(plot_output_list)
  })
  
  # List of genes for the customized profile =========================
  gene_list <- reactive({
    if (!is.null(candidate_genes$plots)) {
      significant_genes <- candidate_genes$plots
      return(significant_genes$geneID)
    } 
  })
  
  return(gene_list)
}

# FUNCTIONS ===================================================================

#' Get the dataframe with the significant genes -------------------------------
#' @export
#' @param in_group list of taxa 
#' @param selected_genes_list list of genes
#' @param rank selected taxonamy rank
#' @param var variable for which to  calculate the significance
#' @param use_common_ancestor boolean if the common anchestor should be used
#' @param reference_taxon taxon which is used as reference
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#'  "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @param ncbi_id_list list of nxbi ids 
#' @param data_full full processed main data
#' @param right_format_features boolean if the features have the right format
#' @param  domains dataframe holding the domain input
#' @return dataframe with the significant genes
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_significant_genes <- function(in_group,
                                  selected_genes_list,
                                  main_rank,
                                  var,
                                  use_common_ancestor,
                                  reference_taxon,
                                  parameters,
                                  ncbi_id_list,
                                  data_full,
                                  right_format_features,
                                  domains,
                                  changed_rank){
  if (is.null(in_group) | length(selected_genes_list) == 0) return()
  
  
  name_list <- get_name_list(TRUE, TRUE) # load name List
  taxa_list <- get_taxa_list(FALSE, ncbi_id_list) # load list of unsorted taxa
  
  if (is.null(changed_rank)) rank <- main_rank
  else rank <- changed_rank
  
  #' Get the rank and the in-group --------------------------------------------
  #' if there is more than one element in the in_group -> use common anchestor
  if (use_common_ancestor == TRUE) {
    ancestor <- get_common_ancestor(in_group, rank,
                                    name_list, taxa_list, ncbi_id_list)
    if (is.null(ancestor)) return("No common anchestor found")
    in_group <- ancestor[1]
    rank <- ancestor[2]
  } else{
    rank <- substring(rank, 4)
  }
  if (is.na(rank)) { return("No common anchestor found")}
  
  #' provide the empty data frame ---------------------------------------------
  significant_genes_df <- data.frame(
    geneID = character(),
    in_group = I(list()),
    out_group = I(list()),
    pvalues = I(list()),
    features = I(list()),
    databases = I(list()))
  
  #' Get the list of genes to look at ----------------------------------------- 
  if (is.element("all", selected_genes_list)) {
    genes <- data_full$geneID
    genes <- genes[!duplicated(genes)]
  } else {
    genes <- selected_genes_list
  }
  genes <- sort(genes)
  
  #' Subset depending on the rank and the in_group ----------------------------
  selected_subset <- get_selected_subset(rank, in_group, name_list, taxa_list)
  selected_subset <- subset(selected_subset,
                            !selected_subset$fullName == reference_taxon)
  
  #' Check for each gene if it is significant ---------------------------------
  for (gene in genes) {
    #' Processing the dataframes for in- and out-group
    selected_gene_df <- subset(data_full, data_full$geneID == gene)
    
    in_group_df <- {
      subset(selected_gene_df,
             selected_gene_df$abbrName %in% selected_subset$abbrName)
    }
    out_group_df <- {
      subset(selected_gene_df,
             !(selected_gene_df$abbrName %in% selected_subset$abbrName))
    }
    out_group_df <- {
      subset(out_group_df, !out_group_df$fullName == reference_taxon)
    }
    
    #' Generate and check the p_values for the gene ---------------------------
    pvalues <- get_p_values(in_group_df, out_group_df, var, gene, parameters)
    
    
    if (!is.null(pvalues)) {
      new_row <- data.frame(geneID = gene,
                            in_group = NA,
                            out_group = NA,
                            pvalues = NA,
                            features = NA)
      new_row$in_group <- list(in_group_df)
      new_row$out_group <- list(out_group_df)
      new_row$pvalues <- list(pvalues)
      features  <- get_features(gene, domains)
      new_row$features <- list(features)
      if (right_format_features) {
        new_row$databases <- list(get_prefix_features(features))
      }
      significant_genes_df <- rbind(significant_genes_df, new_row)
    }
  }
  #' return the significant genes ---------------------------------------------
  if (nrow(significant_genes_df) != 0) {
    significant_genes_df$var <- var
    significant_genes_df$rank <- rank
    return(significant_genes_df)
  }
  else {
    return("No candidate genes found")
  }
}

#' Get the database for each feature in a specific gene-------------------------
#' @export
#' @param data contains "seedID", "orthoID", "feature", "start", "end"
#' @return list of prefixes for the features
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_prefix_features <- function(data){
  features <- data$feature
  choices <- gsub("_.*", "", features)
  choices <- choices[!duplicated(choices)]
  return(choices)
}


#' Get the subset depending on the choosen rank ------------------------------
#' @export
#' @param rank selected rank
#' @param in_group list of taxa
#' @param name_list contains "ncbiID", "fullName", "rank", "parentID"
#' @param taxa_list contains "abbrName, "ncbiID", fullName", "strain", "genus"
#' @return list of prefixes for the features
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_selected_subset <- function(rank, in_group, name_list, taxa_list){
  #' Look if the fullName is in the in_group
  name_list$fullName <- as.character(name_list$fullName)
  name_list_rank <- subset(name_list, name_list$rank == rank)
  in_group_subset <- subset(name_list, name_list$fullName %in% in_group)
  
  #' Look if it has the right rank
  selected_subset <- taxa_list[taxa_list[, rank] %in% in_group_subset$ncbiID, ]
  
  return(selected_subset)
}

#' Generate the in_group ------------------------------------------------------
#' @export
#' @param in_group list of taxa 
#' @param rank selected rank
#' @param name_list contains "ncbiID", "fullName", "rank", "parentID"
#' @param selected_in_group contains "abbrName, "ncbiID", fullName", "strain",
#'  "genus"..
#' @return common anchestor 
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_common_ancestor <- function(in_group,
                                rank,
                                name_list,
                                selected_in_group, 
                                ncbi_id_list){
  
  all_ranks <- get_taxonomy_ranks()
  
  selected_in_group <- {
    selected_in_group[!duplicated(selected_in_group), ]
  }
  
  #' ranks were all elements of the in_group might be in the same taxon
  possible_ranks <- all_ranks[all_ranks >= rank]
  position <-  1
  if (length(in_group) == 1) rank <- substring(rank, 4)
  
  #' find the common ancestor of all taxa in the in_group ---------------------
  while (length(in_group) > 1 & position < length(possible_ranks)) {
    
    current_rank <- substring(possible_ranks[position], 4)
    next_rank <- substring(possible_ranks[position + 1], 4)
    #' dataframe with all elements with fitting rank
    df_in_group <- subset(name_list, name_list$rank == current_rank)
    
    #' subset of df_in_group  with elements that belong to the in_group
    df_in_group <- subset(df_in_group, df_in_group$fullName %in% in_group)
    
    #' get all elements which could belong to the in-group
    possible_in_group <- subset(selected_in_group,
                                select = c(current_rank, next_rank))
    possible_in_group <- {
      possible_in_group[possible_in_group[,current_rank]
                        %in% df_in_group$ncbiID, ]
    }
    possible_in_group <- possible_in_group[!duplicated(possible_in_group), ]
    
    #' only consider elements that have the next higher rank
    subset_next_rank <- taxa_select_gc(next_rank, ncbi_id_list)
    subset_next_rank <- subset_next_rank[!duplicated(subset_next_rank), ]
    subset_next_rank <- {
      subset(subset_next_rank,
             subset_next_rank$ncbiID %in% possible_in_group[, next_rank] )
    }
    in_group <- subset_next_rank$fullName
    position <- position + 1
    rank <- next_rank
  }
  
  #' Return the in-group and the rank -----------------------------------------
  if (position > length(possible_ranks)) return()
  return(c(in_group, rank))
}

#' Decide if the gene is significant ------------------------------------------
#' @export
#' @param in_group contains "supertaxon", "geneID", "ncbiID", "orthoID",
#'  "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#'  "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2" 
#' @param out_group  as in-group but with information containing the out-group
#' @param variable variable(s) to claculate the plots for
#' @param gene gene to calculate the p-values for
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#'  "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @return return the pvalues
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_p_values <- function(in_group, out_group, variable, gene, parameters){
  
  significance_level <- parameters$significance
  
  #' get the p-values for both variables --------------------------------------
  if (variable == "Both") {
    var1 <- parameters$var1
    var2 <- parameters$var2
    
    significant <-  FALSE
    
    pvalues1 <- calculate_p_value(in_group$var1,
                                  out_group$var1,
                                  significance_level)
    pvalues2 <- calculate_p_value(in_group$var2,
                                  out_group$var2,
                                  significance_level)
    
    #' check if the gene has a significant difference -------------------------
    #' in the distribioution of In- and Out-Group
    #' in at least one of the variables
    
    if (is.null(pvalues1)) {
      #' if there is not enough data to calculate a p-value the gene is 
      #' considered as not intereisting
    }
    else if (pvalues1[length(pvalues1)] < significance_level) {
      #' if the last p-value is smaller than the significance level
      #' the gene is considered as significant
      significant <- TRUE
    }
    
    #' analog to pvalues in the first variable
    if (is.null(pvalues2)) {
      
    } else if (pvalues2[length(pvalues2)] < significance_level) {
      significant <-  TRUE
    }
    
    #' Return the conclusion --------------------------------------------------
    #' if the gene is interisting return the p_values 
    if (significant) {
      pvalues <- list(pvalues1, pvalues2)
      return(pvalues)
    }
    #' if the gene is not interisting return NULL
    else return(NULL)
    
    #' get the p-values for one variable --------------------------------------
  } else{
    #' Check which variable is selected and get the p_values
    if (variable == parameters$var1_id) {
      pvalues <- calculate_p_value(in_group$var1,
                                   out_group$var1,
                                   significance_level)
    }
    else {
      pvalues <- calculate_p_value(in_group$var2,
                                   out_group$var2,
                                   significance_level)
    }
    
    #' Analog to getting the significance with both variables
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


#' calculate the p_values -----------------------------------------------------
#' @export
#' @param var_in list of values for the variable concerning the in-group
#' @param var_out list of values for the variable concerning the out-group
#' @param significance_level 
#' @return return the pvalues
#' @author Carla Mölbert (carla.moelbert@gmx.de)
calculate_p_value <- function(var_in, var_out, significance_level){
  
  #' delete all entrys that are NA
  var_in <- var_in[!is.na(var_in)]
  var_out <- var_out[!is.na(var_out)]
  
  #' if there is no data in one of the groups the p-value is NULL
  if (length(var_in) == 0) return(NULL)
  else if (length(var_out) == 0) return(NULL)
  else{
    #' * Kolmogorov-Smirnov Test ----------------------------------------------
    #' H0 : The two samples have the same distribution
    #' package "Matching" is required 
    ks <- ks.boot(var_in, var_out, alternative = "two.sided")
    p_value <- ks$ks.boot.pvalue # probabilitiy to recet H0 if it is correct
    
    if (p_value < significance_level) pvalue <- c(p_value)
    
    else {
      #' * Wilcoxon-Mann-Whitney Test -----------------------------------------
      #' H0: the samples have the same location parameters
      wilcox <- wilcox.test(var_in,
                            var_out,
                            alternative = "two.sided",
                            exact = FALSE)
      p_value_wilcox <- wilcox$p.value
      pvalue <- c(p_value, p_value_wilcox)
    }
    #' return the calculated pvalues ------------------------------------------
    return(pvalue)
  }
}

#' get the list with all the features in the gene -----------------------------
#' @export
#' @param selected_gene gene to get the feartures for
#' @param domains contains "seedID", "orthoID", "feature", "start",   "end"
#' @return dataframe for the specific gene containing "seedID",  "orthoID",
#' "feature", "start",   "end"
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_features <- function(selected_gene, domains){
  subset_domains <- {
    subset(domains,
           substr(domains$seedID,
                  1,
                  nchar(as.character(selected_gene))) == selected_gene)
  }
  subset_domains <- subset_domains[!duplicated(subset_domains), ]
  return(subset_domains)
}


#' Generate the list with all plots -------------------------------------------
#' @export
#' @param genes list of genes
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#' "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @param interesting_features list of databases to take the features from
#' @return list with all plots
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_plot_output_list <- function(genes, parameters, interesting_features) {
  # if we don't have parameters we can not generate plots
  if (is.null(genes)) return()
  # Insert plot output objects the list
  plot_output_list <- lapply(1:nrow(genes), function(i) {
    plotname <- paste(genes[i, 1])
    plot_output_object <- renderPlot(get_multiplot(genes[i, ], parameters,
                                                   interesting_features),
                                     height = 650, width = 700)
  })
  do.call(tagList, plot_output_list) # needed to display properly.
  return(plot_output_list)
}

#' Put the plots for one spicific gene in one multiplot -----------------------
#' @export
#' @param gene_info contains "geneID",  "in_group",  "out_group", "pvalues",
#' "features",  "databases", "var", "rank"
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#' "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @param interesting_features list of databases to take the features from
#' @return grid arrange with the plots that should be shown for this gene
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_multiplot <- function(gene_info, parameters, interesting_features){
  #' Sorting the information to the selected gene ---------------------------- 
  gene <- as.character(gene_info$geneID)
  in_group <- as.data.frame(gene_info$in_group)
  out_group <- as.data.frame(gene_info$out_group)
  features <- as.data.frame(gene_info$features)
  pvalues <- gene_info$pvalues
  var <- gene_info$var
  
  #' Get the barplot ---------------------------------------------------------
  barplot <-  get_barplot_gc(gene,
                             in_group,
                             out_group,
                             features,
                             parameters,
                             interesting_features)
  
  if (is.null(barplot)) {
    barplot <- textGrob("The selected domains are not found in the gene")
  }
  
  #' Get the boxplots  for two variables  -------------------------------------
  if (var == "Both") {
    #' Get information about pvalues
    pvalues <- unlist(pvalues, recursive = FALSE)
    p_value1 <- unlist(pvalues[1])
    p_value2 <- unlist(pvalues[2])
    
    #' Check if the p_values should be printed
    if (parameters$show_p_value == TRUE) {
      info_p_value1 <- get_info_p_values(p_value1)
      info_p_value2 <- get_info_p_values(p_value2)
    }
    else{
      info_p_value1 <- " "
      info_p_value2 <- " "
    }
    
    #' Get information about the plot colour 
    if (parameters$highlight_significant == TRUE) {
      if (is.null(p_value1[1])) colour1 <- "grey"
      else if (p_value1[length(p_value1)] < parameters$significance) {
        colour1 <- "indianred2"
      }
      else colour1 <- "grey"
      
      if (is.null(p_value2[1])) colour2 <- "grey"
      else if (p_value2[length(p_value2)] < parameters$significance) {
        colour2 <- "indianred2"
      }
      else colour2 <- "grey"
    }
    else{
      colour1 <- "grey"
      colour2 <- "grey"
    }
    
    #' Generate the boxplots 
    boxplot1 <- get_boxplot_gc(in_group,
                               out_group,
                               parameters$var1_id,
                               gene,
                               colour1,
                               info_p_value1, parameters)
    
    boxplot2 <- get_boxplot_gc(in_group,
                               out_group,
                               parameters$var2_id,
                               gene,
                               colour2,
                               info_p_value2, parameters)
    
    plots <- grid.arrange(textGrob(gene),
                          arrangeGrob(boxplot1, boxplot2, ncol = 2),
                          barplot,
                          heights = c(0.02, 0.45, 0.458), ncol = 1)
  }else {
    #' get the boxplot if one varibale is selected  ---------------------------
    p <- unlist(pvalues)
    
    #' Check if the p_values should be printed
    if (parameters$show_p_value == TRUE) {
      info <- get_info_p_values(p)
    }else{
      info <- " "
    }
    
    #' Generate the plot 
    boxplot <- get_boxplot_gc(in_group,
                              out_group,
                              var,
                              gene,
                              "grey",
                              info,
                              parameters)
    
    plots <- grid.arrange(textGrob(gene),
                          boxplot,
                          barplot,
                          heights = c(0.02, 0.45, 0.458), ncol = 1)
  }
  
  #' return the plots --------------------------------------------------------
  return(plots)
}

#' Create a Boxplot -----------------------------------------------------------
#' @export
#' @param in_group_df contains "supertaxon", "geneID", "ncbiID", "orthoID",
#'  "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#'  "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2" 
#' @param out_group_df contains "supertaxon", "geneID", "ncbiID", "orthoID",
#'  "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#' "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2" 
#' @param var variable to consider in the boxplot
#' @param  gene gene for which the plot is generated
#' @param  colour colour of the boxes
#' @param  info info about the p-values
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#' "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @return boxplot
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_boxplot_gc <- function(in_group_df,
                           out_group_df,
                           var,
                           gene,
                           colour,
                           info,
                           parameters){
  #' pre-processing the data for the boxplot ---------------------------------
  if (var == parameters$var1_id) {
    in_group <- in_group_df$var1
    out_group <- out_group_df$var1
  }
  else if (var == parameters$var2_id) {
    in_group <- in_group_df$var2
    out_group <- out_group_df$var2
  }
  
  length_in_group <- length(in_group)
  length_out_group <- length(out_group)
  
  in_group <- as.data.frame(in_group)
  names(in_group)[1] <- paste("values")
  in_group$group <- "in_group"
  
  out_group <- as.data.frame(out_group)
  names(out_group)[1] <- paste("values")
  out_group$group <- "Out-Group"
  
  data_boxplot <- rbind(in_group, out_group)
  data_boxplot <- data_boxplot[complete.cases(data_boxplot), ]
  
  names <- c(paste("In-Group \n n=", length_in_group, sep = ""),
             paste("Out-Group \n n=", length_out_group, sep = ""))
  
  #' Generate the boxplot -----------------------------------------------------
  boxplot_gc <- ggplot(data_boxplot, aes(group, values)) +
    geom_boxplot(stat = "boxplot",
                 position = position_dodge(),
                 width = 0.5,
                 fill = colour) +
    labs(x = "", y = var, caption = paste(info)) +
    scale_x_discrete(labels = names) +
    theme_minimal()
  
  boxplot_gc <- boxplot_gc +
    theme(axis.text.x = element_text(size = parameters$x_size_gc, hjust = 1),
          axis.text.y = element_text(size = parameters$y_size_gc),
          axis.title.y = element_text(size = parameters$y_size_gc))
  #' return the boxplot -------------------------------------------------------
  return(boxplot_gc)
}


#' Create a Barplot -----------------------------------------------------------
#' @export
#' @param  selected_gene gene for which the plot is generated
#' @param in_group_df contains "supertaxon", "geneID", "ncbiID", "orthoID",
#'  "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#'  "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2" 
#' @param out_group_df contains "supertaxon", "geneID", "ncbiID", "orthoID",
#'  "var1", "var2", "paralog", "abbrName", "taxonID", "fullname",
#' "supertaxonID", "rank", "category", "presSpec", "mVar1", "mVar2" 
#' @param features contains "seedID",  "orthoID", "feature", "start",   "end"
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#' "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @param interesting_features list of databases for which the features should 
#' be included
#' @return barplot
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_barplot_gc <- function(selected_gene,
                           in_group, out_group,
                           features, parameters, interesting_features){
  
  #' information about the features ------------------------------------------
  features$feature <- as.character(features$feature)
  #' only show features that interest the user
  if (!("all" %in% interesting_features)) {
    features_list <- NULL
    for (feature in interesting_features) {
      subset_features <- subset(features$feature,
                                startsWith(features$feature, feature))
      features_list <- append(features_list, subset_features)
    }
    
    if (is.null(features_list)) return()
    #' only keep rows in which the feature begins with a element out of the
    #' interesing Features
    features_list <- features_list[!duplicated(features_list)]
    features <- subset(features, features$feature %in% features_list)
  }
  
  #' part in in_group and out-group  ------------------------------------------ 
  in_group_domain_df  <-  {
    subset(features, features$orthoID %in% in_group$orthoID)
  }
  out_group_domain_df <- {
    subset(features, features$orthoID %in% out_group$orthoID)
  }
  
  #' get the dataframes for in- and out-group ---------------------------------
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
  
  #' Get the values for the boxplot ------------------ ------------------------ 
  seeds <- unique(features$seedID)
  in_not_empty <- 0
  out_not_empty <- 0
  
  #' Count for each feature how often it is present in each seed --------------
  for (seed in seeds) {
    #' count the features in the in-group  
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
    
    #' count the featueres in the out-group 
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
  
  #' Calculate the average of appearances for each feature --------------------
  if (!is.null(data_in)) {
    data_in$amount <- data_in$amount / in_not_empty
    data_in$type <- "In-Group"
  }
  
  if (!is.null(data_out)) {
    data_out$amount <- data_out$amount / out_not_empty
    data_out$type <- "Out-Group"
  }
  
  #' Get the data for teh barplot --------------------------------------------- 
  if (is.null(data_in) & !is.null(data_out)) {
    data_barplot <- data_out
  }  else if (is.null(data_out) & !is.null(data_in)) {
    data_barplot <- data_in
  } else if (!is.null(data_in) & !is.null(data_out)) {
    data_barplot <- rbind(data_in, data_out)
  } else{
    data_barplot <- NULL
  }
  
  #' generate the barplot -----------------------------------------------------
  if (!is.null(data_barplot)) {
    barplot_gc <- ggplot(data_barplot,
                         aes(x = feature, y = amount, fill = type ),
                         main  = " ") +
      geom_bar(stat = "identity", position = position_dodge(), width = 0.5) +
      scale_fill_grey() +
      labs(x = " ", y = "Average instances per protein", fill = "Group") +
      theme_minimal()
    
    barplot_gc <- barplot_gc +
      theme(axis.text.x = element_text(size = parameters$x_size_gc,
                                       angle = parameters$angle_gc, hjust = 1),
            axis.text.y = element_text(size = parameters$y_size_gc),
            axis.title.y = element_text(size = parameters$y_size_gc),
            legend.position = parameters$legend_gc,
            legend.text = element_text(size = parameters$legend_size_gc ),
            legend.title = element_text(size = parameters$legend_size_gc))
    
    #' return the barplot -----------------------------------------------------
    return(barplot_gc)
  } else (return(NULL))
}

#' Get the p_values to print under the plot -----------------------------------
#' @export
#' @param pvalues list contianing the p-values for a specific gene
#' @return string containing the information about the p-values
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_info_p_values <- function(pvalues) {
  
  if (is.na(pvalues[1])) info_p_values <- "not enough information"
  else if (length(pvalues) == 1) {
    info_p_values <- paste("Kolmogorov-Smirnov-Test:",
                           as.character(pvalues[1]), sep = " " )
  }
  else{
    info_p_values1 <- paste("Kolmogorov-Smirnov-Test:",
                            as.character(pvalues[1]), sep = " " )
    info_p_values2 <- paste("Wilcoxon-Mann-Whitney-Test: ",
                            as.character(round(pvalues[2], 6)), sep = " ")
    info_p_values <- paste(info_p_values1, info_p_values2, sep = "\n")
  }
  
  info_p_values <- paste("p_values:", info_p_values, sep = "\n")
  return(info_p_values)
}

#' get the plots to download --------------------------------------------------
#' @export
#' @param  gene gene for which the plot is generated
#' @param interesting_features list of databases for which the features should 
#' be included
#' @param parameters contains "show_p_value","highlight_significant",
#' "significance", "var1_id", "var2_id", "x_size_gc", "y_size_gc",
#' "interesting_features", "angle_gc", "legend_gc", "legend_size_gc"
#' @return arrange grop containing the plots for this gene
#' @author Carla Mölbert (carla.moelbert@gmx.de)
get_plots_to_download <- function(gene,
                                  parameters,
                                  interesting_features,
                                  significant_genes){
  info_gene <- subset(significant_genes,
                      significant_genes$geneID == gene)
  return(get_multiplot(info_gene, parameters, interesting_features))
}