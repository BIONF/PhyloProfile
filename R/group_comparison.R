#' Group Comparison 
#' 
source("R/functions.R")

group_comparison_ui <- function(id){
  ns <- NS(id)
  fluidPage(
    sidebarPanel(
      uiOutput(ns("get_significant_genes")),
      bsPopover(
        "get_significant_genes",
        "",
        "Select gene to show the plots",
        "right"
      ),
      
      # checkboxInput(
      #   "add_gc_genes_custom_profile",
      #   strong(em("Add to Customized profile")),
      #   value = FALSE,
      #   width = NULL
      # ),
      
      
      
      popify(
        uiOutput(ns("features_of_interest_gc")),
        "",
        "This function is only use full if the features are
              saved in the right format: featuretype_featurename"
      ),
      
      actionButton(ns("gc_downloads"), "Download"),
      width = 3
    ),
    mainPanel(
      tags$style(
        HTML("#plots_gc { height:650px; overflow-y:scroll}")
      ),
      uiOutput(ns("plots_gc")),
      width = 9
    ),
    #* popup for handling the downloads to the Group comparison function -----
    bsModal(
      ns("gc_downloadsBS"),
      "Download",
      "gc_downloads",
      size = "small",

      h5(strong("Download the significant Genes")),
      downloadButton(ns("download_genes_gc"), "Download"),
      h5(""),
      uiOutput(ns("select_plots_to_download ")),
      downloadButton(ns("download_plots_gc"), "Download")
    )
  )
  


}

group_comparison <- function(input, output, session,
                             selected_in_group,
                             list_selected_genes,
                             rank_select,
                             var_name,
                             use_common_anchestor,
                             in_select,
                             demo_data,
                             anno_location,
                             file_domain_input,
                             subset_taxa,
                             get_data_filtered,
                             right_format_features,
                             get_domain_information,
                             plot_gc){
  # Dataframe with Information about the significant Genes --------------------
  # geneID | in_group| out_group | pvalues | features | databases | rank | var
  significant_genes <- NULL
  
  # Parameters for the plots in Group Comparison ------------------------------
  get_parameter_input_gc <- reactive ({
    input_data <- list("show_p_value" = input$show_p_value,
                       "highlight_significant" = input$highlight_significant,
                       "significance" = input$significance,
                       "var1_id" = input$var1_id,
                       "var2_id" = input$var2_id,
                       "x_size_gc" = input$x_size_gc,
                       "y_size_gc" = input$y_size_gc,
                       "interesting_features" = input$interesting_features,
                       "angle_gc" = input$angle_gc,
                       "legend_gc" = input$legend_gc,
                       "legend_size_gc" = input$legend_size_gc)
  })
  
  
  # List with all significant Genes -------------------------------------------
  output$get_significant_genes <- renderUI({
    ns <- session$ns
    plot_gc()

    isolate({
      significant_genes_gc <<- {
        get_significant_genes(selected_in_group(),
                              list_selected_genes(),
                              rank_select(),
                              var_name(),
                              use_common_anchestor(),
                              in_select(),
                              get_parameter_input_gc(),
                              demo_data(),
                              anno_location(),
                              file_domain_input(),
                              subset_taxa(),
                              get_data_filtered(),
                              right_format_features(),
                              get_domain_information())
      }
      if (!is.null(significant_genes_gc)){
        x <- as.vector(significant_genes_gc$geneID)
        choices <- c("all", x)

        # selected Gene
        selectInput(ns("selected_gene"), "Candidate gene(s):",
                    choices,
                    selected = choices[2],
                    multiple = FALSE)

      } else{
        # selected Gene
        selectInput(ns("selected_gene"), "Candidate gene(s):",
                    NULL,
                    selected = NULL,
                    multiple = FALSE)
      }
    })
  })

  # Generate output plots -----------------------------------------------------
  output$plots_gc <- renderUI({
    get_plots_gc()
  })
  
  # Select Feaures you want to see in the barplots (default: All) -------------
  output$features_of_interest_gc <- renderUI({
    ns <- session$ns
    input$selected_gene
    isolate({
      gene <- input$selected_gene
      if (!right_format_features()){
        selectInput(ns("interesting_features"), "Feature type(s) of interest:",
                    NULL,
                    selected = NULL,
                    multiple = TRUE,
                    selectize = FALSE)
      } else if (is.null(gene)){
        selectInput(ns("interesting_features"), "Feature type(s) of interest:",
                    NULL,
                    selected = NULL,
                    multiple = TRUE,
                    selectize = FALSE)
      } else if (gene == ""){
        selectInput(ns("interesting_features"), "Feature type(s) of interest:",
                    NULL,
                    selected = NULL,
                    multiple = TRUE,
                    selectize = FALSE)
      }
      else{
        choices <- c("all")
        if (gene == "all"){
          for (g in significant_genes_gc$geneID){
            x <- subset(significant_genes_gc,
                        significant_genes_gc$geneID == g)
            choices <- append(choices, unlist(x$databases))
          }
          # show each database only once
          choices <- choices[!duplicated(choices)]
        }
        else {
          
          x <- subset(significant_genes_gc,
                      significant_genes_gc$geneID == gene)
          
          choices <- append(choices, unlist(x$databases))
          
        }
        selectInput(ns("interesting_features"), "Feature type(s) of interest:",
                    choices,
                    selected = choices[1],
                    multiple = TRUE,
                    selectize = FALSE)
      }
    })
  })
  
  # Select Plots to download --------------------------------------------------
  output$select_plots_to_download  <- renderUI({
    ns <- session$ns
    plot_gc()
    isolate({
      if (!is.null(significant_genes_gc)){
        x <- as.vector(significant_genes_gc$geneID)
        choice <- c("all", x)
        gene <- subset(choice, choice == input$selected_gene)
        
        selectInput(ns("plots_to_download"), "Select Plots to download:",
                    choice,
                    selected = gene,
                    multiple = TRUE,
                    selectize = FALSE)
        
      } else{
        selectInput(ns("plots_to_download"), "Select Plots to download:",
                    NULL,
                    selected = NULL,
                    multiple = TRUE,
                    selectize = FALSE)
      }
    })
  })
  
  # Downloads for GroupCompairison --------------------------------------------
  # download list of significant genes
  output$download_genes_gc <- downloadHandler(
    filename = function(){
      c("significantGenes.out")
    },
    content = function(file){
      data_out <- significant_genes_gc$geneID
      write.table(data_out, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  # download file with the shown plots
  output$download_plots_gc <- downloadHandler(
    filename = "plotSignificantGenes.zip",
    content = function(file){
      genes <- input$plots_to_download
      
      if ("all" %in% genes){
        genes <- significant_genes_gc$geneID
      }
      
      fs <- c()
      #tmpdir <- tempdir()
      setwd(tempdir())
      
      for (gene in genes){
        path <- paste(gene, ".pdf", sep = "")
        fs <- c(fs, path)
        pdf(path)
        get_multiplot_download_gc(gene, get_parameter_input_gc(),
                                 input$interesting_features)
        dev.off()
      }
      zip(zipfile = file, files = fs)
    },
    contentType = "application/zip"
  )
  
  # observer for the download functions
  observe({
    if (is.null(selected_in_group())
        | length(list_selected_genes()) == 0){
      shinyjs::disable("download_plots_gc")
      shinyjs::disable("download_genes_gc")
    }else if (plot_gc() == FALSE){
      shinyjs::disable("download_plots_gc")
      shinyjs::disable("download_genes_gc")
    }else if (input$selected_gene== "") {
      shinyjs::disable("download_plots_gc")
      shinyjs::disable("download_genes_gc")
    }else{
      shinyjs::enable("download_plots_gc")
      shinyjs::enable("download_genes_gc")
    }
  })
  
  # Deciding which plots should be shown (Group Comparison) -------------------
  get_plots_gc <- reactive({
    gene <- as.character(input$selected_gene)
    plot_gc()
    if (is.null(significant_genes)) return()
    else if (gene == "all"){
      get_plot_output_list(significant_genes,
                           get_parameter_input_gc(),
                           input$interesting_features)
    }else{
      x <- {
        significant_genes_gc[significant_genes$geneID == gene, ]
      }
      if (nrow(x) == 0) return()
      get_plot_output_list(x, get_parameter_input_gc(), input$interesting_features)
    }
  })
}



# Get the list with all the significant genes and the dataset ---------------
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
                                  right_format_features,
                                  domains){
  if (is.null(in_group)
      | length(list_selected_genes_gc) == 0) return()
  
  # load name List
  name_list <- get_name_list(TRUE, TRUE)
  
  # load list of unsorted taxa
  dt <- get_taxa_list(FALSE, subset_taxa)

  # Updateing of the Input ----------------------------------------------------
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
                                domains)
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
get_p_values <- function(var_in, var_out, significance_level){
  
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
get_features <- function(selected_gene, domains){
  
  subdomain_df <- {
    subset(domain_df,
           substr(domain_df$seedID,
                  1,
                  nchar(as.character(selected_gene))) == selected_gene)
  }
  subdomain_df <- subdomain_df[!duplicated(subdomain_df), ]
  subdomain_df
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

get_multiplot_download_gc <- function(gene, input, interesting_features){
  x <- subset(significant_genes_gc,
              significant_genes_gc$geneID == gene)
  get_multiplot(x, input, interesting_features)
}