#' Gene age estimation plot
#'
#' @export
#' @param data for gene age plot (from reactive fn gene_ageDf)
#' @param gene_age_width plot width (from input$gene_age_width)
#' @param gene_age_height plot width (from input$gene_age_height)
#' @param gene_age_text text size (from input$gene_age_text)
#' @return list of genes of a selected age
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

plot_gene_age_ui <- function(id) {
  ns <- NS(id)
  tagList(
    column(
      2,
      downloadButton(ns("gene_age_plot_download"),
                     "Download plot", class = "butDL"),
      tags$head(
        tags$style(HTML(
          ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
      )
    ),
    column(
      10,
      uiOutput(ns("gene_age.ui")),
      br(),
      em(
        h6("01_Species; 02_Family; 03_Class; 04_Phylum;
            05_Kingdom; 06_Superkingdom; 07_Last universal common ancestor;
           Undef_Genes have been filtered out")
      )
    ),
    hr(),
    column(
      4,
      downloadButton(ns("gene_age_table_download"), "Download gene list")
    ),
    column(
      8,
      tableOutput(ns("gene_age.table"))
    )
  )
}

plot_gene_age <- function(input, output, session,
                          data,
                          gene_age_width, gene_age_height, gene_age_text) {

  # render gene age plot -------------------------------------------------------
  output$gene_agePlot <- renderPlot({
    if (is.null(data())) return()
    create_gene_age_plot(gene_age_plotDf(data()), gene_age_text())
  })

  output$gene_age.ui <- renderUI({
    ns <- session$ns
    withSpinner(
      plotOutput(
        ns("gene_agePlot"),
        width = 600 * gene_age_width(),
        height = 150 * gene_age_height(),
        click = ns("plot_click_gene_age")
      )
    )
  })

  output$gene_age_plot_download <- downloadHandler(
    filename = function() {
      "gene_age_plot.pdf"
    },
    content = function(file) {
      ggsave(
        file,
        plot = create_gene_age_plot(gene_age_plotDf(data()), gene_age_text()),
        width = 600 * gene_age_width() * 0.056458333,
        height = 150 * gene_age_height() * 0.056458333,
        units = "cm", dpi = 300, device = "pdf"
      )
    }
  )

  # render genAge.table based on clicked point on gene_agePlot ----------------
  selectedgene_age <- reactive({
    if (is.null(data())) return()
    selected_gene <- get_selected_gene_age(data(), input$plot_click_gene_age$x)
    return(selected_gene)
  })

  output$gene_age.table <- renderTable({
    if (is.null(data())) return()
    if (is.null(input$plot_click_gene_age$x)) return()

    data <- as.data.frame(selectedgene_age())
    data$number <- rownames(data)
    colnames(data) <- c("geneID", "No.")
    data <- data[, c("No.", "geneID")]
    data
  })

  output$gene_age_table_download <- downloadHandler(
    filename = function(){
      c("selectedGeneList.out")
    },
    content = function(file){
      data_out <- selectedgene_age()
      write.table(data_out, file,
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)
    }
  )

  return(selectedgene_age)
}

#' calculate gene age function
#' @param subset_taxa list of input taxa (from reactive fn "subset_taxa")
#' @param filtered_data processed data (from reactive fn "get_data_filtered")
#' @param rank_select selected taxonomy rank (input$rank_select)
#' @param in_select selected taxon name (input$in_select)
#' @param var1_cutoff
#' @param var2_cutoff
#' @param percent_cutoff

estimate_gene_age <- function(subset_taxa, filtered_data,
                              rank_select, in_select,
                              var1_cutoff, var2_cutoff, percent_cutoff){
  rankList <- c("family",
                "class",
                "phylum",
                "kingdom",
                "superkingdom",
                "root")

  # get selected (super)taxon ID
  rankName <- substr(rank_select, 4, nchar(rank_select))
  taxa_list <- get_name_list(FALSE, FALSE)
  superID <- {
    as.numeric(
      taxa_list$ncbiID[taxa_list$fullName == in_select
                       & taxa_list$rank == rankName]
    )
  }

  # full non-duplicated taxonomy data
  Dt <- get_taxa_list(FALSE, subset_taxa)

  # subset of taxonomy data, containing only ranks from rankList
  subDt <- Dt[, c("abbrName", rankList)]

  # get (super)taxa IDs for one of representative species
  # get all taxon info for 1 representative
  first_line <- Dt[Dt[, rankName] == superID, ][1, ]
  sub_first_line <- first_line[, c("abbrName", rankList)]

  # compare each taxon ncbi IDs with selected taxon
  # and create a "category" data frame
  catDf <- data.frame("ncbiID" = character(),
                      "cat" = character(),
                      stringsAsFactors = FALSE)
  for (i in 1:nrow(subDt)) {
    cat <- subDt[i, ] %in% sub_first_line
    cat[cat == FALSE] <- 0
    cat[cat == TRUE] <- 1
    cat <- paste0(cat, collapse = "")
    catDf[i, ] <- c(as.character(subDt[i, ]$abbrName), cat)
  }

  # get main input data
  mdData <- droplevels(filtered_data)
  mdData <- mdData[, c("geneID",
                       "ncbiID",
                       "orthoID",
                       "var1",
                       "var2",
                       "presSpec")]

  ### add "category" into mdData
  mdDataExtended <- merge(mdData,
                          catDf,
                          by = "ncbiID",
                          all.x = TRUE)

  mdDataExtended$var1[mdDataExtended$var1 == "NA"
                      | is.na(mdDataExtended$var1)] <- 0
  mdDataExtended$var2[mdDataExtended$var2 == "NA"
                      | is.na(mdDataExtended$var2)] <- 0

  # remove cat for "NA" orthologs
  # and also for orthologs that do not fit cutoffs
  if (nrow(mdDataExtended[mdDataExtended$orthoID == "NA"
                          | is.na(mdDataExtended$orthoID), ]) > 0) {
    mdDataExtended[mdDataExtended$orthoID == "NA"
                   | is.na(mdDataExtended$orthoID), ]$cat <- NA
  }

  mdDataExtended <- mdDataExtended[complete.cases(mdDataExtended), ]

  # filter by %specpres, var1, var2 ..
  mdDataExtended$cat[mdDataExtended$var1 < var1_cutoff[1]] <- NA
  mdDataExtended$cat[mdDataExtended$var1 > var1_cutoff[2]] <- NA
  mdDataExtended$cat[mdDataExtended$var2 < var2_cutoff[1]] <- NA
  mdDataExtended$cat[mdDataExtended$var2 > var2_cutoff[2]] <- NA
  mdDataExtended$cat[mdDataExtended$presSpec < percent_cutoff[1]] <- NA
  mdDataExtended$cat[mdDataExtended$presSpec > percent_cutoff[2]] <- NA

  mdDataExtended <- mdDataExtended[complete.cases(mdDataExtended), ]

  ### get the furthest common taxon with selected taxon for each gene
  gene_ageDf <- as.data.frame(tapply(mdDataExtended$cat,
                                     mdDataExtended$geneID,
                                     min))

  setDT(gene_ageDf, keep.rownames = TRUE)[]
  setnames(gene_ageDf, 1:2, c("geneID", "cat"))  # rename columns
  row.names(gene_ageDf) <- NULL   # remove row names

  ### convert cat into gene_age
  gene_ageDf$age[gene_ageDf$cat == "0000001"] <- "07_LUCA"
  gene_ageDf$age[gene_ageDf$cat == "0000011" | gene_ageDf$cat == "0000010"] <-
    paste0(
      "06_",
      as.character(
        taxa_list$fullName[taxa_list$ncbiID == sub_first_line$superkingdom
                           & taxa_list$rank == "superkingdom"]
      )
    )
  
  gene_ageDf$age[gene_ageDf$cat == "0000111"] <-
    paste0(
      "05_",
      as.character(
        taxa_list$fullName[taxa_list$ncbiID == sub_first_line$kingdom
                           & taxa_list$rank == "kingdom"]
      )
    )

  gene_ageDf$age[gene_ageDf$cat == "0001111"] <-
    paste0(
      "04_",
      as.character(
        taxa_list$fullName[taxa_list$ncbiID == sub_first_line$phylum
                           & taxa_list$rank == "phylum"]
      )
    )
  
  gene_ageDf$age[gene_ageDf$cat == "0011111"] <-
    paste0(
      "03_",
      as.character(
        taxa_list$fullName[taxa_list$ncbiID == sub_first_line$class
                           & taxa_list$rank == "class"]
      )
    )
  
  gene_ageDf$age[gene_ageDf$cat == "0111111"] <-
    paste0(
      "02_",
      as.character(
        taxa_list$fullName[taxa_list$ncbiID == sub_first_line$family
                           & taxa_list$rank == "family"]
      )
    )

  gene_ageDf$age[gene_ageDf$cat == "1111111"] <-
    paste0(
      "01_",
      as.character(
        taxa_list$fullName[taxa_list$fullName == in_select
                           & taxa_list$rank == rankName]
      )
    )

  # return gene_age data frame
  gene_ageDf <- gene_ageDf[, c("geneID", "cat", "age")]
  gene_ageDf$age[is.na(gene_ageDf$age)] <- "Undef"

  return(gene_ageDf)
}

#' create data for gene age plot (count number of genes in each estimated age)
#' @param gene_ageDf data of estimated gene age (from fn "estimate_gene_age")

gene_age_plotDf <- function(gene_ageDf){
  countDf <- plyr::count(gene_ageDf, c("age"))
  countDf$percentage <- round(countDf$freq / sum(countDf$freq) * 100)
  countDf$pos <- cumsum(countDf$percentage) - (0.5 * countDf$percentage)
  return(countDf)
}

#' get list of gene for a selected age
#' @param gene_ageDf data of estimated gene age (from fn "estimate_gene_age")
#' @param clicked_x x coordinate of selected age

get_selected_gene_age <- function(gene_ageDf, clicked_x){
  data <- gene_ageDf

  # calculate the coordinate range for each age group
  rangeDf <- plyr::count(data, c("age"))

  rangeDf$percentage <- round(rangeDf$freq / sum(rangeDf$freq) * 100)
  rangeDf$rangeStart[1] <- 0
  rangeDf$rangeEnd[1] <- rangeDf$percentage[1]
  if (nrow(rangeDf) > 1) {
    for (i in 2:nrow(rangeDf)) {
      rangeDf$rangeStart[i] <- rangeDf$rangeEnd[i - 1] + 1
      rangeDf$rangeEnd[i] <- rangeDf$percentage[i] + rangeDf$rangeEnd[i - 1]
    }
  }

  # get list of gene for selected age
  if (is.null(clicked_x)) return()
  else{
    corX <- 100 - round(-clicked_x)
    selectAge <- {
      as.character(rangeDf[rangeDf$rangeStart <= corX
                           & rangeDf$rangeEnd >= corX, ]$age)
    }
    subData <- subset(data, age == selectAge)
    data <- data[data$age == selectAge, ]
  }

  # return list of genes
  geneList <- levels(as.factor(subData$geneID))
  return(geneList)
}

#' create gene age plot
#' @param count_df data for gene age plot (from fn "gene_age_plotDf")
#' @param gene_age_text text size (from input$gene_age_text)

create_gene_age_plot <- function(count_df, gene_age_text){
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
  return(p)
}
