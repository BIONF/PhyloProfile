# NOT WORK WITH shinyapp.io ===================================================
# if(!("pacman" %in% installed.packages())) install.packages("pacman")
# library(pacman)
# p_load(shiny,
#        shinyBS,
#        ggplot2,
#        reshape2,
#        plyr,
#        dplyr,
#        tidyr,
#        scales,
#        grid,
#        gridExtra,
#        ape,
#        stringr,
#        gtable,
#        dendextend,
#        ggdendro,
#        gplots,
#        data.table,
#        taxize,
#        install = T)

 # packages <- c("shiny",
 #               "shinyBS",
 #               "ggplot2",
 #               "reshape2",
 #               "plyr",
 #               "dplyr",
 #               "tidyr",
 #               "scales",
 #               "grid",
 #               "gridExtra",
 #               "ape",
 #               "stringr",
 #               "gtable",
 #               "dendextend",
 #               "ggdendro",
 #               "gplots",
 #               "data.table",
 #               "taxize",
 #               "Biostrings",
 #               "zoo",
 #               "RCurl")
# sapply(packages, require, character.only = TRUE)

if (!require("shiny")) install.packages("shiny")
if (!require("shinyBS")) install.packages("shinyBS")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape2")) install.packages("reshape2")
if (!require("plyr")) install.packages("plyr")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("scales")) install.packages("scales")
if (!require("grid")) install.packages("grid")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("ape")) install.packages("ape")
if (!require("stringr")) install.packages("stringr")
if (!require("gtable")) install.packages("gtable")
if (!require("dendextend")) install.packages("dendextend")
if (!require("ggdendro")) install.packages("ggdendro")
if (!require("gplots")) install.packages("gplots")
if (!require("data.table")) install.packages("data.table")
if (!require("Biostrings")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("Biostrings")
}
if (!require("taxize")) install.packages("taxize")
if (!require("zoo")) install.packages("zoo")
if (!require("RCurl")) install.packages("RCurl")
if (!require("shinycssloaders")) {
  if ("devtools" %in% installed.packages() == FALSE){
    install.packages("devtools")
  }
  devtools::install_github("andrewsali/shinycssloaders", force = TRUE)
}
if (!require("Matching")) install.packages("Matching")
source("scripts/taxonomyProcessing.R")
source("scripts/functions.R")
source("scripts/get_oma_browser.R")

# MAIN ========================================================================
options(shiny.maxRequestSize = 99 * 1024 ^ 2)  # size limit for input 99mb

shinyServer(function(input, output, session) {
  # Automatically stop a Shiny app when closing the browser tab ---------------
  # session$onSessionEnded(stopApp)
  session$allowReconnect(TRUE)

  # ===========================================================================
  # PRE-PROCESSING  ===========================================================
  # ===========================================================================

  # check for internet connection ---------------------------------------------
  observe({
    if (has_internet() == FALSE){
      # toggleState("demo")
      toggleState("demo_data")
    }
  })

  # MISSING COMMENT -----------------------------------------------------------
  output$no_internet_msg <- renderUI({
    if (has_internet() == FALSE){
      strong(em("Internet connection is required for using demo data!"),
             style = "color:red")
    } else {
      return()
    }
  })

  # check for the existence of taxonomy info file -----------------------------
  observe({
    if (!file.exists(isolate("data/rankList.txt"))){
      if (has_internet() == TRUE){
        ncol <- max(count.fields("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/rankList.txt",
                                 sep = "\t"))
        df <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/rankList.txt",
                         sep = "\t",
                         quote = "",
                         header = F,
                         fill = T,
                         na.strings = c("", "NA"),
                         col.names = paste0("V", seq_len(ncol)))
        write.table(df, file = "data/rankList.txt",
                    col.names = F,
                    row.names = F,
                    quote = F,
                    sep = "\t") #na = "",
      } else {
        file.create("data/rankList.txt")
      }
    }
  })

  observe({
    if (!file.exists(isolate({
      "data/idList.txt"
      }
      ))){
      if (has_internet() == TRUE){
        ncol <- max(count.fields("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/idList.txt",
                                 comment.char = "",
                                 sep = "\t"))
        df <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/idList.txt",
                         sep = "\t",
                         header = F,
                         fill = T,
                         comment.char = "",
                         na.strings = c("", "NA"),
                         col.names = paste0("V", seq_len(ncol)))
        write.table(df, file = "data/idList.txt",
                    col.names = F,
                    row.names = F,
                    quote = F,
                    sep = "\t") #na = "",
      } else {
        file.create("data/idList.txt")
      }
    }
  })

  observe({
    if (!file.exists(isolate({
      "data/taxonNamesReduced.txt"
      }
      ))){
      if (has_internet() == TRUE){
        ncol <- max(count.fields("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/taxonNamesReduced.txt",
                                 sep = "\t"))
        df <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/taxonNamesReduced.txt",
                         sep = "\t",
                         quote = "",
                         header = F,
                         fill = T,
                         na.strings = c("", "NA"),
                         col.names = paste0("V", seq_len(ncol)))
        write.table(df, file = "data/taxonNamesReduced.txt",
                    col.names = F,
                    row.names = F,
                    quote = F,
                    sep = "\t")
      } else {
        system("cp data/newTaxa.txt data/taxonNamesReduced.txt")
      }
    }
  })

  observe({
    if (!file.exists(isolate({
      "data/taxonomyMatrix.txt"
      }
      ))){
      if (has_internet() == TRUE){
        ncol <- max(count.fields("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/taxonomyMatrix.txt",
                                 sep = "\t"))
        df <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/taxonomyMatrix.txt",
                         sep = "\t",
                         quote = "",
                         header = F,
                         fill = T,
                         na.strings = c("", "NA"),
                         col.names = paste0("V", seq_len(ncol)))

        write.table(df, file = "data/taxonomyMatrix.txt",
                    col.names = F,
                    row.names = F,
                    quote = F,
                    sep = "\t")
      }
    }
  })

  # output mismatched taxa ----------------------------------------------------
  output$not_found_taxa <- renderDataTable(option = list(searching = FALSE), {
    if (input$id_search > 0){
      if (length(taxa_id()) > 0){
        tb <- as.data.frame(taxa_id())
        tb_filtered <- tb[tb$type == "notfound", ]
        not_found_dt <- tb_filtered[, c("name", "new_name", "id")]
        colnames(not_found_dt) <- c("Summitted name",
                                  "Alternative name",
                                  "Alternative ID")
        not_found_dt
      }
    }
  })

  # MISSING COMMENT -----------------------------------------------------------
  output$download_not_found_taxa <- downloadHandler(
    filename = function(){
      c("mismatchedTaxa.txt")
      },
    content = function(file){
      tb <- as.data.frame(taxa_id())
      tb_filtered <- tb[tb$type == "notfound", ]
      not_found_dt <- tb_filtered[, c("name", "new_name", "id")]
      colnames(not_found_dt) <- c("Summitted name",
                                "Alternative name",
                                "Alternative ID")

      write.table(not_found_dt, file,
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)
    }
  )

  # output retrieved taxa IDs -------------------------------------------------
  output$taxa_id <- renderDataTable(option = list(searching = FALSE), {
    if (input$id_search > 0){
      if (length(taxa_id()) > 0){
        tb <- as.data.frame(taxa_id())
        tb_filtered <- tb[tb$type == "retrieved", ]
        retrieved_dt <- tb_filtered[, c("name", "id")]
        colnames(retrieved_dt) <- c("Taxon_name", "Taxon_ID")
        retrieved_dt
      }
    }
  })

  output$download_taxa_id <- downloadHandler(
    filename = function(){
      c("retrievedtaxa_id.txt")
      },
    content = function(file){
      tb <- as.data.frame(taxa_id())
      tb_filtered <- tb[tb$type == "retrieved", ]
      retrieved_dt <- tb_filtered[, c("name", "id")]
      colnames(retrieved_dt) <- c("Taxon name", "Taxon ID")

      write.table(retrieved_dt, file,
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)
    }
  )

  # PARSING VARIABLE 1 AND 2 ==================================================
  # render textinput for variable 1 & 2 ---------------------------------------
  output$var1_id.ui <- renderUI({
    long_dataframe <- get_long_matrix()

    if (is.null(long_dataframe)){
      textInput("var1_id",
                h5("1st variable:"),
                value = "Variable 1",
                width = "100%",
                placeholder = "Name of first variable")
    } else{
      textInput("var1_id", h5("1st variable:"),
                value = colnames(long_dataframe)[4],
                width = "100%",
                placeholder = "Name of first variable")
    }
  })

  output$var2_id.ui <- renderUI({
    long_dataframe <- get_long_matrix()
    if (is.null(long_dataframe)){
      textInput("var2_id",
                h5("1st variable:"),
                value = "Variable 2",
                width = "100%",
                placeholder = "Name of first variable")
    } else{
      textInput("var2_id", h5("1st variable:"),
                value = colnames(long_dataframe)[5],
                width = "100%",
                placeholder = "Name of first variable")
    }
  })

  # variable 1 & 2 cutoff slidebar (main plot) --------------------------------
  output$var1_cutoff.ui <- renderUI({
    if (is.null(input$var1_id)) return()
    if (input$var1_id == ""){
      sliderInput("var1", paste(input$var1_id, "cutoff:"),
                  min = 1,
                  max = 1,
                  step = 0.025,
                  value = c(1.0, 1.0),
                  width = 200)
    } else {
      sliderInput("var1", paste(input$var1_id, "cutoff:"),
                  min = 0,
                  max = 1,
                  step = 0.025,
                  value = c(0.0, 1.0),
                  width = 200)
    }
  })

  output$var2_cutoff.ui <- renderUI({
    if (is.null(input$var2_id)) return()
    if (input$var2_id == ""){
      sliderInput("var2", paste(input$var2_id, "cutoff:"),
                  min = 1,
                  max = 1,
                  step = 0.025,
                  value = c(1.0, 1.0),
                  width = 200)
    } else {
      sliderInput("var2", paste(input$var2_id, "cutoff:"),
                  min = 0,
                  max = 1,
                  step = 0.025,
                  value = c(0.0, 1.0),
                  width = 200)
    }
  })

  # render filter slidebars for Customized plot -------------------------------
  output$var1_filter.ui <- renderUI({
    if (is.null(input$var1_id)) return()
    if (input$var1_id == ""){
      sliderInput("var1cus", paste(input$var1_id, "cutoff:"),
                  min = 1,
                  max = 1,
                  step = 0.025,
                  value = c(1.0, 1.0),
                  width = 200)
    } else {
      sliderInput("var1cus", paste(input$var1_id, "cutoff:"),
                  min = 0,
                  max = 1,
                  step = 0.025,
                  value = c(input$var1[1], input$var1[2]),
                  width = 200)
    }
  })

  output$var2_filter.ui <- renderUI({
    if (is.null(input$var2_id)) return()
    if (input$var2_id == ""){
      sliderInput("var2cus", paste(input$var2_id, "cutoff:"),
                  min = 1,
                  max = 1,
                  step = 0.025,
                  value = c(1.0, 1.0),
                  width = 200)
    } else {
      sliderInput("var2cus", paste(input$var2_id, "cutoff:"),
                  min = 0,
                  max = 1,
                  step = 0.025,
                  value = c(input$var2[1], input$var2[2]),
                  width = 200)
    }
  })

  output$percent_filter.ui <- renderUI({
    sliderInput("percent2",
                "% of present taxa:",
                min = 0,
                max = 1,
                step = 0.025,
                value = input$percent,
                width = 200)
  })

  # render filter slidebars for Distribution plot -----------------------------
  output$var1_dist.ui <- renderUI({
    if (is.null(input$var1_id)) return()
    if (input$var1_id == ""){
      sliderInput("var1_dist",
                  paste(input$var1_id, "cutoff:"),
                  min = 1,
                  max = 1,
                  step = 0.025,
                  value = c(1.0, 1.0),
                  width = 200)
    } else {
      sliderInput("var1_dist", paste(input$var1_id, "cutoff:"),
                  min = 0,
                  max = 1,
                  step = 0.025,
                  value = c(input$var1[1], input$var1[2]),
                  width = 200)
    }
  })

  output$var2_dist.ui <- renderUI({
    if (is.null(input$var2_id)) return()
    if (input$var2_id == ""){
      sliderInput("var2_dist", paste(input$var2_id, "cutoff:"),
                  min = 1,
                  max = 1,
                  step = 0.025,
                  value = c(1.0, 1.0),
                  width = 200)
    } else {
      sliderInput("var2_dist", paste(input$var2_id, "cutoff:"),
                  min = 0,
                  max = 1,
                  step = 0.025,
                  value = c(input$var2[1], input$var2[2]),
                  width = 200)
    }
  })

  output$percent_dist.ui <- renderUI({
    sliderInput("percent_dist",
                "% of present taxa:",
                min = 0,
                max = 1,
                step = 0.025,
                value = input$percent,
                width = 200)
  })

  # render filter slidebars for Gene age estimation plot ----------------------
  output$var1_age.ui <- renderUI({
    if (is.null(input$var1_id)) return()
    if (input$var1_id == ""){
      sliderInput("var1_age", paste(input$var1_id, "cutoff:"),
                  min = 1,
                  max = 1,
                  step = 0.025,
                  value = c(1.0, 1.0),
                  width = 200)
    } else {
      sliderInput("var1_age",
                  paste(input$var1_id, "cutoff:"),
                  min = 0,
                  max = 1,
                  step = 0.025,
                  value = c(input$var1[1], input$var1[2]),
                  width = 200)
    }
  })

  output$var2_age.ui <- renderUI({
    if (is.null(input$var2_id)) return()
    if (input$var2_id == ""){
      sliderInput("var2_age", paste(input$var2_id, "cutoff:"),
                  min = 1,
                  max = 1,
                  step = 0.025,
                  value = c(1.0, 1.0),
                  width = 200)
    } else {
      sliderInput("var2_age", paste(input$var2_id, "cutoff:"),
                  min = 0,
                  max = 1,
                  step = 0.025,
                  value = c(input$var2[1], input$var2[2]),
                  width = 200)
    }
  })

  output$percent_age.ui <- renderUI({
    sliderInput("percent_age", "% of present taxa:",
                min = 0,
                max = 1,
                step = 0.025,
                value = input$percent,
                width = 200)
  })

  # render filter slidebars for Core gene finding function --------------------
  output$var1_cons.ui <- renderUI({
    if (is.null(input$var1_id)) return()
    if (input$var1_id == ""){
      sliderInput("var1_cons", paste(input$var1_id, "cutoff:"),
                  min = 1,
                  max = 1,
                  step = 0.025,
                  value = c(1.0, 1.0),
                  width = 200)
    } else {
      sliderInput("var1_cons",
                  paste(input$var1_id, "cutoff:"),
                  min = 0,
                  max = 1,
                  step = 0.025,
                  value = c(input$var1[1], input$var1[2]),
                  width = 200)
    }
  })

  output$var2_cons.ui <- renderUI({
    if (is.null(input$var2_id)) return()
    if (input$var2_id == ""){
      sliderInput("var2_cons",
                  paste(input$var2_id, "cutoff:"),
                  min = 1,
                  max = 1,
                  step = 0.025,
                  value = c(1.0, 1.0),
                  width = 200)
    } else {
      sliderInput("var2_cons", paste(input$var2_id, "cutoff:"),
                  min = 0,
                  max = 1,
                  step = 0.025,
                  value = c(input$var2[1], input$var2[2]),
                  width = 200)
    }
  })

  output$percent_cons.ui <- renderUI({
    sliderInput("percent_cons", "% of present taxa:",
                min = 0,
                max = 1,
                step = 0.025,
                value = 0.5,
                width = 200)
  })

  # update value for "main" filter slidebars ----------------------------------
  # based on "Customized", "Distribution", "Gene age estimation" slidebars
  observe({
    new_var1 <- input$var1cus

    if (is.null(input$var1_id)) return()
    if (input$var1_id == ""){
      updateSliderInput(session, "var1",
                        value = new_var1,
                        min = 1,
                        max = 1,
                        step = 0.025)
    } else {
      updateSliderInput(session,
                        "var1",
                        value = new_var1,
                        min = 0,
                        max = 1,
                        step = 0.025)
    }
  })
  observe({
    new_var2 <- input$var2cus

    if (is.null(input$var2_id)) return()
    if (input$var2_id == ""){
      updateSliderInput(session, "var2",
                        value = new_var2,
                        min = 1,
                        max = 1,
                        step = 0.025)
    } else {
      updateSliderInput(session, "var2",
                        value = new_var2,
                        min = 0,
                        max = 1,
                        step = 0.025)
    }
  })
  observe({
    new_percent <- input$percent2
    updateSliderInput(session, "percent",
                      value = new_percent,
                      min = 0,
                      max = 1,
                      step = 0.025)
  })

  observe({
    new_var1 <- input$var1_dist

    if (is.null(input$var1_id)) return()
    if (input$var1_id == ""){
      updateSliderInput(session, "var1",
                        value = new_var1,
                        min = 1,
                        max = 1,
                        step = 0.025)
    } else {
      updateSliderInput(session, "var1",
                        value = new_var1,
                        min = 0,
                        max = 1,
                        step = 0.025)
    }
  })
  observe({
    new_var2 <- input$var2_dist

    if (is.null(input$var2_id)) return()
    if (input$var2_id == ""){
      updateSliderInput(session, "var2",
                        value = new_var2,
                        min = 1,
                        max = 1,
                        step = 0.025)
    } else {
      updateSliderInput(session, "var2",
                        value = new_var2,
                        min = 0,
                        max = 1,
                        step = 0.025)
    }
  })
  observe({
    new_percent <- input$percent_dist
    updateSliderInput(session, "percent",
                      value = new_percent,
                      min = 0,
                      max = 1,
                      step = 0.025)
  })

  observe({
    new_var1 <- input$var1_age

    if (is.null(input$var1_id)) return()
    if (input$var1_id == ""){
      updateSliderInput(session, "var1",
                        value = new_var1,
                        min = 1,
                        max = 1,
                        step = 0.025)
    } else {
      updateSliderInput(session, "var1",
                        value = new_var1,
                        min = 0,
                        max = 1,
                        step = 0.025)
    }
  })
  observe({
    new_var2 <- input$var2_age

    if (is.null(input$var2_id)) return()
    if (input$var2_id == ""){
      updateSliderInput(session, "var2",
                        value = new_var2,
                        min = 1,
                        max = 1,
                        step = 0.025)
    } else {
      updateSliderInput(session, "var2",
                        value = new_var2,
                        min = 0,
                        max = 1,
                        step = 0.025)
    }
  })
  observe({
    new_percent <- input$percent_age
    updateSliderInput(session, "percent",
                      value = new_percent,
                      min = 0,
                      max = 1,
                      step = 0.025)
  })

  # render 2. variable relationship according to demo data --------------------
  output$var2_relation.ui <- renderUI({
    if (input$demo_data == "ampk-tor"){
      selectInput("var2_relation", label = h5("Relationship:"),
                  choices = list("Prot-Prot" = "protein",
                                 "Prot-Spec" = "species"),
                  selected = "protein",
                  width = 130)
    } else {
      selectInput("var2_relation", label = h5("Relationship:"),
                  choices = list("Prot-Prot" = "protein",
                                 "Prot-Spec" = "species"),
                  selected = "species",
                  width = 130)
    }
  })

  # ===========================================================================

  # check oneseq fasta file exists --------------------------------------------
  output$one_seq.exist_check <- renderUI({
    #f <- toString(input$oneseq.file)

    if (is.null(input$one_seq_fasta)) return()
    else{
      f <- input$one_seq_fasta$datapath
      if (!file.exists(f)){
        helpText("File not exists!!")
      } else {
        if (length(readLines(f, n = 1)) == 0){
          helpText("is not a fasta file!!")
        } else {
          first_line <- readLines(f, n = 1)
          a <- substr(first_line, 1, 1)
          if (a == ">"){
            HTML('<p><span style="color: #0000ff;"><strong>Please click CLOSE to comfirm!</strong></span></p>')
          } else {
            helpText("is not a fasta file!!")
          }
        }
      }
    }
  })

  # reset cutoffs of main plot ------------------------------------------------
  observeEvent(input$reset_main, {
    shinyjs::reset("var1")
    shinyjs::reset("var2")
    shinyjs::reset("percent")
  })

  # reset config of main plot -------------------------------------------------
  observeEvent(input$reset_main_config, {
    shinyjs::reset("x_size")
    shinyjs::reset("y_size")
    shinyjs::reset("legend_size")
    shinyjs::reset("x_angle")
    shinyjs::reset("dot_zoom")
  })

  # close main config ---------------------------------------------------------
  observeEvent(input$apply_main_config, {
    toggleModal(session, "main_plot_config_bs", toggle = "close")
  })

  # reset cutoffs of Customized plot ------------------------------------------
  observeEvent(input$reset_selected, {
    shinyjs::reset("var1")
    shinyjs::reset("var2")
    shinyjs::reset("percent")
  })

  # reset config of customized plot -------------------------------------------
  observeEvent(input$reset_selected_config, {
    shinyjs::reset("x_size_select")
    shinyjs::reset("y_size_select")
    shinyjs::reset("legend_size_select")
    shinyjs::reset("x_angle_select")
    shinyjs::reset("dot_zoom_select")
  })

  # close customized config ---------------------------------------------------
  observeEvent(input$apply_selected_config, {
    toggleModal(session, "selected_plot_config_bs", toggle = "close")
  })

  # reset colors --------------------------------------------------------------
  observeEvent(input$default_color_var2, {
    shinyjs::reset("low_color_var2")
    shinyjs::reset("high_color_var2")
  })

  observeEvent(input$default_color_var1, {
    shinyjs::reset("low_color_var1")
    shinyjs::reset("high_color_var1")
  })

  observeEvent(input$default_color_para, {
    shinyjs::reset("para_color")
  })

  # ===========================================================================

  # render checkNewick.ui -----------------------------------------------------
  output$checkNewick.ui <- renderUI({
    filein <- input$inputTree
    if (is.null(filein)) return()

    check_newick <- check_newick(filein)
    if (check_newick == 1){
      updateButton(session, "do", disabled = TRUE)
      HTML("<p><em><span style=\"color: #ff0000;\"><strong>ERROR: Parenthesis(-es) missing!</strong></span></em></p>")
    } else if (check_newick == 2){
      updateButton(session, "do", disabled = TRUE)
      HTML("<p><em><span style=\"color: #ff0000;\"><strong>ERROR: Comma(s) missing!</strong></span></em></p>")
    } else if (check_newick == 3){
      updateButton(session, "do", disabled = TRUE)
      HTML("<p><em><span style=\"color: #ff0000;\"><strong>ERROR: Tree contains singleton!</strong></span></em></p>")
    } else if (check_newick == 0){
      return()
    } else {
      updateButton(session, "do", disabled = TRUE)
      strong(em(paste0(check_newick, " not exist in main input file!")),
             style = "color:red")
    }
  })

  # render input_check.ui -----------------------------------------------------
  output$input_check.ui <- renderUI({
    if(is.null(input$main_input)) return()
    input_type <- get_input_type()
    if (input_type == "noGeneID"){
      updateButton(session, "do", disabled = TRUE)
      HTML("<font color=\"red\"><em><strong>ERROR: Unsupported input format. <a href=\"https://github.com/BIONF/PhyloProfile/wiki/Input-Data\" target=\"_blank\">Click here for more info</a></em></strong></font>")
    } else if (input_type == "emptyCell"){
      updateButton(session, "do", disabled = TRUE)
      em(strong("ERROR: Rows have unequal length",
                style = "color:red"))
    }
    else if (input_type == "moreCol"){
      updateButton(session, "do", disabled = TRUE)
      em(strong("ERROR: More columns than column names",
                style = "color:red"))
    } else {
      updateButton(session, "do", disabled = FALSE)
      return()
    }
  })
  
  
  output$select_oma_type <- renderUI({
    if(is.null(input$main_input)) return()
    input_type <- get_input_type()
    if (input_type == "oma"){
      # Options to select the OMA type to generate the output
      selectInput("selected_oma_type", label = "Select type of orthologs:",
                  choices = list("PAIR", "HOG", "OG"),
                  selected = "PAIR")
    } else {
      return()
    }

  })
  
  output$button_oma <- renderUI({
    if(is.null(input$main_input)) return()
    input_type <- get_input_type()
    if (input_type == "oma"){
      shinyBS::bsButton("get_data_oma", "Get data")
    }
  })
  
  output$oma_download <- renderUI({
    if(is.null(input$main_input)) return()
    input_type <- get_input_type()
    if (input_type == "oma"){
      downloadButton("download_files_oma", "Download")
    } 
  })
  
  
  output$download_files_oma <- downloadHandler(
    filenname <- function(){
      "oma_data_to_phyloprofile_input.zip"
      },
    content <- function(file){
      write.table(get_long_matrix(), "long.txt",
                  sep = "\t",
                  row.names = FALSE,
                  col.names = FALSE,
                  quote = FALSE)
      print("long")
      write.table(long_to_fasta(get_long_matrix()), "fasta.txt",
                  sep = "\t",
                  row.names = FALSE,
                  col.names = TRUE,
                  quote = FALSE)
      print("fasta")
      write.table(get_domain_information (), "domain.txt",
                  sep = "\t",
                  row.names = FALSE,
                  col.names = FALSE,
                  quote = FALSE)
      print("domains")
  
      zip(zipfile = file,
          files = c("long.txt", "domain.txt", "fasta.txt")) 
    },
    contentType = "application/zip"
  )

  # render download link for demo online files --------------------------------
  output$main_input_file.ui <- renderUI({
    # if(input$demo == TRUE){
    if (input$demo_data == "demo"){
      strong(a("Download demo input file",
               href = "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.long",
               target = "_blank"))
    } else if (input$demo_data == "ampk-tor"){
      strong(a("Download demo input file",
               href = "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile",
               target = "_blank"))
    } else {
      fileInput("main_input", h5("Upload input file:"))
    }
  })

  output$domain_input_file.ui <- renderUI({
    # if(input$demo == TRUE){
    if (input$demo_data == "demo"){
      strong(a("Download demo domain files",
               href = "https://github.com/BIONF/phyloprofile-data/tree/master/demo/domain_files",
               target = "_blank"))
    } else if (input$demo_data == "ampk-tor"){
      strong(a("Download demo domain file",
               href = "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.domains_F",
               target = "_blank"))
    } else {
      input_type <- get_input_type()
      if (input$anno_choose == "from file"){
        fileInput("file_domain_input", "")
      } else {
        textInput("domainPath", "", "")
      }
    }
  })

  output$download_fastaDemo.ui <- renderUI({
    if (input$demo_data == "demo"){
      strong(a("Download demo fasta file",
               href = "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/fasta_file/concatenatedSeq.fa",
               target = "_blank"))
    } else if (input$demo_data == "ampk-tor"){
      strong(a("Download demo fasta file",
               href = "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.extended.fa",
               target = "_blank"))
    }
  })

  # render description for demo data ------------------------------------------
  output$demo_data_describe <- renderUI({
    if (input$demo_data == "none"){
      return()
    } else if (input$demo_data == "ampk-tor"){
      em(a("Data description",
           href = "https://github.com/BIONF/phyloprofile-data/blob/master/expTestData/ampk-tor/README.md",
           target = "_blank"))
    } else {
      em(a("Data description",
           href = "https://github.com/BIONF/phyloprofile-data/blob/master/demo/README.md",
           target = "_blank"))
    }
  })

  # get status of unkTaxa for conditional panel -------------------------------
  output$unk_taxa_status <- reactive({
    unkTaxa <- unkTaxa()
    length(unkTaxa) > 0
  })

  outputOptions(output, "unk_taxa_status", suspendWhenHidden = FALSE)

  # show full list of unkTaxa -------------------------------------------------
  output$unk_taxa_full <- renderDataTable(option = list(searching = FALSE,
                                                        pageLength = 10), {
    if (length(unkTaxa()) > 0){
      tb <- as.data.frame(unkTaxa())
      names(tb)[1] <- "New taxon"
      tb
    }
  })

  # check if data is loaded and "parse" button is clicked and confirmed -------
  # (get info from input)
  v1 <- reactiveValues(parse = FALSE)
  observeEvent(input$but_parse, {
    toggleModal(session, "parse_confirm", toggle = "close")
    v1$parse <- input$but_parse
    updateButton(session, "but_parse", disabled = TRUE)
    toggleState("new_taxa_ask")
    toggleState("main_input")
  })

  # create rankList.txt, idList.txt, taxonNamesReduced.txt from input file ----
  # (if confirmed by BUTyes).
  # and also create a full taxonomy matrix for sorting taxa
  # (taxonomyMatrix.txt)
  invalidID <- {
    reactiveValues(df = data.frame("Invalid NCBI ID(s)" = as.character(),
                                   stringsAsFactors = F))
  }
  
  observe({
    filein <- input$main_input
    if(is.null(filein)) return()
    input_type <- get_input_type()
    if (input_type == "xml" |
        input_type == "long" |
        input_type == "wide" |
        input_type == "fasta" ){
        inputDf <- as.data.frame(read.table(file = filein$datapath,
                                            sep = "\t",
                                            header = T,
                                            check.names = FALSE,
                                            comment.char = ""))

        if (v1$parse == F) return()
        else{
          # get list of taxa need to be parsed
          #(the one mising taxonomy information)
          if (v1$parse == T){
            unkTaxa <- unkTaxa()
            titleline <- c("geneID", unkTaxa)
          }

          invalidIDtmp <- list()
          # Create 0-row data frame which will be used to store data
          withProgress(message = "Parsing input file", value = 0, {
            allTaxonInfo <- fread("data/taxonNamesFull.txt")
            newTaxaFromFile <- fread("data/newTaxa.txt",
                                     colClasses = c("ncbiID" = "character"))

            allTaxonInfo <- rbind(newTaxaFromFile, allTaxonInfo)

            rankList <- data.frame()
            idList <- data.frame()
            reducedInfoList <- data.frame()

            for (i in 2:length(titleline)){
              ## taxon ID
              refID <- gsub("ncbi", "", titleline[i])

              ## get info for this taxon
              refEntry <- allTaxonInfo[allTaxonInfo$ncbiID == refID, ]

              if (nrow(refEntry) < 1){
                invalidIDtmp <- c(invalidIDtmp, refID)
              } else {
                if (nrow(reducedInfoList[reducedInfoList$X1 == refEntry$ncbiID, ]) == 0){
                  refInfoList <- data.frame(matrix(c(refEntry$ncbiID,
                                                     refEntry$fullName,
                                                     refEntry$rank,
                                                     refEntry$parentID),
                                                   nrow = 1,
                                                   byrow = T),
                                            stringsAsFactors = FALSE)
                  reducedInfoList <- rbind(reducedInfoList, refInfoList)
                }

                # parentID (used to check if hitting last rank, i.e. norank_1)
                lastID <- refEntry$parentID

                # create list of rank for this taxon
                rank <- c(paste0("ncbi", refID), refEntry$fullName)
                if (refEntry$rank == "norank"){
                  rank <- c(rank, paste0("strain"))
                } else {
                  rank <- c(rank, refEntry$rank)
                }

                # create list of IDs for this taxon ---------------------------
                ids <- list(paste0(refEntry$fullName, "#name"))
                if (refEntry$rank == "norank"){
                  ids <- c(ids,
                           paste0(refEntry$ncbiID,
                                  "#",
                                  "strain",
                                  "_",
                                  refEntry$ncbiID))
                } else {
                  ids <- c(ids,
                           paste0(refEntry$ncbiID,
                                  "#",
                                  refEntry$rank))
                }

                # append info into rank and ids -------------------------------
                while (lastID != 1){
                  nextEntry <- allTaxonInfo[allTaxonInfo$ncbiID == lastID, ]

                  if (nrow(reducedInfoList[reducedInfoList$X1 == nextEntry$ncbiID, ]) == 0){
                    nextEntryList <- {
                      data.frame(matrix(c(nextEntry$ncbiID,
                                          nextEntry$fullName,
                                          nextEntry$rank,
                                          nextEntry$parentID),
                                          nrow = 1, byrow = T),
                                 stringsAsFactors = FALSE)
                    }
                    reducedInfoList <- rbind(reducedInfoList,
                                             nextEntryList)
                  }

                  lastID <- nextEntry$parentID

                  if ("norank" %in% nextEntry$rank){
                    rank <- c(rank,
                              paste0(nextEntry$rank,
                                     "_",
                                     nextEntry$ncbiID))
                    ids <- c(ids,
                             paste0(nextEntry$ncbiID,
                                    "#",
                                    nextEntry$rank,
                                    "_",
                                    nextEntry$ncbiID))
                  } else {
                    rank <- c(rank, nextEntry$rank)
                    ids <- c(ids,
                             paste0(nextEntry$ncbiID,
                                    "#",
                                    nextEntry$rank))
                  }
                }

                # last rank and id
                rank <- c(rank, "norank_1")
                ids <- c(ids, "1#norank_1")

                # append into rankList and idList files
                rankListTMP <- data.frame(matrix(unlist(rank),
                                                 nrow = 1, byrow = T),
                                          stringsAsFactors = FALSE)
                rankList <- rbind.fill(rankList, rankListTMP)
                idListTMP <- data.frame(matrix(unlist(ids),
                                               nrow = 1,
                                               byrow = T),
                                        stringsAsFactors = FALSE)
                idList <- rbind.fill(idList, idListTMP)
              }

              # Increment the progress bar, and update the detail text.
              incProgress(1 / (length(titleline) - 1),
                          detail = paste( (i - 1),
                                         "/",
                                         length(titleline) - 1))
            }
          })
          # save invalid IDs to invalidID$df
          invalidID$df <- as.data.frame(unlist(invalidIDtmp))

          if (nrow(invalidID$df) < 1){
            # open existing files (idList.txt, rankList.txt and taxonNamesReduced.txt)
            ncol <- max(count.fields("data/rankList.txt", sep = "\t"))
            # print(ncol)
            oldIDList <- {
              as.data.frame(read.table("data/idList.txt",
                                       sep = "\t",
                                       header = F,
                                       check.names = FALSE,
                                       comment.char = "",
                                       fill = T,
                                       stringsAsFactors = T,
                                       na.strings = c("", "NA"),
                                       col.names = paste0("X", seq_len(ncol))))
            }
            oldRankList <- {
              as.data.frame(read.table("data/rankList.txt",
                                       sep = "\t",
                                       header = F,
                                       check.names = FALSE,
                                       comment.char = "",
                                       fill = T,
                                       stringsAsFactors = T,
                                       na.strings = c("", "NA"),
                                       col.names = paste0("X", seq_len(ncol))))
            }
            oldNameList <- {
              as.data.frame(read.table("data/taxonNamesReduced.txt",
                                       sep = "\t",
                                       header = T,
                                       check.names = FALSE,
                                       comment.char = "",
                                       fill = T,
                                       stringsAsFactors = T))
            }

            # and append new info into those files
            new_idList <- rbind.fill(oldIDList, idList)
            new_rankList <- rbind.fill(oldRankList, rankList)
            colnames(reducedInfoList) <- c("ncbiID",
                                           "fullName",
                                           "rank",
                                           "parentID")
            new_nameList <- rbind.fill(oldNameList,
                                       reducedInfoList)

            write.table(new_idList[!duplicated(new_idList), ],
                        file  = "data/idList.txt",
                        col.names = F,
                        row.names = F,
                        quote = F,
                        sep = "\t")
            write.table(new_rankList[!duplicated(new_rankList), ],
                        file = "data/rankList.txt",
                        col.names = F,
                        row.names = F,
                        quote = F,
                        sep = "\t")
            write.table(new_nameList[!duplicated(new_nameList), ],
                        file = "data/taxonNamesReduced.txt",
                        col.names = T,
                        row.names = F,
                        quote = F,
                        sep = "\t")

            # create taxonomy matrix
            taxMatrix <- taxonomyTableCreator("data/idList.txt",
                                              "data/rankList.txt")
            write.table(taxMatrix,
                        "data/taxonomyMatrix.txt",
                        sep = "\t",
                        eol = "\n",
                        row.names = FALSE,
                        quote = FALSE)
          }
        }
      }
    
  })

  # output invalid NCBI ID ----------------------------------------------------
  output$invalidID.output <- renderTable({
    if (nrow(invalidID$df) < 1) return()
    else{
      outDf <- invalidID$df
      colnames(outDf) <- c("Invalid NCBI ID(s)")
      return(outDf)
    }
  })

  output$end_parsing_msg <- renderUI({
    if (nrow(invalidID$df) < 1) {
      strong(h4("PLEASE RELOAD THIS TOOL AFTER ADDING NEW TAXA!!!"),
             style = "color:red")
    }else{
      HTML('<p><strong><span style="color: #e12525;">SOME INVALID TAXON IDs HAVE BEEN FOUND!!</span><br>Please check the validity of the following IDs in
           <a target="_blank" href="https://www.ncbi.nlm.nih.gov/taxonomy">NCBI taxonomy database</a>!</strong></p>')
    }
    })

  # list of taxonomy ranks for plotting ---------------------------------------
  output$rank_select <- renderUI({
    # if(input$demo == TRUE){
    if (input$demo_data == "demo"){
      selectInput("rank_select", label = "",
                  choices = get_taxonomy_ranks(),
                  selected = "26_phylum")
    } else if (input$demo_data == "ampk-tor"){
      selectInput("rank_select", label = "",
                  choices = get_taxonomy_ranks(),
                  selected = "06_species")
    } else {
      selectInput("rank_select", label = "",
                  choices = get_taxonomy_ranks(),
                  selected = "06_species")
    }
  })

  # then output list of (super)taxa onto UI
  output$select <- renderUI({
    choice <- alltaxa_list()
    choice$fullName <- as.factor(choice$fullName)

    # if(input$demo == TRUE){
    if (input$demo_data == "demo"){
      hellemDf <- data.frame("name" = c("Encephalitozoon hellem",
                                        "Encephalitozoon hellem",
                                        "Encephalitozoon",
                                        "Unikaryonidae",
                                        "Apansporoblastina",
                                        "Apansporoblastina",
                                        "Microsporidia",
                                        "Fungi",
                                        "Eukaryota"),
                             "rank" = c("strain",
                                        "species",
                                        "genus",
                                        "family",
                                        "order",
                                        "class",
                                        "phylum",
                                        "kingdom",
                                        "superkingdom"))
      rank_select <- input$rank_select
      rankName <- substr(rank_select,
                        4,
                        nchar(rank_select))

      selectInput("in_select", "",
                  as.list(levels(choice$fullName)),
                  hellemDf$name[hellemDf$rank == rankName])
    } else if (input$demo_data == "ampk-tor"){
      humanDf <- data.frame("name" = c("Homo sapiens",
                                       "Homo sapiens",
                                       "Homo",
                                       "Hominidae",
                                       "Primates",
                                       "Mammalia",
                                       "Chordata",
                                       "Metazoa",
                                       "Eukaryota"),
                            "rank" = c("strain",
                                       "species",
                                       "genus",
                                       "family",
                                       "order",
                                       "class",
                                       "phylum",
                                       "kingdom",
                                       "superkingdom"))
      rank_select <- input$rank_select
      rankName <- substr(rank_select, 4, nchar(rank_select))

      selectInput("in_select", "",
                  as.list(levels(choice$fullName)),
                  humanDf$name[humanDf$rank == rankName])
    } else {
      selectInput("in_select", "",
                  as.list(levels(choice$fullName)),
                  levels(choice$fullName)[1])
    }
  })

  output$highlight_taxon_ui <- renderUI({
    choice <- alltaxa_list()
    choice$fullName <- as.factor(choice$fullName)

    out <- as.list(levels(choice$fullName))
    out <- append("none", out)

    selectInput("taxon_highlight", "Select (super)taxon to highlight:",
                out, selected = out[1])
  })

  # update highlight_taxon_ui based on double clicked dot ---------------------
  observe({
    choice <- alltaxa_list()
    choice$fullName <- as.factor(choice$fullName)

    out <- as.list(levels(choice$fullName))
    out <- append("none", out)

    if (!is.null(input$plot_dblclick)){
      if (input$x_axis == "genes"){
        corX <- round(input$plot_dblclick$y);
        corY <- round(input$plot_dblclick$x)
      } else {
        corX <- round(input$plot_dblclick$x);
        corY <- round(input$plot_dblclick$y)
      }

      dataHeat <- dataHeat()
      supertaxa <- levels(dataHeat$supertaxon)
      spec <- toString(supertaxa[corX])
      selectedIndex <- match(substr(spec, 6, nchar(spec)), out)

      updateSelectInput(session, "taxon_highlight",
                        label = "Select (super)taxon to highlight:",
                        choices = out,
                        selected = out[selectedIndex])
    }
  })

  # Get list of genes for highlighting ----------------------------------------
  output$highlight_gene_ui <- renderUI({
    geneList <- dataHeat()
    geneList$geneID <- as.factor(geneList$geneID)

    out <- as.list(levels(geneList$geneID))
    out <- append("none", out)

    selectInput("gene_highlight", "Highlight:", out, selected = out[1])
  })

  # update highlight_gene_ui based on double clicked dot ----------------------
  observe({
    if (!is.null(input$plot_dblclick)){
      geneList <- dataHeat()
      if (input$apply_cluster == TRUE){
        geneList <- clusteredDataHeat()
      }

      geneList$geneID <- as.factor(geneList$geneID)

      out <- as.list(levels(geneList$geneID))
      out <- append("none", out)

      clickedInfo <- mainpoint_info()

      if (input$x_axis == "genes"){

        corX <- round(input$plot_dblclick$y);
        corY <- round(input$plot_dblclick$x)
      } else {
        corX <- round(input$plot_dblclick$x);
        corY <- round(input$plot_dblclick$y)
      }
      updateSelectInput(session, "gene_highlight",
                        label = "Highlight:",
                        choices = out,
                        selected = out[corY + 1])
    } else return()
  })

  # print total number of genes -----------------------------------------------
  output$total_gene_number.ui <- renderUI({
    geneList <- preData()
    geneList$geneID <- as.factor(geneList$geneID)
    out <- as.list(levels(geneList$geneID))
    if (length(out) > 0){
      # em(paste0("Total number of genes:  ",length(out)))
      strong(paste0("Total number of genes:  ", length(out)))
    }
  })

  # enable "PLOT" button ------------------------------------------------------
  observeEvent(input$rank_select,  ({
    if (input$rank_select == ""){
      updateButton(session, "do", disabled = TRUE)
    } else{
      unkTaxa <- unkTaxa()
      if (length(unkTaxa) == 0){
        updateButton(session, "do", disabled = FALSE)
      }
    }
  }))

  # move to main tab when "PLOT" button has been clicked ----------------------
  observe({
    # use tabsetPanel "id" argument to change tabs
    if (input$do > 0) {
      updateTabsetPanel(session, "tabs", selected = "Main profile")
    }
  })

  # disable main input, genelist input and initial questions ------------------
  observe({
    # use tabsetPanel "id" argument to change tabs
    if (input$do > 0) {
      toggleState("main_input")
      toggleState("gene_list_selected")
      toggleState("demo_data")
    }
  })

  # disable demo checkbox and update var2_aggregate_by to mean ----------------
  # if using demo data
  observe({
    # if (input$demo == TRUE) {
    if (input$demo_data == "demo") {
      ### disable demo checkbox
      # toggleState("demo")
      ### update var2_aggregate_by to mean
      updateSelectInput(session, "var2_aggregate_by",
                        choices = list("Max" = "max",
                                       "Min" = "min",
                                       "Mean" = "mean",
                                       "Median" = "median"),
                        selected = "mean")
    }
  })

  # ===========================================================================
  # ADD NEW TAXA ==============================================================
  # ===========================================================================
  newTaxa <- reactiveValues()
  newTaxa$Df <- data.frame("ncbiID" = numeric(),
                           "fullName" = character(),
                           "rank" = character(),
                           "parentID" = numeric(),
                           stringsAsFactors = FALSE)
  newIndex <- reactiveValues()
  newIndex$value <- 1

  observeEvent(input$new_add, {
    newTaxa$Df[newIndex$value, ] <- c(input$new_id,
                                     input$new_name,
                                     input$new_rank,
                                     input$new_parent)
    newIndex$value <- newIndex$value + 1
    updateTextInput(session, "new_id",
                    value = as.numeric(input$new_id) + 1)
    updateTextInput(session, "new_name", value = "")
    updateTextInput(session, "new_rank", value = "norank")
    updateTextInput(session, "new_parent", value = "")
  })

  observeEvent(input$new_done, {
    toggleModal(session, "add_taxa_windows", toggle = "close")
    write.table(newTaxa$Df, "data/newTaxa.txt",
                sep = "\t",
                eol = "\n",
                row.names = FALSE,
                quote = FALSE)
  })

  # ===========================================================================
  # PROCESSING INPUT DATA =====================================================
  # ===========================================================================

  # check if data is loaded and "plot" button is clicked ----------------------
  v <- reactiveValues(doPlot = FALSE)
  observeEvent(input$do, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    v$doPlot <- input$do
    filein <- input$main_input
    # if(is.null(filein) & input$demo==FALSE){
    if (is.null(filein) & input$demo_data == "none"){
      v$doPlot <- FALSE
      updateButton(session, "do", disabled = TRUE)
    }
  })

  # PARSING DATA FROM INPUT MATRIX: ===========================================
  # get (super)taxa names (3)
  # calculate percentage of presence (4),
  # max/min/mean/median VAR1 (5) and VAR2 (6)
  # if group input taxa list into higher taxonomy rank

  # check if "no ordering gene IDs" has been checked --------------------------
  output$apply_cluster_check.ui <- renderUI({
    if (input$ordering == FALSE){
      HTML('<p><em>(Check "Ordering sequence IDs" check box in <strong>Input & settings tab</strong>&nbsp;to enable this function)</em></p>')
    }
  })

  observe({
    if (input$ordering == FALSE){
      shinyjs::disable("apply_cluster")
    } else {
      shinyjs::enable("apply_cluster")
    }
  })

  # ===========================================================================
  # DATA & PLOT FOR MAIN PROFILE ==============================================
  # ===========================================================================

  # get list of all sequence IDs for selectize input (customized profile) -----
  output$gene_in <- renderUI({
    filein <- input$main_input
    fileCustom <- input$custom_file

    # if(input$demo == TRUE){
    if (input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      filein <- 1
    }

    if (is.null(filein) & is.null(fileCustom)){
      return(selectInput("in_seq", "", "all"))
      }
    if (v$doPlot == FALSE){
      return(selectInput("in_seq", "", "all"))
      }
    else{
      # full list
      data <- as.data.frame(get_data_filtered())
      data$geneID <- as.character(data$geneID)
      data$geneID <- as.factor(data$geneID)
      outAll <- as.list(levels(data$geneID))
      outAll <- append("all", outAll)
      #selectInput("in_seq","",out,selected=out[1],multiple=TRUE)

      if (input$add_custom_profile == TRUE){
        out <- selectedgene_age()
        if (length(out) > 0){
          selectInput("in_seq",
                      "",
                      out,
                      selected = as.list(out),
                      multiple = TRUE,
                      selectize = FALSE)
        }
        else {
          selectInput("in_seq",
                      "",
                      outAll,
                      selected = outAll[1],
                      multiple = TRUE,
                      selectize = FALSE)
        }
      } else if (input$add_cluster_cutom_profile == TRUE){
        out <- brushed_clusterGene()
        if (length(out) > 0){
          selectInput("in_seq", "",
                      out,
                      selected = as.list(out),
                      multiple = TRUE,
                      selectize = FALSE)
        }
        else {
          selectInput("in_seq", "",
                      outAll,
                      selected = outAll[1],
                      multiple = TRUE,
                      selectize = FALSE)
        }
      }
      else if (input$add_cons_gene_custom_profile == TRUE){
        out <- cons_geneDf()
        if (length(out) > 0){
          selectInput("in_seq", "",
                      out,
                      selected = as.list(out),
                      multiple = TRUE,
                      selectize = FALSE)
        }
        else {
          selectInput("in_seq", "",
                      outAll,
                      selected = outAll[1],
                      multiple = TRUE,
                      selectize = FALSE)
        }
      }
      # else if(input$add_gc_genes_custom_profile == TRUE){
      #   out <- significantGenesGroupCompairison$geneID
      #   if(length(out)>0){
      #     selectInput("in_seq","",out,selected=as.list(out),multiple=TRUE,selectize=FALSE)
      #   }
      #   else {
      #     selectInput("in_seq","",outAll,selected=outAll[1],multiple=TRUE,selectize=FALSE)
      #   }
      #
      # }
      else {
        if (is.null(fileCustom)){
          selectInput("in_seq", "",
                      outAll,
                      selected = outAll[1],
                      multiple = TRUE,
                      selectize = FALSE)
        }
        else {
          customList <- as.data.frame(read.table(file = fileCustom$datapath,
                                                 header = FALSE))
          customList$V1 <- as.factor(customList$V1)
          out <- as.list(levels(customList$V1))
          selectInput("in_seq", "",
                      out,
                      selected = out,
                      multiple = TRUE,
                      selectize = FALSE)
        }
      }
    }
  })

  # get list of all taxa for selectize input ----------------------------------
  output$taxa_in <- renderUI({
    filein <- input$main_input
    # if(input$demo == TRUE){
    if (input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      filein <- 1
    }

    if (is.null(filein)) return(selectInput("in_taxa", "", "all"))
    if (v$doPlot == FALSE) return(selectInput("in_taxa", "", "all"))
    else{
      choice <- alltaxa_list()
      choice$fullName <- as.factor(choice$fullName)

      out <- as.list(levels(choice$fullName))
      out <- append("all", out)

      if (input$apply_cus_taxa == TRUE){
        out <- cus_taxaName()
        selectInput("in_taxa", "",
                    out,
                    selected = out,
                    multiple = TRUE,
                    selectize = FALSE)
      } else {
        selectInput("in_taxa", "",
                    out,
                    selected = out[1],
                    multiple = TRUE,
                    selectize = FALSE)
      }
    }
  })

  # render dot size to dot_size_info ------------------------------------------
  output$dot_size_info <- renderUI({
    if (v$doPlot == FALSE) return()

    dataHeat <- dataHeat()
    dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
    presentVl <- dataHeat$presSpec[!is.na(dataHeat$presSpec)]

    minDot <- (floor(min(presentVl) * 10) / 10 * 5) * (1 + input$dot_zoom)
    maxDot <- (floor(max(presentVl) * 10) / 10 * 5) * (1 + input$dot_zoom)

    em(paste0("current point's size: ", minDot, " - ", maxDot))
  })

  # plot profile into plot.ui -------------------------------------------------
  output$mainPlot <- renderPlot({
    if (input$auto_update == FALSE){
      # Add dependency on the update button
      # (only update when button is clicked)
      input$update_btn

      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        main_plot(v, dataHeat(), clusteredDataHeat(), get_input_main ())
      })
    } else {
      main_plot(v, dataHeat(), clusteredDataHeat(), get_input_main ())
    }
  })

  output$plot.ui <- renderUI({
    # show beschreibung file if no plot present
    if (v$doPlot == FALSE){
      return()
    } else{
      # if auto_update is NOT selected, use update_btn to trigger plot changing
      if (input$auto_update == FALSE){
        # Add dependency on the update button
        # (only update when button is clicked)
        input$update_btn

        # Add all the filters to the data based on the user inputs
        # wrap in an isolate() so that the data won't update
        # every time an input is changed
        isolate({
          withSpinner(
            plotOutput("mainPlot",
                       width = input$width,
                       height = input$height,
                       click = "plot_click",
                       dblclick = "plot_dblclick",
                       hover = hoverOpts(
                         id = "plot_hover",
                         delay = input$hover_delay,
                         delayType = input$hover_policy,
                         nullOutside = input$hover_null_outside
                       )
            )
          )
        })
      }
      # if auto_update is true
      else {
        withSpinner(
          plotOutput("mainPlot",
                     width = input$width,
                     height = input$height,
                     click = "plot_click",
                     dblclick = "plot_dblclick",
                     hover = hoverOpts(
                       id = "plot_hover",
                       delay = input$hover_delay,
                       delayType = input$hover_policy,
                       nullOutside = input$hover_null_outside
                     )
          )
        )
      }
    }
  })

  # download main plot --------------------------------------------------------
  output$plot_download <- downloadHandler(
    filename = function() {
      c("plot.pdf")
      },
    content = function(file) {
      ggsave(file, plot = main_plot(v, dataHeat(),
                                   clusteredDataHeat(),
                                   get_input_main ()),
             width = input$width * 0.056458333,
             height = input$height * 0.056458333,
             units = "cm",
             dpi = 300,
             device = "pdf",
             limitsize = FALSE)
    }
  )

  # # get list of same orthologs (hit_IDs) of a selected point
  # sameOrthoIndex <- reactive({
  #   # check input
  #   if (v$doPlot == FALSE) return()
  #
  #   # info
  #   info <- mainpoint_info()
  #   pos <- info[7]
  # })

  # ===========================================================================
  # PLOT var1/var2 SCORE & % OF PRESENT SPECIES DISTRIBUTION ==================
  # ===========================================================================

  # list of available variables for distribution plot -------------------------
  output$selected.distribution <- renderUI({
    if (nchar(input$var1_id) == 0 & nchar(input$var2_id) == 0){
      varList <- "% present taxa"
    } else if (nchar(input$var1_id) == 0 & nchar(input$var2_id) > 0){
      varList <- as.list(c(input$var2_id, "% present taxa"))
    } else if (nchar(input$var1_id) > 0 & nchar(input$var2_id) == 0){
      varList <- as.list(c(input$var1_id,
                           "% present taxa"))
    } else {
      varList <- as.list(c(input$var1_id,
                           input$var2_id,
                           "% present taxa"))
    }

    selectInput("selected_dist",
                "Choose variable to plot:",
                varList,
                varList[1])
  })

  output$var1DistPlot <- renderPlot(width = 512, height = 356, {
    if (input$auto_update == FALSE){
      # Add dependency on the update button
      # (only update when button is clicked)
      input$update_btn

      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        var1_dist_plot(v, distDf(), input$dist_text_size, input$var1_id)
      })
    } else {
      var1_dist_plot(v, distDf(), input$dist_text_size, input$var1_id)
    }
  })

  output$var2DistPlot <- renderPlot(width = 512, height = 356, {
    if (input$auto_update == FALSE){
      # Add dependency on the update button
      # (only update when button is clicked)
      input$update_btn

      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        var2_dist_plot(v, distDf(), input$dist_text_size, input$var2_id)
      })
    } else {
      var2_dist_plot(v, distDf(), input$dist_text_size, input$var2_id)
    }
  })

  # % present species distribution plot =======================================
  output$presSpecPlot <- renderPlot(width = 512, height = 356, {
    if (input$auto_update == FALSE){
      # Add dependency on the update button
      # (only update when button is clicked)
      input$update_btn

      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        pres_spec_plot(v, presSpecAllDt, input$percent, input$dist_text_size)
      })
    } else {
      pres_spec_plot(v, presSpecAllDt, input$percent, input$dist_text_size)
    }
  })

  # render dist_plot.ui -------------------------------------------------------
  output$dist_plot.ui <- renderUI({
    if (v$doPlot == FALSE){
      return()
    } else{
      if (is.null(input$selected_dist)){
        return()
      } else {
        if (input$selected_dist == "% present taxa"){
          withSpinner(plotOutput("presSpecPlot",
                                 width = input$width,
                                 height = input$height))
        } else{
          if (input$selected_dist == input$var1_id){
            withSpinner(plotOutput("var1DistPlot",
                                   width = input$width,
                                   height = input$height))
          } else if (input$selected_dist == input$var2_id){
            withSpinner(plotOutput("var2DistPlot",
                                   width = input$width,
                                   height = input$height))
          }
        }
      }
    }
  })

  # Download distribution plot ------------------------------------------------
  output$plot_download_dist <- downloadHandler(
    filename = function() {
      paste0("distributionPlot.pdf")
      },
    content = function(file) {
      if (input$selected_dist == input$var1_id){
        ggsave(file, plot = var1_dist_plot(v, distDf(),
                                         input$dist_text_size,
                                         input$var1_id),
               dpi = 300, device = "pdf", limitsize = FALSE)
      }
      if (input$selected_dist == input$var2_id){
        ggsave(file, plot = var2_dist_plot(v, distDf(),
                                         input$dist_text_size,
                                         input$var2_id),
               dpi = 300, device = "pdf", limitsize = FALSE)
      }
      if (input$selected_dist == "% present taxa"){
        ggsave(file, plot = pres_spec_plot(v, presSpecAllDt,
                                         input$percent,
                                         input$dist_text_size),
               dpi = 300, device = "pdf", limitsize = FALSE)
      }
    }
  )

  # ===========================================================================
  # PLOT CUSTOMIZED PROFILE ===================================================
  # ===========================================================================

  # check if all genes and all species are selected ---------------------------
  output$same_profile <- reactive({
    if (v$doPlot == FALSE) return(FALSE)
    if (length(input$in_seq[1]) == 0) return(FALSE)
    else{
      if (input$in_seq[1] == "all" & input$in_taxa[1] == "all") return(TRUE)
    }
  })
  outputOptions(output, "same_profile", suspendWhenHidden = FALSE)

  # change label of plot_custom button if auto_update_selected is unchecked ----
  output$plot_custom_btn <- renderUI({
    if (input$auto_update_selected == FALSE){
      shinyBS::bsButton("plot_custom",
                        "Plot/Update selected sequence(s)/taxa",
                        style = "warning")
    } else {
      shinyBS::bsButton("plot_custom",
                        "Plot selected sequence(s)/taxa",
                        style = "warning")
    }
  })

  # check if button is clicked ------------------------------------------------
  vCt <- reactiveValues(doPlotCustom = FALSE)
  observeEvent(input$plot_custom, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    vCt$doPlotCustom <- input$plot_custom
    filein <- input$main_input
    # if(input$demo == TRUE){ filein = 1 }
    if (input$demo_data == "demo" | input$demo_data == "ampk-tor") filein <- 1
    if (is.null(filein)) vCt$doPlotCustom <- FALSE
  })

  # print list of available customized taxonomy ranks -------------------------
  # (the lowest rank is the same as the chosen main rank)
  output$rank_select_cus <- renderUI({
    mainRank <- input$rank_select
    mainChoices <- get_taxonomy_ranks()
    cusChoices <- mainChoices[mainChoices >= mainRank]

    selectInput("rank_select_cus", label = h5("Select taxonomy rank:"),
                choices = as.list(cusChoices),
                selected = mainRank)
  })

  output$taxa_select_cus <- renderUI({
    choice <- taxa_select_cus()
    choice$fullName <- as.factor(choice$fullName)
    selectInput("taxa_select_cus",
                h5("Choose (super)taxon of interest:"),
                as.list(levels(choice$fullName)),
                levels(choice$fullName)[1])
  })

  output$selected_plot <- renderPlot({
    if (input$auto_update_selected == FALSE){
      # Add dependency on the update button
      # (only update when button is clicked)
      input$plot_custom

      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        selected_plot(vCt,
                      dataHeat(),
                      clusteredDataHeat(),
                      get_input_selected())
      })
    } else {
      selected_plot(vCt,
                    dataHeat(),
                    clusteredDataHeat(),
                    get_input_selected())
    }
  })

  # plot selected sequences heatmap -------------------------------------------
  output$selected_plot.ui <- renderUI({
    if (is.null(input$in_seq[1]) | is.null(input$in_taxa[1]))  return()
    else if (input$in_seq[1] == "all" & input$in_taxa[1] == "all") return()
    else{
      if (input$auto_update_selected == FALSE){
        # Add dependency on the update button
        # (only update when button is clicked)
        input$plot_custom

        # Add all the filters to the data based on the user inputs
        # wrap in an isolate() so that the data won't update every time an input
        # is changed
        isolate({
          withSpinner(
            plotOutput("selected_plot",
                       width = input$selected_width,
                       height = input$selected_height,
                       click = "plot_click_selected",
                       hover = hoverOpts(
                         id = "plot_hover_selected",
                         delay = input$hover_delay,
                         delayType = input$hover_policy,
                         nullOutside = input$hover_null_outside
                       )
            )
          )
        })
      } else {
        withSpinner(
          plotOutput("selected_plot",
                     width = input$selected_width,
                     height = input$selected_height,
                     click = "plot_click_selected",
                     hover = hoverOpts(
                       id = "plot_hover_selected",
                       delay = input$hover_delay,
                       delayType = input$hover_policy,
                       nullOutside = input$hover_null_outside
                     )
          )
        )
      }
    }
  })

  # download selected plot ----------------------------------------------------
  output$selected_download <- downloadHandler(
    filename = function() {
      c("selected_plot.pdf")
      },
    content = function(file) {
      ggsave(file, plot = selected_plot(vCt,
                                        dataHeat(),
                                        clusteredDataHeat(),
                                        get_input_selected()),
             width = input$selected_width * 0.056458333,
             height = input$selected_height * 0.056458333,
             units = "cm", dpi = 300, device = "pdf", limitsize = FALSE)
    }
  )

  # ===========================================================================
  # SHOW CLICKED POINT INFO ===================================================
  # ===========================================================================

  # get value of point_info for activating Detailed Plot button ---------------
  output$point_info_status <- reactive({
    if (input$tabs == "Main profile"){
      # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
      info <- mainpoint_info()
    } else if (input$tabs == "Customized profile"){
      info <- selectedpoint_info()
    } else {
      info <- NULL
    }
    is.null(info)
  })
  outputOptions(output, "point_info_status", suspendWhenHidden = FALSE)

  # show info into "point's info" box -----------------------------------------
  output$point_info <- renderText({
    # GET INFO BASED ON CURRENT TAB
    if (input$tabs == "Main profile"){
      # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
      info <- mainpoint_info()
    } else if (input$tabs == "Customized profile"){
      info <- selectedpoint_info()
    } else {
      return ()
    }

    if (is.null(info)) return()
    else{
      orthoID <- info[2]
      # ## parse orthoID for oneSeq
      # if(input$input_type == "Concatenated fasta file"){
      #   orthoIDTmp <- unlist(strsplit(toString(info[2]),"\\|"))
      #   #orthoID = toString(paste0(orthoIDTmp[2],":",orthoIDTmp[3]))
      #   orthoID = toString(orthoIDTmp[3])
      # }
      if (is.na(orthoID)) return()
      else{
        # if(orthoID=="NA"){orthoID <- info[2]}

        ## print output
        a <- toString(paste("Seed-ID:", info[1]))
        b <- toString(paste0("Hit-ID: ",
                             orthoID,
                             " (",
                             substr(info[3], 6, nchar(info[3])),
                             ")"))
        c <- ""
        if (input$var1_id != ""){
          c <- toString(paste(input$var1_aggregate_by,
                              input$var1_id,
                              ":",
                              info[4]))
        }
        d <- ""
        if (input$var2_id != ""){
          d <- toString(paste(input$var2_aggregate_by,
                              input$var2_id,
                              ":",
                              info[6]))
        }
        e <- toString(paste("% present taxa:", info[5]))
        paste(a, b, c, d, e, sep = "\n")
      }
    }
  })

  # ===========================================================================
  # DETAILED PLOT =============================================================
  # ===========================================================================

  output$detail_plot <- renderPlot({
    p <- detail_plot(v, detail_plotDt(), input$detailed_text,
                     input$var1_id, input$var2_id)
    p
  })

  # plot detailed bar chart ---------------------------------------------------
  output$detail_plot.ui <- renderUI({
    withSpinner(
      plotOutput("detail_plot",
                 width = 800,
                 height = input$detailed_height,
                 click = "plot_click_detail",
                 hover = hoverOpts(
                   id = "plot_hover_2",
                   delay = input$hover_delay,
                   delayType = input$hover_policy,
                   nullOutside = input$hover_null_outside
                 )
      )
    )
  })

  # download detailed plot ----------------------------------------------------
  output$download_detailed <- downloadHandler(
    filename = function() {
      c("detailedPlot.pdf")
      },
    content = function(file) {
      g <- detail_plot(v, detai_plotDt(), input$detailed_text,
                       input$var1_id, input$var2_id)
      ggsave(file,
             plot = g,
             width = 800 * 0.056458333,
             height = input$detailed_height * 0.056458333,
             units = "cm",
             dpi = 300,
             device = "pdf",
             limitsize = FALSE)
    }
  )

  # SHOW info when clicking on detailed plot ----------------------------------
  output$detail_click <- renderText({
    info <- point_infoDetail() # info = seedID, orthoID, var1

    if (is.null(info)) paste("select ortholog")
    else{
      a <- paste0("seedID = ", info[1])
      b <- paste0("hitID = ", info[2])
      c <- ""
      if (input$var1_id != ""){
        c <- paste0(input$var1_id, " = ", info[3])
      }
      d <- ""
      if (input$var2_id != ""){
        d <- paste0(input$var2_id, " = ", info[4])
      }
      paste(a, b, c, d, sep = "\n")
    }
  })

  # FASTA sequence ------------------------------------------------------------
  output$fasta <- renderText({
    if (v$doPlot == FALSE) return()

    info <- point_infoDetail() # info = seedID, orthoID, var1
    if (is.null(info)) return()
    else{
      data <- get_data_filtered()

      seqID <- toString(info[2])
      groupID <- toString(info[1])
      ncbiID <- gsub("ncbi", "", toString(info[5]))

      seqDf <- data.frame("geneID" = groupID,
                          "orthoID" = seqID,
                          "ncbiID" = ncbiID)
      fastaOut <- fasta_out_data(seqDf, input$main_input, input$demo_data,
                               input$input_type, input$one_seq_fasta,
                               input$path,
                               input$dir_format,
                               input$file_ext,
                               input$id_format,
                               get_long_matrix())
      return(paste(fastaOut[1]))
    }
  })

  # ===========================================================================
  # FEATURE ARCHITECTURE PLOT =================================================
  # ===========================================================================

  # check domain file ---------------------------------------------------------
  output$check_domain_files <- renderUI({
    fileDomain <- getDomainFile()
    if (fileDomain == "noFileInput"){
      em("Domain file not provided!!")
    } else if (fileDomain == "noFileInFolder"){
      msg <- paste0(
        "<p><em>Domain file not found!! </em></p>
        <p><em>Please make sure that file name has to be in this format:
        <strong>&lt;seedID&gt;.extension</strong>, where extension is limited to
        <strong>txt</strong>, <strong>csv</strong>, <strong>list</strong>, <strong>domains</strong> or <strong>architecture</strong>.
        </em></p>"
      )
      HTML(msg)
    } else if (fileDomain == "noSelectHit"){
      em("Please select one ortholog sequence!!")
    }
  })

  # check clicked -------------------------------------------------------------
  v3 <- reactiveValues(doPlot3 = FALSE)
  observeEvent(input$do_domain_plot, {
    v3$doPlot3 <- input$do_domain_plot
    filein <- input$main_input
    # if(input$demo == TRUE){ filein = 1 }
    if (input$demo_data == "demo" | input$demo_data == "ampk-tor") filein <- 1
    if (is.null(filein)) v3$doPlot3 <- FALSE
  })

  output$archi_plot <- renderPlot({
    g <- archi_plot(v3,
                    point_infoDetail(),
                    get_domain_information(),
                    input$one_seq_fasta,
                    input$label_archi_size,
                    input$title_archi_size)
    grid.draw(g)
  })

  # render domain architecture plot -------------------------------------------
  output$archi_plot.ui <- renderUI({
    if (v3$doPlot3 == FALSE) {
      domainIN <- unlist(strsplit(toString(input$main_input), ","))
      fileName <- toString(domainIN[1])
      msg <- paste0(
        "<p><strong>No information about domain architecture! Please check:</strong></p>
        <ul>
        <li>if you uploaded the correct domain file/folder using <span style=\"color: #0000ff;\"><em>Upload additional input</em></span> option? (see input example in <span style=\"background-color: #999999; color: #ffffff;\">data/demo/domains/</span> folder); or</li>
        <li>if&nbsp;the selected genes (seed &amp; ortholog) do exist in the uploaded file (please search for the corresponding&nbsp;<span style=\"text-decoration: underline;\"><em>seedID</em></span> and <span style=\"text-decoration: underline;\"><em>hitID</em></span>, which are&nbsp;shown in <span style=\"color: #ffffff; background-color: #999999;\">Detailed plot</span>)</li>
        </ul>"
        )
      HTML(msg)
    } else {
      withSpinner(plotOutput("archi_plot",
                             height = input$archi_height,
                             width = input$archi_width))
    }
  })

  # download architecture plot ------------------------------------------------
  # ***** something strange with archi_plot()
  output$archi_download <- downloadHandler(
    filename = function() {
      c("domains.pdf")
      },
    content = function(file) {
      g <- archi_plot(v3,
                      point_infoDetail(),
                      get_domain_information(),
                      input$one_seq_fasta,
                      input$label_archi_size,
                      input$title_archi_size)
      grid.draw(g)
      ggsave(file, plot = g,
             width = input$selected_width * 0.056458333,
             height = input$selected_height * 0.056458333,
             units = "cm", dpi = 300, device = "pdf", limitsize = FALSE)
    }
  )

  # ===========================================================================
  # GENE AGE ==================================================================
  # ===========================================================================

  output$gene_agePlot <- renderPlot({
    if (input$auto_update == FALSE){
      # Add dependency on the update button
      # (only update when button is clicked)
      input$update_btn

      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        gene_age_plot(gene_ageDfMod(), input$gene_age_text)
      })
    } else {
      gene_age_plot(gene_ageDfMod(), input$gene_age_text)
    }
  })

  output$gene_age.ui <- renderUI({
    if (v$doPlot == FALSE){
      return()
    } else{
      ## if auto_update is NOT selected, use update_btn to trigger plot changing
      if (input$auto_update == FALSE){
        # Add dependency on the update button
        # (only update when button is clicked)
        input$update_btn

        # Add all the filters to the data based on the user inputs
        # wrap in an isolate() so that the data won't update every time an input
        # is changed
        isolate({
          withSpinner(
            plotOutput("gene_agePlot",
                       width = 600 * input$gene_age_width,
                       height = 150 * input$gene_age_height,
                       click = "plot_click_gene_age")
          )
        })
      }
      ## if auto_update is true
      else {
        withSpinner(
          plotOutput("gene_agePlot",
                     width = 600 * input$gene_age_width,
                     height = 150 * input$gene_age_height,
                     click = "plot_click_gene_age")
        )
      }
    }
  })

  # download gene age plot ----------------------------------------------------
  output$gene_age_plot_download <- downloadHandler(
    filename = function() {
      "gene_age_plot.pdf"
      },
    content = function(file) {
      ggsave(file, plot = gene_age_plot(gene_ageDfMod(),
                                       input$gene_age_text),
             width = 600 * input$gene_age_width * 0.056458333,
             height = 150 * input$gene_age_height * 0.056458333,
             units = "cm", dpi = 300, device = "pdf")
    }
  )

  output$gene_age.table <- renderTable({
    if (is.null(input$plot_click_gene_age$x)) return()

    data <- as.data.frame(selectedgene_age())
    data$number <- rownames(data)
    colnames(data) <- c("geneID", "No.")
    data <- data[, c("No.", "geneID")]
    data
  })

  # download gene list from gene_ageTable -------------------------------------
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

  # check if anywhere elese genes are added to the custemized profile ---------
  observe({
    if (input$add_cluster_cutom_profile == TRUE
       | input$add_cons_gene_custom_profile == TRUE
       | input$add_gc_genes_custom_profile == TRUE ){
      shinyjs::disable("add_custom_profile")
    } else {
      shinyjs::enable("add_custom_profile")
    }
  })

  output$add_custom_profile_check.ui <- renderUI({
    if (input$add_cluster_cutom_profile == TRUE
       | input$add_cons_gene_custom_profile == TRUE
       | input$add_gc_genes_custom_profile == TRUE){
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Profile clustering</strong> or <strong>Core genes finding</strong> or <strong>Groupcomparison</strong> &nbsp;to enable this function)</em></p>')
    }
  })

  # reset gene_age_prot_config ------------------------------------------------
  observeEvent(input$reset_gene_age_prot_config, {
    shinyjs::reset("gene_age_width")
    shinyjs::reset("gene_age_height")
    shinyjs::reset("gene_age_text")
  })

  # ===========================================================================
  # CORE GENES ================================================================
  # ===========================================================================

  # render list of available taxa ---------------------------------------------
  output$taxa_list_cons.ui <- renderUI({
    filein <- input$main_input
    # if(input$demo == TRUE){
    if (input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      filein <- 1
    }

    if (is.null(filein)){
      return(selectInput("in_taxa",
                         "Select taxa of interest:",
                         "none"))
      }
    if (v$doPlot == FALSE){
      return(selectInput("in_taxa",
                         "Select taxa of interest:",
                         "none"))
      }
    else{
      choice <- alltaxa_list()
      choice$fullName <- as.factor(choice$fullName)

      out <- as.list(levels(choice$fullName))
      out <- append("none", out)

      if (input$apply_cons_taxa == TRUE){
        out <- consTaxaName()
        selectInput("taxaCons",
                    "Select taxa of interest:",
                    out,
                    selected = out,
                    multiple = TRUE)
      } else {
        selectInput("taxaCons",
                    "Select taxa of interest:",
                    out,
                    selected = out[1],
                    multiple = TRUE)
      }
    }
  })

  output$cons_gene.table <- renderDataTable({
    data <- cons_geneDf()
    if (is.null(data)) return()
    else {
      data <- as.data.frame(data)
      # data$number <- rownames(data)
      # colnames(data) <- c("geneID","No.")
      # data <- data[,c("No.","geneID")]
      data
    }
  })

  # download gene list from cons_gene.table -----------------------------------
  output$cons_gene_table_download <- downloadHandler(
    filename = function(){
      c("consensusGeneList.out")
      },
    content = function(file){
      data_out <- cons_geneDf()
      write.table(data_out, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  # check if anywhere elese genes are added to the custemized profile ---------
  observe({
    if (input$add_cluster_cutom_profile == TRUE
       | input$add_custom_profile == TRUE
       | input$add_gc_genes_custom_profile == TRUE){
      shinyjs::disable("add_cons_gene_custom_profile")
    } else {
      shinyjs::enable("add_cons_gene_custom_profile")
    }
  })

  output$add_cons_gene_custom_profile_check.ui <- renderUI({
    if (input$add_cluster_cutom_profile == TRUE
       | input$add_custom_profile == TRUE
       | input$add_gc_genes_custom_profile == TRUE){
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Profiles clustering</strong> or <strong>Gene age estimating</strong> or r <strong>Group Comparioson</strong>&nbsp;to enable this function)</em></p>')
    }
  })

  # print list of available taxonomy ranks ------------------------------------
  # (the lowest rank is the same as the chosen main rank)
  output$rank_select_cons <- renderUI({
    mainRank <- input$rank_select
    mainChoices <- get_taxonomy_ranks()
    consChoices <- mainChoices[mainChoices >= mainRank]

    selectInput("rank_select_cons", label = h5("Select taxonomy rank:"),
                choices = as.list(consChoices),
                selected = mainRank)
  })

  output$taxa_select_cons <- renderUI({
    choice <- taxa_select_cons()
    choice$fullName <- as.factor(choice$fullName)
    selectInput("taxa_select_cons",
                h5("Choose (super)taxon of interest:"),
                as.list(levels(choice$fullName)),
                levels(choice$fullName)[1])
  })

  # ===========================================================================
  # CLUSTERING PROFILES =======================================================
  # ===========================================================================

  output$dendrogram <- renderPlot({
    if (v$doPlot == FALSE) return()
    dendrogram(clusterDataDend())
  })

  output$cluster.ui <- renderUI({
    withSpinner(
      plotOutput("dendrogram",
                 width = input$cluster_plot.width,
                 height = input$cluster_plot.height,
                 brush = brushOpts(
                   id = "plot_brush",
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

  # check if anywhere elese genes are added to the custemized profile ---------
  observe({
    if (input$add_custom_profile == TRUE
       | input$add_cons_gene_custom_profile == TRUE
       | input$add_gc_genes_custom_profile == TRUE){
      shinyjs::disable("add_cluster_cutom_profile")
    }else{
      shinyjs::enable("add_cluster_cutom_profile")
    }
  })

  output$add_cluster_cutom_profile_check.ui <- renderUI({
    if (input$add_custom_profile == TRUE
       | input$add_cons_gene_custom_profile == TRUE |
       input$add_gc_genes_custom_profile == TRUE ){
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Gene age estimation</strong> or <strong>Core genes finding</strong> or <strong>Group Comparison</strong> &nbsp;to enable this function)</em></p>')
    }
  })

  # ===========================================================================
  # FILTERED DATA FOR DOWNLOADING =============================================
  # ===========================================================================

  # FOR MAIN PROFILE ==========================================================

  # render variable used for identifying representative genes -----------------
  output$ref_var_main.ui <- renderUI({
    if (nchar(input$var2_id) < 1 & nchar(input$var1_id) < 1){
      radioButtons(inputId = "ref_var_main", label = "Reference variable",
                   choices = list(input$var1_id, input$var2_id),
                   selected = input$var1_id)
    } else if (nchar(input$var2_id) < 1){
      radioButtons(inputId = "ref_var_main",
                   label = "Reference variable",
                   choices = list(input$var1_id),
                   selected = input$var1_id)
    } else {
      radioButtons(inputId = "ref_var_main",
                   label = "Reference variable",
                   choices = list(input$var1_id, input$var2_id),
                   selected = input$var1_id)
    }
  })

  # download data -------------------------------------------------------------
  output$download_data <- downloadHandler(
    filename = function(){
      c("filteredData.out")
      },
    content = function(file){
      data_out <- download_data()
      write.table(data_out, file,
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)
    }
  )

  # data table ui tab ---------------------------------------------------------
  output$filtered_main_data <- renderDataTable(rownames = FALSE, {
    if (v$doPlot == FALSE) return()
    #data <- taxa_id()
    #data <- alltaxa_list()
    #data <- sortedtaxa_list()
    #data <- preData()
    #data <- get_data_filtered()
    #data <- dataSupertaxa()
    # data <- dataHeat()
    #data <- detail_plotDt()
    #data <- presSpecAllDt()
    #data <- distDf()
    #data <- gene_ageDf()
    data <- download_data()

    data
  })

  # download FASTA ------------------------------------------------------------
  output$download_fasta.ui <- renderUI({
    # if(input$demo_data == "demo"){
    #   HTML("<p><span style=\"color: #ff0000;\"><em>Be patient! For large number of taxa this can take up to 3 minutes!</em></span></p>")
    # }
  })

  output$download_fasta <- downloadHandler(
    filename = function(){
      c("filteredSeq.fa")
      },
    content = function(file){
      fasta_out_df <- fasta_out_data(as.data.frame(download_data()),
                                 input$main_input, input$demo_data,
                                 input$input_type, input$one_seq_fasta,
                                 input$path,
                                 input$dir_format,
                                 input$file_ext,
                                 input$id_format,
                                 get_long_matrix())
      write.table(fasta_out_df, file,
                  sep = "\t",
                  col.names = FALSE,
                  row.names = FALSE,
                  quote = FALSE)
    }
  )

  # FOR CUSTOMIZED PROFILE ====================================================

  # render variable used for identifying representative genes -----------------
  output$representative_info.ui <- renderUI({
    msg <- paste0("NOTE: According to your choice in [Download filtered data -> Main data], only representative sequences with ",
                  as.character(input$ref_type_main),
                  " ",
                  as.character(input$ref_var_main),
                  "  will be downloaded!")
    strong(em(msg), style = "color:red")
  })


  # download data -------------------------------------------------------------
  output$download_custom_data <- downloadHandler(
    filename = function(){
      c("customFilteredData.out")
      },
    content = function(file){
      data_out <- download_custom_data()
      write.table(data_out, file,
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)
    }
  )

  # data table ui tab ---------------------------------------------------------
  output$filtered_custom_data <- renderDataTable(rownames = FALSE, {
    if (v$doPlot == FALSE) return()
    data <- download_custom_data()
    data
  })

  # download FASTA ------------------------------------------------------------
  output$download_custom_fasta.ui <- renderUI({
    # if(input$demo_data == "demo"){
    #   HTML("<p><span style=\"color: #ff0000;\"><em>Depend on the number of taxa, this might take up to 3 minutes!</em></span></p>")
    # }
  })

  output$download_custom_fasta <- downloadHandler(
    filename = function(){
      c("customFilteredSeq.fa")
    },
    content = function(file){
      fasta_out_df <- fasta_out_data(as.data.frame(download_custom_data()),
                                 input$main_input, input$demo_data,
                                 input$input_type, input$one_seq_fasta,
                                 input$path,
                                 input$dir_format,
                                 input$file_ext,
                                 input$id_format,
                                 get_long_matrix())
      write.table(fasta_out_df,
                  file,
                  sep = "\t",
                  col.names = FALSE,
                  row.names = FALSE,
                  quote = FALSE)
    }
  )

  # ===========================================================================
  # GROUP COMPARISON ==========================================================
  # ===========================================================================

  # Dataframe with Information about the significant Genes --------------------
  # geneID | in_group| out_group | pvalues | features | databases | rank | var
  significant_genes_gc <- NULL

  # Select in_group -----------------------------------------------------------
  output$taxa_list_gc <- renderUI({
    filein <- input$main_input
    if (input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      filein <- 1
    }

    if (is.null(filein)){
      return(selectInput("in_taxa", "Select in_group:", "none"))
    }

    if (v$doPlot == FALSE){
      return(selectInput("in_taxa", "Select in_group:", "none"))
    }else{

      # When the taxonomy rank was changed ------------------------------------
      if (input$apply_taxa_gc == TRUE){
        choice <- taxa_select_gc(input$rank_select_gc, subset_taxa())
        choice$fullName <- as.factor(choice$fullName)
        out <- as.list(levels(choice$fullName))

        x <- out[out == input$taxa_select_gc]
        selectInput("selected_in_group_gc", "Select in_group:",
                    out,
                    selected = x,
                    multiple = TRUE,
                    selectize = FALSE)

      # When the taxonomy is the same as in input and settings ----------------
      } else {

        # check for the rank of the rank in the input
        ranks <- get_taxonomy_ranks()
        pos <- which(ranks == input$rank_select) # position in the list
        higher_rank <- ranks[pos + 1] # take the next higher rank
        higher_rank_name <- substr(higher_rank, 4, nchar(higher_rank))

        name_list <- get_name_list(TRUE, TRUE) # get the taxon names
        dt <- get_taxa_list(FALSE, subset_taxa) # get the taxa

        # get the info for the reference protein from the namelist
        reference <- subset(name_list,
                            name_list$fullName == input$in_select)

        # get the id for every rank for the reference protein
        rank_name <- substr(input$rank_select, 4, nchar(input$rank_select))
        reference_dt <- dt[dt[, rank_name] == reference$ncbiID, ]

        # save the next higher rank
        reference_higher_rank <- reference_dt[higher_rank_name]
        reference_higher_rank <- {
          reference_higher_rank[!duplicated(reference_higher_rank), ]
        }

        # get all the taxa with the same id in the next higher rank
        selected_taxa_dt <- {
          subset(dt, dt[, higher_rank_name] %in% reference_higher_rank)
        }
        selected_taxa_dt <- {
          selected_taxa_dt[!duplicated(selected_taxa_dt[rank_name]), ]
        }

        # get a list with all the ids with reference_higher_rank as parent
        selected_taxa_ids <- selected_taxa_dt[rank_name]
        if (length(selected_taxa_ids[[1]]) >= 1){
          selected_taxa_ids <- selected_taxa_ids[[1]]
        }

        selected_taxa <- subset(name_list, name_list$rank == rank_name)
        selected_taxa <- {
          subset(selected_taxa, selected_taxa$ncbiID %in% selected_taxa_ids)
        }

        default_select <- selected_taxa$fullName

        choice <- taxa_select_gc(input$rank_select, subset_taxa())
        choice$fullName <- as.factor(choice$fullName)
        out <- as.list(levels(choice$fullName))

        selectInput("selected_in_group_gc", "Select in_group:",
                    out,
                    selected = default_select,
                    multiple = TRUE,
                    selectize = FALSE)
      }
    }
  })

  # Functions: Popup Window Select Rank ---------------------------------------
  # print list of available taxonomy ranks
  #(the lowest rank is the same as the chosen main rank)
  output$rank_select_gc <- renderUI({
    main_rank <- input$rank_select
    main_choices <- get_taxonomy_ranks()
    choices_gc <- main_choices[main_choices >= main_rank]

    selectInput("rank_select_gc", label = h5("Select taxonomy rank:"),
                choices = as.list(choices_gc),
                selected = main_rank)
  })

  # Supertaxon of intrest in the popup window for the rank
  output$taxa_select_gc <- renderUI({
    choice <- taxa_select_gc(input$rank_select_gc, subset_taxa())
    choice$fullName <- as.factor(choice$fullName)
    selectInput("taxa_select_gc", h5("Choose (super)taxon of interest:"),
                as.list(levels(choice$fullName)),
                levels(choice$fullName)[1])
  })

  # Select Gene ---------------------------------------------------------------
  output$list_genes_gc <- renderUI({
    filein <- input$main_input

    file_gc <- input$gc_file

    if (input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      filein <- 1
    }

    if (is.null(filein) & is.null(file_gc)){
      return(selectInput("list_selected_genes_gc", "Select sequence:", "none"))
    }else{
      # full list
      data <- as.data.frame(get_data_filtered())
      data$geneID <- as.character(data$geneID)
      data$geneID <- as.factor(data$geneID)
      out_all <- as.list(levels(data$geneID))
      out_all <- append("all", out_all)

      if (is.null(file_gc)){
        selectInput("list_selected_genes_gc", "Select genes:",
                    out_all,
                    selected = out_all[1],
                    multiple = TRUE,
                    selectize = FALSE)
      }else {
        list_gc <- as.data.frame(read.table(file = file_gc$datapath,
                                            header = FALSE))
        list_gc$V1 <- as.factor(list_gc$V1)
        out <- as.list(levels(list_gc$V1))
        selectInput("list_selected_genes_gc", "Select genes:",
                    out,
                    selected = NULL,
                    multiple = FALSE,
                    selectize = FALSE)
      }
    }
  })

  # Buttons to choose the variable --------------------------------------------
  output$variable_button_gc <- renderUI ({
    popify(
      radioButtons(
        inputId = "var_name_gc",
        label = "Select Variable:",
        choices = list(input$var1_id, input$var2_id, "Both"),
        selected = input$var1_id,
        inline = F),
      "",
      "Select variable(s) to generate plots for")
  })

  # Slider to choose significance ---------------------------------------------
  output$significance.ui <- renderUI({
    popify(
      sliderInput("significance", paste("Significance level:"),
                  min = 0,
                  max = 1,
                  step = 0.005,
                  value = c(0.05),
                  width = 200),
      "",
      "maximal probability to reject that in_group and Out-Group have no significant difference by mistake")
  })

  # List with all significant Genes -------------------------------------------
  output$get_significant_genes <- renderUI({
    input$plot_gc

    isolate({
      significant_genes_gc <<- {
        get_significant_genes(input$selected_in_group_gc,
                              input$list_selected_genes_gc,
                              input$rank_select,
                              input$var_name_gc,
                              input$use_common_anchestor,
                              input$in_select,
                              get_input_gc(),
                              input$demo_data,
                              input$anno_choose,
                              input$file_domain_input,
                              subset_taxa(),
                              get_data_filtered(),
                              session,
                              input$right_format_features,
                              get_domain_information())
      }
      if (!is.null(significant_genes_gc)){
        x <- as.vector(significant_genes_gc$geneID)
        choices <- c("all", x)

        # selected Gene
        selectInput("selected_gene_gc", "Candidate gene(s):",
                    choices,
                    selected = choices[2],
                    multiple = FALSE)

      } else{
        # selected Gene
        selectInput("selected_gene_gc", "Candidate gene(s):",
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
    input$selected_gene_gc
    isolate({
      gene <- input$selected_gene_gc
      if (!input$right_format_features){
        selectInput("interesting_features", "Feature type(s) of interest:",
                    NULL,
                    selected = NULL,
                    multiple = TRUE,
                    selectize = FALSE)
      } else if (is.null(gene)){
        selectInput("interesting_features", "Feature type(s) of interest:",
                    NULL,
                    selected = NULL,
                    multiple = TRUE,
                    selectize = FALSE)
      } else if (gene == ""){
        selectInput("interesting_features", "Feature type(s) of interest:",
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
        selectInput("interesting_features", "Feature type(s) of interest:",
                    choices,
                    selected = choices[1],
                    multiple = TRUE,
                    selectize = FALSE)
      }
    })
  })

  # Select Plots to download --------------------------------------------------
  output$select_plots_to_download  <- renderUI({
    input$plot_gc
    isolate({
      if (!is.null(significant_genes_gc)){
        x <- as.vector(significant_genes_gc$geneID)
        choice <- c("all", x)
        gene <- subset(choice, choice == input$selected_gene_gc)

        selectInput("plots_to_download", "Select Plots to download:",
                    choice,
                    selected = gene,
                    multiple = TRUE,
                    selectize = FALSE)

      } else{
        selectInput("plots_to_download", "Select Plots to download:",
                    NULL,
                    selected = NULL,
                    multiple = TRUE,
                    selectize = FALSE)
      }
      })
  })

  # observe Events for the Appearance of the plots ============================
  # reset config of customized plot
  observeEvent(input$reset_config_gc, {
    shinyjs::reset("x_size_gc")
    shinyjs::reset("y_size_gc")
    shinyjs::reset("angle_gc")
    shinyjs::reset("legend_size_gc")
  })

  # close customized config
  observeEvent(input$apply_config_gc, {
    toggleModal(session, "gc_plot_config_bs", toggle = "close")
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
        get_multiplot_download_gc(gene, get_input_gc(),
                                  input$interesting_features)
        dev.off()
      }
      zip(zipfile = file, files = fs)
    },
    contentType = "application/zip"
  )

  # observer for the download functions
  observe({
    if (is.null(input$selected_in_group_gc)
        | length(input$list_selected_genes_gc) == 0){
      shinyjs::disable("download_plots_gc")
      shinyjs::disable("download_genes_gc")
    }else if (input$plot_gc == FALSE){
      shinyjs::disable("download_plots_gc")
      shinyjs::disable("download_genes_gc")
    }else if (input$selected_gene_gc == "") {
      shinyjs::disable("download_plots_gc")
      shinyjs::disable("download_genes_gc")
    }else{
      shinyjs::enable("download_plots_gc")
      shinyjs::enable("download_genes_gc")
    }
  })

  # add_custom_profile --------------------------------------------------------
  # check if add_custom_profile (gene age plot) are being clicked
  observe({
    if (input$add_custom_profile == TRUE |
        input$add_cons_gene_custom_profile == TRUE |
        input$add_cluster_cutom_profile == TRUE){
      shinyjs::disable("add_gc_genes_custom_profile")
    }else{
      shinyjs::enable("add_gc_genes_custom_profile")
    }
  })
  output$add_gc_custom_profile_check <- renderUI({
    if (input$add_custom_profile == TRUE |
        input$add_cons_gene_custom_profile == TRUE |
        input$add_cluster_cutom_profile == TRUE){
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Gene age estimation</strong>or <strong>Profile clustering</strong> or <strong>Core genes finding</strong>&nbsp;to enable this function)</em></p>')
    }
  })

  # ===========================================================================
  # REACTIVE FUCTIONS =========================================================
  # ===========================================================================
  
  # Working with the input file(s) ============================================
  # Get the type of the input file --------------------------------------------
  get_input_type <- reactive({
    filein <- input$main_input
    if (is.null(filein)) return()
    input_type <- check_input_vadility(filein)
    
    return(input_type)
  })
  
  # Turn the input in to a matrix in long format ------------------------------
  get_long_matrix <- reactive({
    if (input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      if (input$demo_data == "demo"){
        long_dataframe <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.long",
                                     sep = "\t",
                                     header = T,
                                     fill = T,
                                     stringsAsFactors = FALSE)
      } else {
        long_dataframe <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile",
                                     sep = "\t",
                                     header = T,
                                     fill = T,
                                     stringsAsFactors = FALSE)
      }
    } else {
      filein <- input$main_input
      if(is.null(filein)) return()
      input_type <- get_input_type()
      
      # XML
      if (input_type == "xml"){
        long_dataframe <- xml_parser(filein$datapath)
        
        # FASTA
      } else if (input_type == "fasta"){
        long_dataframe <- fasta_parser(filein$datapath)
        
        # LONG
      } else if (input_type == "long"){
        long_dataframe <- as.data.frame(read.table(file = filein$datapath,
                                                   sep = "\t",
                                                   header = T,
                                                   check.names = FALSE,
                                                   comment.char = ""))
        # WIDE
      } else if (input_type == "wide"){
        long_dataframe <- wide_to_long(filein$datapath)
        
      } else if (input_type == "oma"){
        # dont load the data before the button is loaded
        if(is.null(input$get_data_oma)) return()
  
        isolate({
          # dont generate data befor the button was clicked
          if (input$get_data_oma[1] == 0) return() 
          oma_type <- input$selected_oma_type
          if(is.null(oma_type))return()
          oma_ids <- as.data.frame(read.table(file = filein$datapath,
                                              sep = "\t",
                                              header = F,
                                              check.names = FALSE,
                                              comment.char = ""))

          oma_ids[, 1] <- as.character(oma_ids[, 1])
          print(oma_type)
          long_dataframe <- oma_ids_to_long(oma_ids[, 1], oma_type)
        })
       
        
      # When there is no  usable data (moreCol, emptyCell, noGeneID)
      }else{
        return (NULL)
      }
    }
    
    # make sure the cells have the same type (factor) 
    # independent of the input file 
    for (i in 1:ncol(long_dataframe)){
      long_dataframe[, i] <- as.factor(long_dataframe[, i])
    }
    return(long_dataframe)
  })
  
  # Save the domain files as a data frame -------------------------------------
  get_domain_information <- reactive({
    domains <- data.frame()
    print(domains)
    files <- c()
  
    if(!is.null(get_input_type())){
      if(get_input_type() == "oma"){
        domains <- long_to_domain(get_long_matrix())
        domains$seedID <- gsub("\\|",":",domains$seedID)
        domains$orthoID <- gsub("\\|",":",domains$orthoID)
        
        return(domains)
      }
    }
    if(!input$demo_data == "ampk-tor" | !input$anno_choose == "from file" ){
        long_df <- get_long_matrix()
        genes <- unlist(long_df$geneID)
        genes <- unique(genes)
        
        for (gene in genes){
          print(gene)
          file <- get_domain_file_gc(gene, 
                                     input$demo_data,
                                     input$anno_choose,
                                     input$file_domain_input,
                                     session,
                                     input$domainPath)
          files <- append(files, file)
        }
      } else{
        file <- get_domain_file_gc(NULL,
                                   input$demo_data,
                                   input$anno_choose,
                                   input$file_domain_input,
                                   session,
                                   input$domainPath)
        files <- c(file)
      }
      
      # parse domain file
      if (input$demo_data == "demo" | input$demo_data == "ampk-tor"){
        
        for(file in files){
          print(file)
          domain_df <- as.data.frame(read.csv(file,
                                              sep = "\t",
                                              header = F,
                                              comment.char = "",
                                              stringsAsFactors = FALSE,
                                              quote = ""))
          domains <- rbind(domains, domain_df)
        }
        
        if (ncol(domains) == 5){
          colnames(domains) <- c("seedID",
                                 "orthoID",
                                 "feature",
                                 "start",
                                 "end")
        } else if (ncol(domains) == 6){
          colnames(domains) <- c("seedID",
                                 "orthoID",
                                 "feature",
                                 "start",
                                 "end",
                                 "weight")
        } else if (ncol(domains) == 7){
          colnames(domains) <- c("seedID",
                                 "orthoID",
                                 "feature",
                                 "start",
                                 "end",
                                 "weight",
                                 "path")
        }
        domains$length <- max(domains$end)
        
      } else {
        
        for(file in files){
          if (file != FALSE){
            exeptions <- c("noFileInput", "noSelectHit",
                           "noSelectHit", "noFileInFolder")
            if(!(file %in% exeptions)){
              domain_df <- as.data.frame(read.table(file,
                                                    sep = "\t",
                                                    header = FALSE,
                                                    comment.char = ""))
              domains <- rbind(domains, domain_df)
            }
          }
        }
        
        if (ncol(domains) == 6){
          colnames(domains) <- c("seedID",
                                 "orthoID",
                                 "length",
                                 "feature",
                                 "start",
                                 "end")
        } else if (ncol(domains) == 7){
          colnames(domains) <- c("seedID",
                                 "orthoID",
                                 "length",
                                 "feature",
                                 "start",
                                 "end",
                                 "weight")
        } else if (ncol(domains) == 8){
          colnames(domains) <- c("seedID",
                                 "orthoID",
                                 "length",
                                 "feature",
                                 "start",
                                 "end",
                                 "weight",
                                 "path")
        }
      }
    print("finish")

    domains$seedID <- gsub("\\|",":",domains$seedID)
    domains$orthoID <- gsub("\\|",":",domains$orthoID)
    print(head(domains))
    return(domains)
   })
  
  # Save the fasta file(s) as a data frame ------------------------------------
  get_fasta_information <- reactive({
    
  })
  
  # Parameters for the plots ==================================================
  # Parameters for the main profile plot --------------------------------------
  get_input_main <- reactive({
    input_data <- list( "x_axis" = input$x_axis,
                        "var1_id" = input$var1_id,
                        "var2_id"  = input$var2_id,
                        "low_color_var1" =  input$low_color_var1,
                        "high_color_var1" = input$high_color_var1,
                        "low_color_var2" = input$low_color_var2,
                        "high_color_var2" = input$high_color_var2,
                        "para_color" = input$para_color,
                        "x_size" = input$x_size,
                        "y_size" = input$y_size,
                        "legend_size" = input$legend_size,
                        "main_legend" = input$main_legend,
                        "dot_zoom" = input$dot_zoom,
                        "x_angle" = input$x_angle,
                        "taxon_highlight" = input$taxon_highlight,
                        "apply_cluster" = input$apply_cluster,
                        "rank_select" = input$rank_select,
                        "gene_highlight" = input$gene_highlight,
                        "auto_update" = input$auto_update,
                        "update_btn" = input$update_btn)

    return (input_data)
  })

  # Parameters for the customized profile plot --------------------------------
  get_input_selected <- reactive ({
    input_data <- list ("x_axis_selected" = input$x_axis_selected,
                        "var1_id" = input$var1_id,
                        "var2_id" = input$var2_id,
                        "low_color_var1" = input$low_color_var1,
                        "high_color_var1" = input$high_color_var1,
                        "low_color_var2" = input$low_color_var2,
                        "high_color_var2" = input$high_color_var2,
                        "para_color" = input$para_color,
                        "x_size_select" = input$x_size_select,
                        "y_size_select" = input$y_size_select,
                        "legend_size_select" = input$legend_size_select,
                        "selected_legend" = input$selected_legend,
                        "dot_zoom_select" = input$dot_zoom_select,
                        "x_angle_select" = input$x_angle_select,
                        "in_taxa" = input$in_taxa,
                        "in_seq" = input$in_seq,
                        "auto_update_selected" = input$auto_update_selected,
                        "plot_custom" = input$plot_custom,
                        "apply_cluster" = input$apply_cluster
    )
  })

  # Parameters for the plots in Group Comparison ------------------------------
  get_input_gc <- reactive ({
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
  
  # get ncbi taxa IDs ---------------------------------------------------------
  # retrieve ID for list of taxa names
  taxa_id <- reactive({
    if (input$id_search > 0){
      taxain <- input$taxa_list
      if (is.null(taxain)) return()
      
      taxa_name_df <- as.data.frame(read.table(file = taxain$datapath,
                                               sep = "\t",
                                               header = F,
                                               check.names = FALSE,
                                               comment.char = ""))
      
      id_df <- data.frame("name" = character(),
                          "new_name" = character(),
                          "id" = character(),
                          "type" = character(),
                          stringsAsFactors = FALSE)
      
      withProgress(message = "Retrieving IDs...", value = 0, {
        for (i in 1:nrow(taxa_name_df)){
          id <- get_uid(sciname = taxa_name_df[i, ])[1]
          if (is.na(id)){
            temp <- gnr_resolve(names = as.character(taxa_name_df[i, ]))
            if (nrow(temp) > 0){
              new_id <- get_uid(sciname = temp[1, 3])[1]
              if (is.na(new_id)){
                id_df[i, ] <- c(as.character(taxa_name_df[i, ]),
                                as.character(temp[1, 3]),
                                paste0("NA"), "notfound")
              } else {
                id_df[i, ] <- c(as.character(taxa_name_df[i, ]),
                                as.character(temp[1, 3]),
                                paste0("ncbi", new_id),
                                "notfound")
              }
            } else {
              id_df[i, ] <- c(as.character(taxa_name_df[i, ]),
                              paste0("no alternative"),
                              paste0("NA"),
                              "notfound")
            }
          } else {
            id_df[i, ] <- c(as.character(taxa_name_df[i, ]),
                            "NA",
                            paste0("ncbi", id),
                            "retrieved")
          }
          # Increment the progress bar, and update the detail text.
          incProgress(1 / nrow(taxa_name_df),
                      detail = paste(i, "/", nrow(taxa_name_df)))
        }
      })
      # return
      id_df
    }
  })
  
  # get input taxa ------------------------------------------------------------
  subset_taxa <- reactive({
    if (input$demo_data == "demo" |
        input$demo_data == "ampk-tor" |
        length(unkTaxa()) == 0){
      long_dataframe <- get_long_matrix()
      if (is.null(long_dataframe)) return()
      inputTaxa <- levels(long_dataframe$ncbiID)
    } else {
      inputTaxa <- readLines(filein$datapath, n = 1)
    }
    
    inputTaxa <- unlist(strsplit(inputTaxa, split = "\t"))
    if (inputTaxa[1] == "geneID"){
      # remove "geneID" element from vector inputTaxa
      inputTaxa <- inputTaxa[-1]
    }
    # return input taxa
    return(inputTaxa)
  })
  
  # check if there is any "unknown" taxon in input matrix ---------------------
  unkTaxa <- reactive({
    long_dataframe <- get_long_matrix ()
    if (is.null(long_dataframe)){
      inputTaxa <- c("NA")
    } else{
      inputTaxa <- levels(long_dataframe$ncbiID)
    }
    
    if (inputTaxa[1] == "NA"){
      return()
    } else {
      inputTaxa <- unlist(strsplit(inputTaxa, split = "\t"))
      if (inputTaxa[1] == "geneID"){
        # remove "geneID" element from vector inputTaxa
        inputTaxa <- inputTaxa[-1]
      }
      
      if (!file.exists(isolate({
        "data/rankList.txt"
      }
      ))){
        return(inputTaxa)
      } else {
        info <- file.info("data/rankList.txt")
        if (info$size == 0){
          return(inputTaxa)
        } else {
          # get list of all available taxon (from /data/rankList.txt)
          pipe_cmd <- paste0("cut -f 1 ", getwd(), "/data/rankList.txt")
          allTaxa <- unlist( (read.table(pipe(pipe_cmd))))
          
          # list of unknown taxa
          unkTaxa <- inputTaxa[!(inputTaxa %in% allTaxa)]
          # return list of unkTaxa
          return(unkTaxa)
        }
      }
    }
  })
  
  # Get list of all (super)taxa -----------------------------------------------
  alltaxa_list <- reactive({
    
    filein <- input$main_input
    
    # if(is.null(filein) & input$demo == FALSE) return()
    if (is.null(filein) & input$demo_data == "none") return()
    
    rank_select <- input$rank_select
    
    if (rank_select == "") return()
    if (length(unkTaxa()) > 0) return()
    
    # load list of unsorted taxa ----------------------------------------------
    Dt <- get_taxa_list(TRUE, subset_taxa())
    
    # load list of taxon name -------------------------------------------------
    nameList <- get_name_list(TRUE, FALSE)
    
    # get rank name from rank_select
    rankName <- substr(rank_select, 4, nchar(rank_select))
    
    # get rank number (number of column in unsorted taxa list - dataframe Dt)
    #    rankNr = 0 + as.numeric(substr(rank_select,1,2))
    
    choice <- as.data.frame
    choice <- rbind(Dt[rankName])
    colnames(choice) <- "ncbiID"
    choice <- merge(choice,
                    nameList,
                    by = "ncbiID",
                    all = FALSE)
  })
  
  sortedtaxa_list <- reactive({
    if (v$doPlot == FALSE) return()
    
    # Get representative taxon ================================================
    
    # load list of unsorted taxa
    Dt <- get_taxa_list(FALSE, subset_taxa())
    
    # load list of taxon name
    nameList <- as.data.frame(read.table("data/taxonNamesReduced.txt",
                                         sep = "\t",
                                         header = T,
                                         fill = TRUE))
    nameList$fullName <- as.character(nameList$fullName)
    
    # input parameters
    rank_select <- input$rank_select
    
    # get rank name from rank_select
    rankName <- substr(rank_select, 4, nchar(rank_select))
    
    # get rank number (number of column in unsorted taxa list - dataframe Dt)
    rankNr <- 0 + as.numeric(substr(rank_select, 1, 2))
    
    # get selected supertaxon ID
    taxa_list <- as.data.frame(read.table("data/taxonNamesReduced.txt",
                                          sep = "\t",
                                          header = T))
    allTaxa <- alltaxa_list()
    rankNameTMP <- allTaxa$rank[allTaxa$fullName == input$in_select]
    if (rankName == "strain"){
      superID <- {
        as.numeric(taxa_list$ncbiID[taxa_list$fullName == input$in_select
                                    & taxa_list$rank == "norank"])
      }
    } else {
      superID <- {
        as.numeric(taxa_list$ncbiID[taxa_list$fullName == input$in_select
                                    & taxa_list$rank == rankNameTMP[1]])
      }
    }
    
    # representative taxon
    repTaxon <- Dt[Dt[, rankName] == superID, ][1, ]
    
    # Sort taxon list based on hclust tree ====================================
    
    # prepare Df for calculating distance matrix ------------------------------
    distDf <- subset(Dt, select = -c(ncbiID, fullName))
    row.names(distDf) <- distDf$abbrName
    distDf <- distDf[, -1]
    
    # get sorted taxon IDs & sort full taxonomy info dataframe ----------------
    treeIn <- input$inputTree
    
    if (is.null(treeIn)){
      taxaTree <- create_rooted_tree(distDf, as.character(repTaxon$abbrName))
    } else {
      taxaTree <- read.tree(file = treeIn$datapath)
    }
    taxonList <- sort_taxa_from_tree(taxaTree)
    sortedDt <- Dt[match(taxonList, Dt$abbrName), ]
    
    # subset to get list of input taxa only -----------------------------------
    inputTaxa <- subset_taxa()
    sortedDt <- subset(sortedDt, abbrName %in% inputTaxa)
    
    # get only taxonIDs list of selected rank and rename columns --------------
    sortedOut <- subset(sortedDt, select = c("abbrName",
                                             "ncbiID",
                                             "fullName",
                                             as.character(rankName)))
    colnames(sortedOut) <- c("abbrName",
                             "species",
                             "fullName",
                             "ncbiID")
    
    # add name of supertaxa into sortedOut list -------------------------------
    sortedOut <- merge(sortedOut, nameList,
                       by = "ncbiID",
                       all.x = TRUE,
                       sort = FALSE)
    sortedOut$species <- as.character(sortedOut$species)
    
    # add order_prefix to supertaxon name -------------------------------------
    # and add prefix "ncbi" to taxon_ncbiID (column "species")
    prefix <- 1001
    
    ## create new column for sorted supertaxon
    sortedOut$sortedSupertaxon <- 0 
    sortedOut$sortedSupertaxon[1] <- paste0(prefix,
                                            "_",
                                            sortedOut$fullName.y[1])
    sortedOut$species[1] <- paste0("ncbi",
                                   sortedOut$species[1])
    
    for (i in 2:nrow(sortedOut)){
      ## increase prefix if changing to another supertaxon
      if (sortedOut$fullName.y[i] != sortedOut$fullName.y[i - 1]){
        prefix <- prefix + 1
      }
      sortedOut$sortedSupertaxon[i] <- paste0(prefix,
                                              "_",
                                              sortedOut$fullName.y[i])
      sortedOut$species[i] <- paste0("ncbi",
                                     sortedOut$species[i])
    }
    
    # final sorted supertaxa list
    sortedOut$taxonID <- 0
    sortedOut$category <- "cat"
    sortedOut <- sortedOut[, c("abbrName",
                               "taxonID",
                               "fullName.x",
                               "species",
                               "ncbiID",
                               "sortedSupertaxon",
                               "rank",
                               "category")]
    colnames(sortedOut) <- c("abbrName",
                             "taxonID",
                             "fullName",
                             "ncbiID",
                             "supertaxonID",
                             "supertaxon",
                             "rank",
                             "category")
    
    sortedOut$taxonID <- as.numeric(sortedOut$taxonID)
    sortedOut$ncbiID <- as.factor(sortedOut$ncbiID)
    sortedOut$supertaxon <- as.factor(sortedOut$supertaxon)
    sortedOut$category <- as.factor(sortedOut$category)
    
    ### return data frame
    return(sortedOut)
  })
  
  # subset data ---------------------------------------------------------------
  preData <- reactive({
    # get list of gene of interest (from a separated file)
    listGene <- list()
    end_index <- input$end_index
    if (is.na(input$end_index)) end_index <- 30
    if (input$gene_list_selected == "from file"){
      listIn <- input$list
      if (!is.null(listIn)){
        list <- as.data.frame(read.table(file = listIn$datapath,
                                         header = FALSE))
        listGeneOri <- list$V1
        if (input$st_index <= length(listGeneOri)){
          listGene <- listGeneOri[listGeneOri[input$st_index:end_index]]
        } else {
          listGene <- listGeneOri
        }
      }
    }
    
    long_dataframe <- get_long_matrix()
    if (is.null(long_dataframe)){
      data <- data.frame("geneID" = character(),
                         "ncbiID" = character(),
                         "orthoID" = character(),
                         "var1" = character(),
                         "var2" = character(),
                         stringsAsFactors = F)
    } else{
      long_dataframe <- unsort_id(long_dataframe, input$ordering)
      
      if (length(listGene) >= 1){
        data <- long_dataframe[long_dataframe$geneID %in% listGene, ]
      } else {
        subsetID <- levels(long_dataframe$geneID)[input$st_index:end_index]
        data <- long_dataframe[long_dataframe$geneID %in% subsetID, ]
      }
      
      if (ncol(data) < 5){
        for (i in 1:(5 - ncol(data))){
          data[paste0("newVar", i)] <- 1
        }
      }
      colnames(data) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
    }
    # return preData
    return(data)
  })
  
  # get all information for input data ----------------------------------------
  get_data_filtered <- reactive({
    mdData <- preData()

    # count number of inparalogs
    paralogCount <- plyr::count(mdData, c("geneID", "ncbiID"))
    mdData <- merge(mdData, paralogCount, by = c("geneID", "ncbiID"))
    colnames(mdData)[ncol(mdData)] <- "paralog"
    
    # (3) GET SORTED TAXONOMY LIST (3)
    taxa_list <- sortedtaxa_list()
    
    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxa_list, "supertaxon")
    
    # merge mdData, mdDataVar2 and taxa_list to get taxonomy info
    taxaMdData <- merge(mdData, taxa_list, by = "ncbiID")
    taxaMdData$var1 <- {
      suppressWarnings(as.numeric(as.character(taxaMdData$var1)))
    }
    taxaMdData$var2 <- {
      suppressWarnings(as.numeric(as.character(taxaMdData$var2)))
    }
    
    # (4) calculate PERCENTAGE of PRESENT SPECIES (4)
    finalPresSpecDt <- calc_pres_spec(taxaMdData, taxaCount)
    
    # (5) calculate max/min/mean/median VAR1 for every supertaxon of each gene (5)
    # remove NA rows from taxaMdData
    taxaMdDataNoNA <- taxaMdData[!is.na(taxaMdData$var1), ]
    # calculate m var1
    mVar1Dt <- aggregate(taxaMdDataNoNA[, "var1"],
                         list(taxaMdDataNoNA$supertaxon,
                              taxaMdDataNoNA$geneID),
                         FUN = input$var1_aggregate_by)
    colnames(mVar1Dt) <- c("supertaxon", "geneID", "mVar1")
    
    # (6) calculate max/min/mean/median VAR2 for each super taxon (6)
    # remove NA rows from taxaMdData
    taxaMdDataNoNA_var2 <- taxaMdData[!is.na(taxaMdData$var2), ]
    # calculate max/min/mean/median VAR2
    if (nrow(taxaMdDataNoNA_var2) > 0){
      mVar2Dt <- aggregate(taxaMdDataNoNA_var2[, "var2"],
                           list(taxaMdDataNoNA_var2$supertaxon,
                                taxaMdDataNoNA_var2$geneID),
                           FUN = input$var2_aggregate_by)
      colnames(mVar2Dt) <- c("supertaxon", "geneID", "mVar2")
    } else {
      mVar2Dt <- taxaMdData[, c("supertaxon", "geneID")]
      mVar2Dt$mVar2 <- 0
    }
    
    # (5+6) & join mVar2 together with mVar1 scores into one df (5+6)
    scoreDf <- merge(mVar1Dt,
                     mVar2Dt,
                     by = c("supertaxon", "geneID"),
                     all = TRUE)
    
    # (4+5+6) add presSpec and mVar1 into taxaMdData (4+5+6)
    presMdData <- merge(taxaMdData,
                        finalPresSpecDt,
                        by = c("geneID", "supertaxon"),
                        all.x = TRUE)
    fullMdData <- merge(presMdData,
                        scoreDf,
                        by = c("geneID", "supertaxon"),
                        all.x = TRUE)
    fullMdData <- merge(fullMdData,
                        taxaCount, by = ("supertaxon"),
                        all.x = TRUE)
    # rename "freq" into "numberSpec"
    names(fullMdData)[names(fullMdData) == "freq"] <- "numberSpec"
    
    fullMdData$fullName <- as.vector(fullMdData$fullName)
    names(fullMdData)[names(fullMdData) == "orthoID.x"] <- "orthoID"
    # parsed input data frame !!!
    fullMdData <- fullMdData[!duplicated(fullMdData), ]
    
    return(fullMdData)
  })
  
  # REDUCE DATA FROM SPECIES LEVEL TO SUPERTAXA LEVEL -------------------------
  # this data set contain only supertaxa
  # and their value (%present, mVar1 & mVar2) for each gene
  dataSupertaxa <- reactive({
    fullMdData <- get_data_filtered()
    
    # to check if working with the lowest taxonomy rank; 1 for NO; 0 for YES
    flag <- 1
    if (length(unique(levels(as.factor(fullMdData$numberSpec)))) == 1){
      if (unique(levels(as.factor(fullMdData$numberSpec))) == 1){
        superDfExt <- fullMdData[, c("geneID",
                                     "supertaxon",
                                     "supertaxonID",
                                     "var1",
                                     "presSpec",
                                     "category",
                                     "orthoID",
                                     "var2",
                                     "paralog")]
        flag <- 0
      }
    }
    
    if (flag == 1){
      # get representative orthoID that has m VAR1 for each supertaxon
      mOrthoID <- fullMdData[, c("geneID",
                                 "supertaxon",
                                 "var1",
                                 "mVar1",
                                 "orthoID")]
      mOrthoID <- subset(mOrthoID,
                         mOrthoID$var1 == mOrthoID$mVar1)
      colnames(mOrthoID) <- c("geneID",
                              "supertaxon",
                              "var1",
                              "mVar1",
                              "orthoID")
      mOrthoID <- mOrthoID[!is.na(mOrthoID$orthoID), ]
      mOrthoID <- mOrthoID[, c("geneID", "supertaxon", "orthoID")]
      mOrthoID <- mOrthoID[!duplicated(mOrthoID[, 1:2]), ]
      
      # get data set for phyloprofile plotting (contains only supertaxa info)
      superDf <- subset(fullMdData, select = c("geneID",
                                               "supertaxon",
                                               "supertaxonID",
                                               "mVar1",
                                               "presSpec",
                                               "category",
                                               "mVar2",
                                               "paralog"))
      superDf$paralog <- 1
      superDf <- superDf[!duplicated(superDf), ]
      
      superDfExt <- merge(superDf, mOrthoID, by = c("geneID", "supertaxon"),
                          all.x = TRUE)
      superDfExt <- superDfExt[, c("geneID",
                                   "supertaxon",
                                   "supertaxonID",
                                   "mVar1",
                                   "presSpec",
                                   "category",
                                   "orthoID",
                                   "mVar2",
                                   "paralog")]
      
      # rename mVar to var
      names(superDfExt)[names(superDfExt) == "mVar1"] <- "var1"
      names(superDfExt)[names(superDfExt) == "mVar2"] <- "var2"
    }
    
    # print(superDfExt[superDfExt$geneID == "ampk_ACACB"
    # & superDfExt$supertaxon == "1001_Chordata",])
    # print("END2222")
    return(superDfExt)
  })
  
  # heatmap data input --------------------------------------------------------
  dataHeat <- reactive({
    
    # get all cutoffs
    percent_cutoff_min <- input$percent[1]
    percent_cutoff_max <- input$percent[2]
    var1_cutoff_min <- input$var1[1]
    var1_cutoff_max <- input$var1[2]
    var2_cutoff_min <- input$var2[1]
    var2_cutoff_max <- input$var2[2]
    
    # check input file
    filein <- input$main_input
    # if(input$demo == TRUE){
    if (input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      filein <- 1
    }
    if (is.null(filein)) return()
    
    dataHeat <- dataSupertaxa()
    
    # get selected supertaxon name
    split <- strsplit(as.character(input$in_select), "_")
    in_select <- as.character(split[[1]][1])
    
    ### replace insufficient values according to the thresholds by NA or 0
    dataHeat$presSpec[dataHeat$supertaxon != in_select
                      & dataHeat$presSpec < percent_cutoff_min] <- 0
    dataHeat$presSpec[dataHeat$supertaxon != in_select
                      & dataHeat$presSpec > percent_cutoff_max] <- 0
    
    dataHeat$presSpec[dataHeat$supertaxon != in_select
                      & dataHeat$var1 < var1_cutoff_min] <- 0
    dataHeat$presSpec[dataHeat$supertaxon != in_select
                      & dataHeat$var1 > var1_cutoff_max] <- 0
    
    if (input$var1_relation == "protein"){
      if (input$var2_relation == "protein"){
        # prot-prot: remove complete cell if one variable not sufficient
        dataHeat$presSpec[dataHeat$supertaxon != in_select
                          & dataHeat$var2 < var2_cutoff_min] <- 0
        dataHeat$presSpec[dataHeat$supertaxon != in_select
                          & dataHeat$var2 > var2_cutoff_max] <- 0
        dataHeat$var2[dataHeat$supertaxon != in_select
                      & dataHeat$var1 < var1_cutoff_min] <- NA
        dataHeat$var2[dataHeat$supertaxon != in_select
                      & dataHeat$var1 > var1_cutoff_max] <- NA
        dataHeat$var1[dataHeat$supertaxon != in_select
                      & dataHeat$var2 < var2_cutoff_min] <- NA
        dataHeat$var1[dataHeat$supertaxon != in_select
                      & dataHeat$var2 > var2_cutoff_max] <- NA
      } else {
        # prot-spec: var1 depend on var2
        dataHeat$presSpec[dataHeat$supertaxon != in_select
                          & dataHeat$var2 < var2_cutoff_min] <- 0
        dataHeat$presSpec[dataHeat$supertaxon != in_select
                          & dataHeat$var2 > var2_cutoff_max] <- 0
      }
    } else {
      if (input$var2_relation == "species"){
        # # spec-spec: remove var1 and var2 independently
        # dataHeat$presSpec[dataHeat$supertaxon != in_select & dataHeat$var1 < var1_cutoff_min] <- 0
        # dataHeat$presSpec[dataHeat$supertaxon != in_select & dataHeat$var1 > var1_cutoff_max] <- 0
      } else {
        # spec-prot: var2 depend on var1
        dataHeat$var2[dataHeat$supertaxon != in_select
                      & dataHeat$var1 < var1_cutoff_min] <- NA
        dataHeat$var2[dataHeat$supertaxon != in_select
                      & dataHeat$var1 > var1_cutoff_max] <- NA
      }
    }
    
    dataHeat$var1[dataHeat$supertaxon != in_select
                  & dataHeat$var1 < var1_cutoff_min] <- NA
    dataHeat$var1[dataHeat$supertaxon != in_select
                  & dataHeat$var1 > var1_cutoff_max] <- NA
    dataHeat$var2[dataHeat$supertaxon != in_select
                  & dataHeat$var2 < var2_cutoff_min] <- NA
    dataHeat$var2[dataHeat$supertaxon != in_select
                  & dataHeat$var2 > var2_cutoff_max] <- NA
    dataHeat <- droplevels(dataHeat)  # delete unused levels
    dataHeat$geneID <- as.factor(dataHeat$geneID)
    return(dataHeat)
  })
  
  # cluster heatmap data -----------------------------------------------------
  clusteredDataHeat <- reactive({
    dataHeat <- dataHeat()
    if (nrow(dataHeat) < 1) return()
    
    # dataframe for calculate distance matrix
    sub_data_heat <- subset(dataHeat, dataHeat$presSpec > 0)
    sub_data_heat <- sub_data_heat[, c("geneID", "supertaxon", "presSpec")]
    sub_data_heat <- sub_data_heat[!duplicated(sub_data_heat), ]
    
    wide_data <- spread(sub_data_heat, supertaxon, presSpec)
    dat <- wide_data[, 2:ncol(wide_data)]  # numerical columns
    rownames(dat) <- wide_data[, 1]
    dat[is.na(dat)] <- 0
    
    # get clustered gene ids
    clusteredGeneIDs <- clustered_gene_list(dat,
                                            input$dist_method,
                                            input$cluster_method)
    
    # sort original data according to clusteredGeneIDs
    dataHeat$geneID <- factor(dataHeat$geneID,
                              levels = clusteredGeneIDs)
    return(dataHeat)
  })
  
  # get info clicked point on main heatmap ------------------------------------
  mainpoint_info <- reactive({
    # check input
    if (v$doPlot == FALSE) return()
    
    # get selected supertaxon name
    taxa_list <- as.data.frame(read.table("data/taxonNamesReduced.txt",
                                          sep = "\t",
                                          header = T))
    rank_select <- input$rank_select
    rankName <- substr(rank_select, 4, nchar(rank_select))
    in_select <- {
      as.numeric(taxa_list$ncbiID[taxa_list$fullName == input$in_select])
    }
    
    dataHeat <- dataHeat()
    if (input$apply_cluster == TRUE){
      dataHeat <- clusteredDataHeat()
    }
    
    # get values
    if (is.null(input$plot_click$x)) return()
    else{
      # get cooridiate point
      if (input$x_axis == "genes"){
        corX <- round(input$plot_click$y);
        corY <- round(input$plot_click$x)
      } else {
        corX <- round(input$plot_click$x);
        corY <- round(input$plot_click$y)
      }
      
      # get geneID
      # genes <- as.matrix(dataHeat[dataHeat$supertaxonID == in_select
      #                             & !is.na(dataHeat$presSpec),])
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
        Percent <- max(na.omit(dataHeat$presSpec[dataHeat$geneID == geneID
                                                 & dataHeat$supertaxon == spec]))
      }
      var2 <- NA
      if (!is.na(dataHeat$var2[dataHeat$geneID == geneID
                               & dataHeat$supertaxon == spec][1])){
        var2 <- max(na.omit(dataHeat$var2[dataHeat$geneID == geneID
                                          & dataHeat$supertaxon == spec]))
      }
      
      # get ortholog ID
      orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID
                                  & dataHeat$supertaxon == spec]
      if (length(orthoID) > 1){
        orthoID <- paste0(orthoID[1], ",...")
      }
      # # get list of all geneID that have the same ortholog
      # geneMatch <- dataHeat$geneID[dataHeat$orthoID == toString(orthoID)]
      # geneMatch <- geneMatch[!is.na(geneMatch)]
      # # list of all available geneID
      # geneList <- preData()
      # geneList$geneID <- as.factor(geneList$geneID)
      # allGenes <- as.list(levels(geneList$geneID))
      # # get index of all matched genes (genes have the same ortholog)
      # pos <- which(allGenes %in% geneMatch)
      # pos <- paste(pos, collapse=",")
      
      # return info of clicked point
      if (is.na(as.numeric(Percent))) return()
      else{
        # info <- c(geneID,
        #           as.character(orthoID),
        #           as.character(spec),
        #           round(as.numeric(var1), 2),
        #           round(as.numeric(Percent), 2),
        #           round(as.numeric(var2), 2), pos)
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
  
  # var1 / var2 distribution data ---------------------------------------------
  distDf <- reactive({
    if (v$doPlot == FALSE) return()
    
    dataOrig <- get_long_matrix()
    if (ncol(dataOrig) < 4){
      colnames(dataOrig) <- c("geneID",
                              "ncbiID",
                              "orthoID")
      splitDt <- dataOrig[, c("orthoID")]
    } else if (ncol(dataOrig) < 5){
      colnames(dataOrig) <- c("geneID",
                              "ncbiID",
                              "orthoID",
                              "var1")
      splitDt <- dataOrig[, c("orthoID",
                              "var1")]
    } else {
      colnames(dataOrig) <- c("geneID",
                              "ncbiID",
                              "orthoID",
                              "var1",
                              "var2")
      splitDt <- dataOrig[, c("orthoID", "var1", "var2")]
    }
    
    splitDt$orthoID[splitDt$orthoID == "NA" | is.na(splitDt$orthoID)] <- NA
    splitDt <- splitDt[complete.cases(splitDt), ]
    
    if (length(levels(as.factor(splitDt$var2))) == 1){
      if (levels(as.factor(splitDt$var2)) == ""){
        splitDt$var2 <- 0
      }
    }
    
    # convert factor into numeric for "var1" & "var2" column
    if ("var1" %in% colnames(splitDt)){
      splitDt$var1 <- suppressWarnings(as.numeric(as.character(splitDt$var1)))
      # filter splitDt based on selected var1 cutoff
      splitDt <- splitDt[splitDt$var1 >= input$var1[1]
                         & splitDt$var1 <= input$var1[2], ]
    }
    if ("var2" %in% colnames(splitDt)){
      splitDt$var2 <- suppressWarnings(as.numeric(as.character(splitDt$var2)))
      # filter splitDt based on selected var2 cutoff
      splitDt <- splitDt[splitDt$var2 >= input$var2[1]
                         & splitDt$var2 <= input$var2[2], ]
    }
    
    # filter data base on customized plot (if chosen)
    if (input$dataset.distribution == "Customized data"){
      # get geneID and supertaxon name for splitDt
      allData <- get_data_filtered()
      splitDtName <- merge(splitDt, allData,
                           by = "orthoID",
                           all.x = TRUE)
      splitDtName$supertaxonMod <- {
        substr(splitDtName$supertaxon,
               6,
               nchar(as.character(splitDtName$supertaxon)))
      }
      splitDtName <- subset(splitDtName,
                            select = c(orthoID,
                                       var1.x,
                                       var2.y,
                                       supertaxonMod,
                                       geneID))
      colnames(splitDtName) <- c("orthoID",
                                 "var1",
                                 "var2",
                                 "supertaxonMod",
                                 "geneID")
      
      # filter
      if (input$in_taxa[1] == "all" & input$in_seq[1] != "all"){
        # select data from dataHeat for selected sequences only
        splitDt <- subset(splitDtName, geneID %in% input$in_seq)
      } else if (input$in_seq[1] == "all" & input$in_taxa[1] != "all"){
        # select data from dataHeat for selected taxa only
        splitDt <- subset(splitDtName, supertaxonMod %in% input$in_taxa)
      } else {
        # select data from dataHeat for selected sequences and taxa
        splitDt <- subset(splitDtName,
                          geneID %in% input$in_seq
                          & supertaxonMod %in% input$in_taxa)
      }
    }
    
    # return dt
    return(splitDt)
  })
  
  # calculate % present species for input file --------------------------------
  presSpecAllDt <- reactive({
    # open main input file
    mdData <- get_long_matrix()
    if (ncol(mdData) < 4){
      colnames(mdData) <- c("geneID",
                            "ncbiID",
                            "orthoID")
    } else if (ncol(mdData) < 5){
      colnames(mdData) <- c("geneID",
                            "ncbiID",
                            "orthoID",
                            "var1")
    } else {
      colnames(mdData) <- c("geneID",
                            "ncbiID",
                            "orthoID",
                            "var1",
                            "var2")
    }
    
    # count number of inparalogs
    paralogCount <- plyr::count(mdData, c("geneID", "ncbiID"))
    mdData <- merge(mdData, paralogCount, by = c("geneID", "ncbiID"))
    colnames(mdData)[ncol(mdData)] <- "paralog"
    
    # (3) GET SORTED TAXONOMY LIST (3) 
    taxa_list <- sortedtaxa_list()
    
    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxa_list, "supertaxon")
    
    # merge mdData, mdDatavar2 and taxa_list to get taxonomy info
    taxaMdData <- merge(mdData, taxa_list, by = "ncbiID")
    if ("var1" %in% colnames(taxaMdData)){
      taxaMdData$var1 <- {
        suppressWarnings(as.numeric(as.character(taxaMdData$var1)))
      }
    }
    if ("var2" %in% colnames(taxaMdData)){
      taxaMdData$var2 <- {
        suppressWarnings(as.numeric(as.character(taxaMdData$var2)))
      }
    }
    # calculate % present species
    finalPresSpecDt <- calc_pres_spec(taxaMdData, taxaCount)
    finalPresSpecDt
  })
  
  # print list of available taxa for customized plot --------------------------
  # (based on rank from rank_select_cus)
  taxa_select_cus <- reactive({
    rank_select_cus <- input$rank_select_cus
    
    if (length(rank_select_cus) == 0) return()
    else{
      ### load list of unsorted taxa
      Dt <- get_taxa_list(TRUE, subset_taxa())
      
      ### load list of taxon name
      nameList <- get_name_list()
      
      # get rank name from rank_select
      rankName <- substr(rank_select_cus, 4, nchar(rank_select_cus))
      choice <- as.data.frame
      choice <- rbind(Dt[rankName])
      colnames(choice) <- "ncbiID"
      choice <- merge(choice, nameList, by = "ncbiID", all = FALSE)
      return(choice)
    }
  })
  
  # get list of taxa based on selected taxa_select_cus ------------------------
  cus_taxaName <- reactive({
    
    taxa_select_cus <- input$taxa_select_cus
    rankName <- substr(input$rank_select_cus, 4, nchar(input$rank_select_cus))
    
    if (taxa_select_cus == "") return()
    
    # load list of unsorted taxa
    Dt <- get_taxa_list(TRUE, subset_taxa())
    
    # get ID of customized (super)taxon
    taxa_list <- get_name_list(FALSE, FALSE)
    superID <- taxa_list$ncbiID[taxa_list$fullName == taxa_select_cus
                                & taxa_list$rank %in% c(rankName, "norank")]
    
    # from that ID, get list of all taxa for main selected taxon
    mainRankName <- substr(input$rank_select, 4, nchar(input$rank_select))
    customizedtaxa_id <- {
      levels(as.factor(Dt[mainRankName][Dt[rankName] == superID, ]))
    }
    
    cus_taxaName <- {
      taxa_list$fullName[taxa_list$rank %in% c(mainRankName, "norank")
                         & taxa_list$ncbiID %in% customizedtaxa_id]
    }
    return(cus_taxaName)
  })
  
  # get info of a clicked point on selected plot ------------------------------
  # (also the same as get info from main plot)
  selectedpoint_info <- reactive({
    # check input
    if (vCt$doPlotCustom == FALSE) return()
    
    # get selected supertaxon name
    taxa_list <- get_name_list(FALSE, FALSE)
    rank_select <- input$rank_select
    rankName <- substr(rank_select, 4, nchar(rank_select))
    in_select <- {
      as.numeric(taxa_list$ncbiID[taxa_list$fullName == input$in_select])
    }
    
    dataHeat <- dataHeat()
    
    if (input$apply_cluster == TRUE){
      dataHeat <- clusteredDataHeat()
    }
    
    # get sub-dataframe of selected taxa and sequences
    dataHeat$supertaxonMod <- substr(dataHeat$supertaxon,
                                     6,
                                     nchar(as.character(dataHeat$supertaxon)))
    if (input$in_taxa[1] == "all" & input$in_seq[1] != "all"){
      # select data from dataHeat for selected sequences only
      dataHeat <- subset(dataHeat, geneID %in% input$in_seq)
    } else if (input$in_seq[1] == "all" & input$in_taxa[1] != "all"){
      # select data from dataHeat for selected taxa only
      dataHeat <- subset(dataHeat, supertaxonMod %in% input$in_taxa)
    } else {
      # select data from dataHeat for selected sequences and taxa
      dataHeat <- subset(dataHeat, geneID %in% input$in_seq
                         & supertaxonMod %in% input$in_taxa)
    }
    
    # drop all other supertaxon that are not in sub-dataframe
    dataHeat$supertaxon <- factor(dataHeat$supertaxon)
    dataHeat$geneID <- factor(dataHeat$geneID)
    
    # get values
    if (is.null(input$plot_click_selected$x)) return()
    else{
      # get cooridiate point
      if (input$x_axis_selected == "genes"){
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
      }
    }
  })
  
  # data for detailed plot ----------------------------------------------------
  detail_plotDt <- reactive({
    if (v$doPlot == FALSE) return()
    
    # GET INFO BASED ON CURRENT TAB
    if (input$tabs == "Main profile"){
      # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
      info <- mainpoint_info()
    } else if (input$tabs == "Customized profile"){
      info <- selectedpoint_info()
    }
    
    if (is.null(info)) return()
    else{
      ### get info for present taxa in selected supertaxon (1)
      plotTaxon <- info[3]
      plotGeneID <- info[1]
      fullDf <- get_data_filtered()
      selDf <- as.data.frame(fullDf[fullDf$geneID == plotGeneID
                                    & fullDf$supertaxon == plotTaxon, ])
      
      ### get all taxa of this supertaxon (2)
      allTaxaDf <- sortedtaxa_list()
      allTaxaDf <- allTaxaDf[allTaxaDf$supertaxon == plotTaxon, ]
      allTaxaDf <- subset(allTaxaDf, select = c("abbrName", "fullName"))
      
      ### merge (1) and (2) together
      joinedDf <- merge(selDf, allTaxaDf,
                        by = c("abbrName"),
                        all.y = TRUE)
      joinedDf <- subset(joinedDf,
                         select = c("abbrName",
                                    "fullName.y",
                                    "geneID",
                                    "orthoID",
                                    "var1",
                                    "var2"))
      names(joinedDf)[names(joinedDf) == "fullName.y"] <- "fullName"
      
      # replace var1/var2 as NA for all "NA orthologs"
      joinedDf$var1[is.na(joinedDf$orthoID)] <- NA
      joinedDf$var2[is.na(joinedDf$orthoID)] <- NA
      
      # remove NA orthologs if required
      if (input$detailed_remove_na == TRUE){
        joinedDf <- joinedDf[!is.na(joinedDf$orthoID), ]
      }
      
      ### return data for detailed plot
      return(joinedDf)
    }
  })
  
  # GET info when clicking on detailed plot -----------------------------------
  point_infoDetail <- reactive({
    selDf <- detail_plotDt()
    # selDf <- selDf[complete.cases(selDf),]
    selDf$orthoID <- as.character(selDf$orthoID)
    # allOrthoID <- sort(selDf$orthoID)
    
    ### get coordinates of plot_click_detail
    if (is.null(input$plot_click_detail$x)) return()
    else{
      corX <- round(input$plot_click_detail$y)
      corY <- round(input$plot_click_detail$x)
    }
    
    ### get pair of sequence IDs & var1
    seedID <- as.character(selDf$geneID[!is.na(selDf$geneID)][1])
    orthoID <- as.character(selDf$orthoID[corX])
    
    var1 <- as.list(selDf$var1[selDf$orthoID == orthoID])
    var1 <- as.character(var1[!is.na(var1)])
    var2 <- as.list(selDf$var2[selDf$orthoID == orthoID])
    var2 <- as.character(var2[!is.na(var2)])
    # ncbiID <- as.character(selDf$abbrName[selDf$orthoID==orthoID])
    ncbiID <- selDf[selDf$orthoID == orthoID, ]$abbrName
    ncbiID <- as.character(ncbiID[!is.na(ncbiID)][1])
    
    ### return info
    if (is.na(orthoID)){
      return(NULL)
    } else {
      if (orthoID != "NA"){
        info <- c(seedID, orthoID, var1, var2, ncbiID)
      }
    }
  })
  
  # get domain file/path ------------------------------------------------------
  getDomainFile <- reactive({
    # click info
    info <- point_infoDetail() # info = seedID, orthoID, var1
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    var1 <- as.character(info[3])
    

    # domain file
    if (input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      if (is.null(info)){
        fileDomain <- "noSelectHit"
        updateButton(session, "do_domain_plot", disabled = TRUE)
      } else {
        updateButton(session, "do_domain_plot", disabled = FALSE)
        if (input$demo_data == "demo"){
          fileDomain <- {
            suppressWarnings(paste0("https://github.com/BIONF/phyloprofile-data/blob/master/demo/domain_files/",
                                    group,
                                    ".domains?raw=true"))
          }
        } else if (input$demo_data == "ampk-tor") {
          fileDomain <- {
            suppressWarnings(paste0("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.domains_F"))
          }
        } 
      }
    } else {
      if (get_input_type() == "oma"){
        fileDomain <- "oma_input"
        updateButton(session, "do_domain_plot", disabled = FALSE)
      } else if (input$anno_choose == "from file"){
        fileDomain <- input$file_domain_input
        #print(fileDomain)
        if (is.null(fileDomain)){
          fileDomain <- "noFileInput"
        } else {
          if (is.null(info)){
            fileDomain <- "noSelectHit"
            updateButton(session, "do_domain_plot", disabled = TRUE)
          } else {
            updateButton(session, "do_domain_plot", disabled = FALSE)
            fileDomain <- fileDomain$datapath
          }
        }
      } else {
        if (is.null(info)){
          fileDomain <- "noSelectHit"
          updateButton(session, "do_domain_plot", disabled = TRUE)
        } else {
          # check file extension
          allExtension <- c("txt", "csv", "list", "domains", "architecture")
          flag <- 0
          for (i in 1:length(allExtension)){
            fileDomain <- paste0(input$domainPath,
                                 "/",
                                 group,
                                 ".",
                                 allExtension[i])
            if (file.exists(fileDomain) == TRUE){
              updateButton(session, "do_domain_plot", disabled = FALSE)
              flag <- 1
              break ()
            }
          }
          
          if (flag == 0){
            fileDomain <- "noFileInFolder"
            updateButton(session, "do_domain_plot", disabled = TRUE)
          }
        }
      }
    }
    return (fileDomain)
  })
  
  # gene age estimation -------------------------------------------------------
  gene_ageDf <- reactive({
    if (v$doPlot == FALSE) return()
    
    rankList <- c("family",
                  "class",
                  "phylum",
                  "kingdom",
                  "superkingdom",
                  "root")
    
    # get selected (super)taxon ID
    rank_select <- input$rank_select
    rankName <- substr(rank_select, 4, nchar(rank_select))
    
    taxa_list <- get_name_list(FALSE, FALSE)
    superID <- {
      as.numeric(taxa_list$ncbiID[taxa_list$fullName == input$in_select
                                  & taxa_list$rank == rankName])
    }
    
    # full non-duplicated taxonomy data
    Dt <- get_taxa_list(FALSE, subset_taxa())
    
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
    for (i in 1:nrow(subDt)){
      cat <- subDt[i, ] %in% sub_first_line
      cat[cat == FALSE] <- 0
      cat[cat == TRUE] <- 1
      cat <- paste0(cat, collapse = "")
      catDf[i, ] <- c(as.character(subDt[i, ]$abbrName), cat)
    }
    
    # get main input data
    mdData <- droplevels(get_data_filtered())
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
                            | is.na(mdDataExtended$orthoID), ]) > 0){
      mdDataExtended[mdDataExtended$orthoID == "NA"
                     | is.na(mdDataExtended$orthoID), ]$cat <- NA
    }
    
    mdDataExtended <- mdDataExtended[complete.cases(mdDataExtended), ]
    
    # filter by %specpres, var1, var2 ..
    mdDataExtended$cat[mdDataExtended$var1 < input$var1[1]] <- NA
    mdDataExtended$cat[mdDataExtended$var1 > input$var1[2]] <- NA
    mdDataExtended$cat[mdDataExtended$var2 < input$var2[1]] <- NA
    mdDataExtended$cat[mdDataExtended$var2 > input$var2[2]] <- NA
    mdDataExtended$cat[mdDataExtended$presSpec < input$percent[1]] <- NA
    mdDataExtended$cat[mdDataExtended$presSpec > input$percent[2]] <- NA
    
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
    gene_ageDf$age[gene_ageDf$cat == "0000011" | gene_ageDf$cat == "0000010"] <- {
      paste0("06_",
             as.character(taxa_list$fullName[taxa_list$ncbiID == sub_first_line$superkingdom
                                             & taxa_list$rank == "superkingdom"]))
    }
    gene_ageDf$age[gene_ageDf$cat == "0000111"] <- {
      paste0("05_",
             as.character(taxa_list$fullName[taxa_list$ncbiID == sub_first_line$kingdom
                                             & taxa_list$rank == "kingdom"]))
    }
    gene_ageDf$age[gene_ageDf$cat == "0001111"] <- {
      paste0("04_",
             as.character(taxa_list$fullName[taxa_list$ncbiID == sub_first_line$phylum
                                             & taxa_list$rank == "phylum"]))
    }
    gene_ageDf$age[gene_ageDf$cat == "0011111"] <- {
      paste0("03_",
             as.character(taxa_list$fullName[taxa_list$ncbiID == sub_first_line$class
                                             & taxa_list$rank == "class"]))
    }
    gene_ageDf$age[gene_ageDf$cat == "0111111"] <- {
      paste0("02_",
             as.character(taxa_list$fullName[taxa_list$ncbiID == sub_first_line$family
                                             & taxa_list$rank == "family"]))
    }
    gene_ageDf$age[gene_ageDf$cat == "1111111"] <- {
      paste0("01_",
             as.character(taxa_list$fullName[taxa_list$fullName == input$in_select
                                             & taxa_list$rank == rankName]))
    }
    
    # return gene_age data frame
    gene_ageDf <- gene_ageDf[, c("geneID", "cat", "age")]
    
    gene_ageDf$age[is.na(gene_ageDf$age)] <- "Undef"
    return(gene_ageDf)
  })
  
  gene_ageDfMod <- reactive({
    gene_ageDf <- gene_ageDf()
    countDf <- plyr::count(gene_ageDf, c("age"))
    countDf$percentage <- round(countDf$freq / sum(countDf$freq) * 100)
    countDf$pos <- cumsum(countDf$percentage) - (0.5 * countDf$percentage)
    return(countDf)
  })
  
  # render genAge.table based on clicked point on gene_agePlot ----------------
  selectedgene_age <- reactive({
    if (v$doPlot == FALSE) return()
    data <- gene_ageDf()
    
    # calculate the coordinate range for each age group
    rangeDf <- plyr::count(data, c("age"))
    
    rangeDf$percentage <- round(rangeDf$freq / sum(rangeDf$freq) * 100)
    rangeDf$rangeStart[1] <- 0
    rangeDf$rangeEnd[1] <- rangeDf$percentage[1]
    if (nrow(rangeDf) > 1){
      for (i in 2:nrow(rangeDf)){
        rangeDf$rangeStart[i] <- rangeDf$rangeEnd[i - 1] + 1
        rangeDf$rangeEnd[i] <- rangeDf$percentage[i] + rangeDf$rangeEnd[i - 1]
      }
    }
    
    # get list of selected age group
    if (is.null(input$plot_click_gene_age$x)) return()
    else{
      corX <- 100 - round(-input$plot_click_gene_age$x)
      selectAge <- {
        as.character(rangeDf[rangeDf$rangeStart <= corX
                             & rangeDf$rangeEnd >= corX, ]$age)
      }
      subData <- subset(data, age == selectAge)
      data <- data[data$age == selectAge, ]
    }
    
    # return list of genes
    geneList <- levels(as.factor(subData$geneID))
    geneList
  })
  
  cons_geneDf <- reactive({
    if (v$doPlot == FALSE) return()
    
    rankName <- substr(input$rank_select, 4, nchar(input$rank_select))
    
    # get ID list of chosen taxa
    taxa_list <- get_name_list(FALSE, FALSE)
    
    if ("none" %in% input$taxaCons) superID <- NA
    else{
      superID <- {
        taxa_list$ncbiID[taxa_list$fullName %in% input$taxaCons
                         & taxa_list$rank %in% c(rankName, "norank")]
      }
    }
    
    # get main input data
    mdData <- get_data_filtered()
    mdData <- mdData[, c("geneID",
                         "ncbiID",
                         "fullName",
                         "supertaxon",
                         "supertaxonID",
                         "rank",
                         "presSpec",
                         "mVar1",
                         "mVar2")]
    
    # filter by selecting taxa
    if (is.na(superID[1])) data <- NULL
    else{
      data <- subset(mdData, supertaxonID %in% superID
                     & presSpec >= input$percent_cons)
      # get supertaxa present in each geneID
      supertaxonCount <- {
        as.data.frame(plyr::count(data,
                                  c("geneID", "supertaxonID")))
      }
      # count number of supertaxa present in each geneID
      # and get only gene that contains all chosen taxa
      count <- as.data.frame(table(supertaxonCount$geneID))
      cons_gene <- subset(count, Freq == length(superID))
      cons_gene$Var1 <- factor(cons_gene$Var1)
      
      return(levels(cons_gene$Var1))
    }
  })
  
  # print list of available taxa for customized plot --------------------------
  # (based on rank from rank_select_cus)
  taxa_select_cons <- reactive({
    rank_select_cons <- input$rank_select_cons
    
    if (length(rank_select_cons) == 0) return()
    else{
      # load list of unsorted taxa
      Dt <- get_taxa_list(TRUE, subset_taxa())
      
      # load list of taxon name
      nameList <- get_name_list(TRUE, FALSE)
      
      # get rank name from rank_select
      rankName <- substr(rank_select_cons, 4, nchar(rank_select_cons))
      choice <- as.data.frame
      choice <- rbind(Dt[rankName])
      colnames(choice) <- "ncbiID"
      choice <- merge(choice, nameList, by = "ncbiID", all = FALSE)
      return(choice)
    }
  })
  
  # get list of taxa based on selected taxa_select_cus ------------------------
  consTaxaName <- reactive({
    
    taxa_select_cons <- input$taxa_select_cons
    rankName <- substr(input$rank_select_cons,
                       4,
                       nchar(input$rank_select_cons))
    
    if (taxa_select_cons == "") return()
    
    # load list of unsorted taxa
    Dt <- get_taxa_list(TRUE, subset_taxa())
    
    # get ID of customized (super)taxon
    taxa_list <- get_name_list(FALSE, FALSE)
    superID <- taxa_list$ncbiID[taxa_list$fullName == taxa_select_cons
                                & taxa_list$rank %in% c(rankName, "norank")]
    
    # from that ID, get list of all taxa for main selected taxon
    mainRankName <- substr(input$rank_select, 4, nchar(input$rank_select))
    constaxa_id <- {
      levels(as.factor(Dt[mainRankName][Dt[rankName] == superID, ]))
    }
    
    consTaxaName <- {
      taxa_list$fullName[taxa_list$rank %in% c(mainRankName, "norank")
                         & taxa_list$ncbiID %in% constaxa_id]
    }
    return(consTaxaName)
  })
  
  # cluster data --------------------------------------------------------------
  clusterDataDend <- reactive({
    if (v$doPlot == FALSE) return()
    # dataframe for calculate distance matrix
    dataHeat <- dataHeat()
    
    sub_data_heat <- subset(dataHeat, dataHeat$presSpec > 0)
    sub_data_heat <- sub_data_heat[, c("geneID", "supertaxon", "presSpec")]
    sub_data_heat <- sub_data_heat[!duplicated(sub_data_heat), ]
    
    wide_data <- spread(sub_data_heat, supertaxon, presSpec)
    dat <- wide_data[, 2:ncol(wide_data)]  # numerical columns
    rownames(dat) <- wide_data[, 1]
    dat[is.na(dat)] <- 0
    
    dd.col <- as.dendrogram(hclust(dist(dat, method = input$dist_method),
                                   method = input$cluster_method))
  })
  
  # render brushed_cluster.table based on clicked point on dendrogram plot ----
  brushed_clusterGene <- reactive({
    if (v$doPlot == FALSE) return()
    
    dd.col <- clusterDataDend()
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
  
  # filtered data for downloading (Main Profile ) -----------------------------
  download_data <- reactive({
    ### check input
    if (v$doPlot == FALSE) return()
    
    ### filtered data
    data_out <- get_data_filtered()
    data_out <- as.data.frame(data_out[data_out$presSpec > 0, ])
    data_out <- data_out[!is.na(data_out$geneID), ]
    
    data_out <- as.data.frame(data_out[data_out$presSpec >= input$percent[1], ])
    data_out <- as.data.frame(data_out[data_out$var1 >= input$var1[1]
                                       & data_out$var1 <= input$var1[2], ])
    
    if (!all(is.na(data_out$var2))){
      data_out <- as.data.frame(data_out[data_out$var2 >= input$var2[1]
                                         & data_out$var2 <= input$var2[2], ])
    } else {
      data_out$var2 <- 0
    }
    
    ### select only representative genes if chosen
    if (input$get_representative_main == TRUE){
      if (is.null(input$ref_var_main)) return()
      else{
        if (input$ref_var_main == input$var1_id){
          data_out_agg <- aggregate(as.numeric(data_out$var1),
                                    by = list(data_out$geneID, data_out$ncbiID),
                                    FUN = input$ref_type_main)
        } else if (input$ref_var_main == input$var2_id){
          data_out_agg <- aggregate(as.numeric(data_out$var2),
                                    by = list(data_out$geneID, data_out$ncbiID),
                                    FUN = input$ref_type_main)
        } else {
          data_out_agg <- data_out[data_out, c("geneID", "ncbiID", "var1")]
        }
        colnames(data_out_agg) <- c("geneID", "ncbiID", "var_best")
        
        data_out_representative <- merge(data_out, data_out_agg,
                                         by = c("geneID", "ncbiID"),
                                         all.x = TRUE)
        
        if (input$ref_var_main == input$var1_id){
          data_out <- {
            data_out_representative[data_out_representative$var1 == data_out_representative$var_best, ]
          }
        } else if (input$ref_var_main == input$var2_id){
          data_out <- {
            data_out_representative[data_out_representative$var2 == data_out_representative$var_best, ]
          }
        } else {
          data_out <- data_out
        }
        # used to select only one ortholog,
        # if there exist more than one "representative"
        data_out$dup <- paste0(data_out$geneID, "#", data_out$ncbiID)
        data_out <- data_out[!duplicated(c(data_out$dup)), ]
      }
    }
    
    # sub select columns of dataout
    data_out <- data_out[, c("geneID",
                             "orthoID",
                             "fullName",
                             "ncbiID",
                             "supertaxon",
                             "var1",
                             "var2",
                             "presSpec")] #,"numberSpec"
    data_out <- data_out[order(data_out$geneID, data_out$supertaxon), ]
    data_out <- data_out[complete.cases(data_out), ]
    
    data_out$geneID <- as.character(data_out$geneID)
    data_out$fullName <- as.character(data_out$fullName)
    data_out$ncbiID <- substr(data_out$ncbiID,
                              5,
                              nchar(as.character(data_out$ncbiID)))
    data_out$supertaxon <- substr(data_out$supertaxon,
                                  6,
                                  nchar(as.character(data_out$supertaxon)))
    data_out$var1 <- as.character(data_out$var1)
    data_out$var2 <- as.character(data_out$var2)
    # data_out$numberSpec <- as.numeric(data_out$numberSpec)
    data_out$presSpec <- as.numeric(data_out$presSpec)
    
    # rename columns
    names(data_out)[names(data_out) == "presSpec"] <- "%Spec"
    # names(data_out)[names(data_out)=="numberSpec"] <- "totalSpec"
    if (nchar(input$var1_id) > 0){
      names(data_out)[names(data_out) == "var1"] <- input$var1_id
    } else {
      data_out <- subset(data_out, select = -c(var1) )
    }
    if (nchar(input$var2_id) > 0){
      names(data_out)[names(data_out) == "var2"] <- input$var2_id
    } else {
      data_out <- subset(data_out, select = -c(var2) )
    }
    
    # return data for downloading
    data_out <- as.matrix(data_out)
    return(data_out)
  })
  
  # filtered data for downloading (Customized Profile) ------------------------
  download_custom_data <- reactive({
    # check input
    if (v$doPlot == FALSE) return()
    
    data <- as.data.frame(download_data())
    
    # get subset of data according to selected genes/taxa
    if (!is.null(input$in_seq) | !is.null(input$in_taxa)){
      if (input$in_seq[1] != "all" & input$in_taxa[1] == "all"){
        # select data for selected sequences only
        custom_data <- subset(data, geneID %in% input$in_seq)
      } else if (input$in_seq[1] == "all" & input$in_taxa[1] != "all"){
        # select data for selected taxa only
        custom_data <- subset(data, supertaxon %in% input$in_taxa)
      } else if (input$in_seq[1] != "all" & input$in_taxa[1] != "all") {
        # select data for selected sequences and taxa
        custom_data <- subset(data, geneID %in% input$in_seq
                              & supertaxon %in% input$in_taxa)
      } else {
        custom_data <- data
      }
    } else {
      custom_data <- data
    }
    # return data
    custom_data <- as.matrix(custom_data)
    custom_data
  })
  
  # Deciding which plots should be shown (Group Comparison) -------------------
  get_plots_gc <- reactive({
    gene <- as.character(input$selected_gene_gc)
    input$plot_gc
    if (is.null(significant_genes_gc)) return()
    else if (gene == "all"){
      get_plot_output_list(significant_genes_gc,
                           get_input_gc(),
                           input$interesting_features)
    }else{
      x <- {
        significant_genes_gc[significant_genes_gc$geneID == gene, ]
      }
      if (nrow(x) == 0) return()
      get_plot_output_list(x, get_input_gc(), input$interesting_features)
    }
  })
  
  # get list of taxa based on selected taxa_select_gc (Group Comparison) ------
  taxa_name_gc <- reactive({
    
    taxa_select_gc <- input$taxa_select_gc
    rank_name <- substr(input$rank_select_gc, 4, nchar(input$rank_select_gc))
    
    if (taxa_select_gc == "") return()
    
    # load list of unsorted taxa
    dt <- get_taxa_list(TRUE, subset_taxa())
    
    # get ID of customized (super)taxon
    taxa_list <- get_name_list(FALSE, FALSE)
    super_id <- {
      taxa_list$ncbiID[taxa_list$fullName == taxa_select_gc
                       & taxa_list$rank %in% c(rank_name, "norank")]
    }
    
    # from that ID, get list of all taxa for main selected taxon
    main_rank_name <- substr(input$rank_select, 4, nchar(input$rank_select))
    taxa_id_gc <- {
      levels(as.factor(dt[main_rank_name][dt[rank_name] == super_id, ]))
    }
    taxa_name_gc <- {
      taxa_list$fullName[taxa_list$rank %in% c(main_rank_name, "norank")
                         & taxa_list$ncbiID %in% taxa_id_gc]
    }
    return(taxa_name_gc)
  })
  
  })