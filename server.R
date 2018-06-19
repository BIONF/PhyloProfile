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

source("scripts/search_taxon_id.R")
source("scripts/parse_main_input.R")
source("scripts/parse_domain_input.R")

source("scripts/get_fasta_seqs.R")
source("scripts/download_filtered_main.R")
source("scripts/download_filtered_customized.R")

source("scripts/parse_phylotree.R")
source("scripts/select_taxon_rank.R")
source("scripts/create_profile_heatmap.R")

source("scripts/identify_core_gene.R")
source("scripts/analyze_distribution.R")
source("scripts/estimate_gene_age.R")
source("scripts/analyze_distribution.R")
source("scripts/cluster_profile.R")

options(shiny.maxRequestSize = 99 * 1024 ^ 2)  # size limit for input 99mb

shinyServer(function(input, output, session) {
  # Automatically stop a Shiny app when closing the browser tab
  # session$onSessionEnded(stopApp)
  session$allowReconnect(TRUE)
  
  # =========================== INITIAL CHECKING  =============================
  
  # * check for internet connection ---------------------------------------------
  observe({
    if (has_internet() == FALSE){
      toggleState("demo_data")
    }
  })
  
  output$no_internet_msg <- renderUI({
    if (has_internet() == FALSE){
      strong(em("Internet connection is required for using demo data!"),
             style = "color:red")
    } else {
      return()
    }
  })
  
  # * check for the existence of taxonomy files ---------------------------------
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
  
  # ======================== INPUT & SETTINGS TAB =============================
  
  # * check the validity of input file and render input_check.ui ----------------
  output$input_check.ui <- renderUI({
    filein <- input$main_input
    if(is.null(filein)) return()
    input_type <- check_input_vadility(filein) #get_input_type()
    
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
  
  # * render download link for Demo online files --------------------------------
  output$main_input_file.ui <- renderUI({
    # if(input$demo == TRUE){
    if (input$demo_data == "lca-micros"){
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
    if (input$demo_data == "lca-micros"){
      strong(a("Download demo domain files",
               href = "https://github.com/BIONF/phyloprofile-data/tree/master/demo/domain_files",
               target = "_blank"))
    } else if (input$demo_data == "ampk-tor"){
      strong(a("Download demo domain file",
               href = "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.domains_F",
               target = "_blank"))
    } else {
      # filein <- input$main_input
      # input_type <- check_input_vadility(filein) #get_input_type()
      
      if (input$anno_location == "from file"){
        fileInput("file_domain_input", "")
      } else {
        textInput("domainPath", "", "")
      }
    }
  })
  
  output$download_fastaDemo.ui <- renderUI({
    if (input$demo_data == "lca-micros"){
      strong(a("Download demo fasta file",
               href = "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/fasta_file/concatenatedSeq.fa",
               target = "_blank"))
    } else if (input$demo_data == "ampk-tor"){
      strong(a("Download demo fasta file",
               href = "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.extended.fa",
               target = "_blank"))
    }
  })
  
  # * render description for Demo data ------------------------------------------
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
  
  # * check OMA input -----------------------------------------------------------
  output$select_oma_type <- renderUI({
    filein <- input$main_input
    if(is.null(filein)) return()
    input_type <- check_input_vadility(filein) #get_input_type()
    
    if (input_type == "oma"){
      # Options to select the OMA type to generate the output
      selectInput("selected_oma_type", label = "Select type of OMA orthologs:",
                  
                  choices = list("PAIR", "HOG", "OG"),
                  selected = "PAIR")
    } else {
      return()
    }
  })
  
  output$button_oma <- renderUI({
    filein <- input$main_input
    if(is.null(filein)) return()
    input_type <- check_input_vadility(filein) #get_input_type()
    
    if (input_type == "oma"){
      shinyBS::bsButton("get_data_oma", "Get data")
    }
  })
  
  # * render link for download OMA files ----------------------------------------
  output$oma_download <- renderUI({
    filein <- input$main_input
    if(is.null(filein)) return()
    input_type <- check_input_vadility(filein) #get_input_type()
    
    if (input_type == "oma"){
      downloadButton("download_files_oma", "Download")
    } 
  })
  
  output$download_files_oma <- downloadHandler(
    filenname <- function(){
      "oma_data_to_phyloprofile_input.zip"
    },
    content <- function(file){
      write.table(get_main_input(), "long.txt",
                  sep = "\t",
                  row.names = FALSE,
                  col.names = TRUE,
                  quote = FALSE)
      
      write.table(long_to_fasta(get_main_input()), "fasta.txt",
                  sep = "\t",
                  row.names = FALSE,
                  col.names = FALSE,
                  quote = FALSE)
      
      write.table(get_domain_information (), "domain.txt",
                  sep = "\t",
                  row.names = FALSE,
                  col.names = FALSE,
                  quote = FALSE)
      
      zip(zipfile = file,
          files = c("long.txt", "domain.txt", "fasta.txt")) 
    },
    contentType = "application/zip"
  )
  
  # * render textinput for Variable 1 & 2 ---------------------------------------
  output$var1_id.ui <- renderUI({
    long_dataframe <- get_main_input()
    
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
    long_dataframe <- get_main_input()
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
  
  # * render 2. variable relationship according to demo data --------------------
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
  
  # * check the existance of the input concatenate fasta file -------------------
  output$concat_fasta.exist_check <- renderUI({
    if (is.null(input$concat_fasta)) return()
    else{
      f <- input$concat_fasta$datapath
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
  
  # * check the validity of input tree file and render checkNewick.ui -----------
  output$checkNewick.ui <- renderUI({
    filein <- input$inputTree
    if (is.null(filein)) return()
    
    check_newick <- check_newick(filein, input$main_input, subset_taxa())
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
  
  # * reset profile plot colors -------------------------------------------------
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
  
  # * render list of taxonomy ranks ------------------------------------------
  output$rank_select <- renderUI({
    # if(input$demo == TRUE){
    if (input$demo_data == "lca-micros"){
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
  
  # * render list of (super)taxa ---------------------------------------------
  output$select <- renderUI({
    choice <- alltaxa_list()
    choice$fullName <- as.factor(choice$fullName)
    
    if (input$demo_data == "lca-micros"){
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
  
  # * enable "PLOT" button ------------------------------------------------------
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
  
  # * move to main tab when "PLOT" button has been clicked ----------------------
  observe({
    # use tabsetPanel "id" argument to change tabs
    if (input$do > 0) {
      updateTabsetPanel(session, "tabs", selected = "Main profile")
    }
  })
  
  # * disable main input, genelist input and demo data checkbox --------------
  observe({
    if (input$do > 0) {
      toggleState("main_input")
      toggleState("gene_list_selected")
      toggleState("demo_data")
    }
  })
  
  # * update var2_aggregate_by to mean if using demo lca-micros data ---------
  observe({
    if (input$demo_data == "lca-micros") {
      ### update var2_aggregate_by to mean
      updateSelectInput(session, "var2_aggregate_by",
                        choices = list("Max" = "max",
                                       "Min" = "min",
                                       "Mean" = "mean",
                                       "Median" = "median"),
                        selected = "mean")
    }
  })
  
  # =========================== RENDER FILTER SLIDEBARS =======================
  
  # * render filter slidebars for Main plot -------------------------------------
  output$var1_cutoff.ui <- renderUI({
    create_slider_cutoff("var1", paste(input$var1_id, "cutoff:"), 0.0, 1.0, input$var1_id)
  })
  
  output$var2_cutoff.ui <- renderUI({
    create_slider_cutoff("var2", paste(input$var2_id, "cutoff:"), 0.0, 1.0, input$var2_id)
  })
  
  output$percent_cutoff.ui <- renderUI({
    create_slider_cutoff("percent", "% of present taxa:", 0.0, 1.0, "percent")
  })
  
  # * render filter slidebars for Customized plot -------------------------------
  output$var1_filter.ui <- renderUI({
    create_slider_cutoff("var1cus", paste(input$var1_id, "cutoff:"), input$var1[1], input$var1[2], input$var1_id)
  })
  
  output$var2_filter.ui <- renderUI({
    create_slider_cutoff("var2cus", paste(input$var2_id, "cutoff:"), input$var2[1], input$var2[2], input$var2_id)
  })
  
  output$percent_filter.ui <- renderUI({
    create_slider_cutoff("percent2", "% of present taxa:", input$percent[1], input$percent[2], "percent")
  })
  
  # * render filter slidebars for Distribution plot -----------------------------
  output$var1_dist.ui <- renderUI({
    create_slider_cutoff("var1_dist", paste(input$var1_id, "cutoff:"), input$var1[1], input$var1[2], input$var1_id)
  })
  
  output$var2_dist.ui <- renderUI({
    create_slider_cutoff("var2_dist", paste(input$var2_id, "cutoff:"), input$var2[1], input$var2[2], input$var2_id)
  })
  
  output$percent_dist.ui <- renderUI({
    create_slider_cutoff("percent_dist", "% of present taxa:", input$percent[1], input$percent[2], "percent")
  })
  
  # * render filter slidebars for Gene age estimation plot ----------------------
  output$var1_age.ui <- renderUI({
    create_slider_cutoff("var1_age", paste(input$var1_id, "cutoff:"), input$var1[1], input$var1[2], input$var1_id)
  })
  
  output$var2_age.ui <- renderUI({
    create_slider_cutoff("var2_age", paste(input$var2_id, "cutoff:"), input$var2[1], input$var2[2], input$var2_id)
  })
  
  output$percent_age.ui <- renderUI({
    create_slider_cutoff("percent_age", "% of present taxa:", input$percent[1], input$percent[2], "percent")
  })
  
  # * render filter slidebars for Core gene finding function --------------------
  output$var1_core.ui <- renderUI({
    create_slider_cutoff("var1_core", paste(input$var1_id, "cutoff:"), input$var1[1], input$var1[2], input$var1_id)
  })
  
  output$var2_core.ui <- renderUI({
    create_slider_cutoff("var2_core", paste(input$var2_id, "cutoff:"), input$var2[1], input$var2[2], input$var2_id)
  })
  
  output$percent_core.ui <- renderUI({
    create_slider_cutoff("percent_core", "% of present taxa:", input$percent[1], input$percent[2], "percent")
  })
  
  # * update value for filter slidebars of Main Plot --------------------------
  # ** based on customized profile
  observe({
    new_var1 <- input$var1cus
    update_slider_cutoff(session, "var1", paste(input$var1_id, "cutoff:"), new_var1, input$var1_id)
  })
  
  observe({
    new_var2 <- input$var2cus
    update_slider_cutoff(session, "var2", paste(input$var2_id, "cutoff:"), new_var2, input$var2_id)
  })
  
  observe({
    new_percent <- input$percent2
    update_slider_cutoff(session, "percent", "% of present taxa:", new_percent, "percent")
  })
  
  # ** based on "Distribution analysis"
  observe({
    new_var1 <- input$var1_dist
    update_slider_cutoff(session, "var1", paste(input$var1_id, "cutoff:"), new_var1, input$var1_id)
  })
  
  observe({
    new_var2 <- input$var2_dist
    update_slider_cutoff(session, "var2", paste(input$var2_id, "cutoff:"), new_var2, input$var2_id)
  })
  
  observe({
    new_percent <- input$percent_dist
    update_slider_cutoff(session, "percent", "% of present taxa:", new_percent, "percent")
  })
  
  # ** based on "Gene age estimation"
  observe({
    new_var1 <- input$var1_age
    update_slider_cutoff(session, "var1", paste(input$var1_id, "cutoff:"), new_var1, input$var1_id)
  })
  
  observe({
    new_var2 <- input$var2_age
    update_slider_cutoff(session, "var2", paste(input$var2_id, "cutoff:"), new_var2, input$var2_id)
  })
  
  observe({
    new_percent <- input$percent_age
    update_slider_cutoff(session, "percent", "% of present taxa:", new_percent, "percent")
  })
  
  # * reset cutoffs of Main plot ------------------------------------------------
  observeEvent(input$reset_main, {
    shinyjs::reset("var1")
    shinyjs::reset("var2")
    shinyjs::reset("percent")
  })
  
  # * reset cutoffs of Customized plot ------------------------------------------
  observeEvent(input$reset_selected, {
    shinyjs::reset("var1")
    shinyjs::reset("var2")
    shinyjs::reset("percent")
  })
  
  # ========================= PARSING UNKNOWN TAXA ============================
  
  # * get list of "unknown" taxa in main input ---------------------
  unkTaxa <- reactive({
    long_dataframe <- get_main_input ()
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
    # return input taxa
    return(inputTaxa)
  })
  
  # * check the status of unkTaxa ---------------------------------------------
  output$unk_taxa_status <- reactive({
    unkTaxa <- unkTaxa()
    length(unkTaxa) > 0
  })
  outputOptions(output, "unk_taxa_status", suspendWhenHidden = FALSE)
  
  # * render list of unkTaxa --------------------------------------------------
  output$unk_taxa_full <- renderDataTable(option = list(searching = FALSE,
                                                        pageLength = 10), {
                                                          if (length(unkTaxa()) > 0){
                                                            tb <- as.data.frame(unkTaxa())
                                                            names(tb)[1] <- "New taxon"
                                                            tb
                                                          }
                                                        })
  
  # * update the form for adding new taxa -------------------------------------
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
  
  # * check if data is loaded and "parse" button is clicked and confirmed -----
  v1 <- reactiveValues(parse = FALSE)
  observeEvent(input$but_parse, {
    toggleModal(session, "parse_confirm", toggle = "close")
    v1$parse <- input$but_parse
    updateButton(session, "but_parse", disabled = TRUE)
    toggleState("new_taxa_ask")
    toggleState("main_input")
  })
  
  # * create rankList.txt, idList.txt, taxonNamesReduced.txt from input file ----
  # * and also create a full taxonomyMatrix.txt for sorting taxa
  invalidID <- {
    reactiveValues(df = data.frame("Invalid NCBI ID(s)" = as.character(),
                                   stringsAsFactors = F))
  }
  
  observe({
    filein <- input$main_input
    if(is.null(filein)) return()
    input_type <- check_input_vadility(filein) #get_input_type()
    
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
              
              # create list of IDs for this taxon
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
              
              # append info into rank and ids
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
  
  # * output invalid NCBI ID --------------------------------------------------
  output$invalidID.output <- renderTable({
    if (nrow(invalidID$df) < 1) return()
    else{
      outDf <- invalidID$df
      colnames(outDf) <- c("Invalid NCBI ID(s)")
      return(outDf)
    }
  })
  
  # * render final msg after taxon parsing ------------------------------------
  output$end_parsing_msg <- renderUI({
    if (nrow(invalidID$df) < 1) {
      strong(h4("PLEASE RELOAD THIS TOOL AFTER ADDING NEW TAXA!!!"),
             style = "color:red")
    }else{
      HTML('<p><strong><span style="color: #e12525;">SOME INVALID TAXON IDs HAVE BEEN FOUND!!</span><br>Please check the validity of the following IDs in
           <a target="_blank" href="https://www.ncbi.nlm.nih.gov/taxonomy">NCBI taxonomy database</a>!</strong></p>')
    }
    })
  
  # * close parsing windows ---------------------------------------------------
  observeEvent(input$new_done, {
    toggleModal(session, "add_taxa_windows", toggle = "close")
    write.table(newTaxa$Df, "data/newTaxa.txt",
                sep = "\t",
                eol = "\n",
                row.names = FALSE,
                quote = FALSE)
  })
  
  
  
  
  # ====================== PROCESSING INPUT DATA ===============================
  
  # * check if data is loaded and "plot" button is clicked --------------------
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
  
  # * check if "no ordering gene IDs" has been checked ------------------------
  output$apply_cluster_check.ui <- renderUI({
    if (input$ordering == FALSE){
      HTML('<p><em>(Check "Ordering sequence IDs" check box in <strong>Input & settings tab</strong>&nbsp;to enable this function)</em></p>')
    }
  })
  # * to enable clustering 
  observe({
    if (input$ordering == FALSE){
      shinyjs::disable("apply_cluster")
    } else {
      shinyjs::enable("apply_cluster")
    }
  })
  
  # * get the type of the input file & return long format dataframe -----------
  get_main_input <- reactive({
    if(input$demo_data == "lca-micros"){
      long_dataframe <- create_long_matrix("lca-micros")
    } else if (input$demo_data == "ampk-tor"){
      long_dataframe <- create_long_matrix("ampk-tor")
    } else {
      filein <- input$main_input
      if(is.null(filein)) return()
      long_dataframe <- create_long_matrix(filein)
    }
    
    return(long_dataframe)
  })
  
  # * parse domain info into data frame -------------------------------------
  get_domain_information <- reactive({
    filein <- input$main_input
    if(is.null(filein)) return()
    input_type <- check_input_vadility(filein)
    main_input <- get_main_input()
    
    domain_df <- parse_domain_input(main_input,
                                    input_type,
                                    input$demo_data,
                                    input$anno_location,
                                    input$file_domain_input,
                                    input$domainPath,
                                    session,
                                    datapath)
    return(domain_df)
  })
  
  # * parse fasta sequences as a data frame -----------------------------------
  # get_fasta_information <- reactive({
  #   
  # })
  
  # * get ID list of input taxa from main input -----------------------------
  subset_taxa <- reactive({
    if (input$demo_data == "lca-micros" |
        input$demo_data == "ampk-tor" |
        length(unkTaxa()) == 0){
      long_dataframe <- get_main_input()
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
  
  # * get NAME list of all (super)taxa -----------------------------------------------
  alltaxa_list <- reactive({
    
    filein <- input$main_input
    
    if (is.null(filein) & input$demo_data == "none") return()
    
    rank_select <- input$rank_select
    
    if (rank_select == "") return()
    if (length(unkTaxa()) > 0) return()
    
    # load list of unsorted taxa
    Dt <- get_taxa_list(TRUE, subset_taxa())
    
    # load list of taxon name
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
  
  # * sort input taxa ---------------------------------------------------------
  sortedtaxa_list <- reactive({
    if (v$doPlot == FALSE) return()
    
    # FIRST, GET REPRESENTATIVE TAXON
    
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
    
    # THEN, SORT TAXON LIST BASED ON HCLUST TREE
    
    # prepare Df for calculating distance matrix
    distDf <- subset(Dt, select = -c(ncbiID, fullName))
    row.names(distDf) <- distDf$abbrName
    distDf <- distDf[, -1]
    
    # get sorted taxon IDs & sort full taxonomy info dataframe
    treeIn <- input$inputTree
    
    if (is.null(treeIn)){
      taxaTree <- create_rooted_tree(distDf, as.character(repTaxon$abbrName))
    } else {
      taxaTree <- read.tree(file = treeIn$datapath)
    }
    taxonList <- sort_taxa_from_tree(taxaTree)
    sortedDt <- Dt[match(taxonList, Dt$abbrName), ]
    
    # subset to get list of input taxa only
    inputTaxa <- subset_taxa()
    sortedDt <- subset(sortedDt, abbrName %in% inputTaxa)
    
    # get only taxonIDs list of selected rank and rename columns
    sortedOut <- subset(sortedDt, select = c("abbrName",
                                             "ncbiID",
                                             "fullName",
                                             as.character(rankName)))
    colnames(sortedOut) <- c("abbrName",
                             "species",
                             "fullName",
                             "ncbiID")
    
    # add name of supertaxa into sortedOut list
    sortedOut <- merge(sortedOut, nameList,
                       by = "ncbiID",
                       all.x = TRUE,
                       sort = FALSE)
    sortedOut$species <- as.character(sortedOut$species)
    
    # add order_prefix to supertaxon name
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
  
  # * get subset data (default: first 30 genes) for plotting ------------------
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
    
    long_dataframe <- get_main_input()
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
  
  # * creating main dataframe for subset taxa (in species/strain level) -------
  # * get (super)taxa names (1)
  # * calculate percentage of presence (2),
  # * max/min/mean/median VAR1 (3) and VAR2 (4)
  get_data_filtered <- reactive({
    mdData <- preData()
    
    # count number of inparalogs
    paralogCount <- plyr::count(mdData, c("geneID", "ncbiID"))
    mdData <- merge(mdData, paralogCount, by = c("geneID", "ncbiID"))
    colnames(mdData)[ncol(mdData)] <- "paralog"
    
    # (1) GET SORTED TAXONOMY LIST (1)
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
    
    # (2) calculate PERCENTAGE of PRESENT SPECIES (2)
    finalPresSpecDt <- calc_pres_spec(taxaMdData, taxaCount)
    
    # (3) calculate max/min/mean/median VAR1 for every supertaxon of each gene (3)
    # remove NA rows from taxaMdData
    taxaMdDataNoNA <- taxaMdData[!is.na(taxaMdData$var1), ]
    # calculate m var1
    mVar1Dt <- aggregate(taxaMdDataNoNA[, "var1"],
                         list(taxaMdDataNoNA$supertaxon,
                              taxaMdDataNoNA$geneID),
                         FUN = input$var1_aggregate_by)
    colnames(mVar1Dt) <- c("supertaxon", "geneID", "mVar1")
    
    # (4) calculate max/min/mean/median VAR2 for each super taxon (4)
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
    
    # (3+4) & join mVar2 together with mVar1 scores into one df (3+4)
    scoreDf <- merge(mVar1Dt,
                     mVar2Dt,
                     by = c("supertaxon", "geneID"),
                     all = TRUE)
    
    # (2+3+4) add presSpec and mVar1 into taxaMdData (2+3+4)
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
  
  # * reduce data from species/strain level to supertaxon (e.g. phylum) level -
  # * This data set contain only supertaxa
  # * and their value (%present, mVar1 & mVar2) for each gene
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
  
  # * heatmap data input --------------------------------------------------------
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
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor"){
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
  
  # * clustered heatmap data -----------------------------------------------------
  clusteredDataHeat <- reactive({
    dataHeat <- dataHeat()
    if (nrow(dataHeat) < 1) return()
    dat <- get_data_clustering(dataHeat, input$dist_method)
    
    # get clustered gene ids
    clusteredGeneIDs <- clustered_gene_list(dat,
                                            input$dist_method,
                                            input$cluster_method)
    
    # sort original data according to clusteredGeneIDs
    dataHeat$geneID <- factor(dataHeat$geneID,
                              levels = clusteredGeneIDs)
    return(dataHeat)
  })
  
  # =========================== MAIN PROFILE TAB ==============================
  
  # * get total number of genes -----------------------------------------------
  output$total_gene_number.ui <- renderUI({
    geneList <- preData()
    geneList$geneID <- as.factor(geneList$geneID)
    out <- as.list(levels(geneList$geneID))
    if (length(out) > 0){
      # em(paste0("Total number of genes:  ",length(out)))
      strong(paste0("Total number of genes:  ", length(out)))
    }
  })
  
  # * get list of taxa for highlighting ---------------------------------------
  output$highlight_taxon_ui <- renderUI({
    choice <- alltaxa_list()
    choice$fullName <- as.factor(choice$fullName)
    
    out <- as.list(levels(choice$fullName))
    out <- append("none", out)
    
    selectInput("taxon_highlight", "Select (super)taxon to highlight:",
                out, selected = out[1])
  })
  
  # * update highlight_taxon_ui based on double clicked dot ---------------------
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
  
  # * get list of genes for highlighting ----------------------------------------
  output$highlight_gene_ui <- renderUI({
    geneList <- dataHeat()
    geneList$geneID <- as.factor(geneList$geneID)
    
    out <- as.list(levels(geneList$geneID))
    out <- append("none", out)
    
    selectInput("gene_highlight", "Highlight:", out, selected = out[1])
  })
  
  # * update highlight_gene_ui based on double clicked dot ----------------------
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
  
  # * reset configuration windows of Main plot ----------------------------------
  observeEvent(input$reset_main_config, {
    shinyjs::reset("x_size")
    shinyjs::reset("y_size")
    shinyjs::reset("legend_size")
    shinyjs::reset("x_angle")
    shinyjs::reset("dot_zoom")
  })
  
  # * close configuration windows of Main plot ----------------------------------
  observeEvent(input$apply_main_config, {
    toggleModal(session, "main_plot_config_bs", toggle = "close")
  })
  
  # * parameters for the main profile plot --------------------------------------
  get_parameter_input_main <- reactive({
    input_para <- list( "x_axis" = input$x_axis,
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
                        # "taxon_highlight" = input$taxon_highlight,
                        # "apply_cluster" = input$apply_cluster,
                        # "rank_select" = input$rank_select,
                        # "gene_highlight" = input$gene_highlight,
                        "guideline" = 1)
    
    return(input_para)
  })
  
  # * render dot size to dot_size_info ------------------------------------------
  output$dot_size_info <- renderUI({
    if (v$doPlot == FALSE) return()
    
    dataHeat <- dataHeat()
    dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
    presentVl <- dataHeat$presSpec[!is.na(dataHeat$presSpec)]
    
    minDot <- (floor(min(presentVl) * 10) / 10 * 5) * (1 + input$dot_zoom)
    maxDot <- (floor(max(presentVl) * 10) / 10 * 5) * (1 + input$dot_zoom)
    
    em(paste0("current point's size: ", minDot, " - ", maxDot))
  })
  
  # * plot main profile -------------------------------------------------
  output$mainPlot <- renderPlot({
    data_heat <- data_main_plot(dataHeat())
    # cluster dataHeat (if selected)
    if (input$apply_cluster == TRUE) {
      data_heat <- data_main_plot(clusteredDataHeat())
    }
    
    profile_plot(data_heat, get_parameter_input_main(), input$taxon_highlight, input$rank_select, input$gene_highlight)
  })
  
  output$plot.ui <- renderUI({
    # show beschreibung file if no plot present
    if (v$doPlot == FALSE){
      return()
    } else{
      # if auto_update is NOT selected, use update_btn to trigger plot changing
      if (input$auto_update == FALSE){
        input$update_btn
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
  
  # * download main profile --------------------------------------------------------
  output$plot_download <- downloadHandler(
    filename = function() {
      c("plot.pdf")
    },
    content = function(file) {
      ggsave(file, plot = main_plot(v, dataHeat(),
                                    clusteredDataHeat(),
                                    get_parameter_input_main ()),
             width = input$width * 0.056458333,
             height = input$height * 0.056458333,
             units = "cm",
             dpi = 300,
             device = "pdf",
             limitsize = FALSE)
    }
  )
  
  # * get info of clicked point on main profile ------------------------------------
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
      
      # return info of clicked point
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
  
  # ======================== CUSTOMIZED PROFILE TAB ===========================
  
  # * get list of all sequence IDs for customized profile -----
  output$gene_in <- renderUI({
    filein <- input$main_input
    fileCustom <- input$custom_file
    
    # if(input$demo == TRUE){
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor"){
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
      else if (input$add_core_gene_custom_profile == TRUE){
        out <- core_geneDf()
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
      else if(input$add_gc_genes_custom_profile == TRUE){
        out <- significantGenesGroupCompairison$geneID
        if(length(out)>0){
          selectInput("in_seq","",out,selected=as.list(out),multiple=TRUE,selectize=FALSE)
        }
        else {
          selectInput("in_seq","",outAll,selected=outAll[1],multiple=TRUE,selectize=FALSE)
        }
        
      }
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
  
  # * render popup for selecting taxon rank and return list of subset taxa ----
  cus_taxaName <- callModule(select_taxon_rank,
                             "select_taxon_rank",
                             rank_select = reactive(input$rank_select),
                             subset_taxa = subset_taxa)
  
  # * get list of all taxa for customized profile ----------------------------------
  output$taxa_in <- renderUI({
    filein <- input$main_input
    # if(input$demo == TRUE){
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor"){
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
  
  # * check if all genes and all species are selected ---------------------------
  output$same_profile <- reactive({
    if (v$doPlot == FALSE) return(FALSE)
    if (length(input$in_seq[1]) == 0) return(FALSE)
    else{
      if (input$in_seq[1] == "all" & input$in_taxa[1] == "all") return(TRUE)
    }
  })
  outputOptions(output, "same_profile", suspendWhenHidden = FALSE)
  
  # * change label of plot_custom button if auto_update_selected is unchecked ----
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
  
  # * check if button (custom)PLOT is clicked ------------------------------------------------
  vCt <- reactiveValues(doPlotCustom = FALSE)
  observeEvent(input$plot_custom, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    vCt$doPlotCustom <- input$plot_custom
    filein <- input$main_input
    # if(input$demo == TRUE){ filein = 1 }
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor") filein <- 1
    if (is.null(filein)) vCt$doPlotCustom <- FALSE
  })
  
  # * reset configuration windows of Customized plot ----------------------------
  observeEvent(input$reset_selected_config, {
    shinyjs::reset("x_size_select")
    shinyjs::reset("y_size_select")
    shinyjs::reset("legend_size_select")
    shinyjs::reset("x_angle_select")
    shinyjs::reset("dot_zoom_select")
  })
  
  # ** close configuration windows of Customized plot ----------------------------
  observeEvent(input$apply_selected_config, {
    toggleModal(session, "selected_plot_config_bs", toggle = "close")
  })
  
  # * parameters for the customized profile plot --------------------------------
  get_parameter_input_customized <- reactive({
    input_para <- list( "x_axis" = input$x_axis_selected,
                        "var1_id" = input$var1_id,
                        "var2_id"  = input$var2_id,
                        "low_color_var1" =  input$low_color_var1,
                        "high_color_var1" = input$high_color_var1,
                        "low_color_var2" = input$low_color_var2,
                        "high_color_var2" = input$high_color_var2,
                        "para_color" = input$para_color,
                        "x_size" = input$x_size_select,
                        "y_size" = input$y_size_select,
                        "legend_size" = input$legend_size_select,
                        "main_legend" = input$selected_legend,
                        "dot_zoom" = input$dot_zoom_select,
                        "x_angle" = input$x_angle_select,
                        "guideline" = 0)
    
    return(input_para)
  })
  
  # * plot customized profile -------------------------------------------
  output$selected_plot <- renderPlot({
    if (vCt$doPlotCustom == FALSE) return()
    if (input$in_seq[1] == "all" & input$in_taxa[1] == "all") return()
    
    data_heat <- data_customized_plot(dataHeat(), input$in_taxa, input$in_seq)
    # cluster dataHeat (if selected)
    if (input$apply_cluster == TRUE) {
      data_heat <- data_customized_plot(clusteredDataHeat(), input$in_taxa, input$in_seq)
    }
    
    profile_plot(data_heat, get_parameter_input_customized(), "none", input$rank_select, "none")
  })
  
  output$selected_plot.ui <- renderUI({
    if (is.null(input$in_seq[1]) | is.null(input$in_taxa[1]))  return()
    else if (input$in_seq[1] == "all" & input$in_taxa[1] == "all") return()
    else{
      if (input$auto_update_selected == FALSE){
        input$plot_custom
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
  
  # * download customized plot ----------------------------------------------------
  output$selected_download <- downloadHandler(
    filename = function() {
      c("selected_plot.pdf")
    },
    content = function(file) {
      ggsave(file, plot = selected_plot(vCt,
                                        dataHeat(),
                                        clusteredDataHeat(),
                                        get_parameter_input_customized()),
             width = input$selected_width * 0.056458333,
             height = input$selected_height * 0.056458333,
             units = "cm", dpi = 300, device = "pdf", limitsize = FALSE)
    }
  )
  
  # * get info of clicked point on customized profile ------------------------------
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
  
  # ============================== POINT INFO =================================
  
  # * get status of point_info for activating Detailed Plot button ---------------
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
  
  # * show info into "point's info" box -----------------------------------------
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
  
  # ============================= DETAILED PLOT ===============================
  
  # * data for detailed plot ----------------------------------------------------
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
  
  # * plot detailed bar chart ---------------------------------------------------
  output$detail_plot <- renderPlot({
    p <- detail_plot(v, detail_plotDt(), input$detailed_text,
                     input$var1_id, input$var2_id)
    p
  })
  
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
  
  # * download detailed plot ----------------------------------------------------
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
  
  # * get info when clicking on detailed plot -----------------------------------
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
    if (length(var2) == 0) var2 = "NA"
    # ncbiID <- as.character(selDf$abbrName[selDf$orthoID==orthoID])
    ncbiID <- selDf[selDf$orthoID == orthoID, ]$abbrName
    ncbiID <- as.character(ncbiID[!is.na(ncbiID)][1])
    
    ### return info
    if (is.na(orthoID)){
      return(NULL)
    } else {
      if (orthoID != "NA"){
        info <- c(seedID, orthoID, var1, var2, ncbiID)
        return(info)
      }
    }
  })
  
  # * show info when clicking on detailed plot ----------------------------------
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
  
  # * render FASTA sequence ------------------------------------------------------------
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
      fastaOut <- get_fasta_seqs(seqDf, input$main_input, input$demo_data,
                                 input$input_type, input$concat_fasta,
                                 input$path,
                                 input$dir_format,
                                 input$file_ext,
                                 input$id_format,
                                 get_main_input())
      return(paste(fastaOut[1]))
    }
  })
  
  # ======================== FEATURE ARCHITECTURE PLOT ========================
  
  # * get domain file/path ------------------------------------------------------
  getDomainFile <- reactive({
    # click info
    info <- point_infoDetail() # info = seedID, orthoID, var1
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    var1 <- as.character(info[3])
    
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor"){
      if (is.null(info)){
        fileDomain <- "noSelectHit"
        updateButton(session, "do_domain_plot", disabled = TRUE)
      } else {
        updateButton(session, "do_domain_plot", disabled = FALSE)
        if (input$demo_data == "lca-micros"){
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
      filein <- input$main_input
      if (check_input_vadility(filein) == "oma"){
        # if (get_input_type() == "oma"){
        fileDomain <- "oma_input"
        updateButton(session, "do_domain_plot", disabled = FALSE)
      } else if (input$anno_location == "from file"){
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
  
  # * check domain file ---------------------------------------------------------
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
  
  # * check if PLOT button is clicked -------------------------------------------------------------
  v3 <- reactiveValues(doPlot3 = FALSE)
  observeEvent(input$do_domain_plot, {
    v3$doPlot3 <- input$do_domain_plot
    filein <- input$main_input
    # if(input$demo == TRUE){ filein = 1 }
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor") filein <- 1
    if (is.null(filein)) v3$doPlot3 <- FALSE
  })
  
  # * render domain architecture plot -------------------------------------------
  output$archi_plot <- renderPlot({
    g <- archi_plot(v3,
                    point_infoDetail(),
                    get_domain_information(),
                    input$concat_fasta,
                    input$label_archi_size,
                    input$title_archi_size)
    grid.draw(g)
  })
  
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
  
  # * download architecture plot ------------------------------------------------
  # *** something strange with archi_plot()
  output$archi_download <- downloadHandler(
    filename = function() {
      c("domains.pdf")
    },
    content = function(file) {
      g <- archi_plot(v3,
                      point_infoDetail(),
                      get_domain_information(),
                      input$concat_fasta,
                      input$label_archi_size,
                      input$title_archi_size)
      grid.draw(g)
      ggsave(file, plot = g,
             width = input$selected_width * 0.056458333,
             height = input$selected_height * 0.056458333,
             units = "cm", dpi = 300, device = "pdf", limitsize = FALSE)
      
    }
  )
  
  # ======================== FILTERED DATA DOWNLOADING ========================
  
  # * for main profile ==========================================================
  main_fasta_download <- reactive({
    main_fasta_out <- get_fasta_seqs(as.data.frame(download_data()),
                                     input$main_input, input$demo_data,
                                     input$input_type, input$concat_fasta,
                                     input$path,
                                     input$dir_format,
                                     input$file_ext,
                                     input$id_format,
                                     get_main_input())
    return(main_fasta_out)
  })
  download_data <- callModule(download_filtered_main,
                              "filtered_main_download", 
                              data = get_data_filtered,
                              fasta = main_fasta_download,
                              var1_id = reactive(input$var1_id),
                              var2_id = reactive(input$var2_id),
                              var1 = reactive(input$var1),
                              var2 = reactive(input$var2),
                              percent = reactive(input$percent))
  
  # * for customized profile ===================================================
  customized_fasta_download <- reactive({
    fasta_out_df <- get_fasta_seqs(as.data.frame(download_custom_data()),
                                   input$main_input, input$demo_data,
                                   input$input_type, input$concat_fasta,
                                   input$path,
                                   input$dir_format,
                                   input$file_ext,
                                   input$id_format,
                                   get_main_input())
  })
  download_custom_data <- callModule(download_filtered_customized,
                                     "filtered_customized_download",
                                     data = download_data,
                                     fasta = customized_fasta_download,
                                     in_seq = reactive(input$in_seq),
                                     in_taxa = reactive(input$in_taxa))
  
  # ============================ ANALYSIS FUNCTIONS ===========================
  
  # * PROFILE CLUSTERING ======================================================
  
  # ** check if anywhere elese genes are added to the custemized profile ---------
  observe({
    if (input$add_custom_profile == TRUE
        | input$add_core_gene_custom_profile == TRUE
        | input$add_gc_genes_custom_profile == TRUE){
      shinyjs::disable("add_cluster_cutom_profile")
    }else{
      shinyjs::enable("add_cluster_cutom_profile")
    }
  })
  
  output$add_cluster_cutom_profile_check.ui <- renderUI({
    if (input$add_custom_profile == TRUE
        | input$add_core_gene_custom_profile == TRUE |
        input$add_gc_genes_custom_profile == TRUE ){
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Gene age estimation</strong> or <strong>Core genes finding</strong> or <strong>Group Comparison</strong> &nbsp;to enable this function)</em></p>')
    }
  })
  
  # ** render cluster tree ----------------------------------------------------
  brushed_clusterGene <- callModule(cluster_profile, "profile_clustering",
                                    data = dataHeat,
                                    dist_method = reactive(input$dist_method),
                                    cluster_method = reactive(input$cluster_method),
                                    cluster_plot.width = reactive(input$cluster_plot.width),
                                    cluster_plot.height = reactive(input$cluster_plot.width),
                                    var1_aggregate_by = reactive(input$var1_aggregate_by))
  
  # * DISTRIBUTION ANALYSIS =====================================================
  
  # ** list of available variables for distribution plot -------------------------
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
  
  # ** var1 / var2 distribution data ---------------------------------------------
  distribution_df <- reactive({
    if (v$doPlot == FALSE) return()
    
    dataOrig <- get_main_input()
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
  
  # ** calculate % present species in supertaxa --------------------------------
  presSpecAllDt <- reactive({
    # open main input file
    mdData <- get_main_input()
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
  
  
  
  # ** render distribution plots ----------------------------------------------
  observe({
    if (v$doPlot == FALSE) return()
    
    if (is.null(input$selected_dist)){
      return()
    } else {
      if (input$selected_dist == "% present taxa"){
        callModule(analyze_distribution, "dist_plot",
                   data = presSpecAllDt, 
                   var_id = reactive(input$selected_dist),
                   var_type = reactive("presSpec"),
                   percent = reactive(input$percent), 
                   dist_text_size = reactive(input$dist_text_size))
      } else{
        if (input$selected_dist == input$var1_id){
          callModule(analyze_distribution, "dist_plot",
                     data = distribution_df,
                     var_id = reactive(input$selected_dist),
                     var_type = reactive("var1"),
                     percent = reactive(input$percent),
                     dist_text_size = reactive(input$dist_text_size))
        } else if (input$selected_dist == input$var2_id){
          callModule(analyze_distribution, "dist_plot",
                     data = distribution_df,
                     var_id = reactive(input$selected_dist),
                     var_type = reactive("var2"),
                     percent = reactive(input$percent),
                     dist_text_size = reactive(input$dist_text_size))
        }
      }
    }
  })
  
  # * GENE AGE ESTIMATION =======================================================
  
  # ** check if anywhere elese genes are added to the custemized profile ---------
  observe({
    if (input$add_cluster_cutom_profile == TRUE
        | input$add_core_gene_custom_profile == TRUE
        | input$add_gc_genes_custom_profile == TRUE ){
      shinyjs::disable("add_custom_profile")
    } else {
      shinyjs::enable("add_custom_profile")
    }
  })
  
  output$add_custom_profile_check.ui <- renderUI({
    if (input$add_cluster_cutom_profile == TRUE
        | input$add_core_gene_custom_profile == TRUE
        | input$add_gc_genes_custom_profile == TRUE){
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Profile clustering</strong> or <strong>Core genes finding</strong> or <strong>Groupcomparison</strong> &nbsp;to enable this function)</em></p>')
      
    }
  })
  
  # ** reset gene_age_prot_config ------------------------------------------------
  observeEvent(input$reset_gene_age_prot_config, {
    shinyjs::reset("gene_age_width")
    shinyjs::reset("gene_age_height")
    shinyjs::reset("gene_age_text")
  })
  
  # ** data for gene age estimation -------------------------------------------
  gene_ageDf <- reactive({
    if (v$doPlot == FALSE) return()
    
    gene_ageDf <- estimate_gene_age(subset_taxa(), get_data_filtered(),
                                    input$rank_select, input$in_select,
                                    input$var1, input$var2, input$percent)
    return(gene_ageDf)
  })
  
  # ** render age distribution plot -------------------------------------------
  selectedgene_age <- callModule(plot_gene_age, "gene_age",
                                 data = gene_ageDf,
                                 gene_age_width = reactive(input$gene_age_width),
                                 gene_age_height = reactive(input$gene_age_height),
                                 gene_age_text = reactive(input$gene_age_text))
  
  # * CORE GENES IDENTIFICATION ===============================================
  
  # ** render list of available taxa ---------------------------------------------
  output$taxa_list_core.ui <- renderUI({
    filein <- input$main_input
    # if(input$demo == TRUE){
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor"){
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
      
      if (input$apply_core_taxa == TRUE){
        out <- core_taxa_name()
        selectInput("taxa_core",
                    "Select taxa of interest:",
                    out,
                    selected = out,
                    multiple = TRUE)
      } else {
        selectInput("taxa_core",
                    "Select taxa of interest:",
                    out,
                    selected = out[1],
                    multiple = TRUE)
      }
    }
  })
  
  # ** render popup for selecting group of taxa to find core genes ------------
  core_taxa_name <- callModule(select_taxon_rank,
                               "select_taxon_rank_core",
                               rank_select = reactive(input$rank_select),
                               subset_taxa = subset_taxa)
  
  # check if anywhere elese genes are added to the custemized profile ---------
  observe({
    if (input$add_cluster_cutom_profile == TRUE
        | input$add_custom_profile == TRUE
        | input$add_gc_genes_custom_profile == TRUE){
      shinyjs::disable("add_core_gene_custom_profile")
    } else {
      shinyjs::enable("add_core_gene_custom_profile")
    }
  })
  
  output$add_core_gene_custom_profile_check.ui <- renderUI({
    if (input$add_cluster_cutom_profile == TRUE
        | input$add_custom_profile == TRUE
        | input$add_gc_genes_custom_profile == TRUE){
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Profiles clustering</strong> or <strong>Gene age estimating</strong> or r <strong>Group Comparioson</strong>&nbsp;to enable this function)</em></p>')
    }
  })
  
  # ** render table contains list of core genes -------------------------------
  core_geneDf <- callModule(identify_core_gene,
                            "core_gene",
                            filtered_data = get_data_filtered,
                            rank_select = reactive(input$rank_select),
                            taxa_core = reactive(input$taxa_core),
                            percent_core = reactive(input$percent_core))
  
  # ** download gene list from core_gene.table -----------------------------------
  output$core_gene_table_download <- downloadHandler(
    filename = function(){
      c("coreGeneList.out")
    },
    content = function(file){
      data_out <- core_geneDf()
      write.table(data_out, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  # * GET NCBI TAXONOMY IDs FROM INPUT LIST OF TAXON NAMES ====================
  callModule(search_taxon_id,"search_taxon_id")
  
  # * GROUP COMPARISON ========================================================
  
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
  
  # Dataframe with Information about the significant Genes --------------------
  # geneID | in_group| out_group | pvalues | features | databases | rank | var
  significant_genes_gc <- NULL
  
  # Select in_group -----------------------------------------------------------
  output$taxa_list_gc <- renderUI({
    filein <- input$main_input
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor"){
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
    
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor"){
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
                              get_parameter_input_gc(),
                              input$demo_data,
                              input$anno_location,
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
        input$add_core_gene_custom_profile == TRUE |
        input$add_cluster_cutom_profile == TRUE){
      shinyjs::disable("add_gc_genes_custom_profile")
    }else{
      shinyjs::enable("add_gc_genes_custom_profile")
    }
  })
  output$add_gc_custom_profile_check <- renderUI({
    if (input$add_custom_profile == TRUE |
        input$add_core_gene_custom_profile == TRUE |
        input$add_cluster_cutom_profile == TRUE){
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Gene age estimation</strong>or <strong>Profile clustering</strong> or <strong>Core genes finding</strong>&nbsp;to enable this function)</em></p>')
    }
  })
  
  # Deciding which plots should be shown (Group Comparison) -------------------
  get_plots_gc <- reactive({
    gene <- as.character(input$selected_gene_gc)
    input$plot_gc
    if (is.null(significant_genes_gc)) return()
    else if (gene == "all"){
      get_plot_output_list(significant_genes_gc,
                           get_parameter_input_gc(),
                           input$interesting_features)
    }else{
      x <- {
        significant_genes_gc[significant_genes_gc$geneID == gene, ]
      }
      if (nrow(x) == 0) return()
      get_plot_output_list(x, get_parameter_input_gc(), input$interesting_features)
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