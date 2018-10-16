#' Load packages
packages <- c("ggplot2", "reshape2",
              "plyr", "dplyr", "tidyr", "scales", "grid",
              "gridExtra", "ape", "stringr", "gtable",
              "dendextend", "ggdendro", "gplots", "data.table",
              "taxize", "zoo", "RCurl", "svMisc",
              "jmuOutlier")
source("R/functions.R")
install_packages(packages)
lapply(packages, library, character.only = TRUE)

#' Install bioconductor packages
bioconductor_pkgs <- c("Biostrings")
install_packages_bioconductor(bioconductor_pkgs)
lapply(bioconductor_pkgs, library, character.only = TRUE)

#' Import function files
source_files = list.files(path = "R",
                          pattern = "*.R$",
                          full.names = TRUE)
lapply(source_files, source, .GlobalEnv)

#' set size limit for input (9999mb)
options(shiny.maxRequestSize = 9999 * 1024 ^ 2)  # size limit for input 9999mb

#' MAIN SERVER ================================================================
shinyServer(function(input, output, session) {
  # Automatically stop a Shiny app when closing the browser tab
  # session$onSessionEnded(stopApp)
  session$allowReconnect(TRUE)
  
  # =========================== INITIAL CHECKING  =============================
  
  # * check for internet connection -------------------------------------------
  observe({
    if (has_internet() == FALSE) {
      toggleState("demo_data")
    }
  })
  
  output$no_internet_msg <- renderUI({
    if (has_internet() == FALSE) {
      strong(em("Internet connection is required for using demo data!"),
             style = "color:red")
    } else {
      return()
    }
  })
  
  # * check for the existence of taxonomy files -------------------------------
  observe({
    if (!file.exists(isolate("data/rankList.txt"))) {
      if (has_internet() == TRUE) {
        ncol <- max(count.fields("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/rankList.txt",
                                 sep = "\t"))
        df <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/rankList.txt",
                         sep = "\t",
                         quote = "",
                         header = FALSE,
                         fill = TRUE,
                         na.strings = c("", "NA"),
                         col.names = paste0("V", seq_len(ncol)))
        write.table(df, file = "data/rankList.txt",
                    col.names = FALSE,
                    row.names = FALSE,
                    quote = FALSE,
                    sep = "\t") #na = "",
      } else {
        file.create("data/rankList.txt")
      }
    }
  })
  
  observe({
    if (!file.exists(isolate("data/idList.txt"))) {
      if (has_internet() == TRUE) {
        ncol <- max(count.fields("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/idList.txt",
                                 comment.char = "",
                                 sep = "\t"))
        df <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/idList.txt",
                         sep = "\t",
                         header = FALSE,
                         fill = TRUE,
                         comment.char = "",
                         na.strings = c("", "NA"),
                         col.names = paste0("V", seq_len(ncol)))
        write.table(df, file = "data/idList.txt",
                    col.names = FALSE,
                    row.names = FALSE,
                    quote = FALSE,
                    sep = "\t") #na = "",
      } else {
        file.create("data/idList.txt")
      }
    }
  })
  
  observe({
    if (!file.exists(isolate("data/taxonNamesReduced.txt"))) {
      if (has_internet() == TRUE) {
        ncol <- max(count.fields("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/taxonNamesReduced.txt",
                                 sep = "\t"))
        df <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/taxonNamesReduced.txt",
                         sep = "\t",
                         quote = "",
                         header = FALSE,
                         fill = TRUE,
                         na.strings = c("", "NA"),
                         col.names = paste0("V", seq_len(ncol)))
        write.table(df, file = "data/taxonNamesReduced.txt",
                    col.names = FALSE,
                    row.names = FALSE,
                    quote = FALSE,
                    sep = "\t")
      } else {
        system("cp data/newTaxa.txt data/taxonNamesReduced.txt")
      }
    }
  })
  
  observe({
    if (!file.exists(isolate("data/taxonomyMatrix.txt"))) {
      if (has_internet() == TRUE) {
        ncol <- max(count.fields("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/taxonomyMatrix.txt",
                                 sep = "\t"))
        df <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/taxonomyMatrix.txt",
                         sep = "\t",
                         quote = "",
                         header = FALSE,
                         fill = TRUE,
                         na.strings = c("", "NA"),
                         col.names = paste0("V", seq_len(ncol)))
        
        write.table(df, file = "data/taxonomyMatrix.txt",
                    col.names = FALSE,
                    row.names = FALSE,
                    quote = FALSE,
                    sep = "\t")
      }
    }
  })
  
  # ======================== INPUT & SETTINGS TAB =============================
  # * check the validity of input file and render input_check.ui --------------
  output$input_check.ui <- renderUI({
    filein <- input$main_input
    if (is.null(filein)) return()
    input_type <- check_input_vadility(filein) #get_input_type()
    
    if (input_type[1] == "noGeneID") {
      updateButton(session, "do", disabled = TRUE)
      HTML(
        "<font color=\"red\"><em><strong>ERROR: Unsupported input format.
        <a href=\"https://github.com/BIONF/PhyloProfile/wiki/Input-Data\"
        target=\"_blank\">Click here for more info</a></em></strong></font>"
      )
    } else if (input_type[1] == "emptyCell") {
      updateButton(session, "do", disabled = TRUE)
      em(strong("ERROR: Rows have unequal length",
                style = "color:red"))
    }
    else if (input_type[1] == "moreCol") {
      updateButton(session, "do", disabled = TRUE)
      em(strong("ERROR: More columns than column names",
                style = "color:red"))
    }
    else {
      valid_type = c("xml", "fasta", "wide", "long", "oma")
      if (!(input_type[1] %in% valid_type)) {
        updateButton(session, "do", disabled = TRUE)
        invalid_oma <- paste(input_type, collapse = "; ")
        msg <- paste0("ERROR: Invalid IDs found! ", invalid_oma)
        em(strong(msg,
                  style = "color:red"))
      } else {
        updateButton(session, "do", disabled = FALSE)
        return()
      }
    }
  })
  
  # * render download link for Demo online files ------------------------------
  output$main_input_file.ui <- renderUI({
    if (input$demo_data == "lca-micros") {
      strong(a("Download demo input file",
               href = "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.long",
               target = "_blank"))
    } else if (input$demo_data == "ampk-tor") {
      strong(a("Download demo input file",
               href = "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile",
               target = "_blank"))
    } else {
      fileInput("main_input", h5("Upload input file:"))
    }
  })
  
  output$domain_input_file.ui <- renderUI({
    if (input$demo_data == "lca-micros") {
      strong(a("Download demo domain files",
               href = "https://github.com/BIONF/phyloprofile-data/tree/master/demo/domain_files",
               target = "_blank"))
    } else if (input$demo_data == "ampk-tor") {
      strong(a("Download demo domain file",
               href = "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.domains_F",
               target = "_blank"))
    } else {
      if (input$anno_location == "from file") {
        fileInput("file_domain_input", "")
      } else {
        textInput("domainPath", "", "")
      }
    }
  })
  
  output$download_fastaDemo.ui <- renderUI({
    if (input$demo_data == "lca-micros") {
      strong(a("Download demo fasta file",
               href = "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/fasta_file/concatenatedSeq.fa",
               target = "_blank"))
    } else if (input$demo_data == "ampk-tor") {
      strong(a("Download demo fasta file",
               href = "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.extended.fa",
               target = "_blank"))
    }
  })
  
  # * render description for Demo data ----------------------------------------
  output$demo_data_describe <- renderUI({
    if (input$demo_data == "none") {
      return()
    } else if (input$demo_data == "ampk-tor") {
      em(a("Data description",
           href = "https://github.com/BIONF/phyloprofile-data/blob/master/expTestData/ampk-tor/README.md",
           target = "_blank"))
    } else {
      em(a("Data description",
           href = "https://github.com/BIONF/phyloprofile-data/blob/master/demo/README.md",
           target = "_blank"))
    }
  })
  
  # * check OMA input ---------------------------------------------------------
  output$check_oma_input <- reactive({
    filein <- input$main_input
    if (is.null(filein)) return()
    input_type <- check_input_vadility(filein)
    input_type == "oma"
  })
  outputOptions(output, "check_oma_input", suspendWhenHidden = FALSE)
  
  # * download OMA data after parsing -----------------------------------------
  output$download_files_oma <- downloadHandler(
    filenname <- function() {
      "oma_data_to_phyloprofile_input.zip"
    },
    content <- function(file) {
      write.table(get_main_input(), "long.txt",
                  sep = "\t",
                  row.names = FALSE,
                  col.names = TRUE,
                  quote = FALSE)

      write.table(get_all_fasta_oma(final_oma_df()), "fasta.txt",
                  sep = "\t",
                  row.names = FALSE,
                  col.names = FALSE,
                  quote = FALSE)
      
      write.table(get_domain_information(), "domain.txt",
                  sep = "\t",
                  row.names = FALSE,
                  col.names = FALSE,
                  quote = FALSE)
      
      zip(zipfile = file,
          files = c("long.txt", "domain.txt", "fasta.txt"))
    },
    contentType = "application/zip"
  )
  
  # * close OMA parsing popup windows -------------------------------------
  observeEvent(input$get_data_oma, {
    toggleModal(session, "get_oma_data_windows", toggle = "close")
    updateButton(session, "get_data_oma", disabled = TRUE)
    toggleState("main_input")
  })
  
  # * render textinput for Variable 1 & 2 -------------------------------------
  output$var1_id.ui <- renderUI({
    long_dataframe <- get_main_input()
    if (is.null(long_dataframe)) {
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
    if (is.null(long_dataframe)) {
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
  
  # * render 2. variable relationship according to demo data ------------------
  output$var2_relation.ui <- renderUI({
    if (input$demo_data == "ampk-tor") {
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
  
  # * check the existance of the input concatenate fasta file -----------------
  output$concat_fasta.exist_check <- renderUI({
    if (is.null(input$concat_fasta)) return()
    else{
      f <- input$concat_fasta$datapath
      if (!file.exists(f)) {
        helpText("File not exists!!")
      } else {
        if (length(readLines(f, n = 1)) == 0) {
          helpText("is not a fasta file!!")
        } else {
          first_line <- readLines(f, n = 1)
          a <- substr(first_line, 1, 1)
          if (a == ">") {
            HTML('<p><span style="color: #0000ff;">
                 <strong>Please click CLOSE to comfirm!</strong></span></p>')
          } else {
            helpText("is not a fasta file!!")
          }
        }
        }
      }
    })
  
  # * check the validity of input tree file and render checkNewick.ui ---------
  check_newick_id <- reactive({
    req(input$inputTree)
    req(input$main_input)
    
    check_newick <- check_newick(input$inputTree, input$main_input, subset_taxa())
    if (check_newick == 0) {
      updateButton(session, "do", disabled = FALSE)
    }
    return(check_newick)
  })
  
  output$checkNewick.ui <- renderUI({
    check_newick <- check_newick_id()
    if (check_newick == 1) {
      updateButton(session, "do", disabled = TRUE)
      HTML("<p><em><span style=\"color: #ff0000;\"><strong>
           ERROR: Parenthesis(-es) missing!</strong></span></em></p>")
    } else if (check_newick == 2) {
      updateButton(session, "do", disabled = TRUE)
      HTML("<p><em><span style=\"color: #ff0000;\"><strong>
           ERROR: Comma(s) missing!</strong></span></em></p>")
    } else if (check_newick == 3) {
      updateButton(session, "do", disabled = TRUE)
      HTML("<p><em><span style=\"color: #ff0000;\"><strong>
           ERROR: Tree contains singleton!</strong></span></em></p>")
    } else if (check_newick == 0) {
      return()
    } else {
      updateButton(session, "do", disabled = TRUE)
      strong(em(paste0(check_newick, " not exist in main input file!")),
             style = "color:red")
    }
  })
  
  # * reset profile plot colors -----------------------------------------------
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
  
  # * render list of taxonomy ranks -------------------------------------------
  output$rank_select <- renderUI({
    if (input$demo_data == "lca-micros") {
      selectInput("rank_select", label = "",
                  choices = get_taxonomy_ranks(),
                  selected = "26_phylum")
    } else if (input$demo_data == "ampk-tor") {
      selectInput("rank_select", label = "",
                  choices = get_taxonomy_ranks(),
                  selected = "06_species")
    } else {
      selectInput("rank_select", label = "",
                  choices = get_taxonomy_ranks(),
                  selected = "06_species")
    }
  })
  
  # * render list of (super)taxa ----------------------------------------------
  output$select <- renderUI({
    choice <- alltaxa_list()
    choice$fullName <- as.factor(choice$fullName)
    
    if (input$demo_data == "lca-micros") {
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
    } else if (input$demo_data == "ampk-tor") {
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
  
  # * enable "PLOT" button ----------------------------------------------------
  observeEvent(input$rank_select,  ({
    if (input$rank_select == "") {
      updateButton(session, "do", disabled = TRUE)
    } else{
      unkTaxa <- unkTaxa()
      if (length(unkTaxa) == 0) {
        updateButton(session, "do", disabled = FALSE)
      }
    }
  }))
  # * move to main tab when "PLOT" button has been clicked --------------------
  observe({
    # use tabsetPanel "id" argument to change tabs
    if (input$do > 0) {
      updateTabsetPanel(session, "tabs", selected = "Main profile")
    }
  })
  
  
  # * disable main input, genelist input and demo data checkbox ---------------
  observe({
    if (input$do > 0) {
      toggleState("main_input")
      toggleState("gene_list_selected")
      toggleState("demo_data")
    }
  })
  
  # * update var2_aggregate_by to mean if using demo lca-micros data ----------
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
  
  # * render filter slidebars for Main plot -----------------------------------
  output$var1_cutoff.ui <- renderUI({
    create_slider_cutoff(
      "var1", paste(input$var1_id, "cutoff:"), 0.0, 1.0, input$var1_id
    )

    # numericInput("var1",
    #              paste(input$var1_id, "cutoff:"),
    #              min = 0,
    #              max = 1,
    #              value = 0)
  })
  
  output$var2_cutoff.ui <- renderUI({
    create_slider_cutoff(
      "var2", paste(input$var2_id, "cutoff:"), 0.0, 1.0, input$var2_id
    )
  })
  
  output$percent_cutoff.ui <- renderUI({
    create_slider_cutoff(
      "percent", "% of present taxa:", 0.0, 1.0, "percent"
    )
  })
  
  # * render filter slidebars for Customized plot -----------------------------
  output$var1_filter.ui <- renderUI({
    create_slider_cutoff(
      "var1cus",
      paste(input$var1_id, "cutoff:"),
      input$var1[1], input$var1[2], input$var1_id
    )
  })
  
  output$var2_filter.ui <- renderUI({
    create_slider_cutoff(
      "var2cus",
      paste(input$var2_id, "cutoff:"),
      input$var2[1], input$var2[2], input$var2_id
    )
  })
  
  output$percent_filter.ui <- renderUI({
    create_slider_cutoff(
      "percent2",
      "% of present taxa:",
      input$percent[1], input$percent[2], "percent"
    )
  })
  
  # * render filter slidebars for Distribution plot ---------------------------
  output$var1_dist.ui <- renderUI({
    create_slider_cutoff(
      "var1_dist",
      paste(input$var1_id, "cutoff:"),
      input$var1[1], input$var1[2], input$var1_id
    )
  })
  
  output$var2_dist.ui <- renderUI({
    create_slider_cutoff(
      "var2_dist",
      paste(input$var2_id, "cutoff:"),
      input$var2[1], input$var2[2], input$var2_id
    )
  })
  
  output$percent_dist.ui <- renderUI({
    create_slider_cutoff(
      "percent_dist",
      "% of present taxa:",
      input$percent[1], input$percent[2], "percent"
    )
  })
  
  # * render filter slidebars for Gene age estimation plot --------------------
  output$var1_age.ui <- renderUI({
    create_slider_cutoff(
      "var1_age",
      paste(input$var1_id, "cutoff:"),
      input$var1[1], input$var1[2], input$var1_id
    )
  })
  
  output$var2_age.ui <- renderUI({
    create_slider_cutoff(
      "var2_age",
      paste(input$var2_id, "cutoff:"),
      input$var2[1], input$var2[2], input$var2_id
    )
  })
  
  output$percent_age.ui <- renderUI({
    create_slider_cutoff(
      "percent_age",
      "% of present taxa:",
      input$percent[1], input$percent[2], "percent"
    )
  })
  
  # * render filter slidebars for Core gene finding function ------------------
  output$var1_core.ui <- renderUI({
    create_slider_cutoff(
      "var1_core",
      paste(input$var1_id, "cutoff:"),
      input$var1[1], input$var1[2], input$var1_id
    )
  })
  
  output$var2_core.ui <- renderUI({
    create_slider_cutoff(
      "var2_core",
      paste(input$var2_id, "cutoff:"),
      input$var2[1], input$var2[2], input$var2_id
    )
  })
  
  output$percent_core.ui <- renderUI({
    create_slider_cutoff(
      "percent_core",
      "% of present taxa:",
      input$percent[1], input$percent[2], "percent"
    )
  })
  
  # * update value for filter slidebars of Main Plot --------------------------
  # ** based on customized profile
  observe({
    new_var1 <- input$var1cus
    update_slider_cutoff(
      session,
      "var1", paste(input$var1_id, "cutoff:"), new_var1, input$var1_id
    )
  })
  
  observe({
    new_var2 <- input$var2cus
    update_slider_cutoff(
      session,
      "var2", paste(input$var2_id, "cutoff:"), new_var2, input$var2_id
    )
  })
  
  observe({
    new_percent <- input$percent2
    update_slider_cutoff(
      session,
      "percent", "% of present taxa:", new_percent, "percent"
    )
  })
  
  # ** based on "Distribution analysis"
  observe({
    new_var1 <- input$var1_dist
    update_slider_cutoff(
      session,
      "var1", paste(input$var1_id, "cutoff:"), new_var1, input$var1_id
    )
  })
  
  observe({
    new_var2 <- input$var2_dist
    update_slider_cutoff(
      session,
      "var2", paste(input$var2_id, "cutoff:"), new_var2, input$var2_id
    )
  })
  
  observe({
    new_percent <- input$percent_dist
    update_slider_cutoff(
      session,
      "percent", "% of present taxa:", new_percent, "percent"
    )
  })
  
  # ** based on "Gene age estimation"
  observe({
    new_var1 <- input$var1_age
    update_slider_cutoff(
      session,
      "var1", paste(input$var1_id, "cutoff:"), new_var1, input$var1_id
    )
  })
  
  observe({
    new_var2 <- input$var2_age
    update_slider_cutoff(
      session,
      "var2", paste(input$var2_id, "cutoff:"), new_var2, input$var2_id
    )
  })
  
  observe({
    new_percent <- input$percent_age
    update_slider_cutoff(
      session,
      "percent", "% of present taxa:", new_percent, "percent"
    )
  })
  
  # * reset cutoffs of Main plot ----------------------------------------------
  observeEvent(input$reset_main, {
    shinyjs::reset("var1")
    shinyjs::reset("var2")
    shinyjs::reset("percent")
  })
  
  # * reset cutoffs of Customized plot ----------------------------------------
  observeEvent(input$reset_selected, {
    shinyjs::reset("var1")
    shinyjs::reset("var2")
    shinyjs::reset("percent")
  })
  
  # ========================= PARSING UNKNOWN TAXA ============================
  
  # * get list of "unknown" taxa in main input --------------------------------
  unkTaxa <- reactive({
    long_dataframe <- get_main_input()
    req(long_dataframe)
    
    if (is.null(long_dataframe)) {
      inputTaxa <- c("NA")
    } else{
      inputTaxa <- levels(long_dataframe$ncbiID)
    }
    
    if (inputTaxa[1] == "NA") {
      return()
    } else {
      inputTaxa <- unlist(strsplit(inputTaxa, split = "\t"))
      if (inputTaxa[1] == "geneID") {
        # remove "geneID" element from vector inputTaxa
        inputTaxa <- inputTaxa[-1]
      }
      
      if (!file.exists(isolate("data/rankList.txt"))) {
        return(inputTaxa)
      } else {
        info <- file.info("data/rankList.txt")
        if (info$size == 0) {
          return(inputTaxa)
        } else {
          # get list of all available taxon (from /data/rankList.txt)
          pipe_cmd <- paste0("cut -f 1 ", getwd(), "/data/rankList.txt")
          allTaxa <- unlist((read.table(pipe(pipe_cmd))))
          
          # list of unknown taxa
          unkTaxa <- inputTaxa[!(inputTaxa %in% allTaxa)]
          if (identical(unkTaxa, character(0))) return()
          
          # get non-ncbi taxa
          unkTaxa <- data.frame("TaxonID" = unkTaxa)
          unkTaxa$id <- substring(unkTaxa$TaxonID, 5)
          unkTaxa$Source <- "ncbi"
          
          pipe_ncbi <- paste0("cut -f 1 ", getwd(), "/data/taxonNamesFull.txt")
          ncbiTaxa <- unlist((read.table(pipe(pipe_ncbi))))
          
          ncbiID <- levels(ncbiTaxa)
          maxNCBI <- max(sort(as.numeric(ncbiID[ncbiID != "ncbiID"])))
          
          if (nrow(unkTaxa[!(unkTaxa$id %in% ncbiTaxa),]) > 0) {
            unk_taxa <- unkTaxa[!(unkTaxa$id %in% ncbiTaxa),]$id
            unkTaxa[unkTaxa$id %in% unk_taxa,]$Source <- "unknown"
            if (any(unk_taxa < maxNCBI)) {
              unkTaxa[unkTaxa$id %in% unk_taxa &
                        unkTaxa$id < maxNCBI,]$Source <- "invalid"
            }
          }
          
          if (nrow(unkTaxa[unkTaxa$id %in% newTaxa,]) > 0) {
            unkTaxa[unkTaxa$id %in% newTaxa,]$Source <- "new"
          }
          
          # check for invalid newly generated IDs in newTaxa.txt file
          pipe_new <- paste0("cut -f 1 ", getwd(), "/data/newTaxa.txt")
          newTaxa <- unlist((read.table(pipe(pipe_new))))

          if (length(newTaxa) > 1) {
            newTaxaList <- levels(newTaxa)
            newTaxaList <- as.integer(newTaxaList[newTaxaList != "ncbiID"])
            
            if (min(newTaxaList) < maxNCBI) {
              invalidList <- as.data.frame(newTaxaList[newTaxaList < maxNCBI])
              colnames(invalidList) <- c("id")
              invalidList$TaxonID <- "newTaxa.txt"
              invalidList$Source <- "invalid"
              unkTaxa <- rbind(invalidList, unkTaxa)
            }
          }
          
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
    if (length(unkTaxa) > 0) {
      if ("invalid" %in% unkTaxa$Source) return("invalid")
      if ("unknown" %in% unkTaxa$Source) return("unknown")
      else return("ncbi")
    } else {
      return(0)
    }
  })
  outputOptions(output, "unk_taxa_status", suspendWhenHidden = FALSE)
  
  # * render list of unkTaxa --------------------------------------------------
  output$unk_taxa_full <-
    renderDataTable(options = list(searching = FALSE, pageLength = 10),{
      if (length(unkTaxa()) > 0) {
        tb <- unkTaxa()
        tb[, c("TaxonID", "Source")]
      }
    })
  
  # * download list of unkTaxa ------------------------------------------------
  output$unk_taxa.download <- downloadHandler(
    filename = function() {
      c("unknown_taxa.txt")
    },
    content = function(file) {
      data_out <- unkTaxa()
      data_out <- data_out[, c("TaxonID", "Source")]
      write.table(data_out, file,
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)
    }
  )
  
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
  
  # * close adding taxa windows -----------------------------------------------
  observeEvent(input$new_done, {
    toggleModal(session, "add_taxa_windows", toggle = "close")
    write.table(newTaxa$Df, "data/newTaxa.txt",
                sep = "\t",
                eol = "\n",
                row.names = FALSE,
                quote = FALSE)
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
  
  # * create rankList.txt, idList.txt, ----------------------------------------
  invalidID <- reactive({
    invalidID <- data.frame("id" = as.character(),
                            "type" = as.character(),
                            stringsAsFactors = FALSE)
    
    filein <- input$main_input
    if (is.null(filein)) return()
    input_type <- check_input_vadility(filein) #get_input_type()
    
    if (input_type == "xml" |
        input_type == "long" |
        input_type == "wide" |
        input_type == "fasta" ) {
      inputDf <- as.data.frame(read.table(file = filein$datapath,
                                          sep = "\t",
                                          header = TRUE,
                                          check.names = FALSE,
                                          comment.char = ""))
      
      if (v1$parse == FALSE) return()
      else {
        # get list of taxa need to be parsed (the one mising taxonomy info)
        if (v1$parse == TRUE) {
          unkTaxaDf <- unkTaxa()
          unkTaxa <- as.character(substring(unkTaxaDf$TaxonID, 5))
          titleline <- c("geneID", unkTaxa)
        }
        # invalidIDtmp <- list()
        
        ncbiTaxonInfo <- fread("data/taxonNamesFull.txt")
        newTaxaFromFile <- fread("data/newTaxa.txt",
                                 colClasses = c("ncbiID" = "character"))
        
        ## join all ncbi taxa and new taxa together
        allTaxonInfo <- rbind(newTaxaFromFile, ncbiTaxonInfo)
        
        ## check missing ids
        if (any(!(unkTaxa %in% allTaxonInfo$ncbiID))) {
          invalid_missing <-
            unkTaxa[!(unkTaxa %in% allTaxonInfo$ncbiID)]
          invalidID_tmp <- data.frame(
            "id" = invalid_missing,
            "type" = rep("missing", length(invalid_missing))
          )
          invalidID <- rbind(invalidID, invalidID_tmp)
        }
        
        ## check IDs & names from newTaxa that are present in taxonNamesFull
        if (nrow(newTaxaFromFile[newTaxaFromFile$ncbiID
                                 %in% ncbiTaxonInfo$ncbiID,]) > 0) {
          invalid_id <- newTaxaFromFile[newTaxaFromFile$ncbiID
                                        %in% ncbiTaxonInfo$ncbiID,]$ncbiID
          invalidID_tmp <- data.frame(
            "id" = invalid_id,
            "type" = rep("id", length(invalid_id))
          )
          invalidID <- rbind(invalidID, invalidID_tmp)
          
          newTaxaFromFile <- newTaxaFromFile[!(newTaxaFromFile$ncbiID
                                               %in% ncbiTaxonInfo$ncbiID),]
        }
        
        if (nrow(newTaxaFromFile[newTaxaFromFile$fullName
                                 %in% ncbiTaxonInfo$fullName,]) > 0) {
          invalid_name <- newTaxaFromFile[newTaxaFromFile$fullName
                                          %in% ncbiTaxonInfo$fullName,]$ncbiID
          invalidID_tmp <- data.frame(
            "id" = invalid_name,
            "type" = rep("name", length(invalid_name))
          )
          invalidID <- rbind(invalidID, invalidID_tmp)
        }
        
        if (nrow(invalidID) > 0) {
          return(invalidID)
        }
        
        ## get unique norank IDs and unique rank
        ## (which are uninformative for sorting taxa)
        allNorankIDs <- list()
        allRanks <- list()
        
        for (i in 2:length(titleline)) {
          # taxon ID
          refID <- titleline[i]
          # get info for this taxon
          refEntry <- allTaxonInfo[allTaxonInfo$ncbiID == refID, ]
          # parent ID
          lastID <- refEntry$parentID
          
          # 
          if (refEntry$rank == "norank") {
            allNorankIDs <- c(allNorankIDs, refEntry$ncbiID)
          }
          allRanks <- c(allRanks, refEntry$rank)
          
          while (lastID != 1) {
            nextEntry <- allTaxonInfo[allTaxonInfo$ncbiID == lastID, ]
            if (nextEntry$rank == "norank") {
              allNorankIDs <- c(allNorankIDs, nextEntry$ncbiID)
            }
            allRanks <- c(allRanks, nextEntry$rank)
            lastID <- nextEntry$parentID
          }
          
          # print progress
          # p <- (i - 1) / (length(titleline) - 1) * 100
          # progress(p)
        }
        uniqueRank <- names(which(table(as.character(allRanks)) == 1))
        uniqueID <- names(which(table(as.character(allNorankIDs)) == 1))
        overID <- names(which(table(as.character(allNorankIDs)) == 
                                (length(titleline) - 1)))
        uniqueID <- c(uniqueID, overID)

        ## parse taxonomy info
        rankList <- data.frame()
        idList <- data.frame()
        reducedInfoList <- data.frame()
        
        # Create 0-row data frame which will be used to store data
        withProgress(message = "Parsing input file", value = 0, {
          for (i in 2:length(titleline)) {
            ## taxon ID
            refID <- titleline[i]
            
            ## get info for this taxon
            refEntry <- allTaxonInfo[allTaxonInfo$ncbiID == refID, ]
            
            if (
              nrow(reducedInfoList[
                reducedInfoList$X1 == refEntry$ncbiID, ]) == 0
            ) {
              refInfoList <- data.frame(matrix(c(refEntry$ncbiID,
                                                 refEntry$fullName,
                                                 refEntry$rank,
                                                 refEntry$parentID),
                                               nrow = 1,
                                               byrow = TRUE),
                                        stringsAsFactors = FALSE)
              reducedInfoList <- rbind(reducedInfoList, refInfoList)
            }
            
            # parentID (used to check if hitting last rank, i.e. norank_1)
            lastID <- refEntry$parentID
            
            # create list of rank for this taxon
            rank <- c(paste0("ncbi", refID), refEntry$fullName)
            if (refEntry$rank == "norank") {
              rank <- c(rank, paste0("strain"))
            } else {
              rank <- c(rank, refEntry$rank)
            }
            
            # create list of IDs for this taxon
            ids <- list(paste0(refEntry$fullName, "#name"))
            if (refEntry$rank == "norank") {
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
            while (lastID != 1) {
              nextEntry <- allTaxonInfo[allTaxonInfo$ncbiID == lastID, ]
              
              if (
                nrow(reducedInfoList[
                  reducedInfoList$X1 == nextEntry$ncbiID, ]) == 0
              ) {
                nextEntryList <-
                  data.frame(matrix(c(nextEntry$ncbiID,
                                      nextEntry$fullName,
                                      nextEntry$rank,
                                      nextEntry$parentID),
                                    nrow = 1, byrow = TRUE),
                             stringsAsFactors = FALSE)
                
                reducedInfoList <- rbind(reducedInfoList,
                                         nextEntryList)
              }
              
              lastID <- nextEntry$parentID
              
              if ("norank" %in% nextEntry$rank) {
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
            
            # change "no_rank" before species into "strain"
            if (rank[3] == "norank" & any(rank == "species")) {
              rank[3] = "strain"
              tmpID <- unlist(strsplit(as.character(ids[2]), split = "#"))
              ids[2] <- paste0(tmpID[1],"#","strain")
            }

            # append into rankList and idList files
            rankListTMP <- data.frame(matrix(unlist(rank),
                                             nrow = 1, byrow = TRUE),
                                      stringsAsFactors = FALSE)
            rankList <- rbind.fill(rankList, rankListTMP)
            idListTMP <- data.frame(matrix(unlist(ids),
                                           nrow = 1,
                                           byrow = TRUE),
                                    stringsAsFactors = FALSE)
            idList <- rbind.fill(idList, idListTMP)
            
            # Increment the progress bar, and update the detail text.
            incProgress(1 / (length(titleline) - 1),
                        detail = paste( (i - 1),
                                        "/",
                                        length(titleline) - 1))
          }
        })
        
        withProgress(message = "Generating taxonomy file...", value = 0, {
          # open existing files (idList, rankList and taxonNamesReduced.txt)
          ncol <- max(count.fields("data/rankList.txt", sep = "\t"))
          oldIDList <-
            as.data.frame(read.table("data/idList.txt",
                                     sep = "\t",
                                     header = FALSE,
                                     check.names = FALSE,
                                     comment.char = "",
                                     fill = TRUE,
                                     stringsAsFactors = TRUE,
                                     na.strings = c("", "NA"),
                                     col.names = paste0("X", seq_len(ncol))))
          
          oldRankList <-
            as.data.frame(read.table("data/rankList.txt",
                                     sep = "\t",
                                     header = FALSE,
                                     check.names = FALSE,
                                     comment.char = "",
                                     fill = TRUE,
                                     stringsAsFactors = TRUE,
                                     na.strings = c("", "NA"),
                                     col.names = paste0("X", seq_len(ncol))))
          
          oldNameList <-
            as.data.frame(read.table("data/taxonNamesReduced.txt",
                                     sep = "\t",
                                     header = TRUE,
                                     check.names = FALSE,
                                     comment.char = "",
                                     fill = TRUE,
                                     stringsAsFactors = TRUE))
          
          
          # and append new info into those files
          new_idList <- rbind.fill(oldIDList, idList)
          new_rankList <- rbind.fill(oldRankList, rankList)
          colnames(reducedInfoList) <- c("ncbiID",
                                         "fullName",
                                         "rank",
                                         "parentID")
          new_nameList <- rbind.fill(oldNameList,
                                     reducedInfoList)

          # write output files (idList, rankList and taxonNamesReduced)
          write.table(new_idList[!duplicated(new_idList), ],
                      file  = "data/idList.txt",
                      col.names = FALSE,
                      row.names = FALSE,
                      quote = FALSE,
                      sep = "\t")
          write.table(new_rankList[!duplicated(new_rankList), ],
                      file = "data/rankList.txt",
                      col.names = FALSE,
                      row.names = FALSE,
                      quote = FALSE,
                      sep = "\t")
          write.table(new_nameList[!duplicated(new_nameList), ],
                      file = "data/taxonNamesReduced.txt",
                      col.names = TRUE,
                      row.names = FALSE,
                      quote = FALSE,
                      sep = "\t")
          

          # create taxonomy matrix (taxonomyMatrix.txt)
          taxMatrix <- taxonomyTableCreator("data/idList.txt",
                                            "data/rankList.txt")
          write.table(taxMatrix,
                      file = "data/taxonomyMatrix.txt",
                      sep = "\t",
                      eol = "\n",
                      row.names = FALSE,
                      quote = FALSE)
        })
      }
    }
    return()
  })
  
  # * output invalid NCBI ID --------------------------------------------------
  output$invalidID.output <- renderTable({
    if (is.null(invalidID())) return()
    else{
      outDf <- invalidID()
      colnames(outDf) <- c("Invalid ID(s)", "Type")
      return(outDf)
    }
  })
  
  # * download list of invalidID ----------------------------------------------
  output$invalidID.download <- downloadHandler(
    filename = function() {
      c("invalid_ids.txt")
    },
    content = function(file) {
      data_out <- invalidID()
      colnames(data_out) <- c("Invalid ID(s)", "Type")
      write.table(data_out, file,
                  sep = "\t",
                  row.names = FALSE,
                  quote = FALSE)
    }
  )
  
  # * render final msg after taxon parsing ------------------------------------
  output$end_parsing_msg <- renderUI({
    if (is.null(invalidID())) {
      strong(h4("PLEASE RELOAD THIS TOOL WHEN FINISHED!!!"),
             style = "color:red")
    } else {
      HTML('<p><strong><span style="color: #e12525;"> SOME INVALID TAXON
           IDs HAVE BEEN FOUND!!</span><br /> </strong></p>
           <p><em>Type="<span style="color: #0000ff;">id</span>"/
           <span style="color: #0000ff;">name</span>:
           IDs/names already exist in NCBI!</em></p>
           <p><em>Type="<span style="color: #0000ff;">missing</span>": IDs
           cannot be found in both NCBI and newTaxa.txt file.</em></p>
           <p>For IDs with type of <em><span style="color: #0000ff;">"id"
           </span></em> and <em><span style="color: #0000ff;">"name"</span>
           </em>, please remove them from newTaxa.txt file or
           renamed their IDs and names.</p>
           <p>For IDs with type of <em><span style="color: #0000ff;">"missing"
           </span></em>, please check the validity of them&nbsp;in
           <a href="https://www.ncbi.nlm.nih.gov/taxonomy" target="_blank"
           rel="noopener"> NCBI taxonomy database</a>!</p>')
    }
    })
  
  # ====================== PROCESSING INPUT DATA ==============================
  
  # * check if data is loaded and "plot" button is clicked --------------------
  v <- reactiveValues(doPlot = FALSE)
  observeEvent(input$do, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    v$doPlot <- input$do
    filein <- input$main_input
    if (is.null(filein) & input$demo_data == "none") {
      v$doPlot <- FALSE
      updateButton(session, "do", disabled = TRUE)
    }
  })
  
  # * check if "no ordering gene IDs" has been checked ------------------------
  output$apply_cluster_check.ui <- renderUI({
    if (input$ordering == FALSE) {
      HTML('<p><em>(Check "Ordering sequence IDs" check box in
           <strong>Input & settings tab</strong>&nbsp;to enable this function)
           </em></p>')
    }
    })
  
  # * to enable clustering ----------------------------------------------------
  observe({
    if (input$ordering == FALSE) {
      shinyjs::disable("apply_cluster")
    } else {
      shinyjs::enable("apply_cluster")
    }
  })
  
  # * get OMA data for input list ---------------------------------------------
  get_oma_browser <- function(id_list, ortho_type) {
    final_omaDf <- data.frame()
    
    withProgress(value = 0, {
      i = 1
      for (seed_id in id_list) {
        incProgress(1 / length(id_list),
                    message = paste0("Retrieving OMA for ", seed_id),
                    detail = paste(i, "/", length(id_list)))
        # get members
        members <- get_members(seed_id, ortho_type)
        oma_seed_id <- OmaDB::getData("protein",seed_id)$omaid
        # get all data
        withProgress(message = "Ortholog ", value = 0, {
          j <- 1
          for (ortho in members) {
            orthoDf <- get_data_for_one_oma(ortho)
            orthoDf$seed <- seed_id
            if (ortho == oma_seed_id) {
              orthoDf$ortho_id <- seed_id
            }
            final_omaDf <- rbind(final_omaDf, orthoDf)
            incProgress(1 / length(members),
                        detail = paste(j, "/", length(members)))
            j <- j + 1
          }
        })
        i <- i + 1
      }
    })
    
    return(final_omaDf)
  }
  
  final_oma_df <- reactive({
    filein <- input$main_input
    if (is.null(filein)) return()
    input_type <- check_input_vadility(filein)
    
    if (input_type == "oma") {
      if (input$get_data_oma[1] == 0) return()
      oma_ids <- as.data.frame(read.table(file = filein$datapath,
                                          sep = "\t",
                                          header = FALSE,
                                          check.names = FALSE,
                                          comment.char = ""))
      oma_ids[,1] <- as.character(oma_ids[,1])
      final_oma_df <- get_oma_browser(oma_ids[,1], input$selected_oma_type)
      return(final_oma_df)
    } else {
      return()
    }
  })
  
  # * convert main input file in any format into long format dataframe --------
  get_main_input <- reactive({
    if (input$demo_data == "lca-micros") {
      long_dataframe <- create_long_matrix("lca-micros")
    } else if (input$demo_data == "ampk-tor") {
      long_dataframe <- create_long_matrix("ampk-tor")
    } else {
      filein <- input$main_input
      if (is.null(filein)) return()
      input_type <- check_input_vadility(filein)
      if (input_type == "oma") {
        if (input$get_data_oma[1] == 0) return()
        long_dataframe <- create_profile_from_oma(final_oma_df())
        for (i in 1:ncol(long_dataframe)) {
          long_dataframe[, i] <- as.factor(long_dataframe[, i])
        }
      } else {
        long_dataframe <- create_long_matrix(filein)
      }
    }
    return(long_dataframe)
  })
  
  # * parse domain info into data frame ---------------------------------------
  get_domain_information <- reactive({
    if (input$demo_data == "none") {
      filein <- input$main_input
      input_type <- check_input_vadility(filein)
    } else {
      input_type <- "demo"
    }
    
    if (input_type == "oma") {
      domain_df <- get_all_domains_oma(final_oma_df())
    } else {
      print("Getting the domains...")
      main_input <- get_main_input()
      domain_df <- parse_domain_input(main_input,
                                      # input_type,
                                      input$demo_data,
                                      input$anno_location,
                                      input$file_domain_input,
                                      input$domainPath,
                                      session,
                                      datapath)
    }
    return(domain_df)
  })
  
  # * get ID list of input taxa from main input -------------------------------
  subset_taxa <- reactive({
    if (input$demo_data == "lca-micros" |
        input$demo_data == "ampk-tor" |
        length(unkTaxa()) == 0) {
      long_dataframe <- get_main_input()
      if (is.null(long_dataframe)) return()
      inputTaxa <- levels(long_dataframe$ncbiID)
    } else {
      inputTaxa <- readLines(filein$datapath, n = 1)
    }
    
    inputTaxa <- unlist(strsplit(inputTaxa, split = "\t"))
    if (inputTaxa[1] == "geneID") {
      # remove "geneID" element from vector inputTaxa
      inputTaxa <- inputTaxa[-1]
    }
    # return input taxa
    return(inputTaxa)
  })
  
  # * get NAME list of all (super)taxa ----------------------------------------
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
                                         header = TRUE,
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
                                          header = TRUE))
    allTaxa <- alltaxa_list()
    rankNameTMP <- allTaxa$rank[allTaxa$fullName == input$in_select]
    if (rankName == "strain") {
      superID <-
        as.numeric(taxa_list$ncbiID[taxa_list$fullName == input$in_select
                                    & taxa_list$rank == "norank"])
    } else {
      superID <-
        as.numeric(taxa_list$ncbiID[taxa_list$fullName == input$in_select
                                    & taxa_list$rank == rankNameTMP[1]])
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
    
    if (is.null(treeIn)) {
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
    
    if (nrow(sortedOut) > 1) {
      for (i in 2:nrow(sortedOut)) {
        ## increase prefix if changing to another supertaxon
        if (sortedOut$fullName.y[i] != sortedOut$fullName.y[i - 1]) {
          prefix <- prefix + 1
        }
        sortedOut$sortedSupertaxon[i] <- paste0(prefix,
                                                "_",
                                                sortedOut$fullName.y[i])
        sortedOut$species[i] <- paste0("ncbi",
                                       sortedOut$species[i])
      }
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
    if (input$gene_list_selected == "from file") {
      listIn <- input$list
      if (!is.null(listIn)) {
        list <- as.data.frame(read.table(file = listIn$datapath,
                                         header = FALSE))
        listGeneOri <- list$V1
        if (input$st_index <= length(listGeneOri)) {
          listGene <- listGeneOri[listGeneOri[input$st_index:end_index]]
        } else {
          listGene <- listGeneOri
        }
      }
    }
    
    long_dataframe <- get_main_input()
    if (is.null(long_dataframe)) return()
    
    if (is.null(long_dataframe)) {
      data <- data.frame("geneID" = character(),
                         "ncbiID" = character(),
                         "orthoID" = character(),
                         "var1" = character(),
                         "var2" = character(),
                         stringsAsFactors = FALSE)
    } else {
      long_dataframe <- unsort_id(long_dataframe, input$ordering)
      
      if (length(listGene) >= 1) {
        data <- long_dataframe[long_dataframe$geneID %in% listGene, ]
      } else {
        subsetID <- levels(long_dataframe$geneID)[input$st_index:end_index]
        data <- long_dataframe[long_dataframe$geneID %in% subsetID, ]
      }
     
      if (ncol(data) < 5) {
        for (i in 1:(5 - ncol(data))) {
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
    if (is.null(preData())) return()

    mdData <- preData()
    
    # count number of inparalogs
    paralogCount <- plyr::count(mdData, c("geneID", "ncbiID"))
    mdData <- merge(mdData, paralogCount, by = c("geneID", "ncbiID"))
    colnames(mdData)[ncol(mdData)] <- "paralog"
    
    # (1) GET SORTED TAXONOMY LIST
    taxa_list <- sortedtaxa_list()
    
    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxa_list, "supertaxon")
    
    # merge mdData, mdDataVar2 and taxa_list to get taxonomy info
    taxaMdData <- merge(mdData, taxa_list, by = "ncbiID")
    taxaMdData$var1 <-
      suppressWarnings(as.numeric(as.character(taxaMdData$var1)))
    taxaMdData$var2 <-
      suppressWarnings(as.numeric(as.character(taxaMdData$var2)))
    
    # (2) calculate PERCENTAGE of PRESENT SPECIES
    finalPresSpecDt <- calc_pres_spec(taxaMdData, taxaCount)
    # (3) calculate max/min/mean/median VAR1 for every supertaxon of each gene
    # remove NA rows from taxaMdData
    taxaMdDataNoNA <- taxaMdData[!is.na(taxaMdData$var1), ]
    # calculate m var1
    mVar1Dt <- aggregate(taxaMdDataNoNA[, "var1"],
                         list(taxaMdDataNoNA$supertaxon,
                              taxaMdDataNoNA$geneID),
                         FUN = input$var1_aggregate_by)
    colnames(mVar1Dt) <- c("supertaxon", "geneID", "mVar1")
    
    # (4) calculate max/min/mean/median VAR2 for each super taxon
    # remove NA rows from taxaMdData
    taxaMdDataNoNA_var2 <- taxaMdData[!is.na(taxaMdData$var2), ]
    # calculate max/min/mean/median VAR2
    if (nrow(taxaMdDataNoNA_var2) > 0) {
      mVar2Dt <- aggregate(taxaMdDataNoNA_var2[, "var2"],
                           list(taxaMdDataNoNA_var2$supertaxon,
                                taxaMdDataNoNA_var2$geneID),
                           FUN = input$var2_aggregate_by)
      colnames(mVar2Dt) <- c("supertaxon", "geneID", "mVar2")
    } else {
      mVar2Dt <- taxaMdData[, c("supertaxon", "geneID")]
      mVar2Dt$mVar2 <- 0
    }
    
    # (3+4) & join mVar2 together with mVar1 scores into one df
    scoreDf <- merge(mVar1Dt,
                     mVar2Dt,
                     by = c("supertaxon", "geneID"),
                     all = TRUE)
    
    # (2+3+4) add presSpec and mVar1 into taxaMdData
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
  
  # * reduce data from lowest level to supertaxon (e.g. phylum) ---------------
  # * This data set contain only supertaxa
  # * and their value (%present, mVar1 & mVar2) for each gene
  dataSupertaxa <- reactive({
    fullMdData <- get_data_filtered()
    
    # to check if working with the lowest taxonomy rank; 1 for NO; 0 for YES
    flag <- 1
    if (length(unique(levels(as.factor(fullMdData$numberSpec)))) == 1) {
      if (unique(levels(as.factor(fullMdData$numberSpec))) == 1) {
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
    
    if (flag == 1) {
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
    return(superDfExt)
  })
  
  # * heatmap data input ------------------------------------------------------
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
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor") {
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
    
    if (input$var1_relation == "protein") {
      if (input$var2_relation == "protein") {
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
      if (input$var2_relation == "species") {
        # # spec-spec: remove var1 and var2 independently
        # dataHeat$presSpec[dataHeat$supertaxon != in_select
        #                   & dataHeat$var1 < var1_cutoff_min] <- 0
        # dataHeat$presSpec[dataHeat$supertaxon != in_select
        #                   & dataHeat$var1 > var1_cutoff_max] <- 0
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
    dataHeat$supertaxon <- as.factor(dataHeat$supertaxon)
    
    ### add gene categories (if provided)
    if (input$color_by_group == TRUE) {
      # get gene category
      gene_category_file <- input$gene_category
      if (is.null(gene_category_file)) {
        cat_dt <- data.frame(
          geneID = levels(dataHeat$geneID)
        )
        cat_dt$group <- "no_category"
      } else {
        
        cat_dt <- as.data.frame(read.table(file = gene_category_file$datapath,
                                           sep = "\t",
                                           header = TRUE,
                                           check.names = FALSE,
                                           comment.char = "",
                                           fill = TRUE))
        colnames(cat_dt) <- c("geneID","group")
      }
      
      # create a dataframe that contain all genes and all taxa
      dataHeat_cat <- data.frame(
        supertaxon = rep(levels(dataHeat$supertaxon), 
                         nlevels(dataHeat$geneID)),
        geneID = rep(levels(dataHeat$geneID), 
                     each = nlevels(dataHeat$supertaxon))
      )
      
      dataHeat_cat <- merge(dataHeat_cat, cat_dt, by = "geneID")
      
      # add categories into dataHeat
      dataHeat <- merge(dataHeat_cat, dataHeat, 
                        by = c("geneID","supertaxon"), 
                        all.x = TRUE)
    }
    return(dataHeat)
  })
  
  # * clustered heatmap data --------------------------------------------------
  clusteredDataHeat <- reactive({
    dataHeat <- dataHeat()
    dat <- get_profiles()
    cutoff <- (input$var1[1])*100

    # do clustering based on distance matrix
   
    row.order <- hclust(get_distance_matrix_profiles(),
                        method = input$cluster_method)$order
   
    # col.order <- hclust(
    #   get_distance_matrix(t(dat),
    #                       method = input$dist_method, 
    #                       (input$var1[1])*100),
    #   method = input$cluster_method)$order
    
    
    # re-order distance matrix accoring to clustering
    dat_new <- dat[row.order, ] #col.order
    #dat_new <- dat[row.order, col.order]
    
    # return clustered gene ID list
    clustered_gene_ids <- as.factor(row.names(dat_new))
    
    # sort original data according to clustered_gene_ids
    dataHeat$geneID <- factor(dataHeat$geneID,
                              levels = clustered_gene_ids)
    
    dataHeat <- dataHeat[!is.na(dataHeat$geneID),]
    # print("Writing the table...")
    # write.table(
    #   clustered_gene_ids,
    #   file = paste0("../gene_order_",
    #                 input$dist_method, "_", (input$var1[1])*100),
    #   col.names = FALSE, row.names = FALSE, quote = FALSE
    # )
    return(dataHeat)
  })
  
  # =========================== MAIN PROFILE TAB ==============================
  
  # * get total number of genes -----------------------------------------------
  output$total_gene_number.ui <- renderUI({
    geneList <- preData()
    geneList$geneID <- as.factor(geneList$geneID)
    out <- as.list(levels(geneList$geneID))
    
    listIn <- input$list
    if (!is.null(listIn)) {
      list <- as.data.frame(read.table(file = listIn$datapath,
                                       header = FALSE))
      out <- as.list(list$V1)
    }
    
    if (length(out) > 0) {
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
  
  # * get list of genes for highlighting --------------------------------------
  output$highlight_gene_ui <- renderUI({
    geneList <- dataHeat()
    geneList$geneID <- as.factor(geneList$geneID)
    
    out <- as.list(levels(geneList$geneID))
    out <- append("none", out)
    
    selectInput("gene_highlight", "Highlight:", out, selected = out[1])
  })
  
  # * reset configuration windows of Main plot --------------------------------
  observeEvent(input$reset_main_config, {
    shinyjs::reset("x_size")
    shinyjs::reset("y_size")
    shinyjs::reset("legend_size")
    shinyjs::reset("x_angle")
    shinyjs::reset("dot_zoom")
  })
  
  # * close configuration windows of Main plot --------------------------------
  observeEvent(input$apply_main_config, {
    toggleModal(session, "main_plot_config_bs", toggle = "close")
  })
  
  # * parameters for the main profile plot ------------------------------------
  get_parameter_input_main <- reactive({
    input_para <- list(
      "x_axis" = input$x_axis,
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
      "guideline" = 1
    )
    return(input_para)
  })
  
  # * render dot size to dot_size_info ----------------------------------------
  output$dot_size_info <- renderUI({
    if (v$doPlot == FALSE) return()
    
    dataHeat <- dataHeat()
    dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
    presentVl <- dataHeat$presSpec[!is.na(dataHeat$presSpec)]
    
    minDot <- (floor(min(presentVl) * 10) / 10 * 5) * (1 + input$dot_zoom)
    maxDot <- (floor(max(presentVl) * 10) / 10 * 5) * (1 + input$dot_zoom)
    
    em(paste0("current point's size: ", minDot, " - ", maxDot))
  })
  
  # * plot main profile -------------------------------------------------------
  mainpoint_info <- callModule(
    create_profile_plot, "main_profile",
    data = dataHeat,
    clusteredDataHeat = clusteredDataHeat,
    apply_cluster = reactive(input$apply_cluster),
    parameters = get_parameter_input_main,
    in_seq = reactive(input$in_seq),
    in_taxa = reactive(input$in_taxa),
    rank_select = reactive(input$rank_select),
    in_select = reactive(input$in_select),
    taxon_highlight = reactive(input$taxon_highlight),
    gene_highlight = reactive(input$gene_highlight),
    width = reactive(input$width),
    height = reactive(input$height),
    x_axis = reactive(input$x_axis),
    type_profile = reactive("main_profile"),
    color_by_group = reactive(input$color_by_group)
  )
  
  
  # ======================== CUSTOMIZED PROFILE TAB ===========================
  
  # * get list of all sequence IDs for customized profile -----
  output$gene_in <- renderUI({
    filein <- input$main_input
    fileCustom <- input$custom_file
    
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor") {
      filein <- 1
    }
    
    if (is.null(filein) & is.null(fileCustom)) {
      return(selectInput("in_seq", "", "all"))
    }
    if (v$doPlot == FALSE) {
      return(selectInput("in_seq", "", "all"))
    }
    else{
      # full list
      data <- as.data.frame(get_data_filtered())
      data$geneID <- as.character(data$geneID)
      data$geneID <- as.factor(data$geneID)
      outAll <- as.list(levels(data$geneID))
      outAll <- append("all", outAll)
      out <- list()
      if (input$add_gene_age_custom_profile == TRUE) {
        out <- as.list(selectedgene_age())
      } else if (input$add_cluster_cutom_profile == TRUE) {
        out <- as.list(brushed_clusterGene())
      }
      else if (input$add_core_gene_custom_profile == TRUE) {
        out <- as.list(core_geneDf())
      }
      else if (input$add_gc_genes_custom_profile == TRUE) {
        out <- as.list(gene_list_gc())
      }
      else {
        if (!is.null(fileCustom)) {
          customList <- as.data.frame(read.table(file = fileCustom$datapath,
                                                 header = FALSE))
          customList$V1 <- as.factor(customList$V1)
          out <- as.list(levels(customList$V1))
        }
      }
      
      if (length(out) > 0) {
        create_select_gene("in_seq", out, out)
      }
      else {
        create_select_gene("in_seq", outAll, outAll[1])
      }
    }
  })
  
  # * render popup for selecting taxon rank and return list of subset taxa ----
  cus_taxaName <- callModule(
    select_taxon_rank,
    "select_taxon_rank",
    rank_select = reactive(input$rank_select),
    subset_taxa = subset_taxa
  )
  
  # * get list of all taxa for customized profile -----------------------------
  output$taxa_in <- renderUI({
    filein <- input$main_input
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor") {
      filein <- 1
    }
    
    if (is.null(filein)) return(selectInput("in_taxa", "", "all"))
    if (v$doPlot == FALSE) return(selectInput("in_taxa", "", "all"))
    else{
      choice <- alltaxa_list()
      choice$fullName <- as.factor(choice$fullName)
      
      out <- as.list(levels(choice$fullName))
      out <- append("all", out)
      if (input$apply_cus_taxa == TRUE) {
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
  
  # * check if all genes and all species are selected -------------------------
  output$same_profile <- reactive({
    if (v$doPlot == FALSE) return(FALSE)
    if (length(input$in_seq[1]) == 0) return(FALSE)
    else{
      if (input$in_seq[1] == "all" & input$in_taxa[1] == "all") return(TRUE)
    }
  })
  outputOptions(output, "same_profile", suspendWhenHidden = FALSE)
  
  # * change label of plot_custom button for not auto_update ------------------
  output$plot_custom_btn <- renderUI({
    if (input$auto_update_selected == FALSE) {
      shinyBS::bsButton("plot_custom",
                        "Plot/Update selected sequence(s)/taxa",
                        style = "warning")
    } else {
      shinyBS::bsButton("plot_custom",
                        "Plot selected sequence(s)/taxa",
                        style = "warning")
    }
  })
  
  # * check if button (custom)PLOT is clicked ---------------------------------
  vCt <- reactiveValues(doPlotCustom = FALSE)
  observeEvent(input$plot_custom, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    vCt$doPlotCustom <- input$plot_custom
    filein <- input$main_input
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor") {
      filein <- 1
    }
    if (is.null(filein)) vCt$doPlotCustom <- FALSE
  })
  
  # * reset configuration windows of Customized plot -------------------------
  observeEvent(input$reset_selected_config, {
    shinyjs::reset("x_size_select")
    shinyjs::reset("y_size_select")
    shinyjs::reset("legend_size_select")
    shinyjs::reset("x_angle_select")
    shinyjs::reset("dot_zoom_select")
  })
  
  # ** close configuration windows of Customized plot -------------------------
  observeEvent(input$apply_selected_config, {
    toggleModal(session, "selected_plot_config_bs", toggle = "close")
  })
  
  # * parameters for the customized profile plot ------------------------------
  get_parameter_input_customized <- reactive({
    input_para <- list(
      "x_axis" = input$x_axis_selected,
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
      "guideline" = 0
    )
    return(input_para)
  })
  
  # * plot customized profile -------------------------------------------------
  selectedpoint_info <- callModule(
    create_profile_plot, "customized_profile",
    data = dataHeat,
    clusteredDataHeat = clusteredDataHeat,
    apply_cluster = reactive(input$apply_cluster),
    parameters = get_parameter_input_customized,
    in_seq = reactive(input$in_seq),
    in_taxa = reactive(input$in_taxa),
    rank_select = reactive(input$rank_select),
    in_select = reactive(input$in_select),
    taxon_highlight = reactive("none"),
    gene_highlight = reactive("none"),
    width = reactive(input$selected_width),
    height = reactive(input$selected_height),
    x_axis = reactive(input$x_axis_selected),
    type_profile = reactive("customized_profile"),
    color_by_group = reactive(input$color_by_group)
  )
  
  # ============================== POINT INFO =================================
  
  # * get status of point_info for activating Detailed Plot button ------------
  output$point_info_status <- reactive({
    if (input$tabs == "Main profile") {
      # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
      info <- mainpoint_info()
    } else if (input$tabs == "Customized profile") {
      info <- selectedpoint_info()
    } else {
      info <- NULL
    }
    is.null(info)
  })
  outputOptions(output, "point_info_status", suspendWhenHidden = FALSE)
  
  # * show info into "point's info" box ---------------------------------------
  output$point_info <- renderText({
    # GET INFO BASED ON CURRENT TAB
    if (input$tabs == "Main profile") {
      # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
      info <- mainpoint_info()
    } else if (input$tabs == "Customized profile") {
      info <- selectedpoint_info()
    } else {
      return()
    }
    
    if (is.null(info)) return()
    else{
      orthoID <- info[2]
      
      if (is.na(orthoID)) return()
      else{
        # if (orthoID=="NA") {orthoID <- info[2]}
        ## print output
        a <- toString(paste("Seed-ID:", info[1]))
        b <- toString(paste0("Hit-ID: ",
                             orthoID,
                             " (",
                             substr(info[3], 6, nchar(info[3])),
                             ")"))
        c <- ""
        if (input$var1_id != "") {
          c <- toString(paste(input$var1_aggregate_by,
                              input$var1_id,
                              ":",
                              info[4]))
        }
        d <- ""
        if (input$var2_id != "") {
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
  
  # * data for detailed plot --------------------------------------------------
  detail_plotDt <- reactive({
    if (v$doPlot == FALSE) return()
    
    # GET INFO BASED ON CURRENT TAB
    if (input$tabs == "Main profile") {
      # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
      info <- mainpoint_info()
    } else if (input$tabs == "Customized profile") {
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
      if (input$detailed_remove_na == TRUE) {
        joinedDf <- joinedDf[!is.na(joinedDf$orthoID), ]
      }
      
      ### return data for detailed plot
      return(joinedDf)
    }
  })
  
  # * render detailed plot ----------------------------------------------------
  
  point_infoDetail <- callModule(
    create_detailed_plot, "detailed_plot",
    data = detail_plotDt,
    var1_id = reactive(input$var1_id),
    var2_id = reactive(input$var2_id),
    detailed_text = reactive(input$detailed_text),
    detailed_height = reactive(input$detailed_height)
  )
  
  # * render FASTA sequence ---------------------------------------------------
  output$fasta <- renderText({
    if (v$doPlot == FALSE) return()
    
    info <- point_infoDetail() # info = seedID, orthoID, var1
    
    if (is.null(info)) return()
    else{
      data <- get_data_filtered()
      
      seqID <- toString(info[2])
      groupID <- toString(info[1])
      ncbiID <- gsub("ncbi", "", toString(info[5]))
      
      if (input$demo_data == "none") {
        filein <- input$main_input
        input_type <- check_input_vadility(filein)
      } else {
        input_type <- "demo"
      }
      
      if (input_type == "oma") {
        fastaOut <- get_selected_fasta_oma(final_oma_df(), seqID)
      } else {
        seqDf <- data.frame("geneID" = groupID,
                            "orthoID" = seqID,
                            "ncbiID" = ncbiID)
        
        fastaOut <- get_fasta_seqs(
          seqDf, input$main_input, input$demo_data,
          input$input_type, input$concat_fasta,
          input$path,
          input$dir_format,
          input$file_ext,
          input$id_format,
          get_main_input()
        )
      }
      
      return(paste(fastaOut[1]))
    }
  })
  
  # ======================== FEATURE ARCHITECTURE PLOT ========================
  
  # * get domain file/path ----------------------------------------------------
  getDomainFile <- reactive({
    # click info
    info <- point_infoDetail() # info = seedID, orthoID, var1
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    var1 <- as.character(info[3])
    
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor") {
      if (is.null(info)) {
        fileDomain <- "noSelectHit"
        updateButton(session, "do_domain_plot", disabled = TRUE)
      } else {
        updateButton(session, "do_domain_plot", disabled = FALSE)
        if (input$demo_data == "lca-micros") {
          fileDomain <-
            suppressWarnings(
              paste0(
                "https://github.com/BIONF/phyloprofile-data/blob/master/demo/domain_files/",
                group,
                ".domains?raw=true"
              )
            )
        } else if (input$demo_data == "ampk-tor") {
          fileDomain <-
            suppressWarnings(
              paste0(
                "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.domains_F"
              )
            )
        }
      }
    } else {
      filein <- input$main_input
      if (check_input_vadility(filein) == "oma") {
        fileDomain <- "oma_input"
        updateButton(session, "do_domain_plot", disabled = FALSE)
      } else if (input$anno_location == "from file") {
        fileDomain <- input$file_domain_input
        if (is.null(fileDomain)) {
          fileDomain <- "noFileInput"
        } else {
          if (is.null(info)) {
            fileDomain <- "noSelectHit"
            updateButton(session, "do_domain_plot", disabled = TRUE)
          } else {
            updateButton(session, "do_domain_plot", disabled = FALSE)
            fileDomain <- fileDomain$datapath
          }
        }
      } else {
        if (is.null(info)) {
          fileDomain <- "noSelectHit"
          updateButton(session, "do_domain_plot", disabled = TRUE)
        } else {
          # check file extension
          allExtension <- c("txt", "csv", "list", "domains", "architecture")
          flag <- 0
          for (i in 1:length(allExtension)) {
            fileDomain <- paste0(input$domainPath,
                                 "/",
                                 group,
                                 ".",
                                 allExtension[i])
            if (file.exists(fileDomain) == TRUE) {
              updateButton(session, "do_domain_plot", disabled = FALSE)
              flag <- 1
              break()
            }
          }
          
          if (flag == 0) {
            fileDomain <- "noFileInFolder"
            updateButton(session, "do_domain_plot", disabled = TRUE)
          }
        }
      }
    }
    return(fileDomain)
  })
  
  # * check domain file -------------------------------------------------------
  output$check_domain_files <- renderUI({
    fileDomain <- getDomainFile()
    if (fileDomain == "noFileInput") {
      em("Domain file not provided!!")
    } else if (fileDomain == "noFileInFolder") {
      msg <- paste0(
        "<p><em>Domain file not found!! </em></p>
        <p><em>Please make sure that file name has to be in this format:
        <strong>&lt;seedID&gt;.extension</strong>, where extension is limited
        to <strong>txt</strong>, <strong>csv</strong>, <strong>list</strong>,
        <strong>domains</strong> or <strong>architecture</strong>.
        </em></p>"
      )
      HTML(msg)
    } else if (fileDomain == "noSelectHit") {
      em("Please select one ortholog sequence!!")
    }
  })
  
  # * render domain plot ------------------------------------------------------
  observeEvent(input$do_domain_plot, {
    callModule(
      create_architecture_plot, "archi_plot",
      point_info = point_infoDetail,
      domain_info = get_domain_information,
      label_archi_size = reactive(input$label_archi_size),
      title_archi_size = reactive(input$title_archi_size),
      archi_height = reactive(input$archi_height),
      archi_width = reactive(input$archi_width)
    )
  })
  
  # ======================== FILTERED DATA DOWNLOADING ========================
  
  # * for main profile ========================================================
  main_fasta_download <- reactive({
    if (input$demo_data == "none") {
      filein <- input$main_input
      input_type <- check_input_vadility(filein)
    } else {
      input_type <- "demo"
    }
    
    if (input_type == "oma") {
      all_oma_df <- final_oma_df()
      filtered_download_df <- as.data.frame(download_data())
      filtered_oma_df <-
        subset(all_oma_df,
               all_oma_df$ortho_id %in% filtered_download_df$orthoID &
                 all_oma_df$seed %in% filtered_download_df$geneID)
      main_fasta_out <- get_all_fasta_oma(filtered_oma_df)
    } else {
      main_fasta_out <- get_fasta_seqs(
        as.data.frame(download_data()),
        input$main_input, input$demo_data,
        input$input_type, input$concat_fasta,
        input$path,
        input$dir_format,
        input$file_ext,
        input$id_format,
        get_main_input()
      )
    }
    
    return(main_fasta_out)
  })
  
  download_data <- callModule(
    download_filtered_main,
    "filtered_main_download",
    data = get_data_filtered,
    fasta = main_fasta_download,
    var1_id = reactive(input$var1_id),
    var2_id = reactive(input$var2_id),
    var1 = reactive(input$var1),
    var2 = reactive(input$var2),
    percent = reactive(input$percent)
  )
  
  # * for customized profile ==================================================
  customized_fasta_download <- reactive({
    if (input$demo_data == "none") {
      filein <- input$main_input
      input_type <- check_input_vadility(filein)
    } else {
      input_type <- "demo"
    }
    
    if (input_type == "oma") {
      all_oma_df <- final_oma_df()
      filtered_download_df <- as.data.frame(download_custom_data())
      filtered_oma_df <-
        subset(all_oma_df,
               all_oma_df$ortho_id %in% filtered_download_df$orthoID &
                 all_oma_df$seed %in% filtered_download_df$geneID)
      fasta_out_df <- get_all_fasta_oma(filtered_oma_df)
    } else {
      fasta_out_df <- get_fasta_seqs(
        as.data.frame(download_custom_data()),
        input$main_input, input$demo_data,
        input$input_type, input$concat_fasta,
        input$path,
        input$dir_format,
        input$file_ext,
        input$id_format,
        get_main_input()
      )
    }
    return(fasta_out_df)
  })
  download_custom_data <- callModule(
    download_filtered_customized,
    "filtered_customized_download",
    data = download_data,
    fasta = customized_fasta_download,
    in_seq = reactive(input$in_seq),
    in_taxa = reactive(input$in_taxa)
  )
  
  # ============================ ANALYSIS FUNCTIONS ===========================
  
  # * PROFILE CLUSTERING ======================================================
  
  # ** check if genes are added anywhere else to the customized profile -------
  observe({
    if (input$add_gene_age_custom_profile == TRUE
        | input$add_core_gene_custom_profile == TRUE
        | input$add_gc_genes_custom_profile == TRUE) {
      shinyjs::disable("add_cluster_cutom_profile")
    }else{
      shinyjs::enable("add_cluster_cutom_profile")
    }
  })
  
  output$add_cluster_cutom_profile_check.ui <- renderUI({
    if (input$add_gene_age_custom_profile == TRUE
        | input$add_core_gene_custom_profile == TRUE |
        input$add_gc_genes_custom_profile == TRUE ) {
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in
           <strong>Gene age estimation</strong> or
           <strong>Core genes finding</strong> or
           <strong>Group comparison</strong>
           &nbsp;to enable this function)</em></p>')
    }
    })
  
  # ** List of possible profile types -----------------------------------------
  output$select_profile_type <- renderUI({
    variable1 <- paste0("profile using ", input$var1_id)
    if (input$var2_id != "") {
      variable2 <- paste0("profile using ", input$var2_id)
      radioButtons(
        "profile_type",
        label = h5("Select the profile type"),
        choiceNames = list(
          "binary profile",
          variable1,
          variable2),
        choiceValues = list(
          "binary", "var1", "var2"
        ),
        selected = "binary",
        inline = FALSE)
    }
    else {
      radioButtons(
        "profile_type",
        label = h5("Select the profile type"),
        choiceNames = list(
          "binary profile",
          variable1),
        choiceValues = list(
          "binary", "var1"
        ),
        selected = "binary",
        inline = FALSE)
    }

  })
  
  # ** List of possible distance methods --------------------------------------
  output$select_dist_method <- renderUI({
    if (is.null(input$profile_type)) return()
    
    if (input$profile_type == "binary") {
    selectInput(
      "dist_method",
      label = h5("Distance measure method:"),
      choices = list("euclidean" = "euclidean",
                     "maximum" = "maximum",
                     "manhattan" = "manhattan",
                     "canberra" = "canberra",
                     "binary" = "binary",
                     "pearson correlation coefficient" = "pearson",
                     "mutual information" = "mutual_information",
                     "distance correlation" = "distance_correlation"
      ),
      selected = "euclidean"
    )
    }
    else {
      selectInput(
        "dist_method",
        label = h5("Distance measure method:"),
        choices = list("mutual information" = "mutual_information",
                       "distance correlation" = "distance_correlation"
        ),
        selected = "mutual_information"
      )
    }
  })
  
  # ** Distance matrix --------------------------------------------------------
  get_distance_matrix_profiles <- reactive({
    if (is.null(input$dist_method)) return()
    profiles <- get_profiles()
    distance_matrix <- get_distance_matrix(profiles, input$dist_method, (input$var1[1])*100)
    # write.table(as.matrix(distance_matrix), file = paste0("../distance_matrix_",
    #                                                       input$dist_method, "_",(input$var1[1])*100 ),
    #             col.names = TRUE, row.names = TRUE, quote = FALSE, sep = " \t")
    return(distance_matrix)
  })
  
  # ** Phylogenetic profiles --------------------------------------------------
  get_profiles <- reactive({
    data_heat <- dataHeat()
    if (nrow(data_heat) < 1) return()
    if (is.null(input$dist_method)) return()
    profiles <- get_data_clustering(data_heat,
                                    input$profile_type,
                                    input$var1_aggregate_by,
                                    input$var2_aggregate_by)
    return(profiles)
  })
  
  
  
  # ** render cluster tree ---------------------------------------------------
  brushed_clusterGene <- callModule(
    cluster_profile, "profile_clustering",
    distance_matrix = get_distance_matrix_profiles,
    cluster_method = reactive(input$cluster_method),
    plot_width = reactive(input$cluster_plot.width),
    plot_height = reactive(input$cluster_plot.height))
  
  
  # * DISTRIBUTION ANALYSIS ===================================================
  
  # ** list of available variables for distribution plot ----------------------
  output$selected.distribution <- renderUI({
    if (nchar(input$var1_id) == 0 & nchar(input$var2_id) == 0) {
      varList <- "% present taxa"
    } else if (nchar(input$var1_id) == 0 & nchar(input$var2_id) > 0) {
      varList <- as.list(c(input$var2_id, "% present taxa"))
    } else if (nchar(input$var1_id) > 0 & nchar(input$var2_id) == 0) {
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
  
  # ** var1 / var2 distribution data ------------------------------------------
  distribution_df <- reactive({
    if (v$doPlot == FALSE) return()
    
    dataOrig <- get_main_input()
    if (ncol(dataOrig) < 4) {
      colnames(dataOrig) <- c("geneID",
                              "ncbiID",
                              "orthoID")
      splitDt <- dataOrig[, c("orthoID")]
    } else if (ncol(dataOrig) < 5) {
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
    
    if (length(levels(as.factor(splitDt$var2))) == 1) {
      if (levels(as.factor(splitDt$var2)) == "") {
        splitDt$var2 <- 0
      }
    }
    
    # convert factor into numeric for "var1" & "var2" column
    if ("var1" %in% colnames(splitDt)) {
      splitDt$var1 <- suppressWarnings(as.numeric(as.character(splitDt$var1)))
      # filter splitDt based on selected var1 cutoff
      splitDt <- splitDt[splitDt$var1 >= input$var1[1]
                         & splitDt$var1 <= input$var1[2], ]
    }
    if ("var2" %in% colnames(splitDt)) {
      splitDt$var2 <- suppressWarnings(as.numeric(as.character(splitDt$var2)))
      # filter splitDt based on selected var2 cutoff
      splitDt <- splitDt[splitDt$var2 >= input$var2[1]
                         & splitDt$var2 <= input$var2[2], ]
    }
    
    # filter data base on customized plot (if chosen)
    if (input$dataset.distribution == "Customized data") {
      # get geneID and supertaxon name for splitDt
      allData <- get_data_filtered()
      splitDtName <- merge(splitDt, allData,
                           by = "orthoID",
                           all.x = TRUE)
      splitDtName$supertaxonMod <-
        substr(splitDtName$supertaxon,
               6,
               nchar(as.character(splitDtName$supertaxon)))
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
      if (input$in_taxa[1] == "all" & input$in_seq[1] != "all") {
        # select data from dataHeat for selected sequences only
        splitDt <- subset(splitDtName, geneID %in% input$in_seq)
      } else if (input$in_seq[1] == "all" & input$in_taxa[1] != "all") {
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
  
  # ** calculate % present species in supertaxa -------------------------------
  presSpecAllDt <- reactive({
    # open main input file
    mdData <- get_main_input()
    if (ncol(mdData) < 4) {
      colnames(mdData) <- c("geneID",
                            "ncbiID",
                            "orthoID")
    } else if (ncol(mdData) < 5) {
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
    if ("var1" %in% colnames(taxaMdData)) {
      taxaMdData$var1 <-
        suppressWarnings(as.numeric(as.character(taxaMdData$var1)))
    }
    if ("var2" %in% colnames(taxaMdData)) {
      taxaMdData$var2 <-
        suppressWarnings(as.numeric(as.character(taxaMdData$var2)))
    }
    # calculate % present species
    finalPresSpecDt <- calc_pres_spec(taxaMdData, taxaCount)
    finalPresSpecDt
  })
  
  # ** render distribution plots ----------------------------------------------
  observe({
    if (v$doPlot == FALSE) return()
    
    if (is.null(input$selected_dist)) {
      return()
    } else {
      if (input$selected_dist == "% present taxa") {
        callModule(
          analyze_distribution, "dist_plot",
          data = presSpecAllDt,
          var_id = reactive(input$selected_dist),
          var_type = reactive("presSpec"),
          percent = reactive(input$percent),
          dist_text_size = reactive(input$dist_text_size),
          dist_width = reactive(input$dist_width)
        )
      } else{
        if (input$selected_dist == input$var1_id) {
          callModule(
            analyze_distribution, "dist_plot",
            data = distribution_df,
            var_id = reactive(input$selected_dist),
            var_type = reactive("var1"),
            percent = reactive(input$percent),
            dist_text_size = reactive(input$dist_text_size),
            dist_width = reactive(input$dist_width)
          )
        } else if (input$selected_dist == input$var2_id) {
          callModule(
            analyze_distribution, "dist_plot",
            data = distribution_df,
            var_id = reactive(input$selected_dist),
            var_type = reactive("var2"),
            percent = reactive(input$percent),
            dist_text_size = reactive(input$dist_text_size),
            dist_width = reactive(input$dist_width)
          )
        }
      }
    }
  })
  
  # * GENE AGE ESTIMATION =====================================================
  
  # ** check if genes are added anywhere else to the customized profile -------
  observe({
    if (input$add_cluster_cutom_profile == TRUE
        | input$add_core_gene_custom_profile == TRUE
        | input$add_gc_genes_custom_profile == TRUE ) {
      shinyjs::disable("add_gene_age_custom_profile")
    } else {
      shinyjs::enable("add_gene_age_custom_profile")
    }
  })
  
  output$add_gene_age_custom_profile_check.ui <- renderUI({
    if (input$add_cluster_cutom_profile == TRUE
        | input$add_core_gene_custom_profile == TRUE
        | input$add_gc_genes_custom_profile == TRUE) {
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in
           <strong>Profile clustering</strong> or
           <strong>Core genes finding</strong> or
           <strong>Group comparison</strong>
           &nbsp;to enable this function)</em></p>')
    }
    })
  
  # ** reset gene_age_prot_config ---------------------------------------------
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
  selectedgene_age <- callModule(
    plot_gene_age, "gene_age",
    data = gene_ageDf,
    gene_age_width = reactive(input$gene_age_width),
    gene_age_height = reactive(input$gene_age_height),
    gene_age_text = reactive(input$gene_age_text)
  )
  
  # * CORE GENES IDENTIFICATION ===============================================
  
  # ** render list of available taxa ------------------------------------------
  output$taxa_list_core.ui <- renderUI({
    filein <- input$main_input
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor") {
      filein <- 1
    }
    if (is.null(filein)) {
      return(selectInput("in_taxa",
                         "Select taxa of interest:",
                         "none"))
    }
    if (v$doPlot == FALSE) {
      return(selectInput("in_taxa",
                         "Select taxa of interest:",
                         "none"))
    }
    else{
      choice <- alltaxa_list()
      choice$fullName <- as.factor(choice$fullName)
      
      out <- as.list(levels(choice$fullName))
      out <- append("none", out)
      
      if (input$apply_core_taxa == TRUE) {
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
  core_taxa_name <- callModule(
    select_taxon_rank,
    "select_taxon_rank_core",
    rank_select = reactive(input$rank_select),
    subset_taxa = subset_taxa
  )
  
  # ** check if genes are added anywhere else to the customized profile -------
  observe({
    if (input$add_cluster_cutom_profile == TRUE
        | input$add_gene_age_custom_profile == TRUE
        | input$add_gc_genes_custom_profile == TRUE) {
      shinyjs::disable("add_core_gene_custom_profile")
    } else {
      shinyjs::enable("add_core_gene_custom_profile")
    }
  })
  
  output$add_core_gene_custom_profile_check.ui <- renderUI({
    if (input$add_cluster_cutom_profile == TRUE
        | input$add_gene_age_custom_profile == TRUE
        | input$add_gc_genes_custom_profile == TRUE) {
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in
           <strong>Profiles clustering</strong> or
           <strong>Gene age estimating</strong> or
           <strong>Group Comparioson</strong>
           &nbsp;to enable this function)</em></p>')
    }
    })
  
  
  
  # ** render table contains list of core genes -------------------------------
  core_geneDf <- callModule(
    identify_core_gene,
    "core_gene",
    filtered_data = get_data_filtered,
    rank_select = reactive(input$rank_select),
    taxa_core = reactive(input$taxa_core),
    percent_core = reactive(input$percent_core),
    core_coverage = reactive(input$core_coverage)
  )
  
  # ** download gene list from core_gene.table --------------------------------
  output$core_gene_table_download <- downloadHandler(
    filename = function() {
      c("coreGeneList.out")
    },
    content = function(file) {
      data_out <- core_geneDf()
      write.table(data_out, file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
  # * GET NCBI TAXONOMY IDs FROM INPUT LIST OF TAXON NAMES ====================
  callModule(search_taxon_id,"search_taxon_id")
  
  # * GROUP COMPARISON ========================================================
  # ** list of all available genes --------------------------------------------
  output$list_genes_gc <- renderUI({
    filein <- input$main_input
    
    file_gc <- input$gc_file
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor") {
      filein <- 1
    }
    
    
    if (v$doPlot == FALSE) {
      return(selectInput("list_selected_genes_gc", "Select sequence(s):", "none"))
    }
    
    if (is.null(filein) & is.null(file_gc)) {
      return(selectInput("list_selected_genes_gc", "Select sequence(s):", "none"))
    } else {
      # full list
      data <- as.data.frame(get_data_filtered())
      data$geneID <- as.character(data$geneID)
      data$geneID <- as.factor(data$geneID)
      out_all <- as.list(levels(data$geneID))
      out_all <- append("all", out_all)
      if (is.null(file_gc)) {
        selectInput("list_selected_genes_gc", "Select sequence(s):",
                    out_all,
                    selected = out_all[1],
                    multiple = TRUE,
                    selectize = FALSE)
      } else {
        list_gc <- as.data.frame(read.table(file = file_gc$datapath,
                                            header = FALSE))
        list_gc$V1 <- as.factor(list_gc$V1)
        out <- as.list(levels(list_gc$V1))
        selectInput("list_selected_genes_gc", "Select sequence(s):",
                    out,
                    selected = NULL,
                    multiple = FALSE,
                    selectize = FALSE)
      }
    }
  })
  
  # ** popup for selecting taxon rank and return list of belonging taxa -------
  gc_taxa_name <- callModule(
    select_taxon_rank,
    "select_taxon_rank_gc",
    rank_select = reactive(input$rank_select),
    subset_taxa = subset_taxa
  )
  
  # ** list of available taxa (for selecting as in_group) ---------------------
  output$taxa_list_gc <- renderUI({
    filein <- input$main_input
    if (input$demo_data == "lca-micros" | input$demo_data == "ampk-tor") {
      filein <- 1
    }
    if (is.null(filein)) {
      return(selectInput("in_taxa", "Select in_group taxa:", "none"))
    }
    if (v$doPlot == FALSE) {
      return(selectInput("in_taxa", "Select in_group taxa:", "none"))
    } else {
      choice <- alltaxa_list()
      choice$fullName <- as.factor(choice$fullName)
      
      out <- as.list(levels(choice$fullName))
      
      #' when the taxonomy rank was changed -----------------------------------
      if (input$apply_taxa_gc == TRUE) {
        out <- gc_taxa_name()
        selectInput("selected_in_group_gc", "Select in_group taxa:",
                    out,
                    selected = out,
                    multiple = TRUE,
                    selectize = FALSE)
      }
      #' when the taxonomy is the same as the initially chosen one ------------
      else {
        #' check for the rank of the rank in the input
        ranks <- get_taxonomy_ranks()
        pos <- which(ranks == input$rank_select) # position in the list
        higher_rank <- ranks[pos + 1] # take the next higher rank
        higher_rank_name <- substr(higher_rank, 4, nchar(higher_rank))
        
        name_list <- get_name_list(TRUE, TRUE) # get the taxon names
        dt <- get_taxa_list(FALSE, subset_taxa) # get the taxa
        
        #' get the info for the reference protein from the namelist
        reference <- subset(name_list, name_list$fullName == input$in_select)
        
        #' get the id for every rank for the reference protein
        rank_name <- substr(input$rank_select, 4, nchar(input$rank_select))
        reference_dt <- dt[dt[, rank_name] == reference$ncbiID, ]
        
        #' save the next higher rank
        reference_higher_rank <- reference_dt[higher_rank_name]
        reference_higher_rank <-
          reference_higher_rank[!duplicated(reference_higher_rank), ]
        
        #' get all the taxa with the same id in the next higher rank
        selected_taxa_dt <-
          subset(dt, dt[, higher_rank_name] %in% reference_higher_rank)
        selected_taxa_dt <-
          selected_taxa_dt[!duplicated(selected_taxa_dt[rank_name]), ]
        
        #' get a list with all the ids with reference_higher_rank as parent
        selected_taxa_ids <- selected_taxa_dt[rank_name]
        if (length(selected_taxa_ids[[1]]) >= 1) {
          selected_taxa_ids <- selected_taxa_ids[[1]]
        }
        
        selected_taxa <- subset(name_list, name_list$rank == rank_name)
        selected_taxa <-
          subset(selected_taxa, selected_taxa$ncbiID %in% selected_taxa_ids)
        
        default_select <- selected_taxa$fullName
        
        selectInput("selected_in_group_gc", "Select in_group taxa:",
                    out,
                    selected = default_select,
                    multiple = TRUE,
                    selectize = FALSE)
      }
    }
  })
  
  # ** buttons to choose the variable -----------------------------------------
  output$variable_button_gc <- renderUI({
    radioButtons(
      inputId = "var_name_gc",
      label = "Select variable(s) to compare:",
      choices = list(input$var1_id, input$var2_id, "Both"),
      selected = input$var1_id,
      inline = FALSE
    )
  })
  
  # ** slider to set significance level ---------------------------------------
  output$significance.ui <- renderUI({
    msg <- paste0(
      "P-value cut-off of the statistic test"
    )
    popify(
      sliderInput(
        "significance",
        paste("Significance level:"),
        min = 0,
        max = 1,
        step = 0.005,
        value = c(0.05),
        width = 200
      ),
      "",
      msg
    )
  })
  
  # ** reset plots config -----------------------------------------------------
  observeEvent(input$reset_config_gc, {
    shinyjs::reset("x_size_gc")
    shinyjs::reset("y_size_gc")
    shinyjs::reset("angle_gc")
    shinyjs::reset("legend_size_gc")
  })
  
  observeEvent(input$apply_config_gc, {
    toggleModal(session, "gc_plot_config_bs", toggle = "close")
  })
  
  # ** check if genes are added anywhere else to the customized profile -------
  observe({
    if (input$add_gene_age_custom_profile == TRUE |
        input$add_core_gene_custom_profile == TRUE |
        input$add_cluster_cutom_profile == TRUE) {
      shinyjs::disable("add_gc_genes_custom_profile")
    }else{
      shinyjs::enable("add_gc_genes_custom_profile")
    }
  })
  
  output$add_gc_custom_profile_check <- renderUI({
    if (input$add_gene_age_custom_profile == TRUE |
        input$add_core_gene_custom_profile == TRUE |
        input$add_cluster_cutom_profile == TRUE) {
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in
           <strong>Gene age estimation</strong> or
           <strong>Profile clustering</strong> or
           <strong>Core genes finding</strong>&nbsp;to enable this function)
           </em></p>')
    }
  })
  
  # ** parameters for the plots in Group Comparison ---------------------------
  get_parameter_input_gc <- reactive({
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
                       "legend_size_gc" = input$legend_size_gc,
                       "p_values_size" = input$p_values_size_gc)
    
  })
  
  # ** render plots for group comparison --------------------------------------
  gene_list_gc <- callModule(
    group_comparison, "group_comparison",
    selected_in_group = reactive(input$selected_in_group_gc),
    selected_genes_list = reactive(input$list_selected_genes_gc),
    main_rank = reactive(input$rank_select),
    selected_variable = reactive(input$var_name_gc),
    use_common_ancestor = reactive(input$use_common_ancestor),
    reference_taxon = reactive(input$in_select),
    ncbi_id_list = subset_taxa,
    filtered_data = get_data_filtered,
    right_format_features = reactive(input$right_format_features),
    domain_information = get_domain_information,
    plot = reactive(input$plot_gc),
    parameter = get_parameter_input_gc,
    selected_point = reactive(input$show_point_gc)
  )
    })
