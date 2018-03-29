# NOT WORK WITH shinyapp.io ===================================================
# if(!("pacman" %in% installed.packages())) install.packages("pacman")
# library(pacman)
# p_load(shiny,shinyBS,ggplot2,reshape2,plyr,dplyr,tidyr,scales,grid,gridExtra,ape,stringr,gtable,dendextend,ggdendro,gplots,data.table,taxize,install=T)

# packages <- c("shiny","shinyBS","ggplot2","reshape2","plyr","dplyr","tidyr","scales","grid","gridExtra","ape","stringr","gtable","dendextend","ggdendro","gplots","data.table","taxize","Biostrings","zoo","RCurl")
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
  if("devtools" %in% installed.packages() == FALSE){
    install.packages("devtools")
  }
  devtools::install_github('andrewsali/shinycssloaders', force = TRUE)
}
if(!require("Matching")) install.packages("Matching")
source("scripts/taxonomyProcessing.R")
source("scripts/functions.R")

# MAIN ========================================================================
options(shiny.maxRequestSize=99*1024^2)  # size limit for input 99mb

shinyServer(function(input, output, session) {
  # Automatically stop a Shiny app when closing the browser tab ---------------
  # session$onSessionEnded(stopApp) 
  session$allowReconnect(TRUE)
  
  # ===========================================================================
  # PRE-PROCESSING  ===========================================================
  # ===========================================================================
  
  # check for internet connection  
  observe({
    if(hasInternet() == FALSE){
      # toggleState("demo")
      toggleState("demo_data")
    }
  })
  output$no_internet_msg <- renderUI({
    if(hasInternet() == FALSE){
      strong(em("Internet connection is required for using demo data!"),style = "color:red")
    } else {
      return()
    }
  })
  
  ####### check for the existence of taxonomy info file #######
  observe({
    if(!file.exists(isolate({"data/rankList.txt"}))){
      if(hasInternet() == TRUE){
        ncol <- max(count.fields("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/rankList.txt", sep = '\t'))
        df <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/rankList.txt", sep="\t", quote='', header=F, fill=T, na.strings=c("","NA"), col.names=paste0('V', seq_len(ncol)))
        write.table(df, file ="data/rankList.txt", col.names = F, row.names = F, quote = F, sep="\t")#na = "",
      } else {
        file.create("data/rankList.txt")
      }
    }
  })
  
  observe({
    if(!file.exists(isolate({"data/idList.txt"}))){
      if(hasInternet() == TRUE){
        ncol <- max(count.fields("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/idList.txt", comment.char="", sep = '\t'))
        df <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/idList.txt", sep="\t", header=F, fill=T,comment.char="", na.strings=c("","NA"), col.names=paste0('V', seq_len(ncol)))
        write.table(df, file ="data/idList.txt", col.names = F, row.names = F, quote = F, sep="\t")#na = "",
      } else {
        file.create("data/idList.txt")
      }
    }
  })
  
  observe({
    if(!file.exists(isolate({"data/taxonNamesReduced.txt"}))){
      if(hasInternet() == TRUE){
        ncol <- max(count.fields("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/taxonNamesReduced.txt", sep = '\t'))
        df <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/taxonNamesReduced.txt", sep="\t", quote='', header=F, fill=T, na.strings=c("","NA"), col.names=paste0('V', seq_len(ncol)))
        write.table(df, file ="data/taxonNamesReduced.txt", col.names = F, row.names = F, quote = F, sep="\t")
      } else {
        system("cp data/newTaxa.txt data/taxonNamesReduced.txt")
      }
    }
  })
  
  observe({
    if(!file.exists(isolate({"data/taxonomyMatrix.txt"}))){
      if(hasInternet() == TRUE){
        ncol <- max(count.fields("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/taxonomyMatrix.txt", sep = '\t'))
        df <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/taxonomyMatrix.txt", sep="\t", quote='', header=F, fill=T, na.strings=c("","NA"), col.names=paste0('V', seq_len(ncol)))
        write.table(df, file ="data/taxonomyMatrix.txt", col.names = F, row.names = F, quote = F, sep="\t")
      }
    }
  })
  
  ################## get ncbi taxa IDs ####################
  ##### retrieve ID for list of taxa names
  taxa_id <- reactive({
    if(input$id_search > 0){
      taxain <- input$taxa_list
      if(is.null(taxain)){return()}
      
      taxaNameDf <- as.data.frame(read.table(file=taxain$datapath, sep='\t',header=F,check.names=FALSE,comment.char=""))
      
      idDf <- data.frame("name"=character(),"new_name"=character(),"id"=character(),"type"=character(),stringsAsFactors=FALSE)
      
      withProgress(message = 'Retrieving IDs...', value = 0,{
        for(i in 1:nrow(taxaNameDf)){
          id <- get_uid(sciname = taxaNameDf[i,])[1]
          if(is.na(id)){
            temp <- gnr_resolve(names = as.character(taxaNameDf[i,]))
            if(nrow(temp) > 0){
              new_id <- get_uid(sciname = temp[1,3])[1]
              if(is.na(new_id)){
                idDf[i,] <- c(as.character(taxaNameDf[i,]),as.character(temp[1,3]),paste0("NA"),"notfound")
              } else {
                idDf[i,] <- c(as.character(taxaNameDf[i,]),as.character(temp[1,3]),paste0("ncbi",new_id),"notfound")
              }
            } else {
              idDf[i,] <- c(as.character(taxaNameDf[i,]),paste0("no alternative"),paste0("NA"),"notfound")
            }
          } else {
            idDf[i,] <- c(as.character(taxaNameDf[i,]),"NA",paste0("ncbi",id),"retrieved")
          }
          # Increment the progress bar, and update the detail text.
          incProgress(1/nrow(taxaNameDf), detail = paste(i,"/",nrow(taxaNameDf)))
        }
      })
      
      ### return
      idDf
    }
  })
  
  ### output mismatched taxa
  output$not_found_taxa <- renderDataTable(option = list(searching = FALSE),{
    if(input$id_search > 0){
      if(length(taxa_id())>0){
        tb <- as.data.frame(taxa_id())
        tbFiltered <- tb[tb$type == "notfound",]
        notFoundDt <- tbFiltered[,c("name","new_name","id")]
        colnames(notFoundDt) <- c("Summitted name","Alternative name","Alternative ID")
        notFoundDt
      }
    }
  })
  
  output$download_not_found_taxa <- downloadHandler(
    filename = function(){c("mismatchedTaxa.txt")},
    content = function(file){
      tb <- as.data.frame(taxa_id())
      tbFiltered <- tb[tb$type == "notfound",]
      notFoundDt <- tbFiltered[,c("name","new_name","id")]
      colnames(notFoundDt) <- c("Summitted name","Alternative name","Alternative ID")
      
      write.table(notFoundDt,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ### output retrieved taxa IDs
  output$taxa_id <- renderDataTable(option = list(searching = FALSE),{
    if(input$id_search > 0){
      if(length(taxa_id())>0){
        tb <- as.data.frame(taxa_id())
        tbFiltered <- tb[tb$type == "retrieved",]
        retrievedDt <- tbFiltered[,c("name","id")]
        colnames(retrievedDt) <- c("Taxon_name","Taxon_ID")
        retrievedDt
      }
    }
  })
  
  output$download_taxa_id <- downloadHandler(
    filename = function(){c("retrievedtaxa_id.txt")},
    content = function(file){
      tb <- as.data.frame(taxa_id())
      tbFiltered <- tb[tb$type == "retrieved",]
      retrievedDt <- tbFiltered[,c("name","id")]
      colnames(retrievedDt) <- c("Taxon name","Taxon ID")
      
      write.table(retrievedDt,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ################# PARSING VARIABLE 1 AND 2 ##################
  ######## render textinput for variable 1 & 2
  output$var1_id.ui <- renderUI({
    # if(input$demo == TRUE){
    #   textInput("var1_id", h5("First variable:"), value = "Domain similarity", width="100%", placeholder="Name of first variable")
    if(input$demo_data == "ampk-tor"){
      textInput("var1_id", h5("1st variable:"), value = "Domain similarity (forward)", width="100%", placeholder="Name of first variable")
    } else if(input$demo_data == "demo"){
      textInput("var1_id", h5("1st variable:"), value = "Domain similarity", width="100%", placeholder="Name of first variable")
    } else {
      filein <- input$mainInput
      if(is.null(filein)){return(textInput("var1_id", h5("1st variable:"), value = "Variable 1", width="100%", placeholder="Name of first variable"))}    # get var1/var2 names based on input (only if input file in long format table)
      
      inputType <- checkInputVadility(filein)
      
      if(inputType == "xml"){
        longDf <- xmlParser(filein$datapath)
        textInput("var1_id", h5("1st variable:"), value = colnames(longDf)[4], width="100%", placeholder="Name of first variable")
      } else if(inputType == "fasta"){
        longDf <- fastaParser(filein$datapath)
        textInput("var1_id", h5("1st variable:"), value = colnames(longDf)[4], width="100%", placeholder="Name of first variable")
      } else if(inputType == "long"){
        headerIn <- readLines(filein$datapath, n = 1)
        headerIn <- unlist(strsplit(headerIn,split = '\t'))
        textInput("var1_id", h5("1st variable:"), value = headerIn[4], width="100%", placeholder="Name of first variable")
      } else{
        textInput("var1_id", h5("1st variable:"), value = "Variable 1", width="100%", placeholder="Name of first variable")
      }
    }
  })
  
  output$var2_id.ui <- renderUI({
    # if(input$demo == TRUE){
    #   textInput("var2_id", h5("Second variable:"), value = "Traceability", width="100%", placeholder="Name of second variable")
    if(input$demo_data == "ampk-tor"){
      textInput("var2_id", h5("2nd variable:"), value = "Domain similarity (backward)", width="100%", placeholder="Name of second variable")
    } else if(input$demo_data == "demo"){
      textInput("var2_id", h5("2nd variable:"), value = "Traceability", width="100%", placeholder="Name of second variable")
    } else {
      filein <- input$mainInput
      if(is.null(filein)){return(textInput("var2_id", h5("2nd variable:"), value = "Variable 2", width="100%", placeholder="Name of second variable"))}    # get var1/var2 names based on input (only if input file in long format table)
      
      inputType <- checkInputVadility(filein)
      
      if(inputType == "xml"){
        longDf <- xmlParser(filein$datapath)
        textInput("var2_id", h5("2nd variable:"), value = colnames(longDf)[5], width="100%", placeholder="Name of second variable")
      } else if(inputType == "fasta"){
        longDf <- fastaParser(filein$datapath)
        textInput("var2_id", h5("2nd variable:"), value = colnames(longDf)[5], width="100%", placeholder="Name of second variable")
      } else if(inputType == "long"){
        headerIn <- readLines(filein$datapath, n = 1)
        headerIn <- unlist(strsplit(headerIn,split = '\t'))
        textInput("var2_id", h5("2nd variable:"), value = headerIn[5], width="100%", placeholder="Name of second variable")
      } else{
        textInput("var2_id", h5("2nd variable:"), value = "Variable 2", width="100%", placeholder="Name of second variable")
      }
    }
  })
  
  ######## variable 1 & 2 cutoff slidebar (main plot)
  output$var1_cutoff.ui <- renderUI({
    if(is.null(input$var1_id)){return()}
    if(input$var1_id == ""){
      sliderInput("var1",paste(input$var1_id,"cutoff:"), min = 1, max = 1, step = 0.025, value = c(1.0,1.0), width = 200)
    } else {
      sliderInput("var1",paste(input$var1_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(0.0,1.0), width = 200)
    }
  })
  
  output$var2_cutoff.ui <- renderUI({
    if(is.null(input$var2_id)){return()}
    if(input$var2_id == ""){
      sliderInput("var2",paste(input$var2_id,"cutoff:"), min = 1, max = 1, step = 0.025, value = c(1.0,1.0), width = 200)
    } else {
      sliderInput("var2",paste(input$var2_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(0.0,1.0), width = 200)
    }
  })
  
  ######## render filter slidebars for Customized plot
  output$var1_filter.ui <- renderUI({
    if(is.null(input$var1_id)){return()}
    if(input$var1_id == ""){
      sliderInput("var1cus",paste(input$var1_id,"cutoff:"), min = 1, max = 1, step = 0.025, value = c(1.0,1.0), width = 200)
    } else {
      sliderInput("var1cus",paste(input$var1_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var1[1],input$var1[2]), width = 200)
    }
  })
  
  output$var2_filter.ui <- renderUI({
    if(is.null(input$var2_id)){return()}
    if(input$var2_id == ""){
      sliderInput("var2cus",paste(input$var2_id,"cutoff:"), min = 1, max = 1, step = 0.025, value = c(1.0,1.0), width = 200)
    } else {
      sliderInput("var2cus",paste(input$var2_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var2[1],input$var2[2]), width = 200)
    }
  })
  
  output$percent_filter.ui <- renderUI({
    sliderInput("percent2",
                "% of present taxa:", min = 0, max = 1, step = 0.025, value = input$percent, width = 200)
  })
  
  ######## render filter slidebars for Distribution plot
  output$var1_dist.ui <- renderUI({
    if(is.null(input$var1_id)){return()}
    if(input$var1_id == ""){
      sliderInput("var1_dist",paste(input$var1_id,"cutoff:"), min = 1, max = 1, step = 0.025, value = c(1.0,1.0), width = 200)
    } else {
      sliderInput("var1_dist",paste(input$var1_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var1[1],input$var1[2]), width = 200)
    }
  })
  
  output$var2_dist.ui <- renderUI({
    if(is.null(input$var2_id)){return()}
    if(input$var2_id == ""){
      sliderInput("var2_dist",paste(input$var2_id,"cutoff:"), min = 1, max = 1, step = 0.025, value = c(1.0,1.0), width = 200)
    } else {
      sliderInput("var2_dist",paste(input$var2_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var2[1],input$var2[2]), width = 200)
    }
  })
  
  output$percent_dist.ui <- renderUI({
    sliderInput("percent_dist",
                "% of present taxa:", min = 0, max = 1, step = 0.025, value = input$percent, width = 200)
  })
  
  ######## render filter slidebars for Gene age estimation plot
  output$var1_age.ui <- renderUI({
    if(is.null(input$var1_id)){return()}
    if(input$var1_id == ""){
      sliderInput("var1_age",paste(input$var1_id,"cutoff:"), min = 1, max = 1, step = 0.025, value = c(1.0,1.0), width = 200)
    } else {
      sliderInput("var1_age",paste(input$var1_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var1[1],input$var1[2]), width = 200)
    }
  })
  
  output$var2_age.ui <- renderUI({
    if(is.null(input$var2_id)){return()}
    if(input$var2_id == ""){
      sliderInput("var2_age",paste(input$var2_id,"cutoff:"), min = 1, max = 1, step = 0.025, value = c(1.0,1.0), width = 200)
    } else {
      sliderInput("var2_age",paste(input$var2_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var2[1],input$var2[2]), width = 200)
    }
  })
  
  output$percent_age.ui <- renderUI({
    sliderInput("percent_age",
                "% of present taxa:", min = 0, max = 1, step = 0.025, value = input$percent, width = 200)
  })
  
  ######## render filter slidebars for Core gene finding function
  output$var1_cons.ui <- renderUI({
    if(is.null(input$var1_id)){return()}
    if(input$var1_id == ""){
      sliderInput("var1_cons",paste(input$var1_id,"cutoff:"), min = 1, max = 1, step = 0.025, value = c(1.0,1.0), width = 200)
    } else {
      sliderInput("var1_cons",paste(input$var1_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var1[1],input$var1[2]), width = 200)
    }
  })
  
  output$var2_cons.ui <- renderUI({
    if(is.null(input$var2_id)){return()}
    if(input$var2_id == ""){
      sliderInput("var2_cons",paste(input$var2_id,"cutoff:"), min = 1, max = 1, step = 0.025, value = c(1.0,1.0), width = 200)
    } else {
      sliderInput("var2_cons",paste(input$var2_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var2[1],input$var2[2]), width = 200)
    }
  })
  
  output$percent_cons.ui <- renderUI({
    sliderInput("percent_cons",
                "% of present taxa:", min = 0, max = 1, step = 0.025, value = 0.5, width = 200)
  })
  
  ######## update value for "main" filter slidebars based on "Customized", "Distribution", "Gene age estimation" slidebars
  observe({
    newVar1 <- input$var1cus
    
    if(is.null(input$var1_id)){return()}
    if(input$var1_id == ""){
      updateSliderInput(session, "var1", value = newVar1,min = 1, max = 1, step = 0.025)
    } else {
      updateSliderInput(session, "var1", value = newVar1,min = 0, max = 1, step = 0.025)
    }
  })
  observe({
    newVar2 <- input$var2cus
    
    if(is.null(input$var2_id)){return()}
    if(input$var2_id == ""){
      updateSliderInput(session, "var2", value = newVar2,min = 1, max = 1, step = 0.025)
    } else {
      updateSliderInput(session, "var2", value = newVar2,min = 0, max = 1, step = 0.025)
    }
  })
  observe({
    newPercent <- input$percent2
    updateSliderInput(session, "percent", value = newPercent,
                      min = 0, max = 1, step = 0.025)
  })
  
  observe({
    newVar1 <- input$var1_dist
    
    if(is.null(input$var1_id)){return()}
    if(input$var1_id == ""){
      updateSliderInput(session, "var1", value = newVar1,min = 1, max = 1, step = 0.025)
    } else {
      updateSliderInput(session, "var1", value = newVar1,min = 0, max = 1, step = 0.025)
    }
  })
  observe({
    newVar2 <- input$var2_dist
    
    if(is.null(input$var2_id)){return()}
    if(input$var2_id == ""){
      updateSliderInput(session, "var2", value = newVar2,min = 1, max = 1, step = 0.025)
    } else {
      updateSliderInput(session, "var2", value = newVar2,min = 0, max = 1, step = 0.025)
    }
  })
  observe({
    newPercent <- input$percent_dist
    updateSliderInput(session, "percent", value = newPercent,
                      min = 0, max = 1, step = 0.025)
  })
  
  observe({
    newVar1 <- input$var1_age
    
    if(is.null(input$var1_id)){return()}
    if(input$var1_id == ""){
      updateSliderInput(session, "var1", value = newVar1,min = 1, max = 1, step = 0.025)
    } else {
      updateSliderInput(session, "var1", value = newVar1,min = 0, max = 1, step = 0.025)
    }
  })
  observe({
    newVar2 <- input$var2_age
    
    if(is.null(input$var2_id)){return()}
    if(input$var2_id == ""){
      updateSliderInput(session, "var2", value = newVar2,min = 1, max = 1, step = 0.025)
    } else {
      updateSliderInput(session, "var2", value = newVar2,min = 0, max = 1, step = 0.025)
    }
  })
  observe({
    newPercent <- input$percent_age
    updateSliderInput(session, "percent", value = newPercent,
                      min = 0, max = 1, step = 0.025)
  })
  
  ####### render 2. variable relationship according to demo data ########
  output$var2_relation.ui <- renderUI({
    if(input$demo_data == "ampk-tor"){
      selectInput("var2_relation", label = h5("Relationship:"),
                  choices = list("Prot-Prot"="protein", "Prot-Spec"="species"),
                  selected = "protein",
                  width = 130)
    } else {
      selectInput("var2_relation", label = h5("Relationship:"),
                  choices = list("Prot-Prot"="protein", "Prot-Spec"="species"),
                  selected = "species",
                  width = 130)
    }
  })
  
  
  ########################################################
  
  ######## check oneseq fasta file exists
  output$oneSeq.existCheck <- renderUI({
    #f <- toString(input$oneseq.file)
    if(is.null(input$oneSeqFasta)){ return()}
    else{
      f <- input$oneSeqFasta$datapath
      if(!file.exists(f)){
        helpText("File not exists!!")
      } else {
        if(length(readLines(f, n=1)) == 0){
          helpText("is not a fasta file!!")
        } else {
          firstLine <- readLines(f, n=1)
          a <- substr(firstLine,1,1)
          if(a == ">"){
            HTML('<p><span style="color: #0000ff;"><strong>Please click CLOSE to comfirm!</strong></span></p>')
          } else {
            helpText("is not a fasta file!!")
          }
        }
      }
    }
  })
  
  ######## reset cutoffs of main plot
  observeEvent(input$reset_main, {
    shinyjs::reset("var1")
    shinyjs::reset("var2")
    shinyjs::reset("percent")
  })
  
  ######## reset config of main plot
  observeEvent(input$reset_main_config, {
    shinyjs::reset("x_size")
    shinyjs::reset("y_size")
    shinyjs::reset("legend_size")
    shinyjs::reset("x_angle")
    shinyjs::reset("dot_zoom")
  })
  
  ######## close main config
  observeEvent(input$apply_main_config, {
    toggleModal(session, "main_plot_config_bs", toggle = "close")
  })
  
  ######## reset cutoffs of Customized plot
  observeEvent(input$reset_selected, {
    shinyjs::reset("var1")
    shinyjs::reset("var2")
    shinyjs::reset("percent")
  })
  
  ######## reset config of customized plot
  observeEvent(input$reset_selected_config, {
    shinyjs::reset("x_size_select")
    shinyjs::reset("y_size_select")
    shinyjs::reset("legend_size_select")
    shinyjs::reset("x_angle_select")
    shinyjs::reset("dot_zoom_select")
  })
  
  ######## close customized config
  observeEvent(input$apply_selected_config, {
    toggleModal(session, "selected_plot_config_bs", toggle = "close")
  })
  
  ######## reset colors
  observeEvent(input$defaultColorVar2, {
    shinyjs::reset("low_color_var2")
    shinyjs::reset("high_color_var2")
  })
  
  observeEvent(input$defaultColorVar1, {
    shinyjs::reset("low_color_var1")
    shinyjs::reset("high_color_var1")
  })
  
  observeEvent(input$default_color_para, {
    shinyjs::reset("para_color")
  })
  
  
  ########################################################
  
  ######## check the validity of input newick tree ########
  checkNewick <- function(filein){
    tree <- read.table(file=filein$datapath, header=F,check.names=FALSE,comment.char="",fill=F)
    
    # get tree structure
    treeStruc <- gsub(regex("\\w"), "", as.character(tree$V1))
    
    open = str_count(treeStruc,"\\(")
    close = str_count(treeStruc,"\\)")
    comma = str_count(treeStruc,"\\,")
    singleton = str_count(treeStruc,"\\(\\)")
    
    # return check
    if(is.null(input$mainInput)){
      return(0) # don't check if main input is absent
    } else {
      if(singleton > 0){
        return(3) # tree contains singleton
      }
      
      if(open != close){
        return(1) # missing parenthesis
      } else {
        if(comma != (open+1)){
          return(2) # missing comma
        } else {
          # get list of tips
          nodeString <- gsub(regex("\\W+"), "#", as.character(tree$V1))
          nodeList <- unlist(strsplit(nodeString, "#"))
          # list of input taxa
          inputTaxa <- subset_taxa()
          
          missingTaxa <- list()
          j <- 1
          for(i in 1:length(nodeList)){
            if(nchar(nodeList[i]) > 0 & !(nodeList[i] %in% inputTaxa)){
              missingTaxa[[j]] <- nodeList[i]
              j <- j+1
            }
          }
          return(paste(missingTaxa, collapse="; ")) # contains taxa that not exist in main input
        }
      }
    }
    return(0)
  }
  
  ######## render checkNewick.ui
  output$checkNewick.ui <- renderUI({
    filein <- input$inputTree
    if(is.null(filein)){return()}
    
    checkNewick = checkNewick(filein)
    if(checkNewick == 1){
      updateButton(session, "do", disabled = TRUE)
      HTML("<p><em><span style=\"color: #ff0000;\"><strong>ERROR: Parenthesis(-es) missing!</strong></span></em></p>")
    } else if(checkNewick == 2){
      updateButton(session, "do", disabled = TRUE)
      HTML("<p><em><span style=\"color: #ff0000;\"><strong>ERROR: Comma(s) missing!</strong></span></em></p>")
    } else if(checkNewick == 3){
      updateButton(session, "do", disabled = TRUE)
      HTML("<p><em><span style=\"color: #ff0000;\"><strong>ERROR: Tree contains singleton!</strong></span></em></p>")
    } else if(checkNewick ==0){
      return()
    } else {
      updateButton(session, "do", disabled = TRUE)
      strong(em(paste0(checkNewick," not exist in main input file!")),style = "color:red")
    }
  })
  
  ######## check validity of main input file
  checkInputVadility <- function(filein){
    inputDt <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=F,check.names=FALSE,comment.char="",fill=T))
    # print(head(inputDt))
    if(is.na(inputDt[1,ncol(inputDt)])){
      return("moreCol")
    } else {
      names(inputDt) <- as.character(unlist(inputDt[1,]))
      
      if(grepl("<?xml",colnames(inputDt)[1])){
        return("xml")
      } else if(grepl(">",colnames(inputDt)[1]) == TRUE){
        return("fasta")
      } else {
        if(grepl("geneID",colnames(inputDt)[1])){
          if(is.na(pmatch("ncbi",colnames(inputDt)[3])) || is.na(pmatch("ncbi",colnames(inputDt)[4])) || is.na(pmatch("ncbi",colnames(inputDt)[5]))){
            return("long")
          } else {
            tmp <- inputDt[inputDt==""][1]
            if( !is.na(tmp) & tmp == ""){
              return("emptyCell")
            } else {
              return("wide")
            }
          }
        } else {
          return("noGeneID")
        }
      }
    }
  }
  
  ######## render input_check.ui
  output$input_check.ui <- renderUI({
    filein <- input$mainInput
    if(is.null(filein)){return()}
    
    inputType <- checkInputVadility(filein)
    if(inputType == "noGeneID"){
      updateButton(session, "do", disabled = TRUE)
      HTML("<font color=\"red\"><em><strong>ERROR: Unsupported input format. <a href=\"https://github.com/BIONF/PhyloProfile/wiki/Input-Data\" target=\"_blank\">Click here for more info</a></em></strong></font>")
    } else if(inputType == "emptyCell"){
      updateButton(session, "do", disabled = TRUE)
      em(strong("ERROR: Rows have unequal length",style = "color:red"))
    }
    else if(inputType == "moreCol"){
      updateButton(session, "do", disabled = TRUE)
      em(strong("ERROR: More columns than column names",style = "color:red"))
    } else {
      updateButton(session, "do", disabled = FALSE)
      return()
    }
  })
  
  
  ######## render download link for demo online files
  output$main_input_file.ui <- renderUI({
    # if(input$demo == TRUE){
    if(input$demo_data == "demo"){
      strong(a("Download demo input file", href="https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.long", target="_blank"))
    } else if(input$demo_data == "ampk-tor"){
      strong(a("Download demo input file", href="https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile", target="_blank"))
    } else {
      fileInput("mainInput",h5("Upload input file:"))
    }
  })
  
  output$domain_input_file.ui <- renderUI({
    # if(input$demo == TRUE){
    if(input$demo_data == "demo"){
      strong(a("Download demo domain files", href="https://github.com/BIONF/phyloprofile-data/tree/master/demo/domain_files", target="_blank"))
    } else if(input$demo_data == "ampk-tor"){
      strong(a("Download demo domain file", href="https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.domains_F", target="_blank"))
    } else {
      if(input$anno_choose == "from file"){
        fileInput("fileDomainInput","")
      } else {
        textInput("domainPath","","")
      }
    }
  })
  
  output$download_fastaDemo.ui <- renderUI({
    if(input$demo_data == "demo"){
      strong(a("Download demo fasta file", href="https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/fasta_file/concatenatedSeq.fa", target="_blank"))
    } else if(input$demo_data == "ampk-tor"){
      strong(a("Download demo fasta file", href="https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.extended.fa", target="_blank"))
    }
  })
  
  
  ####### render description for demo data #######
  output$demo_data_describe <- renderUI({
    if(input$demo_data == "none"){
      return()
    } else if(input$demo_data == "ampk-tor"){
      em(a("Data description", href="https://github.com/BIONF/phyloprofile-data/blob/master/expTestData/ampk-tor/README.md", target="_blank"))
    } else {
      em(a("Data description", href="https://github.com/BIONF/phyloprofile-data/blob/master/demo/README.md", target="_blank"))
    }
  })
  
  ######## get input taxa
  subset_taxa <- reactive({
    # if(input$demo == TRUE){
    if(input$demo_data == "demo"){
      data <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.wide", sep="\t", header=T, fill=T, stringsAsFactors = FALSE)
      inputTaxa <- colnames(data)
    } else if(input$demo_data == "ampk-tor"){
      data <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile", sep="\t", header=T, fill=T, check.names=FALSE,comment.char="")
      inputTaxa <- levels(data$ncbiID)
    } else {
      filein <- input$mainInput
      if(is.null(filein)){return()}
      else{
        
        if(length(unkTaxa()) == 0){
          # get list of input taxa (from main input file)
          inputType <- checkInputVadility(filein)
          if(inputType == "xml"){
            longDf <- xmlParser(filein$datapath)
            inputTaxa <- levels(longDf$ncbiID)
          } else if(inputType == "fasta"){
            longDf <- fastaParser(filein$datapath)
            inputTaxa <- levels(longDf$ncbiID)
          } else if(inputType == "long"){
            inputDf <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
            inputTaxa <- levels(inputDf$ncbiID)
          } else if(inputType == "wide"){
            inputTaxa <- readLines(filein$datapath, n = 1)
          } else {
            inputTaxa <- "NA"
          }
        } else {
          inputTaxa <- readLines(filein$datapath, n = 1)
        }
        
        inputTaxa <- unlist(strsplit(inputTaxa,split = '\t'))
        if(inputTaxa[1] == "geneID"){
          inputTaxa <- inputTaxa[-1]   # remove "geneID" element from vector inputTaxa
        }
        
        # return input taxa
        return(inputTaxa)
      }
    }
  })
  
  ######## check if there is any "unknown" taxon in input matrix
  unkTaxa <- reactive({
    # get list of input taxa (from main input file)
    # if(input$demo == TRUE){
    if(input$demo_data == "demo"){
      data <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.wide", sep="\t", header=T, fill=T, stringsAsFactors = FALSE)
      inputTaxa <- colnames(data)
    } else if(input$demo_data == "ampk-tor"){
      data <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile", sep="\t", header=T, fill=T, check.names=FALSE,comment.char="")
      inputTaxa <- levels(data$ncbiID)
    } else {
      filein <- input$mainInput
      if(is.null(filein)){return()}
      
      inputType <- checkInputVadility(filein)
      if(inputType == "xml"){
        longDf <- xmlParser(filein$datapath)
        inputTaxa <- levels(longDf$ncbiID)
      } else if(inputType == "fasta"){
        longDf <- fastaParser(filein$datapath)
        inputTaxa <- levels(longDf$ncbiID)
      } else if(inputType == "long"){
        inputDf <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
        inputTaxa <- levels(inputDf$ncbiID)
      } else if(inputType == "wide"){
        inputTaxa <- readLines(filein$datapath, n = 1)
      } else {
        inputTaxa <- c("NA")
      }
    }
    
    if(inputTaxa[1] == "NA"){
      return()
    } else {
      inputTaxa <- unlist(strsplit(inputTaxa,split = '\t'))
      if(inputTaxa[1] == "geneID"){
        inputTaxa <- inputTaxa[-1]   # remove "geneID" element from vector inputTaxa
      }
      
      if(!file.exists(isolate({"data/rankList.txt"}))){
        return(inputTaxa)
      } else {
        info = file.info("data/rankList.txt")
        if(info$size == 0){
          return(inputTaxa)
        } else {
          # get list of all available taxon (from /data/rankList.txt)
          pipeCmd <- paste0("cut -f 1 ",getwd(),"/data/rankList.txt")
          allTaxa <- unlist((read.table(pipe(pipeCmd))))
          
          # list of unknown taxa
          unkTaxa <- inputTaxa[!(inputTaxa %in% allTaxa)]
          # return list of unkTaxa
          return(unkTaxa)
        }
      }
    }
  })
  
  ### get status of unkTaxa for conditional panel
  output$unk_taxa_status <- reactive({
    unkTaxa <- unkTaxa()
    length(unkTaxa) > 0
  })
  outputOptions(output, "unk_taxa_status", suspendWhenHidden = FALSE)
  
  ### show full list of unkTaxa
  output$unk_taxa_full <- renderDataTable(option = list(searching = FALSE,pageLength = 10),{
    if(length(unkTaxa())>0){
      tb <- as.data.frame(unkTaxa())
      names(tb)[1] <- "New taxon"
      tb
    }
  })
  
  ######## check if data is loaded and "parse" button (get info from input) is clicked and confirmed
  v1 <- reactiveValues(parse = FALSE)
  observeEvent(input$but_parse, {
    toggleModal(session, "parse_confirm", toggle = "close")
    v1$parse <- input$but_parse
    updateButton(session, "but_parse", disabled = TRUE)
    toggleState("new_taxa_ask")
    toggleState("mainInput")
  })
  
  ######### create rankList.txt, idList.txt, taxonNamesReduced.txt from input file (if confirmed by BUTyes). 
  ######### and also create a full taxonomy matrix for sorting taxa (taxonomyMatrix.txt)
  invalidID <- reactiveValues(df=data.frame("Invalid NCBI ID(s)" = as.character(), stringsAsFactors = F))
  observe({
    filein <- input$mainInput
    if(is.null(filein)){return()}
    else{
      inputType <- checkInputVadility(filein)
      if(inputType == "xml" | inputType == "long" | inputType == "wide" | inputType == "fasta"){
        inputDf <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
        
        if(v1$parse == F){return()}
        else{
          ## get list of taxa need to be parsed (the one mising taxonomy information)
          if(v1$parse == T){
            unkTaxa <- unkTaxa()
            titleline <- c("geneID",unkTaxa)
          }
          
          invalidIDtmp <- list()
          # Create 0-row data frame which will be used to store data
          withProgress(message = 'Parsing input file', value = 0, {
            allTaxonInfo<-fread("data/taxonNamesFull.txt")
            newTaxaFromFile <- fread("data/newTaxa.txt", colClasses=c("ncbiID" = "character"))
            
            allTaxonInfo <- rbind(newTaxaFromFile,allTaxonInfo)
            
            rankList <- data.frame()
            idList <- data.frame()
            reducedInfoList <- data.frame()
            
            for(i in 2:length(titleline)){
              ## taxon ID  
              refID <- gsub("ncbi","",titleline[i])
              
              ## get info for this taxon
              refEntry <- allTaxonInfo[allTaxonInfo$ncbiID == refID,]
              
              if(nrow(refEntry) < 1){
                invalidIDtmp <- c(invalidIDtmp,refID)
              } else {
                if(nrow(reducedInfoList[reducedInfoList$X1 == refEntry$ncbiID,]) == 0){
                  refInfoList <- data.frame(matrix(c(refEntry$ncbiID,refEntry$fullName,refEntry$rank,refEntry$parentID), nrow=1, byrow=T),stringsAsFactors=FALSE)
                  reducedInfoList <- rbind(reducedInfoList,refInfoList)
                }
                
                ## parentID (used to check if hitting last rank, i.e. norank_1)
                lastID <- refEntry$parentID
                
                ## create list of rank for this taxon
                rank <- c(paste0("ncbi",refID),refEntry$fullName)
                if(refEntry$rank == "norank"){
                  rank <- c(rank,paste0("strain"))
                } else {
                  rank <- c(rank,refEntry$rank)
                }
                
                ## create list of IDs for this taxon
                ids <- list(paste0(refEntry$fullName,"#name"))
                if(refEntry$rank == "norank"){
                  ids <- c(ids,paste0(refEntry$ncbiID,"#","strain","_",refEntry$ncbiID))
                } else {
                  ids <- c(ids,paste0(refEntry$ncbiID,"#",refEntry$rank))
                }
                
                ## append info into rank and ids
                while(lastID != 1){
                  nextEntry <- allTaxonInfo[allTaxonInfo$ncbiID == lastID,]
                  
                  if(nrow(reducedInfoList[reducedInfoList$X1 == nextEntry$ncbiID,]) == 0){
                    nextEntryList <- data.frame(matrix(c(nextEntry$ncbiID,nextEntry$fullName,nextEntry$rank,nextEntry$parentID), nrow=1, byrow=T),stringsAsFactors=FALSE)
                    reducedInfoList <- rbind(reducedInfoList,nextEntryList)
                  }
                  
                  lastID <- nextEntry$parentID
                  
                  if("norank" %in% nextEntry$rank){
                    rank <- c(rank,paste0(nextEntry$rank,"_",nextEntry$ncbiID))
                    ids <- c(ids,paste0(nextEntry$ncbiID,"#",nextEntry$rank,"_",nextEntry$ncbiID))
                  } else {
                    rank <- c(rank,nextEntry$rank)
                    ids <- c(ids,paste0(nextEntry$ncbiID,"#",nextEntry$rank))
                  }
                }
                
                ## last rank and id
                rank <- c(rank,"norank_1")
                ids <- c(ids,"1#norank_1")
                
                ## append into rankList and idList files
                rankListTMP <- data.frame(matrix(unlist(rank), nrow=1, byrow=T),stringsAsFactors=FALSE)
                rankList <- rbind.fill(rankList,rankListTMP)
                idListTMP <- data.frame(matrix(unlist(ids), nrow=1, byrow=T),stringsAsFactors=FALSE)
                idList <- rbind.fill(idList,idListTMP)
              }
              
              # Increment the progress bar, and update the detail text.
              incProgress(1/(length(titleline)-1), detail = paste((i-1),"/",length(titleline)-1))
            }
          })
          ### save invalid IDs to invalidID$df
          invalidID$df <- as.data.frame(unlist(invalidIDtmp))
          
          if(nrow(invalidID$df) < 1){
            ### open existing files (idList.txt, rankList.txt and taxonNamesReduced.txt)
            ncol <- max(count.fields("data/rankList.txt", sep = '\t'))
            # print(ncol)
            oldIDList <- as.data.frame(read.table("data/idList.txt", sep='\t', header=F, check.names=FALSE, comment.char="", fill = T, stringsAsFactors=T, na.strings=c("","NA"), col.names=paste0('X', seq_len(ncol))))
            oldRankList <- as.data.frame(read.table("data/rankList.txt", sep='\t', header=F, check.names=FALSE, comment.char="", fill = T, stringsAsFactors=T, na.strings=c("","NA"), col.names=paste0('X', seq_len(ncol))))
            oldNameList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t', header=T, check.names=FALSE, comment.char="", fill = T, stringsAsFactors=T))
            
            ### and append new info into those files
            new_idList <- rbind.fill(oldIDList,idList)
            new_rankList <- rbind.fill(oldRankList,rankList)
            colnames(reducedInfoList) <- c("ncbiID","fullName","rank","parentID")
            new_nameList <- rbind.fill(oldNameList,reducedInfoList)
            
            write.table(new_idList[!duplicated(new_idList),], file ="data/idList.txt", col.names = F, row.names = F, quote = F, sep="\t")
            write.table(new_rankList[!duplicated(new_rankList),], file ="data/rankList.txt", col.names = F, row.names = F, quote = F, sep="\t")
            write.table(new_nameList[!duplicated(new_nameList),], file ="data/taxonNamesReduced.txt", col.names = T, row.names = F, quote = F, sep="\t")
            
            ### create taxonomy matrix
            taxMatrix <- taxonomyTableCreator("data/idList.txt","data/rankList.txt")
            write.table(taxMatrix,"data/taxonomyMatrix.txt",sep="\t",eol="\n",row.names=FALSE,quote = FALSE)
          }
        }
      }
    }
  })
  
  ######## output invalid NCBI ID
  output$invalidID.output <- renderTable({
    if(nrow(invalidID$df) < 1) {return()}
    else{
      outDf <- invalidID$df
      colnames(outDf) <- c("Invalid NCBI ID(s)")
      return(outDf)
    }
  })
  
  output$end_parsing_msg <- renderUI({
    if(nrow(invalidID$df) < 1) {
      strong(h4("PLEASE RELOAD THIS TOOL AFTER ADDING NEW TAXA!!!"),style = "color:red")
    }else{
      HTML('<p><strong><span style="color: #e12525;">SOME INVALID TAXON IDs HAVE BEEN FOUND!!</span><br>Please check the validity of the following IDs in 
           <a target="_blank" href="https://www.ncbi.nlm.nih.gov/taxonomy">NCBI taxonomy database</a>!</strong></p>')
    }
    })
  
  ######## list of taxonomy ranks for plotting
  output$rank_select = renderUI({
    # if(input$demo == TRUE){
    if(input$demo_data == "demo"){
      selectInput("rank_select", label = "",
                  choices = get_taxonomy_ranks(),
                  selected = "26_phylum")
    } else if(input$demo_data == "ampk-tor"){
      selectInput("rank_select", label = "",
                  choices = get_taxonomy_ranks(),
                  selected = "06_species")
    } else {
      selectInput("rank_select", label = "",
                  choices = get_taxonomy_ranks(),
                  selected = "06_species")
    }
  })
  
  ####### GET list of all (super)taxa
  alltaxa_list <- reactive({
    
    filein <- input$mainInput
    # if(is.null(filein) & input$demo == FALSE){return()}
    if(is.null(filein) & input$demo_data == "none"){return()}
    
    rank_select = input$rank_select
    
    if(rank_select == ""){return()}
    
    if(length(unkTaxa()) > 0){return()}
    
    ### load list of unsorted taxa
    Dt <- get_taxa_list(TRUE)
    
    
    ### load list of taxon name
    nameList <- get_name_list(TRUE, FALSE)
    
    rankName = substr(rank_select,4,nchar(rank_select))   # get rank name from rank_select
    #    rankNr = 0 + as.numeric(substr(rank_select,1,2))     # get rank number (number of column in unsorted taxa list - dataframe Dt)
    choice <- as.data.frame
    choice <- rbind(Dt[rankName])
    colnames(choice) <- "ncbiID"
    choice <- merge(choice,nameList,by="ncbiID",all = FALSE)
  })
  
  ### then output list of (super)taxa onto UI
  output$select = renderUI({
    choice <- alltaxa_list()
    choice$fullName <- as.factor(choice$fullName)
    
    # if(input$demo == TRUE){
    if(input$demo_data == "demo"){
      hellemDf <- data.frame("name" = c("Encephalitozoon hellem","Encephalitozoon hellem","Encephalitozoon","Unikaryonidae","Apansporoblastina","Apansporoblastina","Microsporidia","Fungi","Eukaryota"),
                             "rank" = c("strain","species","genus","family","order","class","phylum","kingdom","superkingdom"))
      rank_select = input$rank_select
      rankName = substr(rank_select,4,nchar(rank_select))
      
      selectInput('in_select',"",as.list(levels(choice$fullName)),hellemDf$name[hellemDf$rank == rankName])
    } else if(input$demo_data == "ampk-tor"){
      humanDf <- data.frame("name" = c("Homo sapiens","Homo sapiens","Homo","Hominidae","Primates","Mammalia","Chordata","Metazoa","Eukaryota"),
                            "rank" = c("strain","species","genus","family","order","class","phylum","kingdom","superkingdom"))
      rank_select = input$rank_select
      rankName = substr(rank_select,4,nchar(rank_select))
      
      selectInput('in_select',"",as.list(levels(choice$fullName)),humanDf$name[humanDf$rank == rankName])
    } else {
      selectInput('in_select',"",as.list(levels(choice$fullName)),levels(choice$fullName)[1])
    }
  })
  
  output$highlight_taxon_ui = renderUI({
    choice <- alltaxa_list()
    choice$fullName <- as.factor(choice$fullName)
    
    out <- as.list(levels(choice$fullName))
    out <- append("none",out)
    
    selectInput('taxonHighlight','Select (super)taxon to highlight:',out,selected=out[1])
  })
  
  #### update highlight_taxon_ui based on double clicked dot
  observe({
    choice <- alltaxa_list()
    choice$fullName <- as.factor(choice$fullName)
    
    out <- as.list(levels(choice$fullName))
    out <- append("none",out)
    
    if(!is.null(input$plot_dblclick)){
      if(input$x_axis == "genes"){
        corX = round(input$plot_dblclick$y);
        corY = round(input$plot_dblclick$x)
      } else {
        corX = round(input$plot_dblclick$x);
        corY = round(input$plot_dblclick$y)
      }
      
      dataHeat <- dataHeat()
      supertaxa <- levels(dataHeat$supertaxon)
      spec <- toString(supertaxa[corX])
      selectedIndex = match(substr(spec,6,nchar(spec)),out)
      
      updateSelectInput(session,'taxonHighlight',label = 'Select (super)taxon to highlight:',choices=out,selected=out[selectedIndex])
    }
  })
  
  ######## Get list of genes for highlighting
  output$highlight_gene_ui = renderUI({
    geneList <- dataHeat()
    geneList$geneID <- as.factor(geneList$geneID)
    
    out <- as.list(levels(geneList$geneID))
    out <- append("none",out)
    
    selectInput('geneHighlight','Highlight:',out,selected=out[1])
  })
  
  #### update highlight_gene_ui based on double clicked dot
  observe({
    if(!is.null(input$plot_dblclick)){
      geneList <- dataHeat()
      if(input$apply_cluster == TRUE){
        geneList <- clusteredDataHeat()
      }
      
      geneList$geneID <- as.factor(geneList$geneID)
      
      out <- as.list(levels(geneList$geneID))
      out <- append("none",out)
      
      clickedInfo <- mainpoint_info()
      
      if(input$x_axis == "genes"){
        corX = round(input$plot_dblclick$y);
        corY = round(input$plot_dblclick$x)
      } else {
        corX = round(input$plot_dblclick$x);
        corY = round(input$plot_dblclick$y)
      }
      updateSelectInput(session,'geneHighlight',label = 'Highlight:',choices=out,selected=out[corY+1])
    } else {return()}
  })
  
  ######## print total number of genes
  output$total_gene_number.ui <- renderUI({
    geneList <- preData()
    geneList$geneID <- as.factor(geneList$geneID)
    out <- as.list(levels(geneList$geneID))
    if(length(out) > 0){
      # em(paste0("Total number of genes:  ",length(out)))
      strong(paste0("Total number of genes:  ",length(out)))
    }
  })
  
  ######## enable "PLOT" button
  observeEvent(input$rank_select,({
    if(input$rank_select == ""){
      updateButton(session, "do", disabled = TRUE)
    } else{
      unkTaxa <- unkTaxa()
      if(length(unkTaxa) == 0){
        updateButton(session, "do", disabled = FALSE)
      }
    }
  }))
  
  ######## move to main tab when "PLOT" button has been clicked
  observe({
    # use tabsetPanel 'id' argument to change tabs
    if (input$do > 0) {
      updateTabsetPanel(session, "tabs", selected = "Main profile")
    }
  })
  
  ######## disable main input, genelist input and initial questions
  observe({
    # use tabsetPanel 'id' argument to change tabs
    if (input$do > 0) {
      toggleState("mainInput")
      toggleState("gene_list_selected")
      toggleState("demo_data")
    }
  })
  
  ######## disable demo checkbox and update var2_aggregate_by to mean if using demo data
  observe({
    # if (input$demo == TRUE) {
    if (input$demo_data == "demo") {
      ### disable demo checkbox
      # toggleState("demo")
      ### update var2_aggregate_by to mean
      updateSelectInput(session,"var2_aggregate_by",
                        choices = list("Max"="max", "Min"="min","Mean"="mean","Median"="median"),
                        selected = "mean")
    }
  })
  
  #############################################################
  ######################  ADD NEW TAXA  #######################
  #############################################################
  newTaxa <- reactiveValues()
  newTaxa$Df <- data.frame("ncbiID"=numeric(),"fullName"=character(),"rank"=character(),"parentID"=numeric(),stringsAsFactors=FALSE)
  newIndex <- reactiveValues()
  newIndex$value <- 1
  
  observeEvent(input$new_add, {
    newTaxa$Df[newIndex$value,] <- c(input$new_id,input$new_name,input$new_rank,input$new_parent)
    newIndex$value <- newIndex$value + 1
    updateTextInput(session, "new_id", value=as.numeric(input$new_id)+1)
    updateTextInput(session, "new_name", value="")
    updateTextInput(session, "new_rank", value="norank")
    updateTextInput(session, "new_parent", value="")
  })
  
  observeEvent(input$new_done, {
    toggleModal(session, "add_taxa_windows", toggle = "close")
    write.table(newTaxa$Df,"data/newTaxa.txt",sep="\t",eol="\n",row.names=FALSE,quote = FALSE)
  })
  
  #############################################################
  ##################  PROCESSING INPUT DATA ###################
  #############################################################
  
  ######## check if data is loaded and "plot" button is clicked
  v <- reactiveValues(doPlot = FALSE)
  observeEvent(input$do, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    v$doPlot <- input$do
    filein <- input$mainInput
    # if(is.null(filein) & input$demo==FALSE){
    if(is.null(filein) & input$demo_data == "none"){
      v$doPlot <- FALSE
      updateButton(session, "do", disabled = TRUE)
    }
  })
  
  ######## create rooted tree from a matrix
  createRootedTree <- function(df,rootTaxon){
    ### calculate distance matrix
    taxdis <- tryCatch(taxa2dist(df), error = function(e) e)
    
    ### root tree 
    tree <- as.phylo(hclust(taxdis))
    tree <- root(tree, outgroup = rootTaxon, resolve.root = TRUE)
    
    ### return
    return(tree)
  }
  
  ######## sort supertaxa list based on chosen reference taxon
  sortTaxaFromTree <- function(tree){
    ### and get ordered taxon list
    is_tip <- tree$edge[,2] <= length(tree$tip.label)
    ordered_tips <- tree$edge[is_tip, 2]
    taxonList <- rev(tree$tip.label[ordered_tips])
    
    ### return
    return(taxonList)
  }
  
  sortedtaxa_list <- reactive({
    if(v$doPlot == FALSE){return()}
    
    ####### Get representative taxon ####### 
    
    ### load list of unsorted taxa
    Dt <- get_taxa_list(FALSE)
    
    ### load list of taxon name
    nameList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T,fill = TRUE))
    nameList$fullName <- as.character(nameList$fullName)
    
    ### input parameters
    rank_select = input$rank_select
    rankName = substr(rank_select,4,nchar(rank_select))   # get rank name from rank_select
    rankNr = 0 + as.numeric(substr(rank_select,1,2))     # get rank number (number of column in unsorted taxa list - dataframe Dt)
    
    ### get selected supertaxon ID
    taxa_list <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
    allTaxa <- alltaxa_list()
    rankNameTMP <- allTaxa$rank[allTaxa$fullName == input$in_select]
    if(rankName == "strain"){
      superID <- as.numeric(taxa_list$ncbiID[taxa_list$fullName == input$in_select & taxa_list$rank == "norank"])
    } else {
      superID <- as.numeric(taxa_list$ncbiID[taxa_list$fullName == input$in_select & taxa_list$rank == rankNameTMP[1]])
    }
    
    ### representative taxon
    repTaxon <- Dt[Dt[,rankName]==superID,][1,]
    
    ####### Sort taxon list based on hclust tree #######
    
    ### prepare Df for calculating distance matrix
    distDf <- subset(Dt, select = -c(ncbiID,fullName))
    row.names(distDf) <- distDf$abbrName
    distDf <- distDf[,-1]
    
    ### get sorted taxon IDs & sort full taxonomy info dataframe
    treeIn <- input$inputTree
    if(is.null(treeIn)){
      taxaTree <- createRootedTree(distDf,as.character(repTaxon$abbrName))
    } else {
      taxaTree <- read.tree(file=treeIn$datapath)
    }
    taxonList <- sortTaxaFromTree(taxaTree)
    sortedDt <- Dt[match(taxonList, Dt$abbrName),]
    
    ### subset to get list of input taxa only
    inputTaxa <- subset_taxa()
    sortedDt <- subset(sortedDt, abbrName %in% inputTaxa)
    
    ### get only taxonIDs list of selected rank and rename columns
    sortedOut <- subset(sortedDt,select=c("abbrName","ncbiID","fullName",as.character(rankName)))
    colnames(sortedOut) <- c("abbrName","species","fullName","ncbiID")
    
    ### add name of supertaxa into sortedOut list
    sortedOut <- merge(sortedOut,nameList,by="ncbiID",all.x = TRUE,sort = FALSE)
    sortedOut$species <- as.character(sortedOut$species)
    
    ### add order_prefix to supertaxon name
    ### and add prefix "ncbi" to taxon_ncbiID (column "species")
    prefix = 1001
    
    sortedOut$sortedSupertaxon <- 0   ## create new column for sorted supertaxon
    sortedOut$sortedSupertaxon[1] <- paste0(prefix,"_",sortedOut$fullName.y[1])
    sortedOut$species[1] <- paste0("ncbi",sortedOut$species[1])
    
    for(i in 2:nrow(sortedOut)){
      if(sortedOut$fullName.y[i] != sortedOut$fullName.y[i-1]){     ## increase prefix if changing to another supertaxon
        prefix = prefix + 1
      }
      sortedOut$sortedSupertaxon[i] <- paste0(prefix,"_",sortedOut$fullName.y[i])
      sortedOut$species[i] <- paste0("ncbi",sortedOut$species[i])
    }
    
    ### final sorted supertaxa list
    sortedOut$taxonID <- 0
    sortedOut$category <- "cat"
    sortedOut <- sortedOut[,c("abbrName","taxonID","fullName.x","species","ncbiID","sortedSupertaxon","rank","category")]
    colnames(sortedOut) <- c("abbrName","taxonID","fullName","ncbiID","supertaxonID","supertaxon","rank","category")
    
    sortedOut$taxonID <- as.numeric(sortedOut$taxonID)
    sortedOut$ncbiID <- as.factor(sortedOut$ncbiID)
    sortedOut$supertaxon <- as.factor(sortedOut$supertaxon)
    sortedOut$category <- as.factor(sortedOut$category)
    
    ### return data frame
    return(sortedOut)
  })
  
  ############### PARSING DATA FROM INPUT MATRIX:
  ############### get (super)taxa names (3)
  ############### calculate percentage of presence (4), max/min/mean/median VAR1 (5) and VAR2 (6) if group input taxa list into higher taxonomy rank
  
  ### check if "no ordering gene IDs" has been checked
  output$apply_cluster_check.ui <- renderUI({
    if(input$ordering == FALSE){
      HTML('<p><em>(Check "Ordering sequence IDs" check box in <strong>Input & settings tab</strong>&nbsp;to enable this function)</em></p>')
    }
  })
  
  observe({
    if(input$ordering == FALSE){
      shinyjs::disable('apply_cluster')
    } else {
      shinyjs::enable('apply_cluster')
    }
  })
  
  ### unsorting function to keep user defined geneID order
  unsortID <- function(data,order){
    data$geneID <- as.factor(data$geneID)
    if(order == FALSE){
      data$geneID <- factor(data$geneID, levels = unique(data$geneID))  ######### keep user defined geneID order
    }
    return(data)
  }
  
  ### subset data
  preData <- reactive({
    # get list of gene of interest (from a separated file)
    listGene <- list()
    end_index <- input$end_index
    if(is.na(input$end_index)){end_index <- 30}
    
    if(input$gene_list_selected == 'from file'){
      listIn <- input$list
      if(!is.null(listIn)){
        list <- as.data.frame(read.table(file=listIn$datapath, header=FALSE))
        listGeneOri <- list$V1
        if(input$st_index <= length(listGeneOri)){
          listGene <- listGeneOri[listGeneOri[input$st_index:end_index]]
        } else {
          listGene <- listGeneOri
        }
      }
    }
    
    # if(input$demo == TRUE){
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      if(input$demo_data == "demo"){
        inputDf <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.long", sep="\t", header=T, fill=T, stringsAsFactors = FALSE)
      } else {
        inputDf <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile", sep="\t", header=T, fill=T, stringsAsFactors = FALSE)
      }
      inputDf <- unsortID(inputDf,input$ordering)
      
      if(length(listGene) >= 1){
        data <- inputDf[inputDf$geneID %in% listGene,]
      } else {
        subsetID <- levels(as.factor(inputDf$geneID))[input$st_index:end_index]
        data <- inputDf[inputDf$geneID %in% subsetID,]
      }
      
      if(ncol(data) < 5){
        for(i in 1:(5-ncol(data))){
          data[paste0("newVar",i)] <- 1
        }
      }
      colnames(data) <- c("geneID","ncbiID","orthoID","var1","var2")
    } else {
      filein <- input$mainInput
      if(is.null(filein)){return()}
      
      inputType <- checkInputVadility(filein)
      if(inputType == "xml"){
        longDf <- xmlParser(filein$datapath)
        longDf <- unsortID(longDf,input$ordering)
        
        if(length(listGene) >= 1){
          data <- longDf[longDf$geneID %in% listGene,]
        } else {
          subsetID <- levels(longDf$geneID)[input$st_index:end_index]
          data <- longDf[longDf$geneID %in% subsetID,]
        }
        
        if(ncol(data) < 5){
          for(i in 1:(5-ncol(data))){
            data[paste0("newVar",i)] <- 1
          }
        }
        colnames(data) <- c("geneID","ncbiID","orthoID","var1","var2")
      } else if(inputType == "fasta"){
        longDf <- fastaParser(filein$datapath)
        longDf <- unsortID(longDf,input$ordering)
        
        if(length(listGene) >= 1){
          data <- longDf[longDf$geneID %in% listGene,]
        } else {
          subsetID <- levels(longDf$geneID)[input$st_index:end_index]
          data <- longDf[longDf$geneID %in% subsetID,]
        }
        
        if(ncol(data) < 5){
          for(i in 1:(5-ncol(data))){
            data[paste0("newVar",i)] <- 1
          }
        }
        colnames(data) <- c("geneID","ncbiID","orthoID","var1","var2")
      } else if(inputType == "long"){
        inputDf <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
        inputDf <- unsortID(inputDf,input$ordering)
        
        if(length(listGene) >= 1){
          data <- inputDf[inputDf$geneID %in% listGene,]
        } else {
          subsetID <- levels(inputDf$geneID)[input$st_index:end_index]
          data <- inputDf[inputDf$geneID %in% subsetID,]
        }
        
        if(ncol(data) < 5){
          for(i in 1:(5-ncol(data))){
            data[paste0("newVar",i)] <- 1
          }
        }
        colnames(data) <- c("geneID","ncbiID","orthoID","var1","var2")
      } else if (inputType == "wide"){
        inputDf <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
        mdData <- melt(inputDf,id="geneID")
        mdData <- unsortID(mdData,input$ordering)
        
        if(length(listGene) >= 1){
          mdData <- mdData[mdData$geneID %in% listGene,]
        } else {
          subsetID <- levels(mdData$geneID)[input$st_index:end_index]
          mdData <- mdData[mdData$geneID %in% subsetID,]
        }
        
        splitCol <- data.frame(do.call('rbind', strsplit(as.character(mdData$value), '#', fixed=TRUE)))
        data <- cbind(mdData[,c('geneID','variable')],splitCol)
        
        if(ncol(data) < 5){
          for(i in 1:(5-ncol(data))){
            data[paste0("newVar",i)] <- 0
          }
        }
        colnames(data) <- c("geneID","ncbiID","orthoID","var1","var2")
      } else {
        data <- data.frame("geneID" = character(),"ncbiID" = character(),"orthoID" = character(),"var1" = character(),"var2" = character(), stringsAsFactors = F)
      }
    }
    
    ### return preData
    return(data)
  })
  
  ### get all information for input data
  get_data_filtered <- reactive({
    mdData <- preData()
    ### count number of inparalogs
    paralogCount <- plyr::count(mdData,c('geneID','ncbiID'))
    mdData <- merge(mdData,paralogCount,by=c('geneID','ncbiID'))
    colnames(mdData)[ncol(mdData)] <- "paralog"
    
    ### (3) GET SORTED TAXONOMY LIST (3) ###
    taxa_list <- sortedtaxa_list()
    
    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxa_list,'supertaxon')
    
    # merge mdData, mdDataVar2 and taxa_list to get taxonomy info
    taxaMdData <- merge(mdData,taxa_list,by='ncbiID')
    taxaMdData$var1 <- suppressWarnings(as.numeric(as.character(taxaMdData$var1)))
    taxaMdData$var2 <- suppressWarnings(as.numeric(as.character(taxaMdData$var2)))
    
    # ### (4) calculate PERCENTAGE of PRESENT SPECIES (4) ###
    finalPresSpecDt <- calcPresSpec(taxaMdData, taxaCount)
    
    ### (5) calculate max/min/mean/median VAR1 for every supertaxon of each gene (5) ###
    # remove NA rows from taxaMdData
    taxaMdDataNoNA <- taxaMdData[!is.na(taxaMdData$var1),]
    # calculate m var1
    mVar1Dt <- aggregate(taxaMdDataNoNA[,"var1"],list(taxaMdDataNoNA$supertaxon,taxaMdDataNoNA$geneID),FUN=input$var1_aggregate_by)
    colnames(mVar1Dt) <- c("supertaxon","geneID","mVar1")
    
    ### (6) calculate max/min/mean/median VAR2 for each super taxon (6) ###
    # remove NA rows from taxaMdData
    taxaMdDataNoNA_var2 <- taxaMdData[!is.na(taxaMdData$var2),]
    # calculate max/min/mean/median VAR2
    if(nrow(taxaMdDataNoNA_var2) > 0){
      mVar2Dt <- aggregate(taxaMdDataNoNA_var2[,"var2"],list(taxaMdDataNoNA_var2$supertaxon,taxaMdDataNoNA_var2$geneID),FUN=input$var2_aggregate_by)
      colnames(mVar2Dt) <- c("supertaxon","geneID","mVar2")
    } else {
      mVar2Dt <- taxaMdData[,c("supertaxon","geneID")]
      mVar2Dt$mVar2 <- 0
    }
    
    ### (5+6) & join mVar2 together with mVar1 scores into one df (5+6)
    scoreDf <- merge(mVar1Dt,mVar2Dt, by=c("supertaxon","geneID"), all = TRUE)
    
    ### (4+5+6) add presSpec and mVar1 into taxaMdData (4+5+6)
    presMdData <- merge(taxaMdData,finalPresSpecDt,by=c('geneID','supertaxon'),all.x = TRUE)
    fullMdData <- merge(presMdData,scoreDf,by=c('geneID','supertaxon'), all.x = TRUE)
    fullMdData <- merge(fullMdData,taxaCount,by=('supertaxon'), all.x = TRUE)
    # rename "freq" into "numberSpec"
    names(fullMdData)[names(fullMdData)=="freq"] <- "numberSpec"
    
    fullMdData$fullName <- as.vector(fullMdData$fullName)
    names(fullMdData)[names(fullMdData)=="orthoID.x"] <- "orthoID"
    fullMdData <- fullMdData[!duplicated(fullMdData), ] ### parsed input data frame !!!
    
    return(fullMdData)
  })
  
  #############################################################
  ############## DATA & PLOT FOR MAIN PROFILE #################
  #############################################################
  
  ######## REDUCE DATA FROM SPECIES LEVEL TO SUPERTAXA LEVEL
  ######## this data set contain only supertaxa and their value (%present, mVar1 & mVar2) for each gene
  dataSupertaxa <- reactive({
    fullMdData <- get_data_filtered()
    
    flag <- 1 # to check if working with the lowest taxonomy rank; 1 for NO; 0 for YES
    if(length(unique(levels(as.factor(fullMdData$numberSpec)))) == 1){
      if(unique(levels(as.factor(fullMdData$numberSpec))) == 1){
        superDfExt <- fullMdData[,c("geneID","supertaxon","supertaxonID","var1","presSpec","category","orthoID","var2","paralog")]
        flag <- 0
      }
    }
    
    if(flag == 1){
      ### get representative orthoID that has m VAR1 for each supertaxon
      mOrthoID <- fullMdData[,c('geneID','supertaxon','var1','mVar1','orthoID')]
      mOrthoID <- subset(mOrthoID,mOrthoID$var1 == mOrthoID$mVar1)
      colnames(mOrthoID) <- c('geneID','supertaxon','var1','mVar1','orthoID')
      mOrthoID <- mOrthoID[!is.na(mOrthoID$orthoID),]
      mOrthoID <- mOrthoID[,c('geneID','supertaxon','orthoID')]
      mOrthoID <- mOrthoID[!duplicated(mOrthoID[,1:2]), ]
      
      ### get data set for phyloprofile plotting (contains only supertaxa info)
      superDf <- subset(fullMdData,select=c('geneID','supertaxon','supertaxonID','mVar1','presSpec','category','mVar2','paralog'))
      superDf$paralog <- 1
      superDf <- superDf[!duplicated(superDf), ]
      
      superDfExt <- merge(superDf,mOrthoID, by=c('geneID','supertaxon'),all.x=TRUE)
      superDfExt <- superDfExt[,c("geneID","supertaxon","supertaxonID","mVar1","presSpec","category","orthoID","mVar2",'paralog')]
      
      ### rename mVar to var
      names(superDfExt)[names(superDfExt)=="mVar1"] <- "var1"
      names(superDfExt)[names(superDfExt)=="mVar2"] <- "var2"
    }
    
    # print(superDfExt[superDfExt$geneID == "ampk_ACACB" & superDfExt$supertaxon == "1001_Chordata",])
    # print("END2222")
    
    return(superDfExt)
  })
  
  
  ######## get list of all sequence IDs for selectize input (customized profile)
  output$gene_in = renderUI({
    filein <- input$mainInput
    fileCustom <- input$custom_file
    
    # if(input$demo == TRUE){
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      filein = 1
    }
    
    if(is.null(filein) & is.null(fileCustom)){return(selectInput('inSeq','',"all"))}
    if(v$doPlot == FALSE){return(selectInput('inSeq','',"all"))}
    else{
      ### full list
      data <- as.data.frame(get_data_filtered())
      data$geneID <- as.character(data$geneID)
      data$geneID <- as.factor(data$geneID)
      outAll <- as.list(levels(data$geneID))
      outAll <- append("all",outAll)
      #selectInput('inSeq','',out,selected=out[1],multiple=TRUE)
      
      if (input$add_custom_profile == TRUE){
        out <- selectedgene_age()
        if(length(out)>0){
          selectInput('inSeq','',out,selected=as.list(out),multiple=TRUE,selectize=FALSE)
        }
        else {
          selectInput('inSeq','',outAll,selected=outAll[1],multiple=TRUE,selectize=FALSE)
        }
      } else if(input$add_cluster_cutom_profile == TRUE){
        out <- brushed_clusterGene()
        if(length(out)>0){
          selectInput('inSeq','',out,selected=as.list(out),multiple=TRUE,selectize=FALSE)
        }
        else {
          selectInput('inSeq','',outAll,selected=outAll[1],multiple=TRUE,selectize=FALSE)
        }
      }
      else if(input$add_cons_gene_custom_profile == TRUE){
        out <- cons_geneDf()
        if(length(out)>0){
          selectInput('inSeq','',out,selected=as.list(out),multiple=TRUE,selectize=FALSE)
        }
        else {
          selectInput('inSeq','',outAll,selected=outAll[1],multiple=TRUE,selectize=FALSE)
        }
      }
      # else if(input$addGcGenesCustomProfile == TRUE){
      #   out <- significantGenesGroupCompairison$geneID
      #   if(length(out)>0){
      #     selectInput('inSeq','',out,selected=as.list(out),multiple=TRUE,selectize=FALSE)
      #   }
      #   else {
      #     selectInput('inSeq','',outAll,selected=outAll[1],multiple=TRUE,selectize=FALSE)
      #   }
      # 
      # }
      else {
        if(is.null(fileCustom)){
          selectInput('inSeq','',outAll,selected=outAll[1],multiple=TRUE,selectize=FALSE)
        }
        else {
          customList <- as.data.frame(read.table(file=fileCustom$datapath, header=FALSE))
          customList$V1 <- as.factor(customList$V1)
          out <- as.list(levels(customList$V1))
          selectInput('inSeq','',out,selected=out,multiple=TRUE,selectize=FALSE)
        }
      }
    }
  })
  
  ######## get list of all taxa for selectize input
  output$taxa_in = renderUI({
    filein <- input$mainInput
    # if(input$demo == TRUE){
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      filein = 1
    }
    
    if(is.null(filein)){return(selectInput('inTaxa','',"all"))}
    if(v$doPlot == FALSE){return(selectInput('inTaxa','',"all"))}
    else{
      choice <- alltaxa_list()
      choice$fullName <- as.factor(choice$fullName)
      
      out <- as.list(levels(choice$fullName))
      out <- append("all",out)
      
      if(input$apply_cus_taxa == TRUE){
        out <- cus_taxaName()
        selectInput('inTaxa','',out,selected=out,multiple=TRUE,selectize=FALSE)
      } else {
        selectInput('inTaxa','',out,selected=out[1],multiple=TRUE,selectize=FALSE)
      }
    }
  })
  
  ######## heatmap data input
  dataHeat <- reactive({
    
    ### get all cutoffs
    percent_cutoff_min <- input$percent[1]
    percent_cutoff_max <- input$percent[2]
    var1_cutoff_min <- input$var1[1]
    var1_cutoff_max <- input$var1[2]
    var2_cutoff_min <- input$var2[1]
    var2_cutoff_max <- input$var2[2]
    
    ### check input file
    filein <- input$mainInput
    # if(input$demo == TRUE){
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      filein = 1
    }
    if(is.null(filein)){return()}
    
    dataHeat <- dataSupertaxa()
    
    # get selected supertaxon name
    split <- strsplit(as.character(input$in_select),"_")
    in_select <- as.character(split[[1]][1])
    
    ### replace insufficient values according to the thresholds by NA or 0
    dataHeat$presSpec[dataHeat$supertaxon != in_select & dataHeat$presSpec < percent_cutoff_min] <- 0
    dataHeat$presSpec[dataHeat$supertaxon != in_select & dataHeat$presSpec > percent_cutoff_max] <- 0
    
    dataHeat$presSpec[dataHeat$supertaxon != in_select & dataHeat$var1 < var1_cutoff_min] <- 0
    dataHeat$presSpec[dataHeat$supertaxon != in_select & dataHeat$var1 > var1_cutoff_max] <- 0
    
    if(input$var1_relation == "protein"){
      if(input$var2_relation == "protein"){
        ### prot-prot: remove complete cell if one variable not sufficient
        dataHeat$presSpec[dataHeat$supertaxon != in_select & dataHeat$var2 < var2_cutoff_min] <- 0
        dataHeat$presSpec[dataHeat$supertaxon != in_select & dataHeat$var2 > var2_cutoff_max] <- 0
        
        dataHeat$var2[dataHeat$supertaxon != in_select & dataHeat$var1 < var1_cutoff_min] <- NA
        dataHeat$var2[dataHeat$supertaxon != in_select & dataHeat$var1 > var1_cutoff_max] <- NA
        dataHeat$var1[dataHeat$supertaxon != in_select & dataHeat$var2 < var2_cutoff_min] <- NA
        dataHeat$var1[dataHeat$supertaxon != in_select & dataHeat$var2 > var2_cutoff_max] <- NA
      } else {
        ### prot-spec: var1 depend on var2
        dataHeat$presSpec[dataHeat$supertaxon != in_select & dataHeat$var2 < var2_cutoff_min] <- 0
        dataHeat$presSpec[dataHeat$supertaxon != in_select & dataHeat$var2 > var2_cutoff_max] <- 0
      }
    } else {
      if(input$var2_relation == "species"){
        ### spec-spec: remove var1 and var2 independently
        # dataHeat$presSpec[dataHeat$supertaxon != in_select & dataHeat$var1 < var1_cutoff_min] <- 0
        # dataHeat$presSpec[dataHeat$supertaxon != in_select & dataHeat$var1 > var1_cutoff_max] <- 0
      } else {
        ### spec-prot: var2 depend on var1
        dataHeat$var2[dataHeat$supertaxon != in_select & dataHeat$var1 < var1_cutoff_min] <- NA
        dataHeat$var2[dataHeat$supertaxon != in_select & dataHeat$var1 > var1_cutoff_max] <- NA
      }
    }
    
    dataHeat$var1[dataHeat$supertaxon != in_select & dataHeat$var1 < var1_cutoff_min] <- NA
    dataHeat$var1[dataHeat$supertaxon != in_select & dataHeat$var1 > var1_cutoff_max] <- NA
    dataHeat$var2[dataHeat$supertaxon != in_select & dataHeat$var2 < var2_cutoff_min] <- NA
    dataHeat$var2[dataHeat$supertaxon != in_select & dataHeat$var2 > var2_cutoff_max] <- NA
    
    dataHeat <- droplevels(dataHeat)  ### delete unused levels
    dataHeat$geneID <- as.factor(dataHeat$geneID)
    
    return(dataHeat)
  })
  
  ########### cluster heatmap data
  clusteredDataHeat <- reactive({
    dataHeat <- dataHeat()
    if(nrow(dataHeat) < 1){return()}
    
    # dataframe for calculate distance matrix
    subDataHeat <- subset(dataHeat,dataHeat$presSpec > 0)
    subDataHeat <- subDataHeat[,c('geneID','supertaxon','presSpec')]
    subDataHeat <- subDataHeat[!duplicated(subDataHeat),]
    
    wideData <- spread(subDataHeat, supertaxon, presSpec)
    dat <- wideData[,2:ncol(wideData)]  # numerical columns
    rownames(dat) <- wideData[,1]
    dat[is.na(dat)] <- 0
    
    # get clustered gene ids
    clusteredGeneIDs <- clusteredGeneList(dat,input$dist_method,input$cluster_method)
    
    # sort original data according to clusteredGeneIDs
    dataHeat$geneID <- factor(dataHeat$geneID, levels = clusteredGeneIDs)
    return(dataHeat)
  })
  
  ########### render dot size to dot_size_info
  output$dot_size_info <- renderUI({
    if (v$doPlot == FALSE) return()
    
    dataHeat <- dataHeat()
    dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
    presentVl <- dataHeat$presSpec[!is.na(dataHeat$presSpec)]
    
    minDot <- (floor(min(presentVl)*10)/10*5)*(1+input$dot_zoom)
    maxDot <- (floor(max(presentVl)*10)/10*5)*(1+input$dot_zoom)
    
    em(paste0("current point's size: ",minDot," - ",maxDot))
  })
  
  ########### create profile heatmap
  mainPlot <- function(){
    if (v$doPlot == FALSE) return()
    dataHeat <- dataHeat()
    
    ### cluster dataHeat (if selected)
    if(input$apply_cluster == TRUE){
      dataHeat <- clusteredDataHeat()
    }
    
    ### reduce number of inparalogs based on filtered dataHeat
    dataHeatTB <- data.table(na.omit(dataHeat))
    dataHeatTB[ ,paralogNew := .N, by=c("geneID","supertaxon")]
    dataHeatTB <- data.frame(dataHeatTB[,c("geneID","supertaxon","paralogNew")])
    
    dataHeat <- merge(dataHeat,dataHeatTB,by=c('geneID','supertaxon'), all.x=TRUE)
    dataHeat$paralog <- dataHeat$paralogNew
    dataHeat <- dataHeat[!duplicated(dataHeat),]
    
    ### remove unneeded dots
    dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
    dataHeat$paralog[dataHeat$presSpec < 1] <- NA
    dataHeat$paralog[dataHeat$paralog == 1] <- NA
    
    p <- heatmap.plotting(dataHeat,input$x_axis,input$var1_id,input$var2_id,input$low_color_var1,input$high_color_var1,input$low_color_var2,input$high_color_var2,input$para_color,input$x_size,input$y_size,input$legend_size,input$main_legend,input$dot_zoom,input$x_angle,1)
    
    ### highlight taxon
    if(input$taxonHighlight != "none"){
      ## get selected highlight taxon ID
      rank_select = input$rank_select
      rankName = substr(rank_select,4,nchar(rank_select))   # get rank name from rank_select
      taxa_list <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
      taxonHighlight <- taxa_list$ncbiID[taxa_list$fullName == input$taxonHighlight & taxa_list$rank == rankName]
      if(length(taxonHighlight) == 0L){
        taxonHighlight <- taxa_list$ncbiID[taxa_list$fullName == input$taxonHighlight]
      }
      
      ## get taxonID together with it sorted index
      highlightTaxon <- toString(dataHeat[dataHeat$supertaxonID == taxonHighlight,2][1])
      
      ## get index
      selectedIndex = as.numeric(as.character(substr(highlightTaxon,2,4)))
      
      ## draw a rect to highlight this taxon's column
      if(input$x_axis == "taxa"){
        rect <- data.frame(xmin=selectedIndex-0.5, xmax=selectedIndex+0.5, ymin=-Inf, ymax=Inf)
      } else {
        rect <- data.frame(ymin=selectedIndex-0.5, ymax=selectedIndex+0.5, xmin=-Inf, xmax=Inf)
      }
      
      p = p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                        color="yellow",
                        alpha=0.3,
                        inherit.aes = FALSE)
    }
    
    ### highlight gene
    if(input$geneHighlight != "none"){
      ## get selected highlight gene ID
      geneHighlight <- input$geneHighlight
      
      ## get index
      allGenes <- levels(dataHeat$geneID)
      selectedIndex = match(geneHighlight,allGenes)
      
      ## draw a rect to highlight this taxon's column
      if(input$x_axis == "taxa"){
        rect <- data.frame(ymin=selectedIndex-0.5, ymax=selectedIndex+0.5, xmin=-Inf, xmax=Inf)
      } else {
        rect <- data.frame(xmin=selectedIndex-0.5, xmax=selectedIndex+0.5, ymin=-Inf, ymax=Inf)
      }
      
      p = p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                        color="yellow",
                        alpha=0.3,
                        inherit.aes = FALSE)
    }
    
    ### do plotting
    if(input$auto_update == FALSE){
      # Add dependency on the update button (only update when button is clicked)
      input$updateBtn
      
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
  
  ########### plot profile into plot.ui
  output$mainPlot <- renderPlot({
    if(input$auto_update == FALSE){
      # Add dependency on the update button (only update when button is clicked)
      input$updateBtn
      
      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        mainPlot()
      })
    } else {
      mainPlot()
    }
  })
  
  output$plot.ui <- renderUI({
    # show beschreibung file if no plot present
    if(v$doPlot == FALSE){
      return()
    } else{
      ## if auto_update is NOT selected, use updateBtn to trigger plot changing
      if(input$auto_update == FALSE){
        # Add dependency on the update button (only update when button is clicked)
        input$updateBtn
        
        # Add all the filters to the data based on the user inputs
        # wrap in an isolate() so that the data won't update every time an input
        # is changed
        isolate({
          withSpinner(
            plotOutput("mainPlot",width=input$width,height = input$height,
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
      ## if auto_update is true
      else {
        withSpinner(
          plotOutput("mainPlot",width=input$width,height = input$height,
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
  
  ########### download main plot
  output$plot_download <- downloadHandler(
    filename = function() {c("plot.pdf")},
    content = function(file) {
      ggsave(file, plot = mainPlot(), width = input$width*0.056458333, height = input$height*0.056458333, units="cm", dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
  ######## get info clicked point on main heatmap
  mainpoint_info <- reactive({
    ### check input
    if (v$doPlot == FALSE) return()
    
    # get selected supertaxon name
    taxa_list <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
    rank_select = input$rank_select
    rankName = substr(rank_select,4,nchar(rank_select))
    in_select <- as.numeric(taxa_list$ncbiID[taxa_list$fullName == input$in_select])
    
    dataHeat <- dataHeat()
    if(input$apply_cluster == TRUE){
      dataHeat <- clusteredDataHeat()
    }
    
    ### get values
    if (is.null(input$plot_click$x)) {return()}
    else{
      ### get cooridiate point
      if(input$x_axis == "genes"){
        corX = round(input$plot_click$y);
        corY = round(input$plot_click$x)
      } else {
        corX = round(input$plot_click$x);
        corY = round(input$plot_click$y)
      }
      
      # get geneID
      # genes <- as.matrix(dataHeat[dataHeat$supertaxonID == in_select & !is.na(dataHeat$presSpec),])
      genes <- levels(dataHeat$geneID)
      geneID <- toString(genes[corY])
      
      # get supertaxon (spec)
      supertaxa <- levels(dataHeat$supertaxon)
      spec <- toString(supertaxa[corX])
      
      # get var1, percentage of present species and var2 score
      var1 <- NA
      if(!is.na(dataHeat$var1[dataHeat$geneID == geneID & dataHeat$supertaxon == spec][1])){
        var1 <- max(na.omit(dataHeat$var1[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]))
      }
      Percent <- NA
      if(!is.na(dataHeat$presSpec[dataHeat$geneID == geneID & dataHeat$supertaxon == spec][1])){
        Percent <- max(na.omit(dataHeat$presSpec[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]))
      }
      var2 <- NA
      if(!is.na(dataHeat$var2[dataHeat$geneID == geneID & dataHeat$supertaxon == spec][1])){
        var2 <- max(na.omit(dataHeat$var2[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]))
      }
      
      # get ortholog ID
      orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      if(length(orthoID) > 1){
        orthoID = paste0(orthoID[1],",...")
      }
      # ### get list of all geneID that have the same ortholog
      # geneMatch <- dataHeat$geneID[dataHeat$orthoID == toString(orthoID)]
      # geneMatch <- geneMatch[!is.na(geneMatch)]
      # # list of all available geneID
      # geneList <- preData()
      # geneList$geneID <- as.factor(geneList$geneID)
      # allGenes <- as.list(levels(geneList$geneID))
      # # get index of all matched genes (genes have the same ortholog)
      # pos <- which(allGenes %in% geneMatch)
      # pos <- paste(pos, collapse=',')
      
      ### return info of clicked point
      if(is.na(as.numeric(Percent))){return()}
      else{
        # info <- c(geneID,as.character(orthoID),as.character(spec),round(as.numeric(var1),2),round(as.numeric(Percent),2),round(as.numeric(var2),2),pos)
        info <- c(geneID,as.character(orthoID),as.character(spec),round(as.numeric(var1),2),round(as.numeric(Percent),2),round(as.numeric(var2),2))
        return(info)
      }
    }
  })
  
  # ######## get list of same orthologs (hit_IDs) of a selected point
  # sameOrthoIndex <- reactive({
  #   ### check input
  #   if (v$doPlot == FALSE) return()
  #
  #   ### info
  #   info <- mainpoint_info()
  #   pos <- info[7]
  # })
  
  ###################################################################
  #### PLOT var1/var2 SCORE & % OF PRESENT SPECIES DISTRIBUTION #####
  ###################################################################
  
  ######## list of available variables for distribution plot
  output$selected.distribution = renderUI({
    if(nchar(input$var1_id) == 0 & nchar(input$var2_id) == 0){
      varList <- "% present taxa"
    } else if(nchar(input$var1_id) == 0 & nchar(input$var2_id) > 0){
      varList <- as.list(c(input$var2_id,"% present taxa"))
    } else if(nchar(input$var1_id) > 0 & nchar(input$var2_id) == 0){
      varList <- as.list(c(input$var1_id,"% present taxa"))
    } else {
      varList <- as.list(c(input$var1_id,input$var2_id,"% present taxa"))
    }
    
    selectInput('selected_dist','Choose variable to plot:',varList,varList[1])
  })
  
  ###### var1 / var2 distribution data
  distDf <- reactive({
    if (v$doPlot == FALSE) return()
    
    # open main input file
    # if(input$demo == TRUE){
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      if(input$demo_data == "demo"){
        dataOrig <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.long", sep="\t", header=T, fill=T, stringsAsFactors = FALSE)
      } else {
        dataOrig <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile", sep="\t", header=T, fill=T, stringsAsFactors = FALSE)
      }
      colnames(dataOrig) <- c("geneID","ncbiID","orthoID","var1","var2")
      splitDt <- dataOrig[,c("orthoID","var1","var2")]
    } else {
      filein <- input$mainInput
      
      inputType <- checkInputVadility(filein)
      if(inputType == "xml"){
        dataOrig <- xmlParser(filein$datapath)
        if(ncol(dataOrig) < 4){
          colnames(dataOrig) <- c("geneID","ncbiID","orthoID")
          splitDt <- dataOrig[,c("orthoID")]
        } else if(ncol(dataOrig) < 5){
          colnames(dataOrig) <- c("geneID","ncbiID","orthoID","var1")
          splitDt <- dataOrig[,c("orthoID","var1")]
        } else {
          colnames(dataOrig) <- c("geneID","ncbiID","orthoID","var1","var2")
          splitDt <- dataOrig[,c("orthoID","var1","var2")]
        }
      } else if(inputType == "fasta"){
        dataOrig <- fastaParser(filein$datapath)
        if(ncol(dataOrig) < 4){
          colnames(dataOrig) <- c("geneID","ncbiID","orthoID")
          splitDt <- dataOrig[,c("orthoID")]
        } else if(ncol(dataOrig) < 5){
          colnames(dataOrig) <- c("geneID","ncbiID","orthoID","var1")
          splitDt <- dataOrig[,c("orthoID","var1")]
        } else {
          colnames(dataOrig) <- c("geneID","ncbiID","orthoID","var1","var2")
          splitDt <- dataOrig[,c("orthoID","var1","var2")]
        }
      } else if(inputType == "long"){
        dataOrig <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
        if(ncol(dataOrig) < 4){
          colnames(dataOrig) <- c("geneID","ncbiID","orthoID")
          splitDt <- dataOrig[,c("orthoID")]
        } else if(ncol(dataOrig) < 5){
          colnames(dataOrig) <- c("geneID","ncbiID","orthoID","var1")
          splitDt <- dataOrig[,c("orthoID","var1")]
        } else {
          colnames(dataOrig) <- c("geneID","ncbiID","orthoID","var1","var2")
          splitDt <- dataOrig[,c("orthoID","var1","var2")]
        }
      } else if(inputType == "wide"){
        dataOrig <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
        # convert into paired columns
        mdData <- melt(dataOrig,id="geneID")
        
        # replace NA value with "NA#NA" (otherwise the corresponding orthoID will be empty)
        mdData$value <- as.character(mdData$value)
        mdData$value[is.na(mdData$value)] <- "NA#NA"
        
        # split "orthoID#var1#var2" into 3 columns
        splitDt <- as.data.frame(str_split_fixed(mdData$value, '#', 3))
        colnames(splitDt) <- c("orthoID","var1","var2")
      }
    }
    
    splitDt$orthoID[splitDt$orthoID == "NA" | is.na(splitDt$orthoID)] <- NA
    splitDt <- splitDt[complete.cases(splitDt),]
    
    if(length(levels(as.factor(splitDt$var2))) == 1){
      if(levels(as.factor(splitDt$var2)) == ""){
        splitDt$var2 <- 0
      }
    }
    
    # convert factor into numeric for "var1" & "var2" column
    if("var1" %in% colnames(splitDt)){
      splitDt$var1 <- suppressWarnings(as.numeric(as.character(splitDt$var1)))
      # filter splitDt based on selected var1 cutoff
      splitDt <- splitDt[splitDt$var1 >= input$var1[1] & splitDt$var1 <= input$var1[2],]
    }
    if("var2" %in% colnames(splitDt)){
      splitDt$var2 <- suppressWarnings(as.numeric(as.character(splitDt$var2)))
      # filter splitDt based on selected var2 cutoff
      splitDt <- splitDt[splitDt$var2 >= input$var2[1] & splitDt$var2 <= input$var2[2],]
    }
    
    # filter data base on customized plot (if chosen)
    if(input$dataset.distribution == "Customized data"){
      # get geneID and supertaxon name for splitDt
      allData <- get_data_filtered()
      splitDtName <- merge(splitDt, allData, by = "orthoID", all.x = TRUE)
      splitDtName$supertaxonMod <- substr(splitDtName$supertaxon,6,nchar(as.character(splitDtName$supertaxon)))
      splitDtName <- subset(splitDtName, select=c(orthoID,var1.x,var2.y,supertaxonMod,geneID))
      colnames(splitDtName) <- c("orthoID","var1","var2","supertaxonMod","geneID")
      
      # filter
      if(input$inTaxa[1] == "all" & input$inSeq[1] != "all"){
        splitDt <- subset(splitDtName,geneID %in% input$inSeq) ##### <=== select data from dataHeat for selected sequences only
      } else if(input$inSeq[1] == "all" & input$inTaxa[1] != "all"){
        splitDt <- subset(splitDtName,supertaxonMod %in% input$inTaxa) ##### <=== select data from dataHeat for selected taxa only
      } else {
        splitDt <- subset(splitDtName,geneID %in% input$inSeq & supertaxonMod %in% input$inTaxa) ##### <=== select data from dataHeat for selected sequences and taxa
      }
    }
    
    # return dt
    return(splitDt)
  })
  
  ###### var1 score distribution plot
  var1DistPlot <- function(){
    if (v$doPlot == FALSE) return()
    
    splitDt <- distDf()
    
    if(is.null(levels(as.factor(splitDt$var1)))){return()}
    else{
      splitDt <- splitDt[!is.na(splitDt$var1),]
      splitDt$mean <- mean(splitDt$var1)
      
      # plot var1 score distribution
      p <- ggplot(splitDt, aes(x=var1)) +
        geom_histogram(binwidth=.01, alpha=.5, position="identity") +
        geom_vline(data=splitDt, aes(xintercept=splitDt$mean,colour="red"),
                   linetype="dashed", size=1) +
        #      ggtitle(paste("Mean",input$var1_id,"=",round(mean(splitDt$var1),3))) +
        theme_minimal()
      p <- p + theme(legend.position = "none",
                     #                   plot.title = element_text(hjust = 0.5),
                     axis.title = element_text(size=input$dist_text_size),axis.text = element_text(size=input$dist_text_size)) +
        labs(x = paste0(input$var1_id," (mean = ",round(mean(splitDt$var1),3),")"), y = "Frequency")
      p
    }
  }
  
  output$var1DistPlot <- renderPlot(width = 512, height = 356,{
    if(input$auto_update == FALSE){
      # Add dependency on the update button (only update when button is clicked)
      input$updateBtn
      
      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        var1DistPlot()
      })
    } else {
      var1DistPlot()
    }
  })
  
  ###### var2 score distribution plot
  var2DistPlot <- function(){
    if (v$doPlot == FALSE) return()
    
    splitDt <- distDf()
    if(is.null(levels(as.factor(splitDt$var2)))){return()}
    else{
      splitDt <- splitDt[!is.na(splitDt$var2),]
      splitDt$mean <- mean(splitDt$var2)
      
      # plot var1 score distribution
      p <- ggplot(splitDt, aes(x=var2)) +
        geom_histogram(binwidth=.01, alpha=.5, position="identity") +
        geom_vline(data=splitDt, aes(xintercept=splitDt$mean,colour="red"),
                   linetype="dashed", size=1) +
        #      ggtitle(paste("Mean",input$var2_id,"=",round(mean(splitDt$var2),3))) +
        theme_minimal()
      p <- p + theme(legend.position = "none",
                     #                   plot.title = element_text(size=input$legend_size),#hjust = 0.5,
                     axis.title = element_text(size=input$dist_text_size),axis.text = element_text(size=input$dist_text_size)) +
        labs(x = paste0(input$var2_id," (mean = ",round(mean(splitDt$var2),3),")"), y = "Frequency")
      p
    }
  }
  
  output$var2DistPlot <- renderPlot(width = 512, height = 356,{
    if(input$auto_update == FALSE){
      # Add dependency on the update button (only update when button is clicked)
      input$updateBtn
      
      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        var2DistPlot()
      })
    } else {
      var2DistPlot()
    }
  })
  
  ####### % present species distribution plot
  
  ## calculate % present species for input file
  presSpecAllDt <- reactive({
    # open main input file
    # if(input$demo == TRUE){
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      if(input$demo_data == "demo"){
        mdData <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.long", sep="\t", header=T, fill=T, stringsAsFactors = FALSE)
      } else {
        mdData <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile", sep="\t", header=T, fill=T, stringsAsFactors = FALSE)
      }
      colnames(mdData) <- c("geneID","ncbiID","orthoID","var1","var2")
    } else {
      filein <- input$mainInput
      
      inputType <- checkInputVadility(filein)
      if(inputType == "xml"){
        mdData <- xmlParser(filein$datapath)
        if(ncol(mdData) < 4){
          colnames(mdData) <- c("geneID","ncbiID","orthoID")
        } else if(ncol(mdData) < 5){
          colnames(mdData) <- c("geneID","ncbiID","orthoID","var1")
        } else {
          colnames(mdData) <- c("geneID","ncbiID","orthoID","var1","var2")
        }
      } else if(inputType == "fasta"){
        mdData <- fastaParser(filein$datapath)
        if(ncol(mdData) < 4){
          colnames(mdData) <- c("geneID","ncbiID","orthoID")
        } else if(ncol(mdData) < 5){
          colnames(mdData) <- c("geneID","ncbiID","orthoID","var1")
        } else {
          colnames(mdData) <- c("geneID","ncbiID","orthoID","var1","var2")
        }
      }else if(inputType == "long"){
        mdData <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
        if(ncol(mdData) < 4){
          colnames(mdData) <- c("geneID","ncbiID","orthoID")
        } else if(ncol(mdData) < 5){
          colnames(mdData) <- c("geneID","ncbiID","orthoID","var1")
        } else {
          colnames(mdData) <- c("geneID","ncbiID","orthoID","var1","var2")
        }
      } else if(inputType == "wide"){
        data <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
        # convert into paired columns
        mdData <- melt(data,id="geneID")
        
        # split value column into orthoID, var1 & var2
        splitDt <- (str_split_fixed(mdData$value, '#', 3))
        # then join them back to mdData
        mdData <- cbind(mdData,splitDt)
        # rename columns
        colnames(mdData) <- c("geneID","ncbiID","value","orthoID","var1","var2")
        mdData <- mdData[,c("geneID","ncbiID","orthoID","var1","var2")]
      }
    }
    
    ### count number of inparalogs
    paralogCount <- plyr::count(mdData,c('geneID','ncbiID'))
    mdData <- merge(mdData,paralogCount,by=c('geneID','ncbiID'))
    colnames(mdData)[ncol(mdData)] <- "paralog"
    
    ### (3) GET SORTED TAXONOMY LIST (3) ###
    taxa_list <- sortedtaxa_list()
    
    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxa_list,'supertaxon')
    
    # merge mdData, mdDatavar2 and taxa_list to get taxonomy info
    taxaMdData <- merge(mdData,taxa_list,by='ncbiID')
    if("var1" %in% colnames(taxaMdData)){
      taxaMdData$var1 <- suppressWarnings(as.numeric(as.character(taxaMdData$var1)))
    }
    if("var2" %in% colnames(taxaMdData)){
      taxaMdData$var2 <- suppressWarnings(as.numeric(as.character(taxaMdData$var2)))
    }
    
    # calculate % present species
    finalPresSpecDt <- calcPresSpec(taxaMdData, taxaCount)
    
    finalPresSpecDt
  })
  
  ## % presSpec distribution plot
  presSpecPlot <- function(){
    if (v$doPlot == FALSE) return()
    
    # data
    dt <- presSpecAllDt()
    # remove presSpec < cutoff_min or > cutoff_max
    if(input$percent[1] > 0){
      dt <- dt[dt$presSpec >= input$percent[1] & dt$presSpec <= input$percent[2],]
    } else {
      if(input$percent[2] > 0){
        dt <- dt[dt$presSpec > 0 & dt$presSpec <= input$percent[2],]
      } else {
        dt <- dt[dt$presSpec > 0,]
      }
    }
    
    # calculate mean presSpec score
    dt$mean <- mean(dt$presSpec)
    
    # plot presSpec distribution
    p <- ggplot(dt, aes(x=presSpec)) +
      geom_histogram(binwidth=.01, alpha=.5, position="identity") +
      geom_vline(data=dt, aes(xintercept=dt$mean,colour="red"),
                 linetype="dashed", size=1) +
      #      ggtitle(paste("Mean % present taxa = ",round(mean(dt$presSpec),3))) +
      theme_minimal()
    p <- p + theme(legend.position = "none",
                   #                   plot.title = element_text(hjust = 0.5),
                   axis.title = element_text(size=input$dist_text_size),axis.text = element_text(size=input$dist_text_size)) +
      labs(x = paste0("% present taxa (mean = ",round(mean(dt$presSpec),3),")"), y = "Frequency")
    p
  }
  
  output$presSpecPlot <- renderPlot(width = 512, height = 356,{
    if(input$auto_update == FALSE){
      # Add dependency on the update button (only update when button is clicked)
      input$updateBtn
      
      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        presSpecPlot()
      })
    } else {
      presSpecPlot()
    }
  })
  
  ######## render dist_plot.ui
  output$dist_plot.ui <- renderUI({
    if(v$doPlot == FALSE){
      return()
    } else{
      if(is.null(input$selected_dist)){
        return()
      } else {
        if(input$selected_dist == "% present taxa"){
          withSpinner(plotOutput("presSpecPlot",width=input$width,height = input$height))
        } else{
          if(input$selected_dist == input$var1_id){
            withSpinner(plotOutput("var1DistPlot",width=input$width,height = input$height))
          } else if(input$selected_dist == input$var2_id){
            withSpinner(plotOutput("var2DistPlot",width=input$width,height = input$height))
          }
        }
      }
    }
  })
  
  ######## Download distribution plot
  output$plot_download_dist <- downloadHandler(
    filename = function() {paste0("distributionPlot.pdf")},
    content = function(file) {
      if(input$selected_dist == input$var1_id){
        ggsave(file, plot = var1DistPlot(), dpi=300, device = "pdf", limitsize=FALSE)
      }
      if(input$selected_dist == input$var2_id){
        ggsave(file, plot = var2DistPlot(), dpi=300, device = "pdf", limitsize=FALSE)
      }
      if(input$selected_dist == "% present taxa"){
        ggsave(file, plot = presSpecPlot(), dpi=300, device = "pdf", limitsize=FALSE)
      }
    }
  )
  
  #############################################################
  ################# PLOT CUSTOMIZED PROFILE ###################
  #############################################################
  
  ####### check if all genes and all species are selected
  output$same_profile <- reactive({
    if (v$doPlot == FALSE) return(FALSE)
    if(length(input$inSeq[1]) == 0){ return(FALSE)}
    else{
      if(input$inSeq[1] == "all" & input$inTaxa[1] == "all"){return(TRUE)}
    }
  })
  outputOptions(output, 'same_profile', suspendWhenHidden=FALSE)
  
  ######## change label of plotCustom button if auto_update_selected is unchecked
  output$plot_custom_btn <- renderUI({
    if(input$auto_update_selected == FALSE){
      shinyBS::bsButton("plotCustom", "Plot/Update selected sequence(s)/taxa",style="warning")
    } else {
      shinyBS::bsButton("plotCustom", "Plot selected sequence(s)/taxa",style="warning")
    }
  })
  
  ######## check if button is clicked
  vCt <- reactiveValues(doPlotCustom = FALSE)
  observeEvent(input$plotCustom, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    vCt$doPlotCustom <- input$plotCustom
    filein <- input$mainInput
    # if(input$demo == TRUE){ filein = 1 }
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){ filein = 1 }
    if(is.null(filein)){vCt$doPlotCustom <- FALSE}
  })
  
  ######## print list of available customized taxonomy ranks (the lowest rank is the same as the chosen main rank)
  output$rank_select_cus = renderUI({
    mainRank <- input$rank_select
    mainChoices = get_taxonomy_ranks()
    cusChoices <- mainChoices[mainChoices >= mainRank]
    
    selectInput("rank_select_cus", label = h5("Select taxonomy rank:"),
                choices = as.list(cusChoices),
                selected = mainRank)
  })
  
  ######## print list of available taxa for customized plot (based on rank from rank_select_cus)
  taxa_select_cus <- reactive({
    rank_select_cus = input$rank_select_cus
    
    if(length(rank_select_cus) == 0){return()}
    else{
      ### load list of unsorted taxa
      Dt <- get_taxa_list(TRUE)
      
      ### load list of taxon name
      nameList <- get_name_list()
      
      rankName = substr(rank_select_cus,4,nchar(rank_select_cus))   # get rank name from rank_select
      choice <- as.data.frame
      choice <- rbind(Dt[rankName])
      colnames(choice) <- "ncbiID"
      choice <- merge(choice,nameList,by="ncbiID",all = FALSE)
      return(choice)
    }
  })
  
  output$taxa_select_cus = renderUI({
    choice <- taxa_select_cus()
    choice$fullName <- as.factor(choice$fullName)
    selectInput('taxa_select_cus',h5('Choose (super)taxon of interest:'),as.list(levels(choice$fullName)),levels(choice$fullName)[1])
  })
  
  ######## get list of taxa based on selected taxa_select_cus
  cus_taxaName <- reactive({
    
    taxa_select_cus = input$taxa_select_cus
    rankName = substr(input$rank_select_cus,4,nchar(input$rank_select_cus))
    
    if(taxa_select_cus == ""){return()}
    
    ### load list of unsorted taxa
    Dt <- get_taxa_list(TRUE)

    
    ### get ID of customized (super)taxon
    taxa_list <- get_name_list(FALSE, FALSE)
    superID <- taxa_list$ncbiID[taxa_list$fullName == taxa_select_cus & taxa_list$rank %in% c(rankName,"norank")]
    
    ### from that ID, get list of all taxa for main selected taxon
    mainRankName = substr(input$rank_select,4,nchar(input$rank_select))
    customizedtaxa_id <- levels(as.factor(Dt[mainRankName][Dt[rankName] == superID,]))
    
    cus_taxaName <- taxa_list$fullName[taxa_list$rank %in% c(mainRankName,"norank") & taxa_list$ncbiID %in% customizedtaxa_id]
    
    return(cus_taxaName)
  })
  
  ######## create plot (same as main plot)
  selected_plot <- function(){
    if (vCt$doPlotCustom == FALSE) return()
    if(input$inSeq[1] == "all" & input$inTaxa[1] == "all") return()
    else{
      dataHeat <- dataHeat()
      
      ### cluster dataHeat (if selected)
      if(input$apply_cluster == TRUE){
        dataHeat <- clusteredDataHeat()
      }
      
      ### process data
      dataHeat$supertaxonMod <- substr(dataHeat$supertaxon,6,nchar(as.character(dataHeat$supertaxon)))
      if(input$inTaxa[1] == "all" & input$inSeq[1] != "all"){
        dataHeat <- subset(dataHeat,geneID %in% input$inSeq) ##### <=== select data from dataHeat for selected sequences only
      } else if(input$inSeq[1] == "all" & input$inTaxa[1] != "all"){
        dataHeat <- subset(dataHeat,supertaxonMod %in% input$inTaxa) ##### <=== select data from dataHeat for selected taxa only
      } else {
        dataHeat <- subset(dataHeat,geneID %in% input$inSeq & supertaxonMod %in% input$inTaxa) ##### <=== select data from dataHeat for selected sequences and taxa
      }
      
      ### remove unneeded dots
      dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
      dataHeat$paralog[dataHeat$presSpec < 1] <- NA
      dataHeat$paralog[dataHeat$paralog == 1] <- NA
      
      ### create plot
      p <- heatmap.plotting(dataHeat,input$x_axis_selected,input$var1_id,input$var2_id,input$low_color_var1,input$high_color_var1,input$low_color_var2,input$high_color_var2,input$para_color,input$x_size_select,input$y_size_select,input$legend_size_select,input$selected_legend,input$dot_zoom_select,input$x_angle_select,0)
      
      ### do plotting
      if(input$auto_update_selected == FALSE){
        # Add dependency on the update button (only update when button is clicked)
        input$plotCustom
        
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
  
  output$selected_plot <- renderPlot({
    if(input$auto_update_selected == FALSE){
      # Add dependency on the update button (only update when button is clicked)
      input$plotCustom
      
      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        selected_plot()
      })
    } else {
      selected_plot()
    }
  })
  
  ######## plot selected sequences heatmap
  output$selected_plot.ui <- renderUI({
    if(is.null(input$inSeq[1]) | is.null(input$inTaxa[1])){ return()}
    else if(input$inSeq[1] == "all" & input$inTaxa[1]=="all"){return()}
    else{
      if(input$auto_update_selected == FALSE){
        # Add dependency on the update button (only update when button is clicked)
        input$plotCustom
        
        # Add all the filters to the data based on the user inputs
        # wrap in an isolate() so that the data won't update every time an input
        # is changed
        isolate({
          withSpinner(
            plotOutput("selected_plot",width=input$selected_width,height = input$selected_height,
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
          plotOutput("selected_plot",width=input$selected_width,height = input$selected_height,
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
  
  ######## download selected plot
  output$selected_download <- downloadHandler(
    filename = function() {c("selected_plot.pdf")},
    content = function(file) {
      ggsave(file, plot = selected_plot(), width = input$selected_width*0.056458333, height = input$selected_height*0.056458333, units="cm", dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
  ######## get info of a clicked point on selected plot (also the same as get info from main plot)
  selectedpoint_info <- reactive({
    ### check input
    if (vCt$doPlotCustom == FALSE) return()
    
    # get selected supertaxon name
    taxa_list <- get_name_list(FALSE,FALSE)
    rank_select = input$rank_select
    rankName = substr(rank_select,4,nchar(rank_select))
    in_select <- as.numeric(taxa_list$ncbiID[taxa_list$fullName == input$in_select])
    
    dataHeat <- dataHeat()
    if(input$apply_cluster == TRUE){
      dataHeat <- clusteredDataHeat()
    }
    
    ### get sub-dataframe of selected taxa and sequences
    dataHeat$supertaxonMod <- substr(dataHeat$supertaxon,6,nchar(as.character(dataHeat$supertaxon)))
    if(input$inTaxa[1] == "all" & input$inSeq[1] != "all"){
      dataHeat <- subset(dataHeat,geneID %in% input$inSeq) ##### <=== select data from dataHeat for selected sequences only
    } else if(input$inSeq[1] == "all" & input$inTaxa[1] != "all"){
      dataHeat <- subset(dataHeat,supertaxonMod %in% input$inTaxa) ##### <=== select data from dataHeat for selected taxa only
    } else {
      dataHeat <- subset(dataHeat,geneID %in% input$inSeq & supertaxonMod %in% input$inTaxa) ##### <=== select data from dataHeat for selected sequences and taxa
    }
    
    ### drop all other supertaxon that are not in sub-dataframe
    dataHeat$supertaxon <- factor(dataHeat$supertaxon)
    dataHeat$geneID <- factor(dataHeat$geneID)
    
    ### get values
    if (is.null(input$plot_click_selected$x)) return()
    else{
      ### get cooridiate point
      if(input$x_axis_selected == "genes"){
        corX = round(input$plot_click_selected$y);
        corY = round(input$plot_click_selected$x)
      } else {
        corX = round(input$plot_click_selected$x);
        corY = round(input$plot_click_selected$y)
      }
      
      # get geneID
      genes <- levels(dataHeat$geneID)
      geneID <- toString(genes[corY])
      # get supertaxon (spec)
      supertaxa <- levels(dataHeat$supertaxon)
      spec <- toString(supertaxa[corX])
      # get var1, percentage of present species and var2 score
      var1 <- NA
      if(!is.na(dataHeat$var1[dataHeat$geneID == geneID & dataHeat$supertaxon == spec][1])){
        var1 <- max(na.omit(dataHeat$var1[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]))
      }
      Percent <- NA
      if(!is.na(dataHeat$presSpec[dataHeat$geneID == geneID & dataHeat$supertaxon == spec][1])){
        Percent <- max(na.omit(dataHeat$presSpec[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]))
      }
      var2 <- NA
      if(!is.na(dataHeat$var2[dataHeat$geneID == geneID & dataHeat$supertaxon == spec][1])){
        var2 <- max(na.omit(dataHeat$var2[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]))
      }
      
      # get ortholog ID
      orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      if(length(orthoID) > 1){
        orthoID = paste0(orthoID[1],",...")
      }
      
      if(is.na(as.numeric(Percent))){return()}
      else{
        info <- c(geneID,as.character(orthoID),as.character(spec),round(as.numeric(var1),2),round(as.numeric(Percent),2),round(as.numeric(var2),2))
      }
    }
  })
  
  
  #############################################################
  ################### SHOW CLICKED POINT INFO #################
  #############################################################
  
  ### get value of point_info for activating Detailed Plot button
  output$point_info_status <- reactive({
    if(input$tabs == 'Main profile'){
      info <- mainpoint_info()  # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
    } else if(input$tabs=='Customized profile'){
      info <- selectedpoint_info()
    } else {
      info <- NULL
    }
    is.null(info)
  })
  outputOptions(output, "point_info_status", suspendWhenHidden = FALSE)
  
  ######## show info into "point's info" box
  output$point_info <- renderText({
    ##### GET INFO BASED ON CURRENT TAB
    if(input$tabs == 'Main profile'){
      info <- mainpoint_info()  # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
    } else if(input$tabs=='Customized profile'){
      info <- selectedpoint_info()
    } else {
      return ()
    }
    
    if(is.null(info)){return()}
    else{
      orthoID <- info[2]
      # ## parse orthoID for oneSeq
      # if(input$input_type == 'Concatenated fasta file'){
      #   orthoIDTmp <- unlist(strsplit(toString(info[2]),"\\|"))
      #   #orthoID = toString(paste0(orthoIDTmp[2],":",orthoIDTmp[3]))
      #   orthoID = toString(orthoIDTmp[3])
      # }
      if(is.na(orthoID)){return()}
      else{
        # if(orthoID=="NA"){orthoID <- info[2]}
        
        ## print output
        a <- toString(paste("Seed-ID:",info[1]))
        b <- toString(paste0("Hit-ID: ",orthoID," (",substr(info[3],6,nchar(info[3])),")"))
        c <- ""
        if(input$var1_id != ""){
          c <- toString(paste(input$var1_aggregate_by,input$var1_id,":",info[4]))
        }
        d <- ""
        if(input$var2_id != ""){
          d <- toString(paste(input$var2_aggregate_by,input$var2_id,":",info[6]))
        }
        e <- toString(paste("% present taxa:",info[5]))
        paste(a,b,c,d,e,sep="\n")
      }
    }
  })
  
  #############################################################
  ##################### DETAILED PLOT #########################
  #############################################################
  
  ######## data for detailed plot
  detail_plotDt <- reactive({
    if (v$doPlot == FALSE) return()
    
    ##### GET INFO BASED ON CURRENT TAB
    if(input$tabs == 'Main profile'){
      info <- mainpoint_info()  # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
    } else if(input$tabs=='Customized profile'){
      info <- selectedpoint_info()
    }
    
    if(is.null(info)){return()}
    else{
      ### get info for present taxa in selected supertaxon (1)
      plotTaxon = info[3]
      plotGeneID = info[1]
      fullDf <- get_data_filtered()
      selDf <- as.data.frame(fullDf[fullDf$geneID == plotGeneID & fullDf$supertaxon == plotTaxon,])
      
      ### get all taxa of this supertaxon (2)
      allTaxaDf <- sortedtaxa_list()
      allTaxaDf <- allTaxaDf[allTaxaDf$supertaxon == plotTaxon,]
      allTaxaDf <- subset(allTaxaDf, select = c("abbrName","fullName"))
      
      ### merge (1) and (2) together
      joinedDf <- merge(selDf,allTaxaDf, by= c("abbrName"), all.y=TRUE)
      joinedDf <- subset(joinedDf, select = c("abbrName","fullName.y","geneID","orthoID","var1","var2"))
      names(joinedDf)[names(joinedDf) == 'fullName.y'] <- 'fullName'
      
      # replace var1/var2 as NA for all "NA orthologs"
      joinedDf$var1[is.na(joinedDf$orthoID)] <- NA
      joinedDf$var2[is.na(joinedDf$orthoID)] <- NA
      
      # remove NA orthologs if required
      if(input$detailed_remove_na == TRUE){
        joinedDf <- joinedDf[!is.na(joinedDf$orthoID),]
      }
      
      ### return data for detailed plot
      return(joinedDf)
    }
  })
  
  ######## render detailed plot
  detail_plot <- function(){
    if (v$doPlot == FALSE) return()
    
    selDf <- detail_plotDt()
    selDf$x_label <- paste(selDf$orthoID," (",selDf$fullName,")",sep = "")
    
    # if(input$detailed_remove_na == TRUE){
    #   selDf <- selDf[!is.na(selDf$orthoID),]
    # }
    
    ### create joined DF for plotting var1 next to var2
    var1Df <- subset(selDf, select = c("x_label","var1"))
    var1Df$type <- input$var1_id
    colnames(var1Df) <- c("id","score","var")
    
    var2Df <- subset(selDf, select = c("x_label","var2"))
    var2Df$type <- input$var2_id
    colnames(var2Df) <- c("id","score","var")
    
    detailedDf <- rbind(var1Df,var2Df)
    
    # remove ONE missing variable
    if(nlevels(as.factor(detailedDf$var)) > 1){
      detailedDf <- detailedDf[nchar(detailedDf$var)>0,]
    }
    
    ### keep order of ID (x_label)
    detailedDf$id <- factor(detailedDf$id, levels = unique(detailedDf$id))
    
    ### create plot
    gp = ggplot(detailedDf, aes(y=score,x=id,fill=var)) +
      geom_bar(stat="identity", position=position_dodge(),na.rm = TRUE) +
      coord_flip() +
      labs(x="") +
      labs(fill="") +
      theme_minimal()
    #geom_text(aes(label=var1), vjust=3)
    gp = gp+theme(axis.text.x = element_text(angle=90,hjust=1),
                  axis.text = element_text(size = input$detailed_text),
                  axis.title = element_text(size = input$detailed_text),
                  legend.text = element_text(size = input$detailed_text)
    )
    gp
  }
  
  output$detail_plot <- renderPlot({
    p <- detail_plot()
    p
  })
  
  ######## plot detailed bar chart
  output$detail_plot.ui <- renderUI({
    withSpinner(
      plotOutput("detail_plot",width=800,height = input$detailed_height,
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
  
  ######## download detailed plot
  output$download_detailed <- downloadHandler(
    filename = function() {c("detailedPlot.pdf")},
    content = function(file) {
      g <- detail_plot()
      ggsave(file, plot = g, width = 800*0.056458333, height = input$detailed_height*0.056458333, units="cm", dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
  ######## GET info when clicking on detailed plot
  point_infoDetail <- reactive({
    selDf <- detail_plotDt()
    # selDf <- selDf[complete.cases(selDf),]
    selDf$orthoID <- as.character(selDf$orthoID)
    # allOrthoID <- sort(selDf$orthoID)
    
    ### get coordinates of plot_click_detail
    if (is.null(input$plot_click_detail$x)) return()
    else{
      corX = round(input$plot_click_detail$y)
      corY = round(input$plot_click_detail$x)
    }
    
    ### get pair of sequence IDs & var1
    seedID <- as.character(selDf$geneID[!is.na(selDf$geneID)][1])
    orthoID <- as.character(selDf$orthoID[corX])
    
    var1 <- as.list(selDf$var1[selDf$orthoID==orthoID])
    var1 <- as.character(var1[!is.na(var1)])
    var2 <- as.list(selDf$var2[selDf$orthoID==orthoID])
    var2 <- as.character(var2[!is.na(var2)])
    # ncbiID <- as.character(selDf$abbrName[selDf$orthoID==orthoID])
    ncbiID <- selDf[selDf$orthoID == orthoID,]$abbrName
    ncbiID <- as.character(ncbiID[!is.na(ncbiID)][1])
    
    ### return info
    if(is.na(orthoID)){
      return(NULL)
    } else {
      if(orthoID != "NA"){
        info <- c(seedID,orthoID,var1,var2,ncbiID)
      }
    }
  })
  
  ### SHOW info when clicking on detailed plot
  output$detail_click <- renderText({
    info <- point_infoDetail() # info = seedID, orthoID, var1
    
    if(is.null(info)){paste("select ortholog")}
    else{
      a <- paste0("seedID = ",info[1])
      b <- paste0("hitID = ",info[2])
      c <- ""
      if(input$var1_id != ""){
        c <- paste0(input$var1_id," = ",info[3])
      }
      d <- ""
      if(input$var2_id != ""){
        d <- paste0(input$var2_id," = ",info[4])
      }
      paste(a,b,c,d,sep="\n")
    }
  })
  
  ##### FASTA OUTPUT
  fastaOutData <- function(dataOut){
    # dataOut <- as.data.frame(download_custom_data())
    fastaOutDf <- data.frame()
    
    ### check main input
    filein <- input$mainInput
    if(!is.null(filein)){inputType <- checkInputVadility(filein)}
    else{ inputType <- "NA"}
    
    ### get seqs for AMPK-TOR and microsporidia ONLINE demo data
    if(input$demo_data == "ampk-tor" | input$demo_data == "demo"){
      fastaUrl <- paste0("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/fasta_file/concatenatedSeq.fa")
      if(input$demo_data == "ampk-tor"){
        fastaUrl <- paste0("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.extended.fa")
      }
      
      if(url.exists(fastaUrl)){
        # load fasta file
        faFile <- as.data.frame(read.table(fastaUrl, sep="\t", header=F, fill=T, stringsAsFactors = FALSE, quote = ""))
        faDf <- data.frame("seqID" = faFile$V1[grepl(">",faFile$V1)], "seq" = faFile$V1[!grepl(">",faFile$V1)], stringsAsFactors=FALSE)
        
        # get sequences
        for(j in 1:nrow(dataOut)){
          seqID <- as.character(dataOut$orthoID[j])
          groupID <- as.character(dataOut$geneID[j])
          
          seq <- as.character(faDf$seq[faDf$seqID == paste0(">",seqID)])
          fastaOut <- paste(paste0(">",groupID,"|",seqID),seq,sep="\n")
          fastaOutDf <- rbind(fastaOutDf,as.data.frame(fastaOut))
        }
      } else {
        fastaOut <- paste0(fastaUrl," not found!!!")
        fastaOutDf <- rbind(fastaOutDf,as.data.frame(fastaOut))
      }
    }
    
    ### get seqs for fasta main input
    if(inputType == "fasta"){
      file <- filein$datapath
      fastaFile = readAAStringSet(file)
      
      seq_name = names(fastaFile)
      sequence = paste(fastaFile)
      fa <- data.frame(seq_name, sequence)  # data frame contains all sequences from input file
      
      for(j in 1:nrow(dataOut)){
        seqID <- paste0(as.character(dataOut$geneID[j]),"|ncbi",as.character(dataOut$ncbiID[j]),"|",as.character(dataOut$orthoID[j]))
        
        seq <- fa$sequence[pmatch(seqID,fa$seq_name)]
        
        if(length(seq[1]) < 1){
          fastaOut <- paste0(seqID," not found in ",file,"! Please check again!")
        } else{
          fastaOut <- paste(paste0(">",seqID),seq[1],sep="\n")
        }
        
        fastaOutDf <- rbind(fastaOutDf,as.data.frame(fastaOut))
      }
    }
    
    ### get seqs for extended.fa
    if(input$demo_data == "none" & input$input_type == 'Concatenated fasta file'){
      if(!is.null(input$oneSeqFasta)){
        fasIn <- input$oneSeqFasta
        file <- toString(fasIn$datapath)
        
        # read fasta file and save sequences into dataframe
        fastaFile = readAAStringSet(file)
        
        seq_name = names(fastaFile)
        sequence = paste(fastaFile)
        fa <- data.frame(seq_name, sequence)  # data frame contains all sequences from input file
        
        # get selected sequences
        for(j in 1:nrow(dataOut)){
          seqID <- as.character(dataOut$orthoID[j])
          groupID <- as.character(dataOut$geneID[j])
          seq <- fa$sequence[pmatch(seqID,fa$seq_name)]
          flag <- 1
          if(is.na(seq)){
            seqID <- paste0(as.character(dataOut$geneID[j]),"|ncbi",as.character(dataOut$ncbiID[j]),"|",as.character(dataOut$orthoID[j]))
            seq <- fa$sequence[pmatch(seqID,fa$seq_name)]
            flag <- 0
          }
          
          if(length(seq[1]) < 1){
            fastaOut <- paste0(seqID," not found in ",file,"! Please check the header format in FASTA file!")
          } else{
            if(!is.na(seq[1])){
              if(flag == 1){
                fastaOut <- paste(paste0(">",groupID,"|",seqID),seq[1],sep="\n")
              } else {
                fastaOut <- paste(paste0(">",seqID),seq[1],sep="\n")
              }
              
            } else {
              fastaOut <- paste0(seqID," not found in uploaded FASTA file!!! Please check again!!!")
            }
          }
          fastaOutDf <- rbind(fastaOutDf,as.data.frame(fastaOut))
        }
      } else {
        if(inputType != "fasta"){
          fastaOut <- paste0("Please provide FASTA file(s) in Input & settings page!")
          fastaOutDf <- rbind(fastaOutDf,as.data.frame(fastaOut))
        }
      }
    }
    
    ### get seqs for other cases (input offline fasta files in a folder)
    if(input$demo_data == "none" & input$input_type == "Fasta folder" & inputType != "fasta"){
      if(input$path != ""){
        # path, file format and fasta header format
        path = input$path
        dir_format = input$dir_format
        file_ext = input$file_ext
        id_format = input$id_format
        
        # get list of species IDs
        if(id_format == 1){
          specDf <- as.data.frame(str_split_fixed(strReverse(as.character(dataOut$orthoID)),":",2))
          specDf$specID <- strReverse(as.character(specDf$V2))
        } else if(id_format == 2){
          specDf <- as.data.frame(str_split_fixed(strReverse(as.character(dataOut$orthoID)),"@",2))
          specDf$specID <- strReverse(as.character(specDf$V2))
        } else if(id_format == 3){
          specDf <- as.data.frame(str_split_fixed(strReverse(as.character(dataOut$orthoID)),"|",2))
          specDf$specID <- strReverse(as.character(specDf$V2))
        }
        
        # read all specices FASTA files at once
        fa = data.frame()
        for(i in 1:length(levels(as.factor(specDf$specID)))){
          specID <- as.character(levels(as.factor(specDf$specID))[i])
          
          # full path fasta file
          file <- paste0(path,"/",specID,".",file_ext)
          if(dir_format == 2){
            file <- paste0(path,"/",specID,"/",specID,".",file_ext)
          }
          
          # read fasta file and save sequences into dataframe
          if(file.exists(file)){
            fastaFile = readAAStringSet(file)
            
            seq_name = names(fastaFile)
            sequence = paste(fastaFile)
            fa <- rbind(fa, data.frame(seq_name, sequence))  # data frame contains all sequences from input file
          }
        }
        
        # now get selected sequences
        if(nrow(fa) > 0){
          for(j in 1:nrow(dataOut)){
            seqID <- as.character(dataOut$orthoID[j])
            groupID <- as.character(dataOut$geneID[j])
            
            seq <- fa$sequence[pmatch(seqID,fa$seq_name)]
            
            if(length(seq[1]) < 1){
              fastaOut <- paste0(seqID," not found in ",file,"! Please check id_format in FASTA config again!")
            } else{
              fastaOut <- paste(paste0(">",groupID,"|",seqID),seq[1],sep="\n")
            }
            fastaOutDf <- rbind(fastaOutDf,as.data.frame(fastaOut))
          }
        } else {
          fastaOut <- paste0("No fasta file has been found in ",path,"!!! Please check the full path to FASTA folder and the id_format (header format) in FASTA config again!!!")
          fastaOutDf <- rbind(fastaOutDf,as.data.frame(fastaOut))
        }
        
      } else {
        fastaOut <- paste0("Please provide FASTA files in Input & settings page!")
        fastaOutDf <- rbind(fastaOutDf,as.data.frame(fastaOut))
      }
    }
    
    # remove duplicated sequences
    fastaOutDf <- fastaOutDf[!duplicated(fastaOutDf), ]
    
    return(fastaOutDf)
  }
  
  ######## FASTA sequence
  output$fasta <- renderText({
    if(v$doPlot == FALSE){return()}
    
    info <- point_infoDetail() # info = seedID, orthoID, var1
    if(is.null(info)){return()}
    else{
      data <- get_data_filtered()
      
      seqID <- toString(info[2])
      groupID <- toString(info[1])
      ncbiID <- gsub("ncbi","",toString(info[5]))
      
      seqDf <- data.frame("geneID" = groupID, "orthoID" = seqID, "ncbiID" = ncbiID)
      fastaOut <- fastaOutData(seqDf)
      
      return(paste(fastaOut[1]))
    }
  })
  
  
  #############################################################
  ################ FEATURE ARCHITECTURE PLOT ##################
  #############################################################
  
  ######## get domain file/path
  getDomainFile <- reactive({
    ### click info
    info <- point_infoDetail() # info = seedID, orthoID, var1
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    var1 <- as.character(info[3])
    
    ### domain file
    # if(input$demo == TRUE){
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      if(is.null(info)){
        fileDomain <- "noSelectHit"
        updateButton(session, "doDomainPlot", disabled = TRUE)
      } else {
        updateButton(session, "doDomainPlot", disabled = FALSE)
        if(input$demo_data == "demo"){
          fileDomain <- suppressWarnings(paste0("https://github.com/BIONF/phyloprofile-data/blob/master/demo/domain_files/",group,".domains?raw=true"))
        } else {
          fileDomain <- suppressWarnings(paste0("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.domains_F"))
        }
      }
    } else {
      if(input$anno_choose == "from file"){
        fileDomain <- input$fileDomainInput
        if(is.null(fileDomain)){
          fileDomain <- "noFileInput"
        } else {
          if(is.null(info)){
            fileDomain <- "noSelectHit"
            updateButton(session, "doDomainPlot", disabled = TRUE)
          } else {
            updateButton(session, "doDomainPlot", disabled = FALSE)
            fileDomain <- fileDomain$datapath
          }
        }
      } else {
        if(is.null(info)){
          fileDomain <- "noSelectHit"
          updateButton(session, "doDomainPlot", disabled = TRUE)
        } else {
          ### check file extension
          allExtension <- c("txt","csv","list","domains","architecture")
          flag <- 0
          for(i in 1:length(allExtension)){
            fileDomain <- paste0(input$domainPath,"/",group,".",allExtension[i])
            if(file.exists(fileDomain) == TRUE){
              updateButton(session, "doDomainPlot", disabled = FALSE)
              flag <- 1
              break()
            }
          }
          
          if(flag == 0){
            fileDomain <- "noFileInFolder"
            updateButton(session, "doDomainPlot", disabled = TRUE)
          }
        }
      }
    }
    
    return (fileDomain)
  })
  
  ######## check domain file
  output$check_domain_files <- renderUI({
    fileDomain <- getDomainFile()
    if(fileDomain == "noFileInput"){
      em("Domain file not provided!!")
    } else if(fileDomain == "noFileInFolder"){
      msg <- paste0(
        "<p><em>Domain file not found!! </em></p>
        <p><em>Please make sure that file name has to be in this format:
        <strong>&lt;seedID&gt;.extension</strong>, where extension is limited to
        <strong>txt</strong>, <strong>csv</strong>, <strong>list</strong>, <strong>domains</strong> or <strong>architecture</strong>.
        </em></p>"
      )
      HTML(msg)
    } else if(fileDomain == "noSelectHit"){
      em("Please select one ortholog sequence!!")
    }
  })
  
  ######## check clicked
  v3 <- reactiveValues(doPlot3 = FALSE)
  observeEvent(input$doDomainPlot, {
    v3$doPlot3 <- input$doDomainPlot
    filein <- input$mainInput
    # if(input$demo == TRUE){ filein = 1 }
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){ filein = 1 }
    if(is.null(filein)){v3$doPlot3 <- FALSE}
  })
  
  ######## create domain plot
  archi_plot <- function(){
    if (v3$doPlot3 == FALSE) return()
    
    ### info
    info <- point_infoDetail() # info = seedID, orthoID, var1
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    var1 <- as.character(info[3])
    
    ### parse domain file
    fileDomain <- getDomainFile()
    
    # if(input$demo == TRUE){
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      domainDf <- as.data.frame(read.csv(fileDomain, sep="\t", header=F, comment.char = "", stringsAsFactors = FALSE, quote = ""))
    } else {
      if(fileDomain != FALSE){
        domainDf <- as.data.frame(read.table(fileDomain, sep='\t',header=FALSE,comment.char=""))
      }
    }
    
    if(ncol(domainDf) == 6){
      colnames(domainDf) <- c("seedID","orthoID","length","feature","start","end")
    } else if(ncol(domainDf) == 7){
      colnames(domainDf) <- c("seedID","orthoID","length","feature","start","end","weight")
    } else if(ncol(domainDf) == 8){
      colnames(domainDf) <- c("seedID","orthoID","length","feature","start","end","weight","path")
    }
    
    domainDf$seedID <- gsub("\\|",":",domainDf$seedID)
    domainDf$orthoID <- gsub("\\|",":",domainDf$orthoID)
    

    ### get sub dataframe based on selected groupID and orthoID
    ortho <- gsub("\\|",":",ortho)
    grepID = paste(group,"#",ortho,sep="")
    subDomainDf <- domainDf[grep(grepID,domainDf$seedID),]
    subDomainDf$feature <- as.character(subDomainDf$feature)
    
    if(nrow(subDomainDf) < 1){
      v3$doPlot3 = FALSE
      return()
    } else {
      ### ortho domains df
      orthoDf <- filter(subDomainDf,orthoID==ortho)
      
      ### seed domains df
      seedDf <- filter(subDomainDf,orthoID != ortho)
      if(nrow(seedDf) == 0){seedDf <- orthoDf}
      
      seed = as.character(seedDf$orthoID[1])
      
      ### change order of one dataframe's features based on order of other df's features
      if(length(orthoDf$feature) < length(seedDf$feature)){
        orderedOrthoDf <- orthoDf[order(orthoDf$feature), ]
        orderedSeedDf <- sortDomains(orderedOrthoDf, seedDf)
      } else {
        orderedSeedDf <- seedDf[order(seedDf$feature), ]
        orderedOrthoDf <- sortDomains(orderedSeedDf, orthoDf)
      }
      
      ### join weight values and feature names
      if("weight" %in% colnames(orderedOrthoDf)){
        orderedOrthoDf$yLabel <- paste0(orderedOrthoDf$feature," (",round(orderedOrthoDf$weight,2),")")
        orderedOrthoDf$feature <- orderedOrthoDf$yLabel
      }
      if("weight" %in% colnames(orderedSeedDf)){
        orderedSeedDf$yLabel <- paste0(orderedSeedDf$feature," (",round(orderedSeedDf$weight,2),")")
        orderedSeedDf$feature <- orderedSeedDf$yLabel
      }
      
      ### plotting
      sep = ":"
      if(!is.null(input$oneSeqFasta)){sep="|"}
      
      plot_ortho <- domain.plotting(orderedOrthoDf,ortho,sep,input$label_archi_size,input$title_archi_size,min(subDomainDf$start),max(c(subDomainDf$end,subDomainDf$length)))
      plot_seed <- domain.plotting(orderedSeedDf,seed,sep,input$label_archi_size,input$title_archi_size,min(subDomainDf$start),max(c(subDomainDf$end,subDomainDf$length)))
      
      # grid.arrange(plot_seed,plot_ortho,ncol=1)
      if(ortho == seed){
        arrangeGrob(plot_seed,ncol=1)
      } else {
        seedHeight = length(levels(as.factor(orderedSeedDf$feature)))
        orthoHeight = length(levels(as.factor(orderedOrthoDf$feature)))
        
        arrangeGrob(plot_seed,plot_ortho,ncol=1, heights = c(seedHeight,orthoHeight))
      }
    }
  }
  
  output$archi_plot <- renderPlot({
    g <- archi_plot()
    grid.draw(g)
  })
  
  ######## render domain architecture plot
  output$archi_plot.ui <- renderUI({
    if (v3$doPlot3 == FALSE) {
      domainIN <- unlist(strsplit(toString(input$mainInput),","))
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
      withSpinner(plotOutput("archi_plot",height = input$archiHeight, width = input$archi_width))
    }
  })
  
  ######## download architecture plot ***** something strange with archi_plot()
  output$archi_download <- downloadHandler(
    filename = function() {c("domains.pdf")},
    content = function(file) {
      g <- archi_plot()
      grid.draw(g)
      ggsave(file, plot = g, width = input$selected_width*0.056458333, height = input$selected_height*0.056458333, units="cm", dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
  #############################################################
  ######################## GENE AGE ###########################
  #############################################################
  
  ##### gene age estimation
  gene_ageDf <- reactive({
    if (v$doPlot == FALSE) return()
    
    rankList <- c("family","class","phylum","kingdom","superkingdom","root")
    
    ##### get selected (super)taxon ID
    rank_select <- input$rank_select
    rankName = substr(rank_select,4,nchar(rank_select))
    
    taxa_list <- get_name_list(FALSE, FALSE)
    superID <- as.numeric(taxa_list$ncbiID[taxa_list$fullName == input$in_select & taxa_list$rank == rankName])
    
    ### full non-duplicated taxonomy data
    Dt <- get_taxa_list(FALSE)
    
    ### subset of taxonomy data, containing only ranks from rankList
    subDt <- Dt[,c("abbrName",rankList)]
    
    ### get (super)taxa IDs for one of representative species
    firstLine <- Dt[Dt[,rankName]==superID,][1,]  # get all taxon info for 1 representative
    subFirstLine <- firstLine[,c("abbrName",rankList)]
    
    ### compare each taxon ncbi IDs with selected taxon
    ### and create a "category" data frame
    catDf <- data.frame("ncbiID" = character(), "cat" = character(), stringsAsFactors=FALSE)
    for(i in 1:nrow(subDt)){
      cat <- subDt[i,] %in% subFirstLine
      cat[cat == FALSE] <- 0
      cat[cat == TRUE] <- 1
      cat <- paste0(cat,collapse = "")
      catDf[i,] <- c(as.character(subDt[i,]$abbrName),cat)
    }
    
    ### get main input data
    mdData <- droplevels(get_data_filtered())
    mdData <- mdData[,c("geneID","ncbiID","orthoID","var1","var2","presSpec")]
    
    ### add "category" into mdData
    mdDataExtended <- merge(mdData,catDf,by="ncbiID",all.x = TRUE)
    
    mdDataExtended$var1[mdDataExtended$var1 == "NA" | is.na(mdDataExtended$var1)] <- 0
    mdDataExtended$var2[mdDataExtended$var2 == "NA" | is.na(mdDataExtended$var2)] <- 0
    
    ### remove cat for "NA" orthologs and also for orthologs that do not fit cutoffs
    if(nrow(mdDataExtended[mdDataExtended$orthoID == "NA"| is.na(mdDataExtended$orthoID),]) > 0){
      mdDataExtended[mdDataExtended$orthoID == "NA"| is.na(mdDataExtended$orthoID),]$cat <- NA
    }
    
    mdDataExtended <- mdDataExtended[complete.cases(mdDataExtended),]
    
    # filter by %specpres, var1, var2 ..
    mdDataExtended$cat[mdDataExtended$var1 < input$var1[1]] <- NA
    mdDataExtended$cat[mdDataExtended$var1 > input$var1[2]] <- NA
    mdDataExtended$cat[mdDataExtended$var2 < input$var2[1]] <- NA
    mdDataExtended$cat[mdDataExtended$var2 > input$var2[2]] <- NA
    mdDataExtended$cat[mdDataExtended$presSpec < input$percent[1]] <- NA
    mdDataExtended$cat[mdDataExtended$presSpec > input$percent[2]] <- NA
    
    mdDataExtended <- mdDataExtended[complete.cases(mdDataExtended),]
    
    ### get the furthest common taxon with selected taxon for each gene
    gene_ageDf <- as.data.frame(tapply(mdDataExtended$cat, mdDataExtended$geneID, min))
    
    setDT(gene_ageDf, keep.rownames = TRUE)[]
    setnames(gene_ageDf, 1:2, c("geneID","cat"))  # rename columns
    row.names(gene_ageDf) <- NULL   # remove row names
    
    ### convert cat into gene_age
    gene_ageDf$age[gene_ageDf$cat == "0000001"] <- "07_LUCA"
    gene_ageDf$age[gene_ageDf$cat == "0000011" | gene_ageDf$cat == "0000010"] <- paste0("06_",as.character(taxa_list$fullName[taxa_list$ncbiID == subFirstLine$superkingdom & taxa_list$rank == "superkingdom"]))
    gene_ageDf$age[gene_ageDf$cat == "0000111"] <- paste0("05_",as.character(taxa_list$fullName[taxa_list$ncbiID == subFirstLine$kingdom & taxa_list$rank == "kingdom"]))
    gene_ageDf$age[gene_ageDf$cat == "0001111"] <- paste0("04_",as.character(taxa_list$fullName[taxa_list$ncbiID == subFirstLine$phylum & taxa_list$rank == "phylum"]))
    gene_ageDf$age[gene_ageDf$cat == "0011111"] <- paste0("03_",as.character(taxa_list$fullName[taxa_list$ncbiID == subFirstLine$class & taxa_list$rank == "class"]))
    gene_ageDf$age[gene_ageDf$cat == "0111111"] <- paste0("02_",as.character(taxa_list$fullName[taxa_list$ncbiID == subFirstLine$family & taxa_list$rank == "family"]))
    gene_ageDf$age[gene_ageDf$cat == "1111111"] <- paste0("01_",as.character(taxa_list$fullName[taxa_list$fullName == input$in_select & taxa_list$rank == rankName]))
    
    ### return gene_age data frame
    gene_ageDf <- gene_ageDf[,c("geneID","cat","age")]
    
    gene_ageDf$age[is.na(gene_ageDf$age)] <- "Undef"
    return(gene_ageDf)
  })
  
  gene_ageDfMod <- reactive({
    gene_ageDf <- gene_ageDf()
    countDf <- plyr::count(gene_ageDf,c('age'))
    countDf$percentage <- round(countDf$freq/sum(countDf$freq)*100)
    countDf$pos <- cumsum(countDf$percentage) - (0.5 * countDf$percentage)
    return(countDf)
  })
  
  gene_agePlot <- function(){
    countDf <- gene_ageDfMod()
    p <- ggplot(countDf, aes(fill=age, y=percentage, x=1)) +
      geom_bar(stat="identity") +
      scale_y_reverse() +
      coord_flip() +
      theme_minimal()
    p <- p + geom_text(data=countDf, aes(x = 1, y = 100-pos, label = paste0(freq,"\n",percentage,"%")),size=4*input$gene_age_text)
    p <- p + theme(legend.position="bottom", legend.title = element_blank(), legend.text = element_text(size=12*input$gene_age_text),
                   axis.title = element_blank(), axis.text = element_blank()) +
      scale_fill_brewer(palette="Spectral") +
      guides(fill=guide_legend(nrow=round(nrow(countDf)/3,0),byrow=TRUE))
    p
  }
  
  output$gene_agePlot <- renderPlot({
    if(input$auto_update == FALSE){
      # Add dependency on the update button (only update when button is clicked)
      input$updateBtn
      
      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        gene_agePlot()
      })
    } else {
      gene_agePlot()
    }
  })
  
  output$gene_age.ui <- renderUI({
    if(v$doPlot == FALSE){
      return()
    } else{
      ## if auto_update is NOT selected, use updateBtn to trigger plot changing
      if(input$auto_update == FALSE){
        # Add dependency on the update button (only update when button is clicked)
        input$updateBtn
        
        # Add all the filters to the data based on the user inputs
        # wrap in an isolate() so that the data won't update every time an input
        # is changed
        isolate({
          withSpinner(
            plotOutput("gene_agePlot",width=600*input$gene_age_width,height = 150*input$gene_age_height,
                       click = "plot_click_gene_age")
          )
        })
      }
      ## if auto_update is true
      else {
        withSpinner(
          plotOutput("gene_agePlot",width=600*input$gene_age_width,height = 150*input$gene_age_height,
                     click = "plot_click_gene_age")
        )
      }
    }
  })
  
  ### download gene age plot
  output$gene_age_plot_download <- downloadHandler(
    filename = function() {"gene_age_plot.pdf"},
    content = function(file) { 
      ggsave(file, plot = gene_agePlot(), width=600*input$gene_age_width*0.056458333,height = 150*input$gene_age_height*0.056458333, units="cm", dpi=300, device = "pdf")
    }
  )
  
  ##### render genAge.table based on clicked point on gene_agePlot
  selectedgene_age <- reactive({
    if(v$doPlot == FALSE){return()}
    data <- gene_ageDf()
    
    # calculate the coordinate range for each age group
    rangeDf <- plyr::count(data,c('age'))
    
    rangeDf$percentage <- round(rangeDf$freq/sum(rangeDf$freq)*100)
    rangeDf$rangeStart[1] <- 0
    rangeDf$rangeEnd[1] <- rangeDf$percentage[1]
    if(nrow(rangeDf) > 1){
      for(i in 2:nrow(rangeDf)){
        rangeDf$rangeStart[i] <- rangeDf$rangeEnd[i-1]+1
        rangeDf$rangeEnd[i] <- rangeDf$percentage[i] + rangeDf$rangeEnd[i-1]
      }
    }
    
    # get list of selected age group
    if (is.null(input$plot_click_gene_age$x)) {return()}
    else{
      corX = 100-round(-input$plot_click_gene_age$x)
      selectAge <- as.character(rangeDf[rangeDf$rangeStart <= corX & rangeDf$rangeEnd >= corX,]$age)
      subData <- subset(data, age == selectAge)
      data <- data[data$age == selectAge,]
    }
    
    # return list of genes
    geneList <- levels(as.factor(subData$geneID))
    geneList
  })
  
  output$gene_age.table <- renderTable({
    if (is.null(input$plot_click_gene_age$x)) {return()}
    
    data <- as.data.frame(selectedgene_age())
    data$number <- rownames(data)
    colnames(data) <- c("geneID","No.")
    data <- data[,c("No.","geneID")]
    data
  })
  
  ### download gene list from gene_ageTable
  output$gene_age_table_download <- downloadHandler(
    filename = function(){c("selectedGeneList.out")},
    content = function(file){
      dataOut <- selectedgene_age()
      write.table(dataOut,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ### check if add_cluster_cutom_profile (profile clustering) are being clicked
  observe({
    if(input$add_cluster_cutom_profile == TRUE | input$add_cons_gene_custom_profile == TRUE ){#| input$addGcGenesCustomProfile == TRUE
      shinyjs::disable('add_custom_profile')
    } else {
      shinyjs::enable('add_custom_profile')
    }
  })
  
  output$add_custom_profile_check.ui <- renderUI({
    if(input$add_cluster_cutom_profile == TRUE  | input$add_cons_gene_custom_profile == TRUE ){ #| input$addGcGenesCustomProfile == TRUE
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Profile clustering</strong> or <strong>Core genes finding</strong> or <strong>Groupcomparison</strong> &nbsp;to enable this function)</em></p>')
    }
  })
  
  ### reset gene_age_prot_config
  observeEvent(input$reset_gene_age_prot_config, {
    shinyjs::reset("gene_age_width")
    shinyjs::reset("gene_age_height")
    shinyjs::reset("gene_age_text")
  })
  
  #############################################################
  ##################### CORE GENES ############################
  #############################################################
  
  ### render list of available taxa
  output$taxa_list_cons.ui = renderUI({
    filein <- input$mainInput
    # if(input$demo == TRUE){
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      filein = 1
    }
    
    if(is.null(filein)){return(selectInput('inTaxa','Select taxa of interest:',"none"))}
    if(v$doPlot == FALSE){return(selectInput('inTaxa','Select taxa of interest:',"none"))}
    else{
      choice <- alltaxa_list()
      choice$fullName <- as.factor(choice$fullName)
      
      out <- as.list(levels(choice$fullName))
      out <- append("none",out)
      
      if(input$apply_cons_taxa == TRUE){
        out <- consTaxaName()
        selectInput('taxaCons','Select taxa of interest:',out,selected=out,multiple=TRUE)
      } else {
        selectInput('taxaCons','Select taxa of interest:',out,selected=out[1],multiple=TRUE)
      }
    }
  })
  
  cons_geneDf <- reactive({
    if (v$doPlot == FALSE) return()
    
    rankName = substr(input$rank_select,4,nchar(input$rank_select))
    
    ### get ID list of chosen taxa
    taxa_list <- get_name_list(FALSE, FALSE)
    
    if("none" %in%input$taxaCons){superID = NA}
    else{superID <- taxa_list$ncbiID[taxa_list$fullName %in% input$taxaCons & taxa_list$rank %in% c(rankName,"norank")]}
    
    ### get main input data
    mdData <- get_data_filtered()
    mdData <- mdData[,c("geneID","ncbiID","fullName","supertaxon","supertaxonID","rank","presSpec","mVar1","mVar2")]
    
    ### filter by selecting taxa
    if(is.na(superID[1])){data <- NULL}
    else{
      data <- subset(mdData,supertaxonID %in% superID & presSpec >= input$percent_cons)
      # get supertaxa present in each geneID
      supertaxonCount <- as.data.frame(plyr::count(data,c('geneID','supertaxonID')))
      # count number of supertaxa present in each geneID and get only gene that contains all chosen taxa
      count <- as.data.frame(table(supertaxonCount$geneID))
      cons_gene <- subset(count,Freq == length(superID))
      cons_gene$Var1 <- factor(cons_gene$Var1)
      
      return(levels(cons_gene$Var1))
    }
  })
  
  output$cons_gene.table <- renderDataTable({
    data <- cons_geneDf()
    if(is.null(data)){return()}
    else {
      data <- as.data.frame(data)
      # data$number <- rownames(data)
      # colnames(data) <- c("geneID","No.")
      # data <- data[,c("No.","geneID")]
      data
    }
  })
  
  ### download gene list from cons_gene.table
  output$cons_gene_table_download <- downloadHandler(
    filename = function(){c("consensusGeneList.out")},
    content = function(file){
      dataOut <- cons_geneDf()
      write.table(dataOut,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ### check if add_cluster_cutom_profile (profile clustering) or add_custom_profile (gene age plot) are being clicked
  observe({
    if(input$add_cluster_cutom_profile == TRUE | input$add_custom_profile == TRUE){ # | input$addGcGenesCustomProfile == TRUE
      shinyjs::disable('add_cons_gene_custom_profile')
    } else {
      shinyjs::enable('add_cons_gene_custom_profile')
    }
  })
  
  output$add_cons_gene_custom_profile_check.ui <- renderUI({
    if(input$add_cluster_cutom_profile == TRUE | input$add_custom_profile == TRUE ){ #| input$addGcGenesCustomProfile == TRUE
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Profiles clustering</strong> or <strong>Gene age estimating</strong> or r <strong>Group Comparioson</strong>&nbsp;to enable this function)</em></p>')
    }
  })
  
  ######## print list of available taxonomy ranks (the lowest rank is the same as the chosen main rank)
  output$rank_select_cons = renderUI({
    mainRank <- input$rank_select
    mainChoices = get_taxonomy_ranks()
    consChoices <- mainChoices[mainChoices >= mainRank]
    
    selectInput("rank_select_cons", label = h5("Select taxonomy rank:"),
                choices = as.list(consChoices),
                selected = mainRank)
  })
  
  ######## print list of available taxa for customized plot (based on rank from rank_select_cus)
  taxa_select_cons <- reactive({
    rank_select_cons = input$rank_select_cons
    
    if(length(rank_select_cons) == 0){return()}
    else{
      ### load list of unsorted taxa
      Dt <- get_taxa_list(TRUE)
      
      ### load list of taxon name
      nameList <- get_name_list(TRUE, FALSE)
      
      rankName = substr(rank_select_cons,4,nchar(rank_select_cons))   # get rank name from rank_select
      choice <- as.data.frame
      choice <- rbind(Dt[rankName])
      colnames(choice) <- "ncbiID"
      choice <- merge(choice,nameList,by="ncbiID",all = FALSE)
      return(choice)
    }
  })
  
  output$taxa_select_cons = renderUI({
    choice <- taxa_select_cons()
    choice$fullName <- as.factor(choice$fullName)
    selectInput('taxa_select_cons',h5('Choose (super)taxon of interest:'),as.list(levels(choice$fullName)),levels(choice$fullName)[1])
  })
  
  ######## get list of taxa based on selected taxa_select_cus
  consTaxaName <- reactive({
    
    taxa_select_cons = input$taxa_select_cons
    rankName = substr(input$rank_select_cons,4,nchar(input$rank_select_cons))
    
    if(taxa_select_cons == ""){return()}
    
    ### load list of unsorted taxa
    Dt <- get_taxa_list(TRUE)
    
    ### get ID of customized (super)taxon
    taxa_list <- get_name_list(FALSE, FALSE)
    superID <- taxa_list$ncbiID[taxa_list$fullName == taxa_select_cons & taxa_list$rank %in% c(rankName,"norank")]
    
    ### from that ID, get list of all taxa for main selected taxon
    mainRankName = substr(input$rank_select,4,nchar(input$rank_select))
    constaxa_id <- levels(as.factor(Dt[mainRankName][Dt[rankName] == superID,]))
    
    consTaxaName <- taxa_list$fullName[taxa_list$rank %in% c(mainRankName,"norank") & taxa_list$ncbiID %in% constaxa_id]
    
    return(consTaxaName)
  })
  
  
  #############################################################
  #################### CLUSTERING PROFILES ####################
  #############################################################
  
  ### cluster data
  clusterDataDend <- reactive({
    if(v$doPlot == FALSE){return()}
    # dataframe for calculate distance matrix
    dataHeat <- dataHeat()
    
    subDataHeat <- subset(dataHeat,dataHeat$presSpec > 0)
    subDataHeat <- subDataHeat[,c('geneID','supertaxon','presSpec')]
    subDataHeat <- subDataHeat[!duplicated(subDataHeat),]
    
    wideData <- spread(subDataHeat, supertaxon, presSpec)
    dat <- wideData[,2:ncol(wideData)]  # numerical columns
    rownames(dat) <- wideData[,1]
    dat[is.na(dat)] <- 0
    
    dd.col <- as.dendrogram(hclust(dist(dat, method = input$dist_method), method = input$cluster_method))
  })
  
  ### plot clustered profiles
  dendrogram <- function(dd.col){
    py <- as.ggdend(dd.col)
    p <- ggplot(py, horiz = TRUE, theme=theme_minimal()) +
      theme(axis.title = element_blank(), axis.text.y = element_blank())
    p
  }
  
  output$dendrogram <- renderPlot({
    if(v$doPlot == FALSE){return()}
    dendrogram(clusterDataDend())
  })
  
  output$cluster.ui <- renderUI({
    withSpinner(
      plotOutput("dendrogram",width=input$cluster_plot.width, height=input$cluster_plot.height,
                 brush = brushOpts(
                   id = "plot_brush",
                   delay = input$brush_delay,
                   delayType = input$brush_policy,
                   direction = input$brush_dir,
                   resetOnNew = input$brush_reset)
      )
    )
  })
  
  ### download clustered plot
  output$download_cluster <- downloadHandler(
    filename = function() {"clustered_plot.pdf"},
    content = function(file) {
      ggsave(file, plot = dendrogram(), dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
  #### render brushed_cluster.table based on clicked point on dendrogram plot
  brushed_clusterGene <- reactive({
    if(v$doPlot == FALSE){return()}
    
    dd.col <- clusterDataDend()
    dt <- dendro_data(dd.col)
    dt$labels$label <- levels(dt$labels$label)
    
    # get list of selected gene(s)
    if (is.null(input$plot_brush)) {return()}
    else{
      top = as.numeric(-round(input$plot_brush$ymin))
      bottom = as.numeric(-round(input$plot_brush$ymax))
      
      df <- dt$labels[bottom:top,]
    }
    
    # return list of genes
    df <- df[complete.cases(df),3]
  })
  
  output$brushed_cluster.table <- renderTable({
    if (is.null(input$plot_brush$ymin)) {return()}
    
    data <- as.data.frame(brushed_clusterGene())
    data$number <- rownames(data)
    colnames(data) <- c("geneID","No.")
    data <- data[,c("No.","geneID")]
    data
  })
  
  ### download gene list from brushed_cluster.table
  output$download_cluster_genes <- downloadHandler(
    filename = function(){c("selectedClusteredGeneList.out")},
    content = function(file){
      dataOut <- brushed_clusterGene()
      write.table(dataOut,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ### check if add_custom_profile (gene age plot) are being clicked
  observe({
    if(input$add_custom_profile == TRUE | input$add_cons_gene_custom_profile == TRUE ){ #| input$addGcGenesCustomProfile == TRUE
      shinyjs::disable('add_cluster_cutom_profile')
    }else{
      shinyjs::enable('add_cluster_cutom_profile')
    }
  })
  
  output$add_cluster_cutom_profile_check.ui <- renderUI({
    if(input$add_custom_profile == TRUE | input$add_cons_gene_custom_profile == TRUE ){ #| input$addGcGenesCustomProfile == TRUE
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Gene age estimation</strong> or <strong>Core genes finding</strong> or <strong>Group Comparison</strong> &nbsp;to enable this function)</em></p>')
    }
  })
  
  #############################################################
  ############### FILTERED DATA FOR DOWNLOADING ###############
  #############################################################
  
  ################### FOR MAIN PROFILE ########################
  
  ######## render variable used for identifying representative genes
  output$ref_var_main.ui <- renderUI({
    if(nchar(input$var2_id) < 1 & nchar(input$var1_id) < 1){
      radioButtons(inputId="ref_var_main", label="Reference variable", choices=list(input$var1_id,input$var2_id), selected=input$var1_id)
    } else if(nchar(input$var2_id) < 1){
      radioButtons(inputId="ref_var_main", label="Reference variable", choices=list(input$var1_id), selected=input$var1_id)
    } else {
      radioButtons(inputId="ref_var_main", label="Reference variable", choices=list(input$var1_id,input$var2_id), selected=input$var1_id)
    }
  })
  
  ######## filtered data for downloading
  download_data <- reactive({
    ### check input
    if (v$doPlot == FALSE) return()
    
    ### filtered data
    dataOut <- get_data_filtered()
    dataOut <- as.data.frame(dataOut[dataOut$presSpec > 0,])
    dataOut <- dataOut[!is.na(dataOut$geneID),]
    
    dataOut <- as.data.frame(dataOut[dataOut$presSpec >= input$percent[1],])
    dataOut <- as.data.frame(dataOut[dataOut$var1 >= input$var1[1] & dataOut$var1 <= input$var1[2],])
    
    if(!all(is.na(dataOut$var2))){
      dataOut <- as.data.frame(dataOut[dataOut$var2 >= input$var2[1] & dataOut$var2 <= input$var2[2],])
    } else {
      dataOut$var2 <- 0
    }
    
    ### select only representative genes if chosen
    if(input$get_representative_main == TRUE){
      if(is.null(input$ref_var_main)){return()}
      else{
        if(input$ref_var_main == input$var1_id){
          dataOutAgg <- aggregate(as.numeric(dataOut$var1), by = list(dataOut$geneID,dataOut$ncbiID), FUN = input$ref_type_main)
        } else if (input$ref_var_main == input$var2_id){
          dataOutAgg <- aggregate(as.numeric(dataOut$var2), by = list(dataOut$geneID,dataOut$ncbiID), FUN = input$ref_type_main)
        } else {
          dataOutAgg <- dataOut[dataOut,c("geneID","ncbiID","var1")]
        }
        colnames(dataOutAgg) <- c("geneID","ncbiID","var_best")
        
        dataOutRepresentative <- merge(dataOut,dataOutAgg,by=c("geneID","ncbiID"),all.x=TRUE)
        
        if(input$ref_var_main == input$var1_id){
          dataOut <- dataOutRepresentative[dataOutRepresentative$var1 == dataOutRepresentative$var_best,]
        } else if (input$ref_var_main == input$var2_id){
          dataOut <- dataOutRepresentative[dataOutRepresentative$var2 == dataOutRepresentative$var_best,]
        } else {
          dataOut <- dataOut
        }
        
        dataOut$dup <- paste0(dataOut$geneID,"#",dataOut$ncbiID)  # <- used to select only one ortholog, if there exist more than one "representative"
        dataOut <- dataOut[!duplicated(c(dataOut$dup)),]
      }
    }
    
    ### sub select columns of dataout
    dataOut <- dataOut[,c("geneID","orthoID","fullName","ncbiID","supertaxon","var1","var2","presSpec")] #,"numberSpec"
    dataOut <- dataOut[order(dataOut$geneID,dataOut$supertaxon),]
    dataOut <- dataOut[complete.cases(dataOut),]
    
    dataOut$geneID <- as.character(dataOut$geneID)
    dataOut$fullName <- as.character(dataOut$fullName)
    dataOut$ncbiID <- substr(dataOut$ncbiID,5,nchar(as.character(dataOut$ncbiID)))
    dataOut$supertaxon <- substr(dataOut$supertaxon,6,nchar(as.character(dataOut$supertaxon)))
    dataOut$var1 <- as.character(dataOut$var1)
    dataOut$var2 <- as.character(dataOut$var2)
    # dataOut$numberSpec <- as.numeric(dataOut$numberSpec)
    dataOut$presSpec <- as.numeric(dataOut$presSpec)
    
    ### rename columns
    names(dataOut)[names(dataOut)=="presSpec"] <- "%Spec"
    # names(dataOut)[names(dataOut)=="numberSpec"] <- "totalSpec"
    if(nchar(input$var1_id) > 0){
      names(dataOut)[names(dataOut)=="var1"] <- input$var1_id
    } else {
      dataOut <- subset(dataOut, select = -c(var1) )
    }
    if(nchar(input$var2_id) > 0){
      names(dataOut)[names(dataOut)=="var2"] <- input$var2_id
    } else {
      dataOut <- subset(dataOut, select = -c(var2) )
    }
    
    # return data for downloading
    dataOut <- as.matrix(dataOut)
    return(dataOut)
  })
  
  ######## download data
  output$download_data <- downloadHandler(
    filename = function(){c("filteredData.out")},
    content = function(file){
      dataOut <- download_data()
      write.table(dataOut,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ######## data table ui tab
  output$filtered_main_data <- renderDataTable(rownames = FALSE,{
    if(v$doPlot == FALSE){return()}
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
  
  ######## download FASTA
  output$download_fasta.ui <- renderUI({
    # if(input$demo_data == "demo"){
    #   HTML("<p><span style=\"color: #ff0000;\"><em>Be patient! For large number of taxa this can take up to 3 minutes!</em></span></p>")
    # }
  })
  
  output$download_fasta <- downloadHandler(
    filename = function(){c("filteredSeq.fa")},
    content = function(file){
      fastaOutDf <- fastaOutData(as.data.frame(download_data()))
      write.table(fastaOutDf,file,sep="\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
    }
  )
  
  ################### FOR CUSTOMIZED PROFILE ########################
  
  ######## render variable used for identifying representative genes
  output$representative_info.ui <- renderUI({
    msg <- paste0("NOTE: According to your choice in [Download filtered data -> Main data], only representative sequences with ",as.character(input$ref_type_main)," ",as.character(input$ref_var_main),"  will be downloaded!")
    strong(em(msg),style = "color:red")
  })
  
  ######## filtered data for downloading
  download_custom_data <- reactive({
    ### check input
    if (v$doPlot == FALSE) return()
    
    data <- as.data.frame(download_data())
    
    ### get subset of data according to selected genes/taxa
    if(!is.null(input$inSeq) | !is.null(input$inTaxa)){
      if(input$inSeq[1] != "all" & input$inTaxa[1] == "all"){
        customData <- subset(data,geneID %in% input$inSeq) ##### <=== select data for selected sequences only
      } else if(input$inSeq[1] == "all" & input$inTaxa[1] != "all"){
        customData <- subset(data,supertaxon %in% input$inTaxa) ##### <=== select data for selected taxa only
      } else if(input$inSeq[1] != "all" & input$inTaxa[1] != "all") {
        customData <- subset(data,geneID %in% input$inSeq & supertaxon %in% input$inTaxa) ##### <=== select data for selected sequences and taxa
      } else {
        customData <- data
      }
    } else {
      customData <- data
    }
    ### return data
    customData <- as.matrix(customData)
    customData
  })
  
  ######## download data
  output$download_custom_data <- downloadHandler(
    filename = function(){c("customFilteredData.out")},
    content = function(file){
      dataOut <- download_custom_data()
      write.table(dataOut,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ######## data table ui tab
  output$filtered_custom_data <- renderDataTable(rownames = FALSE,{
    if(v$doPlot == FALSE){return()}
    data <- download_custom_data()
    data
  })
  
  ######## download FASTA
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
      fastaOutDf <- fastaOutData(as.data.frame(download_custom_data()))
      write.table(fastaOutDf,
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
  
  # Select in_group ===========================================================
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
        choice <- taxa_select_gc(input$rank_select_gc)
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
        dt <- get_taxa_list(FALSE) # get the taxa 
        

        # get the info for the reference protein from the namelist
        reference <- subset(name_list,
                            name_list$fullName == input$in_select)
        
        # get the id for every rank for the reference protein
        rank_name <- substr(input$rank_select, 4, nchar(input$rank_select))
        reference_dt <- dt[dt[, rank_name] == reference$ncbiID, ]
        
        # save the next higher rank 
        reference_higher_rank <- reference_dt[higher_rank_name]
        reference_higher_rank <- {
          reference_higher_rank[!duplicated(reference_higher_rank),]
        }


        # get all the taxa with the same id in the next higher rank
        selected_taxa_dt <- {
          subset(dt, dt[,higher_rank_name] %in% reference_higher_rank)
        }
        selected_taxa_dt <- {
          selected_taxa_dt[!duplicated(selected_taxa_dt[rank_name]),]
        }
        
        # get a list with all the ids with reference_higher_rank as parent
        selected_taxa_ids <- selected_taxa_dt[rank_name]
        if(length(selected_taxa_ids[[1]]) >= 1){
          selected_taxa_ids <- selected_taxa_ids[[1]]
        }

        selected_taxa <- subset(name_list, name_list$rank == rank_name)
        selected_taxa <- {
          subset(selected_taxa, selected_taxa$ncbiID %in% selected_taxa_ids)
        }
        
        default_select <- selected_taxa$fullName 
        
        choice <- taxa_select_gc(input$rank_select)
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
  
  # Functions: Popup Window Select Rank =======================================
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
  
  # print list of available taxa
  taxa_select_gc <- function(rank_select_gc){
    
    # if there is no rank set, there can not be any available taxa
    if (length(rank_select_gc) == 0) return()
    else{
      
      # load list of unsorted taxa
      dt <- get_taxa_list(TRUE)
 
      
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
  
  # Supertaxon of intrest in the popup window for the rank
  output$taxa_select_gc <- renderUI({
    choice <- taxa_select_gc(input$rank_select_gc)
    choice$fullName <- as.factor(choice$fullName)
    selectInput("taxa_select_gc", h5("Choose (super)taxon of interest:"),
                as.list(levels(choice$fullName)),
                levels(choice$fullName)[1])
  })
  
  # get list of taxa based on selected taxa_select_gc
  taxa_name_gc <- reactive({
    
    taxa_select_gc <- input$taxa_select_gc
    rank_name <- substr(input$rank_select_gc, 4, nchar(input$rank_select_gc))
    
    if (taxa_select_gc == "") return()
    
    # load list of unsorted taxa
    dt <- get_taxa_list(TRUE)
    
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
  
  # Select Gene ===============================================================
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
  
  # Buttons to choose the variable ============================================
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
  
  # Slider to choose significance =============================================
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
  
  # List with all significant Genes ===========================================
  output$get_significant_genes <- renderUI({
    input$plot_gc
    
    isolate({
      get_significant_genes()
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
  
  # Generate output plots =====================================================
  output$plots_gc <- renderUI(get_plots_gc())
  
  # Select Feaures you want to see in the barplots (default: All)
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
  
  # Select Plots to download ==================================================
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
      }})
  })
  
  
  # Generate the  main Output =================================================
  # Deciding which plots should be shown
  get_plots_gc <- reactive({
    gene <- as.character(input$selected_gene_gc)
    input$plot_gc
    if (is.null(significant_genes_gc))return()
    else if (gene == "all"){
      get_plot_output_list(significant_genes_gc)
    }else{
      x <- {
        significant_genes_gc[significant_genes_gc$geneID == gene, ]
      }
      if (nrow(x) == 0) return()
      get_plot_output_list(x)
    }
  })
  
  # Generate the list with all plots
  get_plot_output_list <- function(genes) {
    # if we dont have parameters we can not generate plots
    if (is.null(genes)) return()
    
    # Insert plot output objects the list
    plot_output_list <- lapply(1:nrow(genes), function(i) {
      plotname <- paste(genes[i, 1])
      plot_output_object <- renderPlot(get_multiplot(genes[i, ]),
                                       height = 650, width = 700)
    })
    do.call(tagList, plot_output_list) # needed to display properly.
    return(plot_output_list)
  }
  
  
  # Generating the plots ======================================================
  # Put the plots for one spicific gene in one multiplot
  get_multiplot <- function(gene_info){
    
    # Sorting the information to the selected gene
    gene <- as.character(gene_info$geneID)
    in_group <- as.data.frame(gene_info$in_group)
    out_group <- as.data.frame(gene_info$out_group)
    features <- as.data.frame(gene_info$features)
    pvalues <- gene_info$pvalues
    var <- gene_info$var
    
    # the barplot does not depent on the selected variable(s)
    barplot <-  get_barplot_gc(gene, in_group, out_group, features)

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
                                 info_p1)
      boxplot2 <- get_boxplot_gc(in_group,
                                 out_group,
                                 input$var2_id,
                                 gene,
                                 c2,
                                 info_p2)
      
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
      
      boxplot <- get_boxplot_gc(in_group, out_group, var, gene, "grey", info)
      
      m <- grid.arrange(textGrob(gene),
                        boxplot,
                        barplot,
                        heights = c(0.02, 0.45, 0.458), ncol = 1)
    }
    return(m)
  }
  
  # Create a Boxplot
  get_boxplot_gc <- function (in_group_df,
                              out_group_df,
                              var,
                              gene,
                              colour,
                              info){
    
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
    
    names <- c(paste("in_group \n n=", a, sep = ""),
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
  
  # Create Barplot
  get_barplot_gc <- function(selected_gene, in_group, out_group, features){
    
    subdomain_df <- features
 
    # only show features that interest the user
    if (!("all" %in% input$interesting_features)){
      ifeatures <- NULL
      for (x in input$interesting_features){
        f <- subset(subdomain_df$feature, startsWith(subdomain_df$feature, x))
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
    
    if (nrow(in_group_domain_df) == 0) data_in <- NULL
    else{
      feature <- unique(in_group_domain_df$feature)
      data_in <- as.data.frame(feature)
      data_in$amount <- 0
    }
    
    if (nrow(out_group_domain_df) == 0)data_out <- NULL
    else{
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
    
    if (is.null(data_in) & !is.null(data_out))data_barplot <- data_out
    else if (is.null(data_out) & !is.null(data_in))data_barplot <- data_in
    else if (!is.null(data_in) & !is.null(data_out)){
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
  
  # Get the p_values to print under the plot
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
  
  # Get the list with all the significant genes and the dataset ===============
  get_significant_genes <- function (){
    if (is.null(input$selected_in_group_gc)
        | length(input$list_selected_genes_gc) == 0)return()
    
    # load name List
    name_list <- get_name_list(TRUE, TRUE)
    
    # load list of unsorted taxa
    dt <- get_taxa_list(FALSE)
    
    # Parameters that are identical for all genes
    in_group <-  input$selected_in_group_gc
    rank <- input$rank_select
    var <- input$var_name_gc
    
    # Updateing of the Input ==================================================
    # if there is more than one element in the in_group
    # we look at the next common anchstor
    ancestor <- get_common_ancestor(in_group, rank, name_list, dt)
    if (is.null(ancestor))return()
    in_group <- ancestor[1]
    rank <- ancestor[2]
    
    # List with all significant genes
    significantGenes <- data.frame(
      geneID = character(),
      in_group = I(list()),
      out_group = I(list()),
      pvalues = I(list()),
      features = I(list()),
      databases = I(list()))
    
    # List of genes to look at
    data_full <- get_data_filtered()
    if (is.element("all", input$list_selected_genes_gc)){
      genes <- data_full$geneID
      genes <- genes[!duplicated(genes)]
    } else {
      genes <- input$list_selected_genes_gc
    }
    genes <- sort(genes)

    # Subset depending on the rank and the in_group
    selected_subset <- get_selected_subset(rank, in_group, name_list, dt)
    selected_subset <- subset(selected_subset,
                              !selected_subset$fullName == input$in_select)
    
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
        subset(out_group_df, !out_group_df$fullName == input$in_select)
      }
      
      # Generate the p_values for the gene
      pvalues <- get_significant(in_group_df, out_group_df, var, gene)
      
      if (!is.null(pvalues)){
        new_row <- data.frame(geneID = gene,
                              in_group = NA,
                              out_group = NA,
                              pvalues = NA,
                              features = NA)
        new_row$in_group <- list(in_group_df)
        new_row$out_group <- list(out_group_df)
        new_row$pvalues <- list(pvalues)
        features  <- get_features(gene)
        new_row$features <- list(features)
        if (input$right_format_features){
          new_row$databases <- list(get_prefix_features(features))
        }
        significantGenes <- rbind(significantGenes, new_row)
      }
    }
    if (nrow(significantGenes) != 0){
      significantGenes$var <- var
      significantGenes$rank <- rank
      significant_genes_gc <<- significantGenes
    }
  }
  
  # Get the database for each feature in a specific gene
  # f is dataframe in $features
  get_prefix_features <- function(f){
    features <- f$feature
    choices <- gsub("_.*", "", features)
    choices <- choices[!duplicated(choices)]
    choices
  }
  
  # Get the Subset depending on the choosen rank
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
  
  # Generate the in_group
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
  
  # Decide if the gene is significant =========================================
  # if the gene is significant return the pvalues
  get_significant <- function(in_g, out_g, var, gene){
    significance_level <- input$significance
    
    if (var == "Both"){
      var1 <- input$var1
      var2 <- input$var2
      
      significant <-  FALSE
      # get the pValues for both variables
      pvalues1 <- get_p_values(in_g$var1, out_g$var1)
      pvalues2 <- get_p_values(in_g$var2, out_g$var2)
      
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
        pvalues <- get_p_values(in_g$var1, out_g$var1)
      }
      else {
        pvalues <- get_p_values(in_g$var2, out_g$var2)
      }
      
      # Analog to getting the significance with both variables
      if (is.null(pvalues)) return(NULL)
      else if (pvalues[length(pvalues)] < significance_level){
        pvalues <-  list(pvalues)
        return(pvalues)
      } else{
        return(NULL)
      }
    }
  }
  
  # calculate the p_values
  get_p_values <- function(var_in, var_out){
    # upper limit for the probability to reject H0 if it is correct
    significance_level <- input$significance
    
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
  
  # get the list with all the features in the gene
  get_features <- function(selected_gene){
    # parse domain file
    file_domain <- get_domain_file_gc(selected_gene)
    
    if (input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      domain_df <- as.data.frame(read.csv(file_domain,
                                          sep = "\t",
                                          header = F,
                                          comment.char = "",
                                          stringsAsFactors = FALSE,
                                          quote = ""))
    } else {
      if (file_domain != FALSE){
        domain_df <- as.data.frame(read.table(file_domain,
                                              sep = "\t",
                                              header = FALSE,
                                              comment.char = ""))
      }
    }
    
    if(ncol(domain_df) == 5){
      colnames(domain_df) <- c("seedID","orthoID","feature","start","end")
    } 
    else if(ncol(domain_df) == 6){
      colnames(domainDf) <- c("seedID","orthoID","feature","start","end","weight")
    } 
    else if(ncol(domain_df) == 7){
      colnames(domain_df) <- c("seedID","orthoID","feature","start","end","weight","path")
    }
    
    subdomain_df <- {
      subset(domain_df,
             substr(domain_df$seedID,
                    1,
                    nchar(as.character(selected_gene))) == selected_gene)
    }
    subdomain_df <- subdomain_df[!duplicated(subdomain_df), ]
    subdomain_df
  }
  
  # get the data where to find the features
  get_domain_file_gc <- function(group){
    # domain file
    if (input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      if (is.null(info)){
        file_domain <- "noSelectHit"
        updateButton(session, "do_domain_plot", disabled = TRUE)
      } else {
        updateButton(session, "do_domain_plot", disabled = FALSE)
        if (input$demo_data == "demo"){
          file_domain <- {
            suppressWarnings(paste0("https://github.com/BIONF/phyloprofile-data/blob/master/demo/domain_files/", group, ".domains?raw=true"))
          }
        } else {
          file_domain <- {
            suppressWarnings(paste0("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.domains_F"))
          }
        }
      }
    }else {
      if (input$anno_choose == "from file"){
        file_domain <- input$file_domain_input
        if (is.null(file_domain)){
          file_domain <- "noFileInput"
        } else {
          if (is.null(info)){
            file_domain <- "noSelectHit"
            updateButton(session, "do_domain_plot", disabled = TRUE)
          } else {
            updateButton(session, "do_domain_plot", disabled = FALSE)
            file_domain <- file_domain$datapath
          }
        }
      } else {
        if (is.null(info)){
          file_domain <- "noSelectHit"
          updateButton(session, "do_domain_plot", disabled = TRUE)
        } else {
          ### check file extension
          all_extension <- c("txt", "csv", "list", "domains", "architecture")
          flag <- 0
          for (i in 1:length(all_extension)){
            file_domain <- paste0(input$domain_path,
                                  "/",
                                  group,
                                  ".",
                                  all_extension[i])
            if (file.exists(file_domain) == TRUE){
              updateButton(session,
                           "do_domain_plot",
                           disabled = FALSE)
              flag <- 1
              break ()
            }
          }
          
          if (flag == 0){
            file_domain <- "noFileInFolder"
            updateButton(session,
                         "do_domain_plot",
                         disabled = TRUE)
          }
        }
      }
    }
    return (file_domain)
  }
  
  
  # Downloads for GroupCompairison ============================================
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
        get_multiplot_download_gc(gene)
        dev.off()
      }
      zip(zipfile = file, files = fs)
    },
    contentType = "application/zip"
  )
  
  get_multiplot_download_gc <- function(gene){
    x <- subset(significant_genes_gc,
                significant_genes_gc$geneID == gene)
    get_multiplot(x)
  }
  
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
  
  # add_custom_profile ========================================================
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
  # Essential functions =======================================================
  # ===========================================================================
  
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
  
  get_name_list <- function (as_character, delete_duplicated){
    name_list <- as.data.frame(read.table("data/taxonNamesReduced.txt",
                                          sep = "\t",
                                          header = T,
                                          fill = TRUE))
    if(as_character) {
      name_list$fullName <- as.character(name_list$fullName)
      name_list$rank <- as.character(name_list$rank)
    }
    if (delete_duplicated){
      name_list <- name_list[!duplicated(name_list), ]
    }
    
    return(name_list)
  }
  
  get_taxa_list <- function(subset_taxa_check){
    dt <- as.data.frame(read.table("data/taxonomyMatrix.txt",
                                   sep = "\t",
                                   header = T,
                                   stringsAsFactors = T))
    if(subset_taxa_check){
      dt <- dt[dt$abbrName  %in% subset_taxa(),]
    }
    return(dt)
  }
  
})
