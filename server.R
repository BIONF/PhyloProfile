############### NOT WORK WITH shinyapp.io #############
# if(!("pacman" %in% installed.packages())) install.packages("pacman")
# library(pacman)
# p_load(shiny,shinyBS,ggplot2,reshape2,plyr,dplyr,tidyr,scales,grid,gridExtra,ape,stringr,gtable,dendextend,ggdendro,gplots,data.table,taxize,install=T)
#######################################################
# packages <- c("shiny","shinyBS","ggplot2","reshape2","plyr","dplyr","tidyr","scales","grid","gridExtra","ape","stringr","gtable","dendextend","ggdendro","gplots","data.table","taxize","Biostrings","zoo","RCurl")
# sapply(packages, require, character.only = TRUE)
#######################################################

if (!require("shiny")) {install.packages("shiny")}
if (!require("shinyBS")) {install.packages("shinyBS")}
if (!require("ggplot2")) {install.packages("ggplot2")}
if (!require("reshape2")) {install.packages("reshape2")}
if (!require("plyr")) {install.packages("plyr")}
if (!require("dplyr")) {install.packages("dplyr")}
if (!require("tidyr")) {install.packages("tidyr")}
if (!require("scales")) {install.packages("scales")}
if (!require("grid")) {install.packages("grid")}
if (!require("gridExtra")) {install.packages("gridExtra")}
if (!require("ape")) {install.packages("ape")}
if (!require("stringr")) {install.packages("stringr")}
if (!require("gtable")) {install.packages("gtable")}
if (!require("dendextend")) {install.packages("dendextend")}
if (!require("ggdendro")) {install.packages("ggdendro")}
if (!require("gplots")) {install.packages("gplots")}
if (!require("data.table")) {install.packages("data.table")}
if (!require("Biostrings")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("Biostrings")
}
if (!require("taxize")) {install.packages("taxize")}
if (!require("zoo")) {install.packages("zoo")}
if (!require("RCurl")) {install.packages("RCurl")}
if (!require("shinycssloaders")) {
  if("devtools" %in% installed.packages() == FALSE){
    install.packages("devtools")
  }
  devtools::install_github('andrewsali/shinycssloaders', force = TRUE)
}

source("scripts/taxonomyProcessing.R")
source("scripts/functions.R")

############################ MAIN ############################
options(shiny.maxRequestSize=99*1024^2)  ## size limit for input 99mb

shinyServer(function(input, output, session) {
  # session$onSessionEnded(stopApp) ### Automatically stop a Shiny app when closing the browser tab
  session$allowReconnect(TRUE)
  #############################################################
  ####################  PRE-PROCESSING  #######################
  #############################################################
  
  ####### check for internet connection  ####### 
  observe({
    if(hasInternet() == FALSE){
      # toggleState("demo")
      toggleState("demo_data")
    }
  })
  output$noInternetMsg <- renderUI({
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
  taxaID <- reactive({
    if(input$idSearch > 0){
      taxain <- input$taxaList
      if(is.null(taxain)){return()}
      
      taxaNameDf <- as.data.frame(read.table(file=taxain$datapath, sep='\t',header=F,check.names=FALSE,comment.char=""))
      
      idDf <- data.frame("name"=character(),"newName"=character(),"id"=character(),"type"=character(),stringsAsFactors=FALSE)
      
      withProgress(message = 'Retrieving IDs...', value = 0,{
        for(i in 1:nrow(taxaNameDf)){
          id <- get_uid(sciname = taxaNameDf[i,])[1]
          if(is.na(id)){
            temp <- gnr_resolve(names = as.character(taxaNameDf[i,]))
            if(nrow(temp) > 0){
              newID <- get_uid(sciname = temp[1,3])[1]
              if(is.na(newID)){
                idDf[i,] <- c(as.character(taxaNameDf[i,]),as.character(temp[1,3]),paste0("NA"),"notfound")
              } else {
                idDf[i,] <- c(as.character(taxaNameDf[i,]),as.character(temp[1,3]),paste0("ncbi",newID),"notfound")
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
  output$notfoundTaxa <- renderDataTable(option = list(searching = FALSE),{
    if(input$idSearch > 0){
      if(length(taxaID())>0){
        tb <- as.data.frame(taxaID())
        tbFiltered <- tb[tb$type == "notfound",]
        notFoundDt <- tbFiltered[,c("name","newName","id")]
        colnames(notFoundDt) <- c("Summitted name","Alternative name","Alternative ID")
        notFoundDt
      }
    }
  })
  
  output$downloadNotFoundTaxa <- downloadHandler(
    filename = function(){c("mismatchedTaxa.txt")},
    content = function(file){
      tb <- as.data.frame(taxaID())
      tbFiltered <- tb[tb$type == "notfound",]
      notFoundDt <- tbFiltered[,c("name","newName","id")]
      colnames(notFoundDt) <- c("Summitted name","Alternative name","Alternative ID")
      
      write.table(notFoundDt,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ### output retrieved taxa IDs
  output$taxaID <- renderDataTable(option = list(searching = FALSE),{
    if(input$idSearch > 0){
      if(length(taxaID())>0){
        tb <- as.data.frame(taxaID())
        tbFiltered <- tb[tb$type == "retrieved",]
        retrievedDt <- tbFiltered[,c("name","id")]
        colnames(retrievedDt) <- c("Taxon_name","Taxon_ID")
        retrievedDt
      }
    }
  })
  
  output$downloadTaxaID <- downloadHandler(
    filename = function(){c("retrievedTaxaID.txt")},
    content = function(file){
      tb <- as.data.frame(taxaID())
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
  output$var1Filter.ui <- renderUI({
    if(is.null(input$var1_id)){return()}
    if(input$var1_id == ""){
      sliderInput("var1cus",paste(input$var1_id,"cutoff:"), min = 1, max = 1, step = 0.025, value = c(1.0,1.0), width = 200)
    } else {
      sliderInput("var1cus",paste(input$var1_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var1[1],input$var1[2]), width = 200)
    }
  })
  
  output$var2Filter.ui <- renderUI({
    if(is.null(input$var2_id)){return()}
    if(input$var2_id == ""){
      sliderInput("var2cus",paste(input$var2_id,"cutoff:"), min = 1, max = 1, step = 0.025, value = c(1.0,1.0), width = 200)
    } else {
      sliderInput("var2cus",paste(input$var2_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var2[1],input$var2[2]), width = 200)
    }
  })
  
  output$percentFilter.ui <- renderUI({
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
  observeEvent(input$resetMain, {
    shinyjs::reset("var1")
    shinyjs::reset("var2")
    shinyjs::reset("percent")
  })
  
  ######## reset config of main plot
  observeEvent(input$resetMainConfig, {
    shinyjs::reset("xSize")
    shinyjs::reset("ySize")
    shinyjs::reset("legendSize")
    shinyjs::reset("xAngle")
    shinyjs::reset("dotZoom")
  })
  
  ######## close main config
  observeEvent(input$applyMainConfig, {
    toggleModal(session, "mainPlotConfigBs", toggle = "close")
  })
  
  ######## reset cutoffs of Customized plot
  observeEvent(input$resetSelected, {
    shinyjs::reset("var1")
    shinyjs::reset("var2")
    shinyjs::reset("percent")
  })
  
  ######## reset config of customized plot
  observeEvent(input$resetSelectedConfig, {
    shinyjs::reset("xSizeSelect")
    shinyjs::reset("ySizeSelect")
    shinyjs::reset("legendSizeSelect")
    shinyjs::reset("xAngleSelect")
    shinyjs::reset("dotZoomSelect")
  })
  
  ######## close customized config
  observeEvent(input$applySelectedConfig, {
    toggleModal(session, "selectedPlotConfigBs", toggle = "close")
  })
  
  ######## reset colors
  observeEvent(input$defaultColorVar2, {
    shinyjs::reset("lowColor_var2")
    shinyjs::reset("highColor_var2")
  })
  
  observeEvent(input$defaultColorVar1, {
    shinyjs::reset("lowColor_var1")
    shinyjs::reset("highColor_var1")
  })
  
  observeEvent(input$defaultColorPara, {
    shinyjs::reset("paraColor")
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
          inputTaxa <- subsetTaxa()
          
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
  
  ######## render inputCheck.ui
  output$inputCheck.ui <- renderUI({
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
  output$mainInputFile.ui <- renderUI({
    # if(input$demo == TRUE){
    if(input$demo_data == "demo"){
      strong(a("Download demo input file", href="https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.long", target="_blank"))
    } else if(input$demo_data == "ampk-tor"){
      strong(a("Download demo input file", href="https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile", target="_blank"))
    } else {
      fileInput("mainInput",h5("Upload input file:"))
    }
  })
  
  output$domainInputFile.ui <- renderUI({
    # if(input$demo == TRUE){
    if(input$demo_data == "demo"){
      strong(a("Download demo domain files", href="https://github.com/BIONF/phyloprofile-data/tree/master/demo/domain_files", target="_blank"))
    } else if(input$demo_data == "ampk-tor"){
      strong(a("Download demo domain file", href="https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.domains_F", target="_blank"))
    } else {
      if(input$annoChoose == "from file"){
        fileInput("fileDomainInput","")
      } else {
        textInput("domainPath","","")
      }
    }
  })
  
  output$downloadFastaDemo.ui <- renderUI({
    if(input$demo_data == "demo"){
      strong(a("Download demo fasta file", href="https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/fasta_file/concatenatedSeq.fa", target="_blank"))
    } else if(input$demo_data == "ampk-tor"){
      strong(a("Download demo fasta file", href="https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.extended.fa", target="_blank"))
    }
  })
  
  
  ####### render description for demo data #######
  output$demoDataDescribe <- renderUI({
    if(input$demo_data == "none"){
      return()
    } else if(input$demo_data == "ampk-tor"){
      em(a("Data description", href="https://github.com/BIONF/phyloprofile-data/blob/master/expTestData/ampk-tor/README.md", target="_blank"))
    } else {
      em(a("Data description", href="https://github.com/BIONF/phyloprofile-data/blob/master/demo/README.md", target="_blank"))
    }
  })
  
  ######## get input taxa
  subsetTaxa <- reactive({
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
  output$unkTaxaStatus <- reactive({
    unkTaxa <- unkTaxa()
    length(unkTaxa) > 0
  })
  outputOptions(output, "unkTaxaStatus", suspendWhenHidden = FALSE)
  
  ### show full list of unkTaxa
  output$unkTaxaFull <- renderDataTable(option = list(searching = FALSE,pageLength = 10),{
    if(length(unkTaxa())>0){
      tb <- as.data.frame(unkTaxa())
      names(tb)[1] <- "New taxon"
      tb
    }
  })
  
  ######## check if data is loaded and "parse" button (get info from input) is clicked and confirmed
  v1 <- reactiveValues(parse = FALSE)
  observeEvent(input$BUTparse, {
    toggleModal(session, "parseConfirm", toggle = "close")
    v1$parse <- input$BUTparse
    updateButton(session, "BUTparse", disabled = TRUE)
    toggleState("newTaxaAsk")
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
            newTaxaFromFile <- fread("data/newTaxa.txt")
            
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
            newIDList <- rbind.fill(oldIDList,idList)
            newRankList <- rbind.fill(oldRankList,rankList)
            colnames(reducedInfoList) <- c("ncbiID","fullName","rank","parentID")
            newNameList <- rbind.fill(oldNameList,reducedInfoList)
            
            write.table(newIDList[!duplicated(newIDList),], file ="data/idList.txt", col.names = F, row.names = F, quote = F, sep="\t")
            write.table(newRankList[!duplicated(newRankList),], file ="data/rankList.txt", col.names = F, row.names = F, quote = F, sep="\t")
            write.table(newNameList[!duplicated(newNameList),], file ="data/taxonNamesReduced.txt", col.names = T, row.names = F, quote = F, sep="\t")
            
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
  
  output$endParsingMsg <- renderUI({
    if(nrow(invalidID$df) < 1) {
      strong(h4("PLEASE RELOAD THIS TOOL AFTER ADDING NEW TAXA!!!"),style = "color:red")
    }else{
      HTML('<p><strong><span style="color: #e12525;">SOME INVALID TAXON IDs HAVE BEEN FOUND!!</span><br>Please check the validity of the following IDs in 
           <a target="_blank" href="https://www.ncbi.nlm.nih.gov/taxonomy">NCBI taxonomy database</a>!</strong></p>')
    }
    })
  
  ######## list of taxonomy ranks for plotting
  output$rankSelect = renderUI({
    # if(input$demo == TRUE){
    if(input$demo_data == "demo"){
      selectInput("rankSelect", label = "",
                  choices = list("Strain"="05_strain","Species" = "06_species","Genus" = "10_genus", "Family" = "14_family", "Order" = "19_order", "Class" = "23_class",
                                 "Phylum" = "26_phylum", "Kingdom" = "28_kingdom", "Superkingdom" = "29_superkingdom","unselected"=""),
                  selected = "26_phylum")
    } else if(input$demo_data == "ampk-tor"){
      selectInput("rankSelect", label = "",
                  choices = list("Strain"="05_strain","Species" = "06_species","Genus" = "10_genus", "Family" = "14_family", "Order" = "19_order", "Class" = "23_class",
                                 "Phylum" = "26_phylum", "Kingdom" = "28_kingdom", "Superkingdom" = "29_superkingdom","unselected"=""),
                  selected = "06_species")
    } else {
      selectInput("rankSelect", label = "",
                  choices = list("Strain"="05_strain","Species" = "06_species","Genus" = "10_genus", "Family" = "14_family", "Order" = "19_order", "Class" = "23_class",
                                 "Phylum" = "26_phylum", "Kingdom" = "28_kingdom", "Superkingdom" = "29_superkingdom","unselected"=""),
                  selected = "06_species")
    }
  })
  
  ####### GET list of all (super)taxa
  allTaxaList <- reactive({
    
    filein <- input$mainInput
    # if(is.null(filein) & input$demo == FALSE){return()}
    if(is.null(filein) & input$demo_data == "none"){return()}
    
    rankSelect = input$rankSelect
    
    if(rankSelect == ""){return()}
    
    if(length(unkTaxa()) > 0){return()}
    
    ### load list of unsorted taxa
    Dt <- as.data.frame(read.table("data/taxonomyMatrix.txt", sep='\t', header=T, stringsAsFactors=T))
    Dt <- Dt[Dt$abbrName  %in% subsetTaxa(),]
    
    ### load list of taxon name
    nameList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T,fill = TRUE))
    nameList$fullName <- as.character(nameList$fullName)
    nameList$rank <- as.character(nameList$rank)
    
    rankName = substr(rankSelect,4,nchar(rankSelect))   # get rank name from rankSelect
    #    rankNr = 0 + as.numeric(substr(rankSelect,1,2))     # get rank number (number of column in unsorted taxa list - dataframe Dt)
    choice <- as.data.frame
    choice <- rbind(Dt[rankName])
    colnames(choice) <- "ncbiID"
    choice <- merge(choice,nameList,by="ncbiID",all = FALSE)
  })
  
  ### then output list of (super)taxa onto UI
  output$select = renderUI({
    choice <- allTaxaList()
    choice$fullName <- as.factor(choice$fullName)
    
    # if(input$demo == TRUE){
    if(input$demo_data == "demo"){
      hellemDf <- data.frame("name" = c("Encephalitozoon hellem","Encephalitozoon hellem","Encephalitozoon","Unikaryonidae","Apansporoblastina","Apansporoblastina","Microsporidia","Fungi","Eukaryota"),
                             "rank" = c("strain","species","genus","family","order","class","phylum","kingdom","superkingdom"))
      rankSelect = input$rankSelect
      rankName = substr(rankSelect,4,nchar(rankSelect))
      
      selectInput('inSelect',"",as.list(levels(choice$fullName)),hellemDf$name[hellemDf$rank == rankName])
    } else if(input$demo_data == "ampk-tor"){
      humanDf <- data.frame("name" = c("Homo sapiens","Homo sapiens","Homo","Hominidae","Primates","Mammalia","Chordata","Metazoa","Eukaryota"),
                            "rank" = c("strain","species","genus","family","order","class","phylum","kingdom","superkingdom"))
      rankSelect = input$rankSelect
      rankName = substr(rankSelect,4,nchar(rankSelect))
      
      selectInput('inSelect',"",as.list(levels(choice$fullName)),humanDf$name[humanDf$rank == rankName])
    } else {
      selectInput('inSelect',"",as.list(levels(choice$fullName)),levels(choice$fullName)[1])
    }
  })
  
  output$highlightTaxonUI = renderUI({
    choice <- allTaxaList()
    choice$fullName <- as.factor(choice$fullName)
    
    out <- as.list(levels(choice$fullName))
    out <- append("none",out)
    
    selectInput('taxonHighlight','Select (super)taxon to highlight:',out,selected=out[1])
  })
  
  #### update highlightTaxonUI based on double clicked dot
  observe({
    choice <- allTaxaList()
    choice$fullName <- as.factor(choice$fullName)
    
    out <- as.list(levels(choice$fullName))
    out <- append("none",out)
    
    if(!is.null(input$plot_dblclick)){
      if(input$xAxis == "genes"){
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
  output$highlightGeneUI = renderUI({
    geneList <- preData()
    geneList$geneID <- as.factor(geneList$geneID)
    
    out <- as.list(levels(geneList$geneID))
    out <- append("none",out)
    
    selectInput('geneHighlight','Highlight:',out,selected=out[1])
  })
  
  #### update highlightGeneUI based on double clicked dot
  observe({
    geneList <- preData()
    geneList$geneID <- as.factor(geneList$geneID)
    
    out <- as.list(levels(geneList$geneID))
    out <- append("none",out)
    
    if(!is.null(input$plot_dblclick)){
      if(input$xAxis == "genes"){
        corX = round(input$plot_dblclick$y);
        corY = round(input$plot_dblclick$x)
      } else {
        corX = round(input$plot_dblclick$x);
        corY = round(input$plot_dblclick$y)
      }
      updateSelectInput(session,'geneHighlight',label = 'Highlight:',choices=out,selected=out[corY+1])
    }
  })
  
  ######## print total number of genes
  output$totalGeneNumber.ui <- renderUI({
    geneList <- preData()
    geneList$geneID <- as.factor(geneList$geneID)
    out <- as.list(levels(geneList$geneID))
    if(length(out) > 0){
      # em(paste0("Total number of genes:  ",length(out)))
      strong(paste0("Total number of genes:  ",length(out)))
    }
  })
  
  ######## enable "PLOT" button
  observeEvent(input$rankSelect,({
    if(input$rankSelect == ""){
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
      toggleState("geneList_selected")
      toggleState("demo_data")
    }
  })
  
  ######## disable demo checkbox and update var2_aggregateBy to mean if using demo data
  observe({
    # if (input$demo == TRUE) {
    if (input$demo_data == "demo") {
      ### disable demo checkbox
      # toggleState("demo")
      ### update var2_aggregateBy to mean
      updateSelectInput(session,"var2_aggregateBy",
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
  
  observeEvent(input$newAdd, {
    newTaxa$Df[newIndex$value,] <- c(input$newID,input$newName,input$newRank,input$newParent)
    newIndex$value <- newIndex$value + 1
    updateTextInput(session, "newID", value=as.numeric(input$newID)+1)
    updateTextInput(session, "newName", value="")
    updateTextInput(session, "newRank", value="norank")
    updateTextInput(session, "newParent", value="")
  })
  
  observeEvent(input$newDone, {
    toggleModal(session, "addTaxaWindows", toggle = "close")
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
  
  sortedTaxaList <- reactive({
    if(v$doPlot == FALSE){return()}
    
    ####### Get representative taxon ####### 
    
    ### load list of unsorted taxa
    Dt <- as.data.frame(read.table("data/taxonomyMatrix.txt", sep='\t', header=T, stringsAsFactors=T))
    
    ### load list of taxon name
    nameList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T,fill = TRUE))
    nameList$fullName <- as.character(nameList$fullName)
    
    ### input parameters
    rankSelect = input$rankSelect
    rankName = substr(rankSelect,4,nchar(rankSelect))   # get rank name from rankSelect
    rankNr = 0 + as.numeric(substr(rankSelect,1,2))     # get rank number (number of column in unsorted taxa list - dataframe Dt)
    
    ### get selected supertaxon ID
    taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
    allTaxa <- allTaxaList()
    rankNameTMP <- allTaxa$rank[allTaxa$fullName == input$inSelect]
    if(rankName == "strain"){
      superID <- as.numeric(taxaList$ncbiID[taxaList$fullName == input$inSelect & taxaList$rank == "norank"])
    } else {
      superID <- as.numeric(taxaList$ncbiID[taxaList$fullName == input$inSelect & taxaList$rank == rankNameTMP[1]])
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
    inputTaxa <- subsetTaxa()
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
  output$applyClusterCheck.ui <- renderUI({
    if(input$ordering == FALSE){
      HTML('<p><em>(Check "Ordering sequence IDs" check box in <strong>Input & settings tab</strong>&nbsp;to enable this function)</em></p>')
    }
  })
  
  observe({
    if(input$ordering == FALSE){
      shinyjs::disable('applyCluster')
    } else {
      shinyjs::enable('applyCluster')
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
    endIndex <- input$endIndex
    if(is.na(input$endIndex)){endIndex <- 30}
    
    if(input$geneList_selected == 'from file'){
      listIn <- input$list
      if(!is.null(listIn)){
        list <- as.data.frame(read.table(file=listIn$datapath, header=FALSE))
        listGeneOri <- list$V1
        if(input$stIndex <= length(listGeneOri)){
          listGene <- listGeneOri[listGeneOri[input$stIndex:endIndex]]
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
        subsetID <- levels(as.factor(inputDf$geneID))[input$stIndex:endIndex]
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
          subsetID <- levels(longDf$geneID)[input$stIndex:endIndex]
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
          subsetID <- levels(longDf$geneID)[input$stIndex:endIndex]
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
          subsetID <- levels(inputDf$geneID)[input$stIndex:endIndex]
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
          subsetID <- levels(mdData$geneID)[input$stIndex:endIndex]
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
  dataFiltered <- reactive({
    mdData <- preData()
    ### count number of inparalogs
    paralogCount <- plyr::count(mdData,c('geneID','ncbiID'))
    mdData <- merge(mdData,paralogCount,by=c('geneID','ncbiID'))
    colnames(mdData)[ncol(mdData)] <- "paralog"
    
    ### (3) GET SORTED TAXONOMY LIST (3) ###
    taxaList <- sortedTaxaList()
    
    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxaList,'supertaxon')
    
    # merge mdData, mdDataVar2 and taxaList to get taxonomy info
    taxaMdData <- merge(mdData,taxaList,by='ncbiID')
    taxaMdData$var1 <- suppressWarnings(as.numeric(as.character(taxaMdData$var1)))
    taxaMdData$var2 <- suppressWarnings(as.numeric(as.character(taxaMdData$var2)))
    
    # ### (4) calculate PERCENTAGE of PRESENT SPECIES (4) ###
    finalPresSpecDt <- calcPresSpec(taxaMdData, taxaCount)
    
    ### (5) calculate max/min/mean/median VAR1 for every supertaxon of each gene (5) ###
    # remove NA rows from taxaMdData
    taxaMdDataNoNA <- taxaMdData[!is.na(taxaMdData$var1),]
    # calculate m var1
    mVar1Dt <- aggregate(taxaMdDataNoNA[,"var1"],list(taxaMdDataNoNA$supertaxon,taxaMdDataNoNA$geneID),FUN=input$var1_aggregateBy)
    colnames(mVar1Dt) <- c("supertaxon","geneID","mVar1")
    
    ### (6) calculate max/min/mean/median VAR2 for each super taxon (6) ###
    # remove NA rows from taxaMdData
    taxaMdDataNoNA_var2 <- taxaMdData[!is.na(taxaMdData$var2),]
    # calculate max/min/mean/median VAR2
    if(nrow(taxaMdDataNoNA_var2) > 0){
      mVar2Dt <- aggregate(taxaMdDataNoNA_var2[,"var2"],list(taxaMdDataNoNA_var2$supertaxon,taxaMdDataNoNA_var2$geneID),FUN=input$var2_aggregateBy)
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
    fullMdData <- dataFiltered()
    
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
  output$geneIn = renderUI({
    filein <- input$mainInput
    fileCustom <- input$customFile
    
    # if(input$demo == TRUE){
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      filein = 1
    }
    
    if(is.null(filein) & is.null(fileCustom)){return(selectInput('inSeq','',"all"))}
    if(v$doPlot == FALSE){return(selectInput('inSeq','',"all"))}
    else{
      ### full list
      data <- as.data.frame(dataFiltered())
      data$geneID <- as.character(data$geneID)
      data$geneID <- as.factor(data$geneID)
      outAll <- as.list(levels(data$geneID))
      outAll <- append("all",outAll)
      #selectInput('inSeq','',out,selected=out[1],multiple=TRUE)
      
      if (input$addCustomProfile == TRUE){
        out <- selectedGeneAge()
        if(length(out)>0){
          selectInput('inSeq','',out,selected=as.list(out),multiple=TRUE,selectize=FALSE)
        }
        else {
          selectInput('inSeq','',outAll,selected=outAll[1],multiple=TRUE,selectize=FALSE)
        }
      } else if(input$addClusterCustomProfile == TRUE){
        out <- brushedClusterGene()
        if(length(out)>0){
          selectInput('inSeq','',out,selected=as.list(out),multiple=TRUE,selectize=FALSE)
        }
        else {
          selectInput('inSeq','',outAll,selected=outAll[1],multiple=TRUE,selectize=FALSE)
        }
      }
      else if(input$addConsGeneCustomProfile == TRUE){
        out <- consGeneDf()
        if(length(out)>0){
          selectInput('inSeq','',out,selected=as.list(out),multiple=TRUE,selectize=FALSE)
        }
        else {
          selectInput('inSeq','',outAll,selected=outAll[1],multiple=TRUE,selectize=FALSE)
        }
      } else {
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
  output$taxaIn = renderUI({
    filein <- input$mainInput
    # if(input$demo == TRUE){
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      filein = 1
    }
    
    if(is.null(filein)){return(selectInput('inTaxa','',"all"))}
    if(v$doPlot == FALSE){return(selectInput('inTaxa','',"all"))}
    else{
      choice <- allTaxaList()
      choice$fullName <- as.factor(choice$fullName)
      
      out <- as.list(levels(choice$fullName))
      out <- append("all",out)
      
      if(input$applyCusTaxa == TRUE){
        out <- cusTaxaName()
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
    split <- strsplit(as.character(input$inSelect),"_")
    inSelect <- as.character(split[[1]][1])
    
    ### replace insufficient values according to the thresholds by NA or 0
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$presSpec < percent_cutoff_min] <- 0
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$presSpec > percent_cutoff_max] <- 0
    
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$var1 < var1_cutoff_min] <- 0
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$var1 > var1_cutoff_max] <- 0
    
    if(input$var1_relation == "protein"){
      if(input$var2_relation == "protein"){
        ### prot-prot: remove complete cell if one variable not sufficient
        dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$var2 < var2_cutoff_min] <- 0
        dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$var2 > var2_cutoff_max] <- 0
        
        dataHeat$var2[dataHeat$supertaxon != inSelect & dataHeat$var1 < var1_cutoff_min] <- NA
        dataHeat$var2[dataHeat$supertaxon != inSelect & dataHeat$var1 > var1_cutoff_max] <- NA
        dataHeat$var1[dataHeat$supertaxon != inSelect & dataHeat$var2 < var2_cutoff_min] <- NA
        dataHeat$var1[dataHeat$supertaxon != inSelect & dataHeat$var2 > var2_cutoff_max] <- NA
      } else {
        ### prot-spec: var1 depend on var2
        dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$var2 < var2_cutoff_min] <- 0
        dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$var2 > var2_cutoff_max] <- 0
      }
    } else {
      if(input$var2_relation == "species"){
        ### spec-spec: remove var1 and var2 independently
        # dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$var1 < var1_cutoff_min] <- 0
        # dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$var1 > var1_cutoff_max] <- 0
      } else {
        ### spec-prot: var2 depend on var1
        dataHeat$var2[dataHeat$supertaxon != inSelect & dataHeat$var1 < var1_cutoff_min] <- NA
        dataHeat$var2[dataHeat$supertaxon != inSelect & dataHeat$var1 > var1_cutoff_max] <- NA
      }
    }
    
    dataHeat$var1[dataHeat$supertaxon != inSelect & dataHeat$var1 < var1_cutoff_min] <- NA
    dataHeat$var1[dataHeat$supertaxon != inSelect & dataHeat$var1 > var1_cutoff_max] <- NA
    dataHeat$var2[dataHeat$supertaxon != inSelect & dataHeat$var2 < var2_cutoff_min] <- NA
    dataHeat$var2[dataHeat$supertaxon != inSelect & dataHeat$var2 > var2_cutoff_max] <- NA
    
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
    clusteredGeneIDs <- clusteredGeneList(dat,input$distMethod,input$clusterMethod)
    
    # sort original data according to clusteredGeneIDs
    dataHeat$geneID <- factor(dataHeat$geneID, levels = clusteredGeneIDs)
    
    return(dataHeat)
  })
  
  ########### render dot size to dotSizeInfo
  output$dotSizeInfo <- renderUI({
    if (v$doPlot == FALSE) return()
    
    dataHeat <- dataHeat()
    dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
    presentVl <- dataHeat$presSpec[!is.na(dataHeat$presSpec)]
    
    minDot <- (floor(min(presentVl)*10)/10*5)*(1+input$dotZoom)
    maxDot <- (floor(max(presentVl)*10)/10*5)*(1+input$dotZoom)
    
    em(paste0("current point's size: ",minDot," - ",maxDot))
  })
  
  ########### create profile heatmap
  mainPlot <- function(){
    if (v$doPlot == FALSE) return()
    dataHeat <- dataHeat()
    
    ### cluster dataHeat (if selected)
    if(input$applyCluster == TRUE){
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
    
    p <- heatmap.plotting(dataHeat,input$xAxis,input$var1_id,input$var2_id,input$lowColor_var1,input$highColor_var1,input$lowColor_var2,input$highColor_var2,input$paraColor,input$xSize,input$ySize,input$legendSize,input$mainLegend,input$dotZoom,input$xAngle,1)
    
    ### highlight taxon
    if(input$taxonHighlight != "none"){
      ## get selected highlight taxon ID
      rankSelect = input$rankSelect
      rankName = substr(rankSelect,4,nchar(rankSelect))   # get rank name from rankSelect
      taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
      taxonHighlight <- taxaList$ncbiID[taxaList$fullName == input$taxonHighlight & taxaList$rank == rankName]
      if(length(taxonHighlight) == 0L){
        taxonHighlight <- taxaList$ncbiID[taxaList$fullName == input$taxonHighlight]
      }
      
      ## get taxonID together with it sorted index
      highlightTaxon <- toString(dataHeat[dataHeat$supertaxonID == taxonHighlight,2][1])
      
      ## get index
      selectedIndex = as.numeric(as.character(substr(highlightTaxon,2,4)))
      
      ## draw a rect to highlight this taxon's column
      if(input$xAxis == "taxa"){
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
      if(input$xAxis == "taxa"){
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
    if(input$autoUpdate == FALSE){
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
    if(input$autoUpdate == FALSE){
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
      ## if autoupdate is NOT selected, use updateBtn to trigger plot changing
      if(input$autoUpdate == FALSE){
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
      ## if autoupdate is true
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
  output$plotDownload <- downloadHandler(
    filename = function() {c("plot.pdf")},
    content = function(file) {
      ggsave(file, plot = mainPlot(), width = input$width*0.056458333, height = input$height*0.056458333, units="cm", dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
  ######## get info clicked point on main heatmap
  mainPointInfo <- reactive({
    ### check input
    if (v$doPlot == FALSE) return()
    
    # get selected supertaxon name
    taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
    rankSelect = input$rankSelect
    rankName = substr(rankSelect,4,nchar(rankSelect))
    inSelect <- as.numeric(taxaList$ncbiID[taxaList$fullName == input$inSelect])
    
    dataHeat <- dataHeat()
    if(input$applyCluster == TRUE){
      dataHeat <- clusteredDataHeat()
    }
    
    ### get values
    if (is.null(input$plot_click$x)) {return()}
    else{
      ### get cooridiate point
      if(input$xAxis == "genes"){
        corX = round(input$plot_click$y);
        corY = round(input$plot_click$x)
      } else {
        corX = round(input$plot_click$x);
        corY = round(input$plot_click$y)
      }
      
      # get geneID
      # genes <- as.matrix(dataHeat[dataHeat$supertaxonID == inSelect & !is.na(dataHeat$presSpec),])
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
  #   info <- mainPointInfo()
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
      allData <- dataFiltered()
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
                     axis.title = element_text(size=input$dist_textSize),axis.text = element_text(size=input$dist_textSize)) +
        labs(x = paste0(input$var1_id," (mean = ",round(mean(splitDt$var1),3),")"), y = "Frequency")
      p
    }
  }
  
  output$var1DistPlot <- renderPlot(width = 512, height = 356,{
    if(input$autoUpdate == FALSE){
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
                     #                   plot.title = element_text(size=input$legendSize),#hjust = 0.5,
                     axis.title = element_text(size=input$dist_textSize),axis.text = element_text(size=input$dist_textSize)) +
        labs(x = paste0(input$var2_id," (mean = ",round(mean(splitDt$var2),3),")"), y = "Frequency")
      p
    }
  }
  
  output$var2DistPlot <- renderPlot(width = 512, height = 356,{
    if(input$autoUpdate == FALSE){
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
    taxaList <- sortedTaxaList()
    
    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxaList,'supertaxon')
    
    # merge mdData, mdDatavar2 and taxaList to get taxonomy info
    taxaMdData <- merge(mdData,taxaList,by='ncbiID')
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
                   axis.title = element_text(size=input$dist_textSize),axis.text = element_text(size=input$dist_textSize)) +
      labs(x = paste0("% present taxa (mean = ",round(mean(dt$presSpec),3),")"), y = "Frequency")
    p
  }
  
  output$presSpecPlot <- renderPlot(width = 512, height = 356,{
    if(input$autoUpdate == FALSE){
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
  output$plotDownload_dist <- downloadHandler(
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
  output$sameProfile <- reactive({
    if (v$doPlot == FALSE) return(FALSE)
    if(length(input$inSeq[1]) == 0){ return(FALSE)}
    else{
      if(input$inSeq[1] == "all" & input$inTaxa[1] == "all"){return(TRUE)}
    }
  })
  outputOptions(output, 'sameProfile', suspendWhenHidden=FALSE)
  
  ######## change label of plotCustom button if autoUpdateSelected is unchecked
  output$plotCustomBtn <- renderUI({
    if(input$autoUpdateSelected == FALSE){
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
  output$rankSelectCus = renderUI({
    mainRank <- input$rankSelect
    mainChoices = list("Strain"="05_strain","Species" = "06_species","Genus" = "10_genus", "Family" = "14_family", "Order" = "19_order", "Class" = "23_class",
                       "Phylum" = "26_phylum", "Kingdom" = "28_kingdom", "Superkingdom" = "29_superkingdom","unselected"="")
    cusChoices <- mainChoices[mainChoices >= mainRank]
    
    selectInput("rankSelectCus", label = h5("Select taxonomy rank:"),
                choices = as.list(cusChoices),
                selected = mainRank)
  })
  
  ######## print list of available taxa for customized plot (based on rank from rankSelectCus)
  taxaSelectCus <- reactive({
    rankSelectCus = input$rankSelectCus
    
    if(length(rankSelectCus) == 0){return()}
    else{
      ### load list of unsorted taxa
      Dt <- as.data.frame(read.table("data/taxonomyMatrix.txt", sep='\t', header=T, stringsAsFactors=T))
      Dt <- Dt[Dt$abbrName  %in% subsetTaxa(),]
      
      ### load list of taxon name
      nameList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T,fill = TRUE))
      nameList$fullName <- as.character(nameList$fullName)
      
      rankName = substr(rankSelectCus,4,nchar(rankSelectCus))   # get rank name from rankSelect
      choice <- as.data.frame
      choice <- rbind(Dt[rankName])
      colnames(choice) <- "ncbiID"
      choice <- merge(choice,nameList,by="ncbiID",all = FALSE)
      return(choice)
    }
  })
  
  output$taxaSelectCus = renderUI({
    choice <- taxaSelectCus()
    choice$fullName <- as.factor(choice$fullName)
    selectInput('taxaSelectCus',h5('Choose (super)taxon of interest:'),as.list(levels(choice$fullName)),levels(choice$fullName)[1])
  })
  
  ######## get list of taxa based on selected taxaSelectCus
  cusTaxaName <- reactive({
    
    taxaSelectCus = input$taxaSelectCus
    rankName = substr(input$rankSelectCus,4,nchar(input$rankSelectCus))
    
    if(taxaSelectCus == ""){return()}
    
    ### load list of unsorted taxa
    Dt <- as.data.frame(read.table("data/taxonomyMatrix.txt", sep='\t', header=T, stringsAsFactors=T))
    Dt <- Dt[Dt$abbrName  %in% subsetTaxa(),]
    
    ### get ID of customized (super)taxon
    taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T,fill = TRUE))
    superID <- taxaList$ncbiID[taxaList$fullName == taxaSelectCus & taxaList$rank %in% c(rankName,"norank")]
    
    ### from that ID, get list of all taxa for main selected taxon
    mainRankName = substr(input$rankSelect,4,nchar(input$rankSelect))
    customizedTaxaID <- levels(as.factor(Dt[mainRankName][Dt[rankName] == superID,]))
    
    cusTaxaName <- taxaList$fullName[taxaList$rank %in% c(mainRankName,"norank") & taxaList$ncbiID %in% customizedTaxaID]
    
    return(cusTaxaName)
  })
  
  ######## create plot (same as main plot)
  selectedPlot <- function(){
    if (vCt$doPlotCustom == FALSE) return()
    if(input$inSeq[1] == "all" & input$inTaxa[1] == "all") {return()}
    else{
      dataHeat <- dataHeat()
      
      ### cluster dataHeat (if selected)
      if(input$applyCluster == TRUE){
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
      p <- heatmap.plotting(dataHeat,input$xAxis_selected,input$var1_id,input$var2_id,input$lowColor_var1,input$highColor_var1,input$lowColor_var2,input$highColor_var2,input$paraColor,input$xSizeSelect,input$ySizeSelect,input$legendSizeSelect,input$selectedLegend,input$dotZoomSelect,input$xAngleSelect,0)
      
      ### do plotting
      if(input$autoUpdateSelected == FALSE){
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
  
  output$selectedPlot <- renderPlot({
    if(input$autoUpdateSelected == FALSE){
      # Add dependency on the update button (only update when button is clicked)
      input$plotCustom
      
      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        selectedPlot()
      })
    } else {
      selectedPlot()
    }
  })
  
  ######## plot selected sequences heatmap
  output$selectedPlot.ui <- renderUI({
    if(is.null(input$inSeq[1]) | is.null(input$inTaxa[1])){ return()}
    else if(input$inSeq[1] == "all" & input$inTaxa[1]=="all"){return()}
    else{
      if(input$autoUpdateSelected == FALSE){
        # Add dependency on the update button (only update when button is clicked)
        input$plotCustom
        
        # Add all the filters to the data based on the user inputs
        # wrap in an isolate() so that the data won't update every time an input
        # is changed
        isolate({
          withSpinner(
            plotOutput("selectedPlot",width=input$selectedWidth,height = input$selectedHeight,
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
          plotOutput("selectedPlot",width=input$selectedWidth,height = input$selectedHeight,
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
  output$selectedDownload <- downloadHandler(
    filename = function() {c("selected_plot.pdf")},
    content = function(file) {
      ggsave(file, plot = selectedPlot(), width = input$selectedWidth*0.056458333, height = input$selectedHeight*0.056458333, units="cm", dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
  ######## get info of a clicked point on selected plot (also the same as get info from main plot)
  selectedPointInfo <- reactive({
    ### check input
    if (vCt$doPlotCustom == FALSE) return()
    
    # get selected supertaxon name
    taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
    rankSelect = input$rankSelect
    rankName = substr(rankSelect,4,nchar(rankSelect))
    inSelect <- as.numeric(taxaList$ncbiID[taxaList$fullName == input$inSelect])
    
    dataHeat <- dataHeat()
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
      if(input$xAxis_selected == "genes"){
        corX = round(input$plot_click_selected$y);
        corY = round(input$plot_click_selected$x)
      } else {
        corX = round(input$plot_click_selected$x);
        corY = round(input$plot_click_selected$y)
      }
      
      # get geneID
      if(input$inSeq[1] == "all"){
        genes <- levels(dataHeat$geneID)
      } else {
        genes <- sort(input$inSeq)
      }
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
  
  ### get value of pointInfo for activating Detailed Plot button
  output$pointInfoStatus <- reactive({
    if(input$tabs == 'Main profile'){
      info <- mainPointInfo()  # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
    } else if(input$tabs=='Customized profile'){
      info <- selectedPointInfo()
    } else {
      info <- NULL
    }
    is.null(info)
  })
  outputOptions(output, "pointInfoStatus", suspendWhenHidden = FALSE)
  
  ######## show info into "point's info" box
  output$pointInfo <- renderText({
    ##### GET INFO BASED ON CURRENT TAB
    if(input$tabs == 'Main profile'){
      info <- mainPointInfo()  # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
    } else if(input$tabs=='Customized profile'){
      info <- selectedPointInfo()
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
          c <- toString(paste(input$var1_aggregateBy,input$var1_id,":",info[4]))
        }
        d <- ""
        if(input$var2_id != ""){
          d <- toString(paste(input$var2_aggregateBy,input$var2_id,":",info[6]))
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
  detailPlotDt <- reactive({
    if (v$doPlot == FALSE) return()
    
    ##### GET INFO BASED ON CURRENT TAB
    if(input$tabs == 'Main profile'){
      info <- mainPointInfo()  # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
    } else if(input$tabs=='Customized profile'){
      info <- selectedPointInfo()
    }
    
    if(is.null(info)){return()}
    else{
      ### get info for present taxa in selected supertaxon (1)
      plotTaxon = info[3]
      plotGeneID = info[1]
      fullDf <- dataFiltered()
      selDf <- as.data.frame(fullDf[fullDf$geneID == plotGeneID & fullDf$supertaxon == plotTaxon,])
      
      ### get all taxa of this supertaxon (2)
      allTaxaDf <- sortedTaxaList()
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
      if(input$detailedremoveNA == TRUE){
        joinedDf <- joinedDf[!is.na(joinedDf$orthoID),]
      }
      
      ### return data for detailed plot
      return(joinedDf)
    }
  })
  
  ######## render detailed plot
  detailPlot <- function(){
    if (v$doPlot == FALSE) return()
    
    selDf <- detailPlotDt()
    selDf$x_label <- paste(selDf$orthoID," (",selDf$fullName,")",sep = "")
    
    # if(input$detailedremoveNA == TRUE){
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
                  axis.text = element_text(size = input$detailedText),
                  axis.title = element_text(size = input$detailedText),
                  legend.text = element_text(size = input$detailedText)
    )
    gp
  }
  
  output$detailPlot <- renderPlot({
    p <- detailPlot()
    p
  })
  
  ######## plot detailed bar chart
  output$detailPlot.ui <- renderUI({
    withSpinner(
      plotOutput("detailPlot",width=800,height = input$detailedHeight,
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
  output$downloadDetailed <- downloadHandler(
    filename = function() {c("detailedPlot.pdf")},
    content = function(file) {
      g <- detailPlot()
      ggsave(file, plot = g, width = 800*0.056458333, height = input$detailedHeight*0.056458333, units="cm", dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
  ######## GET info when clicking on detailed plot
  pointInfoDetail <- reactive({
    selDf <- detailPlotDt()
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
  output$detailClick <- renderText({
    info <- pointInfoDetail() # info = seedID, orthoID, var1
    
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
    # dataOut <- as.data.frame(downloadCustomData())
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
    
    info <- pointInfoDetail() # info = seedID, orthoID, var1
    if(is.null(info)){return()}
    else{
      data <- dataFiltered()
      
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
    info <- pointInfoDetail() # info = seedID, orthoID, var1
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
      if(input$annoChoose == "from file"){
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
  output$checkDomainFiles <- renderUI({
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
  archiPlot <- function(){
    if (v3$doPlot3 == FALSE) return()
    
    ### info
    info <- pointInfoDetail() # info = seedID, orthoID, var1
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
    
    if(ncol(domainDf) == 5){
      colnames(domainDf) <- c("seedID","orthoID","feature","start","end")
    } else if(ncol(domainDf) == 6){
      colnames(domainDf) <- c("seedID","orthoID","feature","start","end","weight")
    } else if(ncol(domainDf) == 7){
      colnames(domainDf) <- c("seedID","orthoID","feature","start","end","weight","path")
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
      
      plot_ortho <- domain.plotting(orderedOrthoDf,ortho,sep,input$labelArchiSize,input$titleArchiSize,min(subDomainDf$start),max(subDomainDf$end))
      plot_seed <- domain.plotting(orderedSeedDf,seed,sep,input$labelArchiSize,input$titleArchiSize,min(subDomainDf$start),max(subDomainDf$end))
      
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
  
  output$archiPlot <- renderPlot({
    g <- archiPlot()
    grid.draw(g)
  })
  
  ######## render domain architecture plot
  output$archiPlot.ui <- renderUI({
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
      withSpinner(plotOutput("archiPlot",height = input$archiHeight, width = input$archiWidth))
    }
  })
  
  ######## download architecture plot ***** something strange with archiPlot()
  output$archiDownload <- downloadHandler(
    filename = function() {c("domains.pdf")},
    content = function(file) {
      g <- archiPlot()
      grid.draw(g)
      ggsave(file, plot = g, width = input$selectedWidth*0.056458333, height = input$selectedHeight*0.056458333, units="cm", dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
  #############################################################
  ######################## GENE AGE ###########################
  #############################################################
  
  ##### gene age estimation
  geneAgeDf <- reactive({
    if (v$doPlot == FALSE) return()
    
    rankList <- c("family","class","phylum","kingdom","superkingdom","root")
    
    ##### get selected (super)taxon ID
    rankSelect <- input$rankSelect
    rankName = substr(rankSelect,4,nchar(rankSelect))
    
    taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
    superID <- as.numeric(taxaList$ncbiID[taxaList$fullName == input$inSelect & taxaList$rank == rankName])
    
    ### full non-duplicated taxonomy data
    Dt <- as.data.frame(read.table("data/taxonomyMatrix.txt", sep='\t', header=T, stringsAsFactors=T))
    
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
    mdData <- droplevels(dataFiltered())
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
    geneAgeDf <- as.data.frame(tapply(mdDataExtended$cat, mdDataExtended$geneID, min))
    
    setDT(geneAgeDf, keep.rownames = TRUE)[]
    setnames(geneAgeDf, 1:2, c("geneID","cat"))  # rename columns
    row.names(geneAgeDf) <- NULL   # remove row names
    
    ### convert cat into geneAge
    geneAgeDf$age[geneAgeDf$cat == "0000001"] <- "07_LUCA"
    geneAgeDf$age[geneAgeDf$cat == "0000011" | geneAgeDf$cat == "0000010"] <- paste0("06_",as.character(taxaList$fullName[taxaList$ncbiID == subFirstLine$superkingdom & taxaList$rank == "superkingdom"]))
    geneAgeDf$age[geneAgeDf$cat == "0000111"] <- paste0("05_",as.character(taxaList$fullName[taxaList$ncbiID == subFirstLine$kingdom & taxaList$rank == "kingdom"]))
    geneAgeDf$age[geneAgeDf$cat == "0001111"] <- paste0("04_",as.character(taxaList$fullName[taxaList$ncbiID == subFirstLine$phylum & taxaList$rank == "phylum"]))
    geneAgeDf$age[geneAgeDf$cat == "0011111"] <- paste0("03_",as.character(taxaList$fullName[taxaList$ncbiID == subFirstLine$class & taxaList$rank == "class"]))
    geneAgeDf$age[geneAgeDf$cat == "0111111"] <- paste0("02_",as.character(taxaList$fullName[taxaList$ncbiID == subFirstLine$family & taxaList$rank == "family"]))
    geneAgeDf$age[geneAgeDf$cat == "1111111"] <- paste0("01_",as.character(taxaList$fullName[taxaList$fullName == input$inSelect & taxaList$rank == rankName]))
    
    ### return geneAge data frame
    geneAgeDf <- geneAgeDf[,c("geneID","cat","age")]
    
    geneAgeDf$age[is.na(geneAgeDf$age)] <- "Undef"
    return(geneAgeDf)
  })
  
  geneAgeDfMod <- reactive({
    geneAgeDf <- geneAgeDf()
    countDf <- plyr::count(geneAgeDf,c('age'))
    countDf$percentage <- round(countDf$freq/sum(countDf$freq)*100)
    countDf$pos <- cumsum(countDf$percentage) - (0.5 * countDf$percentage)
    return(countDf)
  })
  
  geneAgePlot <- function(){
    countDf <- geneAgeDfMod()
    p <- ggplot(countDf, aes(fill=age, y=percentage, x=1)) +
      geom_bar(stat="identity") +
      scale_y_reverse() +
      coord_flip() +
      theme_minimal()
    p <- p + geom_text(data=countDf, aes(x = 1, y = 100-pos, label = paste0(freq,"\n",percentage,"%")),size=4*input$geneAgeText)
    p <- p + theme(legend.position="bottom", legend.title = element_blank(), legend.text = element_text(size=12*input$geneAgeText),
                   axis.title = element_blank(), axis.text = element_blank()) +
      scale_fill_brewer(palette="Spectral") +
      guides(fill=guide_legend(nrow=round(nrow(countDf)/3,0),byrow=TRUE))
    p
  }
  
  output$geneAgePlot <- renderPlot({
    if(input$autoUpdate == FALSE){
      # Add dependency on the update button (only update when button is clicked)
      input$updateBtn
      
      # Add all the filters to the data based on the user inputs
      # wrap in an isolate() so that the data won't update every time an input
      # is changed
      isolate({
        geneAgePlot()
      })
    } else {
      geneAgePlot()
    }
  })
  
  output$geneAge.ui <- renderUI({
    if(v$doPlot == FALSE){
      return()
    } else{
      ## if autoupdate is NOT selected, use updateBtn to trigger plot changing
      if(input$autoUpdate == FALSE){
        # Add dependency on the update button (only update when button is clicked)
        input$updateBtn
        
        # Add all the filters to the data based on the user inputs
        # wrap in an isolate() so that the data won't update every time an input
        # is changed
        isolate({
          withSpinner(
            plotOutput("geneAgePlot",width=600*input$geneAgeWidth,height = 150*input$geneAgeHeight,
                       click = "plot_click_geneAge")
          )
        })
      }
      ## if autoupdate is true
      else {
        withSpinner(
          plotOutput("geneAgePlot",width=600*input$geneAgeWidth,height = 150*input$geneAgeHeight,
                     click = "plot_click_geneAge")
        )
      }
    }
  })
  
  ### download gene age plot
  output$geneAgePlotDownload <- downloadHandler(
    filename = function() {"geneAge_plot.pdf"},
    content = function(file) { 
      ggsave(file, plot = geneAgePlot(), width=600*input$geneAgeWidth*0.056458333,height = 150*input$geneAgeHeight*0.056458333, units="cm", dpi=300, device = "pdf")
    }
  )
  
  ##### render genAge.table based on clicked point on geneAgePlot
  selectedGeneAge <- reactive({
    if(v$doPlot == FALSE){return()}
    data <- geneAgeDf()
    
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
    if (is.null(input$plot_click_geneAge$x)) {return()}
    else{
      corX = 100-round(-input$plot_click_geneAge$x)
      selectAge <- as.character(rangeDf[rangeDf$rangeStart <= corX & rangeDf$rangeEnd >= corX,]$age)
      subData <- subset(data, age == selectAge)
      data <- data[data$age == selectAge,]
    }
    
    # return list of genes
    geneList <- levels(as.factor(subData$geneID))
    geneList
  })
  
  output$geneAge.table <- renderTable({
    if (is.null(input$plot_click_geneAge$x)) {return()}
    
    data <- as.data.frame(selectedGeneAge())
    data$number <- rownames(data)
    colnames(data) <- c("geneID","No.")
    data <- data[,c("No.","geneID")]
    data
  })
  
  ### download gene list from geneAgeTable
  output$geneAgeTableDownload <- downloadHandler(
    filename = function(){c("selectedGeneList.out")},
    content = function(file){
      dataOut <- selectedGeneAge()
      write.table(dataOut,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ### check if addClusterCustomProfile (profile clustering) are being clicked
  observe({
    if(input$addClusterCustomProfile == TRUE | input$addConsGeneCustomProfile == TRUE){
      shinyjs::disable('addCustomProfile')
    } else {
      shinyjs::enable('addCustomProfile')
    }
  })
  
  output$addCustomProfileCheck.ui <- renderUI({
    if(input$addClusterCustomProfile == TRUE  | input$addConsGeneCustomProfile == TRUE){
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Profile clustering</strong> or <strong>Core genes finding</strong>&nbsp;to enable this function)</em></p>')
    }
  })
  
  ### reset geneAgeProtConfig
  observeEvent(input$resetgeneAgeProtConfig, {
    shinyjs::reset("geneAgeWidth")
    shinyjs::reset("geneAgeHeight")
    shinyjs::reset("geneAgeText")
  })
  
  #############################################################
  ##################### CORE GENES ############################
  #############################################################
  
  ### render list of available taxa
  output$taxaList_cons.ui = renderUI({
    filein <- input$mainInput
    # if(input$demo == TRUE){
    if(input$demo_data == "demo" | input$demo_data == "ampk-tor"){
      filein = 1
    }
    
    if(is.null(filein)){return(selectInput('inTaxa','Select taxa of interest:',"none"))}
    if(v$doPlot == FALSE){return(selectInput('inTaxa','Select taxa of interest:',"none"))}
    else{
      choice <- allTaxaList()
      choice$fullName <- as.factor(choice$fullName)
      
      out <- as.list(levels(choice$fullName))
      out <- append("none",out)
      
      if(input$applyConsTaxa == TRUE){
        out <- consTaxaName()
        selectInput('taxaCons','Select taxa of interest:',out,selected=out,multiple=TRUE)
      } else {
        selectInput('taxaCons','Select taxa of interest:',out,selected=out[1],multiple=TRUE)
      }
    }
  })
  
  consGeneDf <- reactive({
    if (v$doPlot == FALSE) return()
    
    rankName = substr(input$rankSelect,4,nchar(input$rankSelect))
    
    ### get ID list of chosen taxa
    taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T,fill = TRUE))
    
    if("none" %in%input$taxaCons){superID = NA}
    else{superID <- taxaList$ncbiID[taxaList$fullName %in% input$taxaCons & taxaList$rank %in% c(rankName,"norank")]}
    
    ### get main input data
    mdData <- dataFiltered()
    mdData <- mdData[,c("geneID","ncbiID","fullName","supertaxon","supertaxonID","rank","presSpec","mVar1","mVar2")]
    
    ### filter by selecting taxa
    if(is.na(superID[1])){data <- NULL}
    else{
      data <- subset(mdData,supertaxonID %in% superID & presSpec >= input$percent_cons)
      # get supertaxa present in each geneID
      supertaxonCount <- as.data.frame(plyr::count(data,c('geneID','supertaxonID')))
      # count number of supertaxa present in each geneID and get only gene that contains all chosen taxa
      count <- as.data.frame(table(supertaxonCount$geneID))
      consGene <- subset(count,Freq == length(superID))
      consGene$Var1 <- factor(consGene$Var1)
      
      return(levels(consGene$Var1))
    }
  })
  
  output$consGene.table <- renderDataTable({
    data <- consGeneDf()
    if(is.null(data)){return()}
    else {
      data <- as.data.frame(data)
      # data$number <- rownames(data)
      # colnames(data) <- c("geneID","No.")
      # data <- data[,c("No.","geneID")]
      data
    }
  })
  
  ### download gene list from consGene.table
  output$consGeneTableDownload <- downloadHandler(
    filename = function(){c("consensusGeneList.out")},
    content = function(file){
      dataOut <- consGeneDf()
      write.table(dataOut,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ### check if addClusterCustomProfile (profile clustering) or addCustomProfile (gene age plot) are being clicked
  observe({
    if(input$addClusterCustomProfile == TRUE | input$addCustomProfile == TRUE){
      shinyjs::disable('addConsGeneCustomProfile')
    } else {
      shinyjs::enable('addConsGeneCustomProfile')
    }
  })
  
  output$addConsGeneCustomProfileCheck.ui <- renderUI({
    if(input$addClusterCustomProfile == TRUE | input$addCustomProfile == TRUE){
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Profiles clustering</strong> or <strong>Gene age estimating</strong>&nbsp;to enable this function)</em></p>')
    }
  })
  
  ######## print list of available taxonomy ranks (the lowest rank is the same as the chosen main rank)
  output$rankSelectCons = renderUI({
    mainRank <- input$rankSelect
    mainChoices = list("Strain"="05_strain","Species" = "06_species","Genus" = "10_genus", "Family" = "14_family", "Order" = "19_order", "Class" = "23_class",
                       "Phylum" = "26_phylum", "Kingdom" = "28_kingdom", "Superkingdom" = "29_superkingdom","unselected"="")
    consChoices <- mainChoices[mainChoices >= mainRank]
    
    selectInput("rankSelectCons", label = h5("Select taxonomy rank:"),
                choices = as.list(consChoices),
                selected = mainRank)
  })
  
  ######## print list of available taxa for customized plot (based on rank from rankSelectCus)
  taxaSelectCons <- reactive({
    rankSelectCons = input$rankSelectCons
    
    if(length(rankSelectCons) == 0){return()}
    else{
      ### load list of unsorted taxa
      Dt <- as.data.frame(read.table("data/taxonomyMatrix.txt", sep='\t', header=T, stringsAsFactors=T))
      Dt <- Dt[Dt$abbrName  %in% subsetTaxa(),]
      
      ### load list of taxon name
      nameList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T,fill = TRUE))
      nameList$fullName <- as.character(nameList$fullName)
      
      rankName = substr(rankSelectCons,4,nchar(rankSelectCons))   # get rank name from rankSelect
      choice <- as.data.frame
      choice <- rbind(Dt[rankName])
      colnames(choice) <- "ncbiID"
      choice <- merge(choice,nameList,by="ncbiID",all = FALSE)
      return(choice)
    }
  })
  
  output$taxaSelectCons = renderUI({
    choice <- taxaSelectCons()
    choice$fullName <- as.factor(choice$fullName)
    selectInput('taxaSelectCons',h5('Choose (super)taxon of interest:'),as.list(levels(choice$fullName)),levels(choice$fullName)[1])
  })
  
  ######## get list of taxa based on selected taxaSelectCus
  consTaxaName <- reactive({
    
    taxaSelectCons = input$taxaSelectCons
    rankName = substr(input$rankSelectCons,4,nchar(input$rankSelectCons))
    
    if(taxaSelectCons == ""){return()}
    
    ### load list of unsorted taxa
    Dt <- as.data.frame(read.table("data/taxonomyMatrix.txt", sep='\t', header=T, stringsAsFactors=T))
    Dt <- Dt[Dt$abbrName  %in% subsetTaxa(),]
    
    ### get ID of customized (super)taxon
    taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T,fill = TRUE))
    superID <- taxaList$ncbiID[taxaList$fullName == taxaSelectCons & taxaList$rank %in% c(rankName,"norank")]
    
    ### from that ID, get list of all taxa for main selected taxon
    mainRankName = substr(input$rankSelect,4,nchar(input$rankSelect))
    consTaxaID <- levels(as.factor(Dt[mainRankName][Dt[rankName] == superID,]))
    
    consTaxaName <- taxaList$fullName[taxaList$rank %in% c(mainRankName,"norank") & taxaList$ncbiID %in% consTaxaID]
    
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
    
    dd.col <- as.dendrogram(hclust(dist(dat, method = input$distMethod), method = input$clusterMethod))
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
      plotOutput("dendrogram",width=input$clusterPlot.width, height=input$clusterPlot.height,
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
  output$downloadCluster <- downloadHandler(
    filename = function() {"clustered_plot.pdf"},
    content = function(file) {
      ggsave(file, plot = dendrogram(), dpi=300, device = "pdf", limitsize=FALSE)
    }
  )
  
  #### render brushedCluster.table based on clicked point on dendrogram plot
  brushedClusterGene <- reactive({
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
  
  output$brushedCluster.table <- renderTable({
    if (is.null(input$plot_brush$ymin)) {return()}
    
    data <- as.data.frame(brushedClusterGene())
    data$number <- rownames(data)
    colnames(data) <- c("geneID","No.")
    data <- data[,c("No.","geneID")]
    data
  })
  
  ### download gene list from brushedCluster.table
  output$downloadClusterGenes <- downloadHandler(
    filename = function(){c("selectedClusteredGeneList.out")},
    content = function(file){
      dataOut <- brushedClusterGene()
      write.table(dataOut,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ### check if addCustomProfile (gene age plot) are being clicked
  observe({
    if(input$addCustomProfile == TRUE | input$addConsGeneCustomProfile == TRUE){
      shinyjs::disable('addClusterCustomProfile')
    }else{
      shinyjs::enable('addClusterCustomProfile')
    }
  })
  
  output$addClusterCustomProfileCheck.ui <- renderUI({
    if(input$addCustomProfile == TRUE | input$addConsGeneCustomProfile == TRUE){
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Gene age estimation</strong> or <strong>Core genes finding</strong>&nbsp;to enable this function)</em></p>')
    }
  })
  
  #############################################################
  ############### FILTERED DATA FOR DOWNLOADING ###############
  #############################################################
  
  ################### FOR MAIN PROFILE ########################
  
  ######## render variable used for identifying representative genes
  output$refVarMain.ui <- renderUI({
    if(nchar(input$var2_id) < 1 & nchar(input$var1_id) < 1){
      radioButtons(inputId="refVarMain", label="Reference variable", choices=list(input$var1_id,input$var2_id), selected=input$var1_id)
    } else if(nchar(input$var2_id) < 1){
      radioButtons(inputId="refVarMain", label="Reference variable", choices=list(input$var1_id), selected=input$var1_id)
    } else {
      radioButtons(inputId="refVarMain", label="Reference variable", choices=list(input$var1_id,input$var2_id), selected=input$var1_id)
    }
  })
  
  ######## filtered data for downloading
  downloadData <- reactive({
    ### check input
    if (v$doPlot == FALSE) return()
    
    ### filtered data
    dataOut <- dataFiltered()
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
    if(input$getRepresentativeMain == TRUE){
      if(is.null(input$refVarMain)){return()}
      else{
        if(input$refVarMain == input$var1_id){
          dataOutAgg <- aggregate(as.numeric(dataOut$var1), by = list(dataOut$geneID,dataOut$ncbiID), FUN = input$refTypeMain)
        } else if (input$refVarMain == input$var2_id){
          dataOutAgg <- aggregate(as.numeric(dataOut$var2), by = list(dataOut$geneID,dataOut$ncbiID), FUN = input$refTypeMain)
        } else {
          dataOutAgg <- dataOut[dataOut,c("geneID","ncbiID","var1")]
        }
        colnames(dataOutAgg) <- c("geneID","ncbiID","var_best")
        
        dataOutRepresentative <- merge(dataOut,dataOutAgg,by=c("geneID","ncbiID"),all.x=TRUE)
        
        if(input$refVarMain == input$var1_id){
          dataOut <- dataOutRepresentative[dataOutRepresentative$var1 == dataOutRepresentative$var_best,]
        } else if (input$refVarMain == input$var2_id){
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
  output$downloadData <- downloadHandler(
    filename = function(){c("filteredData.out")},
    content = function(file){
      dataOut <- downloadData()
      write.table(dataOut,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ######## data table ui tab
  output$filteredMainData <- renderDataTable(rownames = FALSE,{
    if(v$doPlot == FALSE){return()}
    #data <- taxaID()
    #data <- allTaxaList()
    #data <- sortedTaxaList()
    #data <- preData()
    #data <- dataFiltered()
    #data <- dataSupertaxa()
    # data <- dataHeat()
    #data <- detailPlotDt()
    #data <- presSpecAllDt()
    #data <- distDf()
    #data <- geneAgeDf()
    data <- downloadData()
    
    data
  })
  
  ######## download FASTA
  output$downloadFasta.ui <- renderUI({
    # if(input$demo_data == "demo"){
    #   HTML("<p><span style=\"color: #ff0000;\"><em>Be patient! For large number of taxa this can take up to 3 minutes!</em></span></p>")
    # }
  })
  
  output$downloadFasta <- downloadHandler(
    filename = function(){c("filteredSeq.fa")},
    content = function(file){
      fastaOutDf <- fastaOutData(as.data.frame(downloadData()))
      write.table(fastaOutDf,file,sep="\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
    }
  )
  
  ################### FOR CUSTOMIZED PROFILE ########################
  
  ######## render variable used for identifying representative genes
  output$representativeInfo.ui <- renderUI({
    msg <- paste0("NOTE: According to your choice in [Download filtered data -> Main data], only representative sequences with ",as.character(input$refTypeMain)," ",as.character(input$refVarMain),"  will be downloaded!")
    strong(em(msg),style = "color:red")
  })
  
  ######## filtered data for downloading
  downloadCustomData <- reactive({
    ### check input
    if (v$doPlot == FALSE) return()
    
    data <- as.data.frame(downloadData())
    
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
  output$downloadCustomData <- downloadHandler(
    filename = function(){c("customFilteredData.out")},
    content = function(file){
      dataOut <- downloadCustomData()
      write.table(dataOut,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ######## data table ui tab
  output$filteredCustomData <- renderDataTable(rownames = FALSE,{
    if(v$doPlot == FALSE){return()}
    data <- downloadCustomData()
    data
  })
  
  ######## download FASTA
  output$downloadCustomFasta.ui <- renderUI({
    # if(input$demo_data == "demo"){
    #   HTML("<p><span style=\"color: #ff0000;\"><em>Depend on the number of taxa, this might take up to 3 minutes!</em></span></p>")
    # }
  })
  
  output$downloadCustomFasta <- downloadHandler(
    filename = function(){c("customFilteredSeq.fa")},
    content = function(file){
      fastaOutDf <- fastaOutData(as.data.frame(downloadCustomData()))
      write.table(fastaOutDf,file,sep="\t",col.names = FALSE,row.names = FALSE,quote = FALSE)
    }
  )
})
