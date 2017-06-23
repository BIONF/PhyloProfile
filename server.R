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

#############################################################
######################## FUNCTIONS ##########################
#############################################################

########## convert long to wide format ##############
long2wide <- function(inputFile){
  longDf <- as.data.frame(read.table(file=inputFile$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
  
  # rename column names 
  colnames(longDf) <- c("geneID","ncbiID","orthoID","var1","var2")
  longDf$value <- paste0(longDf$orthoID,"#",longDf$var1,"#",longDf$var2)
  longDfmod <- longDf[,c("geneID","ncbiID","value")]
  wideDf <- spread(longDfmod, ncbiID, value)
}



########## calculate percentage of present species ##########
calcPresSpec <- function(taxaMdData, taxaCount){
  ### taxaMdData = df("geneID","ncbiID","orthoID","var1","var2",....,"supertaxon")
  # get geneID and supertaxon
  geneIDsupertaxon <- subset(taxaMdData,select=c('geneID','supertaxon'))
  geneIDsupertaxon <- geneIDsupertaxon[!duplicated(geneIDsupertaxon), ] # remove duplicated rows
  
  # remove NA rows from taxaMdData
  taxaMdDataNoNA <- taxaMdData[taxaMdData$orthoID != "NA",]
  #  taxaMdDataNoNA <- taxaMdData[!is.na(taxaMdData$var1),]
  #  taxaMdDataNoNA <- taxaMdData[!is.na(taxaMdData$var2),]
  
  # count present frequency of supertaxon for each gene
  geneSupertaxonCount <- plyr::count(taxaMdDataNoNA,c('geneID','supertaxon'))
  
  # merge with taxaCount to get total number of species of each supertaxon and calculate presSpec
  presSpecDt <- merge(geneSupertaxonCount,taxaCount,by='supertaxon')
  presSpecDt$presSpec <- presSpecDt$freq.x/presSpecDt$freq.y
  presSpecDt <- presSpecDt[order(presSpecDt$geneID),]
  presSpecDt <- presSpecDt[,c("geneID","supertaxon","presSpec")]
  
  # add absent supertaxon into presSpecDt
  finalPresSpecDt <- merge(presSpecDt,geneIDsupertaxon,by=c('geneID','supertaxon'),all.y = TRUE)
  finalPresSpecDt$presSpec[is.na(finalPresSpecDt$presSpec)] <- 0
  
  # return finalPresSpecDt
  finalPresSpecDt
}

######## function for sorting one domain dataframe (ortho) based on the other domain Df (seed) ########
sortDomains <- function(seedDf, orthoDf){
  # get list of features in seedDf
  featureList <- as.data.frame(levels(as.factor(seedDf$feature)))
  colnames(featureList) <- c("feature")
  # and add order number to each feature
  featureList$orderNo <- seq(length(featureList$feature))
  
  # merge those info to orthoDf
  orderedOrthoDf <- merge(orthoDf,featureList, all.x=TRUE)
  
  # sort orthoDf
  index <- with(orderedOrthoDf, order(orderNo))
  orderedOrthoDf <- orderedOrthoDf[index,]
  
  #turn feature column into a character vector
  orderedOrthoDf$feature <- as.character(orderedOrthoDf$feature)
  #then turn it back into an ordered factor (to keep this order while plotting)
  orderedOrthoDf$feature <- factor(orderedOrthoDf$feature, levels=unique(orderedOrthoDf$feature))
  
  #return sorted df
  orderedOrthoDf
}

######## function for plotting domain architecture ########
domain.plotting <- function(df,geneID,var1,sep,labelSize,titleSize,descSize,minStart,maxEnd){
  gg <- ggplot(df, aes(y=feature, x=end, color = feature)) +
    geom_segment(data=df, aes(y=feature, yend=feature, x=minStart, xend=maxEnd), color="#b2b2b2", size=0.15)
  
  ### draw line and points
  gg <- gg + geom_segment(data=df, aes(x=start, xend=end, y=feature, yend=feature),#, fill=feature),
                          size=1.2)
  gg <- gg + geom_point(data=df, aes(y=feature, x=start), color="#b2b2b2", size=3)
  gg <- gg + geom_point(data=df, aes(y=feature, x=end), color="#edae52", size=3)
  
  ### add text above
  gg <- gg + geom_text(data=df,
                       aes(x=(start+end)/2, y=feature, label=round(weight,2)),
                       color="#9fb059", size=descSize, vjust=-0.75, fontface="bold", family="serif")
  
  ### theme format
  titleMod <- gsub(":",sep,geneID)
  gg <- gg + scale_y_discrete(expand=c(0.075,0))
  gg <- gg + labs(title=paste0(titleMod," - ",var1))
  gg <- gg + theme_minimal()
  gg <- gg + theme(panel.border=element_blank())
  gg <- gg + theme(axis.ticks=element_blank())
  gg <- gg + theme(plot.title=element_text(face="bold",size=titleSize))
  gg <- gg + theme(plot.title=element_text(hjust = 0.5))
  gg <- gg + theme(legend.position="none",axis.title.x=element_blank(),
                   axis.text.y = element_text(size=labelSize),
                   axis.title.y = element_text(size=labelSize))
  
  ### return plot
  return(gg)
}

######## show FASTA sequence in popup windows of selected plot
getFasta <- function(file,seqID){
  fasta <- ""
  ### read file and get sequence
  if(file.exists(file)){
    fastaFile = readAAStringSet(file)
    
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)
    fa <- data.frame(seq_name, sequence)
    
    seq <- fa$sequence[fa$seq_name == seqID]
    if(length(seq) < 1){
      fasta <- paste0(seqID," not found in ",file,"! Please check id_format in FASTA config again!")
    } else{
      fasta <- paste(paste0(">",seqID),seq,sep="\n")
    }
  } else {
    fasta <- paste0(file," not found! Please check the path and dir_format in FASTA config again!")
  }
  ### return
  return(fasta)
}

######## get last n characters from string x
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

############################ MAIN ############################
options(shiny.maxRequestSize=30*1024^2)  ## size limit for input 30mb

shinyServer(function(input, output, session) {
  #  session$onSessionEnded(stopApp) ### Automatically stop a Shiny app when closing the browser tab
  
  ########## uncomment these lines for oneseq version#########
  # observeEvent(input$mainInput,({
  #   updateSelectInput(session,"var2_aggregateBy",
  #                     choices = list("Max"="max", "Min"="min","Mean"="mean","Median"="median"),
  #                     selected = "mean")
  #   updateTextInput(session,"var1_id", value = "FAS")
  #   updateTextInput(session,"var2_id", value = "Traceability")
  # }))
  #############################################################
  
  #############################################################
  ####################  PRE-PROCESSING  #######################
  #############################################################
  
  ################## getting ncbi taxa IDs ####################
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
            temp <- gnr_resolve(names = taxaNameDf[i,])
            newID <- get_uid(sciname = temp[1,3])[1]
            if(is.na(newID)){
              idDf[i,] <- c(as.character(taxaNameDf[i,]),as.character(temp[1,3]),paste0("NA"),"notfound")
            } else {
              idDf[i,] <- c(as.character(taxaNameDf[i,]),as.character(temp[1,3]),paste0("ncbi",newID),"notfound")
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
        colnames(notFoundDt) <- c("Summitted name","Matched name","Matched ID")
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
      colnames(notFoundDt) <- c("Summitted name","Matched name","Matched ID")
      
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
    filein <- input$mainInput
    if(is.null(filein)){return(textInput("var1_id", h5("First variable:"), value = "Variable 1", width="100%", placeholder="Name of first variable"))}    # get var1/var2 names based on input (only if input file in long format table)
    
    if(checkLongFormat() == TRUE){
      headerIn <- readLines(filein$datapath, n = 1)
      headerIn <- unlist(strsplit(headerIn,split = '\t'))
      
      textInput("var1_id", h5("First variable:"), value = headerIn[4], width="100%", placeholder="Name of first variable")
    } else {
      textInput("var1_id", h5("First variable:"), value = "Variable 1", width="100%", placeholder="Name of first variable")
    }
  })
  
  output$var2_id.ui <- renderUI({
    filein <- input$mainInput
    if(is.null(filein)){return(textInput("var2_id", h5("Second variable:"), value = "Variable 2", width="100%", placeholder="Name of second variable"))}    # get var1/var2 names based on input (only if input file in long format table)
    if(checkLongFormat() == TRUE){
      headerIn <- readLines(filein$datapath, n = 1)
      headerIn <- unlist(strsplit(headerIn,split = '\t'))
      
      textInput("var2_id", h5("Second variable:"), value = headerIn[5], width="100%", placeholder="Name of second variable")
    } else {
      textInput("var2_id", h5("Second variable:"), value = "Variable 2", width="100%", placeholder="Name of second variable")
    }
  })
  
  ######## get var1_id and var2_id values
  output$variableID <- reactive({
    ids <- as.list(c(input$var1_id,input$var2_id))
  })
  
  ######## variable 1 & 2 cutoff slidebar (main plot)
  output$var1_cutoff <- renderUI({
    sliderInput("var1",paste(input$var1_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(0.0,1.0), width = 200)
  })
  output$var2_cutoff <- renderUI({
    sliderInput("var2",paste(input$var2_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(0.0,1.0), width = 200)
  })
  
  ######## render filter slidebars for Customized plot
  output$var1Filter.ui <- renderUI({
    sliderInput("var1cus",paste(input$var1_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var1[1],input$var1[2]), width = 200)
  })
  
  output$var2Filter.ui <- renderUI({
    sliderInput("var2cus",paste(input$var2_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var2[1],input$var2[2]), width = 200)
  })
  
  output$percentFilter.ui <- renderUI({
    sliderInput("percent2",
                "% of present taxa:", min = 0, max = 1, step = 0.025, value = input$percent, width = 200)
  })
  
  ######## render filter slidebars for Distribution plot
  output$var1_dist.ui <- renderUI({
    sliderInput("var1_dist",paste(input$var1_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var1[1],input$var1[2]), width = 200)
  })

  output$var2_dist.ui <- renderUI({
    sliderInput("var2_dist",paste(input$var2_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var2[1],input$var2[2]), width = 200)
  })

  output$percent_dist.ui <- renderUI({
    sliderInput("percent_dist",
                "% of present taxa:", min = 0, max = 1, step = 0.025, value = input$percent, width = 200)
  })
  
  ######## render filter slidebars for Gene age estimation plot
  output$var1_age.ui <- renderUI({
    sliderInput("var1_age",paste(input$var1_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var1[1],input$var1[2]), width = 200)
  })
  
  output$var2_age.ui <- renderUI({
    sliderInput("var2_age",paste(input$var2_id,"cutoff:"), min = 0, max = 1, step = 0.025, value = c(input$var2[1],input$var2[2]), width = 200)
  })
  
  output$percent_age.ui <- renderUI({
    sliderInput("percent_age",
                "% of present taxa:", min = 0, max = 1, step = 0.025, value = input$percent, width = 200)
  })
  
  ######## update value for "main" filter slidebars based on "Customized", "Distribution" and "Gene age estimation" slidebars
  observe({
    newVar1 <- input$var1cus
    updateSliderInput(session, "var1", value = newVar1,
                      min = 0, max = 1, step = 0.025)
  })
  observe({
    newVar2 <- input$var2cus
    updateSliderInput(session, "var2", value = newVar2,
                      min = 0, max = 1, step = 0.025)
  })
  observe({
    newPercent <- input$percent2
    updateSliderInput(session, "percent", value = newPercent,
                      min = 0, max = 1, step = 0.025)
  })
  
  observe({
    newVar1 <- input$var1_dist
    updateSliderInput(session, "var1", value = newVar1,
                      min = 0, max = 1, step = 0.025)
  })
  observe({
    newVar2 <- input$var2_dist
    updateSliderInput(session, "var2", value = newVar2,
                      min = 0, max = 1, step = 0.025)
  })
  observe({
    newPercent <- input$percent_dist
    updateSliderInput(session, "percent", value = newPercent,
                      min = 0, max = 1, step = 0.025)
  })
  
  observe({
    newVar1 <- input$var1_age
    updateSliderInput(session, "var1", value = newVar1,
                      min = 0, max = 1, step = 0.025)
  })
  observe({
    newVar2 <- input$var2_age
    updateSliderInput(session, "var2", value = newVar2,
                      min = 0, max = 1, step = 0.025)
  })
  observe({
    newPercent <- input$percent_age
    updateSliderInput(session, "percent", value = newPercent,
                      min = 0, max = 1, step = 0.025)
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
  
  ######## reset all parameters of main plot
  observeEvent(input$resetMain, {
    shinyjs::reset("xSize")
    shinyjs::reset("ySize")
    shinyjs::reset("legendSize")
    shinyjs::reset("var1")
    shinyjs::reset("var2")
    shinyjs::reset("percent")
  })
  
  ######## reset all parameters of Customized plot
  observeEvent(input$resetSelected, {
    shinyjs::reset("xSizeSelect")
    shinyjs::reset("ySizeSelect")
    shinyjs::reset("legendSizeSelect")
    shinyjs::reset("var1")
    shinyjs::reset("var2")
    shinyjs::reset("percent")
  })
  
  ######## reset colors
  observeEvent(input$defaultColorTrace, {
    shinyjs::reset("lowColor_var2")
    shinyjs::reset("highColor_var2")
  })
  
  observeEvent(input$defaultColorVar1, {
    shinyjs::reset("lowColor_var1")
    shinyjs::reset("highColor_var1")
  })
  
  
  ########################################################
  
  ######## check if main input file is in long format
  checkLongFormat <- reactive({
    filein <- input$mainInput
    if(is.null(filein)){return()}
    
    inputDt <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
    if(is.na(pmatch("ncbi",colnames(inputDt)[3])) || is.na(pmatch("ncbi",colnames(inputDt)[4])) || is.na(pmatch("ncbi",colnames(inputDt)[5]))){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  
  ######## check if there is any "unknown" taxon in input matrix
  unkTaxa <- reactive({
    filein <- input$mainInput
    if(is.null(filein)){return()}
    
    # get list of all available taxon (from taxonID.list.fullRankID)
    allTaxa <- unlist(as.list(fread("data/taxonID.list.fullRankID",sep = "\t", select = "abbrName")))
    
    # get list of input taxa (from main input file)
    if(checkLongFormat() == TRUE){
      inputMod <- long2wide(filein)
      inputTaxa <- colnames(inputMod)
    } else {
      inputTaxa <- readLines(filein$datapath, n = 1)
    }
    
    inputTaxa <- unlist(strsplit(inputTaxa,split = '\t'))
    inputTaxa <- inputTaxa[-1]   # remove "geneID" element from vector inputTaxa
    
    # list of unknown taxa
    unkTaxa <- inputTaxa[!(inputTaxa %in% allTaxa)]
    # return list of unkTaxa
    unkTaxa
  })
  
  ### get status of unkTaxa for conditional panel
  output$unkTaxaStatus <- reactive({
    unkTaxa <- unkTaxa()
    length(unkTaxa) > 0
  })
  outputOptions(output, "unkTaxaStatus", suspendWhenHidden = FALSE)
  
  ### show full list of unkTaxa
  output$unkTaxaFull <- renderDataTable(option = list(searching = FALSE),{
    if(length(unkTaxa())>0){
      tb <- as.data.frame(unkTaxa())
      names(tb)[1] <- "New taxon"
      tb
    }
  })
  
  ######## get input taxa (a subset of available taxa in taxonID.list.fullRankID)
  subsetTaxa <- reactive({
    filein <- input$mainInput
    if(is.null(filein)){return()}
    
    if(length(unkTaxa()) == 0){
      # get list of input taxa (from main input file)
      if(checkLongFormat() == TRUE){
        inputMod <- long2wide(filein)
        inputTaxa <- colnames(inputMod)
      } else {
        inputTaxa <- readLines(filein$datapath, n = 1)
      }
      
      inputTaxa <- unlist(strsplit(inputTaxa,split = '\t'))
      inputTaxa <- inputTaxa[-1]   # remove "geneID" element from vector inputTaxa
      
      # return subset of taxonID.list.fullRankID
      inputTaxa
    }
  })
  
  ######## enable "add taxa", "parse" & "upload additional files" button after uploading main file
  observeEvent(input$mainInput, ({
    updateButton(session, "parse", disabled = FALSE)
    updateButton(session, "AddFile", disabled = FALSE)
    updateButton(session, "geneList", disabled = FALSE)
    updateButton(session, "addTaxa", disabled = FALSE)
  }))
  
  ######## check if data is loaded and "parse" button (get info from input) is clicked and confirmed
  v1 <- reactiveValues(parseInput = FALSE)
  observeEvent(input$BUTyes, {
    toggleModal(session, "parseConfirm", toggle = "close")
    v1$parseInput <- input$BUTyes
  })
  observeEvent(input$BUTno, {
    toggleModal(session, "parseConfirm", toggle = "close")
  })
  
  ######### create taxonID.list.fullRankID and taxonNamesReduced.txt from input file (if confirmed by BUTyes)
  observe({
    filein <- input$mainInput
    if(is.null(filein)){return()}
    else{
      if(v1$parseInput == FALSE){return()}
      else{
        if(checkLongFormat() == TRUE){
          inputMod <- long2wide(filein)
          titleline <- toString(paste(colnames(inputMod),collapse = "\t"))
        } else {
          titleline <- readLines(filein$datapath, n=1)
        }
        
        # Create 0-row data frame which will be used to store data
        dat <- data.frame(x = numeric(0), y = numeric(0))   ### use for progess bar
        withProgress(message = 'Parsing input file', value = 0, {
          cmd <- paste("perl ", getwd(),"/data/getTaxonomyInfo.pl",
                       " -i \"", titleline,"\"",
                       " -n ", getwd(),"/data/taxonNamesFull.txt",
                       " -a ", getwd(),"/data/newTaxa.txt",
                       " -o ", getwd(),"/data",
                       sep='')
          system(cmd)
        })
      }
    }
  })
  
  ######## list of taxonomy ranks for plotting
  output$rankSelect = renderUI({
    selectInput("rankSelect", label = h5("Select taxonomy rank:"),
                choices = list("Strain"="05_strain","Species" = "06_species","Genus" = "10_genus", "Family" = "14_family", "Order" = "19_order", "Class" = "23_class",
                               "Phylum" = "26_phylum", "Kingdom" = "28_kingdom", "Superkingdom" = "29_superkingdom","unselected"=""),
                selected = "06_species")
  })
  
  ####### GET list of all (super)taxa
  allTaxaList <- reactive({
    #  output$select = renderUI({
    filein <- input$mainInput
    if(is.null(filein)){return()}
    
    rankSelect = input$rankSelect
    if(rankSelect == ""){return()}
    
    ### load list of unsorted taxa
    Dt <- as.data.frame(read.table("data/taxonID.list.fullRankID", sep='\t',header=T))
    Dt <- Dt[Dt$abbrName  %in% subsetTaxa(),]
    
    ### load list of taxon name
    nameList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T,fill = TRUE))
    nameList$fullName <- as.character(nameList$fullName)
    
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
    selectInput('inSelect',h5('Choose (super)taxon of interest:'),as.list(levels(choice$fullName)),levels(choice$fullName)[1])
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
      toggleState("newTaxaAsk")
      toggleState("parseAsk")
    }
  })
  
  # #### disable clusterGene if input has only 1 gene
  # observe({
  #   filein <- input$mainInput
  #   if(is.null(filein)){return()}
  #   else{
  #     if(checkLongFormat() == TRUE){
  #       dt <- long2wide(filein)
  #     } else {
  #       dt <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
  #     }
  #     
  #     if(nrow(dt) < 2){
  #       updateRadioButtons(session,"ordering","",
  #                          choices = c("alphabetical","none"), selected = "alphabetical",
  #                          inline = F)
  #     } 
  #   }
  # })
  
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
    if(is.null(filein)){
      v$doPlot <- FALSE
      updateButton(session, "do", disabled = TRUE)
    }
  })
  
  ######## sorting supertaxa list based on chosen reference taxon
  nonDupTaxonDf <- reactive({
    if(v$doPlot == FALSE){return()}
    
    ### load list of unsorted taxa
    Dt <- as.data.frame(read.table("data/taxonID.list.fullRankID", sep='\t',header=T))
    names(Dt)[ncol(Dt)] <- "root"
    Dt <- Dt[Dt$abbrName  %in% subsetTaxa(),]
    
    ### reduce the number of columns in the original data frame by removing duplicate columns
    # transpose orig dataframe
    tDt <- as.data.frame(t(Dt))
    # get duplicate rows (columns in original matrix)
    dupDt <- tDt[duplicated(tDt), ]
    # exclude "main" ranks from dupDt
    dupDt_mod <- dupDt[row.names(dupDt) != "strain" & row.names(dupDt) != "species"
                       & row.names(dupDt) != "genus" & row.names(dupDt) != "family"
                       & row.names(dupDt) != "order" & row.names(dupDt) != "class"
                       & row.names(dupDt) != "phylum" & row.names(dupDt) != "kingdom"
                       & row.names(dupDt) != "superkingdom",]
    # transpose again
    tDupDt <- as.data.frame(t(dupDt_mod))
    # list of columns need to be dropped
    drop <- colnames(tDupDt)
    # drop those columns from original Df
    Dt <- Dt[,!(names(Dt) %in% drop)]
  })
  
  sortedTaxaList <- reactive({
    if(v$doPlot == FALSE){return()}
    
    ### load list of unsorted taxa
    Dt <- nonDupTaxonDf()
    
    ### load list of taxon name
    nameList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T,fill = TRUE))
    nameList$fullName <- as.character(nameList$fullName)
    
    ### input parameters
    rankSelect = input$rankSelect
    rankName = substr(rankSelect,4,nchar(rankSelect))   # get rank name from rankSelect
    #rankNr = 0 + as.numeric(substr(rankSelect,1,2))     # get rank number (number of column in unsorted taxa list - dataframe Dt)
    
    # get selected supertaxon ID
    taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
    if(rankName == "strain"){
      superID <- as.integer(taxaList$ncbiID[taxaList$fullName == input$inSelect & taxaList$rank == "norank"])
    } else {
      superID <- as.integer(taxaList$ncbiID[taxaList$fullName == input$inSelect & taxaList$rank == rankName])
    }
    
    ### sort taxa list
    ### first move all species that have the same ID of selected rank (level) to a new data frame
    ### and sort the rest of the data frame based on that rank.
    ### then move species that have different ID of selected rank (level), but have the same ID of the higher level
    ### repeat until reach to last column
    sortedDt <- data.frame()
    
    firstLine <- Dt[Dt[,rankName]==superID,][1,]  # get all taxon info for 1 representative
    sortedDt <- rbind(sortedDt,firstLine) # add to sortedDt
    Dt <- anti_join(Dt, firstLine, by="ncbiID") # & then remove that line from Dt
    
    for(i in 5:ncol(Dt)){
      Dt <- Dt[order(Dt[,i]),]
      matchID <- sortedDt[,i][1]
      subDt <- Dt[Dt[,i]==matchID,]
      if(nrow(subDt) > 0){
        sortedDt <- rbind(sortedDt,subDt)
        Dt <- anti_join(Dt, subDt, by="ncbiID")   # delete already selected lines from Dt dataframe
      }
    }
    
    ### join sortedDt and the rest of Dt list (species of other superkingdom than the one of selected supertaxon)
    sortedDt <- rbind(sortedDt,Dt)
    
    ### get only taxonIDs list of selected rank and rename columns
    sortedOut <- subset(sortedDt,select=c("No.","abbrName","ncbiID","fullName",as.character(rankName)))
    colnames(sortedOut) <- c("No.","abbrName","species","fullName","ncbiID")
    
    ### add name of supertaxa into sortedOut list
    sortedOut <- merge(sortedOut,nameList,by="ncbiID",all.x = TRUE,sort = FALSE)
    
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
    sortedOut <- sortedOut[,c("No.","abbrName","taxonID","fullName.x","species","ncbiID","sortedSupertaxon","rank","category")]
    colnames(sortedOut) <- c("No.","abbrName","taxonID","fullName","ncbiID","supertaxonID","supertaxon","rank","category")
    
    sortedOut$taxonID <- as.integer(sortedOut$taxonID)
    sortedOut$ncbiID <- as.factor(sortedOut$ncbiID)
    sortedOut$supertaxon <- as.factor(sortedOut$supertaxon)
    sortedOut$category <- as.factor(sortedOut$category)
    
    ### return data frame
    sortedOut
  })
  
  ############### FUNCTIONS FOR CLUSTERING PROFILES  ###############
  matrixForDitsCalc <- function(data){
    #### NOTE: input data for this function has to be in wide format
    mdData <- melt(data,id="geneID")
    
    # replace NA value with "NA#NA" (otherwise the corresponding orthoID will be empty)
    mdData$value <- as.character(mdData$value)
    mdData$value[is.na(mdData$value)] <- "NA#NA"
    
    # split value column into orthoID, var1 & var2
    splitDt <- (str_split_fixed(mdData$value, '#', 3))
    # then join them back to mdData
    mdData <- cbind(mdData,splitDt)
    # rename columns
    colnames(mdData) <- c("geneID","ncbiID","value","orthoID","var1","var2")
    clusterMdData <- mdData[,c("geneID","ncbiID","orthoID")]
    
    # replace 1 & 0 for presence / absence genes
    clusterMdData$group[clusterMdData$orthoID != "NA"] <-1
    clusterMdData$group[clusterMdData$orthoID == "NA"] <-0
    clusterMdData <- clusterMdData[,c("geneID","ncbiID","group")]
    
    # convert to wide format
    wideMdData <- spread(clusterMdData, ncbiID, group)
    
    dat <- wideMdData[,2:ncol(wideMdData)]  # numerical columns
    rownames(dat) <- wideMdData[,1]
    
    # return
    dat
  }
  
  clusterData <- function(data){
    #### NOTE: input data for this function has to be in wide format
    dat <- matrixForDitsCalc(data)  # convert into 0/1 matrix
    
    # do clustering
    row.order <- hclust(dist(dat, method = input$distMethod), method = input$clusterMethod)$order # clustering
    #    row.order <- hclust(dist(dat), method = input$clusterMethod)$order # clustering
    col.order <- hclust(dist(t(dat), method = input$distMethod), method = input$clusterMethod)$order
    #    col.order <- hclust(dist(t(dat)), method = input$clusterMethod)$order
    dat_new <- dat[row.order, col.order] # re-order dat accoring to clustering
    
    # get clustered gene IDs
    clusteredGeneIDs <- as.factor(row.names(dat_new))
    
    # sort original data according to clusteredGeneIDs
    data$geneID <- factor(data$geneID, levels = clusteredGeneIDs)
    
    # return clustered data
    data
  }
  
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
  
  ### subset data
  preData <- reactive({
    filein <- input$mainInput
    if(is.null(filein)){return()}
    
    # get rows need to be read
    nrHit <- input$stIndex + input$number - 1
    
    # convert input to wide format (if needed) & get nrHit rows
    if(checkLongFormat() == TRUE){
      inputMod <- long2wide(filein)
      if(input$applyCluster == TRUE){
        inputMod <- clusterData(inputMod)
      }
      #      data <- head(inputMod,nrHit)
      subsetID <- levels(inputMod$geneID)[1:nrHit]
      data <- inputMod[inputMod$geneID %in% subsetID,]
    } else {
      if(input$applyCluster == TRUE){
        oridata <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
        clusteredOridata <- clusterData(oridata)
        
        subsetID <- levels(clusteredOridata$geneID)[1:nrHit]
        data <- clusteredOridata[clusteredOridata$geneID %in% subsetID,]
      } else {
        data <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char="",nrows=nrHit))
      }
    }
    
    # OR just list of gene from a separated input file
    listIn <- input$list
    if(input$geneList_selected == 'from file'){
      if(!is.null(listIn)){
        list <- as.data.frame(read.table(file=listIn$datapath, header=FALSE))
        
        if(checkLongFormat() == TRUE){
          dataOrig <- long2wide(filein)
        } else {
          dataOrig <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
        }
        data <- dataOrig[dataOrig$geneID %in% list$V1,]
        
        if(input$applyCluster == TRUE){
          data <- clusterData(data)
        }
      }
    }
    
    if(input$ordering == FALSE){
      data$geneID <- factor(data$geneID, levels = data$geneID)  ######### keep user defined geneID order
    }
    
    ### return preData
    data
  })
  
  ### get all information for input data
  dataFiltered <- reactive({
    data <- preData()
    # convert into paired columns
    mdData <- melt(data,id="geneID")
    
    # replace NA value with "NA#NA" (otherwise the corresponding orthoID will be empty)
    mdData$value <- as.character(mdData$value)
    mdData$value[is.na(mdData$value)] <- "NA#NA"
    
    # split value column into orthoID, var1 & var2
    splitDt <- (str_split_fixed(mdData$value, '#', 3))
    # then join them back to mdData
    mdData <- cbind(mdData,splitDt)
    # rename columns
    colnames(mdData) <- c("geneID","ncbiID","value","orthoID","var1","var2")
    mdData <- mdData[,c("geneID","ncbiID","var1","orthoID","var2")]
    
    ### (3) GET SORTED TAXONOMY LIST (3) ###
    taxaList <- sortedTaxaList()
    
    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxaList,'supertaxon')
    
    # merge mdData, mdDataTrace and taxaList to get taxonomy info
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
    fullMdData ### parsed input data frame !!!
  })
  
  #############################################################
  ############## DATA & PLOT FOR MAIN PROFILE #################
  #############################################################
  
  ######## REDUCE DATA FROM SPECIES LEVEL TO SUPERTAXA LEVEL
  ######## this data set contain only supertaxa and their value (%present, mVar1 & mVar2) for each gene
  dataSupertaxa <- reactive({
    fullMdData <- dataFiltered()
    
    ### get representative orthoID that has m VAR1 for each supertaxon
    mOrthoID <- fullMdData[,c('geneID','supertaxon','var1','mVar1','orthoID')]
    mOrthoID <- subset(mOrthoID,mOrthoID$var1 == mOrthoID$mVar1)
    colnames(mOrthoID) <- c('geneID','supertaxon','var1','mVar1','orthoID')
    mOrthoID <- mOrthoID[!is.na(mOrthoID$orthoID),]
    mOrthoID <- mOrthoID[,c('geneID','supertaxon','orthoID')]
    mOrthoID <- mOrthoID[!duplicated(mOrthoID[,1:2]), ]
    
    ### get data set for phyloprofile plotting (contains only supertaxa info)
    superDf <- subset(fullMdData,select=c('geneID','supertaxon','supertaxonID','mVar1','presSpec','category','mVar2'))
    superDf <- superDf[!duplicated(superDf), ]
    superDfExt <- merge(superDf,mOrthoID, by=c('geneID','supertaxon'),all.x=TRUE)
    superDfExt <- superDfExt[,c("geneID","supertaxon","supertaxonID","mVar1","presSpec","category","orthoID","mVar2")]
    
    ### output
    names(superDfExt)[names(superDfExt)=="mVar1"] <- "var1"
    names(superDfExt)[names(superDfExt)=="mVar2"] <- "var2"
    superDfExt
  })
  
  
  ######## get list of all sequence IDs for selectize input
  output$geneIn = renderUI({
    filein <- input$mainInput
    fileCustom <- input$customFile
    
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
          selectInput('inSeq','',out,selected=as.list(out),multiple=TRUE)
        } 
        else {
          selectInput('inSeq','',outAll,selected=outAll[1],multiple=TRUE)
        }
      } else if(input$addClusterCustomProfile == TRUE){
        out <- brushedClusterGene()
        if(length(out)>0){
          selectInput('inSeq','',out,selected=as.list(out),multiple=TRUE)
        } 
        else {
          selectInput('inSeq','',outAll,selected=outAll[1],multiple=TRUE)
        }
      } else {
        if(is.null(fileCustom)){
          selectInput('inSeq','',outAll,selected=outAll[1],multiple=TRUE)
        }
        else {
          customList <- as.data.frame(read.table(file=fileCustom$datapath, header=FALSE))
          customList$V1 <- as.factor(customList$V1)
          out <- as.list(levels(customList$V1))
          selectInput('inSeq','',out,selected=out,multiple=TRUE)
        }
      }
    }
  })
  
  ######## get list of all taxa for selectize input
  output$taxaIn = renderUI({
    filein <- input$mainInput
    if(is.null(filein)){return(selectInput('inTaxa','',"all"))}
    if(v$doPlot == FALSE){return(selectInput('inTaxa','',"all"))}
    else{
      choice <- allTaxaList()
      choice$fullName <- as.factor(choice$fullName)
      
      out <- as.list(levels(choice$fullName))
      out <- append("all",out)

      if(input$applyCusTaxa == TRUE){
        out <- cusTaxaName()
        selectInput('inTaxa','',out,selected=out,multiple=TRUE)
      } else {
        selectInput('inTaxa','',out,selected=out[1],multiple=TRUE)
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
    if(is.null(filein)){return()}
    #      data <- read.table(file=filein$datapath, sep='\t',header=T)
    data <- dataSupertaxa()
    
    # get selected supertaxon name
    split <- strsplit(as.character(input$inSelect),"_")
    inSelect <- as.character(split[[1]][1])
    
    ### get sub set of data
    setID <- plyr::count(data,"geneID")
    #    subsetID <- setID[1:300,]
    nrHit <- input$stIndex + input$number - 1
    if(nrHit > nlevels(data$geneID)){nrHit <- nlevels(data$geneID)}
    
    subsetID <- setID[input$stIndex:nrHit,]
    dataHeat <- merge(data,subsetID,by="geneID")
    
    ### replace insufficient values according to the thresholds by NA or 0; and replace var1 0.0 by NA
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$presSpec < percent_cutoff_min] <- 0
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$presSpec > percent_cutoff_max] <- 0
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$var1 < var1_cutoff_min] <- 0
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$var1 > var1_cutoff_max] <- 0
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$var2 < var2_cutoff_min] <- 0
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$var2 > var2_cutoff_max] <- 0
    
    dataHeat$var1[dataHeat$supertaxon != inSelect & dataHeat$var1 < var1_cutoff_min] <- NA
    dataHeat$var1[dataHeat$supertaxon != inSelect & dataHeat$var1 > var1_cutoff_max] <- NA
    dataHeat$var2[dataHeat$supertaxon != inSelect & dataHeat$var2 < var2_cutoff_min] <- NA
    dataHeat$var2[dataHeat$supertaxon != inSelect & dataHeat$var2 > var2_cutoff_max] <- NA
    
    dataHeat <- droplevels(dataHeat)  ### delete unused levels
    dataHeat
  })
  
  ########### create profile heatmap
  mainPlot <- function(){
    if (v$doPlot == FALSE) return()
    dataHeat <- dataHeat()
    
    ### plot format
    if(input$xAxis == "genes"){
      p = ggplot(dataHeat, aes(x = geneID, y = supertaxon))        ## global aes
    } else{
      p = ggplot(dataHeat, aes(y = geneID, x = supertaxon))        ## global aes
    }
    
    if(length(unique(na.omit(dataHeat$var1))) == 1){
      mynewcolor_low <- input$highColor_var1
    } else {
      mynewcolor_low <- input$lowColor_var1
    }

    if(length(unique(na.omit(dataHeat$var2))) != 1){
      p = p + scale_fill_gradient(low = input$lowColor_var2, high = input$highColor_var2, na.value="gray95") +   ## fill color (var2)
        geom_tile(aes(fill = var2))    ## filled rect (var2 score)
    }
    p = p +  geom_point(aes(colour = var1, size = presSpec))  +    ## geom_point for circle illusion (var1 and presence/absence)
      scale_color_gradient(low = input$lowColor_var1,high = input$highColor_var1)#+       ## color of the corresponding aes (var1)
    scale_size(range = c(0,3))             ## to tune the size of circles
    #+ stat_binhex()
    p = p + guides(fill=guide_colourbar(title = input$var2_id), color=guide_colourbar(title = input$var1_id))   # thanks to Arpit Jain :-D
    base_size <- 9
    
    if(input$xAxis == "genes"){
      p = p + labs(y="Taxon")
      p = p+geom_hline(yintercept=0.5,colour="dodgerblue4")
      p = p+geom_hline(yintercept=1.5,colour="dodgerblue4")
    } else{
      p = p + labs(x="Taxon")
      p = p+geom_vline(xintercept=0.5,colour="dodgerblue4")
      p = p+geom_vline(xintercept=1.5,colour="dodgerblue4")
    }
    
    p = p+theme(axis.text.x = element_text(angle=60,hjust=1,size=input$xSize),axis.text.y = element_text(size=input$ySize),
                axis.title.x = element_text(size=input$xSize), axis.title.y = element_text(size=input$ySize),
                legend.title=element_text(size=input$legendSize),legend.text=element_text(size=input$legendSize),legend.position = input$mainLegend)
    
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
          # else, plot profile
          div(id = "plot-container",
              tags$img(src = "spinner.gif",
                       id = "loading-spinner"),
              #uiOutput("plot.ui")
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
        div(id = "plot-container",
            tags$img(src = "spinner.gif",
                     id = "loading-spinner"),
            #uiOutput("plot.ui")
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
  
  
  ########### plot guide axis lines into an absolute panel
  # ##### get margin size of main plot
  # mainMagin <- reactive({
  #   p <- mainPlot()
  #   fullSize <- as.list(par()$fin)
  #   plotSize <- as.list(par()$pin)
  #   print(par()$fin)
  # })
  
  # output$mainAxis <- renderPlot(bg="transparent",{
  #   if(input$autoUpdate == FALSE){
  #     # Add dependency on the update button (only update when button is clicked)
  #     input$updateBtn
  #     isolate({
  #       #list <- allTaxaList()
  #       p <- mainPlot()
  #       g <- ggplotGrob(p)
  # 
  #       if(input$mainXAxisGuide == TRUE & input$mainYAxisGuide == FALSE){
  #         s <- gtable_filter(g, 'axis-b', trim=F)  ### filter to get x-axis
  #       } else if (input$mainXAxisGuide == FALSE & input$mainYAxisGuide == TRUE){
  #         s <- gtable_filter(g, 'axis-l', trim=F)  ### filter to get y-axis
  #       } else if (input$mainXAxisGuide == TRUE & input$mainYAxisGuide == TRUE){
  #         s <- gtable_filter(g, 'axis-b|axis-l', trim=F)  ### filter to get x-axis and y-axis
  #       }
  # 
  #       # draw axis(es)
  #       grid.draw(s)
  #     })
  #   } else {
  #     p <- mainPlot()
  #     g <- ggplotGrob(p)
  # 
  #     if(input$mainXAxisGuide == TRUE & input$mainYAxisGuide == FALSE){
  #       s <- gtable_filter(g, 'axis-b', trim=F)  ### filter to get x-axis
  #     } else if (input$mainXAxisGuide == FALSE & input$mainYAxisGuide == TRUE){
  #       s <- gtable_filter(g, 'axis-l', trim=F)  ### filter to get y-axis
  #     } else if (input$mainXAxisGuide == TRUE & input$mainYAxisGuide == TRUE){
  #       s <- gtable_filter(g, 'axis-b|axis-l', trim=F)  ### filter to get x-axis and y-axis
  #     }
  #     grid.draw(s)
  #   }
  # })
  # 
  # output$mainAxisRender <- renderUI({
  #   if(input$autoUpdate == FALSE){
  #     # Add dependency on the update button (only update when button is clicked)
  #     input$updateBtn
  #     isolate({
  #       plotOutput("mainAxis", width=input$width, height=input$height)
  #     })
  #   } else{
  #     plotOutput("mainAxis", width=input$width, height=input$height)
  #   }
  # })
  
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
      var1 <- dataHeat$var1[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      Percent <- dataHeat$presSpec[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      Trace <- dataHeat$var2[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      # get ortholog ID
      orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      
      ### get list of all geneID that have the same ortholog
      geneMatch <- dataHeat$geneID[dataHeat$orthoID == toString(orthoID)]
      geneMatch <- geneMatch[!is.na(geneMatch)]
      # list of all available geneID
      geneList <- preData()
      geneList$geneID <- as.factor(geneList$geneID)
      allGenes <- as.list(levels(geneList$geneID))
      # get index of all matched genes (genes have the same ortholog)
      pos <- which(allGenes %in% geneMatch)
      pos <- paste(pos, collapse=',')
      
      ### return info of clicked point
      if(is.na(as.numeric(Percent))){return()}
      else{
        info <- c(geneID,as.character(orthoID),as.character(spec),round(as.numeric(var1),2),round(as.numeric(Percent),2),round(as.numeric(Trace),2),pos)
      }
    }
  })
  
  ######## get list of same orthologs (hit_IDs) of a selected point
  sameOrthoIndex <- reactive({
    ### check input
    if (v$doPlot == FALSE) return()
    
    ### info
    info <- mainPointInfo()
    pos <- info[7]
  })
  
  ###################################################################
  #### PLOT var1/var2 SCORE & % OF PRESENT SPECIES DISTRIBUTION #####
  ###################################################################
  
  ######## list of available variables for distribution plot
  output$selected.distribution = renderUI({
    varList <- as.list(c(input$var1_id,input$var2_id,"% present taxa"))
    selectInput('selected_dist','Choose variable to plot:',varList,varList[1])
  })
  
  ###### var1 / var2 distribution data
  distDf <- reactive({
    if (v$doPlot == FALSE) return()
    
    # open main input file
    filein <- input$mainInput
    if(checkLongFormat() == TRUE){
      dataOrig <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
      colnames(dataOrig) <- c("geneID","ncbiID","orthoID","var1","var2")
      splitDt <- dataOrig[,c("orthoID","var1","var2")]
    } else {
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

    splitDt$orthoID[splitDt$orthoID == "NA" | is.na(splitDt$orthoID)] <- NA
    splitDt <- splitDt[complete.cases(splitDt),]

    if(length(levels(as.factor(splitDt$var2))) == 1){
      if(levels(as.factor(splitDt$var2)) == ""){
        splitDt$var2 <- 0
      }
    }

    # convert factor into numeric for "var1" & "var2" column
    splitDt$var1 <- suppressWarnings(as.numeric(as.character(splitDt$var1)))
    splitDt$var2 <- suppressWarnings(as.numeric(as.character(splitDt$var2)))
    
    # filter splitDt based on selected var1 cutoff
    splitDt <- splitDt[splitDt$var1 >= input$var1[1] & splitDt$var1 <= input$var1[2],]
    # filter splitDt based on selected var2 cutoff
    splitDt <- splitDt[splitDt$var2 >= input$var2[1] & splitDt$var2 <= input$var2[2],]
    
    # return dt
    splitDt
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
    filein <- input$mainInput
    if(checkLongFormat() == TRUE){
      mdData <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
      colnames(mdData) <- c("geneID","ncbiID","orthoID","var1","var2")
    } else {
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
    
    ### (3) GET SORTED TAXONOMY LIST (3) ###
    taxaList <- sortedTaxaList()
    
    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxaList,'supertaxon')
    
    # merge mdData, mdDataTrace and taxaList to get taxonomy info
    taxaMdData <- merge(mdData,taxaList,by='ncbiID')
    taxaMdData$var1 <- suppressWarnings(as.numeric(as.character(taxaMdData$var1)))
    taxaMdData$var2 <- suppressWarnings(as.numeric(as.character(taxaMdData$var2)))
    
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
                   axis.title.x = element_text(size=input$dist_textSize),axis.text.x = element_text(size=input$dist_textSize)) +
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
      if(input$selected_dist == input$var1_id){
        plotOutput("var1DistPlot",width=input$width,height = input$height)
      }
      else if(input$selected_dist == input$var2_id){
        plotOutput("var2DistPlot",width=input$width,height = input$height)
      }
      else if(input$selected_dist == "% present taxa"){
        plotOutput("presSpecPlot",width=input$width,height = input$height)
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
  
  ######## change label of plotCustom button if autoUpdateSelected is unchecked
  output$plotCustomBtn <- renderUI({
    if(input$autoUpdateSelected == FALSE){
      bsButton("plotCustom", "Plot/Update selected sequence(s)/taxa",style="warning")
    } else {
      bsButton("plotCustom", "Plot selected sequence(s)/taxa",style="warning")
    }
  })
  
  
  ######## check if button is clicked
  v2 <- reactiveValues(doPlotCustom = FALSE)
  observeEvent(input$plotCustom, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    v2$doPlotCustom <- input$plotCustom
    filein <- input$mainInput
    if(is.null(filein)){v2$doPlotCustom <- FALSE}
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
      Dt <- as.data.frame(read.table("data/taxonID.list.fullRankID", sep='\t',header=T))
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
    Dt <- as.data.frame(read.table("data/taxonID.list.fullRankID", sep='\t',header=T))
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
  
  ######## update 
  
  ######## create plot (same as main plot)
  selectedPlot <- function(){
    if (v2$doPlotCustom == FALSE) return()
    if(input$inSeq[1] == "all" & input$inTaxa[1] == "all") {return()}
    else{
      dataHeat <- dataHeat()
      dataHeat$supertaxonMod <- substr(dataHeat$supertaxon,6,nchar(as.character(dataHeat$supertaxon)))
      if(input$inTaxa[1] == "all" & input$inSeq[1] != "all"){
        dataHeat <- subset(dataHeat,geneID %in% input$inSeq) ##### <=== select data from dataHeat for selected sequences only
      } else if(input$inSeq[1] == "all" & input$inTaxa[1] != "all"){
        dataHeat <- subset(dataHeat,supertaxonMod %in% input$inTaxa) ##### <=== select data from dataHeat for selected taxa only
      } else {
        dataHeat <- subset(dataHeat,geneID %in% input$inSeq & supertaxonMod %in% input$inTaxa) ##### <=== select data from dataHeat for selected sequences and taxa
      }
      
      ### plotting
      if(input$xAxis_selected == "taxa"){
        p = ggplot(dataHeat, aes(y = geneID, x = supertaxon))        ## global aes
      } else {
        p = ggplot(dataHeat, aes(x = geneID, y = supertaxon))
      }
      
      p = p + scale_fill_gradient(low = input$lowColor_var2, high = input$highColor_var2, na.value="gray95") +   ## fill color (var2)
        geom_tile(aes(fill = var2)) +# + scale_fill_gradient(low="gray95", high="red")) +    ## filled rect (var2 score)
        geom_point(aes(colour = var1, size = presSpec))  +    ## geom_point for circle illusion (var1 and presence/absence)
        scale_color_gradient(low = input$lowColor_var1,high = input$highColor_var1)#+       ## color of the corresponding aes (var1)
      scale_size(range = c(0,3))             ## to tune the size of circles
      p = p + guides(fill=guide_colourbar(title = input$var2_id), color=guide_colourbar(title = input$var1_id))
      base_size <- 9
      
      if(input$xAxis_selected == "taxa"){
        p = p + labs(x="Taxon")
        p = p+geom_vline(xintercept=0.5,colour="dodgerblue4")
        p = p+geom_vline(xintercept=1.5,colour="dodgerblue4")
      } else {
        p = p + labs(y="Taxon")
        p = p+geom_hline(yintercept=0.5,colour="dodgerblue4")
        p = p+geom_hline(yintercept=1.5,colour="dodgerblue4")
      }
      
      p = p+theme(axis.text.x = element_text(angle=60,hjust=1,size=input$xSizeSelect),axis.text.y = element_text(size=input$ySizeSelect),
                  axis.title.x = element_text(size=input$xSizeSelect), axis.title.y = element_text(size=input$ySizeSelect),
                  legend.title=element_text(size=input$legendSizeSelect),legend.text=element_text(size=input$legendSizeSelect),legend.position=input$selectedLegend)
      
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
  
  # ######## enable "plot selected sequences" button
  # observeEvent(input$inSeq, ({
  #   if(input$inSeq[1] != "all"){
  #     updateButton(session, "plotCustom", disabled = FALSE)
  #   }
  # }))
  #
  # observeEvent(input$inTaxa, ({
  #   if(input$inTaxa[1] != "all"){
  #     updateButton(session, "plotCustom", disabled = FALSE)
  #   }
  # }))
  
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
          div(id = "plot-container",
              tags$img(src = "spinner.gif",
                       id = "loading-spinner"
              ),
              
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
        div(id = "plot-container",
            tags$img(src = "spinner.gif",
                     id = "loading-spinner"
            ),
            
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
    if (v2$doPlotCustom == FALSE) return()
    
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
      var1 <- dataHeat$var1[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      Percent <- dataHeat$presSpec[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      Trace <- dataHeat$var2[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      # get ortholog ID
      orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      
      if(is.na(as.numeric(Percent))){return()}
      else{
        info <- c(geneID,as.character(orthoID),as.character(spec),round(as.numeric(var1),2),round(as.numeric(Percent),2),round(as.numeric(Trace),2))
      }
    }
  })
  
  
  #############################################################
  ################### SHOW CLICKED POINT INFO #################
  #############################################################
  
  ######## show info into "point's info" box
  output$pointInfo <- renderText({
    ##### GET INFO BASED ON CURRENT TAB
    if(input$tabs == 'Main profile'){
      info <- mainPointInfo()  # info = groupID,orthoID,supertaxon,mVar1,%spec,trace
    } else if(input$tabs=='Customized profile'){
      info <- selectedPointInfo()
    } else {
      return ()
    }
    
    if(is.null(info)){return()}
    else{
      orthoID <- info[2]
      ## parse orthoID for oneSeq
      if(input$input_type == 'oneSeq.extended.fa'){
        orthoIDTmp <- unlist(strsplit(toString(info[2]),"\\|"))
        #orthoID = toString(paste0(orthoIDTmp[2],":",orthoIDTmp[3]))
        orthoID = toString(orthoIDTmp[3])
      }
      if(orthoID=="NA"){orthoID <- info[2]}
      
      ## print output
      #a <- toString(paste(info[1],orthoID, sep = " ; "))
      a <- toString(paste("Seed-ID:",info[1]))
      #b <- toString(paste(substr(info[3],6,nchar(info[3]))))
      b <- toString(paste0("Hit-ID: ",orthoID," (",substr(info[3],6,nchar(info[3])),")"))
      c <- toString(paste(input$var1_aggregateBy,input$var1_id,":",info[4]))
      d <- toString(paste(input$var2_aggregateBy,input$var2_id,":",info[6]))
      e <- toString(paste("% present taxa:",info[5]))
      paste(a,b,c,d,e,sep="\n")
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
      info <- mainPointInfo()  # info = groupID,orthoID,supertaxon,mVar1,%spec,trace
    } else if(input$tabs=='Customized profile'){
      info <- selectedPointInfo()
    }
    
    if(is.null(info)){return()}
    else{
      plotTaxon = info[3]
      plotGeneID = info[1]
      
      fullDf <- dataFiltered()
      selDf <- as.data.frame(fullDf[fullDf$geneID == plotGeneID & fullDf$supertaxon == plotTaxon,])
      selDf
    }
  })
  
  ######## render detailed plot
  output$detailPlot <- renderPlot({
    if (v$doPlot == FALSE) return()
    
    selDf <- detailPlotDt()
    selDf$x_label <- paste(selDf$orthoID,"@",selDf$fullName,sep = "")
    
    ### create joined DF for plotting var1 next to var2
    var1Df <- subset(selDf, select = c("x_label","var1"))
    var1Df$type <- input$var1_id
    colnames(var1Df) <- c("id","score","var")
    var2Df <- subset(selDf, select = c("x_label","var2"))
    var2Df$type <- input$var2_id
    colnames(var2Df) <- c("id","score","var")
    
    detailedDf <- rbind(var1Df,var2Df)
    
    ### create plot
    gp = ggplot(detailedDf, aes(y=score,x=id,fill=var)) +
      geom_bar(stat="identity", position=position_dodge()) +
      coord_flip() +
      labs(x="") +
      theme_minimal()
    #geom_text(aes(label=var1), vjust=3)
    gp = gp+theme(axis.text.x = element_text(angle=90,hjust=1))
    gp
  })
  
  # ######## enable "detailed plot" button
  # observeEvent(pointInfo(), ({
  #   info <- pointInfo()
  #   if(length(info) > 2){
  #     updateButton(session, "go", disabled = FALSE)
  #   }
  # }))
  
  ######## plot detailed bar chart
  output$detailPlot.ui <- renderUI({
    plotOutput("detailPlot",width=800,height = input$detailedHeight,
               click = "plot_click_detail",
               hover = hoverOpts(
                 id = "plot_hover_2",
                 delay = input$hover_delay,
                 delayType = input$hover_policy,
                 nullOutside = input$hover_null_outside
               )
    )
  })
  
  ######## GET info when clicking on detailed plot
  pointInfoDetail <- reactive({
    selDf <- detailPlotDt()
    allOrthoID <- sort(selDf$orthoID)
    
    ### get coordinates of plot_click_detail
    if (is.null(input$plot_click_detail$x)) return()
    else{
      corX = round(input$plot_click_detail$y)
      corY = round(input$plot_click_detail$x)
    }
    
    ### get pair of sequence IDs & var1
    seedID <- toString(selDf$geneID[1])
    orthoID <- toString(allOrthoID[corX])
    var1 <- toString(selDf$var1[selDf$orthoID==orthoID])
    
    ### return info
    if(orthoID != "NA"){
      info <- c(seedID,orthoID,var1)
    }
  })
  
  ### SHOW info when clicking on detailed plot
  output$detailClick <- renderText({
    info <- pointInfoDetail() # info = seedID, orthoID, var1
    if(is.null(info)){return()}
    else{
      a <- paste0("seedID = ",info[1])
      b <- paste0("hitID = ",info[2])
      c <- paste0("var1 = ",info[3])
      paste(a,b,c,sep="\n")
    }
  })
  
  ######## FASTA sequence
  output$fasta <- renderText({
    if(v$doPlot == FALSE){return()}
    
    info <- pointInfoDetail() # info = seedID, orthoID, var1
    if(is.null(info)){return()}
    else{
      seqID <- toString(info[2])
      paste(seqID)
      
      ### fasta path and format
      if(input$input_type == 'oneSeq.extended.fa'){
        #        f <- toString(input$oneseq.file)
        fasIn <- input$oneSeqFasta
        f <- toString(fasIn$datapath)
      } else {
        path = input$path
        dir_format = input$dir_format
        file_ext = input$file_ext
        id_format = input$id_format
        
        ### get species ID and seqID
        specTMP <- unlist(strsplit(seqID,":"))
        specID = specTMP[1]
        if(id_format == 2){
          seqID = specTMP[2]
        }
        
        ### full path fasta file
        f <- paste0(path,"/",specID,".",file_ext)
        if(dir_format == 2){
          f <- paste0(path,"/",specID,"/",specID,".",file_ext)
        }
      }
      
      ### read file and get sequence
      ### get fasta
      paste(getFasta(f,seqID))
    }
  })
  
  
  #############################################################
  ################ FEATURE ARCHITECTURE PLOT ##################
  #############################################################
  
  ######## check clicked
  v3 <- reactiveValues(doPlot3 = FALSE)
  observeEvent(input$do3, {
    v3$doPlot3 <- input$do3
    filein <- input$mainInput
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
    filein3 <- input$fileDomain
    if(is.null(filein3)){
      v3$doPlot3 = FALSE
    } else {
      domainDf <- as.data.frame(read.table(file=filein3$datapath, sep='\t',header=FALSE,comment.char=""))
      colnames(domainDf) <- c("seedID","orthoID","feature","start","end","weight")
      
      if(input$input_type == 'oneSeq.extended.fa'){
        domainDf$seedID <- gsub("\\|",":",domainDf$seedID)
        domainDf$orthoID <- gsub("\\|",":",domainDf$orthoID)
      }
      
      ### get sub dataframe based on selected groupID and orthoID <<<<<============ HOLGER's NEW VERSION OF ONESEQ (14.03.17)
      # orthoNew <- ""
      # if(input$input_type == 'oneSeq.extended.fa'){
      #   last2char <- substrRight(ortho,2)
      #   if(last2char == '|0' | last2char == '|1'){
      #     orthoNew <- substr(ortho, 1, nchar(ortho)-2)
      #     orthoNew <- gsub("\\|",":",orthoNew)
      #   }
      # }
      # if(nchar(orthoNew)>0){ortho <- orthoNew}
      ortho <- gsub("\\|",":",ortho)
      
      grepID = paste(group,"#",ortho,sep="")
      subDomainDf <- domainDf[grep(grepID,domainDf$seedID),]
      
      ### ortho domains df
      orthoDf <- filter(subDomainDf,orthoID==ortho)
      orthoDf$feature <- as.character(orthoDf$feature)
      
      ### seed domains df
      seedDf <- filter(subDomainDf,orthoID != ortho)
      if(nrow(seedDf) == 0){seedDf <- orthoDf}
      
      seedDf$feature <- as.character(seedDf$feature)
      seed = as.character(seedDf$orthoID[1])
      
      ### change order of one dataframe's features based on order of other df's features
      if(length(orthoDf$feature) < length(seedDf$feature)){
        orderedOrthoDf <- orthoDf[order(orthoDf$feature), ]
        orderedSeedDf <- sortDomains(orderedOrthoDf, seedDf)
      } else {
        orderedSeedDf <- seedDf[order(seedDf$feature), ]
        orderedOrthoDf <- sortDomains(orderedSeedDf, orthoDf)
      }
      
      ### plotting
      sep = ":"
      if(!is.null(input$oneSeqFasta)){sep="|"}
      plot_ortho <- domain.plotting(orderedOrthoDf,ortho,var1,sep,input$labelArchiSize,input$titleArchiSize,input$labelDescSize,min(subDomainDf$start),max(subDomainDf$end))
      plot_seed <- domain.plotting(orderedSeedDf,seed,var1,sep,input$labelArchiSize,input$titleArchiSize,input$labelDescSize,min(subDomainDf$start),max(subDomainDf$end))
      
      # grid.arrange(plot_seed,plot_ortho,ncol=1)
      arrangeGrob(plot_seed,plot_ortho,ncol=1)
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
        "<p><span style=\"color: #ff0000;\"><strong>No information about domain architecture! Please check:</strong></span></p>
        <ul style=\"list-style-type: square;\">
        <li>if you selected any sequence in the Detailed plot?</li>
        <li>if you uploaded the domain file using Upload additional file(s) option? (see input example in data/demo/test.domains)</li>
        </ul>"
      )
      HTML(msg)
    } else {
      plotOutput("archiPlot",height = input$archiHeight, width = input$archiWidth)
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
    
    # filename = "domains.pdf",
    # content = function(file) {
    #   g <- archiPlot()
    #   grid.draw(g)
    #   ggsave(file, plot = g, width = input$archiWidth*0.056458333, height = input$archiHeight*0.056458333, units="cm", dpi=300)#, device = "svg")
    # }
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
    superID <- as.integer(taxaList$ncbiID[taxaList$fullName == input$inSelect & taxaList$rank == rankName])

    ### full non-duplicated taxonomy data
    Dt <- nonDupTaxonDf()
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
    mdData <- dataFiltered()
    mdData <- mdData[,c("geneID","ncbiID","orthoID","var1","var2","presSpec")]

    ### add "category" into mdData
    mdDataExtended <- merge(mdData,catDf,by="ncbiID",all.x = TRUE)
    mdDataExtended$var1[mdDataExtended$var1 == "NA" | is.na(mdDataExtended$var1)] <- 0
    mdDataExtended$var2[mdDataExtended$var2 == "NA" | is.na(mdDataExtended$var2)] <- 0
    
    ### remove cat for "NA" orthologs and also for orthologs that do not fit cutoffs
    mdDataExtended[mdDataExtended$orthoID == "NA"| is.na(mdDataExtended$orthoID),]$cat <- NA
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
    geneAgeDf$age[geneAgeDf$cat == "0000011"] <- paste0("06_",as.character(taxaList$fullName[taxaList$ncbiID == subFirstLine$superkingdom & taxaList$rank == "superkingdom"]))
    geneAgeDf$age[geneAgeDf$cat == "0000111"] <- paste0("05_",as.character(taxaList$fullName[taxaList$ncbiID == subFirstLine$kingdom & taxaList$rank == "kingdom"]))
    geneAgeDf$age[geneAgeDf$cat == "0001111"] <- paste0("04_",as.character(taxaList$fullName[taxaList$ncbiID == subFirstLine$phylum & taxaList$rank == "phylum"]))
    geneAgeDf$age[geneAgeDf$cat == "0011111"] <- paste0("03_",as.character(taxaList$fullName[taxaList$ncbiID == subFirstLine$class & taxaList$rank == "class"]))
    geneAgeDf$age[geneAgeDf$cat == "0111111"] <- paste0("02_",as.character(taxaList$fullName[taxaList$ncbiID == subFirstLine$family & taxaList$rank == "family"]))
    geneAgeDf$age[geneAgeDf$cat == "1111111"] <- paste0("01_",as.character(taxaList$fullName[taxaList$fullName == input$inSelect & taxaList$rank == rankName]))
    
    ### return geneAge data frame
    geneAgeDf <- geneAgeDf[,c("geneID","cat","age")]
    
    geneAgeDf$age[is.na(geneAgeDf$age)] <- "Undef"
    geneAgeDf
  })
  
  geneAgePlot <- function(){
    if (v$doPlot == FALSE) return()
    
    geneAgeDf <- geneAgeDf()
    
    countDf <- plyr::count(geneAgeDf,c('age'))
    countDf$percentage <- round(countDf$freq/sum(countDf$freq)*100)
    countDf$pos <- cumsum(countDf$percentage) - (0.5 * countDf$percentage)
    
    p <- ggplot(countDf, aes(fill=age, y=percentage, x=1)) + 
      geom_bar(stat="identity") +
      scale_y_reverse() +
      coord_flip() +
      theme_minimal()
    p <- p + geom_text(data=countDf, aes(x = 1, y = 100-pos, label = paste0(freq,"\n",percentage,"%")),size=4)
    p <- p + theme(legend.position="bottom", legend.title = element_blank(), legend.text = element_text(size=12),
                   axis.title = element_blank(), axis.text = element_blank()) +
      scale_fill_brewer(palette="Spectral") +
      guides(fill=guide_legend(nrow=3,byrow=TRUE))
    
    p
  }
  
  output$geneAgePlot <- renderPlot(height = 150, width = 600, {
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
          plotOutput("geneAgePlot",width=600,height = 150,
                     click = "plot_click_geneAge")
        })
      }
      ## if autoupdate is true
      else {
        plotOutput("geneAgePlot",width=600,height = 150,
                   click = "plot_click_geneAge")
      }
    }
  })
  
  ### download gene age plot
  output$geneAgePlotDownload <- downloadHandler(
    filename = function() {"geneAge_plot.pdf"},
    content = function(file) {
      ggsave(file, plot = geneAgePlot(), dpi=300, device = "pdf", limitsize=FALSE)
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
    if(input$addClusterCustomProfile == TRUE){
      shinyjs::disable('addCustomProfile')
    } else {
      shinyjs::enable('addCustomProfile')
    }
  })
  
  output$addCustomProfileCheck.ui <- renderUI({
    if(input$addClusterCustomProfile == TRUE){
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Cluster profiles</strong>&nbsp;to enable this function)</em></p>')
    }
  })
  
  #############################################################
  #################### CLUSTERING PROFILES ####################
  #############################################################
  
  # output$dist.table <- renderTable({
  #   data <- preData()
  #   dat <- matrixForDitsCalc(data)
  #   d <- dist(dat, method = input$distMethod)
  #   df <- melt(as.matrix(d), varnames = c("row", "col"))
  # })
  
  ### cluster data
  clusterDataDend <- reactive({
    if(v$doPlot == FALSE){return()}
    data <- preData()
    dat <- matrixForDitsCalc(data)
    dd.col <- as.dendrogram(hclust(dist(dat, method = input$distMethod), method = input$clusterMethod))
  })
  
  ### plot clustered profiles
  dendrogram <- function(){
    if(v$doPlot == FALSE){return()}
    
    dd.col <- clusterDataDend()
    
    py <- as.ggdend(dd.col)
    p <- ggplot(py, horiz = TRUE, theme=theme_minimal()) +
      theme(axis.title = element_blank(), axis.text.y = element_blank())
    p
  }
  
  output$dendrogram <- renderPlot({
    dendrogram()
  })
  
  output$cluster.ui <- renderUI({
    plotOutput("dendrogram",width=input$clusterPlot.width, height=input$clusterPlot.height,
               brush = brushOpts(
                 id = "plot_brush",
                 delay = input$brush_delay,
                 delayType = input$brush_policy,
                 direction = input$brush_dir,
                 resetOnNew = input$brush_reset)
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
      # a <- dt$labels$label[bottom]
      # b <- dt$labels$label[top]
      # values <- c(top,b,bottom,a)

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
    if(input$addCustomProfile == TRUE){
      shinyjs::disable('addClusterCustomProfile')
    }else{
      shinyjs::enable('addClusterCustomProfile')
    }
  })
  
  output$addClusterCustomProfileCheck.ui <- renderUI({
    if(input$addCustomProfile == TRUE){
      HTML('<p><em>(Uncheck "Add to Customized profile" check box in <strong>Gene age estimation</strong>&nbsp;to enable this function)</em></p>')
    }
  })
  
  #############################################################
  ############### FILTERED DATA FOR DOWNLOADING ###############
  #############################################################
  
  ################### FOR MAIN PROFILE ########################
  ######## filtered data for downloading
  downloadData <- reactive({
    ### check input
    if (v$doPlot == FALSE) return()
    
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
    
    dataOut <- dataOut[,c("geneID","orthoID","fullName","ncbiID","supertaxon","var1","var2","numberSpec","presSpec")]
    dataOut <- dataOut[order(dataOut$geneID,dataOut$supertaxon),]
    dataOut <- dataOut[complete.cases(dataOut),]
    
    dataOut$geneID <- as.character(dataOut$geneID)
    dataOut$fullName <- as.character(dataOut$fullName)
    dataOut$ncbiID <- substr(dataOut$ncbiID,5,nchar(as.character(dataOut$ncbiID)))
    dataOut$supertaxon <- substr(dataOut$supertaxon,6,nchar(as.character(dataOut$supertaxon)))
    dataOut$var1 <- as.character(dataOut$var1)
    dataOut$var2 <- as.character(dataOut$var2)
    dataOut$numberSpec <- as.integer(dataOut$numberSpec)
    dataOut$presSpec <- as.numeric(dataOut$presSpec)
    
    names(dataOut)[names(dataOut)=="presSpec"] <- "%Spec"
    names(dataOut)[names(dataOut)=="numberSpec"] <- "totalSpec"
    names(dataOut)[names(dataOut)=="var1"] <- input$var1_id
    names(dataOut)[names(dataOut)=="var2"] <- input$var2_id
    dataOut <- as.matrix(dataOut)
    dataOut
  })
  
  ######## download data
  output$downloadData <- downloadHandler(
    filename = function(){c("dataFiltered.out")},
    content = function(file){
      dataOut <- downloadData()
      write.table(dataOut,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ######## data table ui tab
  output$filteredMainData <- renderDataTable({
    if(v$doPlot == FALSE){return()}
    #data <- taxaID()
    #data <- allTaxaList()
    data <- sortedTaxaList()
    #data <- preData()
    #data <- dataFiltered()
    #data <- dataSupertaxa()
    #data <- dataHeat()
    #data <- detailPlotDt()
    #data <- presSpecAllDt()
    #data <- distDf()
    #data <- downloadData()
    data
  })
  
  ################### FOR CUSTOMIZED PROFILE ########################
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
        customData <- subset(data,supertaxonMod %in% input$inTaxa) ##### <=== select data for selected taxa only
      } else if(input$inSeq[1] != "all" & input$inTaxa[1] != "all") {
        customData <- subset(data,geneID %in% input$inSeq & supertaxonMod %in% input$inTaxa) ##### <=== select data for selected sequences and taxa
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
    filename = function(){c("customDataFiltered.out")},
    content = function(file){
      dataOut <- downloadData()
      write.table(dataOut,file,sep="\t",row.names = FALSE,quote = FALSE)
    }
  )
  
  ######## data table ui tab
  output$filteredCustomData <- renderDataTable({
    if(v$doPlot == FALSE){return()}
    data <- downloadCustomData()
    data
  })
  
  #############################################################
  ############### HELP & TEXT OUTPUT for TESTING ##############
  #############################################################
  
  ######## show help
  output$help.ui <- renderUI({
    HTML(
      '
      <h1 style="color: #5e9ca0;">How the input file looks like?</h1>
      <p>The main Input file is a matrix of "values" (e.g. FAS scores, normalized distances, etc.), where rows represent genes and columns represent taxa.</p>
      <p>A gene may be present or absent in some taxa. A present "value" has to be in the range of 0 and 1. An absent "value" is written as NA. Gene ID and its value is concatenated via a "#" symbol.&nbsp;For example:</p>
      <ul style="list-style-type: square;">
      <li>gene0123#0.9837: gene0123 has a value of 0.9837</li>
      <li>gene0999#NA: gene0999 is present but doesn\'t has any value</li>
      <li>NA#NA: there is no ortholog has been found for this taxon</li>
      </ul>
      <p>The header of first column has to be "geneID". The header of each taxon must have this format "ncbi12345", in which 12345 is its NCBI taxon ID.</p>
      <p><em>More detail? Pleas take a look at the example file in /data/demo/test.main :)</em></p>
      <p>&nbsp;</p>
      <h1 style="color: #5e9ca0;">Additional files</h1>
      <p>2 additional input files can be provided are traceability score matrix and feature domain position list.</p>
      <p><strong>Traceability score matrix</strong> must have the same first row (beginning with "geneID" and followed by list of taxa) and first column (list of genes). <span style="color: #ff0000;"><strong>IMPORTANT</strong>: the <span style="text-decoration: underline;">amount</span> and <span style="text-decoration: underline;">order</span> of genes between 2 input matrixes have to be exactly the same!!</span>&nbsp;</p>
      <p><strong>Feature domain position list</strong>&nbsp;has 6 columns separated by tab: (1) pairID = groupID#searchProt_ID#seedID, (2) searchProt_ID, (3) feature name (pfam domain, smart domain,etc.), (4) start position, (5) end position, (6) weight value (only available for seed protein)</p>
      <p><em>Pleas take a look at the example files &nbsp;test.traceability and test.domains in data/demo/ folder for more details :)</em></p>
      <p>&nbsp;</p>
      <h1 style="color: #5e9ca0;">Download function does not work</h1>
      <p>Problem: clicked on the "Download plot" (or Download filtered data) button, entered a file name on to "Download file" window and clicked Save, but the file...was not saved :(</p>
      <p>=&gt; Click on "Open im Browser" to open the app using internet browser. Now the download function should work.</p>
      <p><em>I tested this function using Ubuntu 14.04 LTS and it worked with Firefox web browser.</em>&nbsp;</p>
      <p><em>*** Download function for selected sequences plot does not work :(</em></p>
      <p>&nbsp;</p>
      <h1 style="color: #5e9ca0;">Errors while plotting</h1>
      <p><span style="color: #ff0000;"><strong>Error</strong>: arguments imply differing number of rows: 0, 1</span></p>
      <p>=&gt; does your input matrix contain only 1 line???</p>
      <p><span style="color: #ff0000;"><strong>Error</strong>: id variables not found in data: geneID</span></p>
      <p>=&gt; &nbsp;the header of gene column has to be "geneID".</p>
      <p><span style="color: #ff0000;"><strong>Error</strong>: second argument must be a list</span></p>
      <p>=&gt; please check if any of input genes is absent in all taxa, i.e. the complete row is filled with NA. The tool, unfortunately, does not accept these genes.</p>
      <p>&nbsp;</p>
      <h1 style="color: #5e9ca0;">Errors by&nbsp;showing Point\'s info</h1>
      <p><span style="color: #ff0000;">Error:&nbsp;argument&nbsp;is of length zero</span></p>
      <p>=&gt;&nbsp;re-select (super)taxon to highlight, it will work again :)</p>
      <p>&nbsp;</p>
      <h1 style="color: #5e9ca0;">Errors by plotting&nbsp;domain architecture</h1>
      <p><span style="color: #ff0000;">Error:&nbsp;Aesthetics must be either length 1 or the same as the data (1): x, y, yend, xend</span></p>
      <p>=&gt; please check if the ID of selected sequence, which is shown as orthoID in Detailed plot, is present in domain input file. This error is certainly caused by the missing ID in domain input file.</p>
      <p>&nbsp;</p>
      <h1 style="color: #5e9ca0;">Bug reporting</h1>
      <p>Any bug reports or comments, suggestions are highly appreciated ;-)</p>
      <p>&nbsp;</p>
      <p>&copy; 2016 Vinh Tran</p>
      <p>contact:&nbsp;<a href="mailto:tran@bio.uni-frankfurt.de">tran@bio.uni-frankfurt.de</a></p>
      <p>Please check the latest version at&nbsp;<a href="https://github.com/trvinh/phyloprofile">https://github.com/trvinh/phyloprofile</a></p>
      '
      )
  })

  ############### USED FOR TESTING
  output$testOutput <- renderText({
    # ### print infile
    # filein <- input$mainInput
    # print(toString(filein))
    # filePath <- toString(filein)
    # fileName <- unlist(strsplit(toString(input$mainInput),","))
    # name <- toString(fileName[1])
    # fullPath <- paste0("data/",name,".mDomains")
    # print(fullPath)
    # print(input$plot_dblclick$x)
    # paste(input$var1[1],input$var1[2])
  })
})
