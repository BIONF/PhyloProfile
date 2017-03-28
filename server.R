if (!require("shiny")) {install.packages("shiny")}
if (!require("shinyBS")) {install.packages("shinyBS")}
if (!require("ggplot2")) {install.packages("ggplot2")}
if (!require("reshape")) {install.packages("reshape")}
if (!require("plyr")) {install.packages("plyr")}
if (!require("dplyr")) {install.packages("dplyr")}
if (!require("scales")) {install.packages("scales")}
if (!require("grid")) {install.packages("grid")}
if (!require("gridExtra")) {install.packages("gridExtra")}
if (!require("ape")) {install.packages("ape")}
if (!require("stringr")) {install.packages("stringr")}
if (!require("svglite")) {install.packages("svglite")}
if (!require("gtable")) {install.packages("gtable")}
if (!require("data.table")) {install.packages("data.table")}
if (!require("Biostrings")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("Biostrings")
}

############################################################# 
######################## FUNCTIONS ##########################
#############################################################

######## function for plotting domain architecture ########
plotting <- function(df,geneID,fas,sep,labelSize,titleSize,minStart,maxEnd){
  gg <- ggplot(df, aes(y=feature, x=end, color = feature)) +
    geom_segment(data=df, aes(y=feature, yend=feature, x=minStart, xend=maxEnd), color="#b2b2b2", size=0.15)
  
  ### draw line and points  
  gg <- gg + geom_segment(data=df, aes(x=start, xend=end, y=feature, yend=feature, fill=feature),
                          size=1.2)
  gg <- gg + geom_point(data=df, aes(y=feature, x=start), color="#b2b2b2", size=3)
  gg <- gg + geom_point(data=df, aes(y=feature, x=end), color="#edae52", size=3)
  
  ### add text above
  gg <- gg + geom_text(data=df,
                       aes(x=(start+end)/2, y=feature, label=round(weight,2)),
                       color="#9fb059", size=2.5, vjust=-0.75, fontface="bold", family="Calibri")
  
  ### theme format
  titleMod <- gsub(":",sep,geneID)
  gg <- gg + scale_y_discrete(expand=c(0.075,0))
  gg <- gg + labs(title=paste0(titleMod," - FAS=",fas))
  gg <- gg + theme_bw(base_family="Calibri")
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
  
  ############################################################# 
  ####################  PRE-PROCESSING  #######################
  #############################################################
  
  ######## render filter slidebars for Customized plot
  output$fasFilter.ui <- renderUI({
    sliderInput("fas2",
                "FAS cutoff: ", min = 0, max = 1, step = 0.025, value = input$fas, width = 200)
  })
  
  output$percentFilter.ui <- renderUI({
    sliderInput("percent2",
                "% of present species:", min = 0, max = 1, step = 0.025, value = input$percent, width = 200)
  })
  ######## update value for "main" filter slidebars based on "Customized" slidebars
  observe({
    newFas <- input$fas2
    updateSliderInput(session, "fas", value = newFas,
                      min = 0, max = 1, step = 0.025)
  })
  observe({
    newPercent <- input$percent2
    updateSliderInput(session, "percent", value = newPercent,
                      min = 0, max = 1, step = 0.025)
  })
  
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
    shinyjs::reset("fas")
    shinyjs::reset("percent")
  })
  
  ######## reset all parameters of Customized plot
  observeEvent(input$resetSelected, {
    shinyjs::reset("xSizeSelect")
    shinyjs::reset("ySizeSelect")
    shinyjs::reset("legendSizeSelect")
    shinyjs::reset("fas")
    shinyjs::reset("percent")
  })
  
  ######## reset colors
  observeEvent(input$defaultColorTrace, {
    shinyjs::reset("lowColor_trace")
    shinyjs::reset("highColor_trace")
  })
  
  observeEvent(input$defaultColorFas, {
    shinyjs::reset("lowColor_fas")
    shinyjs::reset("highColor_fas")
  })
  
  ######## check if there is any "unknown" taxon in input matrix
  unkTaxa <- reactive({
    filein <- input$file1
    if(is.null(filein)){return()}
    
    # get list of all available taxon (from taxonID.list.fullrankID)
    allTaxa <- unlist(as.list(fread("data/taxonID.list.fullrankID",sep = "\t", select = "abbrName")))
    
    # get list of input taxa (from main input file)
    inputTaxa <- readLines(filein$datapath, n = 1)
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
  
  ######## enable "add taxa", "parse" & "upload additional files" button after uploading main file
  observeEvent(input$file1, ({
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
    filein <- input$file1
    if(is.null(filein)){return()}
    else{
      if(v1$parseInput == FALSE){return()}
      else{
        titleline <- readLines(filein$datapath, n=1) 
        
        # Create 0-row data frame which will be used to store data
        dat <- data.frame(x = numeric(0), y = numeric(0))   ### use for progess bar
        withProgress(message = 'Parsing input file', value = 0, {
          cmd <- paste("perl ", getwd(),"/data/getTaxonomyInfo.pl", 
                       #                 " -i ", getwd(),"/data/",input$file1,
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
                selected = "")
  })
  
  ####### GET list of all (super)taxa
  allTaxaList <- reactive({
    #  output$select = renderUI({
    filein <- input$file1
    if(is.null(filein)){return()}
    
    rankSelect = input$rankSelect
    if(rankSelect == ""){return()}
    ### load list of unsorted taxa
    Dt <- as.data.frame(read.table("data/taxonID.list.fullRankID", sep='\t',header=T))
    
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
    geneList <- preDataFiltered()
    geneList$geneID <- as.factor(geneList$geneID)
    
    out <- as.list(levels(geneList$geneID))
    out <- append("none",out)
    
    selectInput('geneHighlight','Highlight:',out,selected=out[1])
  })

  #### update highlightGeneUI based on double clicked dot
  observe({
    geneList <- preDataFiltered()
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
  
  ######## disable main input, genelist input and 3 initial questions
  observe({
    # use tabsetPanel 'id' argument to change tabs
    if (input$do > 0) {
      toggleState("file1")
      toggleState("geneList_selected")
      toggleState("sortGene")
      toggleState("newTaxaAsk")
      toggleState("parseAsk")
      #updateTabsetPanel(session, "tabs", selected = "Main profile")
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
    filein <- input$file1
    if(is.null(filein)){
      v$doPlot <- FALSE
      updateButton(session, "do", disabled = TRUE)
    }
  })
  
  ######## sorting supertaxa list based on chosen reference taxon
  sortedTaxaList <- reactive({
    if(v$doPlot == FALSE){return()}

    ### load list of unsorted taxa
    Dt <- as.data.frame(read.table("data/taxonID.list.fullRankID", sep='\t',header=T))
    
    ### reduce the number of columns in the iriginal data frame by removing duplicate columns
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

    ### load list of taxon name
    nameList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T,fill = TRUE))
    nameList$fullName <- as.character(nameList$fullName)

    ### input parameters
    rankSelect = input$rankSelect
    rankName = substr(rankSelect,4,nchar(rankSelect))   # get rank name from rankSelect
    #rankNr = 0 + as.numeric(substr(rankSelect,1,2))     # get rank number (number of column in unsorted taxa list - dataframe Dt)

    # get selected supertaxon ID
    taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
    superID <- as.integer(taxaList$ncbiID[taxaList$fullName == input$inSelect])

    ### sort taxa list
    ### first move all species that have the same ID of selected rank (level) to a new data frame
    ### and sort the rest of the data frame based on that rank.
    ### then move species that have different ID of selected rank (level), but have the same ID of the higher level
    ### repeat until reach to last column

    sortedDt <- data.frame()
    firstLine <- Dt[Dt[,rankName]==superID,][1,]  # get all taxo info for 1 representative
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
  
  ############### PARSING DATA FROM INPUT MATRIX:
  ############### get traceability scores for coressponding orthologs (if trace matrix is provided) (2)
  ############### get (super)taxa names (3)
  ############### calculate percentage of presence (4), max FAS (5) and mean TRACEEBILITY (6) if group input taxa list into higher taxonomy rank
  preDataFiltered <- reactive({
    ### (1) LOADING INPUT MATRIX (1)
    filein <- input$file1
    if(is.null(filein)){return()}
    
    # get rows need to be read
    nrHit <- input$stIndex + input$number - 1
    data <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char="",nrows=nrHit))
    
    # OR just list of gene from a separated input file
    listIn <- input$list
    if(input$geneList_selected == 'from file'){
      if(!is.null(listIn)){
        list <- as.data.frame(read.table(file=listIn$datapath, header=FALSE))
        dataOrig <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
        data <- dataOrig[dataOrig$geneID %in% list$V1,]
      }
    }
    
    if(input$sortGene == "No"){
      data$geneID <- factor(data$geneID, levels = data$geneID)  ######### keep user defined geneID order 
    }
    
    ### return preDataFiltered
    data
  })
  
  dataFiltered <- reactive({
    data <- preDataFiltered()
    # convert into paired columns
    mdData <- melt(data,id="geneID")
    
    # split value column into orthoID and fas
    splitDt <- (str_split_fixed(mdData$value, '#', 2))
    # then join them back to mdData
    mdData <- cbind(mdData,splitDt)
    # rename columns
    colnames(mdData) <- c("geneID","ncbiID","value","orthoID","fas")
    mdData <- mdData[,c("geneID","ncbiID","fas","orthoID")]
    
    ### (2) LOADING TRACEABILITY MATRIX and merge with input mdData (2) ###
    filein2 <- input$file2
    if(is.null(filein2)){
      mdDataTrace <- mdData[,c("geneID","ncbiID")]
      mdDataTrace$traceability <- 0
    } else {
      dataTrace <- as.data.frame(read.table(file=filein2$datapath, sep='\t',header=T,check.names=FALSE,comment.char="",nrows=nrHit))
      
      ## get subset of dataTrace if a list of genes is given
      if(input$geneList_selected == 'from file'){
        if(!is.null(listIn)){
          list <- as.data.frame(read.table(file=listIn$datapath, header=FALSE))
          dataOrig <- as.data.frame(read.table(file=filein2$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
          dataTrace <- dataOrig[dataOrig$geneID %in% list$V1,]
        }
      }
      
      mdDataTrace <- melt(dataTrace,id="geneID")
      colnames(mdDataTrace) <- c("geneID","ncbiID","traceability")
    }
    
    ### (3) GET SORTED TAXONOMY LIST (3) ###
    taxaList <- sortedTaxaList()
    
    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxaList,'supertaxon')
    
    # merge mdData, mdDataTrace and taxaList to get taxonomy info
    taxaMdData <- merge(mdData,taxaList,by='ncbiID')
    taxaMdData$fas <- as.numeric(as.character(taxaMdData$fas))
    taxaMdDataTrace <- merge(mdDataTrace,taxaList,by='ncbiID')  #################### FOR TRACEABILITY SCORES
    
    ### (4) calculate PERCENTAGE of PRESENT SPECIES (4) ###
    # get geneID and supertaxon
    geneIDsupertaxon <- subset(taxaMdData,select=c('geneID','supertaxon'))
    geneIDsupertaxon <- geneIDsupertaxon[!duplicated(geneIDsupertaxon), ] # remove duplicated rows
    
    # remove NA rows from taxaMdData
    taxaMdDataNoNA <- taxaMdData[!is.na(taxaMdData$fas),]
    
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
    
    ### (5) calculate MAX FAS for every supertaxon of each gene (5) ###
    maxFasDt <- aggregate(taxaMdDataNoNA[,"fas"],list(taxaMdDataNoNA$supertaxon,taxaMdDataNoNA$geneID),max)
    colnames(maxFasDt) <- c("supertaxon","geneID","maxFas")
    
    ### (6) calculate mean TRACEABILITY SCORES for each super taxon (6) ###
    meanTraceDt <- aggregate(taxaMdDataTrace[,"traceability"],list(taxaMdDataTrace$supertaxon,taxaMdDataTrace$geneID),mean)
    colnames(meanTraceDt) <- c("supertaxon","geneID","traceability")
    
    ### (5+6) & join mean traceability together with max fas scores into one df (5+6)
    scoreDf <- merge(maxFasDt,meanTraceDt, by=c("supertaxon","geneID"), all = TRUE)
    
    ### (4+5+6) add presSpec and maxFAS into taxaMdData (4+5+6)
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
  ######## this data set contain only supertaxa and their value (%present, max fas & mean tracebility) for each gene
  dataSupertaxa <- reactive({
    fullMdData <- dataFiltered()
    
    ### get representative orthoID that has max FAS for each supertaxon
    maxOrthoID <- fullMdData[,c('geneID','supertaxon','fas','maxFas','orthoID')]
    maxOrthoID <- subset(maxOrthoID,maxOrthoID$fas == maxOrthoID$maxFas)
    colnames(maxOrthoID) <- c('geneID','supertaxon','fas','maxFas','orthoID')
    maxOrthoID <- maxOrthoID[!is.na(maxOrthoID$orthoID),]
    maxOrthoID <- maxOrthoID[,c('geneID','supertaxon','orthoID')]
    maxOrthoID <- maxOrthoID[!duplicated(maxOrthoID[,1:2]), ]
    
    ### get data set for phyloprofile plotting (contains only supertaxa info)
    superDf <- subset(fullMdData,select=c('geneID','supertaxon','supertaxonID','maxFas','presSpec','category','traceability'))
    superDf <- superDf[!duplicated(superDf), ]
    superDfExt <- merge(superDf,maxOrthoID, by=c('geneID','supertaxon'),all.x=TRUE)
    superDfExt <- superDfExt[,c("geneID","supertaxon","supertaxonID","maxFas","presSpec","category","orthoID","traceability")]
    
    ### output
    names(superDfExt)[names(superDfExt)=="maxFas"] <- "fas"
    superDfExt
  })
  

  ######## get list of all sequence IDs for selectize input
  output$geneIn = renderUI({
    filein <- input$file1
    fileCustom <- input$customFile
    
    if(is.null(filein) & is.null(fileCustom)){return(selectInput('inSeq','',"all"))}
    if(v$doPlot == FALSE){return(selectInput('inSeq','',"all"))}
    else{
      if(is.null(fileCustom)){
        data <- as.data.frame(dataFiltered())
        data$geneID <- as.character(data$geneID)
        data$geneID <- as.factor(data$geneID)
        out <- as.list(levels(data$geneID))
        out <- append("all",out)
        selectInput('inSeq','',out,selected=out[1],multiple=TRUE)
       } else {
         customList <- as.data.frame(read.table(file=fileCustom$datapath, header=FALSE))
         customList$V1 <- as.factor(customList$V1)
         out <- as.list(levels(customList$V1))
         selectInput('inSeq','',out,selected=out,multiple=TRUE)
       }
    }
  })
  
  ######## get list of all taxa for selectize input
  output$taxaIn = renderUI({
    filein <- input$file1
    if(is.null(filein)){return(selectInput('inTaxa','Select (super)taxon/(super)taxa of interest:',"all"))}
    if(v$doPlot == FALSE){return(selectInput('inTaxa','Select (super)taxon/(super)taxa of interest:',"all"))}
    else{
      choice <- allTaxaList()
      choice$fullName <- as.factor(choice$fullName)
      
      out <- as.list(levels(choice$fullName))
      out <- append("all",out)
      
      selectInput('inTaxa','Select (super)taxon/(super)taxa of interest:',out,selected=out[1],multiple=TRUE)
    }
  })

  ######## heatmap data input
  dataHeat <- reactive({
    percent_cutoff <- input$percent
    fas_cutoff <- input$fas
    
    ### check input file
    filein <- input$file1
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
    
    # ### check input file
    # filein <- input$file1
    # if(is.null(filein)){return()}
    # #      data <- read.table(file=filein$datapath, sep='\t',header=T)
    # dataHeat <- dataSupertaxa()
    # 
    # # get selected supertaxon name
    # split <- strsplit(as.character(input$inSelect),"_")
    # inSelect <- as.character(split[[1]][1])

    ### replace insufficient values according to the thresholds by NA or 0; and replace FAS 0.0 by NA
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$presSpec < percent_cutoff] <- 0
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$fas < fas_cutoff] <- 0
    dataHeat$fas[dataHeat$supertaxon != inSelect & dataHeat$fas < fas_cutoff] <- NA
    # dataHeat$fas[dataHeat$fas == 0] <- NA
    
    dataHeat <- droplevels(dataHeat)  ### delete unused levels
    dataHeat
  })
  
  ########### create profile heatmap
  mainPlot <- function(){
        if (v$doPlot == FALSE) return()
        dataHeat <- dataHeat()
      
        ### plot format
        if(input$xAxis == "genes"){
          p = ggplot(dataHeat, aes(x = geneID, y = supertaxon)) +        ## global aes
            scale_fill_gradient(low = input$lowColor_trace, high = input$highColor_trace, na.value="gray95") +   ## fill color (traceability)
            geom_tile(aes(fill = traceability)) +    ## filled rect (traceability score)
            geom_point(aes(colour = fas, size = presSpec))  +    ## geom_point for circle illusion (FAS and presence/absence)
            scale_color_gradient(low = input$lowColor_fas,high = input$highColor_fas)#+       ## color of the corresponding aes (FAS)
          scale_size(range = c(0,3))             ## to tune the size of circles
          p = p + labs(y="Taxon")
          
          base_size <- 9
          p = p+geom_hline(yintercept=0.5,colour="dodgerblue4")
          p = p+geom_hline(yintercept=1.5,colour="dodgerblue4")
          p = p+theme(axis.text.x = element_text(angle=60,hjust=1,size=input$xSize),axis.text.y = element_text(size=input$ySize),
                      axis.title.x = element_text(size=input$xSize), axis.title.y = element_text(size=input$ySize),
                      legend.title=element_text(size=input$legendSize),legend.text=element_text(size=input$legendSize),legend.position = input$mainLegend)
        } else {
          p = ggplot(dataHeat, aes(y = geneID, x = supertaxon)) +        ## global aes
            scale_fill_gradient(low = input$lowColor_trace, high = input$highColor_trace, na.value="gray95") +   ## fill color (traceability)
            geom_tile(aes(fill = traceability)) +    ## filled rect (traceability score)
            geom_point(aes(colour = fas, size = presSpec))  +    ## geom_point for circle illusion (FAS and presence/absence
            scale_color_gradient(low = input$lowColor_fas,high = input$highColor_fas)#+       ## color of the corresponding aes
          scale_size(range = c(0,3))             ## to tune the size of circles
          p = p + labs(x="Taxon")
          
          
          base_size <- 9
          p = p+geom_vline(xintercept=0.5,colour="dodgerblue4")
          p = p+geom_vline(xintercept=1.5,colour="dodgerblue4")
          p = p+theme(axis.text.x = element_text(angle=60,hjust=1,size=input$xSize),axis.text.y = element_text(size=input$ySize),
                      axis.title.x = element_text(size=input$xSize), axis.title.y = element_text(size=input$ySize),
                      legend.title=element_text(size=input$legendSize),legend.text=element_text(size=input$legendSize),legend.position = input$mainLegend)
        }
        
        ### highlight taxon
        if(input$taxonHighlight != "none"){
          ## get selected highlight taxon ID
          taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
          taxonHighlight <- as.integer(taxaList$ncbiID[taxaList$fullName == input$taxonHighlight])
          
          ## get taxonID together with it sorted index
          highlightTaxon <- toString(dataHeat[dataHeat$supertaxonID == taxonHighlight,2][1])
          ## get index
          selectedIndex = as.numeric(as.character(substr(highlightTaxon,2,4)))
          ## draw a rect to highlight this taxon's column
          if(input$xAxis == "taxa"){
            rect <- data.frame(xmin=selectedIndex-0.5, xmax=selectedIndex+0.5, ymin=-Inf, ymax=Inf)
            p = p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                              color="yellow",
                              alpha=0.3,
                              inherit.aes = FALSE)
          } else {
            rect <- data.frame(ymin=selectedIndex-0.5, ymax=selectedIndex+0.5, xmin=-Inf, xmax=Inf)
            p = p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                              color="yellow",
                              alpha=0.3,
                              inherit.aes = FALSE)
          }
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
            p = p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                              color="yellow",
                              alpha=0.3,
                              inherit.aes = FALSE)
          } else {
            rect <- data.frame(xmin=selectedIndex-0.5, xmax=selectedIndex+0.5, ymin=-Inf, ymax=Inf)
            p = p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                              color="yellow",
                              alpha=0.3,
                              inherit.aes = FALSE)
          }
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
  
  output$mainAxis <- renderPlot(bg="transparent",{
    if(input$autoUpdate == FALSE){
      # Add dependency on the update button (only update when button is clicked)
      input$updateBtn
      isolate({
        #list <- allTaxaList()
        p <- mainPlot()
        g <- ggplotGrob(p)
        
        if(input$mainXAxisGuide == TRUE & input$mainYAxisGuide == FALSE){
          s <- gtable_filter(g, 'axis-b', trim=F)  ### filter to get x-axis
        } else if (input$mainXAxisGuide == FALSE & input$mainYAxisGuide == TRUE){
          s <- gtable_filter(g, 'axis-l', trim=F)  ### filter to get y-axis
        } else if (input$mainXAxisGuide == TRUE & input$mainYAxisGuide == TRUE){
          s <- gtable_filter(g, 'axis-b|axis-l', trim=F)  ### filter to get x-axis and y-axis
        }
        
        # draw axis(es)
        grid.draw(s)
      })
    } else {
      p <- mainPlot()
      g <- ggplotGrob(p)
      
      if(input$mainXAxisGuide == TRUE & input$mainYAxisGuide == FALSE){
        s <- gtable_filter(g, 'axis-b', trim=F)  ### filter to get x-axis
      } else if (input$mainXAxisGuide == FALSE & input$mainYAxisGuide == TRUE){
        s <- gtable_filter(g, 'axis-l', trim=F)  ### filter to get y-axis
      } else if (input$mainXAxisGuide == TRUE & input$mainYAxisGuide == TRUE){
        s <- gtable_filter(g, 'axis-b|axis-l', trim=F)  ### filter to get x-axis and y-axis
      }
      grid.draw(s)
    }
  })
  
  output$mainAxisRender <- renderUI({
    if(input$autoUpdate == FALSE){
      # Add dependency on the update button (only update when button is clicked)
      input$updateBtn
      isolate({
        plotOutput("mainAxis", width=input$width, height=input$height)
      })
    } else{
      plotOutput("mainAxis", width=input$width, height=input$height)
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
      # get FAS, percentage of present species and traceability score
      FAS <- dataHeat$fas[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      Percent <- dataHeat$presSpec[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      Trace <- dataHeat$traceability[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      # get ortholog ID
      orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      
      ### get list of all geneID that have the same ortholog
      geneMatch <- dataHeat$geneID[dataHeat$orthoID == toString(orthoID)]
      geneMatch <- geneMatch[!is.na(geneMatch)]
      # list of all available geneID
      geneList <- preDataFiltered()
      geneList$geneID <- as.factor(geneList$geneID)
      allGenes <- as.list(levels(geneList$geneID))
      # get index of all matched genes (genes have the same ortholog)
      pos <- which(allGenes %in% geneMatch)
      pos <- paste(pos, collapse=',')

      ### return info of clicked point      
      if(is.na(as.numeric(Percent))){return()}
      else{
        info <- c(geneID,as.character(orthoID),as.character(spec),round(as.numeric(FAS),2),round(as.numeric(Percent),2),round(as.numeric(Trace),2),pos)
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
  
  
  ############################################################# 
  ################# PLOT SELECTED SEQUENCES ###################
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
    filein <- input$file1
    if(is.null(filein)){v2$doPlotCustom <- FALSE}
  })
  
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
        p = ggplot(dataHeat, aes(y = geneID, x = supertaxon)) +        ## global aes
          scale_fill_gradient(low = input$lowColor_trace, high = input$highColor_trace, na.value="gray95") +   ## fill color (traceability)
          geom_tile(aes(fill = traceability)) +# + scale_fill_gradient(low="gray95", high="red")) +    ## filled rect (traceability score)
          geom_point(aes(colour = fas, size = presSpec))  +    ## geom_point for circle illusion (FAS and presence/absence)
          scale_color_gradient(low = input$lowColor_fas,high = input$highColor_fas)#+       ## color of the corresponding aes (FAS)
        scale_size(range = c(0,3))             ## to tune the size of circles
        p = p + labs(x="Taxon")
        
        base_size <- 9
        p = p+geom_vline(xintercept=0.5,colour="dodgerblue4")
        p = p+geom_vline(xintercept=1.5,colour="dodgerblue4")
        p = p+theme(axis.text.x = element_text(angle=60,hjust=1,size=input$xSizeSelect),axis.text.y = element_text(size=input$ySizeSelect),
                    axis.title.x = element_text(size=input$xSizeSelect), axis.title.y = element_text(size=input$ySizeSelect),
                    legend.title=element_text(size=input$legendSizeSelect),legend.text=element_text(size=input$legendSizeSelect),legend.position=input$selectedLegend)
      } else {
        p = ggplot(dataHeat, aes(x = geneID, y = supertaxon)) +        ## global aes
          scale_fill_gradient(low = input$lowColor_trace, high = input$highColor_trace, na.value="gray95") +   ## fill color (traceability)
          geom_tile(aes(fill = traceability)) +    ## filled rect (traceability score)
          geom_point(aes(colour = fas, size = presSpec))  +    ## geom_point for circle illusion (FAS and presence/absence)
          scale_color_gradient(low = input$lowColor_fas,high = input$highColor_fas)#+       ## color of the corresponding aes (FAS)
        scale_size(range = c(0,3))             ## to tune the size of circles
        p = p + labs(y="Taxon")
        
        base_size <- 9
        p = p+geom_vline(xintercept=0.5,colour="dodgerblue4")
        p = p+geom_vline(xintercept=1.5,colour="dodgerblue4")
        p = p+theme(axis.text.x = element_text(angle=60,hjust=1,size=input$xSizeSelect),axis.text.y = element_text(size=input$ySizeSelect),
                    axis.title.x = element_text(size=input$xSizeSelect), axis.title.y = element_text(size=input$ySizeSelect),
                    legend.title=element_text(size=input$legendSizeSelect),legend.text=element_text(size=input$legendSizeSelect),legend.position=input$selectedLegend)
      }
      
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
    if (v$doPlot == FALSE) return()
    
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
      # get FAS, percentage of present species and traceability score
      FAS <- dataHeat$fas[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      Percent <- dataHeat$presSpec[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      Trace <- dataHeat$traceability[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      # get ortholog ID
      orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      
      if(is.na(as.numeric(Percent))){return()}
      else{
        info <- c(geneID,as.character(orthoID),as.character(spec),round(as.numeric(FAS),2),round(as.numeric(Percent),2),round(as.numeric(Trace),2))
      }
    }
  })

  
  ############################################################# 
  ################### SHOW CLICKED POINT INFO #################
  #############################################################
  
  ######## show info into "point's info" box
  output$pointInfo <- renderText({
    ### check input
    if (v$doPlot == FALSE) return()
    
    ##### GET INFO BASED ON CURRENT TAB
    if(input$tabs == 'Main profile'){
      info <- mainPointInfo()  # info = groupID,orthoID,supertaxon,maxFAS,%spec,trace
    } else if(input$tabs=='Customized profile'){
      info <- selectedPointInfo()
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
      c <- toString(paste("maxFas:",info[4],"; %spec:",info[5]))
      d <- toString(paste("traceability:",info[6]))
      paste(a,b,c,d,sep="\n")
    }
  })
  
  
  ############################################################# 
  ##################### DETAILED PLOT #########################
  #############################################################
  
  ######## data for detailed FAS plot
  detailPlotDt <- reactive({
    if (v$doPlot == FALSE) return()
    
    ##### GET INFO BASED ON CURRENT TAB
    if(input$tabs == 'Main profile'){
      info <- mainPointInfo()  # info = groupID,orthoID,supertaxon,maxFAS,%spec,trace
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
  
  ######## render detailed FAS plot
  output$detailPlot <- renderPlot({
    if (v$doPlot == FALSE) return()
    
    selDf <- detailPlotDt()
    selDf$x_label <- paste(selDf$orthoID,"@",selDf$fullName,sep = "")
    
    gp = ggplot(selDf, aes(y=fas,x=x_label)) +
      geom_bar(colour="steelblue", fill="steelblue", stat="identity") +
      coord_flip() +
      labs(x="") #+
    #geom_text(aes(label=fas), vjust=3)
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
    
    ### get pair of sequence IDs & FAS
    seedID <- toString(selDf$geneID[1])
    orthoID <- toString(allOrthoID[corX])
    fas <- toString(selDf$fas[selDf$orthoID==orthoID])
    
    ### return info
    if(orthoID != "NA"){
      info <- c(seedID,orthoID,fas)
    }
  })
  
  ### SHOW info when clicking on detailed plot
  output$detailClick <- renderText({
    info <- pointInfoDetail() # info = seedID, orthoID, FAS
    if(is.null(info)){return()}
    else{
      a <- paste0("seedID = ",info[1])
      b <- paste0("hitID = ",info[2])
      c <- paste0("FAS = ",info[3])
      paste(a,b,c,sep="\n")
    }
  })
  
  ######## FASTA sequence
  output$fasta <- renderText({
    if(v$doPlot == FALSE){return()}
    
    info <- pointInfoDetail() # info = seedID, orthoID, FAS
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
    filein <- input$file1
    if(is.null(filein)){v3$doPlot3 <- FALSE}
  })
  
  ######## create domain plot
  archiPlot <- function(){
    if (v3$doPlot3 == FALSE) return()
    
    ### info
    info <- pointInfoDetail() # info = seedID, orthoID, FAS
    group <- as.character(info[1])
    ortho <- as.character(info[2])
    fas <- as.character(info[3])
    
    ### parse domain file
    filein3 <- input$file3
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
        orthoDf <- orthoDf[order(orthoDf$feature), ]
        seedDf$feature <- factor(seedDf$feature, levels=c(orthoDf$feature,seedDf$feature[seedDf$feature != orthoDf$feature]))
      } else {
        seedDf <- seedDf[order(seedDf$feature), ]
        orthoDf$feature <- factor(orthoDf$feature, levels=c(seedDf$feature,orthoDf$feature[orthoDf$feature != seedDf$feature]))
      }
      
      ### plotting
      sep = ":"
      if(!is.null(input$oneSeqFasta)){sep="|"}
      plot_ortho <- plotting(orthoDf,ortho,fas,sep,input$labelArchiSize,input$titleArchiSize,min(subDomainDf$start),max(subDomainDf$end))
      plot_seed <- plotting(seedDf,seed,fas,sep,input$labelArchiSize,input$titleArchiSize,min(subDomainDf$start),max(subDomainDf$end))
      
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
      domainIN <- unlist(strsplit(toString(input$file1),","))
      fileName <- toString(domainIN[1])
      msg <- paste0(
        "<p><span style=\"color: #ff0000;\"><strong>No information about domain architecture! Please check:</strong></span></p>
        <ul style=\"list-style-type: square;\">
        <li>if you selected any sequence in the Detailed plot?</li>
        <li>if you uploaded the domain file using Upload additional file(s) option? (see input example in data/lca.FASmatrix.mDomains)</li>
        </ul>"
      )
      HTML(msg)
    } else {
      plotOutput("archiPlot",height = input$archiHeight, width = input$archiWidth)
    }
  })
  
  ######## download architecture plot ***** something strange with archiPlot()
  output$archiDownload <- downloadHandler(
    filename = "domains.svg",
    content = function(file) {
      g <- archiPlot()
      grid.draw(g)
      ggsave(file, plot = g, width = input$archiWidth*0.056458333, height = input$archiHeight*0.056458333, units="cm", dpi=300)#, device = "svg")
    }
  )
  
  
  
  ############################################################# 
  ############### FILTERED DATA FOR DOWNLOADING ###############
  #############################################################
  
  ######## filtered data for downloading
  downloadData <- reactive({
    ### check input
    if (v$doPlot == FALSE) return()
    
    dataOut <- dataFiltered()
    dataOut <- as.data.frame(dataOut[dataOut$presSpec > 0,])
    dataOut <- dataOut[!is.na(dataOut$geneID),]
    
    dataOut <- as.data.frame(dataOut[dataOut$presSpec >= input$percent & dataOut$fas >= input$fas,])
    
    dataOut <- dataOut[,c("geneID","orthoID","fullName","ncbiID","supertaxon","fas","numberSpec","presSpec","traceability")]
    dataOut <- dataOut[order(dataOut$geneID,dataOut$supertaxon),]
    dataOut <- dataOut[complete.cases(dataOut),]
    
    dataOut$geneID <- as.character(dataOut$geneID)
    dataOut$fullName <- as.character(dataOut$fullName)
    dataOut$ncbiID <- substr(dataOut$ncbiID,5,nchar(as.character(dataOut$ncbiID)))
    dataOut$supertaxon <- substr(dataOut$supertaxon,6,nchar(as.character(dataOut$supertaxon)))
    dataOut$fas <- as.character(dataOut$fas)
    dataOut$numberSpec <- as.integer(dataOut$numberSpec)
    dataOut$presSpec <- as.numeric(dataOut$presSpec)
    
    names(dataOut)[names(dataOut)=="presSpec"] <- "%Spec"
    names(dataOut)[names(dataOut)=="numberSpec"] <- "totalSpec"
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
  output$dis <- renderDataTable({
    if(v$doPlot == FALSE){return()}
    #data <- allTaxaList()
    #data <- sortedTaxaList()
    #data <- preDataFiltered()
    #data <- dataFiltered()
    #data <- dataSupertaxa()
    #data <- dataHeat()
    #data <- detailPlotDt()
    data <- downloadData()
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
      <p><em>More detail? Pleas take a look at the example file in /data/lca.FASmatrix :)</em></p>
      <p>&nbsp;</p>
      <h1 style="color: #5e9ca0;">Additional files</h1>
      <p>2 additional input files can be provided are traceability score matrix and feature domain position list.</p>
      <p><strong>Traceability score matrix</strong> must have the same first row (beginning with "geneID" and followed by list of taxa) and first column (list of genes). <span style="color: #ff0000;"><strong>IMPORTANT</strong>: the <span style="text-decoration: underline;">amount</span> and <span style="text-decoration: underline;">order</span> of genes between 2 input matrixes have to be exactly the same!!</span>&nbsp;</p>
      <p><strong>Feature domain position list</strong>&nbsp;has 6 columns separated by tab: (1) pairID = groupID#searchProt_ID#seedID, (2) searchProt_ID, (3) feature name (pfam domain, smart domain,etc.), (4) start position, (5) end position, (6) weight value (only available for seed protein)</p>
      <p><em>Pleas take a look at the example files &nbsp;lca.Tracematrix and lca.FASmatrix.mDomains in data/ folder for more details :)</em></p>
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
    # filein <- input$file1
    # print(toString(filein))
    # filePath <- toString(filein)
    # fileName <- unlist(strsplit(toString(input$file1),","))
    # name <- toString(fileName[1])
    # fullPath <- paste0("data/",name,".mDomains")
    # print(fullPath)
    # print(input$plot_dblclick$x)
  })
})