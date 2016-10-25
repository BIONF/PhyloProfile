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

options(shiny.maxRequestSize=30*1024^2)  ## size limit for input 30mb
shinyServer(function(input, output, session) {
  
  ##### check if data is loaded and "plot" button is clicked
  v <- reactiveValues(doPlot = FALSE)
  observeEvent(input$do, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    v$doPlot <- input$do
    filein <- input$file1
    if(is.null(filein)){v$doPlot <- FALSE}
  })
  
  ### check if data is loaded and "parse" button (get info from input) is clicked
  v1 <- reactiveValues(parseInput = FALSE)
  observeEvent(input$parse, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    v1$parseInput <- input$parse
    filein <- input$file1
    if(is.null(filein)){v1$parseInput <- FALSE}
  })
  
  ######### create taxonID.list.fullRankID and taxonNamesReduced.txt from input file (if necessary)
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
                       " -o ", getwd(),"/data",
                       sep='')
          system(cmd) 
        })
      }
    }
  })
  
  ######## list of taxonomy ranks for plotting
  output$rankSelect = renderUI({
    selectInput("rankSelect", label = "Select taxonomy rank:",
                choices = list("Strain"="05_strain","Species" = "06_species","Genus" = "10_genus", "Family" = "14_family", "Order" = "19_order", "Class" = "23_class",
                               "Phylum" = "26_phylum", "Kingdom" = "28_kingdom", "Superkingdom" = "29_superkingdom","unselected"=""), 
                selected = "")
  })
  
  ####### get list of all (super)taxa
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
    rankNr = 0 + as.numeric(substr(rankSelect,1,2))     # get rank number (number of column in unsorted taxa list - dataframe Dt)
    
    choice <- as.data.frame
    choice <- rbind(Dt[rankNr])
    colnames(choice) <- "ncbiID"
    choice <- merge(choice,nameList,by="ncbiID",all = FALSE)
  })
  
  ### then output list of species onto UI
  output$select = renderUI({
    choice <- allTaxaList()
    #    choice$fullName <- paste0(choice$fullName,"_",choice$ncbiID)
    choice$fullName <- as.factor(choice$fullName)
    selectInput('inSelect','Choose (super)taxon of interest:',as.list(levels(choice$fullName)),levels(choice$fullName)[1])
  })
  
  output$highlight = renderUI({
    #    if(input$xAxis == "taxa"){
    choice <- allTaxaList()
    #      choice$fullName <- paste0(choice$fullName,"_",choice$ncbiID)
    choice$fullName <- as.factor(choice$fullName)
    selectInput('inHighlight','Select (super)taxon to highlight:',as.list(levels(choice$fullName)),levels(choice$fullName)[1])
  })
  
  ######## sorting supertaxa list
  sortedTaxaList <- reactive({
    if(v$doPlot == FALSE){return()}
    
    ### load list of unsorted taxa
    Dt <- as.data.frame(read.table("data/taxonID.list.fullRankID", sep='\t',header=T))
    
    ### load list of taxon name
    nameList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T,fill = TRUE))
    nameList$fullName <- as.character(nameList$fullName)
    
    ### input parameters
    rankSelect = input$rankSelect
    rankName = substr(rankSelect,4,nchar(rankSelect))   # get rank name from rankSelect
    rankNr = 0 + as.numeric(substr(rankSelect,1,2))     # get rank number (number of column in unsorted taxa list - dataframe Dt)
    
    # get selected supertaxon ID
    taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
    superID <- as.integer(taxaList$ncbiID[taxaList$fullName == input$inSelect & taxaList$rank == rankName])
    
    ################ sort taxa list using info from pruned common tree
    ### read full common tree
    tree <- read.tree("data/commontree.phy.ids_mod")
    ### prune common tree
    filein <- input$file1
    list <- readLines(filein$datapath, n=1) # get title line of input matrix (which contains list of all taxa IDs)
    listMod <- gsub("ncbi","",list)
    listMod <- gsub("geneID\t","",listMod)
    ids <- unlist(strsplit(listMod,"\t"))
    pruned.tree<-drop.tip(tree, setdiff(tree$tip.label, ids));
    ### get representative ID for rooting (one member of selected supertaxon)
    repID <- Dt[Dt[,rankNr]==superID,][,3][1]
    root <- toString(repID)
    ### reroot the pruned tree
    unrootTree <- unroot(pruned.tree)
    rerootTree <- root(unrootTree,root,resolve.root = TRUE)  ## reroot the tree based on selected ID from input$inSelect
    ### get sorting taxa
    newTree <- read.tree(text=write.tree(rerootTree))
    orderedTaxa <- rev(newTree$tip.label)
    # orderedTaxa
    # length(orderedTaxa)
    # orderedTaxa[2]
    
    ### now sort dataframe Dt
    sortedDt <- data.frame()
    for(i in 1:length(orderedTaxa)){
      subDt <- Dt[Dt[,"ncbiID"]==orderedTaxa[i],]
      #  print(subDt)
      sortedDt <- rbind(sortedDt,subDt)
    }
    # sortedDt
    
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
  
  ######## parsing data from input matrix
  dataFiltered <- reactive({
    ##### matrix input
    filein <- input$file1
    if(is.null(filein)){return()}
    data <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
    # convert into paired columns
    mdData <- melt(data,id="geneID")
    
    # split value column into orthoID and fas
    splitDt <- (str_split_fixed(mdData$value, '#', 2))
    # then join them back to mdData
    mdData <- cbind(mdData,splitDt)
    # rename columns
    colnames(mdData) <- c("geneID","ncbiID","value","orthoID","fas")
    mdData <- mdData[,c("geneID","ncbiID","fas","orthoID")]
    
    #################### traceability matrix input
#    dataTrace <- as.data.frame(read.table("data/lca.Tracematrix", sep='\t',header=T,check.names=FALSE,comment.char = ""))
    filein2 <- input$file2
    if(is.null(filein2)){
      mdDataTrace <- mdData[,c("geneID","ncbiID")]
      mdDataTrace$traceability <- 0
    } else {
      dataTrace <- as.data.frame(read.table(file=filein2$datapath, sep='\t',header=T,check.names=FALSE,comment.char=""))
      mdDataTrace <- melt(dataTrace,id="geneID")
      colnames(mdDataTrace) <- c("geneID","ncbiID","traceability")
    }
    ##### taxonomy file input
    #    taxaList <- as.data.frame(read.table("data/taxonomyList.txt", sep='\t',header=T))
    taxaList <- sortedTaxaList()
    
    # get frequency of all supertaxa
    taxaCount <- plyr::count(taxaList,'supertaxon')
    
    ### merge mdData, mdDataTrace and taxaList to get taxonomy info
    taxaMdData <- merge(mdData,taxaList,by='ncbiID')
    taxaMdData$fas <- as.numeric(as.character(taxaMdData$fas))
    
    taxaMdDataTrace <- merge(mdDataTrace,taxaList,by='ncbiID')  #################### FOR TRACEABILITY SCORES
    
    ############## calculate percent present species ##############
    ### get geneID and supertaxon
    geneIDsupertaxon <- subset(taxaMdData,select=c('geneID','supertaxon'))
    geneIDsupertaxon <- geneIDsupertaxon[!duplicated(geneIDsupertaxon), ] # remove duplicated rows
    
    ### remove NA rows from taxaMdData
    taxaMdDataNoNA <- taxaMdData[!is.na(taxaMdData$fas),]
    
    ### count present frequency of supertaxon for each gene
    geneSupertaxonCount <- plyr::count(taxaMdDataNoNA,c('geneID','supertaxon'))
    
    ### merge with taxaCount to get total number of species of each supertaxon and calculate presSpec
    presSpecDt <- merge(geneSupertaxonCount,taxaCount,by='supertaxon')
    presSpecDt$presSpec <- presSpecDt$freq.x/presSpecDt$freq.y
    presSpecDt <- presSpecDt[order(presSpecDt$geneID),]
    presSpecDt <- presSpecDt[,c("geneID","supertaxon","presSpec")]
    
    ### add absent supertaxon into presSpecDt
    finalPresSpecDt <- merge(presSpecDt,geneIDsupertaxon,by=c('geneID','supertaxon'),all.y = TRUE)
    finalPresSpecDt$presSpec[is.na(finalPresSpecDt$presSpec)] <- 0
    
    ############## calculate max FAS for every supertaxon of each gene ##############
    maxFasDt <- aggregate(taxaMdDataNoNA[,"fas"],list(taxaMdDataNoNA$supertaxon,taxaMdDataNoNA$geneID),max)
    colnames(maxFasDt) <- c("supertaxon","geneID","maxFas")
    
    ############## calculate mean TRACEABILITY SCORES for each super taxon
    meanTraceDt <- aggregate(taxaMdDataTrace[,"traceability"],list(taxaMdDataTrace$supertaxon,taxaMdDataTrace$geneID),mean)
    colnames(meanTraceDt) <- c("supertaxon","geneID","traceability")
    
    ############## & join mean traceability together with max fas scores into one df
    scoreDf <- merge(maxFasDt,meanTraceDt, by=c("supertaxon","geneID"), all = TRUE)

    ############## add presSpec and maxFAS into taxaMdData ##############
    presMdData <- merge(taxaMdData,finalPresSpecDt,by=c('geneID','supertaxon'),all.x = TRUE)
#    fullMdData <- merge(presMdData,maxFasDt,by=c('geneID','supertaxon'), all.x = TRUE)
    fullMdData <- merge(presMdData,scoreDf,by=c('geneID','supertaxon'), all.x = TRUE)
    fullMdData <- merge(fullMdData,taxaCount,by=('supertaxon'), all.x = TRUE)
    
#    names(fullMdData)[names(fullMdData)=="fas.x"] <- "fas"
#    names(fullMdData)[names(fullMdData)=="fas.y"] <- "maxFas"
    names(fullMdData)[names(fullMdData)=="freq"] <- "numberSpec"
    
    fullMdData$fullName <- as.vector(fullMdData$fullName)
    names(fullMdData)[names(fullMdData)=="orthoID.x"] <- "orthoID"
    fullMdData ### parsed input data frame !!!
  })
  
  ### reduce data from species level to supertaxa level
  ### this data set contain only supertaxa and their value (%present and max fas) for each gene
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
  
  ###### get list of all sequence IDs for selectize input
  output$geneIn = renderUI({
    #  output$select = renderUI({
    filein <- input$file1
    if(is.null(filein)){return(selectInput('inSeq','Select sequence IDs of interest:',"all"))}
    if(v$doPlot == FALSE){return(selectInput('inSeq','Select sequence IDs of interest:',"all"))}
    else{
      data <- as.data.frame(dataFiltered())
      data$geneID <- as.character(data$geneID)
      data$geneID <- as.factor(data$geneID)
      out <- as.list(levels(data$geneID))
      out <- append("all",out)
      
      selectInput('inSeq','Select sequence IDs of interest:',out,selected=out[1],multiple=TRUE)
    }
  })
  
  ### text output side panel
  output$totalRows <- renderText({
    if(v$doPlot == FALSE){return()}
    else {
      data <- dataSupertaxa()
      c <- length(unique(data$geneID))
      c("(max =",c,")")
    }
  })
  output$start <- renderText({
    c("start at:")
  })
  
  ### heatmap data input
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
    
    ### replace insufficient values according to the thresholds by NA or 0; and replace FAS 0.0 by NA
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$presSpec < percent_cutoff] <- 0
    dataHeat$presSpec[dataHeat$supertaxon != inSelect & dataHeat$fas < fas_cutoff] <- 0
    dataHeat$fas[dataHeat$supertaxon != inSelect & dataHeat$fas < fas_cutoff] <- NA
    dataHeat$fas[dataHeat$fas == 0] <- NA
    
    dataHeat <- droplevels(dataHeat)  ### delete unused levels
    dataHeat
    # if(input$inSeq != "all"){
    #     dataHeat <- dataHeat[dataHeat$geneID == input$inSeq,]
    # } else {
    #   dataHeat
    # }
  })
  
  ########### plot heatmap
  output$plot2 <- renderPlot(
    {
      if (v$doPlot == FALSE) return()
      
      # ### check input file
      # filein <- input$file1
      # if(is.null(filein)){return()}
      filein2 <- input$file2
      
      dataHeat <- dataHeat()
      ### plotting
      if(input$xAxis == "genes"){
        p = ggplot(dataHeat, aes(x = geneID, y = supertaxon)) +        ## global aes
          scale_fill_gradient(low = "gray95", high = "khaki", na.value="gray95") +   ## fill color
          geom_tile(aes(fill = traceability)) +    ## filled rect (traceability score)
          geom_point(aes(colour = fas, size = presSpec))  +    ## geom_point for circle illusion (FAS and presence/absence)
          scale_color_gradient(low = "darkorange",high = "steelblue")#+       ## color of the corresponding aes
        scale_size(range = c(0,3))             ## to tune the size of circles
        
        base_size <- 9
        p = p+geom_hline(yintercept=0.5,colour="dodgerblue4")
        p = p+geom_hline(yintercept=1.5,colour="dodgerblue4")
        p = p+theme(axis.text.x = element_text(angle=60,hjust=1))
      } else {
        p = ggplot(dataHeat, aes(y = geneID, x = supertaxon)) +        ## global aes
          scale_fill_gradient(low = "gray95", high = "khaki", na.value="gray95") +   ## fill color
          geom_tile(aes(fill = traceability)) +# + scale_fill_gradient(low="gray95", high="red")) +    ## filled rect (traceability score)
          geom_point(aes(colour = fas, size = presSpec))  +    ## geom_point for circle illusion (FAS and presence/absence
          scale_color_gradient(low = "darkorange",high = "steelblue")#+       ## color of the corresponding aes
        scale_size(range = c(0,3))             ## to tune the size of circles
        
        base_size <- 9
        p = p+geom_vline(xintercept=0.5,colour="dodgerblue4")
        p = p+geom_vline(xintercept=1.5,colour="dodgerblue4")
        p = p+theme(axis.text.x = element_text(angle=60,hjust=1))
      }
     
      ## get selected highlight taxon ID
      taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
      inHighlight <- as.integer(taxaList$ncbiID[taxaList$fullName == input$inHighlight])
      
      ## get taxonID together with it sorted index
      highlightTaxon <- toString(dataHeat[dataHeat$supertaxonID == inHighlight,2][1])
      ## get index
      selectedIndex = as.numeric(as.character(substr(highlightTaxon,2,4)))
      ## draw a rect to highlight this taxon's column
      if(input$xAxis == "taxa"){
        rect <- data.frame(xmin=selectedIndex-0.5, xmax=selectedIndex+0.5, ymin=-Inf, ymax=Inf)
        p = p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                          color="yellow",
                          alpha=0.3,
                          inherit.aes = FALSE)
        p
      } else {
        rect <- data.frame(ymin=selectedIndex-0.5, ymax=selectedIndex+0.5, xmin=-Inf, xmax=Inf)
        p = p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                          color="yellow",
                          alpha=0.3,
                          inherit.aes = FALSE)
        
        p
      }
    })
  
  ### show beschreibung file if no plot present
  output$plot.ui <- renderUI({
    if(v$doPlot == FALSE){
      if(!file.exists("www/beschreibung.jpg")){return(paste("WARNING: Cannot load \"beschreibung.jpg\" file in www folder!"))}
      else{return (img(src="beschreibung.jpg", align = "left", height=600, width=800))}
    }
    
    plotOutput("plot2",width=input$width,height = input$height,
               click = "plot_click",
               hover = hoverOpts(
                 id = "plot_hover",
                 delay = input$hover_delay,
                 delayType = input$hover_policy,
                 nullOutside = input$hover_null_outside
               )
    )
  })
  
  ### download plot
  output$plotDownload <- downloadHandler(
    filename = function() {c("plot.png")}, 
    content = function(file) {
      png(file, width = input$width, height = input$height)
      
      dataHeat <- dataHeat()
      ### plotting
      if(input$xAxis == "genes"){
        p = ggplot(dataHeat, aes(x = geneID, y = supertaxon)) +        ## global aes
          scale_fill_gradient(low = "gray95", high = "khaki", na.value="gray95") +   ## fill color
          geom_tile(aes(fill = traceability)) +    ## filled rect (traceability score)
          geom_point(aes(colour = fas, size = presSpec))  +    ## geom_point for circle illusion (FAS and presence/absence)
          scale_color_gradient(low = "darkorange",high = "steelblue")#+       ## color of the corresponding aes
        scale_size(range = c(0,3))             ## to tune the size of circles
        
        base_size <- 9
        p = p+geom_hline(yintercept=0.5,colour="dodgerblue4")
        p = p+geom_hline(yintercept=1.5,colour="dodgerblue4")
        p = p+theme(axis.text.x = element_text(angle=60,hjust=1))
        
        print(p)
        dev.off()
      } else {
        p = ggplot(dataHeat, aes(y = geneID, x = supertaxon)) +        ## global aes
          scale_fill_gradient(low = "gray95", high = "khaki", na.value="gray95") +   ## fill color
          geom_tile(aes(fill = traceability)) +    ## filled rect (traceability score)
          geom_point(aes(colour = fas, size = presSpec))  +    ## geom_point for circle illusion (FAS and presence/absence
          scale_color_gradient(low = "darkorange",high = "steelblue")#+       ## color of the corresponding aes
        scale_size(range = c(0,3))             ## to tune the size of circles
        
        base_size <- 9
        p = p+geom_vline(xintercept=0.5,colour="dodgerblue4")
        p = p+geom_vline(xintercept=1.5,colour="dodgerblue4")
        p = p+theme(axis.text.x = element_text(angle=60,hjust=1))
        
        print(p)
        dev.off()
      }
    }
  )
  
  ######## plot selected sequence
  v2 <- reactiveValues(doPlot2 = FALSE)
  observeEvent(input$do2, {
    # 0 will be coerced to FALSE
    # 1+ will be coerced to TRUE
    v2$doPlot2 <- input$do2
    filein <- input$file1
    if(is.null(filein)){v2$doPlot2 <- FALSE}
  })
  
  output$selectedPlot <- renderPlot({
    if (v2$doPlot2 == FALSE) return()
    if(input$inSeq == "all") {return()}
    else{
      dataHeat <- dataHeat()
      dataHeat <- dataHeat[dataHeat$geneID == input$inSeq,]
      
      ### plotting
      p = ggplot(dataHeat, aes(y = geneID, x = supertaxon)) +        ## global aes
        scale_fill_gradient(low = "gray95", high = "khaki", na.value="gray95") +
        geom_tile(aes(fill = traceability)) +    ## filled rect (traceability score)
        geom_point(aes(colour = fas, size = presSpec))  +    ## geom_point for circle illusion
        scale_color_gradient(low = "darkorange",high = "steelblue")#+       ## color of the corresponding aes
      scale_size(range = c(0,3))             ## to tune the size of circles
      
      base_size <- 9
      p = p+geom_vline(xintercept=0.5,colour="dodgerblue4")
      p = p+geom_vline(xintercept=1.5,colour="dodgerblue4")
      p = p+theme(axis.text.x = element_text(angle=60,hjust=1))
      p
    }
  })
  
  output$selectedPlot.ui <- renderUI({
    if(input$inSeq == "all"){return()}
    plotOutput("selectedPlot",width=input$selectedWidth, height = input$selectedHeight)
  })
  
  ######## get click info and plot detailed chart
  #### get info clicked point
  pointInfo <- reactive({
    ### check input
    if (v$doPlot == FALSE) return()
    
    # get selected supertaxon name
    taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
    rankSelect = input$rankSelect
    rankName = substr(rankSelect,4,nchar(rankSelect))
    inSelect <- as.numeric(taxaList$ncbiID[taxaList$fullName == input$inSelect & taxaList$rank == rankName])
    
    dataHeat <- dataHeat()
    
    ### get values
    if (is.null(input$plot_click$x)) return()
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
      genes <- as.matrix(dataHeat[dataHeat$supertaxonID == inSelect & !is.na(dataHeat$presSpec),])
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
        #substr(spec,6,nchar(as.character(spec)))
      }
    }
  })
  
  ### show info
  output$pointInfo <- renderText({
    ### check input
    if (v$doPlot == FALSE) return()
    
    info <- pointInfo() # info = geneID,supertaxon,maxFAS,%spec
    if(is.null(info)){return()}
    else{
      a <- toString(paste(info[1],info[2], sep = " ; "))
      b <- toString(paste(substr(info[3],6,nchar(info[3]))))
      c <- toString(paste("maxFas:",info[4],"; %spec:",info[5]))
      d <- toString(paste("traceability:",info[6]))
      paste(a,b,c,d,sep="\n")
    }
  })
  
  ### detailed species FAS scores plot
  # data for detailed plot
  detailPlotDt <- reactive({
    if (v$doPlot == FALSE) return()
    
    info <- pointInfo()  # info = geneID,supertaxon,maxFAS,%spec
    if(is.null(info)){return()}
    else{
      plotTaxon = info[3]
      plotGeneID = info[1]
      
      fullDf <- dataFiltered()
      selDf <- as.data.frame(fullDf[fullDf$geneID == plotGeneID & fullDf$supertaxon == plotTaxon,])
      selDf
    }
  })
  
  # render plot
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
  
  # plot detailed bar chart
  output$detailPlot.ui <- renderUI({
    plotOutput("detailPlot",width=400,height = input$detailedHeight,
               click = "plot_click_detail",
               hover = hoverOpts(
                 id = "plot_hover_2",
                 delay = input$hover_delay,
                 delayType = input$hover_policy,
                 nullOutside = input$hover_null_outside
               )
    )
  })
  
  ### filtered data for download
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
  
  ### download data
  output$downloadData <- downloadHandler(
    filename = function(){c("dataFiltered.out")},
    content = function(file){
      dataOut <- downloadData()
      write.table(dataOut,file,sep="\t",row.names = FALSE)
    }
  )
  
  ### data table
  output$dis <- renderDataTable({
    if(v$doPlot == FALSE){return()}
    #data <- allTaxaList()
    #data <- sortedTaxaList()
    #data <- dataFiltered()
    #data <- dataSupertaxa()
    #data <- dataHeat()
    #data <- detailPlotDt()
    data <- downloadData()
    data
  })
  
  ### show help
  output$help.ui <- renderUI({
    #    if(!file.exists("www/beschreibung.jpg")){paste("Cannot load \"beschreibung.jpg\" file in www folder!")}
    #    else{img(src="beschreibung.jpg", align = "left", height=700, width=800)}
    HTML(
      '<h1 style="color: #5e9ca0;">How the input file looks like?</h1>
      <p>Input file is a matrix of "values" (e.g. FAS scores, normalized distances, etc.), where rows represent genes and columns represent taxa.</p>
      <p>A gene may be present or absent in some taxa. A present "value" has to be in the range of 0 and 1. An absent "value" is written as NA.</p>
      <p>The header of first column has to be "geneID". The header of each taxon must have this format "ncbi12345", in which 12345 is its NCBI taxon ID.</p>
      <p><em>Still&nbsp;unclear? Take a look at the example file in /data/lca.FASmatrix :)</em></p>
      <p>&nbsp;</p>
      <h1 style="color: #5e9ca0;">Download function does not work</h1>
      <p>Problem: clicked on the "Download plot" (or Download filtered data) button, entered a file name on to "Download file" window and clicked Save, but the file...was not saved :(</p>
      <p>Solution:&nbsp;</p>
      <p>1) Use macOS to run the app :-P</p>
      <p>2) Click on "Open im Browser" to open the app using internet browser. Now the download function should work.</p>
      <p><em>I tested this function using Ubuntu 14.04 LTS and it worked with Firefox web browser.</em>&nbsp;</p>
      <p>&nbsp;&nbsp;</p>
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
      <h1 style="color: #5e9ca0;">Bug reporting</h1>
      <p>Any bug reports or comments, suggestions are highly appreciated ;-)</p>
      <p>&nbsp;</p>
      <p>&copy; 2016 Vinh Tran</p>
      <p>contact:&nbsp;<a href="mailto:tran@bio.uni-frankfurt.de">tran@bio.uni-frankfurt.de</a></p>
      <p>Please check the latest version at&nbsp;<a href="https://github.com/trvinh/phyloprofile">https://github.com/trvinh/phyloprofile</a></p>'
      )
  })
  
  ############### USED FOR TESTING
  output$testOutput <- renderText({
    # ### print infile
    #  filein <- input$file1
    #  if(is.null(filein)){return()}
    # titleline <- readLines(filein$datapath, n=1)
    # paste("perl ", getwd(),"/data/getTaxonomyInfo.pl", 
    #       #                 " -i ", getwd(),"/data/",input$file1,
    #       " -i \"", titleline,"\"", 
    #       " -n ", getwd(),"/data/taxonNamesFull.txt",
    #       " -o ", getwd(),"/data",
    #       sep='')
    
    # ### print selected supertaxon ID
    # full <- as.character(input$inSelect)
    # split <- strsplit(as.character(input$inSelect),"_")
    # inSelect <- as.integer(split[[1]][2])
    # paste(full,inSelect)
    
    ### print position of highlighted taxa
    # full <- as.character(input$inHighlight)
    # split <- strsplit(as.character(input$inHighlight),"_")
    # inHighlight <- as.integer(split[[1]][2])
    # #paste(full)
    # 
    # dataHeat <- dataHeat()
    # highlightTaxon <- toString(dataHeat[dataHeat$supertaxonID == inHighlight,2][1])
    # selectedIndex = as.numeric(as.character(substr(highlightTaxon,2,4)))
    # paste(highlightTaxon,selectedIndex)
    
    ### print selected taxonomy rank
    #    input$rankSelect
    
    ### print percentage and fas cutoff
    #    paste(input$percent,input$fas)
    #    as.numeric(as.character(substr(input$inHighlight,1,3)))
    
    ### print point info
    #    paste(pointInfo())
    
    ### print taxonName and geneID for detailed plot
    # info <- pointInfo()  # info = geneID,supertaxon,maxFAS,%spec
    # if(is.null(info)){return()}
    # else{
    #   plotTaxon = info[2]
    #   plotGeneID = info[1]
    #   paste(plotTaxon,plotGeneID)
    # }
    
    # ### print value of x and y of plot_click
    #    paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
    
    # ### print value of selected point
    #    taxaList <- as.data.frame(read.table("data/taxonNamesReduced.txt", sep='\t',header=T))
    #    inSelect <- as.numeric(taxaList$ncbiID[taxaList$fullName == input$inSelect])
    
    #    split <- strsplit(as.character(input$inSelect),"_")
    #    inSelect <- as.numeric(split[[1]][2])
    #    dataHeat <- dataHeat()
    
    
    # if (is.null(input$plot_click$x)) return()
    # else{
    # get geneID
    # genes <- as.matrix(dataHeat[dataHeat$supertaxonID == inSelect,])
    # genes[1]
    # geneID <- toString(genes[round(input$plot_click$y)])
    # geneID
    # # get supertaxon (spec)
    # supertaxa <- levels(dataHeat$supertaxon)
    # spec <- toString(supertaxa[round(input$plot_click$x)])
    # paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y, spec)
    # # get FAS and percentage of present species
    # FAS <- dataHeat$fas[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
    # Percent <- dataHeat$presSpec[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
    #
    # if(is.na(as.numeric(Percent))){return()}
    # else{
    #   info <- c(geneID,as.character(spec),round(as.numeric(FAS),2),round(as.numeric(Percent),2))
    #   #substr(spec,6,nchar(as.character(spec)))
    # }
    # }
    
    # ### list of all sequence IDs
    # data <- as.data.frame(dataHeat())
    # data$geneID <- as.character(data$geneID)
    # data$geneID <- as.factor(data$geneID)
    # out <- as.list(levels(data$geneID))
    # paste(out)
  })
  })