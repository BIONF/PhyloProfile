if (!require("shiny")) {install.packages("shiny")}
if (!require("shinyBS")) {install.packages("shinyBS")}
if (!require("ggplot2")) {install.packages("ggplot2")}
if (!require("reshape")) {install.packages("reshape")}
if (!require("plyr")) {install.packages("plyr")}
if (!require("dplyr")) {install.packages("dplyr")}
if (!require("scales")) {install.packages("scales")}
if (!require("grid")) {install.packages("grid")}
if (!require("gridExtra")) {install.packages("gridExtra")}

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
  
  ######### create taxonID.list.fullRankID and taxonNameReduced.txt from input file (if necessary)
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
#    filein <- input$file1
#    if(is.null(filein)){return()}
    selectInput("rankSelect", label = "Select taxonomy rank (for grouping):",
                choices = list("1_strain"="05_strain","2_species" = "06_species","3_genus" = "10_genus", "4_family" = "14_family", "5_order" = "19_order", "6_class" = "23_class",
                               "7_phylum" = "26_phylum", "8_kingdom" = "28_kingdom", "9_superkingdom" = "29_superkingdom","unselected"=""), 
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
    nameList <- as.data.frame(read.table("data/taxonNameReduced.txt", sep='\t',header=T,fill = TRUE))
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
    choice$fullName <- paste0(choice$fullName,"_",choice$ncbiID)
    choice$fullName <- as.factor(choice$fullName)
    selectInput('inSelect','Choose (super)taxon of interest:',as.list(levels(choice$fullName)),levels(choice$fullName)[1])
  })
  
  output$highlight = renderUI({
    choice <- allTaxaList()
    choice$fullName <- paste0(choice$fullName,"_",choice$ncbiID)
    choice$fullName <- as.factor(choice$fullName)
    selectInput('inHighlight','Select (super)taxon to highlight:',as.list(levels(choice$fullName)),levels(choice$fullName)[1])
  })
  
  ######## sorting supertaxa list
  sortedTaxaList <- reactive({
    if(v$doPlot == FALSE){return()}

    ### load list of unsorted taxa
    Dt <- as.data.frame(read.table("data/taxonID.list.fullRankID", sep='\t',header=T))

    ### load list of taxon name
    nameList <- as.data.frame(read.table("data/taxonNameReduced.txt", sep='\t',header=T,fill = TRUE))
    nameList$fullName <- as.character(nameList$fullName)

    ### input parameters
    rankSelect = input$rankSelect
    rankName = substr(rankSelect,4,nchar(rankSelect))   # get rank name from rankSelect
    rankNr = 0 + as.numeric(substr(rankSelect,1,2))     # get rank number (number of column in unsorted taxa list - dataframe Dt)

    # get selected supertaxon ID
    split <- strsplit(as.character(input$inSelect),"_")
    superID <- as.integer(split[[1]][2])

    ### sort taxa list
    ### first move all species that have the same ID of selected rank (level) to a new data frame
    ### then move species that have different ID of selected rank (level), but have the same ID of the higher level
    ### repeat until reach to superkingdom (rankNr = ncol(sortedDt)) or to the end of taxa list
    sortedDt <- data.frame()
    repeat{
      subDt <- Dt[Dt[,rankNr]==superID,]
      if(nrow(subDt) < 1){
        rankNr = rankNr + 1
        if(rankNr > ncol(sortedDt)){break}
        else {superID = sortedDt[rankNr][1,]}
      } else{
        Dt <- anti_join(Dt, subDt, by=rankName)   # delete already removed lines from Dt dataframe
        sortedDt <- rbind(sortedDt,subDt)
        rankNr = rankNr + 1
        if(rankNr > ncol(sortedDt)){break}
        else{superID = sortedDt[rankNr][1,]}
      }
      
      if(nrow(Dt) < 1 | rankNr == ncol(sortedDt)+1){
        break
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

  ######## parsing data from input matrix
  dataFiltered <- reactive({
    ##### matrix input
    filein <- input$file1
    if(is.null(filein)){return()}
    data <- as.data.frame(read.table(file=filein$datapath, sep='\t',header=T))

    # convert into paired columns
    mdData <- melt(data,id="geneID")
    colnames(mdData) <- c("geneID","ncbiID","fas")

    ##### taxonomy file input
#    taxaList <- as.data.frame(read.table("data/taxonomyList.txt", sep='\t',header=T))
    taxaList <- sortedTaxaList()

    # get frequency of all supertaxa
    taxaCount <- plyr::count(taxaList,'supertaxon')

    ### merge mdData and taxaList to get taxonomy info
    taxaMdData <- merge(mdData,taxaList,by='ncbiID')

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
    maxFasDt <- data.frame("geneID"=character(), "supertaxon"=character(),"fas"=numeric(),stringsAsFactors=FALSE)
    allGeneID <- levels(taxaMdData$geneID)

    # Create 0-row data frame which will be used to store data
    dat <- data.frame(x = numeric(0), y = numeric(0))   ### use for progess bar
    withProgress(message = 'please wait....', value = 0, {

        for(i in 1:nlevels(taxaMdData$geneID)){
          # get subset for each gene
          subDt <- taxaMdDataNoNA[taxaMdDataNoNA$geneID == allGeneID[i],]
          # get max FAS for each supertaxon of this gene
          maxSupertaxon <- by(subDt, subDt$supertaxon, function(X) X[which.max(X$fas),])
          maxSupertaxon <- do.call("rbind", maxSupertaxon)
          maxSupertaxon <- maxSupertaxon[,c("geneID","supertaxon","fas")]
          # join into maxFasDt
          maxFasDt <- rbind(maxFasDt,maxSupertaxon)

          # a stand-in for a long-running computation.
          dat <- rbind(dat, data.frame(x = rnorm(1), y = rnorm(1)))
          # Increment the progress bar, and update the detail text.
          incProgress(1/nlevels(taxaMdData$geneID), detail = paste("", percent(i/nlevels(taxaMdData$geneID))))
          # Pause for 0.1 seconds to simulate a long computation.
          Sys.sleep(0.1)
        }
    })

    ############## add presSpec and maxFAS into taxaMdData ##############
    presMdData <- merge(taxaMdData,finalPresSpecDt,by=c('geneID','supertaxon'),all.x = TRUE)
    fullMdData <- merge(presMdData,maxFasDt,by=c('geneID','supertaxon'), all.x = TRUE)

    names(fullMdData)[names(fullMdData)=="fas.x"] <- "fas"
    names(fullMdData)[names(fullMdData)=="fas.y"] <- "maxFas"

    fullMdData ### parsed input data frame !!!
  })

  ### reduce data from species level to supertaxa level
  ### this data set contain only supertaxa and their value (%present and max fas) for each gene
  dataSupertaxa <- reactive({
    fullMdData <- dataFiltered()

    ### get data set for phyloprofile plotting (contains only supertaxa info)
    superDf <- subset(fullMdData,select=c('geneID','supertaxon','supertaxonID','maxFas','presSpec','category'))
    superDf <- superDf[!duplicated(superDf), ]

    ### output
    names(superDf)[names(superDf)=="maxFas"] <- "fas"
    superDf
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
      
      selectInput('inSeq','Select sequence IDs of interest:',out,selected=out[1],multiple=FALSE)
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
  
  # ### dynamic UI (for stIndex and number of rows for profile matrix)
  # observe({
  #   filein <- input$file1
  #   if(is.null(filein)){
  #     return()
  #   }
  #   
  #   # change max number (in input$stIndex & input$number) of plot lines according to input file
  #   # & change step in stIndex according to number of rows plotting
  #   data <- dataSupertaxa()
  #   c <- length(unique(data$geneID))
  # 
  #   updateNumericInput(session, "stIndex",max = c,step = input$number)
  #   updateNumericInput(session, "number",max = c)
  # })
  
#   ### DISTRIBUTION PLOT
#   output$plot1 <- renderPlot(width=500,heigh=500,{
#     # check input file
#     filein <- input$file1
#     if(is.null(filein)){return()}
# 
#     ### get data
#     data <- dataSupertaxa()
# 
#     ### modify data
#     data[is.na(data)] <- 0
#     dataPie <- data[data$presSpec>0.0,]   # dataPie: geneID        supertaxon       fas   presSpec      category
#     dataPie$presSpec <- as.numeric(as.character(dataPie$presSpec))
#     dataPie$category <- as.factor(as.character(dataPie$category))
# 
#     levels(dataPie$category) <- c(levels(dataPie$category), toString(input$inSelect))
#     dataPie$category[dataPie$supertaxon == input$inSelect] <- toString(input$inSelect)
# 
#     ### calculate age for each group
#     dataAge <- data.frame("geneID"=character(0),"Age"=integer(0),"group"=character(0),stringsAsFactors=FALSE)
# 
#     groupList <- as.matrix(dataPie[dataPie[,2] == input$inSelect & !is.na(dataPie[,4]),]) # dataPie[,2] is supertaxon, dataPie[,4] is presSpec
#     for(i in 1:nrow(groupList)){
# #      print (groupList[i,2])
#       # get all rows of this group/gene (based on comparing geneID)
#       sub <- dataPie[dataPie[,1] == groupList[i,1],]  # groupList[i,1] is geneID of row i
#       # count categories exit in this group
#       categories <- plyr::count(sub,"category")
#       # calculate ageType
#       fungi = 0
#       unikonta = 0
#       eukaryota = 0
#       archaea = 0
#       bacteria = 0
#       for(k in 1:nrow(categories)){
#         if(!is.na(categories[k,1])){
# #         print (categories[k,])
#           if(categories[k,1] == "fungi"){fungi = 1}
#           if(categories[k,1] == "unikonta"){unikonta = 1}
#           if(categories[k,1] == "eukaryota"){eukaryota = 1}
#           if(categories[k,1] == "archaea"){archaea = 1}
#           if(categories[k,1] == "bacteria"){bacteria = 1}
#         }
#       }
#       age = as.integer(100000 + fungi*10000 + unikonta*1000 + eukaryota*100 + archaea*10 + bacteria*1)
# #      print (c(groupList[i,2],age))
#       # save to ageData
#       dataAge[i,] <- c(groupList[i,2],age,"")
#     }
# 
#     ### classify proteins into 7 categories
#     dataAge$group<-dataAge$Age
#     dataAge$group[dataAge$Age == "100000"] <- toString(input$inSelect)
#     dataAge$group[dataAge$Age == 110000] <- "2_Fungi"
#     dataAge$group[dataAge$Age == 111000] <- "3_Unikonta"
#     dataAge$group[dataAge$Age == 101000] <- "3_Unikonta"
# 
#     dataAge$group[dataAge$Age == 111100] <- "4_Eukaryota"
#     dataAge$group[dataAge$Age == 101100] <- "4_Eukaryota"
#     dataAge$group[dataAge$Age == 110100] <- "4_Eukaryota"
#     dataAge$group[dataAge$Age == 100100] <- "7_Undef"
# 
#     dataAge$group[dataAge$Age == 111110] <- "5_Archaea"
#     dataAge$group[dataAge$Age == 101110] <- "5_Archaea"
#     dataAge$group[dataAge$Age == 110110] <- "5_Archaea"
#     dataAge$group[dataAge$Age == 111010] <- "5_Archaea"
#     dataAge$group[dataAge$Age == 100110] <- "7_Undef"
#     dataAge$group[dataAge$Age == 101010] <- "7_Undef"
#     dataAge$group[dataAge$Age == 110010] <- "7_Undef"
#     dataAge$group[dataAge$Age == 100010] <- "7_Undef"
# 
#     dataAge$group[dataAge$Age == 111111] <- "6_Bacteria"
#     dataAge$group[dataAge$Age == 101111] <- "6_Bacteria"
#     dataAge$group[dataAge$Age == 110111] <- "6_Bacteria"
#     dataAge$group[dataAge$Age == 111011] <- "6_Bacteria"
#     dataAge$group[dataAge$Age == 111101] <- "6_Bacteria"
#     dataAge$group[dataAge$Age == 100111] <- "6_Bacteria"
#     dataAge$group[dataAge$Age == 101011] <- "6_Bacteria"
#     dataAge$group[dataAge$Age == 101101] <- "6_Bacteria"
#     dataAge$group[dataAge$Age == 110011] <- "7_Undef"
#     dataAge$group[dataAge$Age == 110101] <- "6_Bacteria"
#     dataAge$group[dataAge$Age == 111001] <- "6_Bacteria"
#     dataAge$group[dataAge$Age == 100011] <- "7_Undef"
#     dataAge$group[dataAge$Age == 110001] <- "7_Undef"
#     dataAge$group[dataAge$Age == 101001] <- "7_Undef"
#     dataAge$group[dataAge$Age == 100101] <- "7_Undef"
#     dataAge$group[dataAge$Age == 100001] <- "7_Undef"
# 
#     #head(dataAge)
#     dataAge$group
#     ### plot
#     factor(dataAge$group)
#     w <- plyr::count(dataAge,'group')
#     w
# 
#     # p1 = ggplot(w,aes(x="",y=freq,fill=group)) +
#     #   geom_bar(stat="identity",width=1) +
#     #   coord_polar("y",start=0)
#     # p1 = p1 + scale_fill_brewer(palette = "Set2") +
#     #   theme(axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.background = element_blank()) +
#     #   labs(x="",y="")
#     # #    p1
# 
#     p2 = ggplot(w,aes(x=group,y=freq,fill=group)) +
#       geom_bar(stat="identity") +
#       geom_text(aes(y=freq+5,label=percent(freq/sum(freq))),size=3)
#     p2 = p2 + scale_fill_brewer(palette = "Set2") +
#       theme(axis.text.x=element_text(hjust=1,angle=60),axis.text.y=element_blank(),axis.ticks=element_blank(),panel.background = element_blank()
#             ,legend.position="none") +
#       labs(x="",y="")
#       p2
# 
# #    grid.arrange(p1,p2,ncol=1)
#   })
  
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

      dataHeat <- dataHeat()
      ### plotting
      p = ggplot(dataHeat, aes(y = geneID, x = supertaxon)) +        ## global aes
        scale_fill_gradient(low = "gray95", high = "khaki", na.value="gray95", guide=FALSE) +
        geom_point(aes(colour = fas, size = presSpec))  +    ## geom_point for circle illusion
        scale_color_gradient(low = "darkorange",high = "steelblue")#+       ## color of the corresponding aes
      scale_size(range = c(0,3))             ## to tune the size of circles

      base_size <- 9
      p = p+geom_vline(xintercept=0.5,colour="dodgerblue4")
      p = p+geom_vline(xintercept=1.5,colour="dodgerblue4")
      p = p+theme(axis.text.x = element_text(angle=60,hjust=1))

      ###### highline the selected species
      ## get selected highlight taxon
      full <- as.character(input$inHighlight)
      split <- strsplit(as.character(input$inHighlight),"_")
      inHighlight <- as.integer(split[[1]][2])
      ## get taxonID together with it sorted index
      highlightTaxon <- toString(dataHeat[dataHeat$supertaxonID == inHighlight,2][1])
      ## get index
      selectedIndex = as.numeric(as.character(substr(highlightTaxon,2,4)))
      ## draw a rect to highlight this taxon's column
      rect <- data.frame(xmin=selectedIndex-0.5, xmax=selectedIndex+0.5, ymin=-Inf, ymax=Inf)
      p = p + geom_rect(data=rect, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
                    color="yellow",
                    alpha=0.3,
                    inherit.aes = FALSE)

      p
    })

  ### show help if no plot present
  output$plot.ui <- renderUI({
    if(v$doPlot == FALSE){
      if(!file.exists("www/beschreibung.png")){return(paste("WARNING: Cannot load \"beschreibung.png\" file in www folder!"))}
      else{return (img(src="beschreibung.png", align = "left", height=600, width=800))}}
    
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
      p = ggplot(dataHeat, aes(y = geneID, x = supertaxon)) +        ## global aes
        scale_fill_gradient(low = "gray95", high = "khaki", na.value="gray95", guide=FALSE) +
        geom_point(aes(colour = fas, size = presSpec))  +    ## geom_point for circle illusion
        scale_color_gradient(low = "darkorange",high = "steelblue")#+       ## color of the corresponding aes
      scale_size(range = c(0,3))             ## to tune the size of circles

      base_size <- 9
      p = p+geom_vline(xintercept=0.5,colour="dodgerblue4")
      p = p+geom_vline(xintercept=1.5,colour="dodgerblue4")
      p = p+theme(axis.text.x = element_text(angle=60,hjust=1))
      print(p)

      dev.off()
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
        scale_fill_gradient(low = "gray95", high = "khaki", na.value="gray95", guide=FALSE) +
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
    split <- strsplit(as.character(input$inSelect),"_")
    inSelect <- as.numeric(split[[1]][2])

    dataHeat <- dataHeat()

    if (is.null(input$plot_click$x)) return()
    else{
      # get geneID
      genes <- as.matrix(dataHeat[dataHeat$supertaxonID == inSelect & !is.na(dataHeat$presSpec),])
      geneID <- toString(genes[round(input$plot_click$y)])
      # get supertaxon (spec)
      supertaxa <- levels(dataHeat$supertaxon)
      spec <- toString(supertaxa[round(input$plot_click$x)])
      # get FAS and percentage of present species
      FAS <- dataHeat$fas[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]
      Percent <- dataHeat$presSpec[dataHeat$geneID == geneID & dataHeat$supertaxon == spec]

      if(is.na(as.numeric(Percent))){return()}
      else{
        info <- c(geneID,as.character(spec),round(as.numeric(FAS),2),round(as.numeric(Percent),2))
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
      a <- toString(paste(info[1],substr(info[2],6,nchar(info[2])), sep = " ; "))
      b <- toString(paste("maxFas:",info[3],"; %spec:",info[4]))
      paste(a,b,sep="\n")
#      paste("hello", "world", sep="\n")
    }
  })

  ### detailed species FAS scores plot
  output$detailPlot <- renderPlot({
    if (v$doPlot == FALSE) return()
    
    info <- pointInfo()  # info = geneID,supertaxon,maxFAS,%spec
    if(is.null(info)){return()}
    else{
      plotTaxon = info[2]
      plotGeneID = info[1]

      fullDf <- dataFiltered()
      selDf <- as.data.frame(fullDf[fullDf$geneID == plotGeneID & fullDf$supertaxon == plotTaxon,])
      selDf
      gp = ggplot(selDf, aes(y=fas,x=fullName)) +
        geom_bar(colour="steelblue", fill="steelblue", stat="identity") +
        coord_flip() +
        labs(x="") #+
      #geom_text(aes(label=fas), vjust=3)
      gp = gp+theme(axis.text.x = element_text(angle=90,hjust=1))
      gp
    }
  })

  # plot detailed bar chart
  output$detailPlot.ui <- renderUI({
    plotOutput("detailPlot",width=400,height = input$detailedHeight)
  })

  ### filtered data for download
  downloadData <- reactive({
    ### check input
    if (v$doPlot == FALSE) return()
    
    dataOut <- dataFiltered()
    dataOut <- as.data.frame(dataOut[dataOut$presSpec > 0,])
    dataOut <- dataOut[!is.na(dataOut$geneID),]
    
    dataOut <- as.data.frame(dataOut[dataOut$presSpec >= input$percent & dataOut$fas >= input$fas,])

    dataOut <- dataOut[,c("geneID","abbrName","ncbiID","supertaxon","fas")]
    dataOut <- dataOut[order(dataOut$geneID,dataOut$abbrName),]
    dataOut <- dataOut[complete.cases(dataOut),]
    
    dataOut$geneID <- as.character(dataOut$geneID)
    dataOut$abbrName <- as.character(dataOut$abbrName)
    dataOut$ncbiID <- substr(dataOut$ncbiID,5,nchar(as.character(dataOut$ncbiID)))
    dataOut$supertaxon <- substr(dataOut$supertaxon,6,nchar(as.character(dataOut$supertaxon)))
    dataOut$fas <- as.character(dataOut$fas)
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
    #data <- sortedTaxaList()
    #data <- dataFiltered()
    #data <- dataSupertaxa()
    data <- dataHeat()
    #data <- downloadData()
    data
  })
  
  ### show help
  output$help.ui <- renderUI({
    if(!file.exists("www/beschreibung.png")){paste("Cannot load \"beschreibung.png\" file in www folder!")}
    else{img(src="beschreibung.png", align = "left", height=700, width=800)}
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
    
    # ### print position of highlighted taxa
    # full <- as.character(input$inHighlight)
    # split <- strsplit(as.character(input$inHighlight),"_")
    # inHighlight <- as.integer(split[[1]][2])
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
    # paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
    # ### print value of selected point
    # split <- strsplit(as.character(input$inSelect),"_")
    # inSelect <- as.numeric(split[[1]][2])
    # dataHeat <- dataHeat()
    

    if (is.null(input$plot_click$x)) return()
    else{
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
    }
    
    # ### list of all sequence IDs
    # data <- as.data.frame(dataHeat())
    # data$geneID <- as.character(data$geneID)
    # data$geneID <- as.factor(data$geneID)
    # out <- as.list(levels(data$geneID))
    # paste(out)
 })
})

