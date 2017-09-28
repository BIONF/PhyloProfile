#############################################################
######################## FUNCTIONS ##########################
#############################################################

########## parse orthoXML file ##############
xmlParser <- function(inputFile){
  cmd <- paste("python ", getwd(),"/bin/orthoxmlParser.py",
               " -i ", inputFile,
               sep='')
  dfIN <- as.data.frame(read.table(text = system(cmd,intern=TRUE)))
  
  colnames(dfIN) = as.character(unlist(dfIN[1,])) # the first row will be the header
  dfIN <- subset(dfIN[dfIN$geneID != "geneID",])
  dfIN <- droplevels(dfIN)
  dfIN
}

# ########## convert long to wide format ##############
# long2wide <- function(longDf){
#   # rename column names
#   colnames(longDf) <- c("geneID","ncbiID","orthoID","var1","var2")
#   longDf$value <- paste0(longDf$orthoID,"#",longDf$var1,"#",longDf$var2)
#   longDfmod <- longDf[,c("geneID","ncbiID","value")]
#
#   # count paralogs
#   longDfmod <- data.table(longDfmod)
#   longDfmod[ ,paralog := 1:.N, by=c("geneID","ncbiID")]
#   longDfmod <- data.frame(longDfmod)
#
#   # return wide data frame
#   wideDf <- spread(longDfmod, ncbiID, value)
#   wideDf <- subset(wideDf,paralog == 1)   # remove co-orthologs
#   wideDf <- subset(wideDf, select=-paralog)
# }

############### FUNCTION FOR CLUSTERING PROFILES  ###############
clusteredGeneList <- function(data,distMethod,clusterMethod){
  
  # do clustering
  row.order <- hclust(dist(data, method = distMethod), method = clusterMethod)$order
  col.order <- hclust(dist(t(data), method = distMethod), method = clusterMethod)$order
  
  # re-order data accoring to clustering
  dat_new <- data[row.order, col.order]
  
  # return clustered gene ID list
  clusteredGeneIDs <- as.factor(row.names(dat_new))
  clusteredGeneIDs
}

########## calculate percentage of present species ##########
calcPresSpec <- function(taxaMdData, taxaCount){
  ### taxaMdData = df("geneID","ncbiID","orthoID","var1","var2","paralog",....,"supertaxon")
  taxaMdData <- taxaMdData[taxaMdData$orthoID != "NA",]
  
  # get geneID and supertaxon
  geneIDsupertaxon <- subset(taxaMdData,select=c('geneID','supertaxon','paralog'))
  geneIDsupertaxon <- geneIDsupertaxon[!duplicated(geneIDsupertaxon),] # remove duplicated rows
  
  # remove NA rows from taxaMdData
  taxaMdDataNoNA <- taxaMdData[taxaMdData$orthoID != "NA",]
  
  # count present frequency of supertaxon for each gene
  geneSupertaxonCount <- plyr::count(taxaMdDataNoNA,c('geneID','supertaxon'))
  
  # merge with taxaCount to get total number of species of each supertaxon and calculate presSpec
  presSpecDt <- merge(geneSupertaxonCount,taxaCount,by='supertaxon')
  presSpecDt <- merge(presSpecDt,geneIDsupertaxon,by=c('geneID','supertaxon'))
  presSpecDt$presSpec <- ((presSpecDt$freq.x-presSpecDt$paralog)+1)/presSpecDt$freq.y
  
  presSpecDt <- presSpecDt[presSpecDt$presSpec <= 1,]
  presSpecDt <- presSpecDt[order(presSpecDt$geneID),]
  presSpecDt <- presSpecDt[,c("geneID","supertaxon","presSpec")]
  
  # add absent supertaxon into presSpecDt
  geneIDsupertaxon <- subset(geneIDsupertaxon, select = -c(paralog))
  finalPresSpecDt <- merge(presSpecDt,geneIDsupertaxon,by=c('geneID','supertaxon'),all.y = TRUE)
  finalPresSpecDt$presSpec[is.na(finalPresSpecDt$presSpec)] <- 0
  
  # return finalPresSpecDt
  finalPresSpecDt <- finalPresSpecDt[!duplicated(finalPresSpecDt),] # remove duplicated rows
  return(finalPresSpecDt)
}

######## sort one domain dataframe (ortho) based on the other domain Df (seed) ########
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

######## plot domain architecture ########
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

######## plot profile heatmap ########
heatmap.plotting <- function(data,xAxis,var1_id,var2_id,lowColor_var1,highColor_var1,lowColor_var2,highColor_var2,paraColor,xSize,ySize,legendSize,mainLegend,dotZoom,guideline){
  dataHeat <- data
  
  ### rescale numbers of paralogs
  dataHeat$paralog <- as.numeric(dataHeat$paralog)
  if(length(unique(na.omit(dataHeat$paralog))) > 0){
    maxParalog <- max(na.omit(dataHeat$paralog))
    dataHeat$paralogSize <- (dataHeat$paralog/maxParalog)*3
  }
  
  ### remove prefix number of taxa names but keep the order
  dataHeat$supertaxon <- mapvalues(warn_missing=F,dataHeat$supertaxon,from=as.character(dataHeat$supertaxon),to=substr(as.character(dataHeat$supertaxon),6,nchar(as.character(dataHeat$supertaxon))))
  
  ### format plot
  if(xAxis == "genes"){
    p = ggplot(dataHeat, aes(x = geneID, y = supertaxon))        ## global aes
  } else{
    p = ggplot(dataHeat, aes(y = geneID, x = supertaxon))        ## global aes
  }
  
  if(length(unique(na.omit(dataHeat$var2))) != 1){
    p = p + scale_fill_gradient(low = lowColor_var2, high = highColor_var2, na.value="gray95", limits=c(0,1)) +   ## fill color (var2)
      geom_tile(aes(fill = var2))    ## filled rect (var2 score)
  }
  
  if(length(unique(na.omit(dataHeat$presSpec))) < 3){
    if(length(unique(na.omit(dataHeat$var1))) == 1){
      p = p + geom_point(aes(colour = var1),size = dataHeat$presSpec*5*(1+dotZoom),na.rm = TRUE,show.legend=F)    ## geom_point for circle illusion (var1 and presence/absence)
    } else {
      p = p + geom_point(aes(colour = var1),size = dataHeat$presSpec*5*(1+dotZoom),na.rm = TRUE)    ## geom_point for circle illusion (var1 and presence/absence)
      p = p + scale_color_gradient(low = lowColor_var1,high = highColor_var1, limits=c(0,1)) ## color of the corresponding aes (var1)
    }
  } else {
    if(length(unique(na.omit(dataHeat$var1))) == 1){
      p = p + geom_point(aes(size = presSpec),color = "#336a98",na.rm = TRUE)    ## geom_point for circle illusion (var1 and presence/absence)
    } else {
      p = p + geom_point(aes(colour = var1, size = presSpec),na.rm = TRUE)    ## geom_point for circle illusion (var1 and presence/absence)
      p = p + scale_color_gradient(low = lowColor_var1,high = highColor_var1, limits=c(0,1)) ## color of the corresponding aes (var1)
    }
  }
  
  # plot inparalogs (if available)
  if(length(unique(na.omit(dataHeat$paralog))) > 0){
    p <- p + geom_point(data = dataHeat, aes(size = paralog), color = paraColor, na.rm = TRUE, show.legend = TRUE)
    p <- p + guides(size=guide_legend(title = "# inparalogs"))
    p <- p + scale_size_continuous(range = c(min(na.omit(dataHeat$paralogSize))*(1+dotZoom),max(na.omit(dataHeat$paralogSize))*(1+dotZoom)))  ## to tune the size of circles; "floor(value*10)/10" is used to round "down" the value with one decimal number
  } else {
    # remain the scale of point while filtering
    presentVl <- dataHeat$presSpec[!is.na(dataHeat$presSpec)]
    p <- p + scale_size_continuous(range = c((floor(min(presentVl)*10)/10*5)*(1+dotZoom),(floor(max(presentVl)*10)/10*5)*(1+dotZoom)))  ## to tune the size of circles; "floor(value*10)/10" is used to round "down" the value with one decimal number
  }
  
  #
  p = p + guides(fill=guide_colourbar(title = var2_id), color=guide_colourbar(title = var1_id))
  base_size <- 9
  
  ### guideline for separating ref species
  if(guideline == 1){
    if(xAxis == "genes"){
      p = p + labs(y="Taxon")
      p = p+geom_hline(yintercept=0.5,colour="dodgerblue4")
      p = p+geom_hline(yintercept=1.5,colour="dodgerblue4")
    } else{
      p = p + labs(x="Taxon")
      p = p+geom_vline(xintercept=0.5,colour="dodgerblue4")
      p = p+geom_vline(xintercept=1.5,colour="dodgerblue4")
    }
  }
  
  ### format theme
  p = p + theme_minimal()
  p = p + theme(axis.text.x = element_text(angle=60,hjust=1,size=xSize),axis.text.y = element_text(size=ySize),
                axis.title.x = element_text(size=xSize), axis.title.y = element_text(size=ySize),
                legend.title=element_text(size=legendSize),legend.text=element_text(size=legendSize),legend.position = mainLegend)
  
  ### return plot
  return(p)
}

######## show FASTA sequence in popup windows of selected plot  ########
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

######## get last n characters from string x  ########
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}