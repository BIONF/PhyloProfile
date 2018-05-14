# =============================================================================
# FUNCTIONS ===================================================================
# =============================================================================

# PARSE INPUT FILTE ===========================================================

# parse orthoXML file ---------------------------------------------------------
xmlParser <- function(inputFile){
  cmd <- paste("python ",
               getwd(),
               "/scripts/orthoxmlParser.py",
               " -i ",
               inputFile,
               sep = "")
  
  dfIN <- as.data.frame(read.table(text = system(cmd, intern = TRUE)))
  
  # the first row will be the header
  colnames(dfIN) = as.character(unlist(dfIN[1,])) 
  
  dfIN <- subset(dfIN[dfIN$geneID != "geneID",])
  dfIN <- droplevels(dfIN)
  dfIN
}

# parse fasta input file ------------------------------------------------------
fastaParser <- function(inputFile){
  
  cmd <- paste("python ",
                   getwd(),
                   "/scripts/fastaParser.py",
                   " -i ",
                   inputFile,
                   sep = "")
  
  dfIN <- as.data.frame(read.table(text = system(cmd, intern = TRUE)))
  
  # the first row will be the header
  colnames(dfIN) = as.character(unlist(dfIN[1,])) 
  dfIN <- subset(dfIN[dfIN$geneID != "geneID",])
  dfIN <- droplevels(dfIN)
  
  if(all(is.na(dfIN$var2))){
    dfIN = subset(dfIN, select = -c(var2) )
  }
  if(all(is.na(dfIN$var1))){
    dfIN = subset(dfIN, select = -c(var1) )
  }
  
  dfIN
}

# # parse OMA input file ------------------------------------------------------
# oma_parser <- function(input_file, selected_type ){
#   print("oma parser")
#   
#   # use the python script to get the information from oma
#   cmd <- paste("python ", 
#                getwd(),
#                "/scripts/get_oma_browser.py",
#                " -i ",
#                input_file,
#                " -t ",
#                selected_type,
#                sep = "")
#   print(cmd)
#   dfIN <- as.data.frame(read.table(text = system(cmd, intern = TRUE)))
#   print(dfIN)
#   # the first row will be the header
#   colnames(dfIN) = as.character(unlist(dfIN[1, ]))
#   dfIN <- subset(dfIN[dfIN$geneID != "geneID", ])
#   dfIN <- droplevels(dfIN)
#   dfIN
# }


# PLOTTING ====================================================================
# plot domain architecture ----------------------------------------------------
domain.plotting <- function(df,
                            geneID,
                            sep,
                            labelSize,
                            titleSize,
                            minStart,
                            maxEnd){
  gg <- ggplot(df, aes(y = feature, x = end, color = feature)) +
    geom_segment(data = df,
                 aes(y = feature, yend = feature,
                     x = minStart, xend = maxEnd),
                 color = "white",
                 size = 0)
  
  # draw lines for representing sequence length
  gg <- gg + geom_segment(data=df,
                          aes(x = 0, xend = length,
                              y = feature, yend = feature),
                          size = 1,
                          color = "#b2b2b2")
  
  # draw line and points
  gg <- gg + geom_segment(data=df,
                          aes(x = start, xend = end,
                              y = feature, yend = feature),
                          size = 1.5)
  gg <- gg + geom_point(data = df,
                        aes(y = feature, x=start),
                        color = "#b2b2b2",
                        size = 3,
                        shape = 3)
  gg <- gg + geom_point(data = df,
                        aes(y = feature, x = end),
                        color = "#edae52",
                        size = 3,
                        shape = 5)
  
  # draw dashed line for domain path
  gg <- gg + geom_segment(data = df[df$path == "Y", ], 
                          aes(x = start, xend = end,
                              y = feature, yend = feature),
                          size = 3,
                          linetype = "dashed")
  
  # # add text above
  # gg <- gg + geom_text(data=df,
  #                      aes(x = (start + end) / 2,
  #                          y = feature, label = round(weight, 2)),
  #                      color = "#9fb059",
  #                      size = descSize,
  #                      vjust = -0.75,
  #                      fontface = "bold", family = "serif")
  
  # theme format
  titleMod <- gsub(":",sep,geneID)
  gg <- gg + scale_y_discrete(expand=c(0.075,0))
  gg <- gg + labs(title=paste0(titleMod), y="Feature")
  gg <- gg + theme_minimal()
  gg <- gg + theme(panel.border=element_blank())
  gg <- gg + theme(axis.ticks=element_blank())
  gg <- gg + theme(plot.title=element_text(face="bold",size=titleSize))
  gg <- gg + theme(plot.title=element_text(hjust = 0.5))
  gg <- gg + theme(legend.position="none",axis.title.x=element_blank(),
                   axis.text.y = element_text(size=labelSize),
                   axis.title.y = element_text(size=labelSize),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.x = element_blank())
  
  # return plot
  return(gg)
}

# plot profile heatmap --------------------------------------------------------
heatmap.plotting <- function(data,
                             xAxis,
                             var1_id,
                             var2_id,
                             lowColor_var1,
                             highColor_var1,
                             lowColor_var2,
                             highColor_var2,
                             paraColor,
                             xSize,
                             ySize,
                             legendSize,
                             mainLegend,
                             dotZoom,
                             xAngle,
                             guideline){
  dataHeat <- data
  
  # rescale numbers of paralogs
  dataHeat$paralog <- as.numeric(dataHeat$paralog)
  if(length(unique(na.omit(dataHeat$paralog))) > 0){
    maxParalog <- max(na.omit(dataHeat$paralog))
    dataHeat$paralogSize <- (dataHeat$paralog / maxParalog) * 3
  }
  
  # remove prefix number of taxa names but keep the order
  dataHeat$supertaxon <- {
    mapvalues(warn_missing = F,
              dataHeat$supertaxon,
              from = as.character(dataHeat$supertaxon),  
              to = substr(as.character(dataHeat$supertaxon),
              6,
              nchar(as.character(dataHeat$supertaxon))))
  }
  
  # format plot
  if(xAxis == "genes"){
    p = ggplot(dataHeat, aes(x = geneID, y = supertaxon)) # global aes
  } else{
    p = ggplot(dataHeat, aes(y = geneID, x = supertaxon)) # global aes
  }
  
  if(length(unique(na.omit(dataHeat$var2))) != 1){
    p = p + scale_fill_gradient(low = lowColor_var2,
                                high = highColor_var2,
                                na.value = "gray95",
                                limits = c(0,1)) +  # fill color (var2)
      geom_tile(aes(fill = var2))    # filled rect (var2 score)
  }
  
  if(length(unique(na.omit(dataHeat$presSpec))) < 3){
    if(length(unique(na.omit(dataHeat$var1))) == 1){
      # geom_point for circle illusion (var1 and presence/absence)
      p = p + geom_point(aes(colour = var1),
                         size = dataHeat$presSpec * 5 * (1 + dotZoom),
                         na.rm = TRUE,show.legend=F)  
    } else {
      # geom_point for circle illusion (var1 and presence/absence)
      p = p + geom_point(aes(colour = var1),
                         size = dataHeat$presSpec * 5 * (1 + dotZoom),
                         na.rm = TRUE) 
      # color of the corresponding aes (var1)
      p = p + scale_color_gradient(low = lowColor_var1,
                                   high = highColor_var1,
                                   limits=c(0,1)) 
    }
  } else {
    if(length(unique(na.omit(dataHeat$var1))) == 1){
      # geom_point for circle illusion (var1 and presence/absence)
      p = p + geom_point(aes(size = presSpec),
                         color = "#336a98",
                         na.rm = TRUE)    
    } else {
      # geom_point for circle illusion (var1 and presence/absence)
      p = p + geom_point(aes(colour = var1, size = presSpec),
                         na.rm = TRUE)
      # color of the corresponding aes (var1)
      p = p + scale_color_gradient(low = lowColor_var1,
                                   high = highColor_var1,
                                   limits=c(0,1))
    }
  }
  
  # plot inparalogs (if available)
  if(length(unique(na.omit(dataHeat$paralog))) > 0){
    p <- p + geom_point(data = dataHeat,
                        aes(size = paralog),
                        color = paraColor,
                        na.rm = TRUE,
                        show.legend = TRUE)
    p <- p + guides(size = guide_legend(title = "# of co-orthologs"))
    
    # to tune the size of circles;
    # "floor(value*10)/10" is used to round "down" the value
    # with one decimal number
    min_value <- min(na.omit(dataHeat$paralogSize)) * (1 + dotZoom)
    max_value <- max(na.omit(dataHeat$paralogSize)) * (1 + dotZoom)
  
    p <- p + scale_size_continuous(range = c( min_value, max_value))
  } else {
    # remain the scale of point while filtering
    presentVl <- dataHeat$presSpec[!is.na(dataHeat$presSpec)]
    
    # to tune the size of circles; 
    # "floor(value*10)/10" is used to round "down" the value
    # with one decimal number
    min_value <- (floor(min(presentVl) * 10) / 10 * 5) * (1 + dotZoom)
    max_value <- (floor(max(presentVl) * 10) / 10 * 5) * (1 + dotZoom)
    p <- p + scale_size_continuous(range = c( min_value, max_value))  
  }
  
  p = p + guides(fill = guide_colourbar(title = var2_id),
                 color = guide_colourbar(title = var1_id))
  base_size <- 9
  
  # guideline for separating ref species
  if(guideline == 1){
    if(xAxis == "genes"){
      p = p + labs(y="Taxon")
      p = p + geom_hline(yintercept = 0.5, colour = "dodgerblue4")
      p = p + geom_hline(yintercept = 1.5, colour = "dodgerblue4")
    } else{
      p = p + labs(x="Taxon")
      p = p + geom_vline(xintercept = 0.5,colour = "dodgerblue4")
      p = p + geom_vline(xintercept = 1.5, colour = "dodgerblue4")
    }
  }
  
  # format theme
  p = p + theme_minimal()
  p = p + theme(axis.text.x = element_text(angle = xAngle,
                                           hjust = 1,
                                           size = xSize),
                axis.text.y = element_text(size = ySize),
                axis.title.x = element_text(size = xSize),
                axis.title.y = element_text(size = ySize),
                legend.title = element_text(size = legendSize), 
                legend.text = element_text(size = legendSize),
                legend.position = mainLegend)
  
  # return plot
  return(p)
}

# heatmap data input ----------------------------------------------------------
dataHeat <- function(percent,
                     var1, var2,
                     filein,
                     demo_data,
                     var1_relation, var2_relation,
                     var1_aggregate_by, var2_aggregate_by,
                     end_index,
                     selected,
                     listIn,
                     st_index,
                     ordering,
                     gene_list_selected,
                     v,
                     rank_select,
                     in_select,
                     treeIn) {
  

  # get all cutoffs
  percent_cutoff_min <- percent[1]
  percent_cutoff_max <- percent[2]
  var1_cutoff_min <- var1[1]
  var1_cutoff_max <- var1[2]
  var2_cutoff_min <- var2[1]
  var2_cutoff_max <- var2[2]
  
  
  # if(input$demo == TRUE){
  if(demo_data == "demo" | demo_data == "ampk-tor"){
    filein = 1
  }
  if(is.null(filein)) return()
  dataHeat <- dataSupertaxa(var1_aggregate_by, var2_aggregate_by,
                            end_index,
                            selected,
                            listIn,
                            st_index,
                            demo_data,
                            ordering,
                            filein,
                            gene_list_selected,
                            v,
                            rank_select,
                            in_select,
                            treeIn)
  
  # get selected supertaxon name
  split <- strsplit(as.character(in_select), "_")
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
  
  if(var1_relation == "protein"){
    if(var2_relation == "protein"){
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
    if(var2_relation == "species"){
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
}

# create profile heatmap ------------------------------------------------------
mainPlot <- function(apply_cluster,
                     x_axis,
                     var1_id,
                     var2_id,
                     low_color_var1,
                     high_color_var1,
                     low_color_var2,
                     high_color_var2,
                     para_color,
                     x_size,
                     y_size,
                     legend_size,
                     main_legend,
                     dot_zoom,
                     x_angle,
                     taxon_highlight,
                     rank_select,
                     gene_highlight,
                     auto_update,
                     update_btn,
                     v,
                     percent,
                     var1,
                     var2,
                     mainInput,
                     demo_data,
                     var1_relation,
                     var2_relation,
                     var1_aggregate_by, var2_aggregate_by,
                     end_index,
                     selected,
                     listIn,
                     st_index,
                     ordering,
                     gene_list_selected,
                     in_select,
                     treeIn){
  if (v$doPlot == FALSE) return()
  dataHeat <- dataHeat(percent,
                       var1,
                       var2,
                       mainInput,
                       demo_data,
                       var1_relation,
                       var2_relation,
                       var1_aggregate_by, var2_aggregate_by,
                       end_index,
                       selected,
                       listIn,
                       st_index,
                       ordering,
                       gene_list_selected,
                       v,
                       rank_select,
                       in_select,
                       treeIn)
  
  ### cluster dataHeat (if selected)
  if(apply_cluster == TRUE){
    dataHeat <- clusteredDataHeat()
  }
  
  ### reduce number of inparalogs based on filtered dataHeat
  dataHeatTB <- data.table(na.omit(dataHeat))
  dataHeatTB[, paralogNew := .N, by = c("geneID","supertaxon")]
  dataHeatTB <- data.frame(dataHeatTB[, c("geneID",
                                          "supertaxon",
                                          "paralogNew")])
  
  dataHeat <- merge(dataHeat, dataHeatTB,
                    by = c("geneID", "supertaxon"),
                    all.x = TRUE)
  dataHeat$paralog <- dataHeat$paralogNew
  dataHeat <- dataHeat[!duplicated(dataHeat), ]
  
  ### remove unneeded dots
  dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
  dataHeat$paralog[dataHeat$presSpec < 1] <- NA
  dataHeat$paralog[dataHeat$paralog == 1] <- NA
  
  p <- heatmap.plotting(dataHeat,
                        x_axis,
                        var1_id,
                        var2_id,
                        low_color_var1,
                        high_color_var1,
                        low_color_var2,
                        high_color_var2,
                        para_color,
                        x_size,
                        y_size,
                        legend_size,
                        main_legend,
                        dot_zoom,
                        x_angle,
                        1)
  
  # highlight taxon
  if(taxon_highlight != "none"){
    # get selected highlight taxon ID
    rank_select = rank_select
    rankName = substr(rank_select,
                      4,
                      nchar(rank_select))   # get rank name from rank_select
    taxa_list <- as.data.frame(read.table("data/taxonNamesReduced.txt",
                                          sep = "\t",
                                          header = T))
    taxonHighlight <- {
      taxa_list$ncbiID[taxa_list$fullName == taxon_highlight
                       & taxa_list$rank == rankName]
    }
    if(length(taxonHighlight) == 0L){
      taxonHighlight <- {
        taxa_list$ncbiID[taxa_list$fullName == taxon_highlight]
      }
    }
    
    # get taxonID together with it sorted index
    highlightTaxon <- {
      toString(dataHeat[dataHeat$supertaxonID == taxonHighlight, 2][1])
    }
    
    # get index
    selectedIndex = as.numeric(as.character(substr(highlightTaxon,2,4)))
    
    # draw a rect to highlight this taxon's column
    if(x_axis == "taxa"){
      rect <- data.frame(xmin = selectedIndex - 0.5,
                         xmax = selectedIndex + 0.5,
                         ymin = -Inf,
                         ymax = Inf)
    } else {
      rect <- data.frame(ymin = selectedIndex - 0.5,
                         ymax = selectedIndex + 0.5,
                         xmin = -Inf,
                         xmax = Inf)
    }
    
    p = p + geom_rect(data = rect,
                      aes(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
                      color="yellow",
                      alpha=0.3,
                      inherit.aes = FALSE)
  }
  
  # highlight gene
  if(gene_highlight != "none"){
    # get selected highlight gene ID
    geneHighlight <- gene_highlight
    
    # get index
    allGenes <- levels(dataHeat$geneID)
    selectedIndex = match(geneHighlight, allGenes)
    
    # draw a rect to highlight this taxon's column
    if(x_axis == "taxa"){
      rect <- data.frame(ymin = selectedIndex - 0.5,
                         ymax = selectedIndex + 0.5,
                         xmin = -Inf,
                         xmax = Inf)
    } else {
      rect <- data.frame(xmin = selectedIndex - 0.5,
                         xmax = selectedIndex + 0.5,
                         ymin = -Inf,
                         ymax = Inf)
    }
    
    p = p + geom_rect(data = rect,
                      aes(xmin = xmin, xmax = xmax,
                          ymin = ymin, ymax = ymax),
                      color = "yellow",
                      alpha = 0.3,
                      inherit.aes = FALSE)
  }
  
  # do plotting
  if(auto_update == FALSE){
    # Add dependency on the update button
    # (only update when button is clicked)
    update_btn
    
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


# var1 score distribution plot ------------------------------------------------
var1DistPlot <- function(dist_text_size, var1_id, v){
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
                   axis.title = element_text(size = dist_text_size),
                   axis.text = element_text(size = dist_text_size)) +
      labs(x = paste0(var1_id,
                      " (mean = ",round(mean(splitDt$var1),3),")"),
           y = "Frequency")
    p
  }
}

# var2 score distribution plot ------------------------------------------------
var2DistPlot <- function(dist_text_size, var2_id, v){
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
                   axis.title = element_text(size = dist_text_size),
                   axis.text = element_text(size = dist_text_size)) +
      labs(x = paste0(var2_id," (mean = ",round(mean(splitDt$var2),3),")"), y = "Frequency")
    p
  }
}


# % presSpec distribution plot ------------------------------------------------
presSpecPlot <- function(percent, dist_text_size, v){
  if (v$doPlot == FALSE) return()
  
  # data
  dt <- presSpecAllDt()
  # remove presSpec < cutoff_min or > cutoff_max
  if(percent[1] > 0){
    dt <- dt[dt$presSpec >= percent[1] & dt$presSpec <= percent[2],]
  } else {
    if(percent[2] > 0){
      dt <- dt[dt$presSpec > 0 & dt$presSpec <= percent[2],]
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
                 axis.title = element_text(size = dist_text_size),
                 axis.text = element_text(size = dist_text_size)) +
    labs(x = paste0("% present taxa (mean = ",
                    round(mean(dt$presSpec),3),
                    ")"),
         y = "Frequency")
  p
}

# create plot (same as main plot) ---------------------------------------------
selected_plot <- function(in_seq, in_taxa, apply_cluster,
                          x_axis_selected,
                          var1_id,
                          var2_id,
                          low_color_var1,
                          high_color_var1,
                          low_color_var2,
                          high_color_var2,
                          para_color,
                          x_size_select,
                          y_size_select,
                          legend_size_select,
                          selected_legend,
                          dot_zoom_select,
                          x_angle_select,
                          auto_update_selected,
                          plot_custom,
                          percent,
                          var1,
                          var2,
                          mainInput,
                          demo_data,
                          var1_relation,
                          var2_relation,
                          var1_aggregate_by, var2_aggregate_by,
                          end_index,
                          selected,
                          listIn,
                          st_index,
                          ordering,
                          gene_list_selected,
                          v, rank_select, in_select, treeIn){
  if (vCt$doPlotCustom == FALSE) return()
  if(in_seq[1] == "all" & in_taxa[1] == "all") return()
  else{
    dataHeat <- dataHeat(percent,
                         var1,
                         var2,
                         mainInput,
                         demo_data,
                         var1_relation,
                         var2_relation,
                         var1_aggregate_by, var2_aggregate_by,
                         end_index,
                         selected,
                         listIn,
                         st_index,
                         ordering,
                         gene_list_selected,
                         v,
                         rank_select,
                         in_select,
                         treeIn)
    
    ### cluster dataHeat (if selected)
    if(apply_cluster == TRUE){
      dataHeat <- clusteredDataHeat()
    }
    
    ### process data
    dataHeat$supertaxonMod <- substr(dataHeat$supertaxon,
                                     6,
                                     nchar(as.character(dataHeat$supertaxon)))
    if(in_taxa[1] == "all" & in_seq[1] != "all"){
      # select data from dataHeat for selected sequences only
      dataHeat <- subset(dataHeat,geneID %in% in_seq) 
    } else if(in_seq[1] == "all" & in_taxa[1] != "all"){
      # select data from dataHeat for selected taxa only
      dataHeat <- subset(dataHeat,supertaxonMod %in% in_taxa) 
    } else {
      # select data from dataHeat for selected sequences and taxa
      dataHeat <- subset(dataHeat,geneID %in% in_seq & supertaxonMod %in% in_taxa) 
    }
    
    ### remove unneeded dots
    dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
    dataHeat$paralog[dataHeat$presSpec < 1] <- NA
    dataHeat$paralog[dataHeat$paralog == 1] <- NA
    
    ### create plot
    p <- heatmap.plotting(dataHeat,
                          x_axis_selected,
                          var1_id,
                          var2_id,
                          low_color_var1,
                          high_color_var1,
                          low_color_var2,
                          high_color_var2,
                          para_color,
                          x_size_select,
                          y_size_select,
                          legend_size_select,
                          selected_legend,
                          dot_zoom_select,
                          x_angle_select,0)
    
    ### do plotting
    if(auto_update_selected == FALSE){
      # Add dependency on the update button (only update when button is clicked)
      plot_custom
      
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

# render detailed plot --------------------------------------------------------
detail_plot <- function(var1_id, var2_id, detailed_text, v,
                        input_tabs,
                        var1_aggregate_by,
                        var2_aggregate_by,
                        end_index,
                        selected,
                        listIn,
                        st_index,
                        demo_data,
                        ordering,
                        mainInput,
                        gene_list_selected,
                        rank_select,
                        in_select,
                        inputTree,
                        percent,
                        var1, var2,
                        x_axis){
  if (v$doPlot == FALSE) return()
  
  selDf <- detail_plotDt(v, input_tabs,
                         var1_aggregate_by,
                         var2_aggregate_by,
                         end_index,
                         selected,
                         listIn,
                         st_index,
                         demo_data,
                         ordering,
                         mainInput,
                         gene_list_selected,
                         inputTree,
                         percent,
                         var1, var2,
                         x_axis)
  selDf$x_label <- paste(selDf$orthoID," (",selDf$fullName,")",sep = "")
  
  # if(input$detailed_remove_na == TRUE){
  #   selDf <- selDf[!is.na(selDf$orthoID),]
  # }
  
  ### create joined DF for plotting var1 next to var2
  var1Df <- subset(selDf, select = c("x_label","var1"))
  var1Df$type <- var1_id
  colnames(var1Df) <- c("id","score","var")
  
  var2Df <- subset(selDf, select = c("x_label","var2"))
  var2Df$type <- var2_id
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
                axis.text = element_text(size = detailed_text),
                axis.title = element_text(size = detailed_text),
                legend.text = element_text(size = detailed_text)
  )
  gp
}

# create domain plot ----------------------------------------------------------
archi_plot <- function(demo_data,
                       one_seq_fasta,
                       label_archi_size,
                       title_archi_size){
  if (v3$doPlot3 == FALSE) return()
  ### info
  info <- point_infoDetail() # info = seedID, orthoID, var1
  group <- as.character(info[1])
  ortho <- as.character(info[2])
  var1 <- as.character(info[3])
  
  ### parse domain file
  fileDomain <- getDomainFile()
  # if(input$demo == TRUE){
  if(demo_data == "demo" | demo_data == "ampk-tor"){
    domain_df <- as.data.frame(read.csv(fileDomain, sep="\t", header=F, comment.char = "", stringsAsFactors = FALSE, quote = ""))
    
    if(ncol(domain_df) == 5){
      colnames(domain_df) <- c("seedID",
                               "orthoID",
                               "feature",
                               "start",
                               "end")
      
    } else if(ncol(domain_df) == 6){
      colnames(domain_df) <- c("seedID",
                               "orthoID",
                               "feature",
                               "start",
                               "end",
                               "weight")
      
    } else if(ncol(domain_df) == 7){
      colnames(domain_df) <- c("seedID",
                               "orthoID",
                               "feature",
                               "start",
                               "end",
                               "weight",
                               "path")
    }
    domain_df$length <- max(domain_df$end)
    
  } else {
    if(fileDomain != FALSE){
      domain_df <- as.data.frame(read.table(fileDomain,
                                            sep="\t",
                                            header = FALSE,
                                            comment.char=""))
    }
    
    if(ncol(domain_df) == 6){
      colnames(domain_df) <- c("seedID",
                               "orthoID",
                               "length",
                               "feature",
                               "start",
                               "end")
      
    } else if(ncol(domain_df) == 7){
      colnames(domain_df) <- c("seedID",
                               "orthoID",
                               "length",
                               "feature",
                               "start",
                               "end",
                               "weight")
      
    } else if(ncol(domain_df) == 8){
      colnames(domain_df) <- c("seedID",
                               "orthoID",
                               "length",
                               "feature",
                               "start",
                               "end",
                               "weight",
                               "path")
    }
  }
  
  domain_df$seedID <- gsub("\\|", ":", domain_df$seedID)
  domain_df$orthoID <- gsub("\\|", ":", domain_df$orthoID)
  
  # get sub dataframe based on selected groupID and orthoID
  ortho <- gsub("\\|", ":", ortho)
  grepID = paste(group, "#", ortho, sep = "")
  subdomain_df <- domain_df[grep(grepID, domain_df$seedID), ]
  subdomain_df$feature <- as.character(subdomain_df$feature)
  
  if(nrow(subdomain_df) < 1){
    v3$doPlot3 = FALSE
    return()
  } else {
    
    # ortho domains df
    orthoDf <- filter(subdomain_df, orthoID == ortho)
    
    # seed domains df
    seedDf <- filter(subdomain_df, orthoID != ortho)
    
    if(nrow(seedDf) == 0){seedDf <- orthoDf}
    
    seed = as.character(seedDf$orthoID[1])
    
    # change order of one dataframe's features based on order of other df's features
    if(length(orthoDf$feature) < length(seedDf$feature)){
      orderedOrthoDf <- orthoDf[order(orthoDf$feature), ]
      orderedSeedDf <- sortDomains(orderedOrthoDf, seedDf)
    } else {
      orderedSeedDf <- seedDf[order(seedDf$feature), ]
      orderedOrthoDf <- sortDomains(orderedSeedDf, orthoDf)
    }
    
    # join weight values and feature names
    if("weight" %in% colnames(orderedOrthoDf)){
      orderedOrthoDf$yLabel <- paste0(orderedOrthoDf$feature," (",round(orderedOrthoDf$weight,2),")")
      orderedOrthoDf$feature <- orderedOrthoDf$yLabel
    }
    if("weight" %in% colnames(orderedSeedDf)){
      orderedSeedDf$yLabel <- paste0(orderedSeedDf$feature," (",round(orderedSeedDf$weight,2),")")
      orderedSeedDf$feature <- orderedSeedDf$yLabel
    }
    
    # plotting
    sep = ":"
    if(!is.null(one_seq_fasta)) sep <- "|"
    
    if("length" %in% colnames(subdomain_df)){
      plot_ortho <- domain.plotting(orderedOrthoDf,
                                    ortho,
                                    sep,
                                    label_archi_size,
                                    title_archi_size,
                                    min(subdomain_df$start),
                                    max(c(subdomain_df$end, subdomain_df$length)))
      plot_seed <- domain.plotting(orderedSeedDf,
                                   seed,
                                   sep,
                                   label_archi_size,
                                   title_archi_size,
                                   min(subdomain_df$start),
                                   max(c(subdomain_df$end, subdomain_df$length)))
      
    } else{
      plot_ortho <- domain.plotting(orderedOrthoDf,
                                    ortho,
                                    sep,
                                    label_archi_size,
                                    title_archi_size,
                                    min(subdomain_df$start),
                                    max(subdomain_df$end))
      plot_seed <- domain.plotting(orderedSeedDf,
                                   seed,
                                   sep,
                                   label_archi_size,
                                   title_archi_size,
                                   min(subdomain_df$start),
                                   max(subdomain_df$end))
    }
    
    # grid.arrange(plot_seed,plot_ortho,ncol=1)
    
    if(ortho == seed){
      arrangeGrob(plot_seed, ncol = 1)
    } else {
      seedHeight = length(levels(as.factor(orderedSeedDf$feature)))
      orthoHeight = length(levels(as.factor(orderedOrthoDf$feature)))
      
      arrangeGrob(plot_seed ,plot_ortho, ncol = 1, heights = c(seedHeight, orthoHeight))
    }
  }
}

gene_agePlot <- function(gene_age_text){
  countDf <- gene_ageDfMod()
  p <- ggplot(countDf,
              aes(fill = age, y = percentage, x = 1)) +
    geom_bar(stat = "identity") +
    scale_y_reverse() +
    coord_flip() +
    theme_minimal()
  p <- p + geom_text(data=countDf,
                     aes(x = 1, y = 100 - pos,
                         label = paste0(freq,"\n",percentage,"%")),
                     size = 4 * gene_age_text)
  p <- p + theme(legend.position = "bottom",
                 legend.title = element_blank(),
                 legend.text = element_text(size = 12 * gene_age_text),
                 axis.title = element_blank(),
                 axis.text = element_blank()) +
    scale_fill_brewer(palette="Spectral") +
    guides(fill = guide_legend(nrow = round(nrow(countDf) / 3, 0),
                               byrow = TRUE))
  p
}

# plot clustered profiles -----------------------------------------------------
dendrogram <- function(dd.col){
  py <- as.ggdend(dd.col)
  p <- ggplot(py, horiz = TRUE, theme=theme_minimal()) +
    theme(axis.title = element_blank(), axis.text.y = element_blank())
  p
}

# unused functions ============================================================

# # show FASTA sequence in popup windows of selected plot ---------------------
# getFasta <- function(file,seqID,groupID,faInput){
#   fasta <- ""
#   ### read file and get sequence
#   if(file.exists(file)){
#     fastaFile = readAAStringSet(file)
#     
#     seq_name = names(fastaFile)
#     sequence = paste(fastaFile)
#     fa <- data.frame(seq_name, sequence)
#     seq <- fa$sequence[pmatch(seqID,fa$seq_name)]
#     
#     if(length(seq[1]) < 1){
#       fasta <- paste0(seqID," not found in ",file,"! Please check id_format in FASTA config again!")
#     } else{
#       if(faInput == 1){
#         fasta <- paste(paste0(">",seqID),seq[1],sep="\n")
#       } else {
#         fasta <- paste(paste0(">",groupID,"|",seqID),seq[1],sep="\n")
#       }
#       
#     }
#   } else {
#     fasta <- paste0(file," not found! Please check the path and dir_format in FASTA config again!")
#   }
#   ### return
#   print(paste0(groupID,"|",seqID))
#   return(fasta)
# }

# # convert long to wide format -----------------------------------------------
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


# BASIC FUNCTIONS =============================================================
# reverse string --------------------------------------------------------------
strReverse <- function(x){
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")
}

# get last n characters from string x -----------------------------------------
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# check internet connection ---------------------------------------------------
hasInternet <- function(){
  !is.null(curl::nslookup("r-project.org", error = FALSE))
}


# NEWICK ======================================================================
# used to replace read.tree function of APE package ---------------------------
# when input tree has singletons  
# function to read a Newick string with node labels & (possible) singles
# written by Liam J. Revell 2013
read.newick<-function(file = "", text){
  # check to see if reading from file
  if(file != "") text <- scan(file, sep = "\n", what = "character")
  if(length(text) > 1){
    tree <- lapply(text, newick)
    class(tree) <- "multiPhylo"
  } else tree <- newick(text)
  return(tree)
}

# main Newick string function -------------------------------------------------
# written by Liam J. Revell 2013
newick <- function(text){
  text <- unlist(strsplit(text, NULL))
  tip.label <- vector(mode = "character")
  node.label <- vector(mode = "character") 
  edge <- matrix(c(1, NA), 1, 2) 
  edge.length <- vector()
  currnode <- 1
  Nnode <- currnode
  i <- j <- k <- 1
  while(text[i] != ";"){
    if(text[i] == "("){
      if(j > nrow(edge)) edge <- rbind(edge, c(NA, NA))
      edge[j, 1] <- currnode
      i <- i + 1
      # is the next element a label?
      if(is.na(match(text[i],
                     c("(", ")", ",", ":", ";")))){
        temp <- getLabel(text, i)
        tip.label[k] <- temp$label
        i <- temp$end
        edge[j, 2] <-- k
        k <- k + 1
        # is there a branch length?
        if(text[i] == ":"){
          temp <- getEdgeLength(text, i)
          edge.length[j] <- temp$edge.length
          i <- temp$end
        }	
      } else if(text[i] == "("){
        Nnode <- Nnode + 1 # creating a new internal node
        currnode <- Nnode
        edge[j, 2] <- currnode # move to new internal node
      }
      j <- j + 1
    } else if(text[i] == ")"){
      i <- i + 1
      # is the next element a label?
      if(is.na(match(text[i],
                     c("(", ")", ",", ":", ";")))){
        temp <- getLabel(text, i)
        node.label[currnode] <- temp$label
        i <- temp$end
      }
      # is there a branch length?
      if(text[i] == ":"){
        temp <- getEdgeLength(text, i)
        if(currnode > 1){ 
          ii <- match(currnode, edge[, 2])
          edge.length[ii] <- temp$edge.length
        } else root.edge <- temp$edge.length
        i <- temp$end
      }
      # move down the tree
      if(currnode > 1) currnode <- edge[match(currnode, edge[, 2]), 1]
    } else if(text[i] == ","){
      if(j > nrow(edge)) edge <- rbind(edge, c(NA, NA))
      edge[j, 1] <- currnode
      i <- i + 1
      # is the next element a label?
      if(is.na(match(text[i],
                     c("(", ")", ",", ":", ";")))){
        temp <- getLabel(text, i)
        tip.label[k] <- temp$label
        i <- temp$end
        edge[j, 2] <-- k
        k <- k + 1
        # is there a branch length?
        if(text[i] == ":"){
          temp <- getEdgeLength(text, i)
          edge.length[j] <- temp$edge.length
          i <- temp$end
        }
      } else if(text[i] == "("){
        Nnode <- Nnode + 1 # creating a new internal node
        currnode <- Nnode
        edge[j, 2] <- currnode # move to internal node
      }
      j <- j + 1
    }
  }
  Ntip <- k - 1
  edge[edge > 0] <- edge[edge > 0] + Ntip
  edge[edge < 0] <-- edge[edge < 0]
  edge.length[is.na(edge.length)] <- 0
  if(length(edge.length) == 0) edge.length <- NULL
  node.label[is.na(node.label)] <- ""
  if(length(node.label) == 0) node.label <- NULL
  # assemble into "phylo" object
  tree<-list(edge = edge,
             Nnode = as.integer(Nnode),
             tip.label = tip.label,
             edge.length = edge.length,
             node.label = node.label)
  class(tree) <- "phylo"
  return(tree)
}

# check the validity of input newick tree -------------------------------------
checkNewick <- function(filein, input_demo_data, input_main){
  tree <- read.table(file = filein$datapath,
                     header = F,
                     check.names = FALSE,
                     comment.char = "",
                     fill = F)
  
  # get tree structure
  treeStruc <- gsub(regex("\\w"), "", as.character(tree$V1))
  
  open = str_count(treeStruc, "\\(")
  close = str_count(treeStruc, "\\)")
  comma = str_count(treeStruc, "\\,")
  singleton = str_count(treeStruc, "\\(\\)")
  
  # return check
  if(is.null(input_main)){
    return(0) # don't check if main input is absent
  } else {
    if(singleton > 0){
      return(3) # tree contains singleton
    }
    
    if(open != close){
      return(1) # missing parenthesis
    } else {
      if(comma != (open+1)){
        # return(2) # missing comma
      } else {
        # get list of tips
        nodeString <- gsub(regex("\\W+"), "#", as.character(tree$V1))
        nodeList <- unlist(strsplit(nodeString, "#"))
        # list of input taxa
        inputTaxa <- subset_taxa(input_demo_data, input_main)
        
        missingTaxa <- list()
        j <- 1
        for(i in 1:length(nodeList)){
          if(nchar(nodeList[i]) > 0 & !(nodeList[i] %in% inputTaxa)){
            missingTaxa[[j]] <- nodeList[i]
            j <- j+1
          }
        }
        if(length(missingTaxa) > 0){
          # contains taxa that not exist in main input
          return(paste(missingTaxa, collapse="; ")) 
        } else {
          return(0)
        }
      }
    }
  }
  return(0)
}

# function gets label ---------------------------------------------------------
# written by Liam J. Revell 2011-2013
getLabel <- function(text, start, stop.char = c(",", ":", ")")){
  i <- 0
  label <- vector()
  while(is.na(match(text[i + start], stop.char))){
    label[i + 1] <- text[i + start]
    i <- i + 1
  }
  return(list(label = paste(label, collapse = ""), end = i + start))
}

# function gets branch length -------------------------------------------------
# written by Liam J. Revell 2011-2013
getEdgeLength <- function(text, start){
  i <- start + 1; m <- 1
  temp <- vector()
  while(is.na(match(text[i], c(",", ")", ";")))){
    temp[m] <- text[i]
    i <- i + 1
    m <- m + 1
  }
  return(list(edge.length = as.numeric(paste(temp, collapse="")), end = i))
}


# create rooted tree from a matrix --------------------------------------------
createRootedTree <- function(df,rootTaxon){
  # calculate distance matrix
  taxdis <- tryCatch(taxa2dist(df), error = function(e) e)
  
  # root tree 
  tree <- as.phylo(hclust(taxdis))
  tree <- root(tree, outgroup = rootTaxon, resolve.root = TRUE)
  
  # return
  return(tree)
}

# sort supertaxa list based on chosen reference taxon -----------------------
sortTaxaFromTree <- function(tree){
  # and get ordered taxon list
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  taxonList <- rev(tree$tip.label[ordered_tips])
  
  # return
  return(taxonList)
}

# unsorting function to keep user defined geneID order ----------------------
unsortID <- function(data,order){
  data$geneID <- as.factor(data$geneID)
  if(order == FALSE){
    # keep user defined geneID order
    data$geneID <- factor(data$geneID, levels = unique(data$geneID))  
  }
  return(data)
}



 
# FASTA OUTPUT ----------------------------------------------------------------
fastaOutData <- function(dataOut, filein, demo_data, in_type, one_seq_fasta, path, dir_format, file_ext, id_format){
  # dataOut <- as.data.frame(download_custom_data())
  fastaOutDf <- data.frame()
  
  ### check main input
  if(!is.null(filein)){input_type <- check_input_vadility(filein)}
  else{ input_type <- "NA"}
  
  ### get seqs for AMPK-TOR and microsporidia ONLINE demo data
  if(demo_data == "ampk-tor" | demo_data == "demo"){
    fastaUrl <- paste0("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/fasta_file/concatenatedSeq.fa")
    if(demo_data == "ampk-tor"){
      fastaUrl <- paste0("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.extended.fa")
    }
    
    if(url.exists(fastaUrl)){
      # load fasta file
      faFile <- as.data.frame(read.table(fastaUrl,
                                         sep = "\t",
                                         header = F,
                                         fill = T,
                                         stringsAsFactors = FALSE,
                                         quote = ""))
      faDf <- data.frame("seqID" = faFile$V1[grepl(">",faFile$V1)],
                         "seq" = faFile$V1[!grepl(">",faFile$V1)],
                         stringsAsFactors = FALSE)
      
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
  
  # get seqs for fasta main input -------------------------------------------
  if(input_type == "fasta"){
    file <- filein$datapath
    fastaFile = readAAStringSet(file)
    
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)
    # data frame contains all sequences from input file
    fa <- data.frame(seq_name, sequence)  
    
    for(j in 1:nrow(dataOut)){
      seqID <- paste0(as.character(dataOut$geneID[j]),
                      "|ncbi",
                      as.character(dataOut$ncbiID[j]),
                      "|",
                      as.character(dataOut$orthoID[j]))
      
      seq <- fa$sequence[pmatch(seqID,fa$seq_name)]
      
      if(length(seq[1]) < 1){
        fastaOut <- paste0(seqID,
                           " not found in ",
                           file,
                           "! Please check again!")
      } else{
        fastaOut <- paste(paste0(">", seqID), seq[1], sep = "\n")
      }
      
      fastaOutDf <- rbind(fastaOutDf,as.data.frame(fastaOut))
    }
  }
  
  ### get seqs for extended.fa
  if(demo_data == "none" & in_type == "Concatenated fasta file"){
    if(!is.null(one_seq_fasta)){
      fasIn <- one_seq_fasta
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
      if(input_type != "fasta"){
        fastaOut <- paste0("Please provide FASTA file(s) in Input & settings page!")
        fastaOutDf <- rbind(fastaOutDf,as.data.frame(fastaOut))
      }
    }
  }
  
  ### get seqs for other cases (input offline fasta files in a folder)
  if(demo_data == "none" & in_type == "Fasta folder" & input_type != "fasta"){
    if(path != ""){
     
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
        fastaOut <- paste0("No fasta file has been found in ",
                           path,
                           "!!! Please check the full path to FASTA folder and the id_format (header format) in FASTA config again!!!")
        fastaOutDf <- rbind(fastaOutDf,
                            as.data.frame(fastaOut))
      }
      
    } else {
      fastaOut <- paste0("Please provide FASTA files in Input & settings page!")
      fastaOutDf <- rbind(fastaOutDf,
                          as.data.frame(fastaOut))
    }
  }
  
  # remove duplicated sequences
  fastaOutDf <- fastaOutDf[!duplicated(fastaOutDf), ]
  
  return(fastaOutDf)
}



# GROUP COMPARISON ============================================================
# print list of available taxa ------------------------------------------------
taxa_select_gc <- function(rank_select_gc, input_demo_data, input_main){
  
  # if there is no rank set, there can not be any available taxa
  if (length(rank_select_gc) == 0) return()
  else{
    
    # load list of unsorted taxa
    dt <- get_taxa_list(TRUE, input_demo_data, input_main)
    
    
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

# Deciding which plots should be shown ----------------------------------------
get_plots_gc <- function(gene, plot_gc,
                         significant_genes_gc,
                         show_p_value,
                         highlight_significant,
                         significance,
                         var1_id, 
                         var2_id,
                         x_size_gc,
                         y_size_gc,
                         interesting_features,
                         angle_gc,
                         legend_size_gc,
                         ledgend_gc){
  plot_gc 
  if (is.null(significant_genes_gc)) return()
  else if (gene == "all"){
    get_plot_output_list(significant_genes_gc,
                         show_p_value,
                         highlight_significant,
                         significance,
                         var1_id,
                         var2_id,
                         x_size_gc,
                         y_size_gc,
                         interesting_features,
                         angle_gc,
                         legend_size_gc, 
                         legend)
  }else{
    x <- {
      significant_genes_gc[significant_genes_gc$geneID == gene, ]
    }
    if (nrow(x) == 0) return()
    get_plot_output_list(x, show_p_value,
                         highlight_significant,
                         significance,
                         var1_id,
                         var2_id,
                         x_size_gc,
                         y_size_gc,
                         interesting_features,
                         angle_gc,
                         legend_size_gc,
                         legend)
  }
}




# Generate the list with all plots --------------------------------------------
get_plot_output_list <- function(genes, show_p_value, highlight_significant, significance) {
  # if we dont have parameters we can not generate plots
  if (is.null(genes)) return()
  
  # Insert plot output objects the list
  plot_output_list <- lapply(1:nrow(genes), function(i) {
    plotname <- paste(genes[i, 1])
    plot_output_object <- renderPlot(get_multiplot(genes[i, ],
                                                   show_p_value,
                                                   highlight_significant,
                                                   significance,
                                                   var1_id,
                                                   var2_id,
                                                   x_size_gc,
                                                   y_size_gc,
                                                   interesting_features,
                                                   angle_gc,
                                                   legend_size_gc,
                                                   legend_gc),
                                     height = 650, width = 700)
  })
  do.call(tagList, plot_output_list) # needed to display properly.
  return(plot_output_list)
}


# Generating the plots ========================================================
# Put the plots for one spicific gene in one multiplot ------------------------
get_multiplot <- function(gene_info,
                          show_p_value,
                          highlight_significant,
                          significance,
                          var1_id,
                          var2_id,
                          x_size_gc,
                          y_size_gc,
                          interesting_features,
                          angle_gc,
                          legend_size_gc,
                          legend_gc){

  # Sorting the information to the selected gene
  gene <- as.character(gene_info$geneID)
  in_group <- as.data.frame(gene_info$in_group)
  out_group <- as.data.frame(gene_info$out_group)
  features <- as.data.frame(gene_info$features)
  pvalues <- gene_info$pvalues
  var <- gene_info$var
  
  # the barplot does not depent on the selected variable(s)
  barplot <-  get_barplot_gc(gene,
                             in_group,
                             out_group,
                             features, interesting_features,
                             x_size_gc, angle_gc, y_size_gc,
                             legend_size_gc, legend_gc)
  
  if (is.null(barplot)){
    barplot <- textGrob("The selected domains are not found in the gene")
  }
  # if both variables are selected there are going to be 2 boxplots
  if (var == "Both"){
    pvalues <- unlist(pvalues, recursive = FALSE)
    p1 <- unlist(pvalues[1])
    p2 <- unlist(pvalues[2])
    
    # Check if the p_values should be printed
    if (show_p_value == TRUE){
      info_p1 <- get_info_p_values(p1)
      info_p2 <- get_info_p_values(p2)
    }
    else{
      info_p1 <- " "
      info_p2 <- " "
    }
    
    # check if the significant plot should be highlighted
    if (highlight_significant == TRUE){
      if (is.null (p1[1])) c1 <- "grey"
      else if (p1[length(p1)] < significance) c1 <- "indianred2"
      else c1 <- "grey"
      
      if (is.null (p2[1])) c2 <- "grey"
      else if (p2[length(p2)] < significance) c2 <- "indianred2"
      else c2 <- "grey"
    }
    else{
      c1 <- "grey"
      c2 <- "grey"
    }
    
    boxplot1 <- get_boxplot_gc(in_group,
                               out_group,
                               var1_id,
                               gene,
                               c1,
                               info_p1,
                               var1_id,
                               var2_id,
                               x_size_gc,
                               y_size_gc)
    boxplot2 <- get_boxplot_gc(in_group,
                               out_group,
                               var2_id,
                               gene,
                               c2,
                               info_p2,
                               var1_id,
                               var2_id,
                               x_size_gc,
                               y_size_gc)
    
    m <- grid.arrange(textGrob(gene),
                      arrangeGrob(boxplot1, boxplot2, ncol = 2),
                      barplot,
                      heights = c(0.02, 0.45, 0.458), ncol = 1)
  }else {
    p <- unlist(pvalues)
    
    if (show_p_value == TRUE){
      info <- get_info_p_values(p)
    }else{
      info <- " "
    }
    
    boxplot <- get_boxplot_gc(in_group,
                              out_group,
                              var,
                              gene,
                              "grey",
                              info,
                              var1_id,
                              var2_id,
                              x_size_gc,
                              y_size_gc)
    
    m <- grid.arrange(textGrob(gene),
                      boxplot,
                      barplot,
                      heights = c(0.02, 0.45, 0.458), ncol = 1)
  }
  return(m)
}

# Create a Boxplot ------------------------------------------------------------
get_boxplot_gc <- function (in_group_df,
                            out_group_df,
                            var,
                            gene,
                            colour,
                            info,
                            var1_id,
                            var2_id,
                            x_size_gc,
                            y_size_gc){
   
  if (var == var1_id){
    in_g <- in_group_df$var1
    out_g <- out_group_df$var1
  }
  else if (var == var2_id) {
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
    theme(axis.text.x = element_text(size = x_size_gc, hjust = 1),
          axis.text.y = element_text(size = y_size_gc),
          axis.title.y = element_text(size = y_size_gc))
  boxplot_gc
}

# Create Barplot --------------------------------------------------------------
get_barplot_gc <- function(selected_gene,
                           in_group,
                           out_group,
                           features,
                           interesting_features,
                           x_size_gc,
                           angle_gc, 
                           y_size_gc, 
                           legend_size_gc,
                           legend_gc){
  
  subdomain_df <- features
  subdomain_df$feature <- as.character(subdomain_df$feature)
  # only show features that interest the user
  if (!("all" %in% interesting_features)){
    ifeatures <- NULL
    for (x in interesting_features){
      f <- subset(subdomain_df$feature,
                  startsWith(subdomain_df$feature, x))
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
      theme(axis.text.x = element_text(size = x_size_gc,
                                       angle = angle_gc, hjust = 1),
            axis.text.y = element_text(size = y_size_gc),
            axis.title.y = element_text(size = y_size_gc),
            legend.position = legend_gc,
            legend.text = element_text(size = legend_size_gc),
            legend.title = element_text(size = legend_size_gc))
    barplot_gc
  } else (return(NULL))
}

# Get the p_values to print under the plot ------------------------------------
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

# Get the list with all the significant genes and the dataset -----------------
get_significant_genes <- function (selected_in_group_gc,
                                   list_selected_genes_gc,
                                   rank_select,
                                   var_name_gc,
                                   in_select,
                                   use_common_anchestor,
                                   right_format_features,
                                   significance,
                                   var1,
                                   var2,
                                   var1_id,
                                   var2_id,
                                   demo_data,
                                   anno_choose,
                                   fileDomain_input,
                                   domainPath,
                                   input_demo_data,
                                   input_main,
                                   var1_aggregate_by,
                                   var2_aggregate_by,
                                   end_index,
                                   selected,
                                   listIn,
                                   st_index,
                                   ordering,
                                   gene_list_selected,
                                   v,
                                   treeIn){
  if (is.null(selected_in_group_gc)
      | length(list_selected_genes_gc) == 0) return()
  
  # load name List
  name_list <- get_name_list(TRUE, TRUE)
  
  # load list of unsorted taxa
  dt <- get_taxa_list(FALSE, input_demo_data, input_main)
  
  # Parameters that are identical for all genes
  in_group <-  selected_in_group_gc
  rank <- rank_select
  var <- var_name_gc
  
  # Updateing of the Input ==================================================
  # if there is more than one element in the in_group
  # we look at the next common anchstor
  
  if (use_common_anchestor == TRUE){
    ancestor <- get_common_ancestor(in_group, rank, name_list, dt)
    if (is.null(ancestor))return()
    in_group <- ancestor[1]
    rank <- ancestor[2]
  } else{
    rank <- substring(rank, 4)
  }
  
  
  # List with all significant genes
  significantGenes <- data.frame(
    geneID = character(),
    in_group = I(list()),
    out_group = I(list()),
    pvalues = I(list()),
    features = I(list()),
    databases = I(list()))
  
  # List of genes to look at
  data_full <- get_data_filtered(var1_aggregate_by, var2_aggregate_by,
                                 end_index,
                                 selected,
                                 listIn,
                                 st_index,
                                 demo_data,
                                 ordering,
                                 fileIn,
                                 gene_list_selected,
                                 sortedtaxa_list)
  if (is.element("all", list_selected_genes_gc)){
    genes <- data_full$geneID
    genes <- genes[!duplicated(genes)]
  } else {
    genes <- list_selected_genes_gc
  }
  genes <- sort(genes)
  
  # Subset depending on the rank and the in_group
  selected_subset <- get_selected_subset(rank, in_group, name_list, dt)
  selected_subset <- subset(selected_subset,
                            !selected_subset$fullName == in_select)
  
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
      subset(out_group_df, !out_group_df$fullName == in_select)
    }
    
    # Generate the p_values for the gene
    pvalues <- get_significant(in_group_df, out_group_df, var, gene, significance, var1, var2, var1_id, var2_id)
    
    if (!is.null(pvalues)){
      new_row <- data.frame(geneID = gene,
                            in_group = NA,
                            out_group = NA,
                            pvalues = NA,
                            features = NA)
      new_row$in_group <- list(in_group_df)
      new_row$out_group <- list(out_group_df)
      new_row$pvalues <- list(pvalues)
      features  <- get_features(gene,
                                demo_data,
                                anno_choose,
                                fileDomain_input,
                                domainPath)
      new_row$features <- list(features)
      if (right_format_features){
        new_row$databases <- list(get_prefix_features(features))
      }
      significantGenes <- rbind(significantGenes, new_row)
    }
  }
  if (nrow(significantGenes) != 0){
    significantGenes$var <- var
    significantGenes$rank <- rank
    significantGenes
  }
}

# Get the database for each feature in a specific gene ------------------------
# f is dataframe in $features
get_prefix_features <- function(f){
  features <- f$feature
  choices <- gsub("_.*", "", features)
  choices <- choices[!duplicated(choices)]
  choices
}

# Get the Subset depending on the choosen rank --------------------------------
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

# Generate the in_group -------------------------------------------------------
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


# Decide if the gene is significant -------------------------------------------
# if the gene is significant return the pvalues
get_significant <- function(in_g, out_g, var, gene, significance_level, var1, var2, var1_id, var2_id){

  
  if (var == "Both"){

    significant <-  FALSE
    # get the pValues for both variables
    pvalues1 <- get_p_values(in_g$var1, out_g$var1, significance_level)
    pvalues2 <- get_p_values(in_g$var2, out_g$var2, significance_level)
    
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
    if (var == var1_id){
      pvalues <- get_p_values(in_g$var1, out_g$var1, significance_level)
    }
    else {
      pvalues <- get_p_values(in_g$var2, out_g$var2, significance_level)
    }
    
    # Analog to getting the significance with both variables
    if (is.null(pvalues)) return(NULL)
    else if (is.nan(pvalues[length(pvalues)])){
      return(NULL)
    }else if (pvalues[length(pvalues)] < significance_level){
      pvalues <-  list(pvalues)
      return(pvalues)
    } else{
      return(NULL)
    }
  }
}

# calculate the p_values ------------------------------------------------------
get_p_values <- function(var_in, var_out, significance_level){
  # upper limit for the probability to reject H0 if it is correct
  
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

# get the list with all the features in the gene ------------------------------
get_features <- function(selected_gene, demo_data, anno_choose, file_domain_input, domain_path){
  # parse domain file
  file_domain <- get_domain_file_gc(selected_gene, demo_data)
  
  if (demo_data == "demo" | demo_data == "ampk-tor"){
    domain_df <- as.data.frame(read.csv(file_domain,
                                        sep = "\t",
                                        header = F,
                                        comment.char = "",
                                        stringsAsFactors = FALSE,
                                        quote = ""))
    
    if (ncol(domain_df) == 5){
      colnames(domain_df) <- c("seedID",
                               "orthoID",
                               "feature",
                               "start",
                               "end")
    } else if (ncol(domain_df) == 6){
      colnames(domain_df) <- c("seedID",
                               "orthoID",
                               "feature",
                               "start",
                               "end",
                               "weight")
    } else if (ncol(domain_df) == 7){
      colnames(domain_df) <- c("seedID",
                               "orthoID",
                               "feature",
                               "start",
                               "end",
                               "weight",
                               "path")
    }
    
    domain_df$length <- max(domain_df$end)
    
  } else {
    if (file_domain != FALSE){
      domain_df <- as.data.frame(read.table(file_domain,
                                            sep = "\t",
                                            header = FALSE,
                                            comment.char = ""))
    }
    if (ncol(domain_df) == 6){
      colnames(domain_df) <- c("seedID",
                               "orthoID",
                               "length",
                               "feature",
                               "start",
                               "end")
      
    } else if (ncol(domain_df) == 7){
      colnames(domain_df) <- c("seedID",
                               "orthoID",
                               "length",
                               "feature",
                               "start",
                               "end",
                               "weight")
      
    } else if (ncol(domain_df) == 8){
      colnames(domain_df) <- c("seedID",
                               "orthoID",
                               "length",
                               "feature",
                               "start",
                               "end",
                               "weight",
                               "path")
    }
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

# get the data where to find the features -------------------------------------
get_domain_file_gc <- function(group, demo_data, anno_choose, file_domain_input, domain_path){
  # domain file
  if (demo_data == "demo" | demo_data == "ampk-tor"){
    updateButton(session, "do_domain_plot", disabled = FALSE)
    if (demo_data == "demo"){
      file_domain <- {
        suppressWarnings(paste0("https://github.com/BIONF/phyloprofile-data/blob/master/demo/domain_files/",
                                group,
                                ".domains?raw=true"))
      }
    } else {
      file_domain <- {
        suppressWarnings(paste0("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.domains_F"))
      }
    }
    
  }else {
    if (anno_choose == "from file"){
      file_domain <- file_domain_input
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
          file_domain <- paste0(domain_path,
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

get_multiplot_download_gc <- function(gene,
                                      show_p_value,
                                      highlight_significant,
                                      significance,
                                      var1_id,
                                      var2_id,
                                      x_size_gc,
                                      y_size_gc,
                                      interesting_features,
                                      angle_gc,
                                      legend_size_gc, 
                                      legend_gc){
  x <- subset(significant_genes_gc,
              significant_genes_gc$geneID == gene)
  get_multiplot(x,
                show_p_value,
                highlight_significant,
                significance,
                var1_id, var2_id,
                x_size_gc,
                y_size_gc,
                interesting_features,
                angle_gc,
                legend_size_gc, legend_gc)
}


# ===========================================================================
# Essential functions =======================================================
# ===========================================================================

# get list with all taxanomy ranks --------------------------------------------
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

# Get name list ("data/taxonNamesReduced.txt") --------------------------------
get_name_list <- function (as_character, delete_duplicated){
  name_list <- as.data.frame(read.table("data/taxonNamesReduced.txt",
                                        sep = "\t",
                                        header = T,
                                        fill = TRUE))
  if (as_character) {
    name_list$fullName <- as.character(name_list$fullName)
    name_list$rank <- as.character(name_list$rank)
  }
  if (delete_duplicated){
    name_list <- name_list[!duplicated(name_list), ]
  }
  
  return(name_list)
}

# Get taxa list ("data/taxonomyMatrix.txt") -----------------------------------
get_taxa_list <- function(subset_taxa_check, input_demo_data, input_main ){
  dt <- as.data.frame(read.table("data/taxonomyMatrix.txt",
                                 sep = "\t",
                                 header = T,
                                 stringsAsFactors = T))
  if (subset_taxa_check){
    dt <- dt[dt$abbrName  %in% subset_taxa(input_demo_data, input_main), ]
  }
  return(dt)
}

# check validity of main input file -------------------------------------------
check_input_vadility <- function(filein){
  inputDt <- as.data.frame(read.table(file = filein$datapath,
                                      sep = "\t",
                                      header = F,
                                      check.names = FALSE,
                                      comment.char = "",
                                      fill = T))
  # print(head(inputDt))
  if (is.na(inputDt[1, ncol(inputDt)])){
    return("moreCol")
  } else {
    names(inputDt) <- as.character(unlist(inputDt[1, ]))
    
    # XML -------------------------------------------------------------------
    # a XML file always starts with <?xml
    if (grepl("<?xml", colnames(inputDt)[1])){
      return("xml")
      
      # FASTA ---------------------------------------------------------------
      # a faste file always starts with > (start symbol of the header)
    } else if (grepl(">", colnames(inputDt)[1]) == TRUE){
      return("fasta")
      
    } else {
      # long and wide files always start with geneID
      if (grepl("geneID", colnames(inputDt)[1])){
        # LONG --------------------------------------------------------------
        if (is.na(pmatch("ncbi", colnames(inputDt)[3])) ||
            is.na(pmatch("ncbi", colnames(inputDt)[4])) ||
            is.na(pmatch("ncbi", colnames(inputDt)[5]))){
          return("long")
          
          # WIDE ------------------------------------------------------------
        } else {
          tmp <- inputDt[inputDt == ""][1]
          if ( !is.na(tmp) & tmp == ""){
            return("emptyCell")
          } else {
            return("wide")
          }
        }
      } else {
        # # OMA -------------------------------------------------------------
        # # A list of OMA protein ID's do not start with a gene ID
        # # The list has only one column (is that enough to identify it?)
        # if (ncol(inputDt) == 1) {
        #   return ("oma")
        #
        # # When it has no gene ID and is not a list of oma protein IDs -----
        # } else {
        #   return("noGeneID")
        #   }
        return("noGeneID")
        
      }
    }
  }
}


# REDUCE DATA FROM SPECIES LEVEL TO SUPERTAXA LEVEL -------------------------
# this data set contain only supertaxa
# and their value (%present, mVar1 & mVar2) for each gene
dataSupertaxa <- function(var1_aggregate_by, var2_aggregate_by,
                          end_index,
                          selected,
                          listIn,
                          st_index,
                          demo_data,
                          ordering,
                          fileIn,
                          gene_list_selected,
                          v,
                          rank_select,
                          in_select,
                          treeIn){

  fullMdData <- get_data_filtered(var1_aggregate_by, var2_aggregate_by,
                                  end_index,
                                  selected,
                                  listIn,
                                  st_index,
                                  demo_data,
                                  ordering,
                                  fileIn,
                                  gene_list_selected,
                                  sortedtaxa_list)
  
  # to check if working with the lowest taxonomy rank; 1 for NO; 0 for YES
  flag <- 1 
  if(length(unique(levels(as.factor(fullMdData$numberSpec)))) == 1){
    if(unique(levels(as.factor(fullMdData$numberSpec))) == 1){
      superDfExt <- fullMdData[,c("geneID",
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
  
  if(flag == 1){
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
    
    superDfExt <- merge(superDf,mOrthoID, by=c("geneID", "supertaxon"),
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
    names(superDfExt)[names(superDfExt)=="mVar1"] <- "var1"
    names(superDfExt)[names(superDfExt)=="mVar2"] <- "var2"
  }
  
  # print(superDfExt[superDfExt$geneID == "ampk_ACACB"
  # & superDfExt$supertaxon == "1001_Chordata",])
  # print("END2222")
  
  return(superDfExt)
}

# get all information for input data ----------------------------------------
get_data_filtered <- function(var1_aggregate_by, var2_aggregate_by,
                              end_index,
                              selected,
                              listIn,
                              st_index,
                              demo_data,
                              ordering,
                              fileIn,
                              gene_list_selected,
                              sortedtaxa_list) {

  mdData <- preData(end_index,
                    selected,
                    listIn,
                    st_index,
                    demo_data,
                    ordering,
                    fileIn,
                    gene_list_selected)
  
  
  # count number of inparalogs
  paralogCount <- plyr::count(mdData, c("geneID", "ncbiID"))
  mdData <- merge(mdData,
                  paralogCount,
                  by = c("geneID", "ncbiID"))
  colnames(mdData)[ncol(mdData)] <- "paralog"
  
  # (3) GET SORTED TAXONOMY LIST (3) 
  taxa_list <- sortedtaxa_list
  print(head(taxa_list))
  # calculate frequency of all supertaxa
  taxaCount <- plyr::count(taxa_list,"supertaxon")
  
  # merge mdData, mdDataVar2 and taxa_list to get taxonomy info
  taxaMdData <- merge(mdData, taxa_list, by = "ncbiID")
  taxaMdData$var1 <- {
    suppressWarnings(as.numeric(as.character(taxaMdData$var1)))
  }
  taxaMdData$var2 <- {
    suppressWarnings(as.numeric(as.character(taxaMdData$var2)))
  }
  
  # (4) calculate PERCENTAGE of PRESENT SPECIES (4)
  finalPresSpecDt <- calcPresSpec(taxaMdData, taxaCount)
  
  # (5) calculate max/min/mean/median VAR1 for every supertaxon of each gene (5)
  # remove NA rows from taxaMdData
  taxaMdDataNoNA <- taxaMdData[!is.na(taxaMdData$var1),]
  # calculate m var1
  mVar1Dt <- aggregate(taxaMdDataNoNA[, "var1"],
                       list(taxaMdDataNoNA$supertaxon,
                            taxaMdDataNoNA$geneID),
                       FUN = var1_aggregate_by)
  colnames(mVar1Dt) <- c("supertaxon", "geneID", "mVar1")
  
  # (6) calculate max/min/mean/median VAR2 for each super taxon (6)
  # remove NA rows from taxaMdData
  taxaMdDataNoNA_var2 <- taxaMdData[!is.na(taxaMdData$var2), ]
  # calculate max/min/mean/median VAR2
  if(nrow(taxaMdDataNoNA_var2) > 0){
    mVar2Dt <- aggregate(taxaMdDataNoNA_var2[, "var2"],
                         list(taxaMdDataNoNA_var2$supertaxon,
                              taxaMdDataNoNA_var2$geneID),
                         FUN = var2_aggregate_by)
    colnames(mVar2Dt) <- c("supertaxon","geneID","mVar2")
  } else {
    mVar2Dt <- taxaMdData[,c("supertaxon","geneID")]
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
                      by=c("geneID", "supertaxon"),
                      all.x = TRUE)
  fullMdData <- merge(fullMdData,
                      taxaCount,by=("supertaxon"),
                      all.x = TRUE)
  # rename "freq" into "numberSpec"
  names(fullMdData)[names(fullMdData)=="freq"] <- "numberSpec"
  
  fullMdData$fullName <- as.vector(fullMdData$fullName)
  names(fullMdData)[names(fullMdData)=="orthoID.x"] <- "orthoID"
  # parsed input data frame !!!
  fullMdData <- fullMdData[!duplicated(fullMdData), ] 
  
  return(fullMdData)
}

# subset data ---------------------------------------------------------------
preData <- function(end_index,
                    selected,
                    listIn,
                    st_index,
                    demo_data,
                    ordering,
                    fileIn,
                    gene_list_selected){
  # get list of gene of interest (from a separated file)
  listGene <- list()
  if(is.na(end_index)){end_index <- 30}
  if(gene_list_selected == "from file"){
    if(!is.null(listIn)){
      list <- as.data.frame(read.table(file = listIn$datapath,
                                       header = FALSE))
      listGeneOri <- list$V1
      if(st_index <= length(listGeneOri)){
        listGene <- listGeneOri[listGeneOri[st_index:end_index]]
      } else {
        listGene <- listGeneOri
      }
    }
  }
  
  # if(input$demo == TRUE){
  if(demo_data == "demo" | demo_data == "ampk-tor"){
    if(demo_data == "demo"){
      inputDf <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.long",
                            sep = "\t",
                            header = T,
                            fill = T,
                            stringsAsFactors = FALSE)
    } else {
      inputDf <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile",
                            sep = "\t",
                            header = T,
                            fill = T,
                            stringsAsFactors = FALSE)
    }
    inputDf <- unsortID(inputDf, ordering)
    
    if(length(listGene) >= 1){
      data <- inputDf[inputDf$geneID %in% listGene, ]
    } else {
      subsetID <- levels(as.factor(inputDf$geneID))[st_index:end_index]
      data <- inputDf[inputDf$geneID %in% subsetID, ]
    }
    
    if(ncol(data) < 5){
      for(i in 1:(5-ncol(data))){
        data[paste0("newVar", i)] <- 1
      }
    }
    colnames(data) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
  } else {
    if(is.null(filein)){return()}
    
    input_type <- check_input_vadility(filein)
    # XML -------------------------------------------------------------------
    if(input_type == "xml"){
      longDf <- xmlParser(filein$datapath)
      longDf <- unsortID(longDf,input$ordering)
      
      if(length(listGene) >= 1){
        data <- longDf[longDf$geneID %in% listGene, ]
      } else {
        subsetID <- levels(longDf$geneID)[st_index:end_index]
        data <- longDf[longDf$geneID %in% subsetID, ]
      }
      
      if(ncol(data) < 5){
        for(i in 1:(5-ncol(data))){
          data[paste0("newVar", i)] <- 1
        }
      }
      colnames(data) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
      
      # fasta ----------------------------------------------------------------- 
    } else if(input_type == "fasta"){
      longDf <- fastaParser(filein$datapath)
      longDf <- unsortID(longDf, ordering)
      
      if(length(listGene) >= 1){
        data <- longDf[longDf$geneID %in% listGene, ]
      } else {
        subsetID <- levels(longDf$geneID)[st_index:end_index]
        data <- longDf[longDf$geneID %in% subsetID, ]
      }
      
      if(ncol(data) < 5){
        for(i in 1:(5-ncol(data))){
          data[paste0("newVar", i)] <- 1
        }
      }
      colnames(data) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
      
      # long ------------------------------------------------------------------
    } else if(input_type == "long"){
      inputDf <- as.data.frame(read.table(file = filein$datapath,
                                          sep = "\t",
                                          header = T,
                                          check.names = FALSE, 
                                          comment.char = ""))
      inputDf <- unsortID(inputDf, ordering)
      
      if(length(listGene) >= 1){
        data <- inputDf[inputDf$geneID %in% listGene,]
      } else {
        subsetID <- levels(inputDf$geneID)[st_index:end_index]
        data <- inputDf[inputDf$geneID %in% subsetID,]
      }
      
      if(ncol(data) < 5){
        for(i in 1:(5-ncol(data))){
          data[paste0("newVar",i)] <- 1
        }
      }
      colnames(data) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
      
      # wide ------------------------------------------------------------------
    } else if (input_type == "wide"){
      inputDf <- as.data.frame(read.table(file = filein$datapath,
                                          sep = "\t",
                                          header = T,
                                          check.names = FALSE,
                                          comment.char = ""))
      mdData <- melt(inputDf, id = "geneID")
      mdData <- unsortID(mdData, ordering)
      
      if(length(listGene) >= 1){
        mdData <- mdData[mdData$geneID %in% listGene,]
      } else {
        subsetID <- levels(mdData$geneID)[st_index:end_index]
        mdData <- mdData[mdData$geneID %in% subsetID,]
      }
      
      splitCol <- data.frame(do.call("rbind",
                                     strsplit(as.character(mdData$value),
                                              "#",
                                              fixed = TRUE)))
      data <- cbind(mdData[, c("geneID", "variable")], splitCol)
      
      if(ncol(data) < 5){
        for(i in 1:(5-ncol(data))){
          data[paste0("newVar", i)] <- 0
        }
      }
      colnames(data) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
      
      # # OMA -------------------------------------------------------------------
      # } else if (input_type == "oma"){
      #   # What exactally needs to happen here 
      
      # None of the definded types --------------------------------------------
    } else {
      data <- data.frame("geneID" = character(),
                         "ncbiID" = character(),
                         "orthoID" = character(),
                         "var1" = character(),
                         "var2" = character(),
                         stringsAsFactors = F)
    }
  }
  
  # return preData
  return(data)
}



# Get list of all (super)taxa -------------------------------------------------
alltaxa_list <- function(filein, demo_data, rank_select){
  # input$mainInput, input$demo_data, input$rank_select

  
  # if(is.null(filein) & input$demo == FALSE) return()
  if(is.null(filein) & demo_data == "none") return()

  if(rank_select == "") return()
  
  if(length(unkTaxa(demo_data, filein)) > 0) return()
  
  # load list of unsorted taxa 
  Dt <- get_taxa_list(TRUE, demo_data, filein)
  
  
  # load list of taxon name 
  nameList <- get_name_list(TRUE, FALSE)
  
  # get rank name from rank_select
  rankName = substr(rank_select, 4, nchar(rank_select))
  
  # get rank number (number of column in unsorted taxa list - dataframe Dt)
  #    rankNr = 0 + as.numeric(substr(rank_select,1,2))  
  
  choice <- as.data.frame
  choice <- rbind(Dt[rankName])
  colnames(choice) <- "ncbiID"
  choice <- merge(choice,
                  nameList,
                  by = "ncbiID",
                  all = FALSE)
}


# check if there is any "unknown" taxon in input matrix -----------------------
unkTaxa <- function(demo_data, filein) {

  # get list of input taxa (from main input file)
  # if(input$demo == TRUE){
  if(demo_data == "demo"){
    data <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.wide",
                       sep = "\t",
                       header = T,
                       fill = T,
                       stringsAsFactors = FALSE)
    inputTaxa <- colnames(data)
  } else if(demo_data == "ampk-tor"){
    data <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile",
                       sep = "\t",
                       header = T,
                       fill = T,
                       check.names = FALSE,
                       comment.char = "")
    inputTaxa <- levels(data$ncbiID)
    
  } else {
    if(is.null(filein)) return() 
    
    input_type <- check_input_vadility(filein)
    
    # XML -------------------------------------------------------------------
    if(input_type == "xml"){
      longDf <- xmlParser(filein$datapath)
      inputTaxa <- levels(longDf$ncbiID)
      
      # FASTA -----------------------------------------------------------------
    } else if(input_type == "fasta"){
      longDf <- fastaParser(filein$datapath)
      inputTaxa <- levels(longDf$ncbiID)
      
      # LONG ------------------------------------------------------------------
    } else if(input_type == "long"){
      inputDf <- as.data.frame(read.table(file = filein$datapath,
                                          sep = "\t",
                                          header = T,
                                          check.names = FALSE,
                                          comment.char = ""))
      inputTaxa <- levels(inputDf$ncbiID)
      
      # WIDE ------------------------------------------------------------------
    } else if(input_type == "wide"){
      inputTaxa <- readLines(filein$datapath, n = 1)
      
      # # OMA -------------------------------------------------------------------
      # } else if (input_type == "oma"){
      #   long_df <- oma_parser( filein$datapath, input$selected_oma_type)
      #   print(head(long_df))
      #   inputTaxa <- levels(long_df$ncbiID)
      
      # None of the definded types --------------------------------------------
    }else {
      inputTaxa <- c("NA")
    }
  }
  
  if(inputTaxa[1] == "NA"){
    return()
  } else {
    inputTaxa <- unlist(strsplit(inputTaxa, split = "\t"))
    if(inputTaxa[1] == "geneID"){
      # remove "geneID" element from vector inputTaxa
      inputTaxa <- inputTaxa[-1]  
    }
    
    if(!file.exists(isolate({"data/rankList.txt"}))){
      return(inputTaxa)
    } else {
      info = file.info("data/rankList.txt")
      if(info$size == 0){
        return(inputTaxa)
      } else {
        # get list of all available taxon (from /data/rankList.txt)
        pipeCmd <- paste0("cut -f 1 ", getwd(), "/data/rankList.txt")
        allTaxa <- unlist((read.table(pipe(pipeCmd))))
        
        # list of unknown taxa
        unkTaxa <- inputTaxa[!(inputTaxa %in% allTaxa)]
        # return list of unkTaxa
        return(unkTaxa)
      }
    }
  }
}


# data for detailed plot ------------------------------------------------------
detail_plotDt <- function(v, input_tabs,
                          var1_aggregate_by,
                          var2_aggregate_by,
                          end_index,
                          selected,
                          listIn,
                          st_index,
                          demo_data,
                          ordering,
                          mainInput,
                          gene_list_selected,
                          inputTree,
                          percent,
                          var1, var2,
                          x_axis,
                          sortedtaxa_list){
  if (v$doPlot == FALSE) return()
  
  ##### GET INFO BASED ON CURRENT TAB
  if(input_tabs == "Main profile"){
    info <- mainpoint_info(v, rank_select, percent,
                           var1, var2,
                           mainInput, demo_data,
                           var1_relation, var2_relation,
                           var1_aggregate_by, var2_aggregate_by,
                           end_index,
                           gene_list_selected,
                           list,
                           st_index,
                           ordering,
                           inputTree,
                           x_axis,
                           selected)  # info = groupID,orthoID,supertaxon,mVar1,%spec,var2
  } else if(input_tabs=="Customized profile"){
    info <- selectedpoint_info()
  }
  
  if(is.null(info)) return()
  else{
    ### get info for present taxa in selected supertaxon (1)
    plotTaxon = info[3]
    plotGeneID = info[1]
    fullDf <- get_data_filtered(var1_aggregate_by, var2_aggregate_by,
                                end_index,
                                selected,
                                listIn,
                                st_index,
                                demo_data,
                                ordering,
                                fileIn,
                                gene_list_selected,
                                sortedtaxa_list)
    selDf <- as.data.frame(fullDf[fullDf$geneID == plotGeneID
                                  & fullDf$supertaxon == plotTaxon,])
    
    ### get all taxa of this supertaxon (2)
    allTaxaDf <- sortedtaxa_list
    allTaxaDf <- allTaxaDf[allTaxaDf$supertaxon == plotTaxon,]
    allTaxaDf <- subset(allTaxaDf, select = c("abbrName","fullName"))
    
    ### merge (1) and (2) together
    joinedDf <- merge(selDf,
                      allTaxaDf,
                      by =  c("abbrName"),
                      all.y = TRUE)
    joinedDf <- subset(joinedDf, select = c("abbrName",
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
    if(input$detailed_remove_na == TRUE){
      joinedDf <- joinedDf[!is.na(joinedDf$orthoID),]
    }
    
    ### return data for detailed plot
    return(joinedDf)
  }
}

# get input taxa --------------------------------------------------------------
subset_taxa <- function(input_demo_data, filein){
  
  
  # if(input$demo == TRUE){
  if(input_demo_data == "demo"){
    data <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/demo/test.main.wide",
                       sep = "\t",
                       header = T,
                       fill = T,
                       stringsAsFactors = FALSE)
    inputTaxa <- colnames(data)
  } else if(input_demo_data == "ampk-tor"){
    data <- read.table("https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/expTestData/ampk-tor/ampk-tor.phyloprofile",
                       sep = "\t",
                       header = T,
                       fill = T,
                       check.names = FALSE,
                       comment.char = "")
    inputTaxa <- levels(data$ncbiID)
  } else {
    if(is.null(filein)) return()
    else{
      
      if(length(unkTaxa(input_demo_data, filein)) == 0){
        # get list of input taxa (from main input file)
        input_type <- check_input_vadility(filein)
        
        # XML 
        if(input_type == "xml"){
          longDf <- xmlParser(filein$datapath)
          inputTaxa <- levels(longDf$ncbiID)
          
          # FASTA 
        } else if(input_type == "fasta"){
          longDf <- fastaParser(filein$datapath)
          inputTaxa <- levels(longDf$ncbiID)
          
          # LONG 
        } else if(input_type == "long"){
          inputDf <- as.data.frame(read.table(file = filein$datapath,
                                              sep = "\t",
                                              header = T,
                                              check.names = FALSE,
                                              comment.char = ""))
          inputTaxa <- levels(inputDf$ncbiID)
          
          # WIDE -
        } else if(input_type == "wide"){
          inputTaxa <- readLines(filein$datapath, n = 1)
          
          # # OMA 
          # } else if (input_type == "oma"){
          #   input$get_data_oma
          #   
          #   isolate({
          #     long_df <- oma_parser( filein$datapath, selected_oma_type)
          #     inputTaxa <- levels(longDf$ncbiID)
          #   })
          
          
          # None of the definded types  
        } else {
          inputTaxa <- "NA"
        }
      } else {
        inputTaxa <- readLines(filein$datapath, n = 1)
      }
      
      inputTaxa <- unlist(strsplit(inputTaxa,split = "\t"))
      if(inputTaxa[1] == "geneID"){
        # remove "geneID" element from vector inputTaxa
        inputTaxa <- inputTaxa[-1]   
      }
      
      # return input taxa
      return(inputTaxa)
    }
  }
}

# Function for clustering profiles --------------------------------------------
clusteredGeneList <- function(data,distMethod,clusterMethod){
  
  # do clustering
  row.order <- hclust(dist(data, method = distMethod),
                      method = clusterMethod)$order
  col.order <- hclust(dist(t(data), method = distMethod),
                      method = clusterMethod)$order
  
  # re-order data accoring to clustering
  dat_new <- data[row.order, col.order]
  
  # return clustered gene ID list
  clusteredGeneIDs <- as.factor(row.names(dat_new))
  clusteredGeneIDs
}

# calculate percentage of present species -------------------------------------
calcPresSpec <- function(taxaMdData, taxaCount){
  # taxaMdData = df("geneID",
  #                 "ncbiID",
  #                 "orthoID",
  #                 "var1",
  #                 "var2",
  #                 "paralog",
  #                 ....,
  #                 "supertaxon")
  taxaMdData <- taxaMdData[taxaMdData$orthoID != "NA",]
  
  # get geneID and supertaxon
  geneIDsupertaxon <- subset(taxaMdData, 
                             select = c('geneID',
                                        'supertaxon',
                                        'paralog',
                                        'abbrName'))
  # remove duplicated rows
  geneIDsupertaxon <- geneIDsupertaxon[!duplicated(geneIDsupertaxon), ] 
  
  # remove NA rows from taxaMdData
  taxaMdDataNoNA <- taxaMdData[taxaMdData$orthoID != "NA", ]
  
  # count present frequency of supertaxon for each gene
  geneSupertaxonCount <- plyr::count(taxaMdDataNoNA,
                                     c('geneID', 'supertaxon'))
  
  # merge with taxaCount to get total number of species of each supertaxon
  # and calculate presSpec
  presSpecDt <- merge(geneSupertaxonCount,
                      taxaCount,
                      by = 'supertaxon',
                      all.x = TRUE)
  
  specCount <- plyr::count(geneIDsupertaxon,c('geneID',
                                              'supertaxon'))
  presSpecDt <- merge(presSpecDt,
                      specCount,by = c('geneID',
                                       'supertaxon'))
  
  presSpecDt$presSpec <- presSpecDt$freq / presSpecDt$freq.y
  
  presSpecDt <- presSpecDt[presSpecDt$presSpec <= 1, ]
  presSpecDt <- presSpecDt[order(presSpecDt$geneID), ]
  presSpecDt <- presSpecDt[, c("geneID", "supertaxon", "presSpec")]
  
  # add absent supertaxon into presSpecDt
  geneIDsupertaxon <- subset(geneIDsupertaxon,
                             select = -c(paralog,abbrName))
  finalPresSpecDt <- merge(presSpecDt,
                           geneIDsupertaxon,
                           by = c('geneID', 'supertaxon'),
                           all.y = TRUE)
  finalPresSpecDt$presSpec[is.na(finalPresSpecDt$presSpec)] <- 0
  
  # remove duplicated rows
  finalPresSpecDt <- finalPresSpecDt[!duplicated(finalPresSpecDt), ]
  
  # return finalPresSpecDt
  return(finalPresSpecDt)
}

# sort one domain dataframe (ortho) based on the other domain Df (seed) -------
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
  orderedOrthoDf$feature <- factor(orderedOrthoDf$feature,
                                   levels = unique(orderedOrthoDf$feature))
  
  #return sorted df
  orderedOrthoDf
}


# get info clicked point on main heatmap --------------------------------------
mainpoint_info <- function(v, rank_select, in_select, percent,
                           var1, var2,
                           mainInput, demo_data,
                           var1_relation, var2_relation,
                           var1_aggregate_by, var2_aggregate_by,
                           end_index,
                           gene_list_selected,
                           list,
                           st_index,
                           ordering,
                           inputTree,
                           x_axis, 
                           selected){
  # check input
  if (v$doPlot == FALSE) return()
  # get selected supertaxon name
  taxa_list <- as.data.frame(read.table("data/taxonNamesReduced.txt",
                                        sep = "\t",
                                        header = T))
  rankName = substr(rank_select, 4, nchar(rank_select))
  in_select <- {
    as.numeric(taxa_list$ncbiID[taxa_list$fullName == in_select])
  }
  
  dataHeat <- dataHeat(percent,
                       var1,var2,
                       mainInput, demo_data,
                       var1_relation, var2_relation,
                       var1_aggregate_by, var2_aggregate_by,
                       end_index,
                       selected,
                       list,
                       st_index,
                       ordering,
                       gene_list_selected,
                       v,
                       rank_select,
                       in_select,
                       inputTree)
  if(input$apply_cluster == TRUE){
    dataHeat <- clusteredDataHeat()
  }
  
  # get values
  if (is.null(input$plot_click$x)) return()
  else{
    # get cooridiate point
    if(x_axis == "genes"){
      corX = round(input$plot_click$y);
      corY = round(input$plot_click$x)
    } else {
      corX = round(input$plot_click$x);
      corY = round(input$plot_click$y)
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
    if(!is.na(dataHeat$var1[dataHeat$geneID == geneID
                            & dataHeat$supertaxon == spec][1])){
      var1 <- max(na.omit(dataHeat$var1[dataHeat$geneID == geneID
                                        & dataHeat$supertaxon == spec]))
    }
    Percent <- NA
    if(!is.na(dataHeat$presSpec[dataHeat$geneID == geneID 
                                & dataHeat$supertaxon == spec][1])){
      Percent <- max(na.omit(dataHeat$presSpec[dataHeat$geneID == geneID
                                               & dataHeat$supertaxon == spec]))
    }
    var2 <- NA
    if(!is.na(dataHeat$var2[dataHeat$geneID == geneID 
                            & dataHeat$supertaxon == spec][1])){
      var2 <- max(na.omit(dataHeat$var2[dataHeat$geneID == geneID
                                        & dataHeat$supertaxon == spec]))
    }
    
    # get ortholog ID
    orthoID <- dataHeat$orthoID[dataHeat$geneID == geneID
                                & dataHeat$supertaxon == spec]
    if(length(orthoID) > 1){
      orthoID = paste0(orthoID[1],",...")
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
    if(is.na(as.numeric(Percent))){return()}
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
}