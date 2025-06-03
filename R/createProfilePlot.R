#' Process ortholog IDs
#' @description Process ortholog IDs to identify duplicated IDs
#' @param dataHeat a data frame contains processed profiles (see
#' ?fullProcessedProfile, ?filterProfileData)
#' @return the same dataframe as input, but the ortholog IDs are changed into
#' <taxID:orthoID>. New column "orthoFreq" specifies if the ortholog IDs are
#' single or duplicated
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @export
#' @examples
#' ?processOrthoID
#' data("finalProcessedProfile", package="PhyloProfile")
#' processOrthoID(finalProcessedProfile)

processOrthoID <- function(dataHeat = NULL) {
    if (is.null(dataHeat)) stop("Input data cannot be NULL!")
    orthoID <- geneID <- orthoFreqCount <- orthoFreqNew <- supertaxonID <- NULL
    orthoIDNew <- orthoFreq <- NULL
    # Check if idFormat is "bionf"
    idFormat <- "other"
    firstOrtho <- strsplit(dataHeat$orthoID[1],'\\|')[[1]]
    firstGeneID <- dataHeat$geneID[1]
    firstSupertaxonID <- dataHeat$supertaxonID[1]
    if (
        length(firstOrtho) >= 3 && firstOrtho[1] == firstGeneID && 
        grepl(firstSupertaxonID, firstOrtho[2])
    ) { idFormat <- "bionf" }
    
    # parse orthoID
    if (idFormat == "bionf") {
        dataHeat <- within(
            dataHeat,
            orthoMod <- data.frame(
                do.call('rbind', strsplit(orthoID,'|',fixed=TRUE))
            )
        )
        dataHeat$orthoIDNew <- paste(
            dataHeat$orthoMod$X2, dataHeat$orthoMod$X3, sep = "#"
        )
    } else {
        dataHeat$orthoIDNew <- paste(
            dataHeat$supertaxonID, dataHeat$orthoID, sep = "#"
        )
    }
    dataHeat <- dataHeat[ , !(names(dataHeat) %in% ("orthoMod"))]
    # Count occurrences of orthoIDNew and label as "Single" or "Multiple"
    countOrthoDf <- dataHeat %>%
        count(orthoIDNew, name = "orthoFreq") %>%
        mutate(orthoFreq = ifelse(orthoFreq > 1, "Multiple", "Single"))
    # Merge the frequency labels back to the main data
    dataHeat <- dataHeat %>%
        left_join(countOrthoDf, by = "orthoIDNew")
    # Check for "Multiple" freq in co-orthologs at geneID and supertaxonID level
    freqDt <- dataHeat %>%
        distinct(geneID, supertaxonID, orthoFreq) %>%
        group_by(geneID, supertaxonID) %>%
        mutate(
            orthoFreq = ifelse("Multiple" %in% orthoFreq, "Multiple", orthoFreq)
        ) %>% ungroup()
    # Merge the updated frequency labels back to the main data
    dataHeat <- dataHeat %>% left_join(
        freqDt, by = c("geneID", "supertaxonID"), suffix = c("", "New"), 
        relationship = "many-to-many"
    )
    # Step 3: Update orthoFreq with the new frequency labels and clean up
    dataHeat <- dataHeat %>% mutate(orthoFreq = orthoFreqNew) %>%
        select(-orthoFreqNew)
    return(dataHeat)
}

#' Create data for main profile plot
#' @export
#' @param dataHeat a data frame contains processed profiles (see
#' ?fullProcessedProfile, ?filterProfileData)
#' @return A dataframe for plotting the phylogenetic profile, containing seed
#' protein IDs (geneID), ortholog IDs (orthoID) together with their ncbi
#' taxonomy IDs (ncbiID and abbrName), full names (fullName), indexed supertaxa
#' (supertaxon), values for additional variables (var1, var2) and the aggregated
#' values of those additional variables for each supertaxon (mVar1, mVar2),
#' number of original and filtered co-orthologs in each supertaxon (paralog and
#' paralogNew), number of species in each supertaxon (numberSpec) and the % of
#' species that have orthologs in each supertaxon (presSpec).
#' @importFrom stats na.omit
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{filterProfileData}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' dataMainPlot(finalProcessedProfile)

dataMainPlot <- function(dataHeat = NULL){
    if (is.null(dataHeat)) stop("Input data cannot be NULL!")
    paralogNew <- orthoID <- NULL
    # reduce number of inparalogs based on filtered dataHeat
    dataHeatTb <- data.table(stats::na.omit(dataHeat))
    dataHeatTb[, paralogNew := .N, by = c("geneID", "supertaxon")]
    dataHeatTb <- data.frame(dataHeatTb[, c(
        "geneID", "supertaxon", "paralogNew"
    )])

    dataHeat <- dataHeat %>% left_join(dataHeatTb, by=c("geneID","supertaxon"))
    dataHeat$paralog <- dataHeat$paralogNew
    dataHeat <- dataHeat[!duplicated(dataHeat), ]

    # remove unneeded dots
    dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
    dataHeat$paralog[dataHeat$presSpec < 1] <- NA
    dataHeat$paralog[dataHeat$paralog == 1] <- NA
    # rescale numbers of paralogs
    if (length(unique(stats::na.omit(dataHeat$paralog))) > 0) {
        maxParalog <- max(stats::na.omit(dataHeat$paralog))
        dataHeat$paralogSize <- (dataHeat$paralog / maxParalog) * 3
    }
    # remove prefix number of taxa names but keep the order
    dataHeat$supertaxon <- factor(
        substr(
            as.character(dataHeat$supertaxon), 8 ,
            nchar(as.character(dataHeat$supertaxon))),
        levels = substr(
            levels(as.factor(dataHeat$supertaxon)), 8,
            nchar(levels(as.factor(dataHeat$supertaxon)))))
    # count ortholog IDs for each taxon
    dataHeat <- processOrthoID(dataHeat)
    return(dataHeat)
}

#' Create data for customized profile plot
#' @description Create data for customized profile plot based on a selected
#' list of genes and/or taxa, containing seed protein IDs (geneID), ortholog IDs
#' (orthoID) together with their ncbi taxonomy IDs (ncbiID and abbrName), full
#' names (fullName), indexed supertaxa (supertaxon), values for additional
#' variables (var1, var2) and the aggregated values of those additional
#' variables for each supertaxon (mVar1, mVar2), number of original and filtered
#' co-orthologs in each supertaxon (paralog and paralogNew), number of species
#' in each supertaxon (numberSpec) and the % of species that have orthologs in
#' each supertaxon (presSpec).
#' @export
#' @usage dataCustomizedPlot(dataHeat = NULL, selectedTaxa = "all",
#'     selectedSeq = "all")
#' @param dataHeat a data frame contains processed profiles (see
#' ?fullProcessedProfile, ?filterProfileData)
#' @param selectedTaxa selected subset of taxa. Default = "all".
#' @param selectedSeq selected subset of genes. Default = "all".
#' @return A dataframe contains data for plotting the customized profile.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @importFrom stats na.omit
#' @seealso \code{\link{filterProfileData}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' selectedTaxa <- c("Mammalia", "Saccharomycetes", "Insecta")
#' selectedSeq <- "all"
#' dataCustomizedPlot(finalProcessedProfile, selectedTaxa, selectedSeq)

dataCustomizedPlot <- function(
    dataHeat = NULL, selectedTaxa = "all", selectedSeq = "all"
){
    if (is.null(dataHeat)) stop("Input data cannot be NULL!")
    geneID <- supertaxonMod <- paralogNew <- NULL
    dataHeat$supertaxonMod <- {
        substr(
            dataHeat$supertaxon, 8, nchar(as.character(dataHeat$supertaxon))
        )
    }
    if (selectedTaxa[1] == "all" & selectedSeq[1] != "all") {
        dataHeat <- subset(dataHeat, geneID %in% selectedSeq)
    } else if (selectedSeq[1] == "all" & selectedTaxa[1] != "all") {
        dataHeat <- subset(dataHeat, supertaxonMod %in% selectedTaxa)
    } else {
        dataHeat <- subset(
            dataHeat, geneID %in% selectedSeq & supertaxonMod %in% selectedTaxa
        )
    }
    # reduce number of inparalogs based on filtered dataHeat
    dataHeatTb <- data.table(stats::na.omit(dataHeat))
    dataHeatTb[, paralogNew := .N, by = c("geneID", "supertaxon")]
    dataHeatTb <- data.frame(
        dataHeatTb[, c("geneID", "supertaxon", "paralogNew")]
    )
    dataHeat <- merge(
        dataHeat, dataHeatTb, by = c("geneID", "supertaxon"), all.x = TRUE
    )
    dataHeat$paralog <- dataHeat$paralogNew
    dataHeat <- dataHeat[!duplicated(dataHeat), ]
    # remove unneeded dots
    dataHeat$presSpec[dataHeat$presSpec == 0] <- NA
    dataHeat$paralog[dataHeat$presSpec < 1] <- NA
    dataHeat$paralog[dataHeat$paralog == 1] <- NA
    # rescale numbers of paralogs
    if (length(unique(stats::na.omit(dataHeat$paralog))) > 0) {
        maxParalog <- max(stats::na.omit(dataHeat$paralog))
        dataHeat$paralogSize <- (dataHeat$paralog / maxParalog) * 3
    }
    # remove prefix number of taxa names but keep the order
    dataHeat$supertaxon <- factor(
        substr(
            as.character(dataHeat$supertaxon), 8,
            nchar(as.character(dataHeat$supertaxon))),
        levels = substr(
            levels(dataHeat$supertaxon), 8,
            nchar(levels(dataHeat$supertaxon))))
    # count ortholog IDs for each taxon
    if (nrow(dataHeat) == 0) return(dataHeat)
    dataHeat <- processOrthoID(dataHeat)
    return(dataHeat)
}

#' Create profile heatmap plot
#' @export
#' @param data dataframe for plotting the heatmap phylogentic profile (either
#' full or subset profiles)
#' @param parm plot parameters, including (1) type of x-axis "taxa" or
#' "genes" - default = "taxa"; (2) display gene IDs (default) or gene names;
#' (3+4) names of 2 variables var1ID and var2ID - default = "var1" & "var2"; 
#' (5+6) mid value and color for mid value of var1 -
#' default is 0.5 and #FFFFFF; (7) color for lowest var1 - default = "#FF8C00";
#' (8) color for highest var1 - default = "#4682B4"; (9+10) mid value and color
#' for mid value of var2 - default is 1 and #FFFFFF;(11) color for lowest var2 -
#' default = "#FFFFFF", (12) color for highest var2 - default = "#F0E68C", (13)
#' color of co-orthologs - default = "#07D000"; (14+15+16) text sizes for x, y
#' axis and legend - default = 9 for each; (17) legend position "top", "bottom",
#' "right", "left" or "none" - default = "top"; (18) zoom ratio of the
#' co-ortholog dots from -1 to 3 - default = 0; (19) angle of x-axis from 0 to
#' 90 - default = 60; (20) show/hide separate line for reference taxon 1/0 -
#' default = 0; (21) enable/disable coloring gene categories TRUE/FALSE -
#' default = FALSE; (22) enable/disable coloring duplicated ortholog IDs
#' TRUE/FALSE - default=FALSE). NOTE: Leave blank or NULL to use default values.
#' @return A profile heatmap plot as a ggplot object.
#' @import ggplot2
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{dataMainPlot}}, \code{\link{dataCustomizedPlot}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(finalProcessedProfile)
#' plotParameter <- list(
#'     "xAxis" = "taxa",
#'     "geneIdType" = "geneID",
#'     "var1ID" = "FAS_FW",
#'     "var2ID"  = "FAS_BW",
#'     "midVar1" = 0.5,
#'     "midColorVar1" =  "#FFFFFF",
#'     "lowColorVar1" =  "#FF8C00",
#'     "highColorVar1" = "#4682B4",
#'     "midVar2" = 1,
#'     "midColorVar2" =  "#FFFFFF",
#'     "lowColorVar2" = "#CB4C4E",
#'     "highColorVar2" = "#3E436F",
#'     "paraColor" = "#07D000",
#'     "xSize" = 8,
#'     "ySize" = 8,
#'     "legendSize" = 8,
#'     "mainLegend" = "top",
#'     "dotZoom" = 0,
#'     "xAngle" = 60,
#'     "guideline" = 0,
#'     "colorByGroup" = FALSE,
#'     "catColors" = NULL,
#'     "colorByOrthoID" = FALSE
#' )
#'
#' heatmapPlotting(plotDf, plotParameter)

heatmapPlotting <- function(data = NULL, parm = NULL){
    if (is.null(data)) stop("Input data cannot be NULL!")
    if (is.null(parm))
        parm <- list(
            "xAxis" = "taxa", "geneIdType" = "geneName", 
            "var1ID" = "var1", "var2ID"  = "var2",
            "midVar1" = 0.5, "midColorVar1" =  "#FFFFFF",
            "lowColorVar1" =  "#FF8C00", "highColorVar1" = "#4682B4",
            "midVar2" = 1, "midColorVar2" =  "#FFFFFF",
            "lowColorVar2" = "#CB4C4E", "highColorVar2" = "#3E436F",
            "paraColor" = "#07D000", "xSize" = 8, "ySize" = 8, "legendSize" = 8,
            "mainLegend" = "top", "dotZoom" = 0, "xAngle" = 60, "guideline" = 0,
            "colorByGroup" = FALSE,"catColors" = NULL,"colorByOrthoID" = FALSE)
    geneID<-geneName<-supertaxon<-category<-var1<-var2<- presSpec<-paralog<-NULL
    orthoFreq <- xmin <- xmax <- ymin <- ymax <- NULL
    ### create heatmap plot
    if (!(is.null(parm$geneIdType))) {
        if (parm$geneIdType == "geneName") data$geneID <- data$geneName
    }
    # create a canvas
    if (parm$xAxis == "genes") {
        p <- ggplot(data,aes(x=geneID, y=supertaxon)) +
            labs(x = "Gene ID", y = "Taxon")
    } else {
        p <- ggplot(data, aes(y = geneID, x = supertaxon)) +
            labs(y = "Gene ID", x = "Taxon")
    }
    # create geom_tile & scale_fill_gradient for var2 OR gene category
    if (parm$colorByGroup == TRUE) {
        p <- p + geom_raster(aes(fill = factor(category)), alpha = 0.3)
        if (!is.null(parm$catColors))
            p <- p + scale_fill_manual(values = parm$catColors)
    } else {
        if (length(unique(stats::na.omit(data$var2))) != 1)
            p <- p + scale_fill_gradient2(
                low = parm$lowColorVar2, high = parm$highColorVar2,
                mid = parm$midColorVar2, midpoint = parm$midVar2,
                na.value = "gray95", limits = c(0, 1)) +
                geom_raster(aes(fill = var2))
    }

    # create geom_point for found ortho; coloring by var1 or
    # orthoIDs (only when working on lowest taxonomy rank)
    # and inparalogs (co-orthologs)
    if (length(unique(stats::na.omit(data$presSpec))) < 2) {
        # working on the lowest taxonomy rank
        if (parm$colorByOrthoID == TRUE) {
            # color by ortho IDs
            p <- p + geom_point(
                aes(colour = factor(orthoFreq)), na.rm = TRUE,
                size = data$presSpec*5*(1+parm$dotZoom), show.legend = TRUE
            )
        } else {
            if (length(unique(stats::na.omit(data$var1))) == 1) {
                p <- p + geom_point(
                    aes(colour = var1), na.rm = TRUE,
                    size = data$presSpec*5*(1+parm$dotZoom), show.legend = FALSE
                )
            } else {
                p <- p +
                    geom_point(
                        aes(colour = var1), na.rm = TRUE,
                        size = data$presSpec * 5 * (1 + parm$dotZoom)) +
                    scale_color_gradient2(
                        low = parm$lowColorVar1, high = parm$highColorVar1,
                        mid = parm$midColorVar1, midpoint = parm$midVar1,
                        limits = c(0,1)
                    )
            }
            # plot inparalogs (if available)
            if (length(unique(stats::na.omit(data$paralog))) > 0) {
                # Define min and max for your desired size range
                min_size <- 1
                max_size <- 3*(1+parm$dotZoom)
                # Manually rescale paralog size
                min_val <- min(data$paralog, na.rm = TRUE)
                max_val <- max(data$paralog, na.rm = TRUE)
                # Avoid division by zero (when min == max)
                if (min_val == max_val) {
                    data$paralogSize <- max_size
                } else {
                    data$paralogSize <- 
                        ((data$paralog - min_val) / (max_val - min_val)) * 
                        (max_size - min_size) + min_size
                }
                # create breaks for legend
                paralog_vals <- sort(unique(na.omit(data$paralog)))
                n_vals <- length(paralog_vals)
                custom_breaks <- if (n_vals <= 5) {
                    paralog_vals
                } else {
                    paralog_vals[round(seq(2, n_vals, length.out = 5))]
                }
                p <- p +
                    geom_point(
                        data = data, aes(size = paralog), color=parm$paraColor,
                        na.rm = TRUE, show.legend = TRUE) +
                    guides(size = guide_legend(title = "# of co-orthologs")) +
                scale_size_continuous(
                    breaks = custom_breaks, range = c(min_size, max_size))
            }
        }
    } else {
        # scale dot size based on % preSpec in each super taxon
        if (length(unique(stats::na.omit(data$var1))) == 1) {
            p <- p + geom_point(aes(size=presSpec),color="#336a98",na.rm = TRUE)
        } else {
            p <- p + geom_point(aes(colour=var1, size = presSpec),na.rm = TRUE)+
                scale_color_gradient2(
                    low = parm$lowColorVar1, high = parm$highColorVar1,
                    mid = parm$midColorVar1, midpoint = parm$midVar1,
                    limits = c(0,1)
                )
        }
        # remain the scale of point while filtering % presSpec
        presentVl <- data$presSpec[!is.na(data$presSpec)]
        p <- p + scale_size_continuous(range = c(
            (floor(min(presentVl) * 10) / 10 * 5) * (1 + parm$dotZoom),
            (floor(max(presentVl) * 10) / 10 * 5) * (1 + parm$dotZoom)))
    }

    # create legend
    if (parm$colorByGroup == FALSE) {
        if (parm$colorByOrthoID == FALSE) {
            p <- p + guides(fill = guide_colourbar(title = parm$var2ID),
                            color = guide_colourbar(title = parm$var1ID))
        } else {
            p <- p + guides(fill = guide_colourbar(title = parm$var2ID),
                            color = guide_legend(title = 'OrthoID copy number'))
        }
    } else {
        p <- p + guides(fill = guide_legend("Category"),
                        color = guide_colourbar(title = parm$var1ID))
    }

    # text size, legend position
    p <- p + theme_minimal(base_size = 9)
    vjustValue <- 1
    if (parm$xAngle == 90) vjustValue <- 0.5
    p <- p + theme(
        axis.text.x = element_text(
            angle = parm$xAngle, hjust = 1, size = parm$xSize,
            vjust = vjustValue
        ),
        axis.text.y = element_text(size = parm$ySize),
        axis.title.x = element_text(size = parm$xSize),
        axis.title.y = element_text(size = parm$ySize),
        legend.title = element_text(size = parm$legendSize),
        legend.text = element_text(size = parm$legendSize),
        legend.position = parm$mainLegend)
    return(p)
}


#' Create profile heatmap plot using scattermore
#' @export
#' @param data dataframe for plotting the heatmap phylogentic profile (either
#' full or subset profiles)
#' @param parm plot parameters, including (1) type of x-axis "taxa" or
#' "genes" - default = "taxa"; (2) display gene IDs (default) or gene names;
#' (3+4) names of 2 variables var1ID and var2ID - default = "var1" & "var2"; 
#' (5+6) mid value and color for mid value of var1 -
#' default is 0.5 and #FFFFFF; (7) color for lowest var1 - default = "#FF8C00";
#' (8) color for highest var1 - default = "#4682B4"; (9+10) mid value and color
#' for mid value of var2 - default is 1 and #FFFFFF;(11) color for lowest var2 -
#' default = "#FFFFFF", (12) color for highest var2 - default = "#F0E68C", (13)
#' color of co-orthologs - default = "#07D000"; (14+15+16) text sizes for x, y
#' axis and legend - default = 9 for each; (17) legend position "top", "bottom",
#' "right", "left" or "none" - default = "top"; (18) zoom ratio of the
#' co-ortholog dots from -1 to 3 - default = 0; (19) color dots based on either 
#' "var1" or "var2". NOTE: Leave blank or NULL to 
#' use default values.
#' @return A profile heatmap plot as a ggplot object.
#' @import ggplot2
#' @import scattermore
#' @importFrom grDevices colorRampPalette
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{dataMainPlot}}, \code{\link{dataCustomizedPlot}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(finalProcessedProfile)
#' plotParameter <- list(
#'     "xAxis" = "taxa",
#'     "geneIdType" = "geneID",
#'     "var1ID" = "FAS_FW",
#'     "var2ID"  = "FAS_BW",
#'     "midVar1" = 0.5,
#'     "midColorVar1" =  "#FFFFFF",
#'     "lowColorVar1" =  "#FF8C00",
#'     "highColorVar1" = "#4682B4",
#'     "midVar2" = 1,
#'     "midColorVar2" =  "#FFFFFF",
#'     "lowColorVar2" = "#CB4C4E",
#'     "highColorVar2" = "#3E436F",
#'     "paraColor" = "#07D000",
#'     "xSize" = 8,
#'     "ySize" = 8,
#'     "legendSize" = 8,
#'     "mainLegend" = "top",
#'     "dotZoom" = 0,
#'     "colorVar" = "var1"
#' )
#'
#' heatmapPlottingFast(plotDf, plotParameter)

heatmapPlottingFast <- function(data = NULL, parm = NULL) {
    if (is.null(data)) stop("Input data cannot be NULL!")
    if (is.null(parm))
        parm <- list(
            "xAxis" = "taxa", "geneIdType" = "geneName", 
            "var1ID" = "var1", "var2ID"  = "var2",
            "midVar1" = 0.5, "midColorVar1" =  "#FFFFFF",
            "lowColorVar1" =  "#FF8C00", "highColorVar1" = "#4682B4",
            "midVar2" = 1, "midColorVar2" =  "#FFFFFF",
            "lowColorVar2" = "#CB4C4E", "highColorVar2" = "#3E436F",
            "paraColor" = "#07D000", "xSize" = 8, "ySize" = 8, "legendSize" = 8,
            "mainLegend" = "top", "dotZoom" = 0, "colorVar" = "var1")
    geneID<-geneName<-supertaxon<-category<-var1<-var2<- presSpec<-paralog<-NULL
    orthoFreq <- xmin <- xmax <- ymin <- ymax <- x <- y <- NULL
    if (is.null(parm$colorVar)) parm$colorVar <- "var1"
    ### create heatmap plot
    if (!(is.null(parm$geneIdType))) {
        if (parm$geneIdType == "geneName") data$geneID <- data$geneName
    }   
    # create a canvas
    if (parm$xAxis == "genes") {
        p <- ggplot(data,aes(x = geneID, y = supertaxon)) +
            labs(x = "Gene ID", y = "Taxon")
    } else {
        p <- ggplot(data, aes(y = geneID, x = supertaxon)) +
            labs(y = "Gene ID", x = "Taxon")
    }
    # format theme
    p <- p + theme_minimal(base_size = 9)
    p <- p + theme(
        axis.text = element_blank(),
        axis.title.x = element_text(size = parm$xSize),
        axis.title.y = element_text(size = parm$ySize),
        legend.title = element_text(size = parm$legendSize),
        legend.text = element_text(size = parm$legendSize),
        legend.position = parm$mainLegend)
    # Extract x, y, and color values
    x_values <- as.numeric(data$supertaxon)
    y_values <- as.numeric(data$geneID)
    # Create color vector: green if paralog > 1, else gradient by var1 or var2
    if (parm$colorVar == "var1") {
        cl <- colorRampPalette(c(parm$lowColorVar1, parm$highColorVar1))(100)[
            cut(data$var1, breaks = 100)
        ]
    } else
        cl <- colorRampPalette(c(parm$lowColorVar2, parm$highColorVar2))(100)[
            cut(data$var2, breaks = 100)
        ]
    cl[data$paralog > 1] <- parm$paraColor
    # Add scatter plot
    p <- p + 
        geom_scattermost(
            cbind(x_values, y_values), color = cl,
            pointsize = 2 * (1 + parm$dotZoom), pixels = c(1000, 1000)
        )
    if (parm$colorVar == "var1") {
        p <- p +
            geom_point(data=data.frame(x = double(0)), aes(x, x, color = x)) +
            scale_color_gradientn(
                colors = colorRampPalette(
                    c(parm$lowColorVar1, parm$highColorVar1))(100),
                limits = c(0,1), name = parm$var1ID
            )
    } else {
        p <- p +
            geom_point(data=data.frame(x = double(0)), aes(x, x, color = x)) +
            scale_color_gradientn(
                colors = colorRampPalette(
                    c(parm$lowColorVar2, parm$highColorVar2))(100),
                limits = c(0,1), name = parm$var2ID
            )
    }
    # Add only one real co-ortholog point with fill = green to generate legend
    if(any(data$paralog > 1, na.rm = TRUE)) {
        co_point <- data[data$paralog > 1, ][1, ]
        co_point$x <- as.numeric(co_point$supertaxon)
        co_point$y <- as.numeric(co_point$geneID)
        p <- p + geom_point(
            data = co_point, aes(x = x, y = y, fill = "Co-orthologs"),
            shape = 21, size = 3
        ) + scale_fill_manual(
            values = c("Co-orthologs" = parm$paraColor), name = NULL
        )
    }
    return(p)
}

#' Highlight gene and/or taxon of interest on the phylogenetic profile plot
#' @export
#' @usage highlightProfilePlot(profilePlot = NULL, plotDf = NULL,
#'     taxonHighlight = "none", workingRank = "none", geneHighlight = NULL,
#'     taxDB = NULL, xAxis = "taxa")
#' @param profilePlot initial (highlighted) profile plot
#' @param plotDf dataframe for plotting the heatmap phylogentic profile
#' @param taxonHighlight taxon of interst. Default = "none".
#' @param workingRank working taxonomy rank (needed only for highlight taxon).
#' @param geneHighlight gene of interest. Default = NULL.
#' @param taxDB Path to the taxonomy DB files
#' @param xAxis type of x-axis (either "genes" or "taxa")
#' @return A profile heatmap plot with highlighted gene and/or taxon of interest
#' as ggplot object.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{dataMainPlot}}, \code{\link{dataCustomizedPlot}},
#' \code{\link{heatmapPlotting}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(finalProcessedProfile)
#' plotParameter <- list(
#'     "xAxis" = "taxa",
#'     "geneIdType" = "geneID",
#'     "var1ID" = "FAS_FW",
#'     "var2ID"  = "FAS_BW",
#'     "midVar1" = 0.5,
#'     "midColorVar1" =  "#FFFFFF",
#'     "lowColorVar1" =  "#FF8C00",
#'     "highColorVar1" = "#4682B4",
#'     "midVar2" = 1,
#'     "midColorVar2" =  "#FFFFFF",
#'     "lowColorVar2" = "#CB4C4E",
#'     "highColorVar2" = "#3E436F",
#'     "paraColor" = "#07D000",
#'     "xSize" = 8,
#'     "ySize" = 8,
#'     "legendSize" = 8,
#'     "mainLegend" = "top",
#'     "dotZoom" = 0,
#'     "xAngle" = 60,
#'     "guideline" = 0,
#'     "colorByGroup" = FALSE,
#'     "colorByOrthoID" = FALSE
#' )
#' profilePlot <- heatmapPlotting(plotDf, plotParameter)
#' taxonHighlight <- "none"
#' workingRank <- "class"
#' geneHighlight <- "100265at6656"
#' highlightProfilePlot(
#'     profilePlot, plotDf, taxonHighlight, workingRank, geneHighlight,
#'     NULL, plotParameter$xAxis
#' )

highlightProfilePlot <- function(
        profilePlot = NULL, plotDf = NULL, taxonHighlight = "none",
        workingRank = "none", geneHighlight = NULL, taxDB = NULL, xAxis = "taxa"
){
    if (is.null(plotDf)) stop("Input data cannot be NULL!")
    if (is.null(profilePlot)) stop("Profile plot cannot be NULL!")
    xmin <- xmax <- ymin <- ymax <- NULL
    if (is.null(taxonHighlight) && is.null(geneHighlight)) return(profilePlot)
    # highlight taxon
    if (length(taxonHighlight) > 0 && !("none" %in% taxonHighlight)) {
        # get selected highlight taxon ID
        taxName <- getNameList(taxDB)
        taxonHighlightID <- taxName$ncbiID[
            taxName$fullName %in% taxonHighlight & taxName$rank == workingRank]
        if (length(taxonHighlightID) == 0L) {
            taxonHighlightID <- taxName$ncbiID[
                taxName$fullName %in% taxonHighlight
            ]
        } else if(length(taxonHighlightID) < length(taxonHighlight)) {
            missingID <- taxName$ncbiID[
                taxName$fullName %in% taxonHighlight
            ]
            missingID <- missingID[!(missingID %in% taxonHighlightID)]
            taxonHighlightID <- c(taxonHighlightID, missingID)
        }
        # get taxonID together with it sorted index
        taxonHighlightID <- intersect(taxonHighlightID, plotDf$supertaxonID)
        if (length(taxonHighlightID) > 0) {
            selTaxon <- unique(as.character(
                plotDf$supertaxon[plotDf$supertaxonID %in% taxonHighlightID]
            ))
            selIndex <- match(selTaxon, levels(as.factor(plotDf$supertaxon)))
            if (xAxis == "taxa") {
                rect <- data.frame(
                    xmin=selIndex-0.5, xmax=selIndex+0.5, ymin=-Inf, ymax=Inf)
            } else
                rect <- data.frame(
                    ymin=selIndex-0.5, ymax=selIndex+0.5, xmin=-Inf, xmax=Inf)
            profilePlot <- profilePlot + geom_rect(
                data = rect, color = "yellow", alpha = 0.3, inherit.aes = FALSE,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
        }
    }
    # highlight gene
    if (length(geneHighlight) > 0 && !("none" %in% geneHighlight)) {
        geneHighlight <- intersect(geneHighlight, unique(plotDf$geneID))
        if (length(geneHighlight) > 0) {
            selIndex <- match(geneHighlight, levels(as.factor(plotDf$geneID)))
            if (xAxis == "taxa") {
                rect <- data.frame(
                    ymin=selIndex-0.5, ymax=selIndex+0.5, xmin=-Inf, xmax=Inf)
            } else
                rect <- data.frame(
                    xmin=selIndex-0.5, xmax=selIndex+0.5, ymin=-Inf, ymax=Inf)
            profilePlot <- profilePlot + geom_rect(
                data = rect, color = "yellow", alpha = 0.3, inherit.aes = FALSE,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
        }
    }
    return(profilePlot)
}

#' Add taxonomy rank division lines to the heatmap plot
#' @export
#' @usage addRankDivisionPlot(profilePlot = NULL, plotDf = NULL,
#'     taxDB = NULL, workingRank = NULL, superRank = NULL, xAxis = "taxa",
#'     font = "Arial", groupLabelSize = 14, groupLabelDist = 2,
#'     groupLabelAngle = 90, refLine = TRUE)
#' @param profilePlot initial (highlighted) profile plot
#' @param plotDf dataframe for plotting the heatmap phylogentic profile
#' @param taxDB path to taxonomy database (taxonomyMatrix.txt file required!)
#' @param workingRank working taxonomy rank (e.g. species)
#' @param superRank taxonomy rank for division lines (e.g. superkingdom)
#' @param xAxis type of x-axis (either "genes" or "taxa")
#' @param font font of text. Default = Arial"
#' @param groupLabelSize size of rank labels
#' @param groupLabelDist size of the plot area for rank labels
#' @param groupLabelAngle angle of rank labels
#' @param refLine add vertical line to separate reference taxon
#' @return A profile heatmap plot with highlighted gene and/or taxon of interest
#' as ggplot object.
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{heatmapPlotting}}, \code{\link{highlightProfilePlot}},
#' \code{\link{getTaxonomyMatrix}}
#' @examples
#' data("finalProcessedProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(finalProcessedProfile)
#' plotParameter <- list(
#'     "xAxis" = "taxa",
#'     "geneIdType" = "geneID",
#'     "var1ID" = "FAS_FW",
#'     "var2ID"  = "FAS_BW",
#'     "midVar1" = 0.5,
#'     "midColorVar1" =  "#FFFFFF",
#'     "lowColorVar1" =  "#FF8C00",
#'     "highColorVar1" = "#4682B4",
#'     "midVar2" = 1,
#'     "midColorVar2" =  "#FFFFFF",
#'     "lowColorVar2" = "#CB4C4E",
#'     "highColorVar2" = "#3E436F",
#'     "paraColor" = "#07D000",
#'     "xSize" = 8,
#'     "ySize" = 8,
#'     "legendSize" = 8,
#'     "mainLegend" = "top",
#'     "dotZoom" = 0,
#'     "xAngle" = 60,
#'     "guideline" = 0,
#'     "colorByGroup" = FALSE,
#'     "colorByOrthoID" = FALSE
#' )
#' profilePlot <- heatmapPlotting(plotDf, plotParameter)
#' workingRank <- "class"
#' superRank <- "superkingdom"
#' addRankDivisionPlot(
#'     profilePlot, plotDf, NULL, workingRank, superRank, "taxa", font = "sans"
#' )

addRankDivisionPlot <- function(
        profilePlot = NULL, plotDf = NULL, taxDB = NULL,
        workingRank = NULL, superRank = NULL, xAxis = "taxa", font = "Arial",
        groupLabelSize = 14, groupLabelDist = 2, groupLabelAngle = 90,
        refLine = TRUE
) {
    if (is.null(plotDf)) stop("Input data cannot be NULL!")
    if (is.null(profilePlot)) stop("Profile plot cannot be NULL!")
    profilePlot <- profilePlot + theme(text = element_text(family = font))
    if (is.null(workingRank) || is.null(superRank)) {
        if (refLine == TRUE) {
            # guideline for separating ref species
            if (xAxis == "genes") {
                profilePlot <- profilePlot + labs(y = "Taxon") +
                    geom_hline(yintercept = 0.5, colour = "dodgerblue4") +
                    geom_hline(yintercept = 1.5, colour = "dodgerblue4")
            } else
                profilePlot <- profilePlot + labs(x = "Taxon") +
                    geom_vline(xintercept = 0.5, colour = "dodgerblue4") +
                    geom_vline(xintercept = 1.5, colour = "dodgerblue4")
        }
        return(profilePlot)
    } else {
        # sort supertaxonID based on the sorted supertaxon
        plotDf <- plotDf[with(plotDf,order(supertaxon)),]
        plotDf$supertaxonID <- factor(
            plotDf$supertaxonID, levels = unique(plotDf$supertaxonID)
        )
        # get input (super)taxa
        inputTax <- levels(as.factor(plotDf$supertaxonID))
        # subset taxonomy matrix to contain only superrank for division line
        taxMatrix <- getTaxonomyMatrix(taxDB)
        subTaxMatrix <- subset(
            taxMatrix[taxMatrix[[workingRank]] %in% inputTax,],
            select = c(as.character(workingRank), as.character(superRank))
        )
        subTaxMatrix <- subTaxMatrix[!duplicated(subTaxMatrix),]
        # remove IDs that do not have real NCBI superRank
        nameList <- getNameList(taxDB)
        subNameList <- nameList[
            nameList$ncbiID %in% subTaxMatrix[[superRank]],
            c("ncbiID", "rank", "fullName")
        ]
        colnames(subNameList) <- c(superRank, "rank", "name")
        mergedDf <- merge(subTaxMatrix, subNameList, by = superRank, all.x=TRUE)
        mergedDf <- mergedDf[mergedDf$rank == superRank,]
        subTaxMatrix <- mergedDf
        # group input taxa based on superrank and get their index
        groupedList <- lapply(
            levels(as.factor(subTaxMatrix[[superRank]])),
            function(x) {
                tmp <- subTaxMatrix[subTaxMatrix[[superRank]] == x,]
                tmp$index <- which(
                    levels(as.factor(plotDf$supertaxonID))
                    %in% tmp[[workingRank]]
                )
                return(tmp)
            }
        )
        # add vertical line to divide taxon groups
        max_taxa <- length(unique(as.character(plotDf$geneID)))
        for(i in groupedList) {
            min <- min(i$index)
            max <- max(i$index)
            if (xAxis == "taxa") {
                profilePlot <- profilePlot +
                    geom_vline(xintercept = min - 0.5, colour = "dodgerblue4") +
                    geom_vline(xintercept = max + 0.5, colour = "dodgerblue4") +
                    annotate(
                        geom = "text", angle = groupLabelAngle, hjust = 0,
                        size = groupLabelSize,
                        x = min,
                        y = max_taxa + 0.5,
                        label = unique(as.character(i$name))
                    )

            } else {
                profilePlot <- profilePlot +
                    geom_hline(yintercept = min - 0.5, colour = "dodgerblue4") +
                    geom_hline(yintercept = max + 0.5, colour = "dodgerblue4") +
                    annotate(
                        geom = "text", angle = groupLabelAngle - 90, hjust = 0,
                        size = groupLabelSize,
                        y = min,
                        x = max_taxa + 1,
                        label = unique(as.character(i$name))
                    )
            }
        }
        if (xAxis == "taxa") {
            return(
                profilePlot + coord_cartesian(
                    clip = 'off', ylim = c(1, max_taxa + groupLabelDist)
                )
            )
        } else
            return(
                profilePlot + coord_cartesian(
                    clip = 'off', xlim = c(1, max_taxa + groupLabelDist)
                )
            )
    }
}
