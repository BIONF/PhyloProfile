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
#' @import data.table
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{filterProfileData}}
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' dataMainPlot(fullProcessedProfile)

dataMainPlot <- function(dataHeat = NULL){
    if (is.null(dataHeat)) stop("Input data cannot be NULL!")
    paralogNew <- NULL
    
    # reduce number of inparalogs based on filtered dataHeat
    dataHeatTb <- data.table(na.omit(dataHeat))
    dataHeatTb[, paralogNew := .N, by = c("geneID", "supertaxon")]
    dataHeatTb <- data.frame(dataHeatTb[, c(
        "geneID", "supertaxon", "paralogNew"
    )])
    
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
    if (length(unique(na.omit(dataHeat$paralog))) > 0) {
        maxParalog <- max(na.omit(dataHeat$paralog))
        dataHeat$paralogSize <- (dataHeat$paralog / maxParalog) * 3
    }
    # remove prefix number of taxa names but keep the order
    dataHeat$supertaxon <- factor(
        substr(
            as.character(dataHeat$supertaxon), 6 ,
            nchar(as.character(dataHeat$supertaxon))),
        levels = substr(
            levels(dataHeat$supertaxon), 6, 
            nchar(levels(dataHeat$supertaxon))))
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
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{filterProfileData}}
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' selectedTaxa <- c("Mammalia", "Saccharomycetes", "Insecta")
#' selectedSeq <- "all"
#' dataCustomizedPlot(fullProcessedProfile, selectedTaxa, selectedSeq)

dataCustomizedPlot <- function(
    dataHeat = NULL, selectedTaxa = "all", selectedSeq = "all"
){
    if (is.null(dataHeat)) stop("Input data cannot be NULL!")
    geneID <- supertaxonMod <- paralogNew <- NULL
    dataHeat$supertaxonMod <- {
        substr(
            dataHeat$supertaxon, 6, nchar(as.character(dataHeat$supertaxon))
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
    dataHeatTb <- data.table(na.omit(dataHeat))
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
    if (length(unique(na.omit(dataHeat$paralog))) > 0) {
        maxParalog <- max(na.omit(dataHeat$paralog))
        dataHeat$paralogSize <- (dataHeat$paralog / maxParalog) * 3
    }
    # remove prefix number of taxa names but keep the order
    dataHeat$supertaxon <- factor(
        substr(
            as.character(dataHeat$supertaxon), 6 ,
            nchar(as.character(dataHeat$supertaxon))),
        levels = substr(
            levels(dataHeat$supertaxon), 6, 
            nchar(levels(dataHeat$supertaxon))))
    return(dataHeat)
}


#' Create profile heatmap plot
#' @export
#' @param data dataframe for plotting the heatmap phylogentic profile (either
#' full or subset profiles)
#' @param parm plot parameters, including (1) type of x-axis "taxa" or
#' "genes" - default = "taxa"; (2+3) names of 2 variables var1ID and var2ID -
#' default = "var1" & "var2"; (4) color for lowest var1 - default = "#FF8C00";
#' (5) color for highest var1 - default = "#4682B4"; (6) color for lowest var2 -
#' default = "#FFFFFF", (7) color for highest var2 - default = "#F0E68C", (8)
#' color of co-orthologs - default = "#07D000"; (9+10+11) text sizes for x, y
#' axis and legend - default = 9 for each; (12) legend position "top", "bottom",
#' "right", "left" or "none" - default = "top"; (13) zoom ratio of the
#' co-ortholog dots from -1 to 3 - default = 0; (14) angle of x-axis from 0 to
#' 90 - default = 60; (14) show/hide separate line for reference taxon 1/0 -
#' default = 0; (15) enable/disable coloring gene categories TRUE/FALSE -
#' default = FALSE). NOTE: Leave blank or NULL to use default values.
#' @return A profile heatmap plot as a ggplot object.
#' @import ggplot2
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{dataMainPlot}}, \code{\link{dataCustomizedPlot}}
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(fullProcessedProfile)
#' plotParameter <- list(
#'     "xAxis" = "taxa",
#'     "var1ID" = "FAS_FW",
#'     "var2ID"  = "FAS_BW",
#'     "lowColorVar1" =  "#FF8C00",
#'     "highColorVar1" = "#4682B4",
#'     "lowColorVar2" = "#FFFFFF",
#'     "highColorVar2" = "#F0E68C",
#'     "paraColor" = "#07D000",
#'     "xSize" = 8,
#'     "ySize" = 8,
#'     "legendSize" = 8,
#'     "mainLegend" = "top",
#'     "dotZoom" = 0,
#'     "xAngle" = 60,
#'     "guideline" = 0,
#'     "colorByGroup" = FALSE
#' )
#'
#' heatmapPlotting(plotDf, plotParameter)

heatmapPlotting <- function(data = NULL, parm = NULL){
    if (is.null(data)) stop("Input data cannot be NULL!")
    if (is.null(parm))
        parm <- list(
            "xAxis" = "taxa", "var1ID" = "var1", "var2ID"  = "var2",
            "lowColorVar1" =  "#FF8C00", "highColorVar1" = "#4682B4",
            "lowColorVar2" = "#FFFFFF", "highColorVar2" = "#F0E68C",
            "paraColor" = "#07D000", "xSize" = 8, "ySize" = 8, "legendSize" = 8,
            "mainLegend" = "top", "dotZoom" = 0, "xAngle" = 60, "guideline" = 0,
            "colorByGroup" = FALSE)
    geneID <- supertaxon <- group <- var1 <- var2 <- presSpec <- paralog <- NULL
    xmin <- xmax <- ymin <- ymax <- NULL
    # create heatmap plot with geom_point & scale_color_gradient for present 
    # ortho & var1, geom_tile & scale_fill_gradient for var2
    if (parm$xAxis == "genes") p <- ggplot(data,aes(x = geneID, y = supertaxon))
    else p <- ggplot(data, aes(y = geneID, x = supertaxon))
    if (parm$colorByGroup == TRUE) {
        p <- p + geom_tile(aes(fill = factor(group)), alpha = 0.3)
    } else {
        if (length(unique(na.omit(data$var2))) != 1)
            p <- p + scale_fill_gradient(
                low = parm$lowColorVar2, high = parm$highColorVar2, 
                na.value = "gray95", limits = c(0, 1)) + 
                geom_tile(aes(fill = var2))
    }
    if (length(unique(na.omit(data$presSpec))) < 3) {
        if (length(unique(na.omit(data$var1))) == 1) {
            p <- p + geom_point(aes(colour = var1), na.rm = TRUE, 
                size = data$presSpec*5*(1+parm$dotZoom), show.legend = FALSE)
        } else
            p <- p + 
                geom_point(
                    aes(colour = var1), na.rm = TRUE,
                    size = data$presSpec * 5 * (1 + parm$dotZoom)) +
                scale_color_gradient(
                    low=parm$lowColorVar1,high=parm$highColorVar1,limits=c(0,1))
    } else {
        if (length(unique(na.omit(data$var1))) == 1) {
            p <- p + geom_point(aes(size=presSpec),color="#336a98",na.rm = TRUE)
        } else
            p <- p + geom_point(aes(colour=var1, size = presSpec),na.rm = TRUE)+
                scale_color_gradient(
                    low=parm$lowColorVar1,high=parm$highColorVar1,limits=c(0,1))
    }
    # plot inparalogs (if available)
    if (length(unique(na.omit(data$paralog))) > 0) {
        p <- p + 
            geom_point(
                data = data, aes(size = paralog), color = parm$paraColor,
                na.rm = TRUE, show.legend = TRUE) + 
            guides(size = guide_legend(title = "# of co-orthologs")) +
            scale_size_continuous(range = c(
                min(na.omit(data$paralogSize)) * (1 + parm$dotZoom),
                max(na.omit(data$paralogSize)) * (1 + parm$dotZoom)))
    } else {
        # remain the scale of point while filtering
        presentVl <- data$presSpec[!is.na(data$presSpec)]
        p <- p + scale_size_continuous(range = c(
            (floor(min(presentVl) * 10) / 10 * 5) * (1 + parm$dotZoom),
            (floor(max(presentVl) * 10) / 10 * 5) * (1 + parm$dotZoom)))
    }
    # color gene categories
    if (parm$colorByGroup == FALSE) {
        p <- p + guides(fill = guide_colourbar(title = parm$var2ID),
                        color = guide_colourbar(title = parm$var1ID))
    } else
        p <- p + guides(fill = guide_legend("Category"),
                        color = guide_colourbar(title = parm$var1ID))
    # guideline for separating ref species
    if (parm$guideline == 1) {
        if (parm$xAxis == "genes") {
            p <- p + labs(y = "Taxon") + 
                geom_hline(yintercept = 0.5, colour = "dodgerblue4") + 
                geom_hline(yintercept = 1.5, colour = "dodgerblue4")
        } else
            p <- p + labs(x = "Taxon") + 
                geom_vline(xintercept = 0.5, colour = "dodgerblue4") + 
                geom_vline(xintercept = 1.5, colour = "dodgerblue4")
    }
    p <- p + theme_minimal(base_size = 9)
    p <- p + theme(
        axis.text.x = element_text(angle=parm$xAngle, hjust=1, size=parm$xSize),
        axis.text.y = element_text(size = parm$ySize),
        axis.title.x = element_text(size = parm$xSize),
        axis.title.y = element_text(size = parm$ySize),
        legend.title = element_text(size = parm$legendSize),
        legend.text = element_text(size = parm$legendSize),
        legend.position = parm$mainLegend)
    return(p)
}

#' Highlight gene and/or taxon of interest on the phylogenetic profile plot
#' @export
#' @usage highlightProfilePlot(data, plotParameter = NULL, taxonHighlight =
#'     "none", rankName = "none", geneHighlight = "none")
#' @param data dataframe for plotting the heatmap phylogentic profile (either
#' full or subset profiles)
#' @param plotParameter plot parameters, including (1) type of x-axis "taxa" or
#' "genes" - default = "taxa"; (2+3) names of 2 variables var1ID and var2ID -
#' default = "var1" & "var2"; (4) color for lowest var1 - default = "#FF8C00";
#' (5) color for highest var1 - default = "#4682B4"; (6) color for lowest var2 -
#' default = "#FFFFFF", (7) color for highest var2 - default = "#F0E68C", (8)
#' color of co-orthologs - default = "#07D000"; (9+10+11) text sizes for x, y
#' axis and legend - default = 9 for each; (12) legend position "top", "bottom",
#' "right", "left" or "none" - default = "top"; (13) zoom ratio of the
#' co-ortholog dots from -1 to 3 - default = 0; (14) angle of x-axis from 0 to
#' 90 - default = 60; (14) show/hide separate line for reference taxon 1/0 -
#' default = 0; (15) enable/disable coloring gene categories TRUE/FALSE -
#' default = FALSE). NOTE: Leave blank or NULL to use default values.
#' @param taxonHighlight taxon of interst. Default = "none".
#' @param rankName working taxonomy rank (needed only for highlight taxon).
#' @param geneHighlight gene of interest. Default = "none".
#' @return A profile heatmap plot with highlighted gene and/or taxon of interest
#' as ggplot object.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @importFrom utils data
#' @seealso \code{\link{dataMainPlot}}, \code{\link{dataCustomizedPlot}}
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(fullProcessedProfile)
#' plotParameter <- list(
#'     "xAxis" = "taxa",
#'     "var1ID" = "FAS_FW",
#'     "var2ID"  = "FAS_BW",
#'     "lowColorVar1" =  "#FF8C00",
#'     "highColorVar1" = "#4682B4",
#'     "lowColorVar2" = "#FFFFFF",
#'     "highColorVar2" = "#F0E68C",
#'     "paraColor" = "#07D000",
#'     "xSize" = 8,
#'     "ySize" = 8,
#'     "legendSize" = 8,
#'     "mainLegend" = "top",
#'     "dotZoom" = 0,
#'     "xAngle" = 60,
#'     "guideline" = 0,
#'     "colorByGroup" = FALSE
#' )
#' taxonHighlight <- "none"
#' rankName <- "class"
#' geneHighlight <- "100265at6656"
#' highlightProfilePlot(
#'     plotDf, plotParameter, taxonHighlight, rankName, geneHighlight
#' )

highlightProfilePlot <- function(
    data = NULL, plotParameter = NULL, taxonHighlight = "none",
    rankName = "none", geneHighlight = "none"
){
    if (is.null(data)) stop("Input data cannot be NULL!")
    xmin <- xmax <- ymin <- ymax <- NULL
    p <- heatmapPlotting(data, plotParameter)
    # highlight taxon
    if (taxonHighlight != "none") {
        # get selected highlight taxon ID
        nameReducedFile <- paste(
            system.file(package = "PhyloProfile"),
            "PhyloProfile/data/taxonNamesReduced.txt", sep="/")
        if (!file.exists(nameReducedFile)) {
            taxonNamesReduced <- NULL
            delayedAssign("taxName", taxonNamesReduced)
        } else
            taxName <- read.table(nameReducedFile, sep = "\t", header = TRUE)
        taxonHighlightID <- taxName$ncbiID[
            taxName$fullName == taxonHighlight & taxName$rank == rankName]
        if (length(taxonHighlightID) == 0L)
            taxonHighlightID <- taxName$ncbiID[taxName$fullName==taxonHighlight]
        # get taxonID together with it sorted index
        selTaxon <- toString(data[data$supertaxonID == taxonHighlightID, 2][1])
        selIndex <- grep(selTaxon, levels(as.factor(data$supertaxon)))
        if (plotParameter$xAxis == "taxa") {
            rect <- data.frame(
                xmin=selIndex-0.5, xmax = selIndex+0.5, ymin = -Inf, ymax = Inf)
        } else
            rect <- data.frame(
                ymin=selIndex-0.5, ymax = selIndex+0.5, xmin = -Inf, xmax = Inf)
        p <- heatmapPlotting(data, plotParameter) + geom_rect(
            data = rect, color = "yellow", alpha = 0.3, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
    }
    # highlight gene
    if (geneHighlight != "none") {
        selIndex <- match(geneHighlight, levels(as.factor(data$geneID)))
        if (plotParameter$xAxis == "taxa") {
            rect <- data.frame(
                ymin=selIndex-0.5, ymax = selIndex+0.5, xmin = -Inf, xmax = Inf)
        } else
            rect <- data.frame(
                xmin=selIndex-0.5, xmax = selIndex+0.5, ymin = -Inf, ymax = Inf)
        p <- heatmapPlotting(data, plotParameter) + geom_rect(
            data = rect, color = "yellow", alpha = 0.3, inherit.aes = FALSE,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
    }
    return(p)
}
