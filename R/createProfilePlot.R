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
#' @rawNamespace import(data.table, except = c(set, melt))
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{filterProfileData}}
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' dataMainPlot(fullProcessedProfile)

dataMainPlot <- function(dataHeat = NULL){
    if (is.null(dataHeat)) return()
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
#' selectedTaxa <- c("Mammalia", "Echinoidea", "Gunneridae")
#' selectedSeq <- "all"
#' dataCustomizedPlot(fullProcessedProfile, selectedTaxa, selectedSeq)

dataCustomizedPlot <- function(
    dataHeat = NULL, selectedTaxa = "all", selectedSeq = "all"
){
    if (is.null(dataHeat)) return()
    geneID <- NULL
    supertaxonMod <- NULL
    paralogNew <- NULL

    # process data
    dataHeat$supertaxonMod <- {
        substr(
            dataHeat$supertaxon, 6, nchar(as.character(dataHeat$supertaxon))
        )
    }

    if (selectedTaxa[1] == "all" & selectedSeq[1] != "all") {
        # select data from dataHeat for selected sequences only
        dataHeat <- subset(dataHeat, geneID %in% selectedSeq)
    } else if (selectedSeq[1] == "all" & selectedTaxa[1] != "all") {
        # select data from dataHeat for selected taxa only
        dataHeat <- subset(dataHeat, supertaxonMod %in% selectedTaxa)
    } else {
        # select data from dataHeat for selected sequences and taxa
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

    return(dataHeat)
}


#' Create profile heatmap plot
#' @export
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
#' @return A profile heatmap plot as a ggplot object.
#' @importFrom plyr mapvalues
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 scale_color_gradient
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 guide_colourbar
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_rect
#' @importFrom ggplot2 geom_tile
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{dataMainPlot}}, \code{\link{dataCustomizedPlot}}
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(fullProcessedProfile)
#' plotParameter <- list(
#'     "xAxis" = "taxa",
#'     "var1ID" = "FAS",
#'     "var2ID"  = "Traceability",
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

heatmapPlotting <- function(data = NULL, plotParameter = NULL){
    if (is.null(data)) return()
    if (is.null(plotParameter)) {
        plotParameter <- list(
            "xAxis" = "taxa",
            "var1ID" = "var1", "var2ID"  = "var2",
            "lowColorVar1" =  "#FF8C00", "highColorVar1" = "#4682B4",
            "lowColorVar2" = "#FFFFFF", "highColorVar2" = "#F0E68C",
            "paraColor" = "#07D000",
            "xSize" = 8, "ySize" = 8,
            "legendSize" = 8, "mainLegend" = "top",
            "dotZoom" = 0,
            "xAngle" = 60,
            "guideline" = 0,
            "colorByGroup" = FALSE
        )
    }
    
    geneID <- NULL
    supertaxon <- NULL
    group <- NULL
    var1 <- NULL
    var2 <- NULL
    presSpec <- NULL
    paralog <- NULL
    xmin <- NULL
    xmax <- NULL
    ymin <- NULL
    ymax <- NULL

    # parameters
    xAxis <- plotParameter$xAxis
    var1ID <- plotParameter$var1ID
    var2ID <- plotParameter$var2ID
    lowColorVar1 <- plotParameter$lowColorVar1
    highColorVar1 <- plotParameter$highColorVar1
    lowColorVar2 <- plotParameter$lowColorVar2
    highColorVar2 <- plotParameter$highColorVar2
    paraColor <- plotParameter$paraColor
    xSize <- plotParameter$xSize
    ySize <- plotParameter$ySize
    legendSize <- plotParameter$legendSize
    mainLegend <- plotParameter$mainLegend
    dotZoom <- plotParameter$dotZoom
    xAngle <- plotParameter$xAngle
    guideline <- plotParameter$guideline
    colorByGroup <- plotParameter$colorByGroup

    # rescale numbers of paralogs
    if (length(unique(na.omit(data$paralog))) > 0) {
        maxParalog <- max(na.omit(data$paralog))
        data$paralogSize <- (data$paralog / maxParalog) * 3
    }

    # remove prefix number of taxa names but keep the order
    data$supertaxon <- {
        mapvalues(
            warn_missing = FALSE,
            data$supertaxon,
            from = as.character(data$supertaxon),
            to = substr(
                as.character(data$supertaxon),
                6,
                nchar(as.character(data$supertaxon))
            )
        )
    }

    # format plot
    if (xAxis == "genes") {
        p <- ggplot(data, aes(x = geneID, y = supertaxon))
    } else{
        p <- ggplot(data, aes(y = geneID, x = supertaxon))
    }

    if (colorByGroup == TRUE) {
        p <- p + geom_tile(aes(fill = factor(group)), alpha = 0.3)
    } else {
        if (length(unique(na.omit(data$var2))) != 1) {
            p <- p + scale_fill_gradient(
                low = lowColorVar2,
                high = highColorVar2,
                na.value = "gray95",
                limits = c(0, 1)
            ) +  #fill color (var2)
                geom_tile(aes(fill = var2))    # filled rect (var2 score)
        }
    }

    if (length(unique(na.omit(data$presSpec))) < 3) {
        if (length(unique(na.omit(data$var1))) == 1) {
            # geom_point for circle illusion (var1 and presence/absence)
            p <- p + geom_point(aes(colour = var1),
                                size = data$presSpec * 5 * (1 + dotZoom),
                                na.rm = TRUE, show.legend = FALSE)
        } else {
            # geom_point for circle illusion (var1 and presence/absence)
            p <- p + geom_point(aes(colour = var1),
                                size = data$presSpec * 5 * (1 + dotZoom),
                                na.rm = TRUE)
            # color of the corresponding aes (var1)
            p <- p + scale_color_gradient(
                low = lowColorVar1, high = highColorVar1, limits = c(0, 1)
            )
        }
    } else {
        if (length(unique(na.omit(data$var1))) == 1) {
            # geom_point for circle illusion (var1 and presence/absence)
            p <- p + geom_point(aes(size = presSpec), color = "#336a98",
                                na.rm = TRUE)
        } else {
            # geom_point for circle illusion (var1 and presence/absence)
            p <- p + geom_point(aes(colour = var1, size = presSpec),
                                na.rm = TRUE)
            # color of the corresponding aes (var1)
            p <- p +
                scale_color_gradient(
                    low = lowColorVar1, high = highColorVar1, limits = c(0, 1)
                )
        }
    }

    # plot inparalogs (if available)
    if (length(unique(na.omit(data$paralog))) > 0) {
        p <- p + geom_point(
            data = data, aes(size = paralog), color = paraColor,
            na.rm = TRUE, show.legend = TRUE
        )
        p <- p + guides(size = guide_legend(title = "# of co-orthologs"))

        # to tune the size of circles
        p <- p +
            scale_size_continuous(
                range = c(
                    min(na.omit(data$paralogSize)) * (1 + dotZoom),
                    max(na.omit(data$paralogSize)) * (1 + dotZoom)
                )
            )
    } else {
        # remain the scale of point while filtering
        presentVl <- data$presSpec[!is.na(data$presSpec)]

        # to tune the size of circles;
        # use "floor(value*10)/10" to round "down" the value with one decimal nr
        p <- p +
            scale_size_continuous(
                range = c(
                    (floor(min(presentVl) * 10) / 10 * 5) * (1 + dotZoom),
                    (floor(max(presentVl) * 10) / 10 * 5) * (1 + dotZoom)
                )
            )
    }

    if (colorByGroup == FALSE) {
        p <- p + guides(fill = guide_colourbar(title = var2ID),
                        color = guide_colourbar(title = var1ID))
    } else {
        p <- p + guides(fill = guide_legend("Category"),
                        color = guide_colourbar(title = var1ID))
    }

    baseSize <- 9

    # guideline for separating ref species
    if (guideline == 1) {
        if (xAxis == "genes") {
            p <- p + labs(y = "Taxon")
            p <- p + geom_hline(yintercept = 0.5, colour = "dodgerblue4")
            p <- p + geom_hline(yintercept = 1.5, colour = "dodgerblue4")
        } else{
            p <- p + labs(x = "Taxon")
            p <- p + geom_vline(xintercept = 0.5, colour = "dodgerblue4")
            p <- p + geom_vline(xintercept = 1.5, colour = "dodgerblue4")
        }
    }

    # format theme
    p <- p + theme_minimal()
    p <- p + theme(
        axis.text.x = element_text(angle = xAngle, hjust = 1, size = xSize),
        axis.text.y = element_text(size = ySize),
        axis.title.x = element_text(size = xSize),
        axis.title.y = element_text(size = ySize),
        legend.title = element_text(size = legendSize),
        legend.text = element_text(size = legendSize),
        legend.position = mainLegend
    )

    # return plot
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
#' @seealso \code{\link{dataMainPlot}}, \code{\link{dataCustomizedPlot}}
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' plotDf <- dataMainPlot(fullProcessedProfile)
#' plotParameter <- list(
#'     "xAxis" = "taxa",
#'     "var1ID" = "FAS",
#'     "var2ID"  = "Traceability",
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
#' taxonHighlight <- "Mammalia"
#' rankName <- "class"
#' geneHighlight <- "OG_1019"
#' highlightProfilePlot(
#'     plotDf, plotParameter, taxonHighlight, rankName, geneHighlight
#' )

highlightProfilePlot <- function(
    data = NULL,
    plotParameter = NULL,
    taxonHighlight = "none",
    rankName = "none",
    geneHighlight = "none"
){
    if (is.null(data)) return()
    xmin <- NULL
    xmax <- NULL
    ymin <- NULL
    ymax <- NULL

    # get heatmap
    p <- heatmapPlotting(data, plotParameter)

    # highlight taxon
    if (taxonHighlight != "none") {
        # get selected highlight taxon ID
        nameReducedFile <- paste(
            system.file(package="PhyloProfile"),
            "PhyloProfile/data/taxonNamesReduced.txt",
            sep="/"
        )

        if (!file.exists(nameReducedFile)) {
            fileURL <- paste0(
                "https://raw.githubusercontent.com/BIONF/phyloprofile-data/",
                "master/taxonNamesReduced.txt"
            )
            res <- tryCatch(
                utils::download.file(
                    fileURL, destfile = nameReducedFile, method="auto"
                ),
                error=function(e) 1
            )
        }

        taxaList <- read.table(nameReducedFile, sep = "\t", header = TRUE)

        taxonHighlightID <- {
            taxaList$ncbiID[
                taxaList$fullName == taxonHighlight & taxaList$rank == rankName
            ]
        }

        if (length(taxonHighlightID) == 0L) {
            taxonHighlightID <- {
                taxaList$ncbiID[taxaList$fullName == taxonHighlight]
            }
        }

        # get taxonID together with it sorted index
        highlightTaxon <- {
            toString(
                data[data$supertaxonID == taxonHighlightID, 2][1]
            )
        }

        # get index
        selectedIndex <- as.numeric(
            as.character(substr(highlightTaxon, 2, 4))
        )

        # draw a rect to highlight this taxon's column
        if (plotParameter$xAxis == "taxa") {
            rect <- data.frame(
                xmin = selectedIndex - 0.5,
                xmax = selectedIndex + 0.5,
                ymin = -Inf,
                ymax = Inf
            )
        } else {
            rect <- data.frame(
                ymin = selectedIndex - 0.5,
                ymax = selectedIndex + 0.5,
                xmin = -Inf,
                xmax = Inf
            )
        }

        p <- p + geom_rect(
            data = rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = "yellow",
            alpha = 0.3,
            inherit.aes = FALSE
        )
    }

    # highlight gene
    if (geneHighlight != "none") {
        # get selected highlight gene ID
        geneHighlight <- geneHighlight

        # get index
        allGenes <- levels(data$geneID)
        selectedIndex <- match(geneHighlight, allGenes)

        # draw a rect to highlight this taxon's column
        if (plotParameter$xAxis == "taxa") {
            rect <- data.frame(
                ymin = selectedIndex - 0.5,
                ymax = selectedIndex + 0.5,
                xmin = -Inf,
                xmax = Inf
            )
        } else {
            rect <- data.frame(
                xmin = selectedIndex - 0.5,
                xmax = selectedIndex + 0.5,
                ymin = -Inf,
                ymax = Inf
            )
        }

        p <- p + geom_rect(
            data = rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            color = "yellow",
            alpha = 0.3,
            inherit.aes = FALSE
        )
    }

    return(p)
}
