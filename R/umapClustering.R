#' Prepare data for dimension reduction
#' @export
#' @usage prepareDimRedData(longDf = NULL, taxonRank = NULL, type = "taxa",
#'     taxDB = NULL, filterVar = "both", cutoff = 0, groupLabelsBy = "taxa")
#' @param longDf input phyloprofile file in long format
#' @param taxonRank taxonomy rank for labels (e.g. "phylum")
#' @param type type of clustering, either "taxa" (default) or "genes"
#' @param taxDB path to taxonomy database
#' @param filterVar choose variable (either "var1", "var2" or "both") to filter
#' the data. Default: "both"
#' @param cutoff cutoff to filter data values. Default: 0
#' @param groupLabelsBy group labels by the number of "taxa" (default) or
#' "genes"
#' @return A dataframe in wide format
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @examples
#' rawInput <- system.file(
#'    "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' longDf <- createLongMatrix(rawInput)
#' prepareDimRedData(longDf, "phylum")

prepareDimRedData <- function(
    longDf = NULL, taxonRank = NULL, type = "taxa", taxDB = NULL,
    filterVar = "both", cutoff = 0, groupLabelsBy = "taxa"
) {
    if (is.null(longDf)) stop("Input data cannot be NULL")
    if (is.null(taxonRank)) stop("Taxon rank must be specified!")
    FAS_F <- FAS_B <- geneID <- ncbiID <- abbrName <- fullName <- NULL
    var1 <- var2 <- Freq <- Label <- NULL
    if (type == "genes") groupLabelsBy <- "genes"
    # Handle column renaming and filtering based on input dimensions
    colNames <- c("geneID", "ncbiID", "orthoID", "var1", "var2", "geneName")
    colNames <- colNames[seq_len(ncol(longDf))]
    colnames(longDf) <- colNames
    # Add default values for missing columns
    if (!"var1" %in% colnames(longDf)) longDf$var1 <- 1
    # Apply filters based on cutoff values
    filterFlag <- ifelse("var1" %in% colnames(longDf), 1, 0)
    if (filterFlag) {
        longDf <- longDf %>%
            filter(
                (filterVar == "both" & var1 >= cutoff & var2 >= cutoff) |
                    (filterVar == "var1" & var1 >= cutoff) |
                    (filterVar == "var2" & var2 >= cutoff)
            )
    }
    # calculate mean of var1 and var2 (if var2 available)
    if (all(c("var1", "var2") %in% colnames(longDf)))
        longDfSub <- longDf %>% mutate(var1 = (var1 + var2) / 2)
    longDfSub <- longDfSub %>% select(geneID, ncbiID, var1)

    # get taxon names for input taxa based on a selected supertaxon rank
    taxMatrix <- getTaxonomyMatrix(taxDB)
    nameList <- getNameList(taxDB)
    superTaxonDf <- taxMatrix %>%
        filter(abbrName %in% longDfSub$ncbiID) %>%
        select(abbrName, one_of(taxonRank))
    colnames(superTaxonDf) <- c("ncbiID", "superID")
    superTaxonDf <- dplyr::left_join(
        superTaxonDf, nameList, by = c("superID" = "ncbiID")
    )
    # transform to wide format, add taxon names and ortholog count
    # if co-orthologs present, get max value (e.g. FAS) as the representative
    if (type == "taxa") {
        wideDf <- data.table::dcast(
            setDT(longDfSub), ncbiID ~ geneID, value.var = c("var1"),
            fun.aggregate = max, fill = -1
        )
        wideDf <- dplyr::left_join(
            wideDf, superTaxonDf %>% select(ncbiID, Label = fullName),
            by = "ncbiID"
        )
        # add count (how many seed genes each taxon has, excluding co-orthologs)
        seedWithTax <- longDfSub %>% select(geneID, ncbiID)
        seedWithTax <- seedWithTax[!duplicated(seedWithTax),]
        countDf <- data.frame(table(seedWithTax$ncbiID))
        colnames(countDf) <- c("ncbiID", "Freq")
        wideDf <- dplyr::left_join(
            wideDf, countDf, by = c("ncbiID")
        )
    } else {
        wideDf <- data.table::dcast(
            setDT(longDfSub), geneID ~ ncbiID, value.var = c("var1"),
            fun.aggregate = max, fill = -1
        )
        wideDf$Label <- wideDf$geneID
        # add count (how many taxa has orthologs for each seed gene)
        seedWithTax <- longDfSub %>% select(geneID, ncbiID)
        seedWithTax <- seedWithTax[!duplicated(seedWithTax),]
        countDf <- data.frame(table(seedWithTax$geneID))
        colnames(countDf) <- c("geneID", "Freq")
        wideDf <- dplyr::left_join(
            wideDf, countDf, by = c("geneID")
        )
    }
    # calculate genes / taxa frequency for each label
    if (groupLabelsBy == "genes") {
        countFreqDf <- data.frame(
            wideDf %>% group_by(Label) %>% summarise(n = max(Freq))
        )
    } else {
        countFreqDf <- data.frame(wideDf %>% count(Label))
    }
    wideDf <- left_join(wideDf, countFreqDf, by = "Label")
    wideDf$n <- as.numeric(wideDf$n)
    wideDf$Label <- as.character(wideDf$Label)
    return(data.frame(wideDf, check.names = FALSE))
}

#' Perform dimension reduction 2D
#' @export
#' @usage dimReduction(data4dimRed = NULL, by = "taxa", type = "binary",
#'     randomSeed = 123, reductionTechnique = "umap", dimension = "2d",
#'     tsneIter = 1000)
#' @param data4dimRed data for dimension reduction (from prepareDimRedData)
#' @param by cluster data by "taxa" (default) or "genes"
#' @param type type of data, either "binary" (default) or "non-binary"
#' @param randomSeed random seed. Default: 123
#' @param reductionTechnique dimensionality reduction technique, either "umap"
#' (default) or "tsne"
#' @param dimension either "2d" (default) or "3d"
#' @param tsneIter number of iterations for t-SNE. Default: 1000
#' @import umap tsne
#' @return A table contains coordinates of the 2D dimension reduction
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{prepareDimRedData}}
#' @examples
#' rawInput <- system.file(
#'    "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' longDf <- createLongMatrix(rawInput)
#' data4dimRed <- prepareDimRedData(longDf, "phylum")
#' dimReduction(data4dimRed)

dimReduction <- function(
    data4dimRed = NULL, by = "taxa", type = "binary", randomSeed = 123,
    reductionTechnique = "umap", dimension = "2d", tsneIter = 1000
) {
    if (is.null(data4dimRed)) stop("Input data cannot be NULL!")
    # Determine subset columns based on `by`
    subsetCols <- if (by == "taxa") {
        setdiff(colnames(data4dimRed), c("ncbiID", "Label", "Freq", "n"))
    } else {
        setdiff(colnames(data4dimRed), c("geneID", "Label", "Freq", "n"))
    }
    subsetDt <- data4dimRed[, subsetCols, drop = FALSE]
    # Convert to binary if required
    if (type == "binary") {
        subsetDt <- ifelse(subsetDt >= 0, 1, 0)
    }
    # Define the dimensionality
    dim <- ifelse(dimension == "2d", 2, 3)
    # Perform dimensionality reduction
    outDf <- switch(
        reductionTechnique,
        umap = {performUmap(subsetDt, randomSeed, dim)},
        tsne = {
            tsne::tsne(subsetDt, initial_dims=dim, k=dim, max_iter = tsneIter)
        },
        pca = {performPCA(subsetDt)},
        stop("Unsupported reduction technique: ", reductionTechnique)
    )
    return(outDf)
}

#' Helper function to handle UMAP logic
#' @param umapDt data matrix for UMAP
#' @param randomSeed random seed. Default: 123
#' @param dim dimension, either 2 for 2D (default) or 3 for 3D
#' @import umap
#' @return A table contains coordinates UMAP reduction
#' @author Vinh Tran tran@bio.uni-frankfurt.de

performUmap <- function(umapDt, randomSeed = 123, dim = 2) {
    # Attempt UMAP reduction with fallback for sample size issues
    tryCatch({
        suppressWarnings(
            umap::umap(
                umapDt, random_state = randomSeed, preserve.seed = TRUE,
                n_components = dim
            )$layout
        )
    }, error = function(cond) {
        message("UMAP failed: ", conditionMessage(cond))
        fallbackUmap(umapDt, randomSeed, dim)
    }, warning = function(cond) {
        message(conditionMessage(cond))
        fallbackUmap(umapDt, randomSeed, dim)
    })
}

#' Fallback for UMAP in case of insufficient samples
#' @param umapDt data matrix for UMAP
#' @param randomSeed random seed. Default: 123
#' @param dim dimension, either 2 for 2D (default) or 3 for 3D
#' @import umap
#' @return A table contains coordinates UMAP reduction
#' @author Vinh Tran tran@bio.uni-frankfurt.de

fallbackUmap <- function(umapDt, randomSeed, dim) {
    warning("PROBLEM: Too few samples for UMAP! Adjusting n_neighbors.")
    umap::umap(
        umapDt, random_state = randomSeed, preserve.seed = TRUE,
        n_neighbors = max(1, nrow(umapDt) - 1), n_components = dim
    )$layout
}


#' Helper function to perform PCA
#' @param pcaDt data matrix for PCA
#' @return A table contains coordinates of the first 3 PCs
#' @author Vinh Tran tran@bio.uni-frankfurt.de

performPCA <- function(pcaDt) {
    if (!is.data.frame(pcaDt)) as.data.frame(pcaDt)
    # Perform PCA
    prin_comp <- stats::prcomp(pcaDt)
    
    # Extract the importance (proportion of variance explained by each PC)
    importancePCs <- as.data.frame(prin_comp$sdev^2) / sum(prin_comp$sdev^2)
    colnames(importancePCs) <- "Proportion of Variance"
    # Print the Proportion of Variance for the first 3 PCs
    importancePCs <- importancePCs[seq_len(3), , drop = FALSE]
    # Check if the first 3 PCs explain at least 80% of the variance
    sumVariance <- sum(importancePCs$`Proportion of Variance`)
    if (sumVariance < 0.8) {
        msg <- paste0(
            "The first 3 PCs can only explain ", round(sumVariance*100),
            "% of the data variance (<80%)!"
        )
        warning(msg)
    }
    
    # Extract the principal components scores
    components <- as.data.frame(prin_comp$x)
    # Return only the first 3 principal components
    return(components[, seq_len(3), drop = FALSE])
}


#' Reduce the number of labels for DIM reduction plot based on the
#' gene/taxon frequency
#' @export
#' @usage groupLabelDimRedData(data4dimRed = NULL, freqCutoff = c(0,200))
#' @param data4dimRed data for dimension reduction (from prepareDimRedData)
#' @param freqCutoff gene/taxon frequency cutoff range. Any labels that are
#' outside of this range will be assigned as [Other]
#' @return A dataframe similar to input data4dimRed, but with modified Label
#' column, where less frequent labels are grouped together as "Other"
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{prepareDimRedData}}
#' @examples
#' rawInput <- system.file(
#'    "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' longDf <- createLongMatrix(rawInput)
#' data4dimRed <- prepareDimRedData(longDf, "phylum")
#' groupLabelDimRedData(data4dimRed, freqCutoff = c(3,5))

groupLabelDimRedData <- function(data4dimRed = NULL, freqCutoff = c(0,200)) {
    if (is.null(data4dimRed) || length(data4dimRed) == 0)
        stop("Input data cannot be NULL or EMPTY!")
    data4dimRed$Label[
        data4dimRed$n < freqCutoff[1] | data4dimRed$n > freqCutoff[2]
    ] <- "[Other]"
    return(data4dimRed)
}

#' Generate data for dimension reduction plot
#' @export
#' @usage createDimRedPlotData(dimRedCoord = NULL, data4dimRed = NULL,
#'     freqCutoff = c(0,200), excludeTaxa = "None", currentNCBIinfo = NULL)
#' @param dimRedCoord data contains DIM reduction coordinates (from
#' dimReduction)
#' @param data4dimRed data for dimension reduction (from prepareDimRedData())
#' @param freqCutoff gene/taxon frequency cutoff range. Any labels that are
#' outside of this range will be assigned as [Other]
#' @param excludeTaxa hide taxa from plot. Default: "None"
#' @param currentNCBIinfo table/dataframe of the pre-processed NCBI taxonomy
#' data (/PhyloProfile/data/preProcessedTaxonomy.txt)
#' @importFrom utils tail
#' @return A plot as ggplot object
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{prepareDimRedData}}, \code{\link{dimReduction}}
#' @examples
#' rawInput <- system.file(
#'    "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' longDf <- createLongMatrix(rawInput)
#' data4dimRed <- prepareDimRedData(longDf, "phylum")
#' dimRedCoord <- dimReduction(data4dimRed)
#' createDimRedPlotData(dimRedCoord, data4dimRed)

createDimRedPlotData <- function(
    dimRedCoord = NULL, data4dimRed = NULL, freqCutoff = c(0, 200),
    excludeTaxa = "None", currentNCBIinfo = NULL
) {
    if (is.null(dimRedCoord) || is.null(data4dimRed) ||
        length(dimRedCoord) == 0 || length(data4dimRed) == 0) {
        stop("Input data cannot be NULL or EMPTY!")
    }
    Label <- Freq <- NULL
    data4dimRed$X <- dimRedCoord[,1]
    data4dimRed$Y <- dimRedCoord[,2]
    if (ncol(dimRedCoord) == 3) data4dimRed$Z <- dimRedCoord[,3]
    # join less freq items into "other"
    data4dimRed <- groupLabelDimRedData(data4dimRed, freqCutoff)
    # exclude taxa
    if (length(excludeTaxa) > 0) {
        data4dimRed$X[data4dimRed$Label %in% excludeTaxa] <- NA
        data4dimRed$Y[data4dimRed$Label %in% excludeTaxa] <- NA
    }
    # convert tax IDs into names
    if ("ncbiID" %in% colnames(data4dimRed)) {
        if (is.null(currentNCBIinfo)) {
            dataPath <- system.file(
                "PhyloProfile", "data/",
                package = "PhyloProfile", mustWork = TRUE
            )
            nameFullFile <- paste0(dataPath, "/preProcessedTaxonomy.txt")
            if (file.exists(nameFullFile)) {
                currentNCBIinfo <-as.data.frame(data.table::fread(nameFullFile))
            }
        }
        if (!is.null(currentNCBIinfo)) {
            id2nameDf <- PhyloProfile::id2name(
                gsub("ncbi","",unique(data4dimRed$ncbiID)), currentNCBIinfo
            )
            id2nameDf$ncbiID <- paste0("ncbi", id2nameDf$ncbiID)
            data4dimRed <- merge(
                data4dimRed, id2nameDf, by = "ncbiID", all.x = TRUE
            )
        } else {
            data4dimRed$fullName <- "NA"
        }
    }
    return(data4dimRed)
}


#' Create dimension reduction plot
#' @export
#' @usage plotDimRed(plotDf = NULL, legendPos = "bottom",
#'     colorPalette = "Set2", transparent = 0, textSize = 12, font = "Arial",
#'     highlightTaxa = NULL, dotZoom = 0)
#' @param plotDf data for dimension reduction 2D plot
#' @param legendPos position of legend. Default: "right"
#' @param colorPalette color palette. Default: "Set2"
#' @param transparent transparent level (from 0 to 1). Default: 0
#' @param textSize size of axis and legend text. Default: 12
#' @param font font of text. Default = Arial"
#' @param highlightTaxa list of taxa to be highlighted
#' @param dotZoom dot size zooming factor. Default: 0
#' @return A plot as ggplot object
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{prepareDimRedData}}, \code{\link{dimReduction}},
#' \code{\link{createDimRedPlotData}}
#' @examples
#' rawInput <- system.file(
#'    "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' longDf <- createLongMatrix(rawInput)
#' data4dimRed <- prepareDimRedData(longDf, "phylum")
#' dimRedCoord <- dimReduction(data4dimRed)
#' plotDf <- createDimRedPlotData(dimRedCoord, data4dimRed)
#' plotDimRed(plotDf, font = "sans")

plotDimRed <- function(
    plotDf = NULL, legendPos = "bottom", colorPalette = "Set2",
    transparent = 0, textSize = 12, font = "Arial", highlightTaxa = NULL,
    dotZoom = 0
) {
    if (is.null(plotDf)) stop("Input data cannot be NULL!")
    X <- Y <- Label <- Freq <- NULL
    # add colors
    plotDf <- addDimRedTaxaColors(plotDf, colorPalette, highlightTaxa)
    if ("color" %in% colnames(plotDf)) {
        colorScheme <- structure(
            plotDf$color, .Names = plotDf$Label
        )
    } else {
        colorScheme <- structure(
            head(
                suppressWarnings(
                    RColorBrewer::brewer.pal(
                        nlevels(as.factor(plotDf$Label)), colorPalette
                    )
                ), levels(as.factor(plotDf$Label))
            ), .Names = levels(as.factor(plotDf$Label))
        )
    }
    # generate plot
    plot <- ggplot(plotDf, aes(x = X, y = Y, color = Label)) +
        geom_point(aes(size = Freq), alpha = 1 - transparent) +
        geom_rug(alpha = 1) +
        theme_minimal() +
        labs(x = "", y = "")+
        scale_radius(range = c(1 + dotZoom, 6 + dotZoom))
    # change legend title
    if ("ncbiID" %in% colnames(plotDf)) {
        plot <- plot + guides(
            color = guide_legend(override.aes = list(alpha = 1), ncol = 5),
            size = guide_legend(title = "Number of genes", ncol = 2)
        )
    } else
        plot <- plot + guides(
            color = guide_legend(override.aes = list(alpha = 1), ncol = 5),
            size = guide_legend(title = "Number of taxa", ncol = 1)
        )
    plot <- plot + theme(
            legend.position = legendPos,
            legend.text = element_text(size = textSize),
            legend.title = element_text(size = textSize),
            axis.text = element_text(size = textSize),
            axis.title = element_text(size = textSize),
        )
    plot <- plot + scale_color_manual(values = colorScheme)
    plot <- plot + theme(text = element_text(family = font))
    return(plot)
}


#' Create dimension reduction 3D plot
#' @export
#' @usage plotDimRed3D(plotDf = NULL, legendPos = "bottom",
#'     colorPalette = "Set2", transparent = 0,highlightTaxa = NULL,
#'     dotZoom = 0)
#' @param plotDf data for dimension reduction 3D plot
#' @param legendPos position of legend. Default: "right"
#' @param colorPalette color palette. Default: "Set2"
#' @param transparent transparent level (from 0 to 1). Default: 0
#' @param highlightTaxa list of taxa to be highlighted
#' @param dotZoom dot size zooming factor. Default: 0
#' @return A plot as ggplot object
#' @rawNamespace import(plotly, except = c(last_plot, select))
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{prepareDimRedData}}, \code{\link{dimReduction}},
#' \code{\link{createDimRedPlotData}}
#' @examples
#' rawInput <- system.file(
#'    "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' longDf <- createLongMatrix(rawInput)
#' data4dimRed <- prepareDimRedData(longDf, "phylum")
#' dimRedCoord3d <- dimReduction(data4dimRed, dimension = "3d")
#' plotDf <- createDimRedPlotData(dimRedCoord3d, data4dimRed)
#' plotDimRed3D(plotDf)

plotDimRed3D <- function(
    plotDf = NULL, legendPos = "bottom", colorPalette = "Set2",
    transparent = 0, highlightTaxa = NULL, dotZoom = 0
) {
    if (is.null(plotDf)) stop("Input data cannot be NULL!")
    X <- Y <- Z <- Label <- Freq <- color <- NULL
    if (!("Z" %in% colnames(plotDf))) stop("PlotDf seems to be 2D plot!")
    # add colors
    plotDf <- addDimRedTaxaColors(plotDf, colorPalette, highlightTaxa)
    colorDf <- unique(plotDf %>% select(Label, color))
    # calculate min and max dot size
    minSize <- dotZoom
    maxSize <- (max(plotDf$Freq)*dotZoom/min(plotDf$Freq))
    # generate 3D scatter plot
    if ("ncbiID" %in% colnames(plotDf)) {
        plot <- plot_ly(
            plotDf, x = ~X, y = ~Y, z = ~Z,
            hovertemplate = ~paste0(
                "Taxon ID: ", gsub("ncbi","",ncbiID), "<br>Name: ", fullName,
                "<br>Label: ", Label, "<br>Freq: ", Freq
            )
        )
    } else {
        plot <- plot_ly(
            plotDf, x = ~X, y = ~Y, z = ~Z,
            hovertemplate = ~paste0(
                "Gene ID: ", geneID, "<br>Freq: ", Freq
            )
        )
    }

    if (dotZoom == 0) {
        plot <- plot %>% add_markers(
            color = ~Label, colors = colorDf$color, size = ~Freq,
            alpha = 1 - transparent, span = I(0)
        )
    } else {
        plot <- plot %>% add_markers(
            color = ~Label, colors = colorDf$color, size = ~Freq,
            sizes = c(minSize, maxSize), alpha = 1 - transparent, span = I(0)
        )
    }
    if (legendPos == "none") plot <- plot %>% hide_legend()
    return(plot |> layout(dragmode = "lasso"))
}


#' Add colors for taxa in dimension reduction plot
#' @usage addDimRedTaxaColors(plotDf = NULL, colorPalette = "Set2",
#'     highlightTaxa = NULL)
#' @param plotDf data for dimension reduction plot
#' @param colorPalette color palette. Default: "Set2"
#' @param highlightTaxa list of taxa to be highlighted
#' @return A dataframe for dimension reduction plot with an additional column
#' for the assigned color to each taxon
#' @author Vinh Tran tran@bio.uni-frankfurt.de
#' @seealso \code{\link{prepareDimRedData}}, \code{\link{dimReduction}},
#' \code{\link{createDimRedPlotData}}
#' @examples
#' rawInput <- system.file(
#'    "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' longDf <- createLongMatrix(rawInput)
#' data4dimRed <- prepareDimRedData(longDf, "phylum")
#' dimRedCoord <- dimReduction(data4dimRed)
#' plotDf <- createDimRedPlotData(dimRedCoord, data4dimRed)
#' PhyloProfile:::addDimRedTaxaColors(plotDf, colorPalette = "Set2")

addDimRedTaxaColors <- function(
    plotDf = NULL, colorPalette = "Set2", highlightTaxa = NULL
) {
    if (is.null(plotDf)) stop("plotDf cannot be null!")
    allTaxa <- levels(as.factor(plotDf$Label))
    allColors <- getQualColForVector(allTaxa)
    if (checkColorPalette(allTaxa, colorPalette)) {
        allColors <-
            suppressWarnings(head(
                RColorBrewer::brewer.pal(length(allTaxa), colorPalette),
                length(allTaxa)
            ))
    }
    colorScheme <- data.frame(color = allColors, Label = allTaxa)
    if (!(is.null(highlightTaxa)))
        colorScheme$color[!(colorScheme$Label %in% highlightTaxa)] <- "#d4d6d9"
    plotDf <- merge(plotDf, colorScheme, by = "Label", all.x = TRUE)
    return(plotDf)
}
