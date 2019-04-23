#' Create data for percentage present taxa distribution
#' @param inputData dataframe contains raw input data in long format
#' @param rankName name of the working taxonomy rank (e.g. "species", "family")
#' @return A dataframe ready for analysing the distribution of the percentage of
#' species in the selected supertaxa
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{mainLongRaw}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' createPercentageDistributionData(mainLongRaw, "class")

createPercentageDistributionData <- function(inputData, rankName) {
    mdData <- inputData
    if (ncol(mdData) < 4) {
        colnames(mdData) <- c("geneID", "ncbiID", "orthoID")
    } else if (ncol(mdData) < 5) {
        colnames(mdData) <- c("geneID", "ncbiID", "orthoID", "var1")
    } else {
        colnames(mdData) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
    }

    # count number of inparalogs
    paralogCount <- plyr::count(mdData, c("geneID", "ncbiID"))
    mdData <- merge(mdData, paralogCount, by = c("geneID", "ncbiID"))
    colnames(mdData)[ncol(mdData)] <- "paralog"

    # (3) GET SORTED TAXONOMY LIST (3)
    inputTaxonID <- getInputTaxaID(inputData)
    inputTaxonName <- getInputTaxaName(
        rankName, inputTaxonID
    )
    refTaxon <- inputTaxonName$fullName[1]
    taxaTree <- NULL

    taxaList <- sortInputTaxa(
        inputTaxonID, inputTaxonName, rankName, refTaxon, taxaTree
    )

    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxaList, "supertaxon")

    # merge mdData, mdDatavar2 and taxaList to get taxonomy info
    taxaMdData <- merge(mdData, taxaList, by = "ncbiID")
    
    # calculate % present species
    finalPresSpecDt <- calcPresSpec(taxaMdData, taxaCount)

    finalPresSpecDt[!is.na(finalPresSpecDt$geneID),]
    return(finalPresSpecDt)
}

#' Create data for additional variable distribution
#' @usage createVariableDistributionData(inputData, var1CutoffMin,
#'     var1CutoffMax, var2CutoffMin, var2CutoffMax)
#' @param inputData dataframe contains raw input data in long format
#' @param var1CutoffMin min cutoff for var1
#' @param var1CutoffMax max cutoff for var1
#' @param var2CutoffMin min cutoff for var2
#' @param var2CutoffMax max cutoff for var2
#' @return A dataframe ready for analysing the distribution of the additional
#' variable(s)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{mainLongRaw}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' createVariableDistributionData(
#'     mainLongRaw, 0, 1, 0.5, 1
#' )

createVariableDistributionData <- function(
    inputData,
    var1CutoffMin,
    var1CutoffMax,
    var2CutoffMin,
    var2CutoffMax
) {
    dataOrig <- inputData
    if (ncol(dataOrig) < 4) {
        colnames(dataOrig) <- c("geneID",
                                "ncbiID",
                                "orthoID")
        splitDt <- dataOrig[, c("orthoID")]
    } else if (ncol(dataOrig) < 5) {
        colnames(dataOrig) <- c("geneID",
                                "ncbiID",
                                "orthoID",
                                "var1")
        splitDt <- dataOrig[, c("orthoID",
                                "var1")]
    } else {
        colnames(dataOrig) <- c("geneID",
                                "ncbiID",
                                "orthoID",
                                "var1",
                                "var2")
        splitDt <- dataOrig[, c("orthoID", "var1", "var2")]
    }

    splitDt$orthoID[splitDt$orthoID == "NA" | is.na(splitDt$orthoID)] <- NA
    splitDt <- splitDt[complete.cases(splitDt), ]

    if (length(levels(as.factor(splitDt$var2))) == 1) {
        if (levels(as.factor(splitDt$var2)) == "") {
            splitDt$var2 <- 0
        }
    }

    # Filter based on variable cutoffs
    if ("var1" %in% colnames(splitDt)) {
        # filter splitDt based on selected var1 cutoff
        splitDt <- splitDt[splitDt$var1 >= var1CutoffMin
                            & splitDt$var1 <= var1CutoffMax, ]
    }
    if ("var2" %in% colnames(splitDt)) {
        # filter splitDt based on selected var2 cutoff
        splitDt <- splitDt[splitDt$var2 >= var2CutoffMin
                            & splitDt$var2 <= var2CutoffMax, ]
    }

    return(splitDt)
}

#' Create data for additional variable distribution (for a subset data)
#' @usage createVariableDistributionDataSubset(fullProfileData,
#'     distributionData, selectedGenes, selectedTaxa)
#' @param fullProfileData dataframe contains the full processed profiles
#' @param distributionData dataframe contains the full distribution data
#' @param selectedGenes list of genes of interst
#' @param selectedTaxa list of taxa of interest
#' @return A dataframe ready for analysing the distribution of the additional
#' variable(s) for a subset of genes and/or taxa.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{parseInfoProfile}},
#' \code{\link{createVariableDistributionData}},
#' \code{\link{fullProcessedProfile}}, \code{\link{mainLongRaw}}
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' data("mainLongRaw", package="PhyloProfile")
#' distributionData <- createVariableDistributionData(
#'     mainLongRaw, 0, 1, 0.5, 1
#' )
#' selectedGenes <- "OG_1019"
#' selectedTaxa <- c("Mammalia", "Echinoidea", "Gunneridae", "Mucorales",
#' "Alphaproteobacteria")
#' createVariableDistributionDataSubset(
#'     fullProcessedProfile,
#'     distributionData,
#'     selectedGenes,
#'     selectedTaxa
#' )

createVariableDistributionDataSubset <- function(
    fullProfileData,
    distributionData,
    selectedGenes,
    selectedTaxa
) {
    geneID <- NULL
    orthoID <- NULL
    var1.x <- NULL
    var2.y <- NULL
    supertaxonMod <- NULL

    allData <- fullProfileData
    splitDt <- distributionData

    # get geneID and supertaxon name for splitDt
    splitDtName <- merge(
        splitDt, allData,
        by = "orthoID",
        all.x = TRUE
    )
    splitDtName$supertaxonMod <-
        substr(
            splitDtName$supertaxon,
            6,
            nchar(as.character(splitDtName$supertaxon))
        )
    splitDtName <- subset(
        splitDtName,
        select = c(
            orthoID,
            var1.x,
            var2.y,
            supertaxonMod,
            geneID
        )
    )
    colnames(splitDtName) <- c(
        "orthoID",
        "var1",
        "var2",
        "supertaxonMod",
        "geneID"
    )

    # filter
    if (selectedTaxa[1] == "all" & selectedGenes[1] != "all") {
        # select data from dataHeat for selected sequences only
        splitDt <- subset(splitDtName, geneID %in% selectedGenes)
    } else if (selectedGenes[1] == "all" & selectedTaxa[1] != "all") {
        # select data from dataHeat for selected taxa only
        splitDt <- subset(splitDtName, supertaxonMod %in% selectedTaxa)
    } else {
        # select data from dataHeat for selected sequences and taxa
        splitDt <- subset(
            splitDtName,
            geneID %in% selectedGenes
            & supertaxonMod %in% selectedTaxa
        )
    }

    return(splitDt)
}

#' Create distribution plot
#' @description Create distribution plot for one of the additional variable or
#' the percentage of the species present in the supertaxa.
#' @param data dataframe contains data for plotting
#' @param varName name of the variable that need to be analyzed (either name of
#' variable 1 or variable 2 or "percentage of present taxa")
#' @param varType type of variable (either "var1", "var2" or "presSpec")
#' @param percent range of percentage cutoff
#' @param distTextSize text size of the distribution plot
#' @return A distribution plot as a ggplot object
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 geom_vline
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{mainLongRaw}},
#' \code{\link{createVariableDistributionData}},
#' \code{\link{createVariableDistributionDataSubset}},
#' \code{\link{createPercentageDistributionData}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- createVariableDistributionData(
#'     mainLongRaw, 0, 1, 0.5, 1
#' )
#' varName <- "Variable abc"
#' varType <- "var1"
#' percent <- c(0,1)
#' distTextSize <- 12
#' createVarDistPlot(
#'     data,
#'     varName,
#'     varType,
#'     percent,
#'     distTextSize
#' )

createVarDistPlot <- function(
    data,
    varName,
    varType,
    percent,
    distTextSize
) {
    if (varType == "presSpec") {
        # remove presSpec < cutoffMin or > cutoffMax
        if (percent[1] > 0) {
            data <- data[data$presSpec >= percent[1]
                        & data$presSpec <= percent[2], ]
        } else {
            if (percent[2] > 0) {
                data <- data[data$presSpec > 0 & data$presSpec <= percent[2], ]
            } else {
                data <- data[data$presSpec > 0, ]
            }
        }
    } else {
        data <- data[!is.na(data[,varType]), ]
    }

    data.mean <- mean(data[,varType])

    p <- ggplot(data, aes(x = data[,varType])) +
        geom_histogram(binwidth = .01, alpha = .5, position = "identity") +
        geom_vline(
            data = data,
            aes(xintercept = data.mean, colour = "red"),
            linetype = "dashed",
            size = 1
        ) +
        theme_minimal()
    p <- p +
        theme(
            legend.position = "none",
            axis.title = element_text(size = distTextSize),
            axis.text = element_text(size = distTextSize)
        ) +
        labs(
            x = paste0(
                varName,
                " (mean = ",
                round(mean(data[,varType]), 3),
                ")"
            ),
            y = "Frequency"
        )

    return(p)
}
