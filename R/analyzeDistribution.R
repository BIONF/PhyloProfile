#' Create data for percentage present taxa distribution
#' @usage createPercentageDistributionData(inputData = NULL, rankName = NULL, 
#'     taxDB = NULL)
#' @param inputData dataframe contains raw input data in long format
#' (see ?mainLongRaw)
#' @param rankName name of the working taxonomy rank (e.g. "species", "family")
#' @param taxDB Path to the taxonomy DB files
#' @return A dataframe for analysing the distribution of the percentage of
#' species in the selected supertaxa, containing the seed protein IDs,
#' percentage of their orthologs in each supertaxon and the corresponding
#' supertaxon names.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{mainLongRaw}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' createPercentageDistributionData(mainLongRaw, "class")

createPercentageDistributionData <- function(
    inputData = NULL, rankName = NULL, taxDB = NULL
) {
    if (is.null(inputData) | is.null(rankName))
        stop("Input data or rank name cannot be NULL!")
    allMainRanks <- getTaxonomyRanks()
    if (!(rankName[1] %in% allMainRanks)) stop("Invalid taxonomy rank given!")
    if (ncol(inputData) < 4) {
        colnames(inputData) <- c("geneID", "ncbiID", "orthoID")
    } else if (ncol(inputData) < 5) {
        colnames(inputData) <- c("geneID", "ncbiID", "orthoID", "var1")
    } else {
        colnames(inputData) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
    }
    # count number of inparalogs
    paralogCount <- plyr::count(inputData, c("geneID", "ncbiID"))
    inputData <- merge(inputData, paralogCount, by = c("geneID", "ncbiID"))
    colnames(inputData)[ncol(inputData)] <- "paralog"
    # get sorted taxonomy list
    inputTaxonID <- getInputTaxaID(inputData)
    inputTaxonName <- getInputTaxaName(rankName, inputTaxonID, taxDB)
    refTaxon <- inputTaxonName$fullName[1]
    taxaTree <- NULL
    taxaList <- sortInputTaxa(inputTaxonID, rankName, refTaxon, taxaTree, taxDB)
    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxaList, "supertaxon")
    # merge inputData, inputDatavar2 and taxaList to get taxonomy info
    taxaMdData <- merge(inputData, taxaList, by = "ncbiID")
    # calculate % present species
    finalPresSpecDt <- calcPresSpec(taxaMdData, taxaCount)
    finalPresSpecDt[!is.na(finalPresSpecDt$geneID),]
    return(finalPresSpecDt)
}

#' Create data for additional variable distribution
#' @usage createVariableDistributionData(inputData, var1Cutoff = c(0 ,1),
#'     var2Cutoff = c(0, 1))
#' @param inputData dataframe contains raw input data in long format
#' (see ?mainLongRaw)
#' @param var1Cutoff min and max cutoff for var1. Default = c(0, 1).
#' @param var2Cutoff min and max cutoff for var2. Default = c(0, 1).
#' @return A dataframe for analysing the distribution of the additional
#' variable(s) containing the protein (ortholog) IDs and the values of their
#' variables (var1 and var2).
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{mainLongRaw}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' createVariableDistributionData(
#'     mainLongRaw, c(0, 1), c(0.5, 1)
#' )

createVariableDistributionData <- function(
    inputData, var1Cutoff = c(0, 1), var2Cutoff = c(0, 1)
) {
    if (ncol(inputData) < 4) {
        colnames(inputData) <- c("geneID", "ncbiID", "orthoID")
        splitDt <- inputData[, c("orthoID")]
    } else if (ncol(inputData) < 5) {
        colnames(inputData) <- c("geneID", "ncbiID", "orthoID", "var1")
        splitDt <- inputData[, c("orthoID", "var1")]
    } else {
        colnames(inputData) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
        splitDt <- inputData[, c("orthoID", "var1", "var2")]
    }
    splitDt$orthoID[splitDt$orthoID == "NA" | is.na(splitDt$orthoID)] <- NA
    splitDt <- splitDt[stats::complete.cases(splitDt), ]

    if (length(levels(as.factor(splitDt$var2))) == 1)
        if (levels(as.factor(splitDt$var2)) == "") splitDt$var2 <- 0

    # Filter based on variable cutoffs
    if ("var1" %in% colnames(splitDt))
        splitDt <- splitDt[
            splitDt$var1 >= var1Cutoff[1] & splitDt$var1 <= var1Cutoff[2],]
    if ("var2" %in% colnames(splitDt))
        splitDt <- splitDt[
            splitDt$var2 >= var2Cutoff[1] & splitDt$var2 <= var2Cutoff[2],]
    return(splitDt)
}

#' Create data for additional variable distribution (for a subset data)
#' @usage createVariableDistributionDataSubset(fullProfileData,
#'     distributionData, selectedGenes, selectedTaxa)
#' @param fullProfileData dataframe contains the full processed profiles (see
#' ?fullProcessedProfile, ?filterProfileData or ?fromInputToProfile)
#' @param distributionData dataframe contains the full distribution data (see
#' ?createVariableDistributionData)
#' @param selectedGenes list of genes of interest. Default = "all".
#' @param selectedTaxa list of taxa of interest Default = "all".
#' @return A dataframe for analysing the distribution of the additional
#' variable(s) for a subset of genes and/or taxa containing the protein
#' (ortholog) IDs and the values of their variables (var1 and var2).
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{parseInfoProfile}},
#' \code{\link{createVariableDistributionData}},
#' \code{\link{fullProcessedProfile}}, \code{\link{mainLongRaw}}
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' data("mainLongRaw", package="PhyloProfile")
#' distributionData <- createVariableDistributionData(
#'     mainLongRaw, c(0, 1), c(0.5, 1)
#' )
#' selectedGenes <- "100136at6656"
#' selectedTaxa <- c("Mammalia", "Saccharomycetes", "Insecta")
#' createVariableDistributionDataSubset(
#'     fullProcessedProfile,
#'     distributionData,
#'     selectedGenes,
#'     selectedTaxa
#' )

createVariableDistributionDataSubset <- function(
    fullProfileData, distributionData,
    selectedGenes = "all", selectedTaxa = "all"
) {
    geneID <- orthoID <- var1.x <- var2.y <- supertaxonMod <- NULL
    # check parameters
    if (is.null(fullProfileData) | is.null(distributionData))
        stop("Full processed profiles or distribution data cannot be NULL!")
    # get geneID and supertaxon name for distributionData
    distributionDataName <- merge(
        distributionData, fullProfileData, by = "orthoID", all.x = TRUE
    )
    distributionDataName$supertaxonMod <- substr(
        distributionDataName$supertaxon, 6,
        nchar(as.character(distributionDataName$supertaxon))
    )
    distributionDataName <- subset(
        distributionDataName,
        select = c(orthoID, var1.x, var2.y, supertaxonMod, geneID)
    )
    colnames(distributionDataName) <- c(
        "orthoID", "var1", "var2", "supertaxonMod", "geneID"
    )
    # filter data
    if (selectedTaxa[1] == "all" & selectedGenes[1] != "all") {
        # select data from dataHeat for selected sequences only
        distributionData <- subset(
            distributionDataName, geneID %in% selectedGenes
        )
    } else if (selectedGenes[1] == "all" & selectedTaxa[1] != "all") {
        # select data from dataHeat for selected taxa only
        distributionData <- subset(
            distributionDataName, supertaxonMod %in% selectedTaxa
        )
    } else if (selectedGenes[1] == "all" & selectedTaxa[1] == "all") {
        return(distributionData)
    } else {
        # select data from dataHeat for selected sequences and taxa
        distributionData <- subset(
            distributionDataName,
            geneID %in% selectedGenes & supertaxonMod %in% selectedTaxa
        )
    }
    return(distributionData)
}

#' Create distribution plot
#' @description Create distribution plot for one of the additional variable or
#' the percentage of the species present in the supertaxa.
#' @usage createVarDistPlot(data, varName = "var", varType = "var1",
#'     percent = c(0, 1), textSize = 12)
#' @param data dataframe contains data for plotting (see
#' ?createVariableDistributionData, ?createVariableDistributionDataSubset or
#' ?createPercentageDistributionData)
#' @param varName name of the variable that need to be analyzed (either name of
#' variable 1 or variable 2 or "percentage of present taxa"). Default = "var".
#' @param varType type of variable (either "var1", "var2" or "presSpec").
#' Default = "var1".
#' @param percent range of percentage cutoff (between 0 and 1). Default = c(0,1)
#' @param textSize text size of the distribution plot (in px). Default = 12.
#' @return A distribution plot for the selected variable as a ggplot object
#' @import ggplot2
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{mainLongRaw}},
#' \code{\link{createVariableDistributionData}},
#' \code{\link{createVariableDistributionDataSubset}},
#' \code{\link{createPercentageDistributionData}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' data <- createVariableDistributionData(
#'     mainLongRaw, c(0, 1), c(0.5, 1)
#' )
#' varName <- "Variable abc"
#' varType <- "var1"
#' percent <- c(0,1)
#' textSize <- 12
#' createVarDistPlot(
#'     data,
#'     varName,
#'     varType,
#'     percent,
#'     textSize
#' )

createVarDistPlot <- function(
    data, varName = "var", varType = "var1", percent = c(0, 1), textSize = 12
) {
    if (is.null(data)) stop("Input data cannot be NULL!")
    if (varType == "presSpec") {
        # remove presSpec < cutoffMin or > cutoffMax
        if (percent[1] > 0) {
            data <- data[
                data$presSpec >= percent[1] & data$presSpec <= percent[2],
            ]
        } else {
            if (percent[2] > 0) {
                data <- data[data$presSpec > 0 & data$presSpec <= percent[2], ]
            } else {
                data <- data[data$presSpec > 0, ]
            }
        }
    } else data <- data[!is.na(data[,varType]), ]

    data.mean <- mean(data[,varType])
    p <- ggplot(data, aes(x = data[,varType])) +
        geom_histogram(binwidth = .01, alpha = .5, position = "identity") +
        geom_vline(
            data = data, aes(xintercept = data.mean, colour = "red"),
            linetype = "dashed", linewidth = 1
        ) + theme_minimal()
    p <- p + theme(
            legend.position = "none",
            axis.title = element_text(size = textSize),
            axis.text = element_text(size = textSize)
        ) + labs(
            x = paste0(varName," (mean = ", round(mean(data[,varType]), 3),")"),
            y = "Frequency"
        )
    return(p)
}
