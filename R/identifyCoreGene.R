#' Identify core genes for a list of selected taxa
#' @export
#' @usage getCoreGene(rankName, taxaCore, processedProfileData,
#'     var1Cutoff, var2Cutoff, percentCutoff, coreCoverage)
#' @param rankName taxonomy rank (e.g. "species", "genus", "family")
#' @param taxaCore name list of selected taxa
#' @param processedProfileData dataframe contains the full processed
#' phylogenetic profiles
#' @param var1Cutoff cutoff for var1
#' @param var2Cutoff cutoff for var2
#' @param percentCutoff cutoff for percentage of species present in each
#' supertaxon
#' @param coreCoverage the least number of selected taxa should be considered
#' @return A list of core genes
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{parseInfoProfile}} for creating a full processed
#' profile dataframe
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' rankName <- "class"
#' taxaCore <- c("Mammalia", "Mucorales", "Alphaproteobacteria")
#' processedProfileData <- fullProcessedProfile
#' var1Cutoff <- c(0.75, 1.0)
#' var2Cutoff <- c(0.75, 1.0)
#' percentCutoff <- c(0.0, 1.0)
#' coreCoverage <- 1
#' getCoreGene(
#'     rankName,
#'     taxaCore,
#'     processedProfileData,
#'     var1Cutoff, var2Cutoff,
#'     percentCutoff, coreCoverage
#' )

getCoreGene <- function(
    rankName,
    taxaCore,
    processedProfileData,
    var1Cutoff, var2Cutoff,
    percentCutoff, coreCoverage
) {
    supertaxonID <- NULL
    mVar1 <- NULL
    mVar2 <- NULL
    presSpec <- NULL
    Freq <- NULL

    # get ID list of chosen taxa
    taxaList <- getNameList()

    if ("none" %in% taxaCore) {
        superID <- NA
    } else {
        superID <- taxaList$ncbiID[
            taxaList$fullName %in% taxaCore
            & taxaList$rank %in% c(rankName, "norank")
        ]
    }

    # get main input data
    mdData <- processedProfileData
    if (is.null(mdData)) return()
    mdData <- mdData[, c(
        "geneID",
        "ncbiID",
        "fullName",
        "supertaxon",
        "supertaxonID",
        "rank",
        "presSpec",
        "mVar1",
        "mVar2"
    )]

    # filter by var1 and var2 cutoffs
    var1CutoffMin <- var1Cutoff[1]
    var1CutoffMax <- var1Cutoff[2]
    var2CutoffMin <- var2Cutoff[1]
    var2CutoffMax <- var2Cutoff[2]

    if (!is.null(var1CutoffMax)) {
        if (!is.na(var1CutoffMax)) {
            mdData <- subset(
                mdData, supertaxonID %in% superID
                & mVar1 >= var1CutoffMin
            )
            mdData <- subset(
                mdData, supertaxonID %in% superID
                & mVar1 <= var1CutoffMax
            )
        }
    }

    if (!is.null(var2CutoffMax)) {
        if (!is.na(var2CutoffMax)) {
            mdData <- subset(
                mdData, supertaxonID %in% superID
                & mVar2 >= var2CutoffMin
            )
            mdData <- subset(
                mdData, supertaxonID %in% superID
                & mVar2 <= var2CutoffMax
            )
        }
    }

    # filter by selecting taxa
    if (is.na(superID[1])) return(NULL) #mdData <- NULL
    else {
        data <- subset(
            mdData, supertaxonID %in% superID
            & presSpec >= percentCutoff[1]
        )
        data <- subset(
            data, supertaxonID %in% superID
            & presSpec <= percentCutoff[2]
        )

        # get supertaxa present in each geneID
        supertaxonCount <- as.data.frame(
            plyr::count(data, c("geneID", "supertaxonID"))
        )

        # count number of supertaxa present in each geneID
        # and get min number of supertaxa muss be taken into account
        count <- as.data.frame(table(supertaxonCount$geneID))
        requireCoverage <- length(superID) * (coreCoverage / 100)

        # get only gene that contains orthologs in that coverage # of taxa
        coreGene <- subset(count, Freq >= requireCoverage)
        coreGene$Var1 <- factor(coreGene$Var1)

        return(levels(coreGene$Var1))
    }
}
