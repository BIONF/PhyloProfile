#' Identify core genes for a list of selected taxa
#' @description Identify core genes for a list of selected (super)taxa. The
#' identified core genes must be present in at least a certain proportion of
#' species in each selected (super)taxon (identified via percentCutoff) and
#' that criteria must be fullfilled for a certain percentage of selected taxa
#' or all of them (determined via coreCoverage).
#' @export
#' @usage getCoreGene(rankName, taxaCore = c("none"), processedProfileData,
#'     var1Cutoff = c(0, 1), var2Cutoff = c(0, 1), percentCutoff = c(0, 1),
#'     coreCoverage = 1)
#' @param rankName working taxonomy rank (e.g. "species", "genus", "family")
#' @param taxaCore list of selected taxon names
#' @param processedProfileData dataframe contains the full processed
#' phylogenetic profiles (see ?fullProcessedProfile or ?parseInfoProfile)
#' @param var1Cutoff cutoff for var1. Default = c(0, 1).
#' @param var2Cutoff cutoff for var2. Default = c(0, 1).
#' @param percentCutoff cutoff for percentage of species present in each
#' supertaxon. Default = c(0, 1).
#' @param coreCoverage the least percentage of selected taxa should be
#' considered. Default = 1.
#' @return A list of identified core genes.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{parseInfoProfile}} for creating a full processed
#' profile dataframe
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' rankName <- "class"
#' taxaCore <- c("Mammalia", "Saccharomycetes", "Insecta")
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
    rankName = NULL, taxaCore = c("none"), processedProfileData = NULL,
    var1Cutoff = c(0, 1), var2Cutoff = c(0, 1),
    percentCutoff = c(0, 1), coreCoverage = 1
) {
    if (is.null(processedProfileData)) return()
    if (is.null(rankName)) return()

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
    processedProfileData <- processedProfileData[, c(
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
            processedProfileData <- subset(
                processedProfileData,
                supertaxonID %in% superID & mVar1 >= var1CutoffMin
            )
            processedProfileData <- subset(
                processedProfileData,
                supertaxonID %in% superID & mVar1 <= var1CutoffMax
            )
        }
    }

    if (!is.null(var2CutoffMax)) {
        if (!is.na(var2CutoffMax)) {
            processedProfileData <- subset(
                processedProfileData,
                supertaxonID %in% superID & mVar2 >= var2CutoffMin
            )
            processedProfileData <- subset(
                processedProfileData,
                supertaxonID %in% superID & mVar2 <= var2CutoffMax
            )
        }
    }

    # filter by selecting taxa
    if (is.na(superID[1])) return(NULL)
    else {
        data <- subset(
            processedProfileData,
            supertaxonID %in% superID & presSpec >= percentCutoff[1]
        )
        data <- subset(
            data,
            supertaxonID %in% superID & presSpec <= percentCutoff[2]
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
