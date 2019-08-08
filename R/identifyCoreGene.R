#' Identify core genes for a list of selected taxa
#' @description Identify core genes for a list of selected (super)taxa. The
#' identified core genes must be present in at least a certain proportion of
#' species in each selected (super)taxon (identified via percentCutoff) and
#' that criteria must be fullfilled for a certain percentage of selected taxa
#' or all of them (determined via coreCoverage).
#' @export
#' @usage getCoreGene(rankName, taxaCore = c("none"), profileDt,
#'     var1Cutoff = c(0, 1), var2Cutoff = c(0, 1), percentCutoff = c(0, 1),
#'     coreCoverage = 1)
#' @param rankName working taxonomy rank (e.g. "species", "genus", "family")
#' @param taxaCore list of selected taxon names
#' @param profileDt dataframe contains the full processed
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
#' profileDt <- fullProcessedProfile
#' var1Cutoff <- c(0.75, 1.0)
#' var2Cutoff <- c(0.75, 1.0)
#' percentCutoff <- c(0.0, 1.0)
#' coreCoverage <- 1
#' getCoreGene(
#'     rankName,
#'     taxaCore,
#'     profileDt,
#'     var1Cutoff, var2Cutoff,
#'     percentCutoff, coreCoverage
#' )

getCoreGene <- function(
    rankName = NULL, taxaCore = c("none"), profileDt = NULL, 
    var1Cutoff = c(0, 1), var2Cutoff = c(0, 1),
    percentCutoff = c(0, 1), coreCoverage = 1
) {
    if (is.null(profileDt)) stop("Processed profile cannot be NULL!")
    if (is.null(rankName)) stop("Rank name cannot be NULL!")
    supertaxonID <- mVar1 <- mVar2 <- presSpec <- Freq <- NULL
    # get ID list of chosen taxa & main input profile
    taxaList <- getNameList()
    if ("none" %in% taxaCore) {
        superID <- NA
    } else superID <- taxaList$ncbiID[
        taxaList$fullName%in%taxaCore & taxaList$rank %in% c(rankName,"norank")]
    profileDt <- profileDt[, c(
        "geneID", "ncbiID", "fullName", "supertaxon", "supertaxonID", "rank",
        "presSpec", "mVar1", "mVar2")]
    # filter by var1 and var2 cutoffs
    if (!is.null(var1Cutoff[2])) {
        if (!is.na(var1Cutoff[2])) {
            profileDt <- subset(
                profileDt, supertaxonID %in% superID & mVar1 >= var1Cutoff[1]
                & mVar1 <= var1Cutoff[2])
        }
    }
    if (!is.null(var2Cutoff[2])) {
        if (!is.na(var2Cutoff[2])) {
            profileDt <- subset(
                profileDt, supertaxonID %in% superID & mVar2 >= var2Cutoff[1]
                & mVar2 <= var2Cutoff[2])
        }
    }
    # filter by selecting taxa
    if (is.na(superID[1])) stop("No core gene found!")
    else {
        data <- subset(
            profileDt, supertaxonID %in% superID & presSpec >= percentCutoff[1]
            & presSpec <= percentCutoff[2])
        # get supertaxa present in each geneID
        supertaxonCount <- as.data.frame(
            plyr::count(data, c("geneID", "supertaxonID")))
        # count no. supertaxa for each gene & min no. supertaxa muss be present
        count <- as.data.frame(table(supertaxonCount$geneID))
        requireCoverage <- length(superID) * (coreCoverage / 100)
        # get only gene that contains orthologs in that coverage # of taxa
        coreGene <- subset(count, Freq >= requireCoverage)
        return(levels(factor(coreGene$Var1)))
    }
}
