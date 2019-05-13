# Functions for parsing and pre-processing PhyloProfile input

#' Create a list containing all main taxanomy ranks
#' @export
#' @return A list of all main ranks (from strain to superkingdom)
#' @author Carla Mölbert (carla.moelbert@gmx.de)
#' @examples
#' getTaxonomyRanks()

getTaxonomyRanks <- function(){
    allRanks <- list(
        "Strain " = "strain",
        "Species" = "species",
        "Genus" = "genus",
        "Family" = "family",
        "Order" = "order",
        "Class" = "class",
        "Phylum" = "phylum",
        "Kingdom" = "kingdom",
        "Superkingdom" = "superkingdom",
        "unselected" = ""
    )
    return(allRanks)
}

#' Get list of pre-installed NCBI taxon names
#' @description Get all NCBI taxon names from
#' "PhyloProfile/data/taxonNamesReduced.txt"
#' @export
#' @return List of taxon IDs, their full names, taxonomy ranks and parent IDs
#' obtained from "PhyloProfile/data/taxonNamesReduced.txt"
#' @importFrom utils download.file
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' getNameList()

getNameList <- function() {
    nameReducedFile <- paste(
        system.file(package="PhyloProfile"),
        "PhyloProfile/data/taxonNamesReduced.txt",
        sep="/"
    )

    if (!file.exists(nameReducedFile)) {
        fileURL <- paste0(
            "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/",
            "taxonNamesReduced.txt"
        )
        res <- tryCatch(
            utils::download.file(
                fileURL, destfile = nameReducedFile, method="auto"
            ),
            error=function(e) 1
        )
    }

    nameList <- read.table(
        nameReducedFile, sep = "\t", header = TRUE, fill = TRUE
    )
    nameList$fullName <- as.character(nameList$fullName)
    nameList$rank <- as.character(nameList$rank)
    nameList <- nameList[!duplicated(nameList), ]

    return(nameList)
}

#' Get taxonomy matrix
#' @description Get the (full or subset) taxonomy matrix from
#' "data/taxonomyMatrix.txt" based on an input taxon list
#' @export
#' @param subsetTaxaCheck TRUE/FALSE subset taxonomy matrix based on input taxon
#' IDs. Default = FALSE.
#' @param taxonIDs list of input taxon IDs (e.g. ncbi1234). Default = NULL.
#' @return Data frame contains the (subset of) taxonomy matrix for list of
#' input taxa.
#' @importFrom utils download.file
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' # get full pre-installed taxonomy matrix
#' getTaxonomyMatrix(FALSE, NULL)
#' # get taxonomy matrix for a list of taxon IDs
#' taxonIDs <- c("ncbi10020", "ncbi10090")
#' getTaxonomyMatrix(TRUE, taxonIDs)

getTaxonomyMatrix <- function(subsetTaxaCheck = FALSE, taxonIDs = NULL){
    taxonomyMatrixFile <- paste(
        system.file(package="PhyloProfile"),
        "PhyloProfile/data/taxonomyMatrix.txt",
        sep="/"
    )

    file.exists(taxonomyMatrixFile)
    if (!file.exists(taxonomyMatrixFile)) {
        fileURL <- paste0(
            "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/",
            "taxonomyMatrix.txt"
        )
        res <- tryCatch(
            utils::download.file(
                fileURL, destfile = taxonomyMatrixFile, method="auto"
            ),
            error=function(e) 1
        )
    }

    dt <- read.table(
        taxonomyMatrixFile, sep = "\t", header = TRUE, stringsAsFactors = TRUE
    )
    if (subsetTaxaCheck) {
        if (missing(taxonIDs)) return(dt)
        dt <- dt[dt$abbrName  %in% taxonIDs, ]
    }
    return(dt)
}

#' Get ID list of input taxa from the main input
#' @param rawProfile A dataframe of input phylogenetic profile in long format
#' @return List of all input taxon IDs (e.g. ncbi1234). Default = NULL.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{createLongMatrix}}, \code{\link{mainLongRaw}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' getInputTaxaID(mainLongRaw)

getInputTaxaID <- function(rawProfile = NULL){
    if (is.null(rawProfile)) return()
    inputTaxa <- levels(as.factor(rawProfile$ncbiID))
    inputTaxa <- unlist(strsplit(inputTaxa, split = "\t"))
    # remove "geneID" element from vector inputTaxa
    if (inputTaxa[1] == "geneID") inputTaxa <- inputTaxa[-1]
    # return input taxa
    return(inputTaxa)
}

#' Get NCBI taxon names for a selected list of taxa
#' @description Get NCBI taxon names from
#' "PhyloProfile/data/taxonNamesReduced.txt" for a list of input taxa
#' @param rankName taxonomy rank (e.g. "species","phylum",...)
#' @param taxonIDs list of taxon IDs (e.g. ncbi1234). Default = NULL
#' @return Data frame contains a list of full names, taxonomy ranks and parent
#' IDs for the input taxa.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{getInputTaxaID}} for getting input taxon IDs,
#' \code{\link{getNameList}} for getting the full taxon name list
#' @examples
#' taxonIDs <- c("ncbi10020", "ncbi10090")
#' getInputTaxaName("species", taxonIDs)

getInputTaxaName <- function(rankName, taxonIDs = NULL){
    # check input parameters
    if (missing(rankName)) return("No taxonomy rank name given!")
    allMainRanks <- getTaxonomyRanks()
    if (!(rankName[1] %in% allMainRanks)) return("Invalid taxonomy rank given!")
    # load list of unsorted taxa
    Dt <- getTaxonomyMatrix(TRUE, taxonIDs)
    # load list of taxon name
    nameList <- getNameList()
    # return
    choice <- data.frame(
        "ncbiID" = unlist(Dt[rankName]), stringsAsFactors = FALSE
    )
    choice <- merge(choice, nameList, by = "ncbiID", all = FALSE)
    return(choice)
}

#' Get a subset of input taxa based on a selected taxonomy rank
#' @description Get a subset of taxon ncbi IDs and names from an input list of
#' taxa based on a selected supertaxon (identified by its taxonomy rank and
#' supertaxon name or supertaxon ID).
#' @usage getSelectedTaxonNames(inputTaxonIDs, rank, higherRank, higherID,
#'     higherName)
#' @param inputTaxonIDs list of input taxon IDs (e.g. 876142)
#' @param rank taxonomy rank of input taxa (e.g. "species")
#' @param higherRank selected taxonomy rank (e.g. "phylum")
#' @param higherID supertaxon ID (e.g. 6029). NOTE: either supertaxon ID or
#' name is required, not neccessary to give both.
#' @param higherName supertaxon name (e.g. "Microsporidia"). NOTE: either
#' supertaxon ID or name is required, not neccessary to give both.
#' @export
#' @return A data frame contains ncbi IDs and names of taxa from the input taxon
#' list that belong to the selected supertaxon.
#' @importFrom utils download.file
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' inputTaxonIDs <- c("103449", "1288291", "278021", "27973", "31281", "40302",
#' "586133", "6035", "70536", "876142", "993615")
#' rank <- "species"
#' higherRank <- "genus"
#' higherID <- 6033
#' higherName <- "Encephalitozoon"
#' getSelectedTaxonNames(inputTaxonIDs, rank, higherRank, higherID, higherName)

getSelectedTaxonNames <- function(
    inputTaxonIDs = NULL, rank = NULL,
    higherRank = NULL, higherID = NULL, higherName = NULL
) {
    rankName <- NULL
    if (is.null(inputTaxonIDs) | is.null(rank)) return()
    # load pre-calculated taxonomy daxa
    taxDf <- getTaxonomyMatrix(TRUE, paste0("ncbi", inputTaxonIDs))

    if (is.null(higherID) & is.null(higherName))
        return(
            data.frame(
                ncbiID = taxDf$ncbiID[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                name = taxDf$fullName[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                stringsAsFactors = FALSE
            )
        )

    if (is.null(higherRank)) {
        return(
            data.frame(
                ncbiID = taxDf$ncbiID[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                name = taxDf$fullName[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                stringsAsFactors = FALSE
            )
        )
    } else {
        if (!is.null(higherName) & is.null(higherID)) {
            taxaList <- getNameList()
            superID <- taxaList$ncbiID[
                taxaList$fullName == higherName
                & taxaList$rank %in% c(higherRank, "norank")]
            customizedtaxaID <- levels(
                as.factor(taxDf[rank][taxDf[higherRank] == superID, ])
            )
            return(
                data.frame(
                    ncbiID = taxaList$ncbiID[
                        taxaList$rank %in% c(rank, "norank")
                        & taxaList$ncbiID %in% customizedtaxaID],
                    name = taxaList$fullName[
                        taxaList$rank %in% c(rank, "norank")
                        & taxaList$ncbiID %in% customizedtaxaID],
                    stringsAsFactors = FALSE
                )
            )
        } else if (!is.null(higherID)) {
            return(
                data.frame(
                    ncbiID = taxDf$ncbiID[taxDf[,higherRank] == higherID],
                    name = taxDf$fullName[taxDf[,higherRank] == higherID],
                    stringsAsFactors = FALSE
                )
            )
        }
    }
}

#' Sort list of (super)taxa based on a selected reference (super)taxon
#' @usage sortInputTaxa(taxonIDs = NULL, rankName, refTaxon = NULL,
#'     taxaTree = NULL)
#' @param taxonIDs list of taxon IDs (e.g.: ncbi1234, ncbi9999, ...). Default =
#' NULL.
#' @param rankName working taxonomy rank (e.g. "species", "phylum",...)
#' @param refTaxon selected reference taxon. Default = NULL.
#' @param taxaTree taxonomy tree for the input taxa (optional). Default = NULL.
#' @return A taxonomy matrix for the input taxa ordered by the selected
#' reference taxon. This matrix is sorted either based on the NCBI taxonomy
#' info, or based on an user-defined taxonomy tree (if provided).
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{getNameList}}, \code{\link{getTaxonomyMatrix}},
#' \code{\link{createRootedTree}}, \code{\link{sortTaxaFromTree}},
#' \code{\link{getInputTaxaName}}, \code{\link{getInputTaxaID}},
#' \code{\link{createLongMatrix}}
#' @examples
#' taxonIDs <- c(
#'     "ncbi272557", "ncbi176299", "ncbi3702", "ncbi876142", "ncbi9606"
#' )
#' sortInputTaxa(taxonIDs, "species", "Homo sapiens", NULL)

sortInputTaxa <- function(
    taxonIDs = NULL, rankName, refTaxon = NULL, taxaTree = NULL
){
    ncbiID <- NULL
    fullName <- NULL
    abbrName <- NULL
    # check parameters
    if (missing(rankName)) return("No taxonomy rank name given!")
    allMainRanks <- getTaxonomyRanks()
    if (!(rankName[1] %in% allMainRanks)) return("Invalid taxonomy rank given!")
    if (is.null(refTaxon))  refTaxon <- taxonNames$fullName[1]

    # get list of taxon names
    fullnameList <- getNameList()
    taxonNames <- getInputTaxaName(rankName, taxonIDs)

    # get selected supertaxon ID(s)
    rankNameTMP <- taxonNames$rank[taxonNames$fullName == refTaxon]
    if (rankName == "strain") {
        superID <- fullnameList$ncbiID[
            fullnameList$fullName == refTaxon & fullnameList$rank == "norank"
            ]
    } else {
        superID <- fullnameList$ncbiID[
            fullnameList$fullName == refTaxon
            & fullnameList$rank == rankNameTMP[1]
            ]
    }

    # get full taxonomy data
    Dt <- getTaxonomyMatrix()

    # representative taxon
    repTaxon <- Dt[Dt[, rankName] == superID, ][1, ]

    # THEN, SORT TAXON LIST BASED ON TAXONOMY TREE
    if (is.null(taxaTree)) {
        # prepare Df for calculating distance matrix
        distDf <- subset(Dt, select = -c(ncbiID, fullName))
        row.names(distDf) <- distDf$abbrName
        distDf <- distDf[, -1]
        # create taxonomy tree rooted by refTaxon
        taxaTree <- createRootedTree(distDf, as.character(repTaxon$abbrName))
    } else {
        taxaTree <- ape::root(
            taxaTree,
            outgroup = as.character(repTaxon$abbrName),
            resolve.root = TRUE
        )
    }

    taxonList <- sortTaxaFromTree(taxaTree)
    sortedDt <- Dt[match(taxonList, Dt$abbrName), ]

    # subset to get list of input taxa only
    sortedDt <- subset(sortedDt, abbrName %in% taxonIDs)

    # get only taxonIDs list of selected rank and rename columns
    sortedOut <- subset(
        sortedDt,
        select = c("abbrName", "ncbiID", "fullName", as.character(rankName))
    )
    colnames(sortedOut) <- c("abbrName", "species", "fullName", "ncbiID")

    # add name of supertaxa into sortedOut list
    sortedOut <- merge(
        sortedOut, fullnameList, by = "ncbiID", all.x = TRUE, sort = FALSE
    )
    sortedOut$species <- paste0("ncbi", sortedOut$species)

    ## create new column for sorted supertaxon
    indexSpec <- unlist(lapply(
        seq_len(nlevels(as.factor(sortedOut$fullName.y))),
        function (x) 1000 + x
    ))
    indexSpecDf <- data.frame(
        fullName.y = unique(as.character(sortedOut$fullName.y)),
        sortedSupertaxon = paste0(
            indexSpec, "_", unique(as.character(sortedOut$fullName.y))
        ),
        stringsAsFactors = FALSE
    )
    sortedOut <- plyr::join(indexSpecDf, sortedOut, by = "fullName.y")

    # final sorted supertaxa list
    sortedOut$taxonID <- 0
    sortedOut$category <- "cat"
    sortedOut <- sortedOut[, c(
        "abbrName",
        "taxonID",
        "fullName.x",
        "species",
        "ncbiID",
        "sortedSupertaxon",
        "rank",
        "category"
    )]

    colnames(sortedOut) <- c(
        "abbrName",
        "taxonID",
        "fullName",
        "ncbiID",
        "supertaxonID",
        "supertaxon",
        "rank",
        "category"
    )

    sortedOut$ncbiID <- as.factor(sortedOut$ncbiID)
    sortedOut$supertaxon <- as.factor(sortedOut$supertaxon)
    sortedOut$category <- as.factor(sortedOut$category)

    return(sortedOut)
}

#' Calculate percentage of present species in each super taxon
#' @export
#' @param profileWithTax data frame of main PhyloProfile input together
#' with their taxonomy info (see ?profileWithTaxonomy)
#' @param taxaCount number of species occur in each supertaxon (e.g. phylum
#' or kingdom)
#' @return A data frame with % of present species in each supertaxon
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{profileWithTaxonomy}} for a demo input data
#' @examples
#' # NOTE: for internal testing only - not recommended for outside using
#' data("profileWithTaxonomy", package="PhyloProfile")
#' taxaCount <- plyr::count(profileWithTaxonomy, "supertaxon")
#' taxaCount$freq <- 1
#' calcPresSpec(profileWithTaxonomy, taxaCount)

calcPresSpec <- function(profileWithTax, taxaCount){
    paralog <- NULL
    abbrName <- NULL
    if (missing(profileWithTax)) return ("No input data given")
    if (missing(taxaCount)) return ("No supertaxon count given")

    profileWithTax <- profileWithTax[profileWithTax$orthoID != "NA", ]

    # get geneID and supertaxon
    geneIDSupertaxon <- subset(
        profileWithTax,
        select = c("geneID", "supertaxon", "paralog", "abbrName")
    )
    # remove duplicated rows
    geneIDSupertaxon <- geneIDSupertaxon[!duplicated(geneIDSupertaxon), ]

    # remove NA rows from profileWithTax
    profileWithTaxNoNA <- profileWithTax[profileWithTax$orthoID != "NA", ]

    # count present frequency of supertaxon for each gene
    geneSupertaxonCount <- plyr::count(
        profileWithTaxNoNA, c("geneID", "supertaxon")
    )

    # merge with taxaCount to get total number of species of each supertaxon
    # and calculate presSpec
    presSpecDt <- merge(
        geneSupertaxonCount, taxaCount, by = "supertaxon", all.x = TRUE
    )

    specCount <- plyr::count(geneIDSupertaxon, c("geneID", "supertaxon"))
    presSpecDt <- merge(
        presSpecDt, specCount, by = c("geneID", "supertaxon")
    )

    presSpecDt$presSpec <- presSpecDt$freq / presSpecDt$freq.y

    presSpecDt <- presSpecDt[presSpecDt$presSpec <= 1, ]
    presSpecDt <- presSpecDt[order(presSpecDt$geneID), ]
    presSpecDt <- presSpecDt[, c("geneID", "supertaxon", "presSpec")]

    # add absent supertaxon into presSpecDt
    geneIDSupertaxon <- subset(
        geneIDSupertaxon, select = -c(paralog, abbrName)
    )
    finalPresSpecDt <- merge(
        presSpecDt, geneIDSupertaxon,
        by = c("geneID", "supertaxon"), all.y = TRUE
    )
    finalPresSpecDt$presSpec[is.na(finalPresSpecDt$presSpec)] <- 0

    # remove duplicated and NA rows
    finalPresSpecDt <- finalPresSpecDt[!duplicated(finalPresSpecDt), ]
    finalPresSpecDt <- finalPresSpecDt[complete.cases(finalPresSpecDt), ]

    # return finalPresSpecDt
    return(finalPresSpecDt)
}

#' Parsing info for phylogenetic profiles
#' @description Creating main dataframe for the input phylogenetic profiles
#' based on selected input taxonomy level (e.g. strain, species) and reference
#' taxon. The output contains the number of paralogs, percentage of species
#' presence in each supertaxon, and the max/min/mean/median of VAR1 and VAR2.
#' @usage parseInfoProfile(inputDf, sortedInputTaxa, var1AggregateBy = "max",
#'     var2AggregateBy = "max")
#' @param inputDf input profiles in long format
#' @param sortedInputTaxa sorted taxonomy data for the input taxa
#' (check sortInputTaxa())
#' @param var1AggregateBy aggregate method for VAR1 (max, min, mean
#' or median), applied for calculating var1 of supertaxa. Default = "max".
#' @param var2AggregateBy aggregate method for VAR2 (max, min, mean
#' or median), applied for calculating var2 of supertaxa. Default = "max".
#' @importFrom stats aggregate
#' @return A dataframe contains all info for the input phylogenetic profiles.
#' This full processed profile that is required for several profiling analyses
#' e.g. estimation of gene age (?estimateGeneAge) or identification of core gene
#' (?getCoreGene).
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{createLongMatrix}}, \code{\link{sortInputTaxa}},
#' \code{\link{calcPresSpec}}, \code{\link{mainLongRaw}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' taxonIDs <- getInputTaxaID(mainLongRaw)
#' sortedInputTaxa <- sortInputTaxa(
#'     taxonIDs, "class", "Mammalia", NULL
#' )
#' var1AggregateBy <- "max"
#' var2AggregateBy <- "mean"
#' parseInfoProfile(
#'     mainLongRaw, sortedInputTaxa, var1AggregateBy, var2AggregateBy
#' )

parseInfoProfile <- function(
    inputDf, sortedInputTaxa, var1AggregateBy = "max", var2AggregateBy = "max"
) {
    if (is.null(inputDf) | is.null(sortedInputTaxa)) return()

    # rename columns of 2 additional variables
    if (ncol(inputDf) > 3) {
        if (ncol(inputDf) < 5) colnames(inputDf)[4] <- "var1"
        else colnames(inputDf)[c(4,5)] <- c("var1", "var2")
    }

    # count number of inparalogs
    paralogCount <- plyr::count(inputDf, c("geneID", "ncbiID"))
    inputDf <- merge(inputDf, paralogCount, by = c("geneID", "ncbiID"))
    colnames(inputDf)[ncol(inputDf)] <- "paralog"

    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(sortedInputTaxa, "supertaxon")

    # merge inputDf, mdDataVar2 and sortedInputTaxa to get taxonomy info
    taxaMdData <- merge(inputDf, sortedInputTaxa, by = "ncbiID")

    # (2) calculate PERCENTAGE of PRESENT SPECIES
    finalPresSpecDt <- calcPresSpec(taxaMdData, taxaCount)

    # (3) calculate max/min/mean/median VAR1 for every supertaxon of each gene
    # remove NA rows from taxaMdData
    taxaMdDataNoNA <- taxaMdData[!is.na(taxaMdData$var1), ]
    # calculate m var1
    mVar1Dt <- stats::aggregate(
        taxaMdDataNoNA[, "var1"],
        list(taxaMdDataNoNA$supertaxon, taxaMdDataNoNA$geneID),
        FUN = var1AggregateBy
    )
    colnames(mVar1Dt) <- c("supertaxon", "geneID", "mVar1")

    # (4) calculate max/min/mean/median VAR2 for each super taxon
    # remove NA rows from taxaMdData
    taxaMdDataNoNAVar2 <- taxaMdData[!is.na(taxaMdData$var2), ]
    # calculate max/min/mean/median VAR2
    if (nrow(taxaMdDataNoNAVar2) > 0) {
        mVar2Dt <- aggregate(
            taxaMdDataNoNAVar2[, "var2"],
            list(taxaMdDataNoNAVar2$supertaxon, taxaMdDataNoNAVar2$geneID),
            FUN = var2AggregateBy
        )
        colnames(mVar2Dt) <- c("supertaxon", "geneID", "mVar2")
    } else {
        mVar2Dt <- taxaMdData[, c("supertaxon", "geneID")]
        mVar2Dt$mVar2 <- 0
    }

    # (3+4) & join mVar2 together with mVar1 scores into one df
    scoreDf <- merge(
        mVar1Dt, mVar2Dt, by = c("supertaxon", "geneID"), all = TRUE
    )

    # (2+3+4) add presSpec and mVar1 into taxaMdData
    presMdData <- merge(taxaMdData, finalPresSpecDt,
                        by = c("geneID", "supertaxon"), all.x = TRUE)
    fullMdData <- merge(presMdData, scoreDf,
                        by = c("geneID", "supertaxon"), all.x = TRUE)
    fullMdData <- merge(fullMdData, taxaCount,
                        by = ("supertaxon"), all.x = TRUE)
    # rename "freq" into "numberSpec"
    names(fullMdData)[names(fullMdData) == "freq"] <- "numberSpec"

    fullMdData$fullName <- as.vector(fullMdData$fullName)
    names(fullMdData)[names(fullMdData) == "orthoID.x"] <- "orthoID"

    # parsed input data frame and return
    fullMdData <- fullMdData[!duplicated(fullMdData), ]
    return(fullMdData)
}

#' Reduce the full processed profile data into supertaxon level
#' @description Reduce data of the processed phylogenetic profiles from input
#' taxonomy rank into supertaxon level (e.g. from species to phylum)
#' @param fullProfile dataframe contains the full processed profiles (see
#' ?fullProcessedProfile)
#' @return A reduced dataframe contains only profile data for the selected
#' supertaxon rank. This dataframe contains only supertaxa and their value
#' (%present, mVar1 & mVar2) for each gene.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{parseInfoProfile}} for creating a full processed
#' profile dataframe, \code{\link{fullProcessedProfile}} for a demo full
#' processed profile dataframe
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' reduceProfile(fullProcessedProfile)

reduceProfile <- function(fullProfile) {
    if (is.null(fullProfile)) return()

    # check if working with the lowest taxonomy rank; 1 for NO; 0 for YES
    flag <- 1
    if (length(unique(levels(as.factor(fullProfile$numberSpec)))) == 1) {
        if (unique(levels(as.factor(fullProfile$numberSpec))) == 1) {
            superDfExt <- fullProfile[, c(
                "geneID", "supertaxon", "supertaxonID",
                "var1", "presSpec", "category", "orthoID", "var2", "paralog"
            )]
            flag <- 0
        }
    }
    if (flag == 1) {
        # get representative orthoID that has m VAR1 for each supertaxon
        mOrthoID <- fullProfile[, c(
            "geneID", "supertaxon", "var1", "mVar1", "orthoID"
        )]
        mOrthoID <- subset(mOrthoID, mOrthoID$var1 == mOrthoID$mVar1)
        colnames(mOrthoID) <- c(
            "geneID", "supertaxon", "var1", "mVar1", "orthoID"
        )
        mOrthoID <- mOrthoID[!is.na(mOrthoID$orthoID), ]
        mOrthoID <- mOrthoID[, c("geneID", "supertaxon", "orthoID")]
        mOrthoID <- mOrthoID[!duplicated(mOrthoID[, seq_len(2)]), ]

        # get data set for PhyloProfile plotting (contains only supertaxa info)
        superDf <- subset(fullProfile, select = c(
            "geneID", "supertaxon", "supertaxonID",
            "mVar1", "presSpec", "category", "mVar2", "paralog"
        ))
        superDf$paralog <- 1
        superDf <- superDf[!duplicated(superDf), ]

        superDfExt <- merge(superDf, mOrthoID,
                            by = c("geneID", "supertaxon"), all.x = TRUE)
        superDfExt <- superDfExt[, c(
            "geneID", "supertaxon", "supertaxonID",
            "mVar1", "presSpec", "category", "orthoID", "mVar2", "paralog"
        )]

        # rename mVar to var
        names(superDfExt)[names(superDfExt) == "mVar1"] <- "var1"
        names(superDfExt)[names(superDfExt) == "mVar2"] <- "var2"
    }

    return(superDfExt)
}

#' Filter phylogentic profiles
#' @description Create a filtered data needed for plotting or clustering
#' phylogenetic profiles. NOTE: this function require some intermediate steps
#' using the results from other functions. If you would like to get a full
#' processed data from the raw input, please use the function
#' fromInputToProfile() instead!
#' @usage filterProfileData(superTaxonDf, refTaxon = NULL,
#'     percentCutoff = c(0, 1), coorthologCutoffMax = 9999,
#'     var1Cutoff  = c(0, 1), var2Cutoff = c(0, 1), var1Relation = "protein",
#'     var2Relation = "protein", groupByCat = FALSE, catDt = NULL)
#' @param superTaxonDf a reduced dataframe contains info for all phylogenetic
#' profiles in the selected taxonomy rank.
#' @param refTaxon selected reference taxon. NOTE: This taxon will not be
#' affected by the filtering. If you want to filter all, set refTaxon <- NULL.
#' Default = NULL.
#' @param percentCutoff min and max cutoffs for percentage of species present
#' in a supertaxon. Default = c(0, 1).
#' @param coorthologCutoffMax maximum number of co-orthologs allowed. Default =
#' 9999.
#' @param var1Cutoff min and max cutoffs for var1. Default = c(0, 1).
#' @param var2Cutoff min anc max cutoffs for var2. Default = c(0, 1).
#' @param var1Relation relation of var1 ("protein" for protein-protein or
#' "species" for protein-species). Default = "protein".
#' @param var2Relation relation of var2 ("protein" for protein-protein or
#' "species" for protein-species). Default = "protein".
#' @param groupByCat group genes by their categories (TRUE or FALSE). Default =
#' FALSE.
#' @param catDt dataframe contains gene categories
#' (optional, NULL if groupByCat = FALSE or no info provided). Default = NULL.
#' @return A filtered dataframe for generating profile plot including seed gene
#' IDs (or orthologous group IDs), their ortholog IDs and the corresponding
#' (super)taxa, (super)taxon IDs, number of co-orthologs in each (super)taxon,
#' values for two additional variables var1, var2, % of species present in each
#' supertaxon, and the categories of seed genes (or ortholog groups).
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{parseInfoProfile}} and \code{\link{reduceProfile}}
#' for generating input dataframe, \code{\link{fullProcessedProfile}} for a
#' demo full processed profile dataframe, \code{\link{fromInputToProfile}} for
#' generating fully processed data from raw input.
#' @examples
#' # NOTE: this function require some intermediate steps using the results from
#' # other functions. If you would like to get a full processed data from the
#' # raw input, please use the function fromInputToProfile() instead!
#' data("fullProcessedProfile", package="PhyloProfile")
#' superTaxonDf <- reduceProfile(fullProcessedProfile)
#' refTaxon <- "Mammalia"
#' percentCutoff <- c(0.0, 1.0)
#' coorthologCutoffMax <- 10
#' var1Cutoff <- c(0.75, 1.0)
#' var2Cutoff <- c(0.5, 1.0)
#' var1Relation <- "protein"
#' var2Relation <- "species"
#' groupByCat <- FALSE
#' catDt <- NULL
#' filterProfileData(
#'     superTaxonDf,
#'     refTaxon,
#'     percentCutoff,
#'     coorthologCutoffMax,
#'     var1Cutoff,
#'     var2Cutoff,
#'     var1Relation,
#'     var2Relation,
#'     groupByCat,
#'     catDt
#' )

filterProfileData <- function(
    superTaxonDf,
    refTaxon = NULL,
    percentCutoff = c(0, 1),
    coorthologCutoffMax = 9999,
    var1Cutoff = c(0, 1),
    var2Cutoff = c(0, 1),
    var1Relation = "protein",
    var2Relation = "protein",
    groupByCat = FALSE,
    catDt = NULL
) {
    if (is.null(superTaxonDf)) return()
    if (is.null(refTaxon)) refTaxon = "NA"

    ### remove index from supertaxon names
    superTaxonDf$taxonMod <- gsub("^[[:digit:]]*_", "", superTaxonDf$supertaxon)

    ### replace insufficient values according to the thresholds by NA or 0
    # based on presSpec or # of co-orthologs
    numberCoortholog <- levels(as.factor(superTaxonDf$paralog))
    if (length(numberCoortholog) > 1) {
        superTaxonDf$presSpec[
            superTaxonDf$taxonMod != refTaxon
            & superTaxonDf$paralog > coorthologCutoffMax
            ] <- 0
        if (var2Relation == "protein") {
            superTaxonDf$var2[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$paralog > coorthologCutoffMax
                ] <- NA
        }
    } else {
        if (length(levels(as.factor(superTaxonDf$presSpec))) > 1)
            superTaxonDf$presSpec[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$presSpec < percentCutoff[1]
                ] <- 0
        superTaxonDf$presSpec[
            superTaxonDf$taxonMod != refTaxon
            & superTaxonDf$presSpec > percentCutoff[2]
            ] <- 0
        if (var2Relation == "protein") {
            superTaxonDf$var2[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$presSpec < percentCutoff[1]
                ] <- NA
            superTaxonDf$var2[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$presSpec > percentCutoff[2]
                ] <- NA
        }
    }

    superTaxonDf$presSpec[
        superTaxonDf$taxonMod != refTaxon & superTaxonDf$var1 < var1Cutoff[1]
        ] <- 0
    superTaxonDf$presSpec[
        superTaxonDf$taxonMod != refTaxon & superTaxonDf$var1 > var1Cutoff[2]
        ] <- 0
    if (var1Relation == "protein") {
        if (var2Relation == "protein") {
            # prot-prot: remove complete cell if one variable not sufficient
            superTaxonDf$presSpec[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$var2 < var2Cutoff[1]
                ] <- 0
            superTaxonDf$presSpec[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$var2 > var2Cutoff[2]
                ] <- 0
            superTaxonDf$var2[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$var1 < var1Cutoff[1]
                ] <- NA
            superTaxonDf$var2[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$var1 > var1Cutoff[2]
                ] <- NA
            superTaxonDf$var1[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$var2 < var2Cutoff[1]
                ] <- NA
            superTaxonDf$var1[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$var2 > var2Cutoff[2]
                ] <- NA
        } else {
            # prot-spec: var1 depend on var2
            superTaxonDf$presSpec[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$var2 < var2Cutoff[1]
                ] <- 0
            superTaxonDf$presSpec[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$var2 > var2Cutoff[2]
                ] <- 0
        }
    } else {
        if (var2Relation == "species") {
            # spec-spec: remove var1 and var2 independently
            superTaxonDf$presSpec[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$var1 < var1Cutoff[1]
                ] <- 0
            superTaxonDf$presSpec[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$var1 > var1Cutoff[2]
                ] <- 0
        } else {
            # spec-prot: var2 depend on var1
            superTaxonDf$var2[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$var1 < var1Cutoff[1]
                ] <- NA
            superTaxonDf$var2[
                superTaxonDf$taxonMod != refTaxon
                & superTaxonDf$var1 > var1Cutoff[2]
                ] <- NA
        }
    }

    superTaxonDf$var1[
        superTaxonDf$taxonMod != refTaxon & superTaxonDf$var1 < var1Cutoff[1]
        ] <- NA
    superTaxonDf$var1[
        superTaxonDf$taxonMod != refTaxon & superTaxonDf$var1 > var1Cutoff[2]
        ] <- NA
    superTaxonDf$var2[
        superTaxonDf$taxonMod != refTaxon & superTaxonDf$var2 < var2Cutoff[1]
        ] <- NA
    superTaxonDf$var2[
        superTaxonDf$taxonMod != refTaxon & superTaxonDf$var2 > var2Cutoff[2]
        ] <- NA

    superTaxonDf <- droplevels(superTaxonDf)  # delete unused levels
    superTaxonDf$geneID <- as.factor(superTaxonDf$geneID)
    superTaxonDf$supertaxon <- as.factor(superTaxonDf$supertaxon)

    ### add gene categories (if provided)
    if (groupByCat == TRUE) {
        if (is.null(catDt)) {
            catDt <- data.frame( geneID = levels(superTaxonDf$geneID))
            catDt$group <- "noCategory"
        }

        # create a dataframe that contain all genes and all taxa
        dataHeatCat <- data.frame(
            supertaxon = rep(
                levels(superTaxonDf$supertaxon), nlevels(superTaxonDf$geneID)
            ),
            geneID = rep(
                levels(superTaxonDf$geneID),
                each = nlevels(superTaxonDf$supertaxon)
            )
        )

        dataHeatCat <- merge(dataHeatCat, catDt, by = "geneID")

        # add categories into superTaxonDf
        superTaxonDf <- merge(
            dataHeatCat, superTaxonDf,
            by = c("geneID","supertaxon"), all.x = TRUE
        )
    }

    return(superTaxonDf)
}

#' Complete processing of raw input phylogenetic profiles
#' @description Create a processed and filtered data for plotting or analysing
#' phylogenetic profiles from raw input file (from raw input to final filtered
#' dataframe)
#' @usage fromInputToProfile(rawInput, rankName, refTaxon = NULL,
#'     taxaTree = NULL, var1AggregateBy = "max", var2AggregateBy = "max",
#'     percentCutoff = c(0, 1), coorthologCutoffMax = 9999,
#'     var1Cutoff = c(0, 1), var2Cutoff = c(0, 1), var1Relation = "protein",
#'     var2Relation = "protein", groupByCat = FALSE, catDt = NULL)
#' @param rawInput input file (in long, wide, multi-fasta or orthoxml format)
#' @param rankName taxonomy rank (e.g. "species","phylum",...)
#' @param refTaxon selected reference taxon name (used for sorting and will be
#' protected from filtering). Default = NULL.
#' @param taxaTree input taxonomy tree for taxa in input profiles (optional).
#' Default = NULL.
#' @param var1AggregateBy aggregate method for var1 (min, max, mean or median).
#' Default = "max".
#' @param var2AggregateBy aggregate method for VAR2 (min, max, mean or median).
#' Default = "max".
#' @param percentCutoff min and max cutoffs for percentage of species present
#' in a supertaxon. Default = c(0, 1).
#' @param coorthologCutoffMax maximum number of co-orthologs allowed. Default =
#' 9999.
#' @param var1Cutoff min and max cutoffs for var1. Default = c(0, 1).
#' @param var2Cutoff min and max cutoffs for var2. Default = c(0, 1).
#' @param var1Relation relation of var1 ("protein" for protein-protein or
#' "species" for protein-species). Default = "protein".
#' @param var2Relation relation of var2 ("protein" for protein-protein or
#' "species" for protein-species). Default = "protein".
#' @param groupByCat group genes by their categories (TRUE or FALSE). Default =
#' FALSE.
#' @param catDt dataframe contains gene categories. Default = NULL.
#' @return Dataframe required for generating phylogenetic profile plot or
#' clustering analysis. It contains seed gene IDs (or orthologous group IDs),
#' their ortholog IDs and the corresponding (super)taxa, (super)taxon IDs,
#' number of co-orthologs in each (super)taxon, values for two additional
#' variables var1, var2, % of species present in each supertaxon, and the
#' categories of seed genes (or ortholog groups).
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{createLongMatrix}}, \code{\link{getInputTaxaID}},
#' \code{\link{getInputTaxaName}}, \code{\link{sortInputTaxa}},
#' \code{\link{parseInfoProfile}}, \code{\link{reduceProfile}},
#' \code{\link{filterProfileData}}
#' @examples
#' rawInput <- system.file(
#'     "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' rankName <- "class"
#' refTaxon <- "Mammalia"
#' taxaTree <- NULL
#' var1AggregateBy <- "max"
#' var2AggregateBy <- "mean"
#' percentCutoff <- c(0.0, 1.0)
#' coorthologCutoffMax <- 10
#' var1Cutoff <- c(0.75, 1.0)
#' var2Cutoff <- c(0.5, 1.0)
#' var1Relation <- "protein"
#' var2Relation <- "species"
#' groupByCat <- FALSE
#' catDt <- NULL
#' fromInputToProfile(
#'     rawInput,
#'     rankName,
#'     refTaxon,
#'     taxaTree,
#'     var1AggregateBy,
#'     var2AggregateBy,
#'     percentCutoff,
#'     coorthologCutoffMax,
#'     var1Cutoff,
#'     var2Cutoff,
#'     var1Relation,
#'     var2Relation,
#'     groupByCat,
#'     catDt
#' )

fromInputToProfile <- function(
    rawInput,
    rankName,
    refTaxon = NULL,
    taxaTree = NULL,
    var1AggregateBy = "max",
    var2AggregateBy = "max",
    percentCutoff = c(0, 1),
    coorthologCutoffMax = 9999,
    var1Cutoff = c(0, 1),
    var2Cutoff = c(0, 1),
    var1Relation = "protein",
    var2Relation = "protein",
    groupByCat = FALSE,
    catDt = NULL
) {
    if (missing(rawInput) | missing(rankName)) return("Missing input")
    if (is.null(rawInput) | is.null(rankName)) return("Missing input")
    # convert raw input into long format
    inputDf <- createLongMatrix(rawInput)

    # get input taxon IDs and names
    inputTaxonID <- getInputTaxaID(inputDf)

    # sort input taxa based on selected reference taxon or input taxonomy tree
    sortedInputTaxa <- sortInputTaxa(inputTaxonID, rankName, refTaxon, taxaTree)

    # parse info (additional values...) into profile df
    fullMdData <- parseInfoProfile(
        inputDf, sortedInputTaxa, var1AggregateBy, var2AggregateBy
    )

    # reduce profile df into supertaxon level
    superTaxonDf <- reduceProfile(fullMdData)

    # create final df
    dataHeat <- filterProfileData(
        superTaxonDf,
        refTaxon,
        percentCutoff,
        coorthologCutoffMax,
        var1Cutoff,
        var2Cutoff,
        var1Relation,
        var2Relation,
        groupByCat = FALSE,
        catDt = NULL
    )

    return(dataHeat)
}
