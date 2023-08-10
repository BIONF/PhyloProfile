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
#' @param taxDB Path to the taxonomy DB files
#' @return List of taxon IDs, their full names, taxonomy ranks and parent IDs
#' obtained from "PhyloProfile/data/taxonNamesReduced.txt"
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' getNameList()

getNameList <- function(taxDB = NULL) {
    if (is.null(taxDB)) {
        nameReducedFile <- paste(
            system.file(package = "PhyloProfile"),
            "PhyloProfile/data/taxonNamesReduced.txt", sep = "/"
        )
    } else nameReducedFile <- paste(taxDB, "taxonNamesReduced.txt", sep = "/")

    if (!file.exists(nameReducedFile)) {
        utils::data(taxonNamesReduced)
    } else {
        taxonNamesReduced <- utils::read.table(
            nameReducedFile, sep = "\t", header = TRUE, fill = TRUE,
            comment.char = ""
        )
    }

    taxonNamesReduced$fullName <- as.character(taxonNamesReduced$fullName)
    taxonNamesReduced$rank <- as.character(taxonNamesReduced$rank)
    taxonNamesReduced <- taxonNamesReduced[!duplicated(taxonNamesReduced), ]

    return(taxonNamesReduced)
}

#' Get taxonomy matrix
#' @description Get the (full or subset) taxonomy matrix from
#' "data/taxonomyMatrix.txt" based on an input taxon list
#' @export
#' @param taxDB Path to the taxonomy DB files
#' @param subsetTaxaCheck TRUE/FALSE subset taxonomy matrix based on input taxon
#' IDs. Default = FALSE
#' @param taxonIDs list of input taxon IDs (e.g. ncbi1234). Default = NULL
#' @return Data frame contains the (subset of) taxonomy matrix for list of
#' input taxa.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' # get full pre-installed taxonomy matrix
#' getTaxonomyMatrix()
#' # get taxonomy matrix for a list of taxon IDs
#' taxonIDs <- c("ncbi9606", "ncbi10116")
#' getTaxonomyMatrix(NULL, TRUE, taxonIDs)

getTaxonomyMatrix <- function(
        taxDB = NULL, subsetTaxaCheck = FALSE, taxonIDs = NULL
){
    if (is.null(taxDB)) {
        taxonomyMatrixFile <- paste(
            system.file(package = "PhyloProfile"),
            "PhyloProfile/data/taxonomyMatrix.txt", sep = "/"
        )
    } else taxonomyMatrixFile <- paste(taxDB, "taxonomyMatrix.txt", sep = "/")


    if (!file.exists(taxonomyMatrixFile)) {
        utils::data(taxonomyMatrix)
    } else {
        taxonomyMatrix <- utils::read.table(
            taxonomyMatrixFile, sep = "\t", header = TRUE,
            stringsAsFactors = TRUE
        )
    }

    if (subsetTaxaCheck) {
        if (missing(taxonIDs)) return(taxonomyMatrix)
        taxonomyMatrix <- taxonomyMatrix[
            taxonomyMatrix$abbrName  %in% taxonIDs, ]
    }
    return(taxonomyMatrix)
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
    if (is.null(rawProfile)) stop("Input profile data cannot be NULL!")
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
#' @param taxDB Path to the taxonomy DB files
#' @return Data frame contains a list of full names, taxonomy ranks and parent
#' IDs for the input taxa.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{getInputTaxaID}} for getting input taxon IDs,
#' \code{\link{getNameList}} for getting the full taxon name list
#' @examples
#' taxonIDs <- c("ncbi9606", "ncbi10116")
#' getInputTaxaName("species", taxonIDs)

getInputTaxaName <- function(rankName, taxonIDs = NULL, taxDB = NULL){
    # check input parameters
    if (missing(rankName)) return("No taxonomy rank name given!")
    allMainRanks <- getTaxonomyRanks()
    if (!(rankName[1] %in% allMainRanks)) return("Invalid taxonomy rank given!")
    # load list of unsorted taxa
    Dt <- getTaxonomyMatrix(taxDB, TRUE, taxonIDs)
    # load list of taxon name
    nameList <- getNameList(taxDB)
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
#' @usage getSelectedTaxonNames(inputTaxonIDs = NULL, rank = NULL,
#'     higherRank = NULL, higherID = NULL, higherName = NULL, taxDB = NULL)
#' @param inputTaxonIDs list of input taxon IDs (e.g. c("10116", "122586"))
#' @param rank taxonomy rank of input taxa (e.g. "species")
#' @param higherRank selected taxonomy rank (e.g. "phylum")
#' @param higherID supertaxon ID (e.g. 7711). NOTE: either supertaxon ID or
#' name is required, not neccessary to give both
#' @param higherName supertaxon name (e.g. "Chordata"). NOTE: either
#' supertaxon ID or name is required, not neccessary to give both
#' @param taxDB Path to the taxonomy DB files
#' @export
#' @return A data frame contains ncbi IDs and names of taxa from the input taxon
#' list that belong to the selected supertaxon.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' inputTaxonIDs <- c("10116", "122586", "123851", "13616", "188937", "189518",
#' "208964", "224129", "224324", "237631", "243230")
#' rank <- "species"
#' higherRank <- "phylum"
#' higherID <- 7711
#' getSelectedTaxonNames(inputTaxonIDs, rank, higherRank, higherID, NULL)
#' higherName <- "Chordata"
#' getSelectedTaxonNames(inputTaxonIDs, rank, higherRank, NULL, higherName,NULL)

getSelectedTaxonNames <- function(
    inputTaxonIDs = NULL, rank = NULL,
    higherRank = NULL, higherID = NULL, higherName = NULL, taxDB = NULL
) {
    rankName <- NULL
    if (is.null(inputTaxonIDs) | is.null(rank))
        stop("Input taxa and taxonomy rank cannot be NULL!")
    taxDf <- getTaxonomyMatrix(taxDB, TRUE, paste0("ncbi", inputTaxonIDs))
    if (is.null(higherID) & is.null(higherName))
        return(
            data.frame(
                ncbiID = taxDf$ncbiID[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                name = taxDf$fullName[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                stringsAsFactors = FALSE))
    if (is.null(higherRank)) {
        return(
            data.frame(
                ncbiID = taxDf$ncbiID[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                name = taxDf$fullName[
                    taxDf$abbrName %in% paste0("ncbi", inputTaxonIDs)],
                stringsAsFactors = FALSE))
    } else {
        if (!is.null(higherName) & is.null(higherID)) {
            taxaList <- getNameList(taxDB)
            superID <- taxaList$ncbiID[
                taxaList$fullName == higherName
                & taxaList$rank %in% c(higherRank, "norank")]
            customizedtaxaID <- levels(
                as.factor(taxDf[rank][taxDf[higherRank] == superID, ]))
            return(
                data.frame(
                    ncbiID = taxaList$ncbiID[
                        taxaList$rank %in% c(rank, "norank")
                        & taxaList$ncbiID %in% customizedtaxaID],
                    name = taxaList$fullName[
                        taxaList$rank %in% c(rank, "norank")
                        & taxaList$ncbiID %in% customizedtaxaID],
                    stringsAsFactors = FALSE))
        } else if (!is.null(higherID)) {
            return(
                data.frame(
                    ncbiID = taxDf$ncbiID[taxDf[,higherRank] == higherID],
                    name = taxDf$fullName[taxDf[,higherRank] == higherID],
                    stringsAsFactors = FALSE))
        }
    }
}

#' Sort list of (super)taxa based on a selected reference (super)taxon
#' @usage sortInputTaxa(taxonIDs = NULL, rankName, refTaxon = NULL,
#'     taxaTree = NULL, sortedTaxonList = NULL, taxDB = NULL)
#' @param taxonIDs list of taxon IDs (e.g.: ncbi1234, ncbi9999, ...). Default =
#' NULL
#' @param rankName working taxonomy rank (e.g. "species", "phylum",...)
#' @param refTaxon selected reference taxon. Default = NULL
#' @param taxaTree taxonomy tree for the input taxa (optional). Default = NULL
#' @param sortedTaxonList list of sorted taxa (optional). Default = NULL
#' @param taxDB Path to the taxonomy DB files
#' @return A taxonomy matrix for the input taxa ordered by the selected
#' reference taxon. This matrix is sorted either based on the NCBI taxonomy
#' info, or based on an user-defined taxonomy tree (if provided).
#' @importFrom ape read.tree
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{getNameList}}, \code{\link{getTaxonomyMatrix}},
#' \code{\link{createUnrootedTree}}, \code{\link{sortTaxaFromTree}},
#' \code{\link{getInputTaxaName}}, \code{\link{getInputTaxaID}},
#' \code{\link{createLongMatrix}}
#' @examples
#' taxonIDs <- c(
#'     "ncbi10116", "ncbi123851", "ncbi3702", "ncbi13616", "ncbi9606"
#' )
#' sortInputTaxa(taxonIDs, "species", "Homo sapiens", NULL, NULL)

sortInputTaxa <- function(
    taxonIDs = NULL, rankName, refTaxon = NULL, taxaTree = NULL,
    sortedTaxonList = NULL, taxDB = NULL
){
    ncbiID <- fullName <- abbrName <- NULL
    if (missing(rankName)) return("No taxonomy rank name given!")
    allMainRanks <- getTaxonomyRanks()
    if (!(rankName[1] %in% allMainRanks)) return("Invalid taxonomy rank given!")
    # get list of taxon names
    fullnameList <- getNameList(taxDB)
    taxonNames <- getInputTaxaName(rankName, taxonIDs, taxDB)
    if (is.null(refTaxon))  refTaxon <- taxonNames$fullName[1]
    # get selected supertaxon ID(s)
    rankNameTMP <- taxonNames$rank[taxonNames$fullName == refTaxon]
    if (rankName == "strain") {
        superID <- fullnameList$ncbiID[fullnameList$fullName == refTaxon]
    } else
        superID <- fullnameList$ncbiID[
            fullnameList$fullName == refTaxon
            & fullnameList$rank == rankNameTMP[1]]
    # get full taxonomy data & representative taxon
    Dt <- getTaxonomyMatrix(taxDB)
    repTaxon <- Dt[Dt[, rankName] == superID, ][1, ]
    # THEN, SORT TAXON LIST BASED ON TAXONOMY TREE or SORTED TAXON LIST
    if (is.null(sortedTaxonList)) {
        if (is.null(taxaTree)) {
            if (is.null(taxDB)) {
                preCalcTreeFile <- paste(
                    system.file(package = "PhyloProfile"),
                    "PhyloProfile/data/preCalcTree.nw", sep = "/"
                )
            } else preCalcTreeFile <- paste(taxDB, "preCalcTree.nw", sep = "/")
            
            if (file.exists(preCalcTreeFile)) {
                preTree <- ape::read.tree(preCalcTreeFile)
            } else preTree <- createUnrootedTree(Dt)
            
            if (!(repTaxon$abbrName %in% preTree$tip.label))
                message(c(repTaxon$abbrName, " not found in ", preCalcTreeFile))
            taxaTree <- ape::root(
                preTree, outgroup=as.character(repTaxon$abbrName),
                resolve.root=TRUE
            )
        } else
            taxaTree <- ape::root(
                taxaTree, outgroup=as.character(repTaxon$abbrName),
                resolve.root=TRUE
            )
        taxonList <- sortTaxaFromTree(taxaTree)
    } else {
        taxonList <- sortedTaxonList
    }
    sortedDt <- Dt[match(taxonList, Dt$abbrName), ]
    # subset to get list of input taxa only
    sortedDt <- subset(sortedDt, abbrName %in% taxonIDs)
    # get only taxonIDs list of selected rank and rename columns
    sortedOut <- subset(
        sortedDt,
        select = c("abbrName", "ncbiID", "fullName", as.character(rankName)))
    colnames(sortedOut) <- c("abbrName", "species", "fullName", "ncbiID")
    # add name of supertaxa into sortedOut list
    sortedOut <- merge(
        sortedOut, fullnameList, by = "ncbiID", all.x = TRUE, sort = FALSE)
    sortedOut$species <- paste0("ncbi", sortedOut$species)
    ## create new column for sorted supertaxon
    indexSpec <- unlist(lapply(
        seq_len(nlevels(as.factor(sortedOut$fullName.y))),function (x) 100000+x)
    )
    indexSpecDf <- data.frame(
        fullName.y = unique(as.character(sortedOut$fullName.y)),
        sortedSupertaxon = paste0(
            indexSpec, "_", unique(as.character(sortedOut$fullName.y))
        ), stringsAsFactors = FALSE)
    sortedOut <- merge(indexSpecDf, sortedOut, by = "fullName.y")
    # final sorted supertaxa list
    sortedOut$taxonID <- 0
    sortedOut$category <- "cat"
    sortedOut <- sortedOut[, c(
        "abbrName", "taxonID", "fullName.x", "species", "ncbiID",
        "sortedSupertaxon", "rank", "category")]
    colnames(sortedOut) <- c(
        "abbrName", "taxonID", "fullName", "ncbiID", "supertaxonID",
        "supertaxon", "rank", "category")
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
    paralog <- abbrName <- NULL
    if (missing(profileWithTax)) return("No input data given")
    if (missing(taxaCount)) return("No supertaxon count given")
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
    presSpecDt <- presSpecDt[
        , c("geneID", "supertaxon", "presSpec", "freq", "freq.y")
    ]
    colnames(presSpecDt) <- c(
        "geneID", "supertaxon", "presSpec", "presentTaxa", "totalTaxa"
    )
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
    finalPresSpecDt <- finalPresSpecDt[stats::complete.cases(finalPresSpecDt), ]
    # return finalPresSpecDt
    return(finalPresSpecDt)
}

#' Parsing info for phylogenetic profiles
#' @description Creating main dataframe for the input phylogenetic profiles
#' based on selected input taxonomy level (e.g. strain, species) and reference
#' taxon. The output contains the number of paralogs, the max/min/mean/median
#' of VAR1 and VAR2.
#' @usage parseInfoProfile(inputDf, sortedInputTaxa, taxaCount, coorthoCOMax)
#' @param inputDf input profiles in long format
#' @param sortedInputTaxa sorted taxonomy data for the input taxa
#' (check sortInputTaxa())
#' @param taxaCount dataframe counting present taxa in each supertaxon
#' @param coorthoCOMax maximum number of co-orthologs allowed
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
#'     taxonIDs, "class", "Mammalia", NULL, NULL
#' )
#' taxaCount <- plyr::count(sortedInputTaxa, "supertaxon")
#' coorthoCOMax <- 999
#' parseInfoProfile(
#'     mainLongRaw, sortedInputTaxa, taxaCount, coorthoCOMax
#' )

parseInfoProfile <- function(
    inputDf, sortedInputTaxa, taxaCount, coorthoCOMax = 9999
) {
    if (is.null(inputDf) | is.null(sortedInputTaxa))
        stop("Input profiles and sorted taxonomy data cannot be NULL!")
    if (ncol(inputDf) > 3) {
        if (ncol(inputDf) < 5) colnames(inputDf)[4] <- "var1"
        else colnames(inputDf)[c(4,5)] <- c("var1", "var2")
    }
    if (coorthoCOMax < 1) coorthoCOMax <- 1
    # count number of inparalogs & calculate frequency of all supertaxa
    paralogCount <- plyr::count(inputDf, c("geneID", "ncbiID"))
    inputDf <- merge(inputDf, paralogCount, by = c("geneID", "ncbiID"))
    colnames(inputDf)[ncol(inputDf)] <- "paralog"
    # filter by number of coorthologs
    inputDf <- inputDf[!(inputDf$paralog > coorthoCOMax),]
    # merge inputDf and sortedInputTaxa to get taxonomy info
    taxaMdData <- merge(inputDf, sortedInputTaxa, by = "ncbiID")
    # merge with taxaCount to get number of species
    fullMdData <- merge(taxaMdData, taxaCount, by = ("supertaxon"), all.x=TRUE)
    names(fullMdData)[names(fullMdData) == "freq"] <- "numberSpec"
    fullMdData$fullName <- as.vector(fullMdData$fullName)
    names(fullMdData)[names(fullMdData) == "orthoID.x"] <- "orthoID"
    fullMdData <- fullMdData[!duplicated(fullMdData), ]
    return(fullMdData)
}

#' Filter phylogentic profiles
#' @description Create a filtered data needed for plotting or clustering
#' phylogenetic profiles. NOTE: this function require some intermediate steps
#' using the results from other functions. If you would like to get a full
#' processed data from the raw input, please use the function
#' fromInputToProfile() instead!
#' @usage filterProfileData(DF, taxaCount, refTaxon = NULL,
#'     percentCO = c(0, 1), coorthoCOMax = 9999,
#'     var1CO  = c(0, 1), var2CO = c(0, 1), var1Rel = "protein",
#'     var2Rel = "protein", groupByCat = FALSE, catDt = NULL,
#'     var1AggregateBy = "max", var2AggregateBy = "max")
#' @param DF a reduced dataframe contains info for all phylogenetic
#' profiles in the selected taxonomy rank.
#' @param taxaCount dataframe counting present taxa in each supertaxon
#' @param refTaxon selected reference taxon. NOTE: This taxon will not be
#' affected by the filtering. If you want to filter all, set refTaxon <- NULL.
#' Default = NULL.
#' @param percentCO min and max cutoffs for percentage of species present
#' in a supertaxon. Default = c(0, 1).
#' @param coorthoCOMax maximum number of co-orthologs allowed. Default =
#' 9999.
#' @param var1CO min and max cutoffs for var1. Default = c(0, 1).
#' @param var2CO min anc max cutoffs for var2. Default = c(0, 1).
#' @param var1Rel relation of var1 ("protein" for protein-protein or
#' "species" for protein-species). Default = "protein".
#' @param var2Rel relation of var2 ("protein" for protein-protein or
#' "species" for protein-species). Default = "protein".
#' @param groupByCat group genes by their categories (TRUE or FALSE). Default =
#' FALSE.
#' @param catDt dataframe contains gene categories
#' (optional, NULL if groupByCat = FALSE or no info provided). Default = NULL.
#' @param var1AggregateBy aggregate method for VAR1 (max, min, mean
#' or median), applied for calculating var1 of supertaxa. Default = "max".
#' @param var2AggregateBy aggregate method for VAR2 (max, min, mean
#' or median), applied for calculating var2 of supertaxa. Default = "max".
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
#' rankName <- "class"
#' refTaxon <- "Mammalia"
#' percentCutoff <- c(0.0, 1.0)
#' coorthologCutoffMax <- 10
#' var1Cutoff <- c(0.75, 1.0)
#' var2Cutoff <- c(0.5, 1.0)
#' var1Relation <- "protein"
#' var2Relation <- "species"
#' groupByCat <- FALSE
#' catDt <- NULL
#' var1AggregateBy <- "max"
#' var2AggregateBy <- "max"
#' taxonIDs <- levels(as.factor(fullProcessedProfile$ncbiID))
#' sortedInputTaxa <- sortInputTaxa(
#'     taxonIDs, rankName, refTaxon, NULL, NULL
#' )
#' taxaCount <- plyr::count(sortedInputTaxa, "supertaxon")
#' filterProfileData(
#'     fullProcessedProfile,
#'     taxaCount,
#'     refTaxon,
#'     percentCutoff,
#'     coorthologCutoffMax,
#'     var1Cutoff,
#'     var2Cutoff,
#'     var1Relation,
#'     var2Relation,
#'     groupByCat,
#'     catDt,
#'     var1AggregateBy,
#'     var2AggregateBy
#' )

filterProfileData <- function(
    DF, taxaCount, refTaxon = NULL, percentCO = c(0, 1), coorthoCOMax = 9999,
    var1CO = c(0, 1), var2CO = c(0, 1), var1Rel = "protein",
    var2Rel = "protein", groupByCat = FALSE, catDt = NULL,
    var1AggregateBy = "max", var2AggregateBy = "max"
) {
    if (is.null(DF)) stop("Profile data cannot be NULL!")
    if (is.null(refTaxon)) refTaxon <- "NA"
    ### check if working with lowest rank (species/strain), e.g. flag == 0
    flag <- 1
    if (length(unique(levels(as.factor(DF$numberSpec)))) == 1) {
        if (unique(levels(as.factor(DF$numberSpec))) == 1) {
            flag <- 0
        }
    }
    ### remove index from supertaxon names
    DF$taxonMod <- gsub("^[[:digit:]]*_", "", DF$supertaxon)
    ### filter by var1 and var2
    if (flag == 0) {
        if (var1Rel == "protein") {
            DF$var1[DF$taxonMod != refTaxon & DF$var1 < var1CO[1]] <- NA
            DF$var1[DF$taxonMod != refTaxon & DF$var1 > var1CO[2]] <- NA
        } else {
            DF$var1[DF$taxonMod != refTaxon & DF$var1 < var1CO[1]] <- 0
            DF$var1[DF$taxonMod != refTaxon & DF$var1 > var1CO[2]] <- 0
        }
        if (var2Rel == "protein") {
            DF$var2[DF$taxonMod != refTaxon & DF$var2 < var2CO[1]] <- NA
            DF$var2[DF$taxonMod != refTaxon & DF$var2 > var2CO[2]] <- NA
        } else {
            DF$var2[DF$taxonMod != refTaxon & DF$var2 < var2CO[1]] <- 0
            DF$var2[DF$taxonMod != refTaxon & DF$var2 > var2CO[2]] <- 0
        }
    } else {
        if (var1Rel == "protein") {
            DF$var1[DF$var1 < var1CO[1]] <- NA
            DF$var1[DF$var1 > var1CO[2]] <- NA
        } else {
            DF$var1[DF$var1 < var1CO[1]] <- 0
            DF$var1[DF$var1 > var1CO[2]] <- 0
        }
        if (var2Rel == "protein") {
            DF$var2[DF$var2 < var2CO[1]] <- NA
            DF$var2[DF$var2 > var2CO[2]] <- NA
        } else {
            DF$var2[DF$var2 < var2CO[1]] <- 0
            DF$var2[DF$var2 > var2CO[2]] <- 0
        }
    }
    ### calculate % present taxa and filter by percentCO
    DFtmp <- DF[stats::complete.cases(DF), ]
    finalPresSpecDt <- calcPresSpec(DFtmp, taxaCount)
    DF <- Reduce(
        function(x, y) merge(x, y, by = c("geneID", "supertaxon"), all.x=TRUE),
        list(DF, finalPresSpecDt))
    DF$presSpec[DF$presSpec < percentCO[[1]] | DF$presSpec > percentCO[[2]]] <-0
    DF$orthoID[DF$presSpec == 0] <- NA
    if (var1Rel == "protein") {
        DF$orthoID[is.na(DF$var1)] <- NA
        DF$var1[is.na(DF$orthoID)] <- NA
    }
    if (var2Rel == "protein") {
        DF$orthoID[is.na(DF$var2)] <- NA
        DF$var2[is.na(DF$orthoID)] <- NA
    }
    DF$presSpec[is.na(DF$orthoID)] <- 0
    if (var1Rel == "protein")
        DF$var1[DF$presSpec == 0] <- NA
    if (var2Rel == "protein")
        DF$var2[DF$presSpec == 0] <- NA

    ### remove paralog count if NOT working with lowest rank (species/strain)
    if (flag == 1) DF$paralog <- 1

    DF <- droplevels(DF)  # delete unused levels
    DF$geneID <- as.factor(DF$geneID)
    DF$supertaxon <- as.factor(DF$supertaxon)

    # calculate max/min/mean/median VAR1 for every supertaxon of each gene
    DFNoNA <- DF[!is.na(DF$var1), ]
    mVar1Dt <- stats::aggregate(
        DFNoNA[, "var1"],
        list(DFNoNA$supertaxon, DFNoNA$geneID),
        FUN = var1AggregateBy)
    colnames(mVar1Dt) <- c("supertaxon", "geneID", "mVar1")
    # calculate max/min/mean/median VAR2 for each supertaxon
    DFNoNAVar2 <- DF[!is.na(DF$var2), ]
    if (nrow(DFNoNAVar2) > 0) {
        mVar2Dt <- stats::aggregate(
            DFNoNAVar2[, "var2"],
            list(DFNoNAVar2$supertaxon, DFNoNAVar2$geneID),
            FUN = var2AggregateBy)
        colnames(mVar2Dt) <- c("supertaxon", "geneID", "mVar2")
    } else {
        mVar2Dt <- DF[, c("supertaxon", "geneID")]
        mVar2Dt$mVar2 <- 0
    }
    # join mVar2 together with mVar1 scores into one df
    scoreDf <- merge(mVar1Dt, mVar2Dt, by = c("supertaxon","geneID"), all=TRUE)
    # add into DF
    DF <- Reduce(
        function(x, y) merge(x, y, by = c("geneID", "supertaxon"), all.x=TRUE),
        list(DF, scoreDf))

    ### add gene categories (if provided)
    originalOrder <- levels(as.factor(DF$geneID))
    if (groupByCat == TRUE) {
        if (is.null(catDt)) {
            catDt <- data.frame( geneID = levels(DF$geneID))
            catDt$group <- "noCategory"
        }
        dfCat <- data.frame(
            supertaxon = rep(levels(DF$supertaxon), nlevels(DF$geneID)),
            geneID = rep(levels(DF$geneID), each = nlevels(DF$supertaxon)))
        dfCat <- merge(dfCat, catDt, by = "geneID")
        DF <- merge(dfCat, DF, by = c("geneID","supertaxon"), all.x = TRUE)
        DF$category <- DF$group
        DF$geneID <- factor(DF$geneID, levels = originalOrder)
    }
    return(DF) #[!is.na(DF$orthoID),])
}

#' Reduce the filtered profile data into supertaxon level
#' @description Reduce data of the processed phylogenetic profiles from input
#' taxonomy rank into supertaxon level (e.g. from species to phylum)
#' @param filteredProfile dataframe contains the filtered profiles (see
#' ?parseInfoProfile, ?filterProfileData and ?filteredProfile)
#' @return A reduced dataframe contains only profile data for the selected
#' supertaxon rank. This dataframe contains only supertaxa and their value
#' (mVar1 & mVar2) for each gene.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{parseInfoProfile}} for creating a full processed
#' profile dataframe, \code{\link{filterProfileData}} for filter processed
#' profile and \code{\link{filteredProfile}} for a demo filtered
#' profile dataframe
#' @examples
#' data("filteredProfile", package="PhyloProfile")
#' reduceProfile(filteredProfile)

reduceProfile <- function(filteredProfile) {
    if (is.null(filteredProfile)) stop("Profile data cannot be NULL!")
    # check if working with the lowest taxonomy rank; 1 for NO; 0 for YES
    flag <- 1
    if (length(unique(levels(as.factor(filteredProfile$numberSpec)))) == 1) {
        if (unique(levels(as.factor(filteredProfile$numberSpec))) == 1) {
            superDfExt <- filteredProfile[, c(
                "geneID", "supertaxon", "supertaxonID",
                "var1", "presSpec", "category", "orthoID", "var2", "paralog"
            )]
            superDfExt$presentTaxa <- 1
            superDfExt$totalTaxa <- 1
            flag <- 0
        }
    }
    if (flag == 1) {
        # get representative orthoID that has m VAR1 for each supertaxon
        mOrthoID <- filteredProfile[, c(
            "geneID", "supertaxon", "var1", "mVar1", "orthoID", "presSpec"
        )]
        mOrthoID <- subset(mOrthoID, mOrthoID$var1 == mOrthoID$mVar1)
        colnames(mOrthoID) <- c(
            "geneID", "supertaxon", "var1", "mVar1", "orthoID", "presSpec"
        )
        mOrthoID <- mOrthoID[!is.na(mOrthoID$orthoID), ]
        mOrthoID <- mOrthoID[, c("geneID", "supertaxon", "orthoID", "presSpec")]
        mOrthoID <- mOrthoID[!duplicated(mOrthoID[, seq_len(2)]), ]
        # get data set for PhyloProfile plotting (contains only supertaxa info)
        superDf <- subset(filteredProfile, select = c(
            "geneID", "supertaxon", "supertaxonID", "mVar1", "category",
            "mVar2", "paralog", "presentTaxa", "totalTaxa"
        ))
        superDf <- superDf[!duplicated(superDf), ]
        superDfExt <- merge(
            superDf, mOrthoID, by = c("geneID", "supertaxon"), all.x = TRUE
        )
        superDfExt <- superDfExt[, c(
            "geneID", "supertaxon", "supertaxonID",
            "mVar1", "presSpec", "category", "orthoID", "mVar2", "paralog",
            "presentTaxa", "totalTaxa"
        )]
        # rename mVar to var
        names(superDfExt)[names(superDfExt) == "mVar1"] <- "var1"
        names(superDfExt)[names(superDfExt) == "mVar2"] <- "var2"
    }
    return(superDfExt)
}

#' Complete processing of raw input phylogenetic profiles
#' @description Create a processed and filtered data for plotting or analysing
#' phylogenetic profiles from raw input file (from raw input to final filtered
#' dataframe)
#' @usage fromInputToProfile(rawInput, rankName, refTaxon = NULL,
#'     taxaTree = NULL, sortedTaxonList = NULL, var1AggregateBy = "max",
#'     var2AggregateBy = "max", percentCutoff = c(0, 1),
#'     coorthologCutoffMax = 9999, var1Cutoff = c(0, 1), var2Cutoff = c(0, 1),
#'     var1Relation = "protein", var2Relation = "protein", groupByCat = FALSE,
#'     catDt = NULL, taxDB = NULL)
#' @param rawInput input file (in long, wide, multi-fasta or orthoxml format)
#' @param rankName taxonomy rank (e.g. "species","phylum",...)
#' @param refTaxon selected reference taxon name (used for sorting and will be
#' protected from filtering). Default = NULL.
#' @param taxaTree input taxonomy tree for taxa in input profiles (optional).
#' Default = NULL.
#' @param sortedTaxonList list of sorted taxa (optional). Default = NULL.
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
#' @param catDt dataframe contains gene categories. Default = NULL
#' @param taxDB Path to the taxonomy DB files
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
#' sortedTaxonList <- NULL
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
#'     sortedTaxonList,
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
    sortedTaxonList = NULL,
    var1AggregateBy = "max",
    var2AggregateBy = "max",
    percentCutoff = c(0, 1),
    coorthologCutoffMax = 9999,
    var1Cutoff = c(0, 1),
    var2Cutoff = c(0, 1),
    var1Relation = "protein",
    var2Relation = "protein",
    groupByCat = FALSE,
    catDt = NULL,
    taxDB = NULL
) {
    if (missing(rawInput) | missing(rankName)) return("Missing input")
    if (is.null(rawInput) | is.null(rankName)) return("Missing input")
    # convert raw input into long format
    inputDf <- createLongMatrix(rawInput)
    # get input taxon IDs and names
    inputTaxonID <- getInputTaxaID(inputDf)
    # sort input taxa based on selected reference taxon or input taxonomy tree
    sortedInputTaxa <- sortInputTaxa(
        inputTaxonID, rankName, refTaxon, taxaTree, sortedTaxonList, taxDB
    )
    # count present taxa in each supertaxon
    taxaCount <- plyr::count(sortedInputTaxa, "supertaxon")
    # parse info (additional values...) into profile df
    fullMdData <- parseInfoProfile(inputDf, sortedInputTaxa, taxaCount)
    # filter profile
    filteredDf <- filterProfileData(
        fullMdData,
        taxaCount,
        refTaxon,
        percentCutoff,
        coorthologCutoffMax,
        var1Cutoff,
        var2Cutoff,
        var1Relation,
        var2Relation,
        groupByCat = FALSE,
        catDt = NULL,
        var1AggregateBy, var2AggregateBy
    )
    # reduce profile df into supertaxon level
    dataHeat <- reduceProfile(filteredDf)
    return(dataHeat)
}
