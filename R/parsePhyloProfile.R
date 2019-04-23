# Functions for parsing and pre-processing PhyloProfile input

#' Get list of pre-installed NCBI taxon names
#' @description Get all NCBI taxon names from "data/taxonNamesReduced.txt"
#' @export
#' @return List of taxon IDs, their full names, taxonomy ranks and parent IDs
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
                fileURL,
                destfile = nameReducedFile,
                method="auto"
            ),
            error=function(e) 1
        )
    }

    nameList <- as.data.frame(read.table(
        nameReducedFile,
        sep = "\t",
        header = TRUE,
        fill = TRUE
    ))
    nameList$fullName <- as.character(nameList$fullName)
    nameList$rank <- as.character(nameList$rank)
    nameList <- nameList[!duplicated(nameList), ]

    return(nameList)
}

#' Get taxonomy matrix
#' @description Get full taxonomy matrix from "data/taxonomyMatrix.txt" or
#' only a subset of matrix based on an input taxon list
#' @export
#' @param subsetTaxaCheck subset taxonomy matrix based on input taxon IDs
#' (TRUE/FALSE)
#' @param taxonIDs list of input taxon IDs
#' @return Data frame contains the (subset of) taxonomy matrix
#' @importFrom utils download.file
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' # get full pre-installed taxonomy matrix
#' getTaxonomyMatrix(FALSE, NULL)
#' # get taxonomy matrix for a list of taxon IDs
#' taxonIDs <- c("ncbi10020", "ncbi10090")
#' getTaxonomyMatrix(TRUE, taxonIDs)

getTaxonomyMatrix <- function(subsetTaxaCheck, taxonIDs){
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
                fileURL,
                destfile = taxonomyMatrixFile,
                method="auto"
            ),
            error=function(e) 1
        )
    }

    dt <- as.data.frame(read.table(
        taxonomyMatrixFile,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = TRUE
    ))
    if (subsetTaxaCheck) {
        dt <- dt[dt$abbrName  %in% taxonIDs, ]
    }
    return(dt)
}

#' Get ID list of input taxa from the main input
#' @param rawProfile long dataframe of input profile
#' @return List of all input taxon IDs
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{createLongMatrix}}, \code{\link{mainLongRaw}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' getInputTaxaID(mainLongRaw)

getInputTaxaID <- function(rawProfile){
    if (is.null(rawProfile)) return()
    inputTaxa <- levels(as.factor(rawProfile$ncbiID))

    inputTaxa <- unlist(strsplit(inputTaxa, split = "\t"))
    if (inputTaxa[1] == "geneID") {
        # remove "geneID" element from vector inputTaxa
        inputTaxa <- inputTaxa[-1]
    }
    # return input taxa
    return(inputTaxa)
}

#' Get NCBI taxon names for a selected list of taxa
#' @description Get NCBI taxon names from "data/taxonNamesReduced.txt" for
#' a selected list of taxon
#' @param rankName taxonomy rank (e.g. "species","phylum",...)
#' @param taxonIDs list of taxon IDs (check getInputTaxaID())
#' @return List of full names, taxonomy ranks and parent IDs for the input taxa
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{getInputTaxaID}} for getting input taxon IDs,
#' \code{\link{getNameList}} for getting the full taxon name list
#' @examples
#' taxonIDs <- c("ncbi10020", "ncbi10090")
#' getInputTaxaName("species", taxonIDs)

getInputTaxaName <- function(rankName, taxonIDs){
    # load list of unsorted taxa
    Dt <- getTaxonomyMatrix(TRUE, taxonIDs)

    # load list of taxon name
    nameList <- getNameList()

    choice <- as.data.frame
    choice <- rbind(Dt[rankName])
    colnames(choice) <- "ncbiID"
    choice <- merge(choice,
                    nameList,
                    by = "ncbiID",
                    all = FALSE)
    return(choice)
}

#' Sort list of (super)taxa based on a selected reference (super)taxon
#' @param taxonIDs list of taxon IDs (e.g.: ncbi1234, ncbi9999, ...)
#' @param taxonNames list of taxon names and their corresponding taxonomy
#' ranks, IDs,...
#' @param rankName working taxonomy rank (e.g. "species", "phylum",...)
#' @param refTaxon selected reference taxon
#' @param taxaTree input taxonomy tree (optional)
#' @return Taxonomy matrix for the input taxa ordered by the selected reference
#' taxon. This matrix is sorted either based on the NCBI taxonomy info, or
#' based on an user-defined taxonomy tree (if provided).
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
#' taxonNames <- getInputTaxaName("species", taxonIDs)
#' sortInputTaxa(taxonIDs, taxonNames, "species", "Homo sapiens", NULL)

sortInputTaxa <- function(taxonIDs,
                            taxonNames,
                            rankName,
                            refTaxon,
                            taxaTree){
    ncbiID <- NULL
    fullName <- NULL
    abbrName <- NULL

    fullnameList <- getNameList()

    # get selected supertaxon ID(s)
    rankNameTMP <- taxonNames$rank[taxonNames$fullName == refTaxon]
    if (rankName == "strain") {
        superID <- fullnameList$ncbiID[
            fullnameList$fullName == refTaxon
            & fullnameList$rank == "norank"
            ]
    } else {
        superID <- fullnameList$ncbiID[
            fullnameList$fullName == refTaxon 
            & fullnameList$rank == rankNameTMP[1]
            ]
    }

    # get full taxonomy data
    Dt <- getTaxonomyMatrix(FALSE, NULL)

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
        select = c(
            "abbrName", "ncbiID", "fullName", as.character(rankName)
        )
    )
    colnames(sortedOut) <- c("abbrName", "species", "fullName", "ncbiID")

    # add name of supertaxa into sortedOut list
    sortedOut <- merge(
        sortedOut, fullnameList,
        by = "ncbiID",
        all.x = TRUE,
        sort = FALSE
    )
    sortedOut$species <- as.character(sortedOut$species)

    # add orderPrefix to supertaxon name
    # and add prefix "ncbi" to taxonNcbiID (column "species")
    prefix <- 1001

    ## create new column for sorted supertaxon
    sortedOut$sortedSupertaxon <- 0
    sortedOut$sortedSupertaxon[1] <- paste0(prefix,
                                            "_",
                                            sortedOut$fullName.y[1])
    sortedOut$species[1] <- paste0("ncbi", sortedOut$species[1])

    if (nrow(sortedOut) > 1) {
        for (i in 2:nrow(sortedOut)) {
            ## increase prefix if changing to another supertaxon
            if (sortedOut$fullName.y[i] != sortedOut$fullName.y[i - 1]) {
                prefix <- prefix + 1
            }
            sortedOut$sortedSupertaxon[i] <- paste0(prefix,
                                                    "_",
                                                    sortedOut$fullName.y[i])
            sortedOut$species[i] <- paste0("ncbi", sortedOut$species[i])
        }
    }

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

#' Calculate percentage of present species in each ortholog group
#' @export
#' @param profileWithTax long data frame of main PhyloProfile input together
#' with their taxonomy info
#' @param taxaCount number of species occur in each supertaxon
#' @return A data frame with % of present species for each seed protein in
#' each supertaxon
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
    profileWithTax <- profileWithTax[profileWithTax$orthoID != "NA", ]

    # get geneID and supertaxon
    geneIDSupertaxon <- subset(
        profileWithTax,
        select = c("geneID", "supertaxon", "paralog", "abbrName")
    )
    # remove duplicated rows
    geneIDSupertaxon <- geneIDSupertaxon[!duplicated(geneIDSupertaxon), ]

    # remove NA rows from profileWithTax
    profileWithTaxNoNA <-
        profileWithTax[profileWithTax$orthoID != "NA", ]

    # count present frequency of supertaxon for each gene
    geneSupertaxonCount <- plyr::count(
        profileWithTaxNoNA,
        c("geneID", "supertaxon")
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
    finalPresSpecDt <- merge(presSpecDt,
                                geneIDSupertaxon,
                                by = c("geneID", "supertaxon"),
                                all.y = TRUE)
    finalPresSpecDt$presSpec[is.na(finalPresSpecDt$presSpec)] <- 0

    # remove duplicated rows
    finalPresSpecDt <- finalPresSpecDt[!duplicated(finalPresSpecDt), ]

    # return finalPresSpecDt
    return(finalPresSpecDt)
}

#' Parsing info for phylogenetic profiles
#' @description Creating main dataframe for the input phylogenetic profiles with
#' the selected input taxonomy level (e.g. strain, species) and reference taxon.
#' The output contains the number of paralogs, percentage of species presence
#' in each supertaxon, and the max/min/mean/median of VAR1 and VAR2.
#' @usage parseInfoProfile( inputDf, sortedInputTaxa, var1AggregateBy,
#'     var2AggregateBy)
#' @param inputDf input profiles in long format
#' @param sortedInputTaxa sorted taxonomy data for the input taxa
#' (check sortInputTaxa())
#' @param var1AggregateBy aggregate method for VAR1 (min, max, mean or median)
#' @param var2AggregateBy aggregate method for VAR2 (min, max, mean or median)
#' @importFrom stats aggregate
#' @return Dataframe contains all info for input profiles (a full processed
#' profile)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{createLongMatrix}}, \code{\link{sortInputTaxa}},
#' \code{\link{calcPresSpec}}, \code{\link{mainLongRaw}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' inputDf <- mainLongRaw
#' taxonIDs <- getInputTaxaID(inputDf)
#' taxonNames <- getInputTaxaName("class", taxonIDs)
#' sortedInputTaxa <- sortInputTaxa(
#'     taxonIDs, taxonNames, "class", "Mammalia", NULL
#' )
#' var1AggregateBy <- "max"
#' var2AggregateBy <- "mean"
#' parseInfoProfile(
#'     inputDf, sortedInputTaxa, var1AggregateBy, var2AggregateBy
#' )

parseInfoProfile <- function(
    inputDf,
    sortedInputTaxa,
    var1AggregateBy,
    var2AggregateBy
) {
    mdData <- inputDf
    # rename columns of 2 additional variables
    if (ncol(mdData) > 3) {
        if (ncol(mdData) < 5) {
            colnames(mdData)[4] <- "var1"
        } else {
            colnames(mdData)[c(4,5)] <- c("var1", "var2")
        }
    }

    # count number of inparalogs
    paralogCount <- plyr::count(mdData, c("geneID", "ncbiID"))
    mdData <- merge(mdData, paralogCount, by = c("geneID", "ncbiID"))
    colnames(mdData)[ncol(mdData)] <- "paralog"

    # (1) GET SORTED TAXONOMY LIST
    taxaList <- sortedInputTaxa

    # calculate frequency of all supertaxa
    taxaCount <- plyr::count(taxaList, "supertaxon")

    # merge mdData, mdDataVar2 and taxaList to get taxonomy info
    taxaMdData <- merge(mdData, taxaList, by = "ncbiID")

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
    presMdData <- merge(taxaMdData,
                        finalPresSpecDt,
                        by = c("geneID", "supertaxon"),
                        all.x = TRUE)
    fullMdData <- merge(presMdData,
                        scoreDf,
                        by = c("geneID", "supertaxon"),
                        all.x = TRUE)
    fullMdData <- merge(fullMdData,
                        taxaCount, by = ("supertaxon"),
                        all.x = TRUE)
    # rename "freq" into "numberSpec"
    names(fullMdData)[names(fullMdData) == "freq"] <- "numberSpec"

    fullMdData$fullName <- as.vector(fullMdData$fullName)
    names(fullMdData)[names(fullMdData) == "orthoID.x"] <- "orthoID"

    # parsed input data frame and return
    fullMdData <- fullMdData[!duplicated(fullMdData), ]
    return(fullMdData)
}

#' Reduce the full processed profile into supertaxon level
#' @description Reduce data of the phylogenetic profiles from input taxonomy
#' rank into supertaxon level (e.g. from species to phylum)
#' @param fullProfile dataframe contains the full processed profiles
#' @return A reduced dataframe contains only profiles for the selected rank.
#' This dataframe contains only supertaxa and their value (%present, mVar1 &
#' mVar2) for each gene.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{parseInfoProfile}} for creating a full processed
#' profile dataframe, \code{\link{fullProcessedProfile}} for a demo full
#' processed profile dataframe
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' reduceProfile(fullProcessedProfile)

reduceProfile <- function(fullProfile) {
    fullMdData <- fullProfile

    # to check if working with the lowest taxonomy rank; 1 for NO; 0 for YES
    flag <- 1
    if (length(unique(levels(as.factor(fullMdData$numberSpec)))) == 1) {
        if (unique(levels(as.factor(fullMdData$numberSpec))) == 1) {
            superDfExt <- fullMdData[, c(
                "geneID",
                "supertaxon",
                "supertaxonID",
                "var1",
                "presSpec",
                "category",
                "orthoID",
                "var2",
                "paralog"
            )]
            flag <- 0
        }
    }

    if (flag == 1) {
        # get representative orthoID that has m VAR1 for each supertaxon
        mOrthoID <- fullMdData[, c(
            "geneID",
            "supertaxon",
            "var1",
            "mVar1",
            "orthoID"
        )]
        mOrthoID <- subset(mOrthoID, mOrthoID$var1 == mOrthoID$mVar1)
        colnames(mOrthoID) <- c(
            "geneID",
            "supertaxon",
            "var1",
            "mVar1",
            "orthoID"
        )
        mOrthoID <- mOrthoID[!is.na(mOrthoID$orthoID), ]
        mOrthoID <- mOrthoID[, c("geneID", "supertaxon", "orthoID")]
        mOrthoID <- mOrthoID[!duplicated(mOrthoID[, seq_len(2)]), ]

        # get data set for PhyloProfile plotting (contains only supertaxa info)
        superDf <- subset(fullMdData, select = c(
            "geneID",
            "supertaxon",
            "supertaxonID",
            "mVar1",
            "presSpec",
            "category",
            "mVar2",
            "paralog"
        ))
        superDf$paralog <- 1
        superDf <- superDf[!duplicated(superDf), ]

        superDfExt <- merge(superDf, mOrthoID, by = c("geneID", "supertaxon"),
                            all.x = TRUE)
        superDfExt <- superDfExt[, c(
            "geneID",
            "supertaxon",
            "supertaxonID",
            "mVar1",
            "presSpec",
            "category",
            "orthoID",
            "mVar2",
            "paralog"
        )]

        # rename mVar to var
        names(superDfExt)[names(superDfExt) == "mVar1"] <- "var1"
        names(superDfExt)[names(superDfExt) == "mVar2"] <- "var2"
    }

    return(superDfExt)
}

#' Create data for plotting the phylogentic profiles
#' @usage createProfileData(superTaxonData, refTaxon, percentCutoff,
#'     coorthologCutoffMax, var1Cutoff, var2Cutoff, var1Relation,
#'     var2Relation, groupByCat, catDt)
#' @param superTaxonData a reduced dataframe contains info for all profiles in
#' the selected taxonomy rank.
#' @param refTaxon selected reference taxon
#' @param percentCutoff min and max cutoffs for percentage of species present
#' in a supertaxon
#' @param coorthologCutoffMax maximum number of co-orthologs allowed
#' @param var1Cutoff min and max cutoffs for var1
#' @param var2Cutoff min anc max cutoffs for var2
#' @param var1Relation relation of var1 ("protein" for protein-protein or
#' "species" for protein-species)
#' @param var2Relation relation of var2 ("protein" for protein-protein or
#' "species" for protein-species)
#' @param groupByCat group genes by their categories (TRUE or FALSE)
#' @param catDt dataframe contains gene categories
#' (optional, NULL if groupByCat = FALSE or no info provided)
#' @return A dataframe ready for generating profile plot.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{parseInfoProfile}} and \code{\link{reduceProfile}}
#' for generating input dataframe, \code{\link{fullProcessedProfile}} for a
#' demo full processed profile dataframe
#' @examples
#' data("fullProcessedProfile", package="PhyloProfile")
#' superTaxonData <- reduceProfile(fullProcessedProfile)
#' refTaxon <- "Mammalia"
#' percentCutoff <- c(0.0, 1.0)
#' coorthologCutoffMax <- 10
#' var1Cutoff <- c(0.75, 1.0)
#' var2Cutoff <- c(0.5, 1.0)
#' var1Relation <- "protein"
#' var2Relation <- "species"
#' groupByCat <- FALSE
#' catDt <- NULL
#' createProfileData(
#'     superTaxonData,
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

createProfileData <- function(
    superTaxonData,
    refTaxon,
    percentCutoff,
    coorthologCutoffMax,
    var1Cutoff,
    var2Cutoff,
    var1Relation,
    var2Relation,
    groupByCat,
    catDt
) {
    dataHeat <- superTaxonData

    # cutoffs
    percentCutoffMin <- percentCutoff[1]
    percentCutoffMax <- percentCutoff[2]
    var1CutoffMin <- var1Cutoff[1]
    var1CutoffMax <- var1Cutoff[2]
    var2CutoffMin <- var2Cutoff[1]
    var2CutoffMax <- var2Cutoff[2]

    # get selected supertaxon name
    inSelect <- refTaxon

    ### replace insufficient values according to the thresholds by NA or 0
    # based on presSpec or # of co-orthologs
    numberCoortholog <- levels(as.factor(dataHeat$paralog))
    if (length(numberCoortholog) > 1) {
        dataHeat$presSpec[
            dataHeat$supertaxon != inSelect
            & dataHeat$paralog > coorthologCutoffMax
        ] <- 0
        if (var2Relation == "protein") {
            dataHeat$var2[
                dataHeat$supertaxon != inSelect
                & dataHeat$paralog > coorthologCutoffMax
            ] <- NA
        }
    } else {
        if (length(levels(as.factor(dataHeat$presSpec))) > 1)
            dataHeat$presSpec[
                dataHeat$supertaxon != inSelect
                & dataHeat$presSpec < percentCutoffMin
            ] <- 0
        dataHeat$presSpec[
            dataHeat$supertaxon != inSelect
            & dataHeat$presSpec > percentCutoffMax
        ] <- 0
        if (var2Relation == "protein") {
            dataHeat$var2[
                dataHeat$supertaxon != inSelect
                & dataHeat$presSpec < percentCutoffMin
            ] <- NA
            dataHeat$var2[
                dataHeat$supertaxon != inSelect
                & dataHeat$presSpec > percentCutoffMax
            ] <- NA
        }
    }

    dataHeat$presSpec[
        dataHeat$supertaxon != inSelect
        & dataHeat$var1 < var1CutoffMin
    ] <- 0
    dataHeat$presSpec[
        dataHeat$supertaxon != inSelect
        & dataHeat$var1 > var1CutoffMax
    ] <- 0
    if (var1Relation == "protein") {
        if (var2Relation == "protein") {
            # prot-prot: remove complete cell if one variable not sufficient
            dataHeat$presSpec[
                dataHeat$supertaxon != inSelect
                & dataHeat$var2 < var2CutoffMin
            ] <- 0
            dataHeat$presSpec[
                dataHeat$supertaxon != inSelect
                & dataHeat$var2 > var2CutoffMax
            ] <- 0
            dataHeat$var2[
                dataHeat$supertaxon != inSelect
                & dataHeat$var1 < var1CutoffMin
            ] <- NA
            dataHeat$var2[
                dataHeat$supertaxon != inSelect
                & dataHeat$var1 > var1CutoffMax
            ] <- NA
            dataHeat$var1[
                dataHeat$supertaxon != inSelect
                & dataHeat$var2 < var2CutoffMin
            ] <- NA
            dataHeat$var1[
                dataHeat$supertaxon != inSelect
                & dataHeat$var2 > var2CutoffMax
            ] <- NA
        } else {
            # prot-spec: var1 depend on var2
            dataHeat$presSpec[
                dataHeat$supertaxon != inSelect
                & dataHeat$var2 < var2CutoffMin
            ] <- 0
            dataHeat$presSpec[
                dataHeat$supertaxon != inSelect
                & dataHeat$var2 > var2CutoffMax
            ] <- 0
        }
    } else {
        if (var2Relation == "species") {
            # spec-spec: remove var1 and var2 independently
            dataHeat$presSpec[
                dataHeat$supertaxon != inSelect
                & dataHeat$var1 < var1CutoffMin
            ] <- 0
            dataHeat$presSpec[
                dataHeat$supertaxon != inSelect
                & dataHeat$var1 > var1CutoffMax
            ] <- 0
        } else {
            # spec-prot: var2 depend on var1
            dataHeat$var2[
                dataHeat$supertaxon != inSelect
                & dataHeat$var1 < var1CutoffMin
            ] <- NA
            dataHeat$var2[
                dataHeat$supertaxon != inSelect
                & dataHeat$var1 > var1CutoffMax
            ] <- NA
        }
    }

    dataHeat$var1[
        dataHeat$supertaxon != inSelect & dataHeat$var1 < var1CutoffMin
    ] <- NA
    dataHeat$var1[
        dataHeat$supertaxon != inSelect & dataHeat$var1 > var1CutoffMax
    ] <- NA
    dataHeat$var2[
        dataHeat$supertaxon != inSelect & dataHeat$var2 < var2CutoffMin
    ] <- NA
    dataHeat$var2[
        dataHeat$supertaxon != inSelect & dataHeat$var2 > var2CutoffMax
    ] <- NA

    dataHeat <- droplevels(dataHeat)  # delete unused levels
    dataHeat$geneID <- as.factor(dataHeat$geneID)
    dataHeat$supertaxon <- as.factor(dataHeat$supertaxon)

    ### add gene categories (if provided)
    if (groupByCat == TRUE) {
        if (is.null(catDt)) {
            catDt <- data.frame( geneID = levels(dataHeat$geneID))
            catDt$group <- "noCategory"
        }

        # create a dataframe that contain all genes and all taxa
        dataHeatCat <- data.frame(
            supertaxon = rep(
                levels(dataHeat$supertaxon), nlevels(dataHeat$geneID)
            ),
            geneID = rep(
                levels(dataHeat$geneID), each = nlevels(dataHeat$supertaxon)
            )
        )

        dataHeatCat <- merge(dataHeatCat, catDt, by = "geneID")

        # add categories into dataHeat
        dataHeat <- merge(
            dataHeatCat, dataHeat, by = c("geneID","supertaxon"), all.x = TRUE
        )
    }

    return(dataHeat)
}

#' Create data for plotting profiles (from raw input to final dataframe)
#' @description Create data needed for plotting phylogenetic profile
#' from raw input file.
#' @usage fromInputToProfile(rawInput, rankName, refTaxon, taxaTree,
#'     var1AggregateBy, var2AggregateBy, percentCutoff,
#'     coorthologCutoffMax, var1Cutoff, var2Cutoff, var1Relation,
#'     var2Relation, groupByCat, catDt)
#' @param rawInput input file (in long, wide, multi-fasta or orthoxml format)
#' @param rankName taxonomy rank (e.g. "species","phylum",...)
#' @param refTaxon selected reference taxon
#' @param taxaTree input taxonomy tree (optional)
#' @param var1AggregateBy aggregate method for VAR1 (min, max, mean or median)
#' @param var2AggregateBy aggregate method for VAR2 (min, max, mean or median)
#' @param percentCutoff min and max cutoffs for percentage of species present
#' in a supertaxon
#' @param coorthologCutoffMax maximum number of co-orthologs allowed
#' @param var1Cutoff min and max cutoffs for var1
#' @param var2Cutoff min and max cutoffs for var2
#' @param var1Relation relation of var1 ("protein" for protein-protein or
#' "species" for protein-species)
#' @param var2Relation relation of var2 ("protein" for protein-protein or
#' "species" for protein-species)
#' @param groupByCat group genes by their categories (TRUE or FALSE)
#' @param catDt dataframe contains gene categories
#' @return dataframe for generating profile plot
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @seealso \code{\link{createLongMatrix}}, \code{\link{getInputTaxaID}},
#' \code{\link{getInputTaxaName}}, \code{\link{sortInputTaxa}},
#' \code{\link{parseInfoProfile}}, \code{\link{reduceProfile}},
#' \code{\link{createProfileData}}
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
    refTaxon,
    taxaTree,
    var1AggregateBy,
    var2AggregateBy,
    percentCutoff,
    coorthologCutoffMax,
    var1Cutoff,
    var2Cutoff,
    var1Relation,
    var2Relation,
    groupByCat,
    catDt
) {
    # convert raw input into long format
    inputDf <- createLongMatrix(rawInput)

    # get input taxon IDs and names
    inputTaxonID <- getInputTaxaID(inputDf)
    inputTaxonName <- getInputTaxaName(rankName, inputTaxonID)

    # sort input taxa based on selected reference taxon or input taxonomy tree
    sortedtaxaList <- sortInputTaxa(
        inputTaxonID, inputTaxonName, rankName, refTaxon, taxaTree
    )

    # parse info (additional values...) into profile df
    fullMdData <- parseInfoProfile(
        inputDf,
        sortedInputTaxa = sortedtaxaList,
        var1AggregateBy, var2AggregateBy
    )

    # reduce profile df into supertaxon level
    dataSupertaxa <- reduceProfile(fullMdData)

    # cutoffs
    percentCutoffMin <- percentCutoff[1]
    percentCutoffMax <- percentCutoff[2]
    var1CutoffMin <- var1Cutoff[1]
    var1CutoffMax <- var1Cutoff[2]
    var2CutoffMin <- var2Cutoff[1]
    var2CutoffMax <- var2Cutoff[2]

    # create final df
    dataHeat <- createProfileData(superTaxonData = dataSupertaxa,
                                    refTaxon = refTaxon,
                                    percentCutoff,
                                    coorthologCutoffMax,
                                    var1Cutoff,
                                    var2Cutoff,
                                    var1Relation,
                                    var2Relation,
                                    groupByCat = FALSE,
                                    catDt = NULL)

    return(dataHeat)
}
