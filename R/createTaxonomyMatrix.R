#' Get taxonomy info for a list of taxa
#' @description Get NCBI taxonomy IDs, ranks and names for an input
#' taxon list
#' @param inputTaxa input taxon list
#' @param currentNCBIinfo table/dataframe of the current taxonomy data
#' (taxonNamesFull.txt)
#' @return A list of 3 dataframes: idList, rankList and reducedInfoList
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @examples
#' inputTaxa <- c("272557", "176299")
#' ncbiFilein <- system.file(
#'     "extdata", "data/taxonNamesFull.txt",
#'     package = "phyloprofile", mustWork = TRUE
#' )
#' currentNCBIinfo <- as.data.frame(data.table::fread(ncbiFilein))
#' getIDsRank(inputTaxa, currentNCBIinfo)

getIDsRank <- function(inputTaxa, currentNCBIinfo){
    titleline <- inputTaxa
    allTaxonInfo <- currentNCBIinfo

    ## get unique norank IDs and unique rank
    ## (which are uninformative for sorting taxa)
    allNorankIDs <- list()
    allRanks <- list()

    for (i in seq_len(length(titleline))) {
        # taxon ID
        refID <- titleline[i]
        # get info for this taxon
        refEntry <- allTaxonInfo[allTaxonInfo$ncbiID == refID, ]
        # parent ID
        lastID <- refEntry$parentID

        if (refEntry$rank == "norank") {
            allNorankIDs <- c(allNorankIDs, refEntry$ncbiID)
        }
        allRanks <- c(allRanks, refEntry$rank)

        while (lastID != 1) {
            nextEntry <- allTaxonInfo[allTaxonInfo$ncbiID == lastID, ]
            if (nextEntry$rank == "norank") {
                allNorankIDs <- c(allNorankIDs, nextEntry$ncbiID)
            }
            allRanks <- c(allRanks, nextEntry$rank)
            lastID <- nextEntry$parentID
        }

        p <- i / length(titleline) * 100
        svMisc::progress(p)
        if (p == 100) message("Getting taxonomy info done!")
    }
    uniqueRank <- names(which(table(as.character(allRanks)) == 1))
    uniqueID <- names(which(table(as.character(allNorankIDs)) == 1))
    overID <- names(
        which(table(as.character(allNorankIDs)) == length(titleline))
    )
    uniqueID <- c(uniqueID, overID)

    ## parse taxonomy info
    rankList <- data.frame()
    idList <- data.frame()
    reducedInfoList <- data.frame()

    for (i in seq_len(length(titleline))) {
        ## taxon ID
        refID <- titleline[i]

        ## get info for this taxon
        refEntry <- allTaxonInfo[allTaxonInfo$ncbiID == refID, ]

        if (
            nrow(reducedInfoList[
                reducedInfoList$X1 == refEntry$ncbiID, ]) == 0
        ) {
            refInfoList <- data.frame(
                matrix(
                    c(
                        refEntry$ncbiID,
                        refEntry$fullName,
                        refEntry$rank,
                        refEntry$parentID
                    ),
                    nrow = 1,
                    byrow = TRUE
                ),
                stringsAsFactors = FALSE
            )
            reducedInfoList <- rbind(reducedInfoList, refInfoList)
        }

        # parentID (used to check if hitting last rank, i.e. norank_1)
        lastID <- refEntry$parentID

        # create list of rank for this taxon
        rank <- c(paste0("ncbi", refID), refEntry$fullName)
        if (refEntry$rank == "norank") {
            rank <- c(rank, paste0("strain"))
        } else {
            rank <- c(rank, refEntry$rank)
        }

        # create list of IDs for this taxon
        ids <- list(paste0(refEntry$fullName, "#name"))
        if (refEntry$rank == "norank") {
            ids <- c(
                ids,
                paste0(
                    refEntry$ncbiID,
                    "#",
                    "strain",
                    "_",
                    refEntry$ncbiID
                )
            )
        } else {
            ids <- c(
                ids,
                paste0(
                    refEntry$ncbiID,
                    "#",
                    refEntry$rank
                )
            )
        }

        # append info into rank and ids
        while (lastID != 1) {
            nextEntry <- allTaxonInfo[allTaxonInfo$ncbiID == lastID, ]

            if (
                nrow(reducedInfoList[
                    reducedInfoList$X1 == nextEntry$ncbiID, ]) == 0
            ) {
                nextEntryList <-
                    data.frame(
                        matrix(
                            c(
                                nextEntry$ncbiID,
                                nextEntry$fullName,
                                nextEntry$rank,
                                nextEntry$parentID
                            ),
                            nrow = 1, byrow = TRUE
                        ),
                        stringsAsFactors = FALSE
                    )

                reducedInfoList <- rbind(
                    reducedInfoList,
                    nextEntryList
                )
            }

            lastID <- nextEntry$parentID

            if ("norank" %in% nextEntry$rank) {
                rank <- c(
                    rank,
                    paste0(
                        nextEntry$rank,
                        "_",
                        nextEntry$ncbiID
                    )
                )
                ids <- c(
                    ids,
                    paste0(
                        nextEntry$ncbiID,
                        "#",
                        nextEntry$rank,
                        "_",
                        nextEntry$ncbiID
                    )
                )
            } else {
                rank <- c(rank, nextEntry$rank)
                ids <- c(
                    ids,
                    paste0(
                        nextEntry$ncbiID,
                        "#",
                        nextEntry$rank
                    )
                )
            }
        }

        # last rank and id
        rank <- c(rank, "norank_1")
        ids <- c(ids, "1#norank_1")

        # change "noRank" before species into "strain"
        if (rank[3] == "norank" & any(rank == "species")) {
            rank[3] = "strain"
            tmpID <- unlist(strsplit(as.character(ids[2]), split = "#"))
            ids[2] <- paste0(tmpID[1],"#","strain")
        }

        # append into rankList and idList files
        rankListTMP <- data.frame(
            matrix(
                unlist(rank), nrow = 1, byrow = TRUE
            ),
            stringsAsFactors = FALSE
        )
        rankList <- plyr::rbind.fill(rankList, rankListTMP)
        idListTMP <- data.frame(
            matrix(unlist(ids), nrow = 1, byrow = TRUE),
            stringsAsFactors = FALSE
        )
        idList <- plyr::rbind.fill(idList, idListTMP)

        # Increment the progress bar, and update the detail text.
        p <- i / length(titleline) * 100
        svMisc::progress(p)
        if (p == 100) message("Creating dataframe done!")
    }

    # return idList, rankList and reducedInfoList
    return(list(idList, rankList, reducedInfoList))
}

#' Indexing all available ranks (including norank)
#' @param rankListFile Input file, where each row is a rank list of a taxon
#' @importFrom utils read.table
#' @return A dataframe containing a list of indexed ranks
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' rankListFile <- system.file(
#'     "extdata", "data/rankList.txt", package = "phyloprofile", mustWork = TRUE
#' )
#' rankIndexing(rankListFile)

rankIndexing <- function(rankListFile){
    ### input is a dataframe, where each row is a rank list of a taxon
    rankList <- as.data.frame(
        read.table(
            rankListFile,
            sep = '\t', header = FALSE, fill = TRUE,
            stringsAsFactors = TRUE, na.strings = c("","NA")
        )
    )

    ### get all available ranks from input rankList
    uList <- unlist(rankList)
    # remove unique rank by replacing with NA (they are useless for sorting)
    uList[!duplicated(uList)] <- NA
    # get final list of available ranks (remove NA items)
    allInputRank <- as.character(unique(uList))
    allInputRank <- allInputRank[!is.na(allInputRank)]

    ### initial index for main ranks
    mainRank <- c(
        "strain","forma","subspecies","varietas",
        "subspecies","species","species subgroup","species group",
        "subgenus","genus","subtribe","tribe",
        "subfamily","family","superfamily",
        "parvorder","infraorder","suborder","order","superorder",
        "cohort","infraclass","subclass","class","superclass",
        "subphylum","phylum","superphylum",
        "subkingdom","kingdom","superkingdom"
    )
    rank2Index <- new.env()
    for (i in seq_len(length(mainRank))) {
        rank2Index[[mainRank[i]]] = i
    }

    ### the magic happens here
    for (k in seq_len(nrow(rankList))) {
        ## get subset of rank list for current taxon
        ## which contains only ranks existing in allInputRank
        subList <- rankList[k,][!is.na(rankList[k,])]
        filter <- vapply(
            subList, function(x) x %in% allInputRank,
            FUN.VALUE = logical(1)
        )
        subList <- subList[filter]

        ## now go to each rank and check...
        for (i in seq_len(length(subList))) {
            # if it has no index (probably only for norank), ...
            if (is.null(rank2Index[[subList[i]]])) {
                # then set index for this rank = the [available] index
                # of prev rank + 1
                for (j in seq_len(length(subList))) {
                    # print(subList[i - j])
                    if (j < i) {
                        if (!is.null(rank2Index[[subList[i - j]]])) {
                            rank2Index[[subList[i]]] =
                                rank2Index[[subList[i - j]]] + 1
                            break
                        } else {
                            j = j - 1
                        }
                    }
                }
            }
            # else, check if the current index is smaller than one of prev rank,
            else {
                if (i > 1) {
                    preRank <- subList[i - 1]
                    # if so, increase index of this current rank by
                    # (index of previous rank + 1)
                    if (rank2Index[[subList[i]]] <= rank2Index[[preRank]]) {
                        rank2Index[[subList[i]]] = rank2Index[[preRank]] + 1
                    }
                }
            }
        }
    }

    ### output a list of indexed ranks
    index2RankDf <- data.frame(
        "index" = character(), "rank" = character(), stringsAsFactors = FALSE
    )

    for (i in seq_len(length(allInputRank))) {
        index2RankDf[i,] = c(rank2Index[[allInputRank[i]]], allInputRank[i])
    }

    index2RankDf$index <- as.numeric(index2RankDf$index)
    index2RankDf <- index2RankDf[with(index2RankDf, order(index2RankDf$index)),]

    return(index2RankDf)
}

#' Align ncbi taxonomy IDs of all taxa into a sorted rank list
#' @export
#' @param idListFile a text file whose each row is a rank+ID list of a taxon
#' @param rankListFile a text file whose each row is a rank list of a taxon
#' @importFrom data.table transpose
#' @importFrom reshape2 melt
#' @importFrom stats complete.cases
#' @importFrom utils count.fields
#' @importFrom utils read.table
#' @return An aligned taxonomy dataframe
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{rankIndexing}}
#' @examples
#' idListFile <- system.file(
#'     "extdata", "data/idList.txt", package = "phyloprofile", mustWork = TRUE
#' )
#' rankListFile <- system.file(
#'     "extdata", "data/rankList.txt", package = "phyloprofile", mustWork = TRUE
#' )
#' taxonomyTableCreator(idListFile, rankListFile)

taxonomyTableCreator <- function(idListFile, rankListFile){
    index <- NULL
    ### get indexed rank list
    index2RankDf <- rankIndexing(rankListFile)

    ### load idList file
    ncol <- max(count.fields(rankListFile, sep = '\t'))
    idList <- as.data.frame(
        read.table(
            idListFile,
            sep = '\t', header = FALSE, check.names = FALSE, comment.char = "",
            fill = TRUE,
            stringsAsFactors = TRUE, na.strings = c("","NA"),
            col.names = paste0('X', seq_len(ncol))
        )
    )

    colnames(idList)[1] <- "tip"

    ### get ordered rank list
    orderedRank <- factor(index2RankDf$rank, levels = index2RankDf$rank)

    ### create a dataframe containing ordered ranks
    fullRankIDdf <- data.frame(
        "rank" = matrix(unlist(orderedRank),
                        nrow = length(orderedRank),
                        byrow = TRUE),
        stringsAsFactors = FALSE
    )
    fullRankIDdf$index <- as.numeric(rownames(fullRankIDdf))

    for (i in seq_len(nrow(idList))) {
        ### get list of all IDs for this taxon
        taxonDf <- data.frame(idList[i,])
        taxonName <- unlist(
            strsplit(as.character(idList[i,]$tip), "#", fixed = TRUE)
        )
        ### convert into long format
        mTaxonDf <- suppressWarnings(melt(taxonDf, id = "tip"))

        ### get rank names and corresponding IDs
        splitCol <- data.frame(
            do.call(
                'rbind',
                strsplit(
                    as.character(mTaxonDf$value), '#',
                    fixed = TRUE
                )
            )
        )
        mTaxonDf <- cbind(mTaxonDf, splitCol)

        ### remove NA cases
        mTaxonDf <- mTaxonDf[complete.cases(mTaxonDf),]

        ### subselect mTaxonDf to keep only 2 column rank id and rank name
        mTaxonDf <- mTaxonDf[, c("X1","X2")]
        if (mTaxonDf$X2[1] != index2RankDf$rank[1]) {
            mTaxonDf <- rbind(
                data.frame("X1" = mTaxonDf$X1[1], "X2" = index2RankDf$rank[1]),
                mTaxonDf
            )
        }

        ### rename columns
        colnames(mTaxonDf) <- c(taxonName[1], "rank")

        ### merge with index2RankDf (contains all available ranks from input)
        fullRankIDdf <- merge(
            fullRankIDdf, mTaxonDf, by = c("rank"), all.x = TRUE
        )

        ### reorder ranks
        fullRankIDdf <- fullRankIDdf[order(fullRankIDdf$index),]

        ### replace NA id by id of previous rank
        fullRankIDdf <- zoo::na.locf(fullRankIDdf)

        ### print process
        p <- i / nrow(idList) * 100
        svMisc::progress(p)
        if (p == 100) message("Creating taxonomy matrix done!")
    }

    ### remove index column
    fullRankIDdf <- subset(fullRankIDdf, select = -c(index))

    ### transpose into wide format
    tFullRankIDdf <- data.table::transpose(fullRankIDdf)

    ### set first row to column names
    colnames(tFullRankIDdf) = as.character(unlist(tFullRankIDdf[1,]))
    tFullRankIDdf <- tFullRankIDdf[-1,]

    ### add "abbrName  ncbiID  fullName" columns
    abbrName <- paste0("ncbi", tFullRankIDdf[,1])
    ncbiID <- tFullRankIDdf[,1]
    fullName <- colnames(fullRankIDdf)[-c(1)]

    fullRankIDdf <- cbind(abbrName, ncbiID)
    fullRankIDdf <- cbind(fullRankIDdf, fullName)
    fullRankIDdf <- cbind(fullRankIDdf, tFullRankIDdf)

    ### rename last column to "root"
    names(fullRankIDdf)[ncol(fullRankIDdf)] <- "root"

    ### return
    return(fullRankIDdf)
}

#' TAXA2DIST
#' @importFrom stats as.dist
#' @param x taxa matrix
#' @param varstep var-step
#' @param check check
#' @param labels labels
#' @return a distance matrix
#' @author function from taxize library

taxa2dist <- function(x, varstep = FALSE, check = TRUE, labels) {
    rich <- apply(x, 2, function(taxa) length(unique(taxa)))
    S <- nrow(x)
    if (check) {
        keep <- rich < S & rich > 1
        rich <- rich[keep]
        x <- x[, keep]
    }
    i <- rev(order(rich))
    x <- x[, i]
    rich <- rich[i]
    if (varstep) {
        add <- -diff(c(nrow(x), rich, 1))
        add <- add/c(S, rich)
        add <- add/sum(add) * 100
    }
    else {
        add <- rep(100/(ncol(x) + check), ncol(x) + check)
    }
    if (!is.null(names(add)))
        names(add) <- c("Base", names(add)[-length(add)])
    if (!check)
        add <- c(0, add)
    out <- matrix(add[1], nrow(x), nrow(x))
    for (i in seq_len(ncol(x))) {
        out <- out + add[i + 1] * outer(x[, i], x[, i], "!=")
    }
    out <- as.dist(out)
    attr(out, "method") <- "taxa2dist"
    attr(out, "steps") <- add
    if (missing(labels)) {
        attr(out, "Labels") <- rownames(x)
    }
    else {
        if (length(labels) != nrow(x))
            warning("Labels are wrong: needed ", nrow(x), " got ",
                    length(labels))
        attr(out, "Labels") <- as.character(labels)
    }
    if (!check && any(out <= 0))
        warning("you used 'check=FALSE' and some distances are zero
                        -- was this intended?")
    out
}
