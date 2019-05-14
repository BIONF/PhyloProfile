#' Get taxonomy info for a list of taxa
#' @description Get NCBI taxonomy IDs, ranks and names for an input taxon list.
#' @param inputTaxa NCBI ID list of input taxa.
#' @param currentNCBIinfo table/dataframe of the pre-processed NCBI taxonomy
#' data (/PhyloProfile/data/taxonNamesFull.txt)
#' @return A list of 3 dataframes: idList, rankList and reducedInfoList. The
#' "rankList" contains taxon names and all taxonomy ranks of the input taxa
#' including also the noranks from the input rank to the taxonomy root. The
#' "idList" contains input taxon IDs, taxon names, all the ranks from current
#' rank to the taxonomy root together with their IDs (with the format
#' "id#rank"). The reducedInfoList is a subset of taxonNamesFull.txt file,
#' containing the NCBI IDs, taxon fullnames, their current rank and their
#' direct parent ID.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @examples
#' inputTaxa <- c("272557", "176299")
#' ncbiFilein <- system.file(
#'     "extdata", "data/taxonNamesFull.txt",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' currentNCBIinfo <- as.data.frame(data.table::fread(ncbiFilein))
#' getIDsRank(inputTaxa, currentNCBIinfo)

getIDsRank <- function(inputTaxa = NULL, currentNCBIinfo = NULL){
    if (is.null(currentNCBIinfo)) return()
    ## get all taxonomy info for input taxa
    inputTaxaInfo <- pbapply::pblapply(
        seq_len(length(inputTaxa)),
        function (x) {
            refID <- inputTaxa[x]
            # get info for this taxon
            refEntry <- currentNCBIinfo[currentNCBIinfo$ncbiID == refID, ]
            lastID <- refEntry$parentID
            inputTaxaInfo <- refEntry

            while (lastID != 1) {
                nextEntry <- currentNCBIinfo[currentNCBIinfo$ncbiID == lastID, ]
                inputTaxaInfo <- rbindlist(
                    list(inputTaxaInfo, nextEntry),
                    use.names = TRUE,
                    fill = TRUE,
                    idcol = NULL
                )
                lastID <- nextEntry$parentID
            }
            return(inputTaxaInfo)
        }
    )

    ## get reduced taxonomy info (subset of taxonNamesFull.txt)
    reducedInfoDf <- unique(rbindlist(inputTaxaInfo))

    ## get list of all ranks and rank IDs
    rankMod <- NULL
    ncbiID <- NULL
    inputRankIDDf <- lapply(
        seq_len(length(inputTaxaInfo)),
        function (x) {
            inputTaxaInfo[[x]]$rankMod <- inputTaxaInfo[[x]]$rank
            if (inputTaxaInfo[[x]]$rank[1] == "norank") {
                inputTaxaInfo[[x]]$rankMod[1] <-
                    paste0("strain_", inputTaxaInfo[[x]]$ncbiID[1])
            }
            inputTaxaInfo[[x]]$rankMod <- with(
                inputTaxaInfo[[x]],
                ifelse(rankMod == "norank", paste0("norank_", ncbiID), rankMod)
            )
            inputTaxaInfo[[x]]$id <- paste0(
                inputTaxaInfo[[x]]$ncbiID, "#", inputTaxaInfo[[x]]$rankMod
            )
            return(inputTaxaInfo[[x]][ , c("rankMod", "id")])
        }
    )

    ## create matrix of all ranks
    inputRankList <- lapply(
        seq_len(length(inputRankIDDf)),
        function (x) {
            ll <- c(
                paste0("ncbi", inputTaxa[x]),
                reducedInfoDf$fullName[reducedInfoDf$ncbiID == inputTaxa[x]],
                inputRankIDDf[[x]]$rank,
                "norank_1"
            )
            ll <- gsub("strain_[[:digit:]]+", "strain", ll)
            return(data.frame(
                matrix(ll, nrow = 1, byrow = TRUE), stringsAsFactors = FALSE
            ))
        }
    )
    inputRankDf <- do.call(plyr::rbind.fill, inputRankList)

    ## create matrix of all ID#ranks
    inputIDList <- lapply(
        seq_len(length(inputRankIDDf)),
        function (x) {
            ll <- c(
                paste0(
                    reducedInfoDf$fullName[reducedInfoDf$ncbiID ==inputTaxa[x]],
                    "#name"
                ),
                inputRankIDDf[[x]]$id,
                "1#norank_1"
            )
            return(data.frame(
                matrix(ll, nrow = 1, byrow = TRUE), stringsAsFactors = FALSE
            ))
        }
    )
    inputIDDf <- do.call(plyr::rbind.fill, inputIDList)

    # return
    newCol <- seq(ncol(inputIDDf) + 1, ncol(inputRankDf))
    inputIDDf[paste0("X", newCol)] <- NA
    return(list(inputIDDf, inputRankDf, as.data.frame(reducedInfoDf)))
}

#' Indexing all available ranks (including norank)
#' @param rankListFile Input file, where each row is a rank list of a taxon
#' (see rankListFile in example)
#' @importFrom utils read.table
#' @return A dataframe containing a list of all possible ranks and their indexed
#' values.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' rankListFile <- system.file(
#'     "extdata", "data/rankList.txt", package = "PhyloProfile", mustWork = TRUE
#' )
#' rankIndexing(rankListFile)

rankIndexing <- function(rankListFile = NULL){
    if (is.null(rankListFile)) return()
    ### input is a dataframe, where each row is a rank list of a taxon
    rankList <- read.table(
        rankListFile,
        sep = '\t', header = FALSE, fill = TRUE,
        stringsAsFactors = TRUE, na.strings = c("","NA")
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
            subList, function(x) x %in% allInputRank, FUN.VALUE = logical(1)
        )
        subList <- subList[filter]

        ## now go to each rank and check...
        for (i in seq_len(length(subList))) {
            # if it has no index (probably only for norank), ...
            if (is.null(rank2Index[[subList[i]]])) {
                # then set index for this rank = the [available] index
                # of prev rank + 1
                for (j in seq_len(length(subList))) {
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
    index2RankList <- lapply(
        seq_len(length(allInputRank)),
        function (x) {
            data.frame(
                index = rank2Index[[allInputRank[x]]],
                rank = allInputRank[x],
                stringsAsFactors = FALSE
            )
        }
    )
    index2RankDf <- do.call(rbind, index2RankList)
    index2RankDf <- index2RankDf[with(index2RankDf, order(index2RankDf$index)),]

    return(index2RankDf)
}

#' Align NCBI taxonomy IDs of list of taxa into a sorted rank list.
#' @export
#' @param idListFile a text file whose each row is a rank+ID list of a taxon
#' (see idListFile in example)
#' @param rankListFile a text file whose each row is a rank list of a taxon
#' (see rankListFile in example)
#' @importFrom data.table transpose
#' @importFrom reshape2 melt
#' @importFrom stats complete.cases
#' @importFrom utils count.fields
#' @importFrom utils read.table
#' @return An aligned taxonomy dataframe which contains all the available
#' taxonomy ranks from the id and rank list file. This dataframe can be used for
#' creating a well resolved taxonomy tree (see ?createRootedTree) and sorting
#' taxa based on a selected reference taxon (see ?sortInputTaxa).
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{rankIndexing}}, \code{\link{createRootedTree}},
#' \code{\link{sortInputTaxa}}
#' @examples
#' idListFile <- system.file(
#'     "extdata", "data/idList.txt", package = "PhyloProfile", mustWork = TRUE
#' )
#' rankListFile <- system.file(
#'     "extdata", "data/rankList.txt", package = "PhyloProfile", mustWork = TRUE
#' )
#' taxonomyTableCreator(idListFile, rankListFile)

taxonomyTableCreator <- function(idListFile = NULL, rankListFile = NULL){
    if (is.null(idListFile) | is.null(rankListFile)) return()
    index <- NULL
    ### get indexed rank list
    index2RankDf <- rankIndexing(rankListFile)

    ### load idList file
    ncol <- max(count.fields(rankListFile, sep = '\t'))
    idList <- read.table(
        idListFile,
        sep = '\t', header = FALSE, check.names = FALSE, comment.char = "",
        fill = TRUE,
        stringsAsFactors = TRUE, na.strings = c("","NA"),
        col.names = paste0('X', seq_len(ncol))
    )
    colnames(idList)[1] <- "tip"

    ### get ordered rank list
    orderedRank <- factor(index2RankDf$rank, levels = index2RankDf$rank)

    ### create a dataframe containing ordered ranks
    fullRankIDdf <- data.frame(
        rank = matrix(
            unlist(orderedRank), nrow = length(orderedRank), byrow = TRUE
        ),
        stringsAsFactors = FALSE
    )
    fullRankIDdf$index <- as.numeric(rownames(fullRankIDdf))
    fullRankIDdf <- data.table(fullRankIDdf)
    setkey(fullRankIDdf, rank)

    mTaxonDf <- pbapply::pblapply(
        seq_len(nrow(idList)),
        function (x) {
            ### get list of all IDs for this taxon
            taxonDf <- data.frame(idList[x,])
            taxonName <- unlist(
                strsplit(as.character(idList[x,]$tip), "#", fixed = TRUE)
            )
            ### convert into long format
            mTaxonDf <- suppressWarnings(melt(taxonDf, id = "tip"))

            ### get rank names and corresponding IDs
            splitCol <- data.frame(
                do.call(
                    'rbind',
                    strsplit(as.character(mTaxonDf$value), '#', fixed = TRUE)
                )
            )
            mTaxonDf <- cbind(mTaxonDf, splitCol)

            ### remove NA cases
            mTaxonDf <- mTaxonDf[complete.cases(mTaxonDf),]

            ### subselect mTaxonDf to keep only 2 column rank id and rank name
            mTaxonDf <- mTaxonDf[, c("X1","X2")]
            if (mTaxonDf$X2[1] != index2RankDf$rank[1]) {
                mTaxonDf <- rbind(
                    data.frame(
                        "X1" = mTaxonDf$X1[1], "X2" = index2RankDf$rank[1]
                    ),
                    mTaxonDf
                )
            }

            ### rename columns
            colnames(mTaxonDf) <- c(taxonName[1], "rank")

            # return
            mTaxonDf <- data.table(mTaxonDf)
            setkey(mTaxonDf, rank)
            return(mTaxonDf)
        }
    )

    ### merge into data frame contains all available ranks from input
    mTaxonDfFull <- c(list(fullRankIDdf), mTaxonDf)
    fullRankIDdf <- Reduce(
        function (x, y) merge(x, y, all.x = TRUE), mTaxonDfFull
    )

    ### reorder ranks
    fullRankIDdf <- fullRankIDdf[order(fullRankIDdf$index),]

    ### replace NA id by id of previous rank
    fullRankIDdf <- zoo::na.locf(fullRankIDdf)

    ### remove index column
    fullRankIDdf <- subset(fullRankIDdf, select = -c(index))

    ### transpose into wide format
    tFullRankIDdf <- data.table::transpose(fullRankIDdf)

    ### set first row to column names
    colnames(tFullRankIDdf) = as.character(unlist(tFullRankIDdf[1,]))
    tFullRankIDdf <- tFullRankIDdf[-1,]

    ### add "abbrName  ncbiID  fullName" columns
    fullRankIDdfOut <- data.frame(
        abbrName = paste0("ncbi", unlist(tFullRankIDdf[,1])),
        ncbiID = unlist(tFullRankIDdf[,1]),
        fullName = colnames(fullRankIDdf)[-c(1)],
        stringsAsFactors = FALSE
    )
    fullRankIDdfOut <- cbind(fullRankIDdfOut, tFullRankIDdf)

    ### rename last column to "root"
    names(fullRankIDdfOut)[ncol(fullRankIDdfOut)] <- "root"

    ### return
    return(fullRankIDdfOut)
}
