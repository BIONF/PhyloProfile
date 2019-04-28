# Get OMA information functions
# using OmaDB package (https://github.com/klarakaleb/OmaDB)

#' Check OMA IDs
#' @description Check if input IDs are valid IDs for OMA Browser
#' (either OMA IDs or UniProt IDs)
#' @export
#' @param ids list of ids needs to be checked
#' @return list of invalid IDs (not readable for OMA)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' checkOmaID("HUMAN29398")

checkOmaID <- function(ids) {
    ids <- as.character(ids)
    getOma <- OmaDB::getProtein(ids)
    return(setdiff(ids, names(getOma)))
}

#' Get OMA members
#' @description Get OMA orthologs for a seed protein from OMA Browser
#' @export
#' @param id ID of the seed protein (OMA or UniProt ID)
#' @param orthoType type of OMA orthologs: either "HOG", "OG"
#' (orthologous group) or "PAIR" (orthologous pair - CURRENTLY NOT WORKING)
#' @return list of ortholog members
#' @author Carla Mölbert {carla.moelbert@gmx.de}
#' @examples
#' getOmaMembers("HUMAN29397", "OG")

getOmaMembers <- function(id, orthoType) {
    # get the members of the Hierarchical Orthologous Group
    if (orthoType == "HOG") {
        members <- suppressWarnings(
            OmaDB::getHOG(id = id, level = "root", members = TRUE)$members$omaid
        )
    }
    # get the members of the Ortholoug group
    else if (orthoType == "OG") {
        members <- suppressWarnings(
            OmaDB::getOMAGroup(id = id)$members$omaid
        )
    }
    # # get the members of the Orthologous Pair
    # else if (orthoType == "PAIR") {
    #   members <- OmaDB::resolveURL(OmaDB::getData(type = "protein",
    #                                               id = id)$orthologs)$omaid
    #   # add query ID into output list
    #   seed <- OmaDB::getData("protein",id)$omaid
    #   members <- c(seed,members)
    # }

    return(members)
}

#' Get domain annotation from OMA Browser
#' @export
#' @description Get domain annotation from OMA Browser based on a URL or a
#' raw data frame contains annotation info from OMA
#' @param domainURL URL address for domain annotation of ONE OMA id
#' @return data frame contains feature names, start and end positions
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' getOmaDomainFromURL("https://omabrowser.org/api/protein/7916808/domains/")

getOmaDomainFromURL <- function(domainURL) {
    if (grepl("https://", domainURL[1])) {
        domains <- OmaDB::resolveURL(domainURL)$regions
    } else {
        domains <- domainURL$regions
    }

    pos <- strsplit(domains$location, ":")
    allPos <- unlist(pos)
    start <- allPos[c(TRUE, FALSE)]
    end <- allPos[c(FALSE, TRUE)]

    feature <- ifelse(
        nchar(domains$name) > 0,
        paste0(domains$source, "_", domains$domainid,
               " (",domains$name,")"),
        paste0(domains$source, "_", domains$domainid)
    )
    feature <- gsub("#", "-", feature)
    featureRep <- vapply(pos, FUN.VALUE = numeric(1), function (x) length(x) / 2)

    return(data.frame(
        feature = rep(feature, featureRep),
        start,
        end,
        stringsAsFactors = FALSE
    ))
}

#' Get taxonomy ID, sequence and annotation for one OMA sequence
#' @export
#' @param id oma ID of a ortholog
#' @return data frame
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' getOmaDataForOneOrtholog("HUMAN29397")

getOmaDataForOneOrtholog <- function(id) {
    # get ncbi taxonomy id
    specName <- substr(id, 1, 5)
    taxonID <- paste0(
        "ncbi", OmaDB::getTaxonomy(members = specName, newick = FALSE)$id
    )
    # get raw data
    raw <- OmaDB::getProtein(id)

    # get sequence
    seq <- as.character(raw$sequence)

    # get sequence length
    length <- raw$sequence_length

    # get annotation
    rawDomains <- suppressWarnings(raw$domains)
    domainDf <- suppressWarnings(getOmaDomainFromURL(rawDomains))
    domainDfJoin <- c(domainDf, sep = "#")
    domains <- paste(unlist(do.call(paste, domainDfJoin)), collapse = "\t")

    # return data frame contains all info
    omaDf <- data.frame(
        orthoID = id,
        taxonID = taxonID,
        seq = seq,
        length = length,
        domains = domains,
        stringsAsFactors = FALSE
    )
    return(omaDf)
}

#' Get OMA info for a query protein and its orthologs
#' @description Get taxonomy IDs, sequences and annotations for an OMA
#' orthologous group (or OMA HOG).
#' @export
#' @param seedID protein query ID in OMA or UniProt format
#' @param orthoType type of OMA orthologs
#' @return data frame contains info for all sequences of that OMA group (or HOG)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' getDataForOneOma("HUMAN29397", "OG")

getDataForOneOma <- function(seedID, orthoType){
    # get members
    idList <- getOmaMembers(seedID, orthoType)
    specName <- substr(idList, 1, 5)

    # get taxonomy IDs
    taxonID <- paste0(
        "ncbi",
        lapply(
            specName,
            function (x) OmaDB::getTaxonomy(members = x, newick = FALSE)$id
        )
    )

    # get sequences, protein lengths and their domain annotations
    raw <- OmaDB::getProtein(idList)

    seq <- lapply(
        seq(length(raw)),
        function (x) as.character(raw[[x]]$sequence)
    )

    length <- lapply(
        seq(length(raw)),
        function (x) raw[[x]]$sequence_length
    )

    domains <- lapply(
        seq(length(raw)),
        function (x) {
            rawDomains <- suppressWarnings(raw[[x]]$domains)
            domainDf <- suppressWarnings(getOmaDomainFromURL(rawDomains))
            domainDfJoin <- c(domainDf, sep = "#")
            paste(unlist(do.call(paste, domainDfJoin)), collapse = "\t")
        }
    )

    # return data frame
    return(data.frame(
        orthoID = idList,
        taxonID = taxonID,
        seq = unlist(seq),
        length = unlist(length),
        domains = unlist(domains),
        seed = seedID,
        stringsAsFactors = FALSE
    ))
}

#' Create phylogenetic profile from a raw OMA dataframe
#' @export
#' @param finalOmaDf raw OMA data for a list of input IDs
#' @return phylogenetic profiles in long format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' omaData <- getDataForOneOma("HUMAN29397", "OG")
#' createProfileFromOma(omaData)
#' @seealso \code{\link{getDataForOneOma}}

createProfileFromOma <- function(finalOmaDf) {
    profileDf <- finalOmaDf[, c("seed", "taxonID", "orthoID")]
    colnames(profileDf) <- c("geneID", "ncbiID", "orthoID")
    return(profileDf[!duplicated(profileDf), ])
}

#' Create domain annotation dataframe from a raw OMA dataframe
#' @export
#' @param finalOmaDf raw OMA data for a list of input IDs
#' @return domain annotation in a dataframe to input into PhyloProfile, which
#' contains seed ID, ortholog ID, ortholog length, annotated feature, start
#' and end position of that feature.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' omaData <- getDataForOneOma("HUMAN29397", "OG")
#' getAllDomainsOma(omaData)
#' @seealso \code{\link{getDataForOneOma}}

getAllDomainsOma <- function(finalOmaDf) {
    seedID <- finalOmaDf$seed
    orthoID <- finalOmaDf$orthoID
    length <- finalOmaDf$length
    feature <- finalOmaDf$domains

    featureCount <- sapply(
        feature,
        function (x) length(unlist(strsplit(x, "\t")))
    )

    domainDf <- as.data.frame(
        stringr::str_split_fixed((unlist(strsplit(feature, "\t"))), "#", 3)
    )
    colnames(domainDf) <- c("feature", "start", "end")
    domainDf$seedID <- rep(paste0(seedID, "#", orthoID), featureCount)
    domainDf$orthoID <- rep(orthoID, featureCount)
    domainDf$length <- rep(length, featureCount)

    seedDf <- domainDf[domainDf$orthoID %in% levels(as.factor(seedID)),]
    seedDfFull <- data.frame(
        feature = rep(seedDf$feature, nlevels(as.factor(domainDf$seedID))),
        start = rep(seedDf$start, nlevels(as.factor(domainDf$seedID))),
        end = rep(seedDf$end, nlevels(as.factor(domainDf$seedID))),
        seedID = rep(levels(as.factor(domainDf$seedID)), each = nrow(seedDf)),
        orthoID = rep(seedDf$orthoID, nlevels(as.factor(domainDf$seedID))),
        length = rep(seedDf$length, nlevels(as.factor(domainDf$seedID))),
        stringsAsFactors = FALSE
    )

    return(
        rbind(
            domainDf, seedDfFull
        )[, c("seedID", "orthoID", "length", "feature", "start", "end")]
    )
}

#' Get all fasta sequences from a raw OMA dataframe
#' @export
#' @param finalOmaDf raw OMA data for a list of input IDs
#' @return A list contains all fasta sequences
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' omaData <- getDataForOneOma("HUMAN29397", "OG")
#' getAllFastaOma(omaData)
#' @seealso \code{\link{getDataForOneOma}}

getAllFastaOma <- function(finalOmaDf) {
    fastaDf <- finalOmaDf[, c("orthoID", "seq")]
    fastaOut <- paste(paste0(">", fastaDf$orthoID), fastaDf$seq, sep = "\n")
    return(unique(fastaOut))
}

#' Get selected fasta sequences from a raw OMA dataframe
#' @export
#' @param finalOmaDf raw OMA data for a list of input IDs
#' @param seqID sequence need to be returned
#' @return required sequence in fasta format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' omaData <- getDataForOneOma("HUMAN29397", "OG")
#' getSelectedFastaOma(omaData, "HUMAN29397")
#' @seealso \code{\link{getDataForOneOma}}

getSelectedFastaOma <- function(finalOmaDf, seqID) {
    selectedDf <- subset(
        finalOmaDf[, c("orthoID", "seq")],
        finalOmaDf$orthoID == seqID
    )
    header <- paste0(">", selectedDf$orthoID[1])
    return(paste(header, selectedDf$seq[1], sep = "\n"))
}
