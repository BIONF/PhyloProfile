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
    invalid <- list()
    for (id in ids) {
        id <- as.character(id)
        data <- OmaDB::getProtein("protein", id)
        if (is.null(data$entry_nr)) invalid <- c(invalid, id)
    }
    return(invalid)
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
            OmaDB::getOMAGroup(type = "group", id = id)$members$omaid
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
    
    domains$feature <- NA
    domains$start <- NA
    domains$end <- NA
    for (i in seq_len(nrow(domains))) {
        pos <- unlist(strsplit(domains$location[i], ":"))
        domains[i,]$start <- pos[1]
        domains[i,]$end <- pos[2]
        
        if (nchar(domains$name[i]) > 0) {
            domains[i,]$feature <- paste0(
                domains$source[i],"_",domains$domainid[i],
                " (",domains$name[i],")"
            )
        } else {
            domains[i,]$feature <- paste0(
                domains$source[i],"_",domains$domainid[i]
            )
        }
        domains[i,]$feature <- gsub("#", "-", domains[i,]$feature)
    }
    return(domains[, c("feature","start","end")])
}

#' Get taxonomy ID, sequence and annotation for one OMA sequence
#' @export
#' @param id oma ID of a ortholog
#' @return data frame
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' getOmaDataForOneOrtholog("HUMAN29397")

getOmaDataForOneOrtholog <- function(id) {
    omaDf <- data.frame(
        "orthoID" = character(),
        "taxonID" = character(),
        "seq" = character(),
        "length" = numeric(),
        "domains" = character(),
        stringsAsFactors = FALSE
    )
    
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
    omaDf[1,] <- c(id, taxonID, seq, length, domains)
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
    finalOmaDf <- data.frame()
    
    # get members
    members <- getOmaMembers(seedID, orthoType)
    omaSeedID <- OmaDB::getProtein(seedID)$omaid
    
    # get all data
    j <- 1
    for (ortho in members) {
        orthoDf <- getOmaDataForOneOrtholog(ortho)
        orthoDf$seed <- seedID
        if (ortho == omaSeedID) {
            orthoDf$orthoID <- seedID
        }
        finalOmaDf <- rbind(finalOmaDf, orthoDf)
        
        p <- j / length(members) * 100
        svMisc::progress(p)
        j <- j + 1
    }
    
    return(finalOmaDf)
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

#' Create domain annotation dataframe for one OMA protein
#' @param domainID protein domain ID
#' @param orthoID protein ID
#' @param length protein length
#' @param domainList list of all domains and their positions for this protein
#' @return domain annotation in a dataframe
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

createDomainDf <- function(domainID, orthoID, length, domainList) {
    domainDf = data.frame(
        seedID = character(),
        orthoID = character(),
        length = numeric(),
        feature = character(),
        start = numeric(),
        end = numeric(),
        stringsAsFactors = FALSE
    )
    
    for (i in seq_len(length(domainList))) {
        annoInfo <- strsplit(domainList[i], "#")[[1]]
        domainDf[i,] <- c( domainID, orthoID, length, annoInfo)
    }
    
    return(domainDf)
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
    omaDomainDf <- data.frame()
    for (i in seq_len(nrow(finalOmaDf))) {
        domainID <- paste0(
            finalOmaDf[i,]$seed, "#", finalOmaDf[i,]$orthoID
        )
        
        seedLine <-
            finalOmaDf[finalOmaDf$orthoID == finalOmaDf[i,]$seed, ]
        seedDomains <- strsplit(as.character(seedLine$domains), "\t")[[1]]
        seedDomainDf <- createDomainDf(
            domainID, seedLine$orthoID, seedLine$length, seedDomains
        )
        omaDomainDf <- rbind(omaDomainDf, seedDomainDf)
        
        orthoDomains <-
            strsplit(as.character(finalOmaDf[i,]$domains), "\t")[[1]]
        orthoDomainDf <- createDomainDf(
            domainID,
            finalOmaDf[i,]$orthoID,
            finalOmaDf[i,]$length,
            orthoDomains
        )
        omaDomainDf <- rbind(omaDomainDf, orthoDomainDf)
    }
    
    omaDomainDf$length <- as.numeric(omaDomainDf$length)
    omaDomainDf$start <- as.numeric(omaDomainDf$start)
    omaDomainDf$end <- as.numeric(omaDomainDf$end)
    return(omaDomainDf)
}

#' Get all fasta sequences from a raw OMA dataframe
#' @export
#' @param finalOmaDf raw OMA data for a list of input IDs
#' @return dataframe contains all fasta sequences
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' omaData <- getDataForOneOma("HUMAN29397", "OG")
#' getAllFastaOma(omaData)
#' @seealso \code{\link{getDataForOneOma}}

getAllFastaOma <- function(finalOmaDf) {
    omaFastaDf <- data.frame()
    
    fastaDf <- finalOmaDf[, c("orthoID", "seq")]
    for (i in seq_len(nrow(fastaDf))) {
        seqID <- as.character(fastaDf$orthoID[i])
        seq <- as.character(fastaDf$seq[i])
        fastaOut <- paste(paste0(">", seqID), seq, sep = "\n")
        omaFastaDf <- rbind(omaFastaDf, as.data.frame(fastaOut))
    }
    
    return(omaFastaDf[!duplicated(omaFastaDf), ])
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
