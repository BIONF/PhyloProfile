#' Parse domain input file
#' @export
#' @param seed seed ID(s)
#' @param inputFile name of input data (demo data, file name or path to folder)
#' @param type type of data (demo, file or folder)
#' @importFrom utils read.csv
#' @return A dataframe containing protein domain info
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{getDomainOnline}}, \code{\link{getDomainFile}},
#' \code{\link{getDomainFolder}}
#' @examples
#' seed <- "OG_1009"
#' inputFile <- system.file(
#'     "extdata", "domainFiles/OG_1009.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' type <- "file"
#' parseDomainInput(seed, inputFile, type)

parseDomainInput <- function(seed, inputFile, type) {
    domains <- data.frame()
    file <- NULL

    # get domain file(s) from online data
    if (type == "demo") {
        file <- getDomainOnline(inputFile)
    }
    # or from single file
    else if (type == "file") {
        file <- getDomainFile(inputFile)
    }
    # or from a domain folder
    else {
        file <- getDomainFolder(seed, inputFile)
        if (file == "noSelectHit") return("noSelectHit")
        else if (file == "noFileInFolder") return("noFileInFolder")
    }

    # parse domain file
    # for demo data
    if (type == "demo") {
        domainDf <- as.data.frame(read.csv(file,
                                            sep = "\t",
                                            header = FALSE,
                                            comment.char = "",
                                            stringsAsFactors = FALSE,
                                            quote = ""))
        domains <- rbind(domains, domainDf)

        if (ncol(domains) == 5) {
            colnames(domains) <- c(
                "seedID", "orthoID", "feature", "start", "end"
            )
        } else if (ncol(domains) == 6) {
            colnames(domains) <- c(
                "seedID", "orthoID", "feature", "start", "end", "weight"
            )
        } else if (ncol(domains) == 7) {
            colnames(domains) <- c(
                "seedID", "orthoID", "feature", "start", "end", "weight", "path"
            )
        }
    }
    # for user input data
    else {
        if (file != FALSE) {
            exeptions <- c("noFileInput", "noSelectHit", "noFileInFolder")
            if (!(file %in% exeptions)) {
                domainDf <- as.data.frame(read.table(
                    file,
                    sep = "\t",
                    header = FALSE,
                    comment.char = ""
                ))
                domains <- rbind(domains, domainDf)
            }
        }

        if (ncol(domains) == 5) {
            colnames(domains) <- c(
                "seedID", "orthoID", "feature", "start", "end"
            )
        } else if (ncol(domains) == 6) {
            colnames(domains) <- c(
                "seedID", "orthoID", "length", "feature", "start", "end"
            )
        } else if (ncol(domains) == 7) {
            colnames(domains) <- c(
                "seedID", "orthoID", "length", "feature", "start", "end",
                "weight"
            )
        } else if (ncol(domains) == 8) {
            colnames(domains) <- c(
                "seedID", "orthoID", "length", "feature", "start", "end",
                "weight", "path"
            )
        } else {
            return("ERR")
        }
    }

    if (nrow(domains) == 0) return("ERR-0")

    domains$seedID <- as.character(domains$seedID)
    domains$orthoID <- as.character(domains$orthoID)
    domains$seedID <- gsub("\\|",":",domains$seedID)
    domains$orthoID <- gsub("\\|",":",domains$orthoID)

    return(domains)
}

#' Get domain file(s) for online data set
#' @param demoData demo data name (either lca-micros or ampk-tor)
#' @return Domain file and its complete directory path
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' getDomainOnline("lca-micors")
#' }

getDomainOnline <- function(demoData) {
    if (demoData == "lca-micros") {
        fileDomain <- {
            suppressWarnings(
                paste0(
                    "https://github.com/BIONF/phyloprofile-data/blob/master/",
                    "demo/domainFiles/concatenate.domains?raw=true"
                )
            )
        }
    } else {
        fileDomain <- {
            suppressWarnings(
                paste0(
                    "https://raw.githubusercontent.com/BIONF/phyloprofile-data",
                    "/master/expTestData/ampk-tor/ampk-tor.domains_F?raw=true"
                )
            )
        }
    }

    return(fileDomain)
}

#' Get domain file from a single (concatenate) file
#' @param inputFile concatenate domain file
#' @return Domain file and its complete directory path
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' inputFile <- system.file(
#'     "extdata", "domainFiles/OG_1009.domains",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' getDomainFile("lca-micors")
#' }

getDomainFile <- function(inputFile) {
    fileDomain <- inputFile
    return(fileDomain)
}

#' Get domain files from a folder
#' @param seed seed ID
#' @param domainPath path to domain folder
#' @return Domain file and its complete directory path
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' domainPath <- paste0(
#'     path.package("PhyloProfile", quiet = FALSE), "extdata/domainFiles"
#' )
#' getDomainOnline("OG_1009", domainPath)
#' }

getDomainFolder <- function(seed, domainPath){
    if (is.null(seed)) {
        fileDomain <- "noSelectHit"
    } else {
        # check file extension
        allExtension <- c("txt", "csv", "list", "domains", "architecture")
        flag <- 0

        for (i in seq_len(length(allExtension))) {
            fileDomain <- paste0(
                domainPath, "/", seed, ".", allExtension[i]
            )
            if (file.exists(fileDomain) == TRUE) {
                flag <- 1
                break()
            }
        }

        if (flag == 0) {
            fileDomain <- "noFileInFolder"
        }
    }

    return(fileDomain)
}
