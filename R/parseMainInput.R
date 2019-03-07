# Parse main input functions
# These functions check the validity of the main input file
# and convert different input format into long format.

#' Check the validity of the main input file
#' @description Check if input file has one of the following format: orthoXML,
#' multiple FASTA, wide or long matrix, or a list of OMA IDs.
#' @export
#' @param filein input file
#' @return The format of the input file format, or type of error
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{checkOmaID}}
#' @examples
#' filein <- system.file(
#'     "extdata", "test.main.wide", package = "PhyloProfile", mustWork = TRUE
#' )
#' checkInputValidity(filein)

checkInputValidity <- function(filein) {
    inputDt <- as.data.frame(read.table(
        file = filein,
        sep = "\t",
        header = FALSE,
        check.names = FALSE,
        comment.char = "",
        fill = TRUE
    ))

    if (is.na(inputDt[1, ncol(inputDt)])) {
        return("moreCol")
    } else {
        names(inputDt) <- as.character(unlist(inputDt[1, ]))

        # XML format (starts with <?xml)
        if (grepl("<?xml", colnames(inputDt)[1])) {
            return("xml")
        }
        # FASTA format (starts with ">" )
        else if (grepl(">", colnames(inputDt)[1]) == TRUE) {
            return("fasta")
        }
        # LONG or WIDE format (starts with "geneID")
        else {
            if (grepl("geneID", colnames(inputDt)[1])) {
                # LONG format
                if (is.na(pmatch("ncbi", colnames(inputDt)[3])) ||
                    is.na(pmatch("ncbi", colnames(inputDt)[4])) ||
                    is.na(pmatch("ncbi", colnames(inputDt)[5]))) {
                    return("long")
                }
                # WIDE format
                else {
                    tmp <- inputDt[inputDt == ""][1]
                    if (!is.na(tmp) & tmp == "") {
                        return("emptyCell")
                    } else {
                        return("wide")
                    }
                }
            }
            # OMA ids
            else {
                invalidOma <- checkOmaID(levels(inputDt[,1]))
                if (length(invalidOma) == 0) {
                    return("oma")
                } else {
                    return(invalidOma)
                }
            }
        }
    }
}

#' Parse orthoXML input file
#' @param inputFile input file in xml format
#' @return A data frame containing input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' inputFile <- system.file(
#'     "extdata", "test.main.xml", package = "PhyloProfile", mustWork = TRUE
#' )
#' xmlParser(inputFile)
#' }

xmlParser <- function(inputFile){
    path <- paste(
        system.file(package="PhyloProfile"),
        "PhyloProfile/scripts/orthoxmlParser.py",
        sep="/"
    )
    cmd <- paste("python ", path, " -i ", inputFile, sep = "")

    dfIn <- as.data.frame(read.table(text = system(cmd, intern = TRUE)))

    # the first row will be the header
    colnames(dfIn) <- as.character(unlist(dfIn[1, ]))

    dfIn <- subset(dfIn[dfIn$geneID != "geneID", ])
    dfIn <- droplevels(dfIn)

    return(dfIn)
}

#' Parse multi-fasta input file
#' @param inputFile input multiple fasta file
#' @return A data frame containing input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' inputFile <- system.file(
#'     "extdata", "test.main.fasta", package = "PhyloProfile", mustWork = TRUE
#' )
#' fastaParser(inputFile)
#' }

fastaParser <- function(inputFile){
    path <- paste(
        system.file(package="PhyloProfile"),
        "PhyloProfile/scripts/fastaParser.py",
        sep="/"
    )
    cmd <- paste("python ", path, " -i ", inputFile, sep = "")

    var1 <- NULL
    var2 <- NULL

    dfIn <- as.data.frame(read.table(text = system(cmd, intern = TRUE)))

    # the first row will be the header
    colnames(dfIn) <- as.character(unlist(dfIn[1, ]))
    dfIn <- subset(dfIn[dfIn$geneID != "geneID", ])
    dfIn <- droplevels(dfIn)

    # remove var1 and var2 columns if they are all NAs
    if (all(is.na(dfIn$var2))) {
        dfIn <- subset(dfIn, select = -c(var2) )
    }
    if (all(is.na(dfIn$var1))) {
        dfIn <- subset(dfIn, select = -c(var1) )
    }

    return(dfIn)
}

#' Transform input file in wide matrix into long matrix format
#' @param inputFile input file in wide matrix format
#' @return A data frame containing input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' \dontrun{
#' inputFile <- system.file(
#'     "extdata", "test.main.wide", package = "PhyloProfile", mustWork = TRUE
#' )
#' wideToLong(inputFile)
#' }

wideToLong <- function(inputFile){
    wideDataframe <- as.data.frame(read.table(
        file = inputFile,
        sep = "\t",
        header = TRUE,
        check.names = FALSE,
        comment.char = ""
    ))
    longDataframe <- data.frame()
    rowNrLong <- 0
    ncbiIDs <- colnames(wideDataframe)

    for (rowNr in seq_len(nrow(wideDataframe))) {
        geneID <- wideDataframe[rowNr, 1]
        for (columnNr in 2:ncol(wideDataframe)) {
            currentCell <- as.character(wideDataframe[rowNr, columnNr])
            newRowInfo <- unlist(strsplit(currentCell, "#"))
            rowNrLong <- rowNrLong + 1
            longDataframe[rowNrLong, 1] <- geneID
            longDataframe[rowNrLong, 2] <- ncbiIDs[columnNr]
            longDataframe[rowNrLong, 3] <- newRowInfo[1]
            longDataframe[rowNrLong, 4] <- newRowInfo[2]
            longDataframe[rowNrLong, 5] <- newRowInfo[3]
        }
    }

    colnames(longDataframe) <- c("geneID", "ncbiID", "orthoID", "var1", "var2")
    return(longDataframe)
}

#' Create a long matrix format for all kinds of input file
#' @export
#' @param inputFile input file in orthoXML, multiple FASTA, wide or long
#' matrix format.
#' @return A data frame that contains input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{xmlParser}}, \code{\link{fastaParser}},
#' \code{\link{wideToLong}}
#' @examples
#' inputFile <- system.file(
#'     "extdata", "test.main.wide", package = "PhyloProfile", mustWork = TRUE
#' )
#' createLongMatrix(inputFile)

createLongMatrix <- function(inputFile){
    if (inputFile[1] == "lca-micros") {
        inputURL <- paste0(
            "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/",
            "demo/test.main.long"
        )
        longDataframe <-
            read.table(
                inputURL,
                sep = "\t",
                header = TRUE,
                fill = TRUE,
                stringsAsFactors = FALSE
            )
    } else if (inputFile[1] == "ampk-tor") {
        inputURL <- paste0(
            "https://raw.githubusercontent.com/BIONF/phyloprofile-data/master/",
            "expTestData/ampk-tor/ampk-tor.phyloprofile"
        )
        longDataframe <-
            read.table(
                inputURL,
                sep = "\t",
                header = TRUE,
                fill = TRUE,
                stringsAsFactors = FALSE
            )
    } else {
        filein <- inputFile
        if (is.null(filein)) return()
        inputType <- checkInputValidity(filein)

        # XML
        if (inputType == "xml") {
            longDataframe <- xmlParser(filein)
        }
        # FASTA
        else if (inputType == "fasta") {
            longDataframe <- fastaParser(filein)
        }
        # LONG
        else if (inputType == "long") {
            longDataframe <- as.data.frame(read.table(
                file = filein,
                sep = "\t",
                header = TRUE,
                check.names = FALSE,
                comment.char = ""
            ))
        }
        # WIDE
        else if (inputType == "wide") {
            longDataframe <- wideToLong(filein)
        }
        else {
            return(NULL)
        }
    }

    # make sure all columns have the same type (factor)
    for (i in seq_len(ncol(longDataframe))) {
        longDataframe[, i] <- as.factor(longDataframe[, i])
    }
    
    # remove duplicated lines
    longDataframe <- longDataframe[!duplicated(longDataframe),]

    # longDataframe$orthoID <- gsub("\\|",":",longDataframe$orthoID)
    return(longDataframe)
}
