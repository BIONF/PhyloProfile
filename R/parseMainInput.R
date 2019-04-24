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
    inputDt <- read.table(
        file = filein,
        sep = "\t",
        header = FALSE,
        check.names = FALSE,
        comment.char = "",
        fill = TRUE
    )

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
#' @export
#' @param inputFile input file in xml format
#' @return A data frame containing input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @importFrom xml2 read_xml
#' @importFrom xml2 xml_find_all
#' @importFrom xml2 xml_attr
#' @examples
#' inputFile <- system.file(
#'     "extdata", "test.main.xml", package = "PhyloProfile", mustWork = TRUE
#' )
#' xmlParser(inputFile)

xmlParser <- function(inputFile){
    data <- read_xml(inputFile)
    
    # get all genes for each taxon
    species <- xml_find_all(data, ".//species")
    
    speciesDf <- data.frame(
        "ncbiID" = character(),
        "refGeneID" = character(),
        "orthoID" = character(),
        stringsAsFactors=FALSE
    )
    c <- 1
    
    for (i in seq_len(length(species))) {
        speciesName <- xml_attr(species[i], "name")
        speciesID <- xml_attr(species[i], "NCBITaxId")
        
        gene <- xml_find_all(species[i], ".//gene")
        for (j in seq_len(length(gene))) {
            if (grepl("protId", gene[j])) {
                refGeneID <- xml_attr(gene[j], "id")
                orthoIDTmp <- unlist(
                    strsplit(
                        xml_attr(gene[j], "protId"), split = '|', fixed = TRUE
                    )
                )
                orthoID <- orthoIDTmp[1]
            } else if (grepl("geneId", gene[j])) {
                refGeneID <- xml_attr(gene[j], "id")
                orthoID <- xml_attr(gene[j], "protId")
            }
            speciesDf[c, ] <- c(paste0("ncbi", speciesID), refGeneID, orthoID)
            c <- c + 1
        }
    }

    # get orthologs and their scores
    orthoDf <- data.frame(
        "geneID" = character(),
        "refGeneID" = character(),
        stringsAsFactors=FALSE
    )
    orthoDf[1,] <- c(NA,NA)
    
    scoreType <- xml_find_all(data, ".//scoreDef")
    for(s in seq_len(length(scoreType))){
        scoreID <- xml_attr(scoreType[s], "id")
        orthoDf[scoreID] <- NA
    }
    
    c <- 1
    orthoGroup <- xml_find_all(data, ".//orthologGroup")
    for (i in seq_len(length(orthoGroup))) {
        groupID <- xml_attr(orthoGroup[i], "id")
        gene <- xml_find_all(orthoGroup[i], ".//geneRef")
        for (j in seq_len(length(gene))) {
            refGeneID <- xml_attr(gene[j], "id")
            orthoDf[c, "geneID"] <- groupID
            orthoDf[c, "refGeneID"] <- refGeneID
            scores <- list()
            score <- xml_find_all(gene[j], ".//score")
            for (k in seq_len(length(score))) {
                id <- xml_attr(score[k], "id")
                value <- xml_attr(score[k], "value")
                orthoDf[c, id] <- value
            }
            c <- c + 1
        }
    }
    
    # merge into final dataframe
    finalDf <- merge(speciesDf, orthoDf, all.y = TRUE, by = "refGeneID")
    
    # remove refGeneID column and reorder columns
    finalDf <- finalDf[, -1]
    refcols <- c("geneID", "ncbiID", "orthoID")
    finalDf <- finalDf[, c(refcols, setdiff(names(finalDf), refcols))]
    
    # remove columns that contains only NA
    finalDf <- finalDf[, colSums(is.na(finalDf)) < nrow(finalDf)]
    
    # return final long dataframe
    return(finalDf)
}

#' Parse multi-fasta input file
#' @export
#' @param inputFile input multiple fasta file
#' @return A data frame containing input data in long-format
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @importFrom Biostrings readAAStringSet
#' @examples
#' inputFile <- system.file(
#'     "extdata", "test.main.fasta", package = "PhyloProfile", mustWork = TRUE
#' )
#' fastaParser(inputFile)

fastaParser <- function(inputFile){
    # split sequence IDs into columns
    fastaFile <- Biostrings::readAAStringSet(inputFile)
    seqID <- names(fastaFile)
    faDf <- data.frame(seqID)

    finalDf <- data.frame(
        do.call('rbind', strsplit(as.character(faDf$seqID),'|',fixed=TRUE))
    )
    
    # rename columns
    colnames(finalDf) <- c("geneID", "ncbiID", "orthoID")
    
    cidx <- seq_len(ncol(finalDf) - 3)
    colnames(finalDf)[cidx + 3] <- paste0("value_", cidx)
    
    return(finalDf)
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
    wideDataframe <- read.table(
        file = inputFile,
        sep = "\t",
        header = TRUE,
        check.names = FALSE,
        comment.char = ""
    )
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
            longDataframe <- read.table(
                file = filein,
                sep = "\t",
                header = TRUE,
                check.names = FALSE,
                comment.char = "",
                stringsAsFactors = FALSE
            )
        }
        # WIDE
        else if (inputType == "wide") {
            longDataframe <- wideToLong(filein)
        }
        else {
            return(NULL)
        }
    }

    # convert geneID, ncbiID and orthoID into factor and var1, var2 into numeric
    for (i in seq_len(3)) {
        longDataframe[, i] <- as.factor(longDataframe[, i])
    }
    if (ncol(longDataframe) > 3) {
        for (j in seq(4, ncol(longDataframe))){
            longDataframe[,j] <- suppressWarnings(as.numeric(longDataframe[,j]))
        }
    }
    
    # remove duplicated lines
    longDataframe <- longDataframe[!duplicated(longDataframe),]

    # longDataframe$orthoID <- gsub("\\|",":",longDataframe$orthoID)
    return(longDataframe)
}
