#' Get fasta sequences from main input file in multi-fasta format
#' @export
#' @param seqIDs list of sequences IDs. Set seqIDs = "all" if you want to get
#' all fasta sequences from the input file.
#' @param file raw phylogenetic profile input file in multi-fasta format.
#' @return A dataframe with one column contains sequences in fasta format.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' file <- system.file(
#'     "extdata", "test.main.fasta",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' getFastaFromFasInput("all", file)

getFastaFromFasInput <- function(seqIDs = NULL, file = NULL) {
    if (is.null(seqIDs)) stop("No sequence ID given!")
    if (is.null(seqIDs) | is.null(file)) 
        stop("Sequence ID and input fasta file cannot be NULL!")
    faFile <- Biostrings::readAAStringSet(file)

    if (length(seqIDs) == 1 & seqIDs[1] == "all") seqIDs <- names(faFile)
    return(data.frame(
        fasta = paste(
            paste0(">", seqIDs),
            lapply(
                pmatch(seqIDs, names(faFile)),
                function (x) as.character(faFile[[x]])
            ),
            sep = "\n"
        ), stringsAsFactors = FALSE
    ))
}

#' Get fasta sequences from main input file in multi-fasta format
#' @export
#' @param seqIDs list of sequences IDs. Set seqIDs = "all" if you want to get
#' all fasta sequences from the concatenated input fasta file.
#' @param concatFasta input concatenated fasta file.
#' @return A dataframe with one column contains sequences in fasta format.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @examples
#' concatFasta <- system.file(
#'     "extdata", "fastaFiles/concatenatedFile.fa",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' getFastaFromFasInput("all", concatFasta)

getFastaFromFile <- function(seqIDs = NULL, concatFasta = NULL) {
    if (is.null(seqIDs)) stop("No sequence ID given!")
    if (!is.null(concatFasta)) {
        # read fasta file and save sequences into dataframe
        faFile <- Biostrings::readAAStringSet(toString(concatFasta))

        if (length(seqIDs) == 1 & seqIDs[1] == "all") seqIDs <- names(faFile)
        if (
            length(unique(pmatch(seqIDs, names(faFile)))) == 1 &
            is.na(unique(pmatch(seqIDs, names(faFile)))[1])
        ) {
            return(data.frame(
                fasta = paste0("Please check FASTA header format!"),
                stringsAsFactors = FALSE
            ))
        } else {
            return(data.frame(
                fasta = paste(
                    paste0(">", seqIDs),
                    lapply(
                        pmatch(seqIDs, names(faFile)),
                        function (x) as.character(faFile[[x]])
                    ),
                    sep = "\n"
                ), stringsAsFactors = FALSE
            ))
        }
    } else {
        return(data.frame(
            fasta = paste0("Please provide FASTA file(s)!"),
            stringsAsFactors = FALSE
        ))
    }
}

#' Get fasta sequences
#' @description Get fasta sequences for the input phylogenetic profiles.
#' @export
#' @usage getFastaFromFolder(seqIDs = NULL, path = NULL, dirFormat = NULL,
#'     fileExt = NULL, idFormat = NULL)
#' @param seqIDs list of sequences IDs.
#' @param path path to fasta folder.
#' @param dirFormat directory format (either 1 for "path/speciesID.fa*" or 2 for
#' "path/speciesID/speciesID.fa*")
#' @param fileExt fasta file extension ("fa", "fasta", "fas" or "txt")
#' @param idFormat fasta header format (1 for ">speciesID:seqID", 2 for
#' ">speciesID@seqID", 3 for ">speciesID|seqID" or 4 for "seqID")
#' @return A dataframe with one column contains sequences in fasta format.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{mainLongRaw}}
#' @examples
#' seqIDs <- "RAT@10116@1|D3ZUE4"
#' path <- system.file(
#'     "extdata", "fastaFiles", package = "PhyloProfile", mustWork = TRUE
#' )
#' dirFormat <- 1
#' fileExt <- "fa"
#' idFormat <- 3
#' getFastaFromFolder(seqIDs, path, dirFormat, fileExt, idFormat)

getFastaFromFolder <- function(
    seqIDs = NULL, path = NULL, dirFormat = NULL, fileExt = NULL, idFormat=NULL
) {
    if (is.null(seqIDs)) stop("No sequence ID given!")
    if (is.null(path)) return(
        data.frame(
            fasta = paste0("Please provide path to FASTA file(s)!"),
            stringsAsFactors = FALSE))
    if (path != "") {
        # get list of species IDs
        if (idFormat == 1) {
            specID <- unlist(strsplit(seqIDs, "\\:"))[1]
        } else if (idFormat == 2) {
            specID <- unlist(strsplit(seqIDs, "@"))[1]
        } else if (idFormat == 3) {
            specID <- unlist(strsplit(seqIDs, "\\|"))[1]
        }
        # read all specices FASTA files at once
        if (idFormat == 4) {
            fileList <- list.files(path, pattern = fileExt)
            if(length(fileList) == 0) {
                return(data.frame(
                    fasta = paste0(
                        "No FASTA files with ", fileExt, "found in ", path
                    ), stringsAsFactors = FALSE))
            } else {
                fileWithPath <- paste0(path, "/", fileList)
                faFile <- Biostrings::readAAStringSet(fileWithPath)
            }
        } else {
            specID <- as.character(levels(as.factor(specID)))
            specID <- specID[specID != ""]
            file <- paste0(path, "/", specID, ".", fileExt)
            if (dirFormat == 2)
                file <- paste0(path, "/", specID, "/", specID, ".", fileExt)
            file <- file[file.exists(file)]
            if (length(file) == 0) {
                return(data.frame(
                    fasta = paste0(
                        "No FASTA files found in ", path, "! Please check ",
                        "settings for file extension and ID format again."
                    ), stringsAsFactors = FALSE
                ))
            } else faFile <- Biostrings::readAAStringSet(file)
        }
        # now get selected sequences
        if (length(unique(pmatch(seqIDs, names(faFile)))) == 1 
            & is.na(unique(pmatch(seqIDs, names(faFile)))[1])) {
            return(data.frame(
                fasta = paste0("Please check FASTA header format!"),
                stringsAsFactors = FALSE))
        } else {
            return(data.frame(
                fasta = paste(
                    paste0(">", seqIDs),
                    lapply(
                        pmatch(seqIDs, names(faFile)),
                        function (x) as.character(faFile[[x]])
                    ), sep = "\n"), stringsAsFactors = FALSE))
        }
    } else {
        return(data.frame(
            fasta = paste0("Please provide path to FASTA file(s)!"),
            stringsAsFactors = FALSE))
    }
}
