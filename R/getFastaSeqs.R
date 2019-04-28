#' Get fasta sequences
#' @export
#' @usage getFastaSeqs(dataIn, filein, demoData, inputTypeFasta,
#'     concatFasta, path, dirFormat, fileExt, idFormat)
#' @param dataIn a dataframe of the input phylogenetic profiles, that contains
#' at least 3 columns: geneIDs, orthoIDs and ncbiIDs
#' @param filein main input file (for checking input type)
#' @param demoData name of demo data set (either "ampk-tor", "lca-micros",
#' or "none")
#' @param inputTypeFasta source of fasta sequences ("Concatenated fasta file"
#' or "Fasta folder")
#' @param concatFasta input concatenated fasta file
#' @param path path to fasta folder
#' @param dirFormat directory format (either 1 for "path/speciesID.fa*" or 2 for
#' "path/speciesID/speciesID.fa*")
#' @param fileExt fasta file extension ("fa", "fasta", "fas" or "txt")
#' @param idFormat fasta header format (1 for ">speciesID:seqID", 2 for
#' ">speciesID@seqID", 3 for ">speciesID|seqID" or 4 for "seqID")
#' @return A dataframe contains fasta sequences
#' @importFrom IRanges reverse
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{checkInputValidity}}, \code{\link{mainLongRaw}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' dataIn <- mainLongRaw
#' filein <- system.file(
#'     "extdata", "test.main.long", package = "PhyloProfile", mustWork = TRUE
#' )
#' demoData <- "none"
#' inputTypeFasta <- "Concatenated fasta file"
#' concatFasta <- system.file(
#'     "extdata", "fastaFiles/concatenatedFile.fa",
#'     package = "PhyloProfile", mustWork = TRUE
#' )
#' # the following parameters are conly required if upload a fasta folder!
#' path <- NULL
#' dirFormat <- NULL
#' fileExt <- NULL
#' idFormat <- NULL
#' # get fasta sequences
#' getFastaSeqs(
#'     dataIn, filein, demoData,
#'     inputTypeFasta,
#'     concatFasta,
#'     path,
#'     dirFormat,
#'     fileExt,
#'     idFormat
#' )

getFastaSeqs <- function(
    dataIn, filein, demoData,
    inputTypeFasta,
    concatFasta,
    path,
    dirFormat,
    fileExt,
    idFormat
){
    # check main input
    if (!is.null(filein)) {
        inputType <- checkInputValidity(filein)
    } else{
        inputType <- "NA"
    }

    # get seqs for AMPK-TOR and microsporidia ONLINE demo data -----------------
    if (demoData == "ampk-tor" | demoData == "lca-micros") {
        fastaURL <-
            paste0(
                "https://raw.githubusercontent.com/BIONF/phyloprofile-data/",
                "master/demo/fastaFile/concatenatedSeq.fa"
            )
        if (demoData == "ampk-tor") {
            fastaURL <-
                paste0(
                    "https://raw.githubusercontent.com/BIONF/",
                    "phyloprofile-data/master/expTestData/ampk-tor/",
                    "ampk-tor.extended.fa"
                )
        }

        if (RCurl::url.exists(fastaURL)) {
            # load fasta file
            faFile <- Biostrings::readAAStringSet(fastaURL)
            
            # get sequences
            seqIDs <- as.character(dataIn$orthoID)
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
        } else {
            return(data.frame(
                fasta = paste0(fastaURL, " not found!!!"),
                stringsAsFactors = FALSE
            ))
        }
    }

    # get seqs for fasta main input --------------------------------------------
    if (inputType == "fasta") {
        file <- filein
        faFile <- Biostrings::readAAStringSet(file)
        
        seqIDs <- paste0(
            as.character(dataIn$geneID), "|ncbi",
            as.character(dataIn$ncbiID), "|",
            as.character(dataIn$orthoID)
        )
        
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

    # get seqs for main input in other formats ---------------------------------
    else {
        # * get seqs from concatenated fasta file ------------------------------
        if (demoData == "none" &
            inputTypeFasta == "Concatenated fasta file") {
            if (!is.null(concatFasta)) {
                fasIn <- concatFasta
                file <- toString(fasIn)

                # read fasta file and save sequences into dataframe
                faFile <- Biostrings::readAAStringSet(file)

                seqIDs <- as.character(dataIn$orthoID)
                if (
                    length(unique(pmatch(seqIDs, names(faFile)))) == 1 & 
                    is.na(unique(pmatch(seqIDs, names(faFile)))[1])
                ) {
                    seqIDs <- paste0(
                        as.character(dataIn$geneID), "|ncbi",
                        as.character(dataIn$ncbiID), "|",
                        as.character(dataIn$orthoID)
                    )
                } 
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
                if (inputType != "fasta") {
                    return(data.frame(
                        fasta = paste0(
                            "Please provide FASTA file(s) in ",
                            "Input & settings page!"
                        ),
                        stringsAsFactors = FALSE
                    ))
                }
            }
        }

        # * get seqs from an folder --------------------------------------------
        if (demoData == "none"
            & inputTypeFasta == "Fasta folder"
            & inputType != "fasta") {
            if (path != "") {
                # get list of species IDs
                if (idFormat == 1) {
                    specDf <-
                        as.data.frame(
                            stringr::str_split_fixed(
                                reverse(as.character(dataIn$orthoID)),
                                ":", 2
                            )
                        )
                    specDf$specID <- reverse(as.character(specDf$V2))
                } else if (idFormat == 2) {
                    specDf <-
                        as.data.frame(
                            stringr::str_split_fixed(
                                reverse(as.character(dataIn$orthoID)),
                                "@", 2
                            )
                        )
                    specDf$specID <- reverse(as.character(specDf$V2))
                } else if (idFormat == 3) {
                    specDf <-
                        as.data.frame(
                            stringr::str_split_fixed(
                                reverse(as.character(dataIn$orthoID)),
                                "|", 2
                            )
                        )
                    specDf$specID <- reverse(as.character(specDf$V2))
                }

                # read all specices FASTA files at once
                if (idFormat == 4) {
                    fileList <- list.files(path, pattern = fileExt)
                    if(length(fileList) == 0) {
                        return(data.frame(
                            fasta = paste0("No FASTA files with ", fileExt, 
                                           "found in ", path),
                            stringsAsFactors = FALSE
                        ))
                    } else {
                        fileWithPath <- paste0(path, "/", fileList)
                        faFile <- Biostrings::readAAStringSet(fileWithPath)
                    }
                } else {
                    specID <- as.character(
                        levels(as.factor(specDf$specID))
                    )
                    specID <- specID[specID != ""]
                    
                    file <- paste0(path, "/", specID, ".", fileExt)
                    if (dirFormat == 2) {
                        file <- paste0(
                            path, "/", specID, "/", specID, ".", fileExt
                        )
                    }
                    file <- file[file.exists(file)]
                    if(length(file) == 0) {
                        return(data.frame(
                            fasta = paste0("No FASTA files found in ", path,
                                           "! Please check settings for file",
                                           " extension and ID format again."),
                            stringsAsFactors = FALSE
                        ))
                    } else {
                        faFile <- Biostrings::readAAStringSet(file)
                    }
                }

                # now get selected sequences
                seqIDs <- as.character(dataIn$orthoID)
                
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
                    fasta = paste0(
                        "Please provide FASTA file(s) in Input & settings page!"
                    ),
                    stringsAsFactors = FALSE
                ))
            }
        }
    }
}
