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
#' @param dirFormat directory format (either "path/speciesID.fa*" or
#' "path/speciesID/speciesID.fa*")
#' @param fileExt fasta file extension ("fa", "fasta", "fas" or "txt")
#' @param idFormat fasta header format (">speciesID:seqID",
#' ">speciesID@seqID", ">speciesID|seqID" or only "seqID")
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

    fastaOutDf <- data.frame()

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
            for (j in seq_len(nrow(dataIn))) {
                seqID <- as.character(dataIn$orthoID[j])
                groupID <- as.character(dataIn$geneID[j])

                # seq <- as.character(
                #     faDf$seq[faDf$seqID == paste0(">", seqID)]
                # )
                seq <- as.character(faFile[[match(seqID, names(faFile))]])
                fastaOut <- paste(paste0(">", seqID), seq, sep = "\n")
                fastaOutDf <- rbind(fastaOutDf, as.data.frame(fastaOut))
            }
        } else {
            fastaOut <- paste0(fastaURL, " not found!!!")
            fastaOutDf <- rbind(fastaOutDf, as.data.frame(fastaOut))
        }
    }

    # get seqs for fasta main input --------------------------------------------
    if (inputType == "fasta") {
        file <- filein
        fastaFile <- Biostrings::readAAStringSet(file)

        seqName <- names(fastaFile)
        sequence <- paste(fastaFile)
        # data frame contains all sequences from input file
        fa <- data.frame(seqName, sequence)

        for (j in seq_len(nrow(dataIn))) {
            seqID <- paste0(
                as.character(dataIn$geneID[j]), "|ncbi",
                as.character(dataIn$ncbiID[j]), "|",
                as.character(dataIn$orthoID[j])
            )

            seq <- fa$sequence[pmatch(seqID, fa$seqName)]

            if (length(seq[1]) < 1) {
                fastaOut <- paste0(seqID,
                                    " not found in ",
                                    file,
                                    "! Please check again!")
            } else{
                fastaOut <- paste(paste0(">", seqID), seq[1], sep = "\n")
            }
            fastaOutDf <- rbind(fastaOutDf, as.data.frame(fastaOut))
        }
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
                fastaFile <- Biostrings::readAAStringSet(file)

                seqName <- names(fastaFile)
                sequence <- paste(fastaFile)
                # data frame contains all sequences from input file
                fa <- data.frame(seqName, sequence)

                # get selected sequences
                for (j in seq_len(nrow(dataIn))) {
                    seqID <- as.character(dataIn$orthoID[j])
                    groupID <- as.character(dataIn$geneID[j])
                    seq <- fa$sequence[pmatch(seqID, fa$seqName)]
                    flag <- 1
                    if (is.na(seq)) {
                        seqID <- paste0(
                            as.character(dataIn$geneID[j]), "|ncbi",
                            as.character(dataIn$ncbiID[j]), "|",
                            as.character(dataIn$orthoID[j])
                        )
                        seq <- fa$sequence[pmatch(seqID, fa$seqName)]
                        flag <- 0
                    }
                    if (length(seq[1]) < 1) {
                        fastaOut <-
                            paste0(
                                seqID,
                                " not found in ",
                                file,
                                "! Please check the header format in
                                FASTA file!"
                            )
                    } else {
                        if (!is.na(seq[1])) {
                            if (flag == 1) {
                                fastaOut <- paste(
                                    paste0(">", seqID), seq[1], sep = "\n"
                                )
                            } else {
                                fastaOut <- paste(
                                    paste0(">", seqID), seq[1], sep = "\n"
                                )
                            }

                        } else {
                            fastaOut <-
                                paste0(
                                    seqID,
                                    " not found in uploaded FASTA file!!! ",
                                    "Please check again!!!"
                                )
                        }
                    }
                    fastaOutDf <- rbind(
                        fastaOutDf, as.data.frame(fastaOut)
                    )
                }
            } else {
                if (inputType != "fasta") {
                    fastaOut <- {
                        paste0(
                            "Please provide FASTA file(s) in ",
                            "Input & settings page!"
                        )
                    }
                    fastaOutDf <- rbind(
                        fastaOutDf, as.data.frame(fastaOut)
                    )
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
                fa <- data.frame()
                if (idFormat == 4) {
                    fileList <- list.files(path, pattern = fileExt)
                    for (file in fileList){
                        fileWithPath = paste0(path, "/", file)
                        fastaFile <-
                            Biostrings::readAAStringSet(fileWithPath)

                        seqName <- names(fastaFile)
                        sequence <- paste(fastaFile)
                        # data frame contains all sequences from input file
                        fa <- rbind(fa, data.frame(seqName, sequence))
                    }
                } else {
                    for (
                        i in
                        seq_len(length(levels(as.factor(specDf$specID))))
                    ) {
                        specID <- as.character(
                            levels(as.factor(specDf$specID))[i]
                        )

                        # full path fasta file
                        file <- paste0(path, "/", specID, ".", fileExt)
                        if (dirFormat == 2) {
                            file <- paste0(
                                path, "/", specID, "/", specID, ".", fileExt
                            )
                        }

                        # read fasta file and save sequences into dataframe
                        if (file.exists(file)) {
                            fastaFile <- Biostrings::readAAStringSet(file)
                            seqName <- names(fastaFile)
                            sequence <- paste(fastaFile)
                            # data frame contains all sequences from input file
                            fa <- rbind(fa, data.frame(seqName, sequence))
                        }
                    }
                }

                # now get selected sequences
                if (nrow(fa) > 0) {
                    for (j in seq_len(nrow(dataIn))) {
                        seqID <- as.character(dataIn$orthoID[j])
                        groupID <- as.character(dataIn$geneID[j])

                        seq <- fa$sequence[pmatch(seqID, fa$seqName)]

                        if (length(seq[1]) < 1) {
                            fastaOut <- {
                                paste0(
                                    seqID,
                                    " not found in ",
                                    file,
                                    "! Please check idFormat in FASTA config ",
                                    "again!"
                                )
                            }
                        } else {
                            fastaOut <- paste(
                                paste0(">", seqID), seq[1], sep = "\n"
                            )
                        }
                        fastaOutDf <- rbind(
                            fastaOutDf, as.data.frame(fastaOut)
                        )
                    }
                } else {
                    fastaOut <-
                        paste0(
                            "No fasta file has been found in ",
                            path,
                            "!!! Please check the full path to FASTA folder ",
                            "and the idFormat (header format) in FASTA config",
                            " again!!!"
                        )
                    fastaOutDf <- rbind(
                        fastaOutDf, as.data.frame(fastaOut)
                    )
                }
            } else {
                fastaOut <- {
                    paste0(
                        "Please provide FASTA files in Input & settings page!"
                    )
                }
                fastaOutDf <- rbind(fastaOutDf, as.data.frame(fastaOut))
            }
        }
    }
    # remove duplicated sequences
    fastaOutDf <- fastaOutDf[!duplicated(fastaOutDf), ]

    return(fastaOutDf)
}