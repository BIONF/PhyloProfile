#' Get fasta sequences for AMPK-TOR and BUSCO Arthropoda online demo data
#' @param seqIDs list of sequences IDs. Set seqIDs = "all" if you want to get
#' all fasta sequences from the data set.
#' @param demoData name of demo data set (either "ampk-tor" or "arthropoda").
#' Default = "arthropoda".
#' @return A dataframe with one column contains sequences in fasta format.
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @seealso \code{\link{mainLongRaw}}
#' @examples
#' data("mainLongRaw", package="PhyloProfile")
#' seqIDs <- unique(as.character(mainLongRaw$orthoID))
#' getFastaDemo(seqIDs, "arthropoda")

getFastaDemo <- function(seqIDs = NULL, demoData = "arthropoda") {
    if (is.null(seqIDs)) return()
    if (demoData == "ampk-tor" | demoData == "arthropoda") {
        fastaURL <- paste0(
            "https://raw.githubusercontent.com/BIONF/phyloprofile-data/",
            "master/expTestData/arthropoda/arthropoda.extended.fa"
        )
        if (demoData == "ampk-tor") {
            fastaURL <- paste0(
                "https://raw.githubusercontent.com/BIONF/",
                "phyloprofile-data/master/expTestData/ampk-tor/",
                "ampk-tor.extended.fa"
            )
        }
        
        if (RCurl::url.exists(fastaURL)) {
            # load fasta file
            faFile <- Biostrings::readAAStringSet(fastaURL)
            
            # get sequences
            if (length(seqIDs) == 1 & seqIDs[1] == "all")
                seqIDs <- names(faFile)
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
}