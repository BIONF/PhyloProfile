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
    if (is.null(seqIDs)) stop("No sequence ID given!")
    demoBaseUrl<-"https://github.com/BIONF/phyloprofile-data/raw/master/rdata/"
    if (demoData == "ampk-tor" | demoData == "arthropoda") {
        if (demoData == "ampk-tor") {
            load(url(paste0(demoBaseUrl, "ampkTorFasta.RData")))
            faFile <- ampkTorFasta
        } else {
            load(url(paste0(demoBaseUrl, "arthropodaFasta.RData")))
            faFile <- arthropodaFasta
        }
        
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
    }
}
