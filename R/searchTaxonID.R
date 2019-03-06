#' @title Search NCBI taxonomy IDs for a list of taxon names
#' @param taxa list of taxon names
#' @importFrom taxize get_uid
#' @importFrom taxize gnr_resolve
#' @return dataframe contains NCBI taxonomy IDs for all input taxa
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
#' @export
#' @examples
#' taxa <- c("homo sapiens")
#' searchTaxonIDOnline(taxa)

searchTaxonIDOnline <- function(taxa){
    idDf <- data.frame("name" = character(),
                        "newName" = character(),
                        "id" = character(),
                        "type" = character(),
                        stringsAsFactors = FALSE)
    for (i in seq_len(length(taxa))) {
        idDf[i,] <- getTaxonID(taxa[i])
    }
    return(idDf)
}

getTaxonID <- function(taxName){
    id <- taxize::get_uid(sciname = taxName)[1]

    idDf <- data.frame("name" = character(),
                        "newName" = character(),
                        "id" = character(),
                        "type" = character(),
                        stringsAsFactors = FALSE)

    if (is.na(id)) {
        temp <- taxize::gnr_resolve(names = as.character(taxName))
        if (nrow(temp) > 0) {
            newID <- taxize::get_uid(sciname = temp[1, 3])[1]
            if (is.na(newID)) {
                idDf[1, ] <- c(as.character(taxName),
                                as.character(temp[1, 3]),
                                paste0("NA"), "notfound")
            } else {
                idDf[1, ] <- c(as.character(taxName),
                                as.character(temp[1, 3]),
                                paste0("ncbi", newID),
                                "notfound")
            }
        } else {
            idDf[1, ] <- c(as.character(taxName),
                            paste0("no alternative"),
                            paste0("NA"),
                            "notfound")
        }
    } else {
        idDf[1, ] <- c(as.character(taxName),
                        "NA",
                        paste0("ncbi", id),
                        "retrieved")
    }
    return(idDf)
}
