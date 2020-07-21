#' Reset taxonomy data
#' @return none
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
resetTaxData <- function() {
    packgePath <- find.package("PhyloProfile")
    data <- c(
        "idList.txt", "rankList.txt", "taxonomyMatrix.txt",
        "taxonNamesReduced.txt"
    )
    files <- paste0(packgePath, "/PhyloProfile/data/", data)
    # delete old files
    message("1) Deleting old taxonomy data...")
    for(f in files) {
        if (file.exists(f)) file.remove(f)
    }
    # rewrite precalculated data again
    message("2) Re-writing taxonomy data files...")
    data(rankList)
    write.table(
        rankList, file = paste0(packgePath, "/PhyloProfile/data/rankList.txt"),
        col.names = FALSE,
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )
    data(idList)
    write.table(
        idList, file = paste0(packgePath, "/PhyloProfile/data/idList.txt"),
        col.names = FALSE,
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )
    data(taxonNamesReduced)
    write.table(
        taxonNamesReduced,
        file = paste0(packgePath, "/PhyloProfile/data/taxonNamesReduced.txt"),
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )
    data(taxonomyMatrix)
    write.table(
        taxonomyMatrix,
        file = paste0(packgePath, "/PhyloProfile/data/taxonomyMatrix.txt"),
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )
    message("2) FINISHED!<br>")
    message("Now you can restart PhyloProfile or re-upload your input file.")
}
