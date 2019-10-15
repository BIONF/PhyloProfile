#' Update NCBI taxonomy database
#' @return none (updated files will be saved in PhyloProfile package folder)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
updateNcbiTax <- function() {
    # get latest NCBI taxonomy database and processing
    message("1) Getting latest NCBI taxonomy database and processing...")
    preProcessedTaxonomy <- processNcbiTaxonomy()
    # save to preProcessedTaxonomy.txt file within PhyloProfile pkg folder
    message("2) Saving new data into preProcessedTaxanomy.txt file...")
    packgePath <- find.package("PhyloProfile")
    write.table(
        preProcessedTaxonomy,
        file = paste0(packgePath,"/PhyloProfile/data/preProcessedTaxonomy.txt"),
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )
    message("3) FINISHED! Your NCBI taxonomy database has been updated!<br>")
    message("Now you can re-upload your input or restart PhyloProfile.")
}
