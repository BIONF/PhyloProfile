#' Import user defined NCBI taxonomy database
#' @return none (imported files will be saved in PhyloProfile package folder)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
exportNcbiTax <- function(outDir, taxDB) {
    taxFiles <- c(
        "newTaxa.txt",
        "taxonNamesReduced.txt",
        "idList.txt",
        "rankList.txt",
        "taxonomyMatrix.txt"
    )
    # Check required files
    message("1) Checking taxonomy files in ", outDir,"...")
    flag = 1
    for (file in taxFiles) {
        if (file.exists(paste0(outDir, "/", file))) flag = 0
    }

    if (flag == 0) {
        message(
            "<p><span style=\"color: #ff0000;\"><strong>WARNING</strong></span>: Some of the taxonomy files already exist in the selected folder! Please rename or move them to another place!</p>"
        )
    } else {
        message("2) Exporting data from ", taxDB, "...")
        if (file.exists(paste0(taxDB, "/preCalcTree.nw")))
            taxFiles <- c(taxFiles, "preCalcTree.nw")
        for (file in taxFiles) {
            system(
                paste(
                    "cp",
                    paste0(taxDB, "/", file),
                    paste0(outDir, "/", file)
                )
            )
        }
        message(
            "3) FINISHED! The current NCBI taxonomy database has been saved in ", outDir, "!<br>"
        )
    }
}
