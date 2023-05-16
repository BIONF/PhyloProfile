#' Import user defined NCBI taxonomy database
#' @return none (imported files will be saved in PhyloProfile package folder)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
importNcbiTax <- function(inTaxDir, taxDB) {
    taxFiles <- c(
        "newTaxa.txt",
        "taxonNamesReduced.txt",
        "idList.txt",
        "rankList.txt",
        "taxonomyMatrix.txt"
    )
    # Check required files
    message("1) Checking taxonomy files in ", inTaxDir,"...")
    flag = 1
    for (file in taxFiles) {
        if (!file.exists(paste0(inTaxDir, "/", file))) flag = 0
    }

    if (!file.exists(paste0(inTaxDir, "/preCalcTree.nw"))) {
        if (file.exists(paste0(taxDB, "/preCalcTree.nw")))
            unlink(paste0(taxDB, "/preCalcTree.nw"))
    } else taxFiles <- c(taxFiles, "preCalcTree.nw")

    if (flag == 0) {
        message(
            "<p><span style=\"color: #ff0000;\"><strong>ERROR</strong></span>: Some of the taxonomy files cannot be found! Please check <a href=\"https://github.com/BIONF/PhyloProfile/wiki/PhyloProfile-and-the-NCBI-taxonomy-database\">this link</a> for more info.</p>"
        )
    } else {
        message("2) Importing data into ", taxDB, "...")
        for (file in taxFiles) {
            system(
                paste(
                    "cp",
                    paste0(inTaxDir, "/", file),
                    paste0(taxDB, "/", file)
                )
            )
        }

        message(
            "3) FINISHED! Your NCBI taxonomy database has been imported!<br>"
        )
        message("Now you can re-upload your input or restart PhyloProfile.")
        message(
            "To get the original taxonomy data, plase use the Reset function!"
        )
    }
}
