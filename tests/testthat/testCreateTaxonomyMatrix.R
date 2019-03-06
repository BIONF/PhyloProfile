context("test creating taxonomy matrix function")

test_that("get taxonomy IDs for list of taxa", {
    allTaxonInfo <- data.table::fread("taxonNamesTest.txt")
    inputTaxa <- c("272557")

    taxonomyInfo <- getIDsRank(inputTaxa, allTaxonInfo)
    reducedInfoList <- as.data.frame(taxonomyInfo[3])

    expect_that(nrow(reducedInfoList), equals(10))
})

test_that("create taxonomy matrix file", {
    taxMatrix <- taxonomyTableCreator("idListTest.txt", "rankListTest.txt")
    expect_that(taxMatrix, is_a("data.frame"))
})
