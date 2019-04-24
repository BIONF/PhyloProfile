context("test parsing and pre-processing PhyloProfile input")

test_that("test connection to taxonomy files and getting input taxa", {
    rankName <- "species"
    inputDf <- createLongMatrix("mainWideTest.txt")

    inputTaxonID <- getInputTaxaID(inputDf)
    expect_true(length(inputTaxonID) == 4)

    inputTaxonName <- getInputTaxaName(rankName, inputTaxonID)
    expect_true(nrow(inputTaxonName) == 4)
})
