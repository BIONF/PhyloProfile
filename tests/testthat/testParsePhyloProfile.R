context("test parsing and pre-processing phyloprofile input")

test_that("test connection to taxonomy files and getting input taxa", {
    rankName <- "species"
    inputDf <- phyloprofile::createLongMatrix("mainWideTest.txt")

    inputTaxonID <- phyloprofile::getInputTaxaID(inputDf)
    expect_true(length(inputTaxonID) == 4)

    inputTaxonName <- phyloprofile::getInputTaxaName(rankName, inputTaxonID)
    expect_true(nrow(inputTaxonName) == 4)
})
