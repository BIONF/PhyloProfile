context("test functions for getting data from OMA browser")

test_that("test for a wrong OMA ID", {
    id <- "RATNOBLABLA"
    expect_error(checkOmaID(id))
})

test_that("test for getting HOG orthologs", {
    id <- "HUMAN29397"
    orthologs <- getOmaMembers(id, "HOG")
    expect_true(length(orthologs) > 1)
})

test_that("test for getting protein annotations and fasta for one OMA id", {
    id <- "HUMAN29397"
    omaDf <- getOmaDataForOneOrtholog(id)
    expect_true(ncol(omaDf) == 5)
})

test_that("test for parsing annotation", {
    id <- "HUMAN29397"
    omaDf <- getDataForOneOma(id, "HOG")
    domainDf <- getAllDomainsOma(omaDf)
    expect_true(ncol(domainDf) == 6)
})


