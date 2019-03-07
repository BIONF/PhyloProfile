context("test parsing main input file into long dataframe")

test_that("test fasta input", {
    a <- PhyloProfile::createLongMatrix("mainFastaTest.txt")
    expect_true(ncol(a) == 5)
})

test_that("test wide input", {
    a <- PhyloProfile::createLongMatrix("mainWideTest.txt")
    expect_true(nrow(a) == 8)
})

test_that("test xml input", {
    a <- PhyloProfile::createLongMatrix("mainOrthoxmlTest.xml")
    expect_true(ncol(a) == 6)
})
