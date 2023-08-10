context("test parsing domain input file(s) into dataframe")

test_that("test with single file", {
    a <- parseDomainInput(NULL, "domains/101621at6656.domains", "file")
    expect_true(ncol(a) == 10)
})

test_that("test with domain folder", {
    a <- parseDomainInput("101621at6656", "domains", "folder")
    expect_true(ncol(a) == 10)

    b <- parseDomainInput("OG_1020", "domains", "folder")
    expect_true(b == "noFileInFolder")
})
