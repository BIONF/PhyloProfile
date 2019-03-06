context("test searching taxonomy ID function")

test_that("return ID for a taxon", {
    taxon <- "homo sapiens"
    a <- searchTaxonIDOnline(taxon)
    expect_that(a, is_a("data.frame"))
})
