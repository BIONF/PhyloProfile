context("test searching taxonomy ID function")

pkg_path <- find.package("phyloprofile")
r_source <- paste0(pkg_path,"/phyloprofile/R")

source(paste0(r_source,"/functions.R"))

test_that("check internet connection", {
  a <- has_internet()
  expect_true(a)
})

test_that("return ID for a taxon", {
  taxon <- "homo sapiens"
  a <- search_taxonID_online(taxon)
  expect_that(a, is_a("data.frame"))
})
