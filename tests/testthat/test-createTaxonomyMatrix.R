context("test creating taxonomy matrix function")

test_that("get taxonomy IDs for list of taxa", {
  allTaxonInfo <- data.table::fread("taxonNames_test.txt")
  input_taxa <- c("245562")

  taxonomy_info <- get_ids_rank(input_taxa, allTaxonInfo)
  reducedInfoList <- as.data.frame(taxonomy_info[3])

  expect_that(nrow(reducedInfoList), equals(14))
})

test_that("create taxonomy matrix file", {
  tax_matrix <- taxonomy_table_creator("idList_test.txt", "rankList_test.txt")
  expect_that(tax_matrix, is_a("data.frame"))
})
