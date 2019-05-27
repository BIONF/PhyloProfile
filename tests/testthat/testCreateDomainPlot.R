context("test creating domain plot for a pair seed and ortholog proteins")

test_that("test domain plot", {
    info <- c("101621at6656", "101621at6656|AGRPL@224129@0|224129_0:001955|1")
    domainDf <- parseDomainInput(
        "101621at6656","domains/101621at6656.domains","file"
    )
    expect_true(nrow(domainDf) == 3)
})
