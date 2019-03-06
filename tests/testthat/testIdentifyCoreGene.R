context("test identification of core genes for a selected set of taxa")

test_that("test core gene estimation", {
    # load full processed data
    data("fullProcessedProfile", package="phyloprofile")

    # calculate gene ages
    rankName <- "class"
    taxaCore <- c("Mammalia", "Mucorales", "Alphaproteobacteria")
    var1Cutoff <- c(0,1)
    var2Cutoff <- c(0,1)
    percentCutoff <- c(0,1)
    coreCoverage <- 100
    coreGene <- getCoreGene(rankName,
                               taxaCore,
                               fullProcessedProfile,
                               var1Cutoff, var2Cutoff,
                               percentCutoff, coreCoverage)

    expect_true(length(coreGene) == 1)
})
