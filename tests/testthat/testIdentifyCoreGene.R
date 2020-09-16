context("test identification of core genes for a selected set of taxa")

test_that("test core gene estimation", {
    data("fullProcessedProfile", package="PhyloProfile")
    rankName <- "class"
    refTaxon <- "Mammalia"
    taxaCore <- c("Mammalia", "Saccharomycetes", "Insecta")
    var1Cutoff <- c(0,1)
    var2Cutoff <- c(0,1)
    percentCutoff <- c(0,1)
    coreCoverage <- 100
    taxonIDs <- levels(as.factor(fullProcessedProfile$ncbiID))
    sortedInputTaxa <- sortInputTaxa(
        taxonIDs, rankName, refTaxon, NULL
    )
    taxaCount <- plyr::count(sortedInputTaxa, "supertaxon")
    coreGene <- getCoreGene(
        rankName,
        taxaCore,
        fullProcessedProfile,
        taxaCount,
        var1Cutoff, var2Cutoff,
        percentCutoff, coreCoverage
    )
    expect_true(length(coreGene) == 3)
})
