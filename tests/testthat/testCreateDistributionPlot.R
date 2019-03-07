context("test creating data and plots for distritbution analysis function")

test_that("test distritbution analyze functions", {
    # raw input
    inputData <- PhyloProfile::createLongMatrix("mainWideTest.txt")

    # selected rank name
    rankName <- "species"

    # variable thresholds
    var1CutoffMin <- 0.0
    var1CutoffMax <- 1.0
    var2CutoffMin <- 0.0
    var2CutoffMax <- 1.0

    # distribution data for percentage of species present in supertaxa
    percentDistributionData <- createPercentageDistributionData(
        inputData, rankName
    )
    expect_true(nrow(percentDistributionData) == 6)

    # distribution data for 2 additional variables
    distributionData <- createVariableDistributionData(
        inputData,
        var1CutoffMin,
        var1CutoffMax,
        var2CutoffMin,
        var2CutoffMax
    )
    expect_true(nrow(distributionData) == 6)

    # distribution data for 2 additional variables of a subset of taxa
    inputTaxonID <- PhyloProfile::getInputTaxaID(inputData)
    inputTaxonName <- PhyloProfile::getInputTaxaName(rankName,
                                                         inputTaxonID)
    refTaxon <- inputTaxonName$fullName[1]
    taxaTree <- NULL

    sortedTaxa <- PhyloProfile::sortInputTaxa(inputTaxonID,
                                                 inputTaxonName,
                                                 rankName,
                                                 refTaxon,
                                                 taxaTree)

    fullProfileData <- PhyloProfile::parseInfoProfile(inputData,
                                                          sortedTaxa,
                                                          "max",
                                                          "mean")
    selectedGenes <- c("OG_1017", "OG_1019")
    selectedTaxa <- c("Arabidopsis thaliana", "Encephalitozoon intestinalis")
    subsetDistributionData <- createVariableDistributionDataSubset(
        fullProfileData,
        distributionData,
        selectedGenes,
        selectedTaxa
    )
    expect_true(ncol(subsetDistributionData) == 5)

    # plot distribution of var1
    p <- createVarDistPlot(distributionData,
                              "variable 1 name",
                              "var1",
                              NULL,
                              12)
    expect_true(nrow(p$data) == 6)
})
