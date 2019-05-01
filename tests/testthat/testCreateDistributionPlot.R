context("test creating data and plots for distritbution analysis function")

test_that("test distritbution analyze functions", {
    # raw input
    inputData <- createLongMatrix("mainWideTest.txt")

    # selected rank name
    rankName <- "species"

    # distribution data for percentage of species present in supertaxa
    percentDistributionData <- createPercentageDistributionData(
        inputData, rankName
    )
    expect_true(nrow(percentDistributionData) == 6)

    # distribution data for 2 additional variables
    distributionData <- createVariableDistributionData(
        inputData, c(0, 1), c(0.5, 1)
    )
    expect_true(nrow(distributionData) == 4)

    # distribution data for 2 additional variables of a subset of taxa
    inputTaxonID <- getInputTaxaID(inputData)
    inputTaxonName <- getInputTaxaName(rankName, inputTaxonID)
    refTaxon <- inputTaxonName$fullName[1]
    taxaTree <- NULL

    sortedTaxa <- sortInputTaxa(
        inputTaxonID, rankName, refTaxon, taxaTree
    )

    fullProfileData <- parseInfoProfile(
        inputData,
        sortedTaxa,
        "max", "mean"
    )
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
    p <- createVarDistPlot(
        distributionData, "variable 1 name", "var1", NULL, 12
    )
    expect_true(nrow(p$data) == 4)
})
