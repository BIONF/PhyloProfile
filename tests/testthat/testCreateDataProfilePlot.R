context("test creating data for plotting profiles")

test_that("test profile plot data generation", {
    data("superTaxonProfile", package="PhyloProfile")

    plotDf <- dataMainPlot(superTaxonProfile)
    expect_true(nrow(plotDf) == 91)

    selectedTaxa <- c("Mammalia", "Saccharomycetes", "Insecta")
    selectedSeq <- "all"
    customizedPlotDf <- dataCustomizedPlot(
        superTaxonProfile, selectedTaxa, selectedSeq
    )
    expect_true(nrow(customizedPlotDf) == 11)
})
