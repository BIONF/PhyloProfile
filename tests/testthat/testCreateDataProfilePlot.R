context("test creating data for plotting profiles")

test_that("test profile plot data generation", {
    data("fullProcessedProfile", package="PhyloProfile")

    plotDf <- dataMainPlot(fullProcessedProfile)
    expect_true(nrow(plotDf) == 168)

    selectedTaxa <- c("Mammalia", "Saccharomycetes", "Insecta")
    selectedSeq <- "all"
    customizedPlotDf <- dataCustomizedPlot(
        fullProcessedProfile, selectedTaxa, selectedSeq
    )
    expect_true(nrow(customizedPlotDf) == 67)
})
