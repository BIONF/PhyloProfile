context("test creating data for plotting profiles")

test_that("test profile plot data generation", {
    data("finalProcessedProfile", package="PhyloProfile")

    plotDf <- dataMainPlot(finalProcessedProfile)
    expect_true(nrow(plotDf) == 88)

    selectedTaxa <- c("Mammalia", "Saccharomycetes", "Insecta")
    selectedSeq <- "all"
    customizedPlotDf <- dataCustomizedPlot(
        finalProcessedProfile, selectedTaxa, selectedSeq
    )
    expect_true(nrow(customizedPlotDf) == 11)
})
