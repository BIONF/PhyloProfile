context("test creating data for plotting profiles")

test_that("test profile plot data generation", {
    data("fullProcessedProfile", package="phyloprofile")
    
    plotDf <- dataMainPlot(fullProcessedProfile)
    expect_true(nrow(plotDf) == 20)
    
    selectedTaxa <- c("Mammalia", "Echinoidea", "Gunneridae")
    selectedSeq <- "all"
    customizedPlotDf <- dataCustomizedPlot(
        fullProcessedProfile, selectedTaxa, selectedSeq
    )
    expect_true(nrow(customizedPlotDf) == 10)
})
