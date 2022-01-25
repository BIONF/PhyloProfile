#' Download plot seetings
#' @return a yaml file or rscript file
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}
downloadPlotSettings <- function (settings, outFile, type) {
    if (file.exists(outFile)) {
        msg <- paste(
            "Output file exists! Please choose another directory or rename",
            "the output file to continue!"
        )
        message(msg)
    } else {
        if (type == "list") {
            yaml::write_yaml(settings, outFile, indent = 4)
            message("FINISHED! The output file is ", outFile)
        } else {
            yaml::write_yaml(settings, paste0(outFile, ".yml"), indent = 4)
            lapply(
                writePlottingScript(settings), 
                write, paste0(outFile, ".rscript"), append = TRUE
            )
            msg1 <- paste0(
                "FINISHED! The output files are ", outFile,
                ".yml and ", outFile, ".rscript"
            )
            msg2 <- paste0(
                "You can use this command to generate the plot:<br>",
                "Rscript ", outFile, ".rscript ", outFile, ".yml ",
                "input.phyloprofile output_plot.pdf" 
            )
            message(msg1, "<br>", msg2)
        }
    }
}

writePlottingScript <- function(settingsFile) {
    out <- c(
        "library(PhyloProfile)",
        "library(ggplot2)",
        
        "args = commandArgs(trailingOnly = TRUE)",
        "if (is.na(args[1])) stop('Setting file missing')",
        "if (is.na(args[2])) stop('Input file missing')",
        "if (is.na(args[3])) stop('Output file missing')",
        
        "##### Reading setting file",
        paste(
            "if (!file.exists(args[1]))",
            "stop(paste(args[1], 'not found!'))"
            
        ),
        "settings <- yaml::read_yaml(args[1])",
        "print('Reading setting file...')",
        "rawInput <- args[2]",
        paste(
            "if (!file.exists(rawInput))",
            "stop(paste(rawInput, 'not found!'))"
        ),
        "rankName <- settings$rank",
        "refTaxon <- settings$refspec",
        "taxaTree <- NULL",
        "var1AggregateBy <- settings$var1AggregateBy",
        "var2AggregateBy <- settings$var2AggregateBy",
        "percentCutoff <- settings$percentCutoff",
        "coorthologCutoffMax <- 999",
        "var1Cutoff <- settings$var1Cutoff",
        "var2Cutoff <- settings$var2Cutoff",
        "var1Relation <- settings$var1Relation",
        "var2Relation <- settings$var2Relation",
        "groupByCat <- FALSE",
        "catDt <- NULL",
        
        "##### Processing input",
        "print('Processing input...')",
        "inputDf <- createLongMatrix(rawInput)",
        "inputDf <- inputDf[!duplicated(inputDf),]",
        "colnames(inputDf) <- c('geneID', 'ncbiID', 'orthoID', 'var1', 'var2')",
        
        "sortedTaxa <- sortInputTaxa(
            taxonIDs = getInputTaxaID(inputDf),
            rankName = rankName,
            refTaxon = refTaxon,
            taxaTree = NULL
        )",
        
        "taxaCount <- plyr::count(sortedTaxa, 'supertaxon')",
        
        "fullData <- parseInfoProfile(
            inputDf = inputDf,
            sortedInputTaxa = sortedTaxa,
            taxaCount = taxaCount,
            coorthoCOMax = coorthologCutoffMax
        )",
        
        "filteredDf <- filterProfileData(
            fullData,
            taxaCount,
            refTaxon,
            percentCutoff,
            coorthologCutoffMax,
            var1Cutoff,
            var2Cutoff,
            var1Relation,
            var2Relation,
            groupByCat,
            catDt,
            var1AggregateBy,
            var2AggregateBy
        )",
        "dataHeat <- reduceProfile(filteredDf)",
        
        "##### Clustering",
        "if (settings$clusterProfile == TRUE) {",
        "print('Clustering...')",
        "clusteredDataHeat <- dataHeat",
        "clusterMethod <- settings$clusterMethod",
        "distMethod <- settings$distMethodClustering",
        paste0(
            "profiles4cluster <- getDataClustering(
            dataHeat, settings$profileTypeClustering,
            var1AggregateBy,
            var2AggregateBy
            )"
        ),
        "distanceMatrix <- getDistanceMatrix(profiles4cluster, distMethod)",
        "row.order <- hclust(distanceMatrix, method = clusterMethod)$order",
        "datNew <- profiles4cluster[row.order, ]",
        "clusteredGeneIDs <- as.factor(row.names(datNew))",
        paste0(
            "clusteredDataHeat$geneID <- factor(clusteredDataHeat$geneID,",
            "levels=clusteredGeneIDs)"
        ),
        paste0(
            "clusteredDataHeat <- clusteredDataHeat",
            "[!is.na(clusteredDataHeat$geneID),]"
        ),
        "}",
        
        "print('Generating plot...')",
        "plotDf <- dataMainPlot(dataHeat)",
        paste(
            "if (settings$clusterProfile == TRUE)",
            "plotDf <- dataMainPlot(clusteredDataHeat)"
        ),
        "plotParameter <- list(
            'xAxis' = settings$xAxis,
            'var1ID' = settings$var1ID,
            'var2ID' = settings$var2ID,
            'midVar1' = settings$midVar1,
            'midColorVar1' = settings$midColorVar1,
            'lowColorVar1' = settings$lowColorVar1,
            'highColorVar1' = settings$highColorVar1,
            'midVar2' = settings$midVar2,
            'midColorVar2' = settings$midColorVar2,
            'lowColorVar2' = settings$lowColorVar2,
            'highColorVar2' = settings$highColorVar2,
            'paraColor' = settings$paraColor,
            'xSize' = settings$xSize,
            'ySize' = settings$ySize,
            'legendSize' = settings$legendSize,
            'mainLegend' = settings$mainLegend,
            'dotZoom' = settings$dotZoom,
            'xAngle' = settings$xAngle,
            'guideline' = 0,
            'colorByGroup' = settings$colorByGroup
        )",
        "p <- heatmapPlotting(plotDf, plotParameter)",
        "ggsave(args[3], plot = p, width = settings$width * 0.056458333, 
            height = settings$height * 0.056458333,
            units = 'cm', dpi = 300, device = 'pdf', limitsize = FALSE)",
        "print('DONE! Your plot is saved in', args[3])"
    )
    return(out)
}