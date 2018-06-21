#' Detailed plot
#' 
#' @export
#' @param data (from reactive fn "detail_plotDt")
#' @param var1_id name of variable 1 (from input$var1_id)
#' @param var2_id name of variable 2 (from input$var2_id)
#' @param detailed_text text size (from input$detailed_text)
#' @param detailed_height plot height (from input$detailed_height)
#' @return information of selected protein on detailed plot
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

create_detailed_plot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("detail_plot.ui")),
    downloadButton(ns("download_detailed"), "Download plot", class = "butDL"),
    tags$head(
      tags$style(HTML(
        ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
    ),
    hr(),
    verbatimTextOutput(ns("detail_click"))
  )
}

create_detailed_plot <- function(input, output, session, data,
                                 var1_id, var2_id,
                                 detailed_text, detailed_height){
  
  # render detailed plot -------------------------------------------------------
  output$detail_plot <- renderPlot({
    detail_plot(data(), detailed_text(), var1_id(), var2_id())
  })
  
  output$detail_plot.ui <- renderUI({
    ns <- session$ns
    withSpinner(
      plotOutput(
        ns("detail_plot"),
        width = 800,
        height = detailed_height(),
        click = ns("plot_click_detail")
      )
    )
  })

  output$download_detailed <- downloadHandler(
    filename = function() {
      c("detailedPlot.pdf")
    },
    content = function(file) {
      g <- detail_plot(detai_plotDt(), detailed_text(), var1_id(), var2_id())
      ggsave(
        file,
        plot = g,
        width = 800 * 0.056458333,
        height = detailed_height() * 0.056458333,
        units = "cm",
        dpi = 300,
        device = "pdf",
        limitsize = FALSE
      )
    }
  )

  # get info when clicking on detailed plot ------------------------------------
  point_infoDetail <- reactive({
    selDf <- data()
    selDf$orthoID <- as.character(selDf$orthoID)
    
    # get coordinates of plot_click_detail
    if (is.null(input$plot_click_detail$x)) return()
    else{
      corX <- round(input$plot_click_detail$y)
      corY <- round(input$plot_click_detail$x)
    }
    
    # get pair of sequence IDs & var1
    seedID <- as.character(selDf$geneID[!is.na(selDf$geneID)][1])
    orthoID <- as.character(selDf$orthoID[corX])
    
    var1 <- as.list(selDf$var1[selDf$orthoID == orthoID])
    var1 <- as.character(var1[!is.na(var1)])
    var2 <- as.list(selDf$var2[selDf$orthoID == orthoID])
    var2 <- as.character(var2[!is.na(var2)])
    if (length(var2) == 0) var2 = "NA"

    ncbiID <- selDf[selDf$orthoID == orthoID, ]$abbrName
    ncbiID <- as.character(ncbiID[!is.na(ncbiID)][1])

    # return info
    if (is.na(orthoID)) {
      return(NULL)
    } else {
      if (orthoID != "NA") {
        info <- c(seedID, orthoID, var1, var2, ncbiID)
        return(info)
      }
    }
  })
  
  # * show info when clicking on detailed plot ----------------------------------
  output$detail_click <- renderText({
    info <- point_infoDetail() # info = seedID, orthoID, var1

    if (is.null(info)) paste("select ortholog")
    else{
      a <- paste0("seedID = ", info[1])
      b <- paste0("hitID = ", info[2])
      c <- ""
      if (var1_id() != "") {
        c <- paste0(var1_id(), " = ", info[3])
      }
      d <- ""
      if (var2_id() != "") {
        d <- paste0(var2_id(), " = ", info[4])
      }
      paste(a, b, c, d, sep = "\n")
    }
  })
  
  return(point_infoDetail)
}


#' create detailed plot
#' @export
#' @param sel_df data for plotting  (from reactive fn "detail_plotDt")
#' @param detailed_text text size (input$detailed_text)
#' @param var1_id name of variable 1 (input$var1_id)
#' @param var2_id name of variable 2 (input$var2_id)
#' @return detailed plot (ggplot object)
#' @author  Vinh Tran {tran@bio.uni-frankfurt.de}

detail_plot <- function(sel_df, detailed_text, var1_id, var2_id){
  sel_df$x_label <- paste(sel_df$orthoID,
                         " (",
                         sel_df$fullName,
                         ")",
                         sep = "")

  # create joined DF for plotting var1 next to var2
  var1Df <- subset(sel_df, select = c("x_label", "var1"))
  var1Df$type <- var1_id
  colnames(var1Df) <- c("id", "score", "var")

  var2Df <- subset(sel_df, select = c("x_label", "var2"))
  var2Df$type <- var2_id
  colnames(var2Df) <- c("id", "score", "var")

  detailed_df <- rbind(var1Df, var2Df)

  # remove ONE missing variable
  if (nlevels(as.factor(detailed_df$var)) > 1) {
    detailed_df <- detailed_df[nchar(detailed_df$var) > 0, ]
  }

  # keep order of ID (x_label)
  detailed_df$id <- factor(detailed_df$id, levels = unique(detailed_df$id))

  # create plot
  gp <- ggplot(detailed_df, aes(y = score, x = id, fill = var)) +
    geom_bar(stat = "identity", position = position_dodge(), na.rm = TRUE) +
    coord_flip() +
    labs(x = "") +
    labs(fill = "") +
    theme_minimal()
  gp <- gp + theme(axis.text.x = element_text(angle = 90, hjust = 1),
                axis.text = element_text(size = detailed_text),
                axis.title = element_text(size = detailed_text),
                legend.text = element_text(size = detailed_text)
  )
  
  return(gp)
}
