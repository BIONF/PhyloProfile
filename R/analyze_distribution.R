#' Distribution plots
#'
#' @export
#' @param data data for plotting (from reactive fn "presSpecAllDt")
#' @param var_id name of variable (either input$var1_id, input$var2_id or
#' "% present taxa"; from input$selected_dist)
#' @param var_type type of variable (either var1, var2 or presSpec)
#' @param percent percentage cutoff (from input$percent)
#' @param dist_text_size text size of distribution plot
#' (from input$dist_text_size)
#' @return
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

analyze_distribution_ui <- function(id) {
  ns <- NS(id)
  tagList(
    column(
      2,
      tags$head(
        tags$style(HTML(
          "#plot_download_dist{background-color:#A9E2F3}"))
      ),
      downloadButton("plot_download_dist", "Download plot")
    ),
    column(
      10,
      uiOutput(ns("dist_plot.ui"))
    )
  )
}

analyze_distribution <- function(input, output, session,
                                 data,
                                 var_id, var_type,
                                 percent,
                                 dist_text_size){

  # render dist_plot.ui -------------------------------------------------------
  output$dist_plot.ui <- renderUI({
    ns <- session$ns
    withSpinner(plotOutput(ns("distribution_plot")))
  })

  output$distribution_plot <- renderPlot(width = 512, height = 356, {
    var_dist_plot(data(), var_id(), var_type(), percent(), dist_text_size())
  })

  output$plot_download_dist <- downloadHandler(
    filename = function() {
      paste0("distributionPlot.pdf")
    },
    content = function(file) {
      ggsave(
        file,
        plot = var_dist_plot(
          data(),
          var_id(), var_type(),
          percent(), dist_text_size()
        ),
        dpi = 300, device = "pdf", limitsize = FALSE)
    }
  )
}


#' Create distribution plot function
var_dist_plot <- function(data, var_id, var_type, percent, dist_text_size){
  if (var_type == "presSpec") {
    # remove presSpec < cutoff_min or > cutoff_max
    if (percent[1] > 0) {
      data <- data[data$presSpec >= percent[1] & data$presSpec <= percent[2], ]
    } else {
      if (percent[2] > 0) {
        data <- data[data$presSpec > 0 & data$presSpec <= percent[2], ]
      } else {
        data <- data[data$presSpec > 0, ]
      }
    }
  } else {
    data <- data[!is.na(data[,var_type]), ]
  }

  data.mean <- mean(data[,var_type])

  p <- ggplot(data, aes(x = data[,var_type])) +
    geom_histogram(binwidth = .01, alpha = .5, position = "identity") +
    geom_vline(
      data = data,
      aes(xintercept = data.mean, colour = "red"),
      linetype = "dashed",
      size = 1
    ) +
    theme_minimal()
  p <- p + 
    theme(
      legend.position = "none",
      axis.title = element_text(size = dist_text_size),
      axis.text = element_text(size = dist_text_size)
    ) +
    labs(
      x = paste0(var_id,
                 " (mean = ",
                 round(mean(data[,var_type]), 3),
                 ")"),
      y = "Frequency"
    )
  return(p)
}
