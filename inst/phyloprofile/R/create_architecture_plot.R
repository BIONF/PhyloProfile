#' Protein domain architecture plot
#'
#' @export
#' @param point_info() info of clicked point
#' (from reactive fn "point_infoDetail")
#' @param domain_info() domain information
#' (from reactive fn "get_domain_information")
#' @param label_archi_size lable size (from input$label_archi_size)
#' @param title_archi_size title size (from input$title_archi_size)
#' @param archi_height plot height (from input$archi_height)
#' @param archi_width plot width (from input$archi_width)
#' @return
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

source("R/functions.R")

create_architecture_plot_ui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("archi_plot.ui")),
    downloadButton(ns("archi_download"), "Download plot", class = "butDL"),
    tags$head(
      tags$style(HTML(
        ".butDL{background-color:#476ba3;} .butDL{color: white;}"))
    ),
    textOutput(ns("selected_domain"))
  )
}

create_architecture_plot <- function(input, output, session,
                                     point_info,
                                     domain_info,
                                     label_archi_size, title_archi_size,
                                     archi_height, archi_width){
  output$archi_plot <- renderPlot({
    if (is.null(nrow(domain_info()))) return()
    g <- archi_plot(point_info(),
                    domain_info(),
                    label_archi_size(),
                    title_archi_size())
    if (any(g == "ERR_0")) {
      msg_plot()
    } else {
      grid.draw(g)
    }
  })

  output$archi_plot.ui <- renderUI({
    ns <- session$ns
    if (is.null(nrow(domain_info()))) {
      msg <- paste0(
        "<p><em>Wrong domain file has been uploaded!
        Please check the correct format in
        <a href=\"https://github.com/BIONF/PhyloProfile/wiki/Input-Data#ortholog-annotations-eg-domains\"
        target=\"_blank\" rel=\"noopener\">our PhyloProfile wiki</a>.</em></p>"
      )
      HTML(msg)
    } else {
      withSpinner(
        plotOutput(
          ns("archi_plot"),
          height = archi_height(),
          width = archi_width(),
          click = ns("archi_click")
        )
      )
    }
  })

  output$archi_download <- downloadHandler(
    filename = function() {
      c("domains.pdf")
    },
    content = function(file) {
      g <- archi_plot(point_info(),
                      domain_info(),
                      label_archi_size(),
                      title_archi_size())
      grid.draw(g)
      ggsave(
        file, plot = g,
        width = archi_width() * 0.056458333,
        height = archi_height() * 0.056458333,
        units = "cm", dpi = 300, device = "pdf", limitsize = FALSE
      )
    }
  )

  output$selected_domain <- renderText({
    if (is.null(input$archi_click$y)) return()
    convertY(unit(input$archi_click$y, "npc"), "native")
  })
}


#' Create architecure plot for both seed and ortho protein
#' @export
#' @param info info of clicked point (from reactive fn "point_infoDetail")
#' @param domain_df domain info (from reactive fn "get_domain_information")
#' @param label_archi_size lable size (from input$label_archi_size)
#' @param title_archi_size title size (from input$title_archi_size)
#' @return plot as arrangeGrob object
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

archi_plot <- function(info,
                       domain_df,
                       label_archi_size, title_archi_size){
  # info
  group <- as.character(info[1])
  ortho <- as.character(info[2])
  var1 <- as.character(info[3])

  # get sub dataframe based on selected group_id and orthoID
  ortho <- gsub("\\|", ":", ortho)
  grepID <- paste(group, "#", ortho, sep = "")

  subdomain_df <- domain_df[grep(grepID, domain_df$seedID), ]
  subdomain_df$feature <- as.character(subdomain_df$feature)

  if (nrow(subdomain_df) < 1) {
    return(paste0("ERR_0"))
  } else {

    # ortho domains df
    ortho_df <- filter(subdomain_df, orthoID == ortho)

    # seed domains df
    seed_df <- filter(subdomain_df, orthoID != ortho)

    if (nrow(seed_df) == 0) seed_df <- ortho_df

    seed <- as.character(seed_df$orthoID[1])

    # return ERR_0 if seed_df and ortho_df are empty
    if (nrow(seed_df) == 0) return(paste0("ERR_0"))

    # change order of one dataframe's features
    # based on order of other df's features
    if (length(ortho_df$feature) < length(seed_df$feature)) {
      ordered_ortho_df <- ortho_df[order(ortho_df$feature), ]
      ordered_seed_df <- sort_domains(ordered_ortho_df, seed_df)
    } else {
      ordered_seed_df <- seed_df[order(seed_df$feature), ]
      ordered_ortho_df <- sort_domains(ordered_seed_df, ortho_df)
    }

    # join weight values and feature names
    if ("weight" %in% colnames(ordered_ortho_df)) {
      ordered_ortho_df$yLabel <- paste0(ordered_ortho_df$feature,
                                        " (",
                                        round(ordered_ortho_df$weight, 2),
                                        ")")
    } else {
      ordered_ortho_df$yLabel <- ordered_ortho_df$feature
    }
    if ("weight" %in% colnames(ordered_seed_df)) {
      ordered_seed_df$yLabel <- paste0(ordered_seed_df$feature,
                                       " (",
                                       round(ordered_seed_df$weight, 2),
                                       ")")
    } else {
      ordered_seed_df$yLabel <- ordered_seed_df$feature
    }

    # create color scheme for all features
    # the same features in seed & ortholog will have the same colors
    feature_seed <- levels(as.factor(ordered_seed_df$feature))
    feature_ortho <- levels(as.factor(ordered_ortho_df$feature))
    all_features <- c(feature_seed, feature_ortho)
    all_colors <- get_qual_col_for_vector(
      all_features,
      length(all_features)
    )

    color_scheme <- structure(
      all_colors,
      .Names = all_features
    )

    # plotting
    sep <- ":"

    if ("length" %in% colnames(subdomain_df)) {
      plot_ortho <- domain_plotting(ordered_ortho_df,
                                    ortho,
                                    sep,
                                    label_archi_size,
                                    title_archi_size,
                                    min(subdomain_df$start),
                                    max(c(subdomain_df$end,
                                          subdomain_df$length)),
                                    color_scheme)
      plot_seed <- domain_plotting(ordered_seed_df,
                                   seed,
                                   sep,
                                   label_archi_size,
                                   title_archi_size,
                                   min(subdomain_df$start),
                                   max(c(subdomain_df$end,
                                         subdomain_df$length)),
                                   color_scheme)

    } else{
      plot_ortho <- domain_plotting(ordered_ortho_df,
                                    ortho,
                                    sep,
                                    label_archi_size,
                                    title_archi_size,
                                    min(subdomain_df$start),
                                    max(subdomain_df$end),
                                    color_scheme)
      plot_seed <- domain_plotting(ordered_seed_df,
                                   seed,
                                   sep,
                                   label_archi_size,
                                   title_archi_size,
                                   min(subdomain_df$start),
                                   max(subdomain_df$end),
                                   color_scheme)
    }

    # grid.arrange(plot_seed,plot_ortho,ncol=1)
    if (ortho == seed) {
      arrangeGrob(plot_seed, ncol = 1)
    } else {
      seed_height <- length(levels(as.factor(ordered_seed_df$feature)))
      ortho_height <- length(levels(as.factor(ordered_ortho_df$feature)))

      arrangeGrob(plot_seed, plot_ortho, ncol = 1,
                  heights = c(seed_height, ortho_height))
    }
  }
}

#' Create architecure plot for single protein (seed/ortho)
#' @export
#' @param df data (seed/ortho) for ploting
#' @param geneID ID of seed or ortho protein
#' @param sep separate indicator for title
#' @param label_size lable size (from input$label_archi_size)
#' @param title_size title size (from input$title_archi_size)
#' @param min_start the smallest start position of all domains
#' @param max_end the highest stop position of all domains
#' @return plot as ggplot object
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

domain_plotting <- function(df,
                            geneID,
                            sep,
                            label_size, title_size,
                            min_start, max_end,
                            color_scheme){
  gg <- ggplot(df, aes(y = feature, x = end, color = as.factor(feature))) +
    geom_segment(data = df,
                 aes(y = feature, yend = feature,
                     x = min_start, xend = max_end),
                 color = "white",
                 size = 0) +
    scale_color_manual(values = color_scheme)

  # draw lines for representing sequence length
  if ("length" %in% colnames(df)) {
    gg <- gg + geom_segment(data = df,
                            aes(x = 0, xend = length,
                                y = feature, yend = feature),
                            size = 1,
                            color = "#b2b2b2")
  }

  # draw line and points
  gg <- gg + geom_segment(data = df,
                          aes(x = start, xend = end,
                              y = feature, yend = feature),
                          size = 1.5)
  gg <- gg + geom_point(data = df,
                        aes(y = feature, x = start),
                        color = "#b2b2b2",
                        size = 3,
                        shape = 3)
  gg <- gg + geom_point(data = df,
                        aes(y = feature, x = end),
                        color = "#edae52",
                        size = 3,
                        shape = 5)

  # draw dashed line for domain path
  gg <- gg + geom_segment(data = df[df$path == "Y", ],
                          aes(x = start, xend = end,
                              y = feature, yend = feature),
                          size = 3,
                          linetype = "dashed")

  # # add text above
  # gg <- gg + geom_text(data = df,
  #                      aes(x = (start + end) / 2,
  #                          y = feature, label = round(weight,2)),
  #                        color = "#9fb059",
  #                        size = descSize,
  #                        vjust = -0.75,
  #                        fontface = "bold",
  #                        family = "serif")

  # theme format
  title_mod <- gsub(":", sep, geneID)
  gg <- gg + scale_y_discrete(expand = c(0.075, 0),
                              breaks = df$feature,
                              labels = df$yLabel)
  gg <- gg + labs(title = paste0(title_mod), y = "Feature")
  gg <- gg + theme_minimal()
  gg <- gg + theme(panel.border = element_blank())
  gg <- gg + theme(axis.ticks = element_blank())
  gg <- gg + theme(plot.title = element_text(face = "bold", size = title_size))
  gg <- gg + theme(plot.title = element_text(hjust = 0.5))
  gg <- gg + theme(legend.position = "none", axis.title.x = element_blank(),
                   axis.text.y = element_text(size = label_size),
                   axis.title.y = element_text(size = label_size),
                   panel.grid.minor.x = element_blank(),
                   panel.grid.major.x = element_blank())
  # return plot
  return(gg)
}

#' sort one domain dataframe (ortho) based on the other domain Df (seed)
#' @export
#' @param seed_df data of seed protein
#' @param ortho_df data of ortholog protein
#' @return sorted domain list (dataframe)
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

sort_domains <- function(seed_df, ortho_df){
  # get list of features in seed_df
  feature_list <- as.data.frame(levels(as.factor(seed_df$feature)))
  colnames(feature_list) <- c("feature")
  # and add order number to each feature
  feature_list$orderNo <- seq(length(feature_list$feature))

  # merge those info to ortho_df
  ordered_ortho_df <- merge(ortho_df, feature_list, all.x = TRUE)

  # sort ortho_df
  index <- with(ordered_ortho_df, order(orderNo))
  ordered_ortho_df <- ordered_ortho_df[index, ]

  #turn feature column into a character vector
  ordered_ortho_df$feature <- as.character(ordered_ortho_df$feature)
  #then turn it back into an ordered factor (to keep this order while plotting)
  ordered_ortho_df$feature <- factor(ordered_ortho_df$feature,
                                     levels = unique(ordered_ortho_df$feature))
  #return sorted df
  return(ordered_ortho_df)
}

#' plot error message
#' @export
#' @param
#' @return error message in a ggplot object
#' @author Vinh Tran {tran@bio.uni-frankfurt.de}

msg_plot <- function() {
  msg <- paste(
    "No information about domain architecture!",
    "Please check:","if you uploaded the correct domain file/folder; or ",
    "if the selected genes (seed & ortholog) do exist in the uploaded file",
    "(please search for the corresponding seedID and hitID)",
    sep = "\n"
  )
  x <- c(1,2,3,4,5)
  y <- c(1,2,3,4,5)
  g <- ggplot(data.frame(x, y), aes(x,y)) +
    geom_point(color = "white") +
    annotate("text", label = msg, x = 3.5, y = 0.5, size = 5, colour = "red") +
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_blank()) +
    ylim(0,1)
  return(g)
}
