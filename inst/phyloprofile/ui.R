#' Load packages
packages <- c("shiny", "shinyBS", "shinyjs", "colourpicker", "DT", "devtools")

source("R/functions.R")
install_packages(packages)
lapply(packages, library, character.only = TRUE)

if (!require("shinycssloaders")) {
  devtools::install_github("andrewsali/shinycssloaders")
  library(shinycssloaders)
}

#' Import function files
source_files = list.files(path = "R",
                          pattern = "*.R$",
                          full.names = TRUE)
lapply(source_files, source, .GlobalEnv)

#' MAIN UI ====================================================================
shinyUI(
  fluidPage(

    tags$style(type = "text/css", "body {padding-top: 80px;}"),
    useShinyjs(),

    # Application title
    titlePanel("", windowTitle = "PhyloProfile"),

    # TOP WELLPANEL FOR PLOT CONFIGURATION ------------------------------------
    conditionalPanel(
      condition = "input.tabs=='Main profile'",
      wellPanel(
        fluidRow(
          column(
            2,
            radioButtons(
              inputId = "x_axis",
              label = "Choose type of x-axis:",
              choices = list("taxa", "genes"),
              selected = "taxa",
              inline = TRUE
            ),
            hr(),
            checkboxInput(
              "auto_update",
              strong(em("Auto update plot")),
              value = FALSE,
              width = NULL
            )
          ),
          column(
            1,
            create_plot_size("width", "Width (px)", 600),
            actionButton("main_plot_config", "Appearance")
          ),
          column(
            1,
            create_plot_size("height", "Height (px)", 600)
          ),

          column(
            2,
            uiOutput("var1_cutoff.ui")
          ),
          column(
            2,
            uiOutput("var2_cutoff.ui")
          ),
          column(
            2,
            uiOutput("percent_cutoff.ui")
          ),
          column(
            2,
            numericInput(
              "coortholog",
              "Min co-orthologs",
              min = 1,
              max = 999,
              step = 1,
              value = 1,
              width = 150
            ),
            shinyBS::bsButton(
              "reset_main",
              "Reset cutoffs",
              style = "danger"
            )
          )
        )
      )
    ),

    conditionalPanel(
      condition = "input.tabs=='Customized profile'",
      wellPanel(
        fluidRow(
          column(
            2,
            radioButtons(
              inputId = "x_axis_selected",
              label = "Choose type of x-axis:",
              choices = list("taxa", "genes"),
              selected = "taxa",
              inline = TRUE
            ),
            hr(),
            checkboxInput(
              "auto_update_selected",
              strong(em("Auto update plot")),
              value = FALSE,
              width = NULL
            )
          ),
          column(
            1,
            create_plot_size("selected_width", "Width (px)", 600),
            actionButton("selected_plot_config", "Appearance")
          ),
          column(
            1,
            create_plot_size("selected_height", "Height (px)", 600)
          ),

          column(
            2,
            uiOutput("var1_filter.ui")
          ),
          column(
            2,
            uiOutput("var2_filter.ui")
          ),
          column(
            2,
            uiOutput("percent_filter.ui")
          ),
          column(
            2,
            uiOutput("coortholog_filter.ui"),
            shinyBS::bsButton(
              "reset_selected",
              "Reset cutoffs",
              style = "danger"
            )
          )
        )
      )
    ),

    # MAIN NARVARPAGE TABS ----------------------------------------------------
    navbarPage(
      em(strong("PhyloProfile v0.3.2")),
      id = "tabs",
      collapsible = TRUE,
      inverse = TRUE,
      fluid = TRUE,
      position = "fixed-top",

      # INPUT TAB -------------------------------------------------------------
      tabPanel(
        "Input & settings",
        # * 1st column --------------------------------------------------------
        column(
          4,
          # ** Main input -----------------------------------------------------
          strong(h4("Main input:")),
          conditionalPanel(
            condition = "input.do",
            em(
              strong("RELOAD THIS TOOL TO UPLOAD A NEW INPUT FILE!!!",
                     style = "color:red")
            )
          ),

          selectInput(
            "demo_data", label = h5("Use online demo data:"),
            choices = list("None" = "none",
                           "AMPK-TOR" = "ampk-tor",
                           "LCA Microsporidia" = "lca-micros"),
            selected = "none",
            width = "80%"
          ),

          uiOutput("no_internet_msg"),
          uiOutput("demo_data_describe"),
          uiOutput("main_input_file.ui"),
          uiOutput("input_check.ui"),

          fluidRow(
            column(
              6,
              conditionalPanel(
                condition = "output.check_oma_input",
                shinyBS::bsButton("open_oma_windows", "Get data from OMA"),
                br()
              )
            )
          ),

          # ** Variable 1 -----------------------------------------------------
          fluidRow(
            column(
              4,
              uiOutput("var1_id.ui")
            ),
            column(
              4,
              selectInput(
                "var1_aggregate_by",
                label = h5("Aggregate by:"),
                choices = list("Max" = "max",
                               "Min" = "min",
                               "Mean" = "mean",
                               "Median" = "median"),
                selected = "max",
                width = 130
              )
            ),
            column(
              4,
              selectInput(
                "var1_relation", label = h5("Relationship:"),
                choices = list("Prot-Prot" = "protein",
                               "Prot-Spec" = "species"),
                selected = "protein",
                width = 130
              ),
              bsPopover(
                "var1_relation",
                "",
                "select if variable is the comparison between
                *seed protein - ortholog protein* or
                *seed protein - search taxon*",
                "top"
              )
            )
          ),

          # ** Variable 2 -----------------------------------------------------
          fluidRow(
            column(
              4,
              uiOutput("var2_id.ui")
            ),
            column(
              4,
              selectInput(
                "var2_aggregate_by",
                label = h5("Aggregate by:"),
                choices = list("Max" = "max",
                               "Min" = "min",
                               "Mean" = "mean",
                               "Median" = "median"),
                selected = "max",
                width = 130
              )
            ),
            column(
              4,
              uiOutput("var2_relation.ui")
            )
          ),

          hr(),

          # ** Domain input ---------------------------------------------------
          strong(h4("Additional annotation input:")),
          radioButtons(
            inputId = "anno_location", label = "",
            choices = list("from file", "from folder"),
            selected = "from file",
            inline = TRUE
          ),

          uiOutput("domain_input_file.ui"),

          hr(),
          em(
            a("Click here to download demo files",
              href = "https://github.com/BIONF/phyloprofile-data",
              target = "_blank")
          )
        ),

        # * 2nd column --------------------------------------------------------
        column(
          3,
          bsAlert("input_msg_ui"),

          # ** List of new taxa -----------------------------------------------
          conditionalPanel(
            condition = "output.unk_taxa_status == 'unknown' ||
                        output.unk_taxa_status == 'ncbi' ||
                        output.unk_taxa_status == 'invalid'",
            strong(h4("New taxa were found:")),
            dataTableOutput("unk_taxa_full"),
            br(),
            downloadButton("unk_taxa.download", "Download ID list")
          ),

          # ** Other input options --------------------------------------------
          conditionalPanel(
            condition = "output.unk_taxa_status == 0",
            strong(h4("Choose genes of interest:")),
            radioButtons(
              inputId = "gene_list_selected",
              label = "",
              choices = list("all", "from file"),
              selected = "all",
              inline = TRUE
            ),

            conditionalPanel(
              condition = "input.gene_list_selected == 'from file'",
              fileInput("list", "")
            ),

            hr(),

            checkboxInput(
              "ordering",
              strong("Order sequence IDs"),
              value = TRUE
            ),

            hr(),

            HTML("<b>Order taxa</b>"),

            radioButtons(
              inputId = "order_taxa",
              label = "",
              choices = list("automatically",
                             "by user defined tree"),
              selected = "automatically",
              inline = TRUE
            ),

            bsPopover("order_taxa", "", "in newick format", "bottom"),

            conditionalPanel(
              condition = "input.order_taxa == 'by user defined tree'",
              fileInput("inputTree", "")
            ),

            uiOutput("checkNewick.ui"),
            hr(),

            strong(h4("Other optional input:")),

            shinyBS::bsButton("fasta_upload", "FASTA file(s)"),
            h5(""),

            shinyBS::bsButton("upload_gene_category", "Gene categories"),
            h5(""),
            hr(),

            strong(h4("Color configuration:")),
            actionButton(
              "set_color",
              "Change colors",
              style = "padding:4px; font-size:100%"
            ),

            hr()
          )
        ),

        # * 3rd column --------------------------------------------------------
        column(
          4,

          # ** Msg for parsing new taxa ---------------------------------------
          conditionalPanel(
            condition = "output.unk_taxa_status == 'unknown' ||
                        output.unk_taxa_status == 'ncbi' ||
                        output.unk_taxa_status == 'invalid'",

            conditionalPanel(
              condition = "output.unk_taxa_status == 'invalid'",
              HTML(
                "<p><em>Some new taxa have <span style=\"color: #ff0000;\">
                invalid IDs</span> (either in newTaxa.txt or in the main
                profile input or both). IDs of non-NCBI taxa have to be greater
                than 2268208.</em></p>
                <p><em>Please replace those IDs before continuing!</em></p>"
              )
            ),

            conditionalPanel(
              condition = "output.unk_taxa_status == 'unknown'",
              HTML(
                '<p><em>NCBI taxonomy information of some taxa can neither
                </em></p>
                <ul>
                <li><em>be retrieved from NCBI (<span style="color: #0000ff;
                ">Source="ncbi"</span>) nor </em></li>
                <li><em>be found in <span style="color: #ff0000;">
                phyloprofile/data/newTaxa.txt</span>&nbsp;(<span
                style="color: #0000ff;">Source="new"</span>) file</em></li>
                </ul>
                <p><strong><em>Please add taxonomy information for those
                unknown taxa and <span style="color: #ff0000;">
                reload the tool</span> to continue!</em></strong></p>'
              ),
              h5(""),
              shinyBS::bsButton(
                "add_taxa",
                "Add taxonomy info",
                disabled = FALSE,
                style = "warning"
              )
            ),

            conditionalPanel(
              condition = "output.unk_taxa_status == 'ncbi'",
              HTML(
                '<p><em>NCBI taxonomy information of some taxa can either
                </em></p>
                <ul>
                <li><em>be retrieved from NCBI (<span style="color: #0000ff;
                ">Source="ncbi"</span>) or </em></li>
                <li><em>be found in <span style="color: #ff0000;">
                phyloprofile/data/newTaxa.txt</span>&nbsp;(<span
                style="color: #0000ff;">Source="new"</span>) file</em></li>
                </ul>
                <p><strong><em>Click here to get required taxonomy information
                for those taxa!</em></strong></p>'
              ),
              h5(""),
              shinyBS::bsButton(
                "but_parse",
                "Get taxonomy info",
                disabled = FALSE,
                style = "warning"
              ),

              hr(),
              uiOutput("end_parsing_msg"),
              tableOutput("invalidID.output"),
              hr(),
              conditionalPanel(
                condition = "output.unk_taxa_status == 'invalid'",
                downloadButton("invalidID.download", "Download invalid IDs")
              )
            )
          ),

          # ** List of ranks & available taxa ---------------------------------
          conditionalPanel(
            condition = "output.unk_taxa_status == 0",
            strong(h4("Seed (super)taxon:")),
            br(),

            strong(h5("Select taxonomy rank:")),
            withSpinner(
              uiOutput("rank_select"),
              proxy.height = "50px",
              type = 7,
              size = 0.5
            ),
            br(),

            strong(h5("Choose (super)taxon of interest:")),
            withSpinner(
              uiOutput("select"),
              proxy.height = "50px",
              type = 7,
              size = 0.5
            ),
            br(),

            shinyBS::bsButton(
              "do",
              "PLOT",
              type = "action",
              style = "danger",
              size = "large",
              disabled = FALSE
            ),
            h5("")
          )
        )
      ),

      # MAIN PROFILE TAB ======================================================
      tabPanel(
        "Main profile",
        sidebarLayout(
          # * sidebar panel for profile highlight -----------------------------
          sidebarPanel(
            uiOutput("total_gene_number.ui"),

            column(
              4,
              numericInput(
                "st_index",
                "Show from:",
                min = 1,
                max = 1600,
                value = 1,
                width = 100
              ),
              style = "padding:0px;"
            ),

            column(
              4,
              numericInput(
                "end_index",
                "...to:",
                min = 1,
                max = 1600,
                value = 30,
                width = 100
              ),
              style = "padding:0px;"
            ),

            column(
              4,
              uiOutput("highlight_gene_ui")
            ),

            bsPopover(
              "highlight_gene_ui",
              "",
              "Select gene to highlight",
              "bottom"
            ),

            bsPopover(
              "st_index",
              "",
              "Set start index for sequence range",
              "bottom"
            ),

            bsPopover(
              "end_index",
              "",
              "Set end index for sequence range",
              "bottom"
            ),

            br(),

            uiOutput("highlight_taxon_ui"),

            checkboxInput(
              "color_by_group",
              strong("Highlight genes by categories"),
              value = FALSE
            ),

            conditionalPanel(
              condition = "input.auto_update == false",
              shinyBS::bsButton(
                "update_btn",
                "Update plot",
                style = "warning"
              )
            )
          ),
          # * main panel for profile plot -------------------------------------
          mainPanel(
            conditionalPanel(
              condition = "input.do > 0",
              create_profile_plot_ui("main_profile")
            )
            # ,
            # conditionalPanel(
            #   condition = "input.main_x_axis_guide == true |
            #   input.main_y_axis_guide == true",
            #   absolutePanel(
            #     id = "absAxis",
            #     bottom = 0, left = 0,
            #     heigh = NULL, width = NULL,
            #     fixed = TRUE,
            #     draggable = TRUE,
            #     style = "opacity: 0.80",
            #
            #     uiOutput("mainAxisRender")
            #   )
            # )
          )
        )
      ),

      # CUSTOMIZED PROFILE TAB ================================================
      tabPanel(
        "Customized profile",
        sidebarLayout(
          # * sidebar panel for subseting data --------------------------------
          sidebarPanel(
            width = 4,
            column(
              12,
              style = "padding:0px;",
              strong("Select sequence(s) of interest:")
            ),

            column(
              12,
              fluidRow(
                column(
                  8,
                  style = "padding:0px;",
                  uiOutput("gene_in")
                ),
                column(
                  4,
                  fileInput("custom_file", "", width = "100%")
                )
              )
            ),

            column(
              12,
              style = "padding:0px;",
              strong("Select (super)taxon/(super)taxa of interest:")
            ),
            column(
              12,
              fluidRow(
                column(
                  8,
                  style = "padding:0px;",
                  uiOutput("taxa_in")
                ),
                column(
                  4,
                  h3(""),
                  shinyBS::bsButton("cus_taxa", "Browse...")
                )
              )
            ),

            h5(""),
            shinyBS::bsButton(
              "plot_custom",
              "Update selected sequence(s)/taxa",
              style = "warning"
            )
          ),

          # * main panel for customized profile plot --------------------------
          mainPanel(
            conditionalPanel(
              condition = "output.same_profile == true",
              h4(
                "Please select subset of genes and/
                or taxa for customized profile!"
              )
            ),

            conditionalPanel(
              condition = "input.do > 0",
              create_profile_plot_ui("customized_profile")
            )
          )
        )
      ),

      # FUNCTION TAB ==========================================================
      navbarMenu(
        "Function",
        # * Profiles clustering -----------------------------------------------
        tabPanel(
          "Profiles clustering",
          h4(strong("Profiles clustering")),
          bsAlert("desc_clustering_ui"),

          wellPanel(
            fluidRow(
              column(
                2,
                uiOutput("select_profile_type")
              ),
              column(
                3,
                uiOutput("select_dist_method")
              ),

              column(
                3,
                selectInput(
                  "cluster_method", label = h5("Cluster method:"),
                  choices = list("single" = "single",
                                 "complete" = "complete",
                                 "average (UPGMA)" = "average",
                                 "mcquitty (WPGMA)" = "mcquitty",
                                 "median (WPGMC)" = "median",
                                 "centroid (UPGMC)" = "centroid"),
                  selected = "complete"
                )
              ),

              column(
                1,
                create_plot_size("cluster_plot.width", "Width (px)", 600)
              ),
              column(

                1,
                create_plot_size("cluster_plot.height", "Height (px)", 600)
              ),
              column(
                2,
                checkboxInput(
                  "apply_cluster",
                  em(strong("Apply clustering to profile plot",
                            style = "color:darkblue")),
                  value = FALSE
                ),

                uiOutput("apply_cluster_check.ui"),

                checkboxInput(
                  "add_cluster_cutom_profile",
                  strong(em("Add selected genes to Customized profile",
                            style = "color:red")),
                  value = FALSE,
                  width = NULL
                ),
                uiOutput("add_cluster_cutom_profile_check.ui")
              )
            )
          ),

          cluster_profile_ui("profile_clustering")
        ),

        # * Distribution analysis ---------------------------------------------
        tabPanel(
          "Distribution analysis",
          h4(strong("Distribution analysis")),
          bsAlert("desc_distribution_ui"),

          wellPanel(
            fluidRow(
              column(
                2,
                selectInput(
                  "dataset.distribution", "Select data",
                  choices = c("Main data", "Customized data"),
                  selected = "Main data"
                ),
                uiOutput("selected.distribution")
              ),
              column(
                2,
                uiOutput("var1_dist.ui")
              ),
              column(
                2,
                uiOutput("var2_dist.ui")
              ),
              column(
                2,
                uiOutput("percent_dist.ui")
              ),
              column(
                2,
                create_text_size("dist_text_size", "Label size", 12, 100)
              ),
              column(
                2,
                create_plot_size("dist_width", "Width (px)", 600)
              )
            )
          ),
          analyze_distribution_ui("dist_plot")
        ),


        # * Gene age estimation -----------------------------------------------
        tabPanel(
          "Gene age estimation",
          h4(strong("Gene age estimation")),
          bsAlert("desc_gene_age_ui"),

          wellPanel(
            fluidRow(
              column(
                2,
                uiOutput("var1_age.ui")
              ),
              column(
                2,
                uiOutput("var2_age.ui")
              ),
              column(
                2,
                uiOutput("percent_age.ui")
              ),
              column(
                2,
                strong("Appearance"),
                bsButton("gene_age_prot_config", "Plot config")
              ),
              column(
                4,
                checkboxInput(
                  "add_gene_age_custom_profile",
                  strong(em("Add selected genes to Customized profile",
                            style = "color:red")),
                  value = FALSE,
                  width = NULL
                ),
                uiOutput("add_gene_age_custom_profile_check.ui")
              )
            )
          ),
          plot_gene_age_ui("gene_age")
        ),

        # * Core gene identification  -----------------------------------------
        tabPanel(
          "Core gene identification",
          h4(strong("Core gene identification")),
          bsAlert("desc_core_gene_ui"),

          wellPanel(
            fluidRow(
              column(
                3,
                uiOutput("var1_core.ui")
              ),
              column(
                3,
                uiOutput("var2_core.ui")
              ),
              column(
                3,
                uiOutput("percent_core.ui")
              ),
              column(
                3,
                sliderInput(
                  "core_coverage",
                  "Core taxa coverage",
                  min = 0,
                  max = 100,
                  value = 100,
                  step = 5
                )
              ),
              column(
                12,
                uiOutput("taxa_list_core.ui"),
                shinyBS::bsButton("browse_taxa_core", "Browse")
              )
            )
          ),
          hr(),

          column(
            4,
            downloadButton("core_gene_table_download", "Download gene list"),
            checkboxInput(
              "add_core_gene_custom_profile",
              strong(em("Add core genes to Customized profile",
                        style = "color:red")),
              value = FALSE,
              width = NULL
            ),
            uiOutput("add_core_gene_custom_profile_check.ui")
          ),
          identify_core_gene_ui("core_gene")
        ),

        # * Group Comparison  -------------------------------------------------
        tabPanel(
          "Group comparison",
          h4(strong("Group comparison")),
          bsAlert("desc_gc_ui"),
          wellPanel(
            fluidRow(
              column(
                3,
                uiOutput("variable_button_gc"),
                popify(
                  checkboxInput(
                    "right_format_features",
                    "Annotation format:
                    ’Type_Name’",
                    value = TRUE,
                    width = NULL
                  ),
                  "",
                  "E.g.: pfam_ApbA, smart_SRP54"
                )
              ),
              column(
                2,
                uiOutput("list_genes_gc"),
                popify(
                  fileInput("gc_file", NULL, width = "100%"),
                  "",
                  "Upload list of genes of interest"
                )
              ),
              column(
                2,
                uiOutput("taxa_list_gc"), # Select In-Group
                shinyBS::bsButton("taxa_gc", "Browse"),
                checkboxInput(
                  "use_common_ancestor",
                  "Use common ancestor",
                  value = TRUE,
                  width = NULL
                ),
                bsPopover(
                  "use_common_ancestor",
                  "",
                  "All taxa that have the same common ancestor with
                  the selected taxa above will be considered as the in-group",
                  "top"
                )
              ),
              column(
                3,
                uiOutput("significance.ui"),
                checkboxInput(
                  "add_gc_genes_custom_profile",
                  strong(em("Add candidate gene(s) to Customized profile",
                            style = "color:red")),
                  value = FALSE,
                  width = NULL
                ),
                uiOutput("add_gc_custom_profile_check")
              ),
              column(
                2,
                popify(
                  actionButton("gc_plot_config", "Plot config"),
                  "",
                  "Change the appearance of the plots"
                ),
                hr(),
                bsButton("plot_gc", "COMPARE!", style = "warning")
              )
            )
          ),
          group_comparison_ui("group_comparison")
        ),

        # * Search for NCBI taxonomy IDs  -------------------------------------
        search_taxon_id_ui("search_taxon_id")
      ),

      # DATA DOWNLOAD TAB =====================================================
      navbarMenu(
        "Download filtered data",
        download_filtered_main_ui("filtered_main_download"),
        download_filtered_customized_ui("filtered_customized_download")
      ),

      # HELP TAB ==============================================================
      navbarMenu(
        "Help",
        tabPanel(
          a(
            "Wiki",
            href = "https://github.com/BIONF/PhyloProfile/wiki",
            target = "_blank"
          )
        ),
        tabPanel(
          a(
            "About",
            href = "https://BIONF.github.io/PhyloProfile/",
            target = "_blank"
          )
        )
      )
    ),

    # LIST OF POP-UP WINDOWS ==================================================

    # * popup for getting taxa from OMA browser -------------------------------
    bsModal(
      "get_oma_data_windows",
      "Get OMA data",
      "open_oma_windows",
      size = "small",

      selectInput(
        "selected_oma_type",
        label = "Select type of OMA orthologs:",
        choices = list("HOG", "OG", "PAIR"),
        selected = "HOG"
      ),
      shinyBS::bsButton("get_data_oma", "Get data", style = "danger"),
      downloadButton("download_files_oma", "Save data"),
      br(),
      em("This windows will close automatically when eveything
           is done!", style = "color:red")
    ),

    # * popup for adding new taxa from input file -----------------------------
    bsModal(
      "add_taxa_windows",
      "Add new taxa",
      "add_taxa",
      size = "medium",

      HTML(
        "<p><em>Use this form to add taxon that does not exist in NCBI taxonomy
         database (or alternatively you can manually prepare the
        <span style=\"text-decoration: underline;\">
        <span style=\"color: #ff0000; text-decoration: underline;\">
        phyloprofile/data/newTaxa.txt file with the following description
        for each field).</em></p>
        <p><span style=\"color: #ff0000;\"><em><strong>
        NOTE: ID and name of new taxon must be
        <span style=\"text-decoration: underline;\">
        different</span> from any existing NCBI taxa.</strong></em></span></p>"
      ),

      textInput(
        "new_id",
        "ID (must be a number and greater than 2268208,
        e.g. 9000001)",
        9000001,
        width = 500
      ),
      textInput(
        "new_name",
        "Name (e.g. Saccharomyces cerevisiae strain ABC)",
        "",
        width = 500
      ),
      textInput(
        "new_rank",
        "Rank (e.g. \"norank\" (for strain), species, order, etc.)",
        "norank",
        width = 500
      ),
      textInput(
        "new_parent",
        "Parent ID (NCBI taxonomy ID of the next higher rank,
      e.g. 4932 (S.cerevisiae species))",
        4932,
        width = 500
      ),
      actionButton("new_add", "Add new taxon"),

      hr(),
      fileInput("new_taxa_file","Or upload file contains IDs for new taxa"),
      HTML(
        "<p><em>Taxonomy file for new taxa has to be a tab-delimited text file
        and has the following header (please follow the rule above):</em></p>
        <p>ncbiID &nbsp;fullName &nbsp;rank &nbsp;parentID</p>"
      ),
      bsAlert("wrong_new_taxa"),

      hr(),
      shinyBS::bsButton("new_done", "Finish adding", style = "warning",
                        disabled = TRUE)
    ),

    # * popup for confirming parsing taxa from input file ---------------------
    bsModal(
      "parse_confirm",
      "Get taxonomy info",
      "but_parse",
      size = "small",

      HTML(
        '<p>Fetching Missing Taxonomy Information and Post-processing.</p>
        <p><em>This windows will close automatically when eveything is done.
        Please wait...</em></p>
        <p><strong><span style="color: #ff0000;">PLEASE RELOAD THIS TOOL WHEN
        FINISHED!!!</span></strong></p>'
      )
    ),

    # * popup for plotting detailed plot --------------------------------------
    bsModal(
      "modal_bs",
      "Detailed plot",
      "detailed_btn",
      size = "large",

      fluidRow(
        column(
          2,
          create_plot_size("detailed_height", "Height (px)", 100)
        ),
        column(
          3,
          create_text_size("detailed_text", "Text size (px)", 12, 150)
        ),
        column(
          7,
          checkboxInput(
            "detailed_remove_na",
            strong("Hide taxa that have no ortholog (NAs)",
                   style = "color:red"),
            value = FALSE
          )
        )
      ),
      hr(),

      create_detailed_plot_ui("detailed_plot"),

      shinyBS::bsButton("do_domain_plot",
                        "Show domain architecture",
                        disabled = TRUE),
      uiOutput("check_domain_files"),
      br(),

      h4("Sequence:"),
      verbatimTextOutput("fasta")
    ),

    # * popup for plotting domain architecture plot ---------------------------
    bsModal(
      "plot_archi",
      "Domain architecture",
      "do_domain_plot",
      size = "large",

      fluidRow(
        column(
          2,
          create_plot_size("archi_height", "Plot height(px)", 400)
        ),
        column(
          2,
          create_plot_size("archi_width", "Plot width(px)", 800)
        ),
        column(
          2,
          create_text_size("title_archi_size", "Title size(px)", 11, 150)
        ),
        column(
          2,
          create_text_size("label_archi_size", "SeqID size(px)", 11, 150)
        )
      ),
      uiOutput("test.ui"),
      create_architecture_plot_ui("archi_plot")
    ),

    # * popup for setting plot colors (profiles) ------------------------------
    bsModal(
      "color",
      "Set colors for profile",
      "set_color",
      size = "small",

      colourpicker::colourInput(
        "low_color_var1",
        "Low variable 1 (dot)",
        value = "darkorange"
      ),
      colourpicker::colourInput(
        "high_color_var1",
        "High variable 1 (dot)",
        value = "steelblue"
      ),
      actionButton(
        "default_color_var1",
        "Default",
        style = "padding:4px; font-size:100%"
      ),
      hr(),

      colourpicker::colourInput(
        "low_color_var2",
        "Low variable 2 (background)",
        value = "grey95"
      ),
      colourpicker::colourInput(
        "high_color_var2",
        "High variable 2 (background)",
        value = "khaki"
      ),
      actionButton(
        "default_color_var2",
        "Default",
        style = "padding:4px; font-size:100%"
      ),
      hr(),

      colourpicker::colourInput(
        "para_color",
        "Color for inparalogs",
        value = "#07d000"
      ),
      actionButton(
        "default_color_para",
        "Default",
        style = "padding:4px; font-size:100%"
      )
    ),

    # * popup for FASTA upload ------------------------------------------------
    bsModal(
      "fasta_upload_bs",
      "FASTA upload",
      "fasta_upload",
      size = "small",

      selectInput(
        "input_type", "Choose location for:",
        c("Concatenated fasta file", "Fasta folder")
      ),
      hr(),

      uiOutput("default_color_para.ui"),

      conditionalPanel(
        condition = "input.input_type == 'Concatenated fasta file'",
        fileInput("concat_fasta", ""),
        uiOutput("concat_fasta.exist_check")
      ),
      conditionalPanel(
        condition = "input.input_type == 'Fasta folder'",
        textInput("path", "Main FULL path:", ""),
        selectInput(
          "dir_format", "Directory format:",
          choices = list("path/speciesID.fa*" = 1,
                         "path/speciesID/speciesID.fa*" = 2),
          selected = "Path/speciesID.fasta"
        ),
        selectInput(
          "file_ext", "File extension:",
          choices = list("fa" = "fa",
                         "fasta" = "fasta",
                         "fas" = "fas",
                         "txt" = "txt"),
          selected = "fa"
        ),
        selectInput(
          "id_format",
          "ID format:",
          choices = list(">speciesID:seqID" = 1,
                         ">speciesID@seqID" = 2,
                         ">speciesID|seqID" = 3,
                         ">seqID" = 4),
          selected = 4
        )
      )
    ),

    # * popup for upload gene category ----------------------------------------
    bsModal(
      "upload_gene_category_bs",
      "Upload gene categories",
      "upload_gene_category",
      size = "small",
      fileInput("gene_category", "")
    ),

    # * popup for setting Main plot configurations ----------------------------
    bsModal(
      "main_plot_config_bs",
      "Plot appearance configuration",
      "main_plot_config",
      size = "small",

      column(
        6,
        create_text_size("x_size", "X-axis label size (px)", 8, 100)
      ),
      column(
        6,
        create_text_size("y_size", "Y-axis label size (px)", 8, 100)
      ),

      column(
        6,
        create_text_size("legend_size", "Legend label size (px)", 8, 150)
      ),
      column(
        6,
        selectInput(
          "main_legend", label = "Legend position:",
          choices = list("Right" = "right",
                         "Left" = "left",
                         "Top" = "top",
                         "Bottom" = "bottom",
                         "Hide" = "none"),
          selected = "right",
          width = 150
        )
      ),
      column(
        12,
        HTML("<strong>Angle for x-axis label</strong>:<br>"),
        sliderInput(
          "x_angle",
          "",
          min = 0,
          max = 90,
          step = 10,
          value = 60,
          width = 250
        ),
        br()
      ),

      column(
        12,
        HTML("<strong>Zooming factor (α) for dots on
             profile</strong>:<br>"),
        sliderInput(
          "dot_zoom", "",
          min = -1,
          max = 3,
          step = 0.1,
          value = 0,
          width = 250
        ),
        HTML("<em>dot size = (1+α)*default_size<br>default_size
             =[0:5]</em>"),
        uiOutput("dot_size_info"),
        br()
      ),

      br(),
      hr(),
      shinyBS::bsButton("reset_main_config", "Reset", style = "danger"),
      shinyBS::bsButton("applyMainConfig", "Done", style = "warning")
    ),

    # * popup for setting Customized plot configurations ----------------------
    bsModal(
      "selected_plot_config_bs",
      "Plot appearance configuration",
      "selected_plot_config",
      size = "small",

      column(
        6,
        create_text_size("x_size_select", "X-axis label size (px)", 8, 100)
      ),
      column(
        6,
        create_text_size("y_size_select", "Y-axis label size (px)", 8, 100)
      ),

      column(
        6,
        create_text_size("legend_size_select", "Legend label size (px)", 8, 150)
      ),
      column(
        6,
        selectInput(
          "selected_legend", label = "Legend position:",
          choices = list("Right" = "right",
                         "Left" = "left",
                         "Top" = "top",
                         "Bottom" = "bottom",
                         "Hide" = "none"),
          selected = "right",
          width = 150
        )
      ),
      column(
        12,
        HTML("<strong>Angle for x-axis label</strong>:<br>"),
        sliderInput(
          "x_angle_select", "",
          min = 0,
          max = 90,
          step = 10,
          value = 60,
          width = 250
        ),
        br()
      ),

      column(
        12,
        HTML("<strong>Zooming factor (α) for dots on profile</strong>:<br>"),
        sliderInput(
          "dot_zoom_select", "",
          min = -1,
          max = 3,
          step = 0.1,
          value = 0,
          width = 250
        ),
        HTML("<em>dot size = (1+α)*default_size<br>default_size=[0:5]</em>"),
        uiOutput("dot_size_infoSelect"),
        br()
      ),

      br(),
      hr(),
      shinyBS::bsButton("reset_selected_config", "Reset", style = "danger"),
      shinyBS::bsButton("apply_selected_config", "Done", style = "warning")
    ),

    # * popup for setting Gene age plot configurations ------------------------
    bsModal(
      "gene_age_prot_config_bs",
      "Plot appearance configuration",
      "gene_age_prot_config",
      size = "small",

      sliderInput(
        "gene_age_width",
        "Width zoom (*600px)",
        min = 0,
        max = 5,
        step = 0.1,
        value = 1,
        width = "100%"
      ),
      sliderInput(
        "gene_age_height", "Height zoom (*150px)",
        min = 0,
        max = 5,
        step = 0.1,
        value = 1,
        width = "100%"
      ),
      sliderInput(
        "gene_age_text", "Text size zoom",
        min = 0,
        max = 5,
        step = 0.1,
        value = 1,
        width = "100%"
      ),
      br(),
      hr(),
      shinyBS::bsButton(
        "reset_gene_age_prot_config",
        "Reset",
        style = "danger"
      )
    ),

    # * popup for setting Group compariosn plot configurations ----------------
    bsModal(
      "gc_plot_config_bs",
      "Plot appearance configuration",
      "gc_plot_config",
      size = "small",

      column(
        6,
        create_text_size("x_size_gc", "X-axis label size (px)", 10, 100)
      ),
      column(
        6,
        create_text_size("y_size_gc", "Y-axis label size (px)", 10, 100)
      ),
      column(
        6,
        create_text_size("legend_size_gc", "Legend label size (px)", 10, 150)
      ),
      column(
        6,
        selectInput(
          "legend_gc", label = "Legend position:",
          choices = list("Right" = "right",
                         "Left" = "left",
                         "Top" = "top",
                         "Bottom" = "bottom",
                         "Hide" = "none"),
          selected = "right",
          width = 150
        )
      ),
      column(
        6,
        create_text_size("p_values_size_gc", "P-value label size (px)", 10, 100)
      ),
      column(
        6,
        selectInput(
          "show_point_gc", label = "Show location parameter:",
          choices = list("Mean" = "mean",
                         "Median" = "median"),
          selected = "mean",
          width = 150)
        ),
      column(
        12,
        sliderInput(
          "angle_gc", "Angle of the X-axis label",
          min = 0,
          max = 180,
          step = 1,
          value = 90,
          width = 250
        )
      ),
      column(
        12,
        checkboxInput(
          "show_p_value",
          strong("Show P-Values"),
          value = TRUE,
          width = 250
        )
      ),
      column(
        12,
        popify(
          checkboxInput(
            "highlight_significant",
            strong("Highlight significant plots"),
            value = TRUE,
            width = 250
          ),
          "",
          "If both variables are selected the significant Plot is colored"
        )
      ),

      br(),
      hr(),
      shinyBS::bsButton("reset_config_gc", "Reset", style = "danger"),
      shinyBS::bsButton("apply_config_gc", "Done", style = "warning")
    ),

    # * popup for select taxa on Customized Profile ---------------------------
    bsModal(
      "cus_taxa_bs",
      "Select taxon/taxa of interest",
      "cus_taxa",
      size = "small",

      select_taxon_rank_ui("select_taxon_rank"),
      checkboxInput(
        "apply_cus_taxa",
        strong("Apply to customized profile",
               style = "color:red"),
        value = FALSE
      )
    ),

    # * popup for select taxa on Core gene finding ----------------------------
    bsModal(
      "browse_taxa_core_bs",
      "Select taxon/taxa of interest",
      "browse_taxa_core",
      size = "small",

      select_taxon_rank_ui("select_taxon_rank_core"),
      checkboxInput(
        "apply_core_taxa",
        strong("Apply", style = "color:red"),
        value = FALSE
      )
    ),

    # * popup for select taxa on Group comparison -----------------------------
    bsModal(
      "taxa_gc_bs",
      "Select taxon/taxa of interest",
      "taxa_gc",
      size = "small",

      select_taxon_rank_ui("select_taxon_rank_gc"),
      checkboxInput(
        "apply_taxa_gc",
        strong("Apply",
               style = "color:red"),
        value = FALSE
      )
    ),

    # POINT INFO BOX ==========================================================
    conditionalPanel(
      condition =
        "input.tabs=='Main profile' || input.tabs=='Customized profile'",

      absolutePanel(
        bottom = 5, left = 30,
        fixed = TRUE,
        draggable = TRUE,
        h5("Point's info:"),
        verbatimTextOutput("point_info"),
        conditionalPanel(
          condition = "output.point_info_status == 0",
          shinyBS::bsButton(
            "detailed_btn",
            "Detailed plot",
            style = "success",
            disabled = FALSE
          )
        ),
        style = "opacity: 0.80"
      )
    )
  )
)
