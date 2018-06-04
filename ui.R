# pacman NOT YET WORK WITH shinyapp.io ========================================
#if(!("pacman" %in% installed.packages())) install.packages("pacman")
# library(pacman)
# p_load(shiny, shinyBS, shinyjs, DT, colourpicker, install = T)
# -----------------------------------------------------------------------------
# packages <- c("shiny", "shinyBS", "shinyjs", "DT", "colourpicker")
# sapply(packages, require, character.only = TRUE)
# -----------------------------------------------------------------------------

if (!require("shiny")) install.packages("shiny")
if (!require("shinyBS")) install.packages("shinyBS")
if (!require("DT")) install.packages("DT")
if (!require("colourpicker")) install.packages("colourpicker")
if (!require("shinyjs")) install.packages("shinyjs")
if (!require("shinycssloaders")) {
  if ("devtools" %in% installed.packages() == FALSE){
    install.packages("devtools")
  }
  devtools::install_github("andrewsali/shinycssloaders")
}

source("scripts/search_taxon_id.R")
source("scripts/download_filtered_main.R")
source("scripts/download_filtered_customized.R")

source("scripts/select_taxon_rank.R")

source("scripts/identify_core_gene.R")
source("scripts/analyze_distribution.R")

# MAIN UI =====================================================================

shinyUI(fluidPage(
  tags$style(type = "text/css", "body {padding-top: 80px;}"),
  
  # Application title
  titlePanel("", windowTitle = "PhyloProfile"),
  useShinyjs(),
  
  # TOP wellpanel for plot configuration --------------------------------------
  conditionalPanel(
    condition = "input.tabs=='Main profile'",
    wellPanel(
      fluidRow(
        column(2,
               radioButtons(
                 inputId = "x_axis",
                 label = "Choose type of x-axis:",
                 choices = list("taxa", "genes"),
                 selected = "taxa",
                 inline = T
               ),
               hr(),
               checkboxInput("auto_update",
                             strong(em("Auto update plot")),
                             value = TRUE,
                             width = NULL)
        ),
        column(1,
               numericInput("width",
                            "Width (px)",
                            min = 600,
                            max = 3200,
                            step = 50,
                            value = 600,
                            width = 100),
               actionButton("main_plot_config", "Appearance")
        ),
        column(1,
               numericInput("height",
                            "Height (px)",
                            min = 600,
                            max = 1600,
                            step = 50,
                            value = 600,
                            width = 100)
        ),
        column(2,
               uiOutput("var1_cutoff.ui")
        ), column(2,
                  uiOutput("var2_cutoff.ui")
        ),
        column(2,
               sliderInput("percent",
                           "% of present taxa:",
                           min = 0,
                           max = 1,
                           step = 0.025,
                           value = c(0.0, 1.0),
                           width = 200)
        ),
        column(2,
               shinyBS::bsButton("reset_main",
                                 "Reset cutoffs",
                                 style = "danger"),
               hr(),
               downloadButton("plot_download", "Download profile"),
               tags$head(
                 tags$style(HTML("#plot_download{background-color:#A9E2F3}"))
               )
        )
      )
    )
  ),
  
  conditionalPanel(
    condition = "input.tabs=='Customized profile'",
    wellPanel(
      fluidRow(
        column(2,
               radioButtons(
                 inputId = "x_axis_selected",
                 label = "Choose type of x-axis:",
                 choices = list("taxa", "genes"),
                 selected = "taxa",
                 inline = T
               ),
               hr(),
               checkboxInput("auto_update_selected",
                             strong(em("Auto update plot")),
                             value = TRUE,
                             width = NULL)
        ),
        column(1,
               numericInput("selected_width",
                            "Width (px)",
                            min = 100,
                            max = 1000,
                            step = 50,
                            value = 600,
                            width = 100),
               actionButton("selected_plot_config", "Appearance")
        ),
        column(1,
               numericInput("selected_height",
                            "Height (px)",
                            min = 100,
                            max = 1600,
                            step = 50,
                            value = 600,
                            width = 100)
        ),
        column(2,
               uiOutput("var1_filter.ui")
        ),
        column(2,
               uiOutput("var2_filter.ui")
        ),
        column(2,
               uiOutput("percent_filter.ui")
        ),
        column(2,
               shinyBS::bsButton("reset_selected",
                                 "Reset cutoffs",
                                 style = "danger"),
               hr(),
               downloadButton("selected_download", "Download profile"),
               tags$head(
                 tags$style(
                   HTML("#selected_download{background-color:#A9E2F3}"))
               )
        )
      )
    )
  ),
  
  # main narvarpage tabs ------------------------------------------------------
  navbarPage(
    em(strong("PhyloProfile v0.3.0-beta")),
    id = "tabs",
    collapsible = TRUE,
    inverse = TRUE,
    fluid = TRUE,
    position = "fixed-top",
    
    # INPUT TAB ---------------------------------------------------------------
    tabPanel(
      "Input & settings",
      column(4,
             strong(h4("Main input:")),
             conditionalPanel(
               condition = "input.do",
               em(strong("RELOAD THIS TOOL TO UPLOAD A NEW INPUT FILE!!!",
                         style = "color:red"))
             ),
             
             # checkboxInput("demo",
             #               em(strong("Use demo files"),
             #               style = "color:darkblue")),
             selectInput("demo_data", label = h5("Use online demo data:"),
                         choices = list("None" = "none",
                                        "AMPK-TOR" = "ampk-tor",
                                        "LCA Microsporidia" = "lca-micros"),
                         selected = "none",
                         width = "80%"),
             
             uiOutput("no_internet_msg"),
             uiOutput("demo_data_describe"),
             uiOutput("main_input_file.ui"),
             uiOutput("input_check.ui"),
             
             fluidRow(
               column(6, uiOutput("select_oma_type")),
               column(6, uiOutput("button_oma"),
                      uiOutput("oma_download"))
             ),
             
             fluidRow(
               column(4,
                      uiOutput("var1_id.ui")
               ),
               column(4,
                      selectInput("var1_aggregate_by",
                                  label = h5("Aggregate by:"),
                                  choices = list("Max" = "max",
                                                 "Min" = "min",
                                                 "Mean" = "mean",
                                                 "Median" = "median"),
                                  selected = "max",
                                  width = 130)
               ),
               column(4,
                      selectInput("var1_relation", label = h5("Relationship:"),
                                  choices = list("Prot-Prot" = "protein",
                                                 "Prot-Spec" = "species"),
                                  selected = "protein",
                                  width = 130),
                      bsPopover("var1_relation",
                                "",
                                "select if variable is the comparison between
                                *seed protein - ortholog protein* or
                                *seed protein - search taxon*",
                                "top")
               )
               ),
             fluidRow(
               column(4,
                      uiOutput("var2_id.ui")
               ),
               column(4,
                      selectInput("var2_aggregate_by",
                                  label = h5("Aggregate by:"),
                                  choices = list("Max" = "max",
                                                 "Min" = "min",
                                                 "Mean" = "mean",
                                                 "Median" = "median"),
                                  selected = "max",
                                  width = 130)
               ),
               column(4,
                      uiOutput("var2_relation.ui")
                      # selectInput("var2_relation",
                      #              label = h5("Relationship:"),
                      #              choices = list("Prot-Prot"="protein",
                      #                       "Prot-Spec"="species"),
                      #             selected = "species",
                      #             width = 130)
               )
             ),
             
             hr(),
             strong(h4("Additional annotation input:")),
             radioButtons(inputId = "anno_location", label = "",
                          choices = list("from file", "from folder"),
                          selected = "from file",
                          inline = T),
             uiOutput("domain_input_file.ui"),
             
             hr(),
             em(a("Click here to download demo files",
                  href = "https://github.com/BIONF/phyloprofile-data",
                  target = "_blank"))
             # em("Click here to download demo files:"),
             # em(a("(1) Main inputs,",
             #       href={ "https://github.com/BIONF/phyloprofile-data/
             #           tree/master/demo",
             #      target="_blank")),
             # em(a("(2) Domain annotations (optional),",
             #      href="https://github.com/BIONF/phyloprofile-data/
             #            tree/master/demo/domain_files",
             #      target="_blank")),
             # em(a("(3) FASTA sequence files (optional)",
             #      href="https://github.com/BIONF/phyloprofile-data/
             #            tree/master/demo/fasta_files",
             #      target="_blank"))
               ),
      column(3,
             conditionalPanel(
               condition = "output.unk_taxa_status == 1",
               strong(h4("New taxa were found:")),
               dataTableOutput("unk_taxa_full")
             ),
             
             conditionalPanel(
               condition = "output.unk_taxa_status == 0",
               strong(h4("Choose genes of interest:")),
               radioButtons(inputId = "gene_list_selected",
                            label = "",
                            choices = list("all", "from file"),
                            selected = "all",
                            inline = T),
               conditionalPanel(
                 condition = "input.gene_list_selected == 'from file'",
                 fileInput("list", "")
               ),
               hr(),
               
               checkboxInput("ordering",
                             strong("Order sequence IDs"),
                             value = TRUE),
               hr(),
               
               HTML("<b>Order taxa</b>"),
               radioButtons(inputId = "order_taxa",
                            label = "",
                            choices = list("automatically",
                                           "by user defined tree"),
                            selected = "automatically",
                            inline = T),
               bsPopover("order_taxa", "", "in newick format", "bottom"),
               
               conditionalPanel(
                 condition = "input.order_taxa == 'by user defined tree'",
                 fileInput("inputTree", "")
               ),
               uiOutput("checkNewick.ui"),
               hr(),
               
               shinyBS::bsButton("get_config", "FASTA config"),
               h5(""),
               actionButton("set_color",
                            "COLORS config",
                            style = "padding:4px; font-size:100%"),
               hr()
             )
      ),
      column(4,
             conditionalPanel(
               condition = "output.unk_taxa_status",
               strong(h4("PLEASE CHECK:")),
               
               em("Do you have any taxon,
                  which doesn't exist in the NCBI taxonomy database?"),
               radioButtons("new_taxa_ask",
                            "",
                            c("Yes" = "Yes", "No" = "No"),
                            inline = T,
                            selected = "No"),
               conditionalPanel(condition = "input.new_taxa_ask == 'Yes'",
                                shinyBS::bsButton("add_taxa",
                                                  "Add info for new taxa",
                                                  disabled = FALSE,
                                                  style = "warning")
               ),
               conditionalPanel(condition = "input.new_taxa_ask == 'No'",
                                shinyBS::bsButton("but_parse",
                                                  "Get taxonomy info
                                                  from NCBI *",
                                                  disabled = FALSE,
                                                  style = "warning"),
                                helpText(em(
                                  "(*) Taxonomy information for a given taxa
                                  list contains all taxonomy ranks and their
                                  correspoding NCBI IDs"))
                                ),
               
               hr(),
               uiOutput("end_parsing_msg"),
               tableOutput("invalidID.output")
               # strong(h4("PLEASE RELOAD THIS TOOL AFTER ADDING NEW TAXA!!!"),
               #        style = "color:red")
                                ),
             
             conditionalPanel(
               condition = "output.unk_taxa_status == 0",
               strong(h4("Seed (super)taxon:")),
               br(),
               
               strong(h5("Select taxonomy rank:")),
               withSpinner(uiOutput("rank_select"),
                           proxy.height = "50px",
                           type = 7,
                           size = 0.5),
               br(),
               strong(h5("Choose (super)taxon of interest:")),
               withSpinner(uiOutput("select"),
                           proxy.height = "50px",
                           type = 7,
                           size = 0.5),
               br(),
               shinyBS::bsButton("do",
                                 "PLOT",
                                 type = "action",
                                 style = "danger",
                                 size = "large",
                                 disabled = TRUE),
               h5("")
             )
               )
      ),
    
    # MAIN PROFILE TAB ========================================================
    tabPanel(
      "Main profile",
      sidebarLayout(
        sidebarPanel(
          uiOutput("total_gene_number.ui"),
          
          column(4, numericInput("st_index",
                                 "Show from:",
                                 min = 1,
                                 max = 1600,
                                 value = 1,
                                 width = 100),
                 style = "padding:0px;"),
          column(4, numericInput("end_index",
                                 "...to:",
                                 min = 1,
                                 max = 1600,
                                 value = 30,
                                 width = 100),
                 style = "padding:0px;"),
          column(4, uiOutput("highlight_gene_ui")),
          bsPopover("st_index",
                    "",
                    "Set start index for sequence range",
                    "bottom"),
          bsPopover("end_index",
                    "",
                    "Set end index for sequence range",
                    "bottom"),
          
          br(),
          
          uiOutput("highlight_taxon_ui"),
          bsPopover("highlight_taxon_ui",
                    "",
                    "OR double click on heatmap",
                    "right"),
          
          conditionalPanel(
            condition = "input.auto_update == false",
            shinyBS::bsButton("update_btn",
                              "Update plot",
                              style = "warning")
          )
        ),
        
        mainPanel(
          uiOutput("plot.ui"),
          
          conditionalPanel(
            condition = "input.main_x_axis_guide == true |
            input.main_y_axis_guide == true",
            absolutePanel(
              id = "absAxis",
              bottom = 0, left = 0,
              heigh = NULL, width = NULL,
              fixed = TRUE,
              draggable = TRUE,
              style = "opacity: 0.80",
              
              uiOutput("mainAxisRender")
            )
          )
      )
      )
    ),
    
    # CUSTOMIZED PROFILE TAB ==================================================
    tabPanel(
      "Customized profile",
      sidebarLayout(
        sidebarPanel(
          width = 4,
          column(12,
                 style = "padding:0px;",
                 strong("Select sequence(s) of interest:")
          ),
          
          column(12,
                 fluidRow(
                   column(8,
                          style = "padding:0px;",
                          uiOutput("gene_in")
                   ),
                   column(4,
                          fileInput("custom_file", "", width = "100%")
                   )
                 )
          ),
          column(12,
                 style = "padding:0px;",
                 strong("Select (super)taxon/(super)taxa of interest:")
          ),
          column(12,
                 fluidRow(
                   column(8,
                          style = "padding:0px;",
                          uiOutput("taxa_in")
                   ),
                   column(4,
                          h3(""),
                          shinyBS::bsButton("cus_taxa", "Browse...")
                          # select_taxon_rank_ui("cus_taxa")
                   )
                 )
          ),
          
          h5(""),
          uiOutput("plot_custom_btn")
        ),
        mainPanel(
          conditionalPanel(
            condition = "output.same_profile == true",
            h4("Please select subset of genes and/
               or taxa for customized profile!")
            ),
          uiOutput("selected_plot.ui")
          )
      )
      ),
    
    # FUNCTION TAB ============================================================
    navbarMenu(
      "Function",
      
      # Profiles clustering ---------------------------------------------------
      tabPanel(
        "Profiles clustering",
        h4(strong("Profiles clustering")),
        
        wellPanel(
          fluidRow(
            column(3,
                   selectInput("dist_method",
                               label = h5("Distance measure method:"),
                               choices = list("euclidean" = "euclidean",
                                              "maximum" = "maximum",
                                              "manhattan" = "manhattan",
                                              "canberra" = "canberra",
                                              "binary" = "binary"),
                               selected = "euclidean")
            ),
            column(3,
                   selectInput("cluster_method", label = h5("Cluster method:"),
                               choices = list("single" = "single",
                                              "complete" = "complete",
                                              "average (UPGMA)" = "average",
                                              "mcquitty (WPGMA)" = "mcquitty",
                                              "median (WPGMC)" = "median",
                                              "centroid (UPGMC)" = "centroid"),
                               selected = "complete")
            ),
            column(1,
                   numericInput("cluster_plot.width",
                                h5("Width (px)"),
                                min = 200,
                                max = 3200,
                                step = 50,
                                value = 600,
                                width = 100)
            ),
            column(2,
                   numericInput("cluster_plot.height",
                                h5("Height (px)"),
                                min = 200,
                                max = 3200,
                                step = 50,
                                value = 400,
                                width = 100)
            ),
            column(3,
                   checkboxInput("apply_cluster",
                                 em(strong("Apply clustering to profile plot",
                                           style = "color:red")),
                                 value = FALSE),
                   uiOutput("apply_cluster_check.ui"),
                   
                   tags$head(
                     tags$style(HTML(
                       "#download_cluster{background-color:#A9E2F3}"))
                   ),
                   downloadButton("download_cluster", "Download plot")
            )
          )
        ),
        
        column(8,
               uiOutput("cluster.ui")
        ),
        column(4,
               downloadButton("download_cluster_genes", "Download gene list"),
               checkboxInput("add_cluster_cutom_profile",
                             strong(em("Add to Customized profile")),
                             value = FALSE,
                             width = NULL),
               uiOutput("add_cluster_cutom_profile_check.ui"),
               tableOutput("brushed_cluster.table")
        )
      ),
      
      # Distribution analysis -------------------------------------------------
      tabPanel(
        "Distribution analysis",
        h4(strong("Distribution analysis")),
        
        wellPanel(
          fluidRow(
            column(2,
                   selectInput("dataset.distribution", "Select data",
                               choices = c("Main data", "Customized data"),
                               selected = "Main data"),
                   uiOutput("selected.distribution")
            ),
            column(2,
                   uiOutput("var1_dist.ui")
            ),
            column(2,
                   uiOutput("var2_dist.ui")
            ),
            column(2,
                   uiOutput("percent_dist.ui")
            ),
            column(1,
                   numericInput("dist_text_size", "Label size",
                                min = 2,
                                max = 99,
                                step = 1,
                                value = 12,
                                width = 100)
            ),
            column(2,
                   strong("Download"),
                   tags$head(
                     tags$style(HTML(
                       "#plot_download_dist{background-color:#A9E2F3}"))
                   ),
                   downloadButton("plot_download_dist", "Download plot")
            )
          )
        ),
        
        uiOutput("dist_plot.ui")
        # analyze_distribution_ui("dist_plot")
      ),
      
      # Gene age estimation ---------------------------------------------------
      tabPanel(
        "Gene age estimation",
        h4(strong("Gene age estimation")),
        
        wellPanel(
          fluidRow(
            column(2,
                   uiOutput("var1_age.ui")
            ),
            column(2,
                   uiOutput("var2_age.ui")
            ),
            column(2,
                   uiOutput("percent_age.ui")
            ),
            column(2,
                   strong("Appearance"),
                   bsButton("gene_age_prot_config", "Plot config")
            ),
            column(2,
                   strong("Download"),
                   tags$head(
                     tags$style(HTML(
                       "#gene_age_plot_download{background-color:#A9E2F3}"))
                   ),
                   downloadButton("gene_age_plot_download", "Download plot")
            )
          )
        ),
        
        uiOutput("gene_age.ui"),
        conditionalPanel(
          condition = "input.do",
          em(h6("01_Species; 02_Family; 03_Class; 04_Phylum;
                05_Kingdom; 06_Superkingdom; 07_Last universal common ancestor;
                Undef_Genes have been filtered out"))
          ),
        hr(),
        column(4,
               downloadButton("gene_age_table_download", "Download gene list"),
               checkboxInput("add_custom_profile",
                             strong(em("Add to Customized profile")),
                             value = FALSE,
                             width = NULL),
               uiOutput("add_custom_profile_check.ui")
        ),
        tableOutput("gene_age.table")
          ),
      
      # Core gene identification  -----------------------------------------
      tabPanel(
        "Core gene identification",
        h4(strong("Core gene identification")),
        
        wellPanel(
          fluidRow(
            column(2,
                   uiOutput("var1_cons.ui")
            ),
            column(2,
                   uiOutput("var2_cons.ui")
            ),
            column(2,
                   uiOutput("percent_cons.ui")
            ),
            column(6,
                   uiOutput("taxa_list_cons.ui"),
                   shinyBS::bsButton("browse_taxa_cons", "Browse")
            )
          )
        ),
        hr(),
        column(4,
               downloadButton("cons_gene_table_download", "Download gene list"),
               checkboxInput("add_cons_gene_custom_profile",
                             strong(em("Add to Customized profile")),
                             value = FALSE,
                             width = NULL),
               uiOutput("add_cons_gene_custom_profile_check.ui")
        ),
        identify_core_gene_ui("cons_gene")
      ),
      
      # Search for NCBI taxonomy IDs  -----------------------------------------
      search_taxon_id_ui("search_taxon_id"),
      
      # Group Comparison  -----------------------------------------------------
      tabPanel("Group Comparison",
               h4(strong("Group Comparison")),
               wellPanel(
                 fluidRow(
                   column(3,
                          uiOutput("variable_button_gc")
                   ),
                   column(3,
                          uiOutput("list_genes_gc"),
                          popify(fileInput("gc_file", NULL, width = "100%"),
                                 "",
                                 "Select file of sequences")
                   ),
                   column(3,
                          uiOutput("taxa_list_gc"), # Select In-Group
                          checkboxInput("use_common_anchestor",
                                        "Use common anchestor",
                                        value = TRUE,
                                        width = NULL),
                          bsPopover("use_common_anchestor",
                                    "",
                                    "The common anchestor of the selected taxa is used as the in-group",
                                    "top"),
                          shinyBS::bsButton("taxa_gc", "Change rank")
                          
                   ),
                   column(3,
                          uiOutput("significance.ui"),
                          checkboxInput("right_format_features",
                                        "feature format:
                                        ’featuretype_featurename’",
                                        value = TRUE,
                                        width = NULL),
                          actionButton("plot_gc", "Plot"),
                          popify(actionButton("gc_plot_config", "Appearance"),
                                 "",
                                 "Change the appearance of the plots")
                          
                   )
                 )
                 ),
               sidebarPanel(
                 uiOutput("get_significant_genes"),
                 bsPopover("get_significant_genes",
                           "",
                           "Select gene to show the plots",
                           "right"),
                 
                 
                 checkboxInput("add_gc_genes_custom_profile",
                               strong(em("Add to Customized profile")),
                               value = FALSE,
                               width = NULL),
                 uiOutput("add_gc_custom_profile_check"),
                 
                 popify(uiOutput("features_of_interest_gc"),
                        "",
                        "This function is only use full if the features are
                        saved in the right format: featuretype_featurename"),
                 
                 actionButton("gc_downloads",
                              "Download"),
                 
                 width = 3
                 ),
               mainPanel(
                 tags$style(HTML("
                                 #plots_gc {
                                 height:650px;
                                 overflow-y:scroll
                                 }
                                 ")),
                 uiOutput("plots_gc"),
                 width = 9
                 )
                 )
               ),
    
    # DATA DOWNLOAD TAB ===================================================
    navbarMenu(
      "Download filtered data",
      download_filtered_main_ui("filtered_main_download"),
      download_filtered_customized_ui("filtered_customized_download")
    ),
    
    
    # HELP TAB ============================================================
    navbarMenu(
      "Help",
      tabPanel(
        a("Wiki",
          href = {
            "https://github.com/BIONF/PhyloProfile/wiki"
          },
          target = "_blank")
      ),
      tabPanel(
        a("About",
          href = "https://BIONF.github.io/PhyloProfile/",
          target = "_blank")
      )
    )
               ),
  
  # LIST OF POP-UP WINDOWS =================================================
  
  # popup to confirm parsing data from input file ---------------------------
  bsModal("add_taxa_windows", "Add new taxa", "add_taxa", size = "medium",
          helpText(em(
            "Use this form to add taxon that does not exist in NCBI taxonomy
            database (or alternatively you can prepare the data/newTaxa.txt
            file with the following description for each field).")),
          textInput("new_id",
                    "ID (must be a number and greater than 2077091,
                    e.g. 9000001)",
                    9000001,
                    width = 500),
          textInput("new_name",
                    "Name (e.g. Saccharomyces cerevisiae strain ABC)",
                    "",
                    width = 500),
          textInput("new_rank",
                    "Rank (e.g. \"norank\" (for strain),species,
                    order, etc.)",
                    "norank",
                    width = 500),
          textInput("new_parent",
                    "Parent ID (NCBI taxonomy ID of the next higher rank,
                    e.g. 4932 (S.cerevisiae species))",
                    4932,
                    width = 500),
          actionButton("new_add", "Add"),
          actionButton("new_done", "Done")
          ),
  
  # popup to confirm parsing data from input file ---------------------------
  bsModal("parse_confirm", "Get info from input",
          "but_parse",
          size = "small",
          HTML("Fetching Missing Taxonomy Information from NCBI and
               Post-processing...<br><br>"),
          em("This windows will close automatically when eveything
             is done!"),
          br(),
          strong("PLEASE RELOAD THIS TOOL WHEN FINISHED!!!",
                 style = "color:red"),
          # conditionalPanel(
          # condition = 'output.taxonomyParseStatus == 1',
          # strong(h4("Invalid NCBI IDs:")),
          dataTableOutput("invalidIDout")
          # )
          ),
  
  # popup windows for setting plot colors -----------------------------------
  bsModal("color", "Set colors for profile",
          "set_color",
          size = "small",
          colourpicker::colourInput("low_color_var1",
                                    "Low variable 1",
                                    value = "darkorange"),
          colourpicker::colourInput("high_color_var1",
                                    "High variable 1",
                                    value = "steelblue"),
          actionButton("default_color_var1",
                       "Default",
                       style = "padding:4px; font-size:100%"),
          hr(),
          colourpicker::colourInput("low_color_var2",
                                    "Low variable 2",
                                    value = "grey95"),
          colourpicker::colourInput("high_color_var2",
                                    "High variable 2",
                                    value = "khaki"),
          actionButton("default_color_var2",
                       "Default",
                       style = "padding:4px; font-size:100%"),
          hr(),
          colourpicker::colourInput("para_color",
                                    "Color for inparalogs",
                                    value = "#07d000"),
          actionButton("default_color_para",
                       "Default",
                       style = "padding:4px; font-size:100%")
  ),
  
  # popup windows for FASTA configurations ----------------------------------
  bsModal("config",
          "FASTA config",
          "get_config",
          size = "small",
          selectInput("input_type", "Choose location for:",
                      c("Concatenated fasta file", "Fasta folder")
          ),
          hr(),
          uiOutput("default_color_para.ui"),
          conditionalPanel(
            condition = "input.input_type == 'Concatenated fasta file'",
            #            textInput("oneseq.file","Path:",""),
            fileInput("concat_fasta", ""),
            uiOutput("one_seq.exist_check")
          ),
          conditionalPanel(
            condition = "input.input_type == 'Fasta folder'",
            textInput("path", "Main FULL path:", ""),
            selectInput("dir_format", "Directory format:",
                        choices = list("path/speciesID.fa*" = 1,
                                       "path/speciesID/speciesID.fa*" = 2),
                        selected = "Path/speciesID.fasta"),
            selectInput("file_ext", "File extension:",
                        choices = list("fa" = "fa",
                                       "fasta" = "fasta",
                                       "fas" = "fas",
                                       "txt" = "txt"),
                        selected = "fa"),
            selectInput("id_format",
                        "ID format:",
                        choices = list(">speciesID:seqID" = 1,
                                       ">speciesID@seqID" = 2,
                                       ">speciesID|seqID" = 3),
                        selected = 1)
          )
  ),
  
  # popup windows for setting main plot configurations ----------------------
  bsModal("main_plot_config_bs",
          "Plot appearance configuration",
          "main_plot_config",
          size = "small",
          column(6,
                 numericInput("x_size",
                              "X-axis label size (px)",
                              min = 8,
                              max = 99,
                              step = 1,
                              value = 8,
                              width = 100)
          ),
          column(6,
                 
                 numericInput("y_size",
                              "Y-axis label size (px)",
                              min = 8,
                              max = 99,
                              step = 1,
                              value = 8,
                              width = 100)
          ),
          
          column(6,
                 numericInput("legend_size",
                              "Legend label size (px)",
                              min = 8,
                              max = 99,
                              step = 1,
                              value = 8,
                              width = 150)
          ),
          column(6,
                 selectInput("main_legend", label = "Legend position:",
                             choices = list("Right" = "right",
                                            "Left" = "left",
                                            "Top" = "top",
                                            "Bottom" = "bottom",
                                            "Hide" = "none"),
                             selected = "right",
                             width = 150)
          ),
          column(12,
                 HTML("<strong>Angle for x-axis label</strong>:<br>"),
                 sliderInput("x_angle",
                             "",
                             min = 0,
                             max = 90,
                             step = 10,
                             value = 60,
                             width = 250),
                 br()
          ),
          
          column(12,
                 HTML("<strong>Zooming factor (α) for dots on
                      profile</strong>:<br>"),
                 sliderInput("dot_zoom", "",
                             min = -1,
                             max = 3,
                             step = 0.1,
                             value = 0,
                             width = 250),
                 HTML("<em>size = (1+α)*default_size<br>default_size
                      =[0:5]</em>"),
                 uiOutput("dot_size_info"),
                 br()
                 ),
          
          br(),
          hr(),
          shinyBS::bsButton("reset_main_config", "Reset", style = "danger"),
          shinyBS::bsButton("applyMainConfig", "Done", style = "warning")
          ),
  
  # popup windows for setting gene age plot configurations --------------------
  bsModal("gene_age_prot_config_bs", "Plot appearance configuration",
          "gene_age_prot_config",
          size = "small",
          sliderInput("gene_age_width",
                      "Width zoom (*600px)",
                      min = 0,
                      max = 5,
                      step = 0.1,
                      value = 1,
                      width = "100%"),
          sliderInput("gene_age_height", "Height zoom (*150px)",
                      min = 0,
                      max = 5,
                      step = 0.1,
                      value = 1,
                      width = "100%"),
          sliderInput("gene_age_text", "Text size zoom",
                      min = 0,
                      max = 5,
                      step = 0.1,
                      value = 1,
                      width = "100%"),
          br(),
          hr(),
          shinyBS::bsButton("reset_gene_age_prot_config",
                            "Reset",
                            style = "danger")
  ),
  
  # popup windows for setting main plot configurations ------------------------
  bsModal("selected_plot_config_bs", "Plot properties configuration",
          "selected_plot_config",
          size = "small",
          column(6,
                 numericInput("x_size_select", "X-axis label size(px)",
                              min = 8,
                              max = 99,
                              step = 1,
                              value = 8,
                              width = 150)
          ),
          column(6,
                 
                 numericInput("y_size_select", "Y-axis label size (px)",
                              min = 8,
                              max = 99,
                              step = 1,
                              value = 8,
                              width = 100)
          ),
          
          column(6,
                 numericInput("legend_size_select", "Legend label size (px)",
                              min = 8,
                              max = 99,
                              step = 1,
                              value = 8,
                              width = 150)
          ),
          column(6,
                 selectInput("selected_legend", label = "Legend position:",
                             choices = list("Right" = "right",
                                            "Left" = "left",
                                            "Top" = "top",
                                            "Bottom" = "bottom",
                                            "Hide" = "none"),
                             selected = "right",
                             width = 150)
          ),
          column(12,
                 HTML("<strong>Angle for x-axis label</strong>:<br>"),
                 sliderInput("x_angle_select", "",
                             min = 0,
                             max = 90,
                             step = 10,
                             value = 60,
                             width = 250),
                 br()
          ),
          
          column(12,
                 HTML("<strong>Zooming factor (α) for dots on
                      profile</strong>:<br>"),
                 sliderInput("dot_zoom_select", "",
                             min = -1,
                             max = 3,
                             step = 0.1,
                             value = 0,
                             width = 250),
                 HTML("<em>size = (1+α)*default_size<br>default
                      _size=[0:5]</em>"),
                 uiOutput("dot_size_infoSelect"),
                 br()
                 ),
          
          br(),
          hr(),
          shinyBS::bsButton("reset_selected_config", "Reset", style = "danger"),
          shinyBS::bsButton("apply_selected_config", "Done", style = "warning")
          ),
  
  # popup windows for setting group compariosn plot configurations ------------
  bsModal("gc_plot_config_bs", "Plot appearance configuration",
          "gc_plot_config",
          size = "small",
          column(6,
                 numericInput("x_size_gc", "X-axis label size (px)",
                              min = 8,
                              max = 99,
                              step = 1,
                              value = 10,
                              width = 100)
          ),
          column(6,
                 numericInput("y_size_gc", "Y-axis label size (px)",
                              min = 8,
                              max = 99,
                              step = 1,
                              value = 10,
                              width = 100)
          ),
          
          
          column(6,
                 numericInput("legend_size_gc", "Legend label size (px)",
                              min = 8,
                              max = 99,
                              step = 1,
                              value = 10,
                              width = 150)
          ),
          column(6,
                 selectInput("legend_gc", label = "Legend position:",
                             choices = list("Right" = "right",
                                            "Left" = "left",
                                            "Top" = "top",
                                            "Bottom" = "bottom",
                                            "Hide" = "none"),
                             selected = "right",
                             width = 150)
          ),
          column(12,
                 sliderInput("angle_gc", "Angle of the X-axis label",
                             min = 0,
                             max = 180,
                             step = 1,
                             value = 90,
                             width = 250)
          ),
          column(12,
                 checkboxInput("show_p_value",
                               strong("Show P-Values"),
                               value = FALSE,
                               width = 250)
          ),
          column(12,
                 popify(checkboxInput("highlight_significant",
                                      strong("Highlight significant plots"),
                                      value = TRUE,
                                      width = 250),
                        "",
                        "If both variables are selected
                        the significant Plot is colored")
                 ),
          
          br(),
          hr(),
          shinyBS::bsButton("reset_config_gc", "Reset", style = "danger"),
          shinyBS::bsButton("apply_config_gc", "Done", style = "warning")
          ),
  
  
  
  
  # popup windows for select taxa on Customized Profile -----------------------
  bsModal("cus_taxa_bs",
          "Select taxon/taxa of interest",
          "cus_taxa", size = "small",
          select_taxon_rank_ui("select_taxon_rank"),
          checkboxInput("apply_cus_taxa",
                        strong("Apply to customized profile",
                               style = "color:red"),
                        value = FALSE)
  ),
  
  # popup windows for select taxa on Consensus gene finding -------------------
  bsModal("browse_taxa_cons_bs",
          "Select taxon/taxa of interest",
          "browse_taxa_cons",
          size = "small",
          select_taxon_rank_ui("select_taxon_rank_cons"),
          checkboxInput("apply_cons_taxa",
                        strong("Apply", style = "color:red"),
                        value = FALSE)
  ),
  
  # popup windows for detailed plot -------------------------------------------
  bsModal("modal_bs",
          "Detailed plot",
          "detailed_btn",
          size = "large",
          uiOutput("detail_plot.ui"),
          checkboxInput("detailed_remove_na",
                        strong("Hide taxa that have no ortholog (NAs)",
                               style = "color:red"),
                        value = FALSE),
          fluidRow(
            column(
              3,
              numericInput("detailed_height",
                           "Plot height (px)",
                           min = 100,
                           max = 1600,
                           step = 50,
                           value = 100,
                           width = 150)
            ),
            column(
              3,
              numericInput("detailed_text",
                           "Text size (px)",
                           min = 3,
                           max = 30,
                           step = 1,
                           value = 12,
                           width = 150)
            ),
            column(
              3,
              strong("Download"),
              tags$head(
                tags$style(HTML(
                  "#plot_download_dist{background-color:#A9E2F3}"))
              ),
              downloadButton("download_detailed", "Download plot")
            )
          ),
          hr(),
          verbatimTextOutput("detail_click"),
          shinyBS::bsButton("do_domain_plot",
                            "Show domain architecture",
                            disabled = TRUE),
          uiOutput("check_domain_files"),
          br(),
          h4("Sequence:"),
          verbatimTextOutput("fasta")
  ),
  
  # popup windows for domain architecture plot --------------------------------
  bsModal("plot_archi",
          "Domain architecture",
          "do_domain_plot",
          size = "large",
          fluidRow(
            column(2,
                   numericInput("archi_height",
                                "plot height(px)",
                                min = 100,
                                max = 1600,
                                step = 50,
                                value = 400,
                                width = 100)
            ),
            column(2,
                   numericInput("archi_width",
                                "plot width(px)",
                                min = 100,
                                max = 1600,
                                step = 50,
                                value = 800,
                                width = 100)
            ),
            column(2,
                   numericInput("title_archi_size",
                                "Title size(px)",
                                min = 8,
                                max = 99,
                                step = 1,
                                value = 11,
                                width = 150)
            ),
            column(2,
                   numericInput("label_archi_size",
                                "SeqID size(px)",
                                min = 8,
                                max = 99,
                                step = 1,
                                value = 11,
                                width = 150)
            )
          ),
          uiOutput("archi_plot.ui"),
          downloadButton("archi_download", "Download plot")
  ),
  
  # popup window to change the rank in the group comparison function ----------
  bsModal("taxa_gc_bs", "Select taxon/taxa of interest", "taxa_gc",
          size = "small",
          uiOutput("rank_select_gc"),
          uiOutput("taxa_select_gc"),
          checkboxInput("apply_taxa_gc",
                        strong("Apply to group comparison",
                               style = "color:red"),
                        value = FALSE)
  ),
  
  # popup window to handle the downloads to the group comparison function -----
  bsModal("gc_downloadsBS", "Download", "gc_downloads", size = "small",
          h5(strong("Download the significant Genes")),
          downloadButton("download_genes_gc", "Download"),
          h5(""),
          uiOutput("select_plots_to_download "),
          downloadButton("download_plots_gc", "Download")
  ),
  
  # POINT INFO BOX ============================================================
  conditionalPanel(
    condition = "input.tabs=='Main profile' ||
    input.tabs=='Customized profile'",
    # PONIT's INFO BOX --------------------------------------------------------
    absolutePanel(
      bottom = 5, left = 30,
      fixed = TRUE,
      draggable = TRUE,
      h5("Point's info:"),
      verbatimTextOutput("point_info"),
      conditionalPanel(
        condition = "output.point_info_status == 0",
        shinyBS::bsButton("detailed_btn",
                          "Detailed plot",
                          style = "success",
                          disabled = FALSE)
      ),
      style = "opacity: 0.80"
    )
  )
)
)