# ui for marApp

header <- dashboardHeader(
    title = "marApp",
    titleWidth = sidewidth
)

sidebar <- dashboardSidebar(
    width = sidewidth,
    sidebarMenu(
        menuItem("Upload data", tabName = "data"),
        menuItem("Site frequency spectrum", tabName = "sfs"),
        menuItem("Mutations-area relationship", tabName = "mar"),
        menuItem("Extinction simulation", tabName = "ext"),
        menuItem("Conservation estimator", tabName = "con")
    )
)

body <- dashboardBody(
    tags$head(tags$style(type = "text/css", ".slider-animate-button { font-size: 20pt !important; }")),
    chooseSliderSkin("Flat", color = "green"),
    tabItems(
        # Zero tab for decisions
        tabItem(
            tabName = "con",
            fluidRow(
                box(
                    width = 12,
                    title = "Conservation scenarios", status = "info",
                    radioButtons("mode0",
                        "Estimate genetic diversity loss or build habitat protection goals?",
                        choices = mode0_choices, selected = "estimate loss"
                    ),
                    sliderInput("structure_slider", "Adjust population structure (low to high):",
                        min = 0, max = 1, value = 0.3, step = 0.01
                    ),
                    numericInput("habitat_loss", "Proportion of habitat lost (0 to 1):",
                        value = 0.3, min = 0, max = 1, step = 0.01
                    ),
                    numericInput("gd_target", "Proportion of genetic diversity to protect (0 to 1):",
                        value = 0.9, min = 0, max = 1, step = 0.01
                    ),
                    actionButton("go0", "Estimate", width = 120)
                )
            ),
            conditionalPanel(
                condition = "input.go0",
                fluidRow(
                    box(
                        width = 12,
                        title = "Genetic diversity estimation",
                        uiOutput("reportUI"),
                        downloadButton("downloadReport", "Download report")
                    ),
                )
            )
        ),
        # First tab content
        tabItem(
            tabName = "data",
            fluidRow(
                box(
                    width = 12, collapsible = TRUE,
                    title = "Input validation", status = "info",
                    radioButtons("mode",
                        label = "Run custom dataset or demo?",
                        choices = c("Demo", "Custom"), selected = "Demo"
                    ),
                    uiOutput("uploadNotes"),
                    actionButton("go1", "Load data", width = 120)
                )
            ),
            conditionalPanel(
                condition = "input.go1",
                fluidRow(
                    box(
                        width = 12,
                        title = "Coordinate file preview",
                        verbatimTextOutput("print_mapsdata")
                    ),
                ),
                fluidRow(
                    box(
                        width = 12,
                        title = "Genotype file preview",
                        verbatimTextOutput("print_genodata")
                    )
                ),
                fluidRow(
                    box(
                        width = 12,
                        title = "Genomaps object",
                        leafletOutput("map_genomaps"),
                        # h5("`Sample Raster` layer shows number of samples in grouped cells.\n
                        #    Given Leaflet's automatic reprojection, some cells might be displaced."),
                        # h5("`Sample Points` layer shows inputted sample geo-locations. \n
                        #    Sample IDs are available when selected."),
                    )
                )
            )
        ),
        # Second tab content
        tabItem(
            tabName = "sfs",
            fluidRow(
                box(
                    width = 12, collapsible = TRUE,
                    title = "SFS options", status = "info",
                    radioButtons(
                        inputId = "sfs_plottype",
                        label = "SFS plot style:",
                        choices = c("Bars" = "bar", "Density curves" = "density", "Both" = "both"),
                        selected = "bar", inline = TRUE
                    ),
                    checkboxInput("log_sfs_x", "Log scale x-axis (inspect rare variants)", value = FALSE),
                    checkboxInput("log_sfs_y", "Log scale y-axis (inspect common variants)", value = FALSE),
                    actionButton("go2", "Compute SFS", width = 150),
                )
            ),
            fluidRow(
                box(
                    title = "The site-frequency spectrum",
                    width = 12,
                    withMathJax(includeMarkdown("docs/sfs_explanation.md")),
                )
            ),
            conditionalPanel(
                condition = "input.go2",
                fluidRow(
                    box(
                        width = 12,
                        title = "Site frequency spectrum",
                        plotlyOutput("plot_sfsdf")
                    )
                )
            )
        ),
        # Third tab content
        tabItem(
            tabName = "mar",
            fluidRow(
                box(
                    width = 12, collapsible = TRUE,
                    title = "Mutations-area relationship (MAR) options", status = "info",
                    selectInput(
                        inputId = "scheme",
                        label = "Select the MARsampling scheme:",
                        choices = mar:::.MARsampling_schemes,
                        selected = "random"
                    ),
                    selectInput(
                        inputId = "Mtype",
                        label = "Select the genetic diversity metrics: ",
                        choices = Mchoices,
                        selected = c("M", "thetapi"),
                        multiple = TRUE
                    ),
                    selectInput(
                        inputId = "Atype",
                        label = "Select the area metrics: ",
                        choices = Achoices,
                        selected = "Asq"
                    ),
                    numericInput(
                        inputId = "nrep",
                        label = "Number of replicates:",
                        value = 5,
                        min = 1,
                        max = 20
                    ),
                    numericInput(
                        inputId = "seed",
                        label = "Random seed (optional; leave blank to draw one at random):",
                        value = NA,
                        min = 1,
                        step = 1
                    ),
                    actionButton("go3", "Calculate MAR/GDAR", width = 150)
                )
            ),
            fluidRow(
                box(
                    title = "The mutations (genetic diversity) area relationship",
                    width = 12,
                    withMathJax(includeMarkdown("docs/mar_explanation.md"))
                )
            ),
            conditionalPanel(
                condition = "input.go3",
                fluidRow(
                    box(
                        title = "Summary of MAR/GDAR",
                        width = 12,
                        DT::dataTableOutput("print_marres"),
                        helpText("2.5% and 97.5% columns give the 95% confidence interval of each parameter."),
                        helpText(textOutput("seed_used_mar", inline = TRUE)),
                        downloadButton("download_mardf", "Download data")
                    )
                ),
                fluidRow(
                    box(
                        width = 6,
                        title = "MAR sampling process",
                        uiOutput("slider_mar"),
                        plotOutput("anim_mardf")
                    ),
                    box(
                        width = 6,
                        title = "MAR/GDAR plots",
                        checkboxInput("log_mar", "Plot MAR/GDAR on log scale", value = FALSE),
                        uiOutput("mtype_plot_mar_ui"),
                        plotlyOutput("plot_mardf")
                    )
                )
            )
        ),
        tabItem(
            tabName = "ext",
            fluidRow(
                box(
                    width = 12, collapsible = TRUE,
                    title = "MAR extinction simulations", status = "info",
                    selectInput(
                        inputId = "scheme_ext",
                        label = "Select the MARsampling scheme:",
                        choices = mar:::.MARsampling_schemes,
                        selected = "random"
                    ),
                    selectInput(
                        inputId = "Mtype_ext",
                        label = "Select the genetic diversity metrics: ",
                        choices = Mchoices,
                        selected = c("M", "thetapi"),
                        multiple = TRUE
                    ),
                    selectInput(
                        inputId = "Atype_ext",
                        label = "Select the area metrics: ",
                        choices = Achoices_ext,
                        selected = "A"
                    ),
                    numericInput(
                        inputId = "nrep_ext",
                        label = "Number of replicates:",
                        value = 5,
                        min = 1,
                        max = 20
                    ),
                    numericInput(
                        inputId = "seed_ext",
                        label = "Random seed (optional; leave blank to draw one at random):",
                        value = NA,
                        min = 1,
                        step = 1
                    ),
                    actionButton("go4", "Simulate extinction", width = 150)
                )
            ),
            fluidRow(
                box(
                    title = "Extinction prediction using MAR",
                    width = 12,
                    withMathJax(includeMarkdown("docs/ext_explanation.md"))
                )
            ),
            conditionalPanel(
                condition = "input.go4",
                fluidRow(
                    box(
                        title = "Summary of MAR/GDAR extinction",
                        width = 12,
                        DT::dataTableOutput("print_extres"),
                        helpText("2.5% and 97.5% columns give the 95% confidence interval of each parameter."),
                        helpText(textOutput("seed_used_ext", inline = TRUE)),
                        downloadButton("download_extdf", "Download data")
                    )
                ),
                fluidRow(
                    box(
                        width = 6,
                        title = "MAR extinction process",
                        uiOutput("select_ext"),
                        uiOutput("slider_ext"),
                        plotOutput("anim_extdf")
                    ),
                    box(
                        width = 6,
                        title = "MAR/GDAR extinction plots",
                        uiOutput("mtype_plot_ext_ui"),
                        plotlyOutput("plot_extdf")
                    )
                )
            )
        )
    )
)

tagList(
    dashboardPage(
        header = header,
        sidebar = sidebar,
        body = body,
        skin = "green"
    ),
    tags$footer(
        "© 2024 MOI LAB. Developed by Meixi Lin.",
        style = "width:300px; padding:10px; background-color: #222D32; color: white"
    )
)
