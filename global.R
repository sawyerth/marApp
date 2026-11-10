library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(markdown) # to fix deployment issues

library(leaflet)
library(DT)
library(plotly)

library(ggplot2)
library(dplyr)
library(reshape2)

# library(sads)
library(sf)
library(terra)
library(mar)

sidewidth <- 300
mycrs <- "+proj=longlat +datum=WGS84"

mode0_choices <- c("estimate loss", "build goals")

Mchoices <- mar:::.Mtype
names(Mchoices) <- c("Segregating sites", "Endemic segregating sites", "Watterson's theta", "Nucleotide diversity")
Achoices <- mar:::.Atype
names(Achoices) <- c("Area of cells (km^2)", "Area of squares (degree^2)", "Number of individuals sampled")
Achoices_ext <- Achoices[c("Area of cells (km^2)", "Number of individuals sampled")]

options(shiny.maxRequestSize = 30 * 1024^2) # maximum 30 MB upload

theme_set(theme_bw(base_size = 10))

quiet <- function(x) {
    sink(tempfile())
    on.exit(sink())
    invisible(force(x))
}

# bold text
bt <- function(x) {
    shiny::p(shiny::strong(x))
}

# get name by value
get_name <- function(vec, x) {
    names(vec[vec == x])
}

# run a file-parsing expression and surface a friendly validation message
# (instead of a raw crash) if the file does not match the expected format
safe_parse <- function(expr, filelabel) {
    result <- tryCatch(expr, error = function(e) e)
    shiny::validate(shiny::need(
        !inherits(result, "error"),
        paste0(
            "Failed to parse the ", filelabel, ". Please check that it follows the ",
            "format described above. Details: ",
            if (inherits(result, "error")) conditionMessage(result) else ""
        )
    ))
    return(result)
}

SAMPLEMAP_PAL <- rev(grDevices::terrain.colors(255))


rowcol_to_extent <- function(samplemap, bbox) {
    stopifnot(length(bbox) == 4)
    rr <- if (inherits(samplemap, "Raster")) terra::rast(samplemap) else samplemap
    e <- as.vector(terra::ext(rr)) # xmin, xmax, ymin, ymax
    xres <- terra::xres(rr)
    yres <- terra::yres(rr)
    terra::ext(unname(c(
        e["xmin"] + (bbox[3] - 1) * xres, # xmin from first column
        e["xmin"] + bbox[4] * xres, # xmax from last column
        e["ymax"] - bbox[2] * yres, # ymin from last row (rows count from north)
        e["ymax"] - (bbox[1] - 1) * yres # ymax from first row
    )))
}

#######################################################################
# Conservation estimator
# The report is built from regular shiny tags so nothing is written to disk
# while the app runs; a standalone HTML copy is only produced when the user
# clicks "Download report".

# Solve M = (1 - A)^z for the unknown given the chosen mode.
# A: proportion of habitat lost; M: proportion of genetic diversity retained.
con_estimate <- function(mode, z, aloss, gtarg) {
    if (mode == "estimate loss") {
        Atar <- aloss
        Mtar <- (1 - Atar)^z
    } else {
        Mtar <- gtarg
        Atar <- 1 - Mtar^(1 / z)
    }
    list(mode = mode, z = z, aloss = aloss, gtarg = gtarg, Atar = Atar, Mtar = Mtar)
}

con_plot <- function(est) {
    preddf <- data.frame(A = seq(0, 1, by = 0.005))
    preddf$M <- (1 - preddf$A)^est$z
    forplot <- data.frame(A = est$Atar, M = est$Mtar)
    ggplot(data = forplot, aes(x = A, y = M)) +
        geom_line(data = preddf, color = "lightgreen") +
        geom_point(size = 3, color = "darkblue") +
        geom_segment(aes(x = A, xend = A, y = 0, yend = M), linetype = "dashed", color = "lightgray") +
        geom_segment(aes(x = 0, xend = A, y = M, yend = M), linetype = "dashed", color = "lightgray") +
        annotate("text", x = est$Atar, y = 0, label = scales::percent(est$Atar), vjust = 1.5, size = 3.5) +
        annotate("text", x = 0, y = est$Mtar, label = scales::percent(est$Mtar), hjust = -0.2, size = 3.5) +
        scale_x_continuous(labels = scales::percent) +
        scale_y_continuous(labels = scales::percent) +
        labs(x = "% of Area lost", y = "% of Segregating sites remained")
}

# Report body as shiny tags. `plot` is the tag that shows the figure: a
# plotOutput() in the app, an <img> in the downloaded file.
con_report_tags <- function(est, plot) {
    pct <- scales::percent
    b <- function(x) strong(x, .noWS = "outside")
    isloss <- est$mode == "estimate loss"
    tagList(
        h3(if (isloss) "Genetic diversity loss based on MAR" else "Habitat conservation goal based on MAR"),
        h4("Overview"),
        p(
            "This report summarizes predictions based on the ", b("mutations-area relationship (MAR)"),
            ", which describes how genetic diversity scales with geographic area. The MAR follows a power law:"
        ),
        p("$$M = cA^z$$"),
        p("where:"),
        tags$ul(
            tags$li("\\(M\\) is the genetic diversity measured as the number of segregating mutations,"),
            tags$li("\\(A\\) is the area sampled for genetic data,"),
            tags$li("\\(z\\) is the scaling exponent, influenced by population structure, demographic history, and species biological traits.")
        ),
        p(if (isloss) {
            "This can be rearranged to estimate the genetic diversity retained given habitat loss:"
        } else {
            "This can be rearranged to estimate habitat conservation goals given genetic diversity goals:"
        }),
        p(if (isloss) {
            "$$\\text{Proportion of genetic diversity retained} = (1 - \\text{Proportion of habitat lost})^z$$"
        } else {
            "$$\\text{Proportion of genetic diversity to preserve} = (1 - \\text{Maximum proportion of habitat lost})^z$$"
        }),
        h4("User input"),
        tags$ul(
            tags$li(b("Mode: "), est$mode),
            tags$li(b("Population structure (z): "), est$z),
            if (isloss) {
                tags$li(b("Percentage of habitat already lost: "), pct(est$aloss))
            } else {
                tags$li(b("Percentage of genetic diversity to preserve: "), pct(est$gtarg))
            }
        ),
        h4(if (isloss) "Estimate genetic diversity loss" else "Build habitat protection goals"),
        plot,
        h4("Interpretation"),
        if (isloss) {
            p(
                "Given the current loss of ", b(pct(est$Atar)), " of habitat, and a MAR scaling factor z of ",
                b(est$z), ", approximately ", b(pct(est$Mtar)), " of genetic diversity is still intact, and ",
                b(pct(1 - est$Mtar)), " of genetic diversity is already lost."
            )
        } else {
            p(
                "To preserve at least ", b(pct(est$Mtar)), " of genetic diversity intact, at most ",
                b(pct(est$Atar)), " of habitat can be lost, and at least ", b(pct(1 - est$Atar)),
                " of habitat needs to be protected, given a MAR scaling factor z of ", b(est$z), "."
            )
        },
        p(
            "Genetic diversity is measured as allelic richness in this report. ",
            "Higher population structure (higher z) results in faster genetic diversity loss with habitat loss. ",
            "To retain more genetic diversity, more area must be protected."
        )
    )
}

# Write a standalone HTML report to `file` (the temp path from downloadHandler).
con_report_save <- function(est, file) {
    png <- tempfile(fileext = ".png")
    on.exit(unlink(png))
    ggsave(png, con_plot(est), width = 7, height = 5, dpi = 120)
    img <- tags$img(src = base64enc::dataURI(file = png, mime = "image/png"), style = "max-width:100%;")
    page <- tagList(
        tags$head(
            tags$title("MAR conservation report"),
            tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.9/MathJax.js?config=TeX-AMS-MML_HTMLorMML")
        ),
        div(
            style = "font-family: sans-serif; max-width: 800px; margin: 2em auto; padding: 0 1em;",
            con_report_tags(est, img)
        )
    )
    htmltools::save_html(page, file)
}
