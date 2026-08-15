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
con_reports <- c("report_loss.Rmd", "report_goal.Rmd")
names(con_reports) <- mode0_choices

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
