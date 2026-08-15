# `marApp`: Mutations-Area Relationship shinyApp

## Overview

`marApp` is an web-based application that enables users to perform MAR analysis without requirements for coding experiences. Users can upload their genotype and coordinate data, select the desired analysis steps, and visualize the results.

## Installation

The hosted version of `marApp` is available on [shinyapps.io](https://www.shinyapps.io), but it can also be run locally. Running locally avoids uploading data to a third-party server and avoids shinyapps.io's usage limits, at the cost of a heavier local install (the app pulls in `shinydashboard`, `leaflet`, `plotly`, and `sf`/`terra`, among others) — this is why `marApp` is distributed as a standalone app rather than bundled into the `mar` package itself.

1. Install R (>= 4.1) and the `mar` package and its dependencies:
   ```r
   install.packages(c("shiny", "shinydashboard", "shinyWidgets", "markdown",
                       "leaflet", "DT", "plotly", "ggplot2", "dplyr", "reshape2",
                       "sf", "raster"))
   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   remotes::install_github("meixilin/mar")
   ```
2. Clone or download this repository, then launch the app from within the `marApp` directory:
   ```r
   shiny::runApp(".")
   ```

## Troubleshooting

### 1. No output is generated in the `Mutations-area relationship` tab.

1. Reload the web page.
2. Navigate to the `Upload data` tab.
3. Click the `Load data` button to load the example data. Or upload your own data and click the `Load data` button.
4. Navigate to the `Mutations-area relationship` tab.
5. Click the `Calculate MAR/GDAR` button.
