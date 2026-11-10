function(input, output, session) {
    #######################################################################
    # Conservation estimator. The report is shown with regular shiny tags (no
    # file is written); a standalone HTML copy is generated only when the user
    # clicks "Download report", straight into the download's temp path.
    con_est <- eventReactive(input$go0, {
        req(input$mode0)
        isloss <- input$mode0 == "estimate loss"
        prop <- if (isloss) input$habitat_loss else input$gd_target
        validate(
            need(is.numeric(prop) && !is.na(prop) && prop >= 0 && prop <= 1, "Please enter a proportion between 0 and 1."),
            need(isloss || input$structure_slider > 0, "Population structure (z) must be above 0 to build goals.")
        )
        con_estimate(input$mode0, input$structure_slider, input$habitat_loss, input$gd_target)
    })

    output$con_plot <- renderPlot(con_plot(con_est()))

    output$reportUI <- renderUI({
        withMathJax(con_report_tags(con_est(), plotOutput("con_plot")))
    })

    output$downloadReport <- downloadHandler(
        filename = function() {
            paste0("Report_", Sys.Date(), ".html")
        },
        content = function(file) {
            con_report_save(con_est(), file)
        }
    )

    #######################################################################
    # renderUI for Upload data input selection
    output$uploadNotes <- renderUI({
        req(input$mode)
        if (input$mode == "Custom") {
            tagList(
                includeMarkdown("docs/upload_req.md"),
                fileInput(
                    inputId = "in_coords",
                    label = "Coordinate file",
                    accept = c(".txt", ".txt.gz", ".csv", ".csv.gz", ".tsv", ".tsv.gz")
                ),
                fileInput(
                    inputId = "in_genomes",
                    label = "Genotype file",
                    accept = c(".txt", ".txt.gz", ".tsv", ".tsv.gz", ".vcf", ".vcf.gz")
                )
            )
        }
    })

    #######################################################################
    # Load and build genomaps object. Use gm1001g as default.
    # Custom uploads are validated: parsing errors surface an informative
    # message instead of crashing the app.
    mapsdata <- reactive({
        req(input$go1)
        if (input$mode == "Custom") {
            req(input$in_coords)
            mypath <- input$in_coords$datapath
            obj <- safe_parse(mar::lonlat_parser(mypath, mapres = NULL, mapcrs = mycrs), "coordinate file")
        } else {
            obj <- mar::gm1001g$maps
        }
        return(obj)
    })

    genodata <- reactive({
        req(input$go1)
        if (input$mode == "Custom") {
            req(input$in_genomes)
            mypath <- input$in_genomes$datapath
            if (grepl(".vcf", input$in_genomes$name)) {
                obj <- safe_parse(quiet(mar::vcf_parser(mypath)), "genotype file")
            } else {
                obj <- safe_parse(mar::text_parser(mypath), "genotype file")
            }
        } else {
            obj <- mar::gm1001g$geno
        }
        return(obj)
    })

    # Build genomaps
    gm <- reactive({
        req(mapsdata(), genodata())
        safe_parse(
            mar::genomaps(geno = genodata(), maps = mapsdata()),
            "uploaded data (sample IDs in the coordinate and genotype files must match and be in the same order)"
        )
    })

    #######################################################################
    # Print data
    output$print_mapsdata <- renderPrint({
        mar:::print.marmaps(mapsdata())
    })
    output$print_genodata <- renderPrint({
        mar:::print.margeno(genodata())
    })

    #######################################################################
    # Plot leaflet
    # TODO: leaflet automatically reprojects, so the raster is not exactly what used
    output$map_genomaps <- leaflet::renderLeaflet({
        # convert gm() to a data frame
        mapdf <- cbind(
            data.frame(ID = gm()$maps$sample.id),
            as.data.frame(gm()$maps$lonlat)
        )
        leaflet() %>%
            addTiles() %>%
            addRasterImage(
                x = mar:::.get_samplemap(gm()$maps),
                group = "Sample Raster"
            ) %>%
            addCircleMarkers(
                data = mapdf,
                clusterOptions = markerClusterOptions(),
                popup = ~ID, label = ~ID, lng = ~LON, lat = ~LAT,
                group = "Sample Points"
            ) %>%
            addLayersControl(
                overlayGroups = c("Sample Raster", "Sample Points"),
                options = layersControlOptions(collapsed = FALSE)
            )
    })

    #########################################################################
    sfslist <- eventReactive(input$go2, {
        req(gm())
        genosfs <- mar::sfs(gm = gm(), folded = FALSE)
        neutralsfs <- mar::expsfs(gm = gm(), folded = FALSE)
        allsfs <- list(observed = genosfs, neutral = neutralsfs)
        return(allsfs)
    })

    # plot sfs (bars, density curves, or both)
    output$plot_sfsdf <- renderPlotly({
        forplot <- lapply(names(sfslist()), function(model) {
            sfsvec <- sfslist()[[model]]
            data.frame(AC = as.numeric(names(sfsvec)), model = model, value = as.numeric(sfsvec))
        }) %>% dplyr::bind_rows()
        # The neutral expectation theta/k is an expected number of sites, so it
        # falls below 1 at high allele counts (and is 0 at k = 0). On a log y-axis
        # those bins would plot as negative values (or -Inf), implying fewer than
        # one site; observed counts are integers, so nothing below 1 can be seen
        # in the data either. Omit them, as the plot.sfs() method in the mar
        # package effectively does.
        if (isTRUE(input$log_sfs_y)) {
            forplot <- forplot[forplot$value >= 1, ]
        }
        pp <- ggplot(data = forplot, mapping = aes(x = AC, y = value, fill = model, color = model, group = model))
        if (input$sfs_plottype %in% c("bar", "both")) {
            pp <- pp + geom_col(position = "dodge", alpha = if (input$sfs_plottype == "both") 0.5 else 1)
        }
        if (input$sfs_plottype %in% c("density", "both")) {
            pp <- pp + geom_line(linewidth = 1)
        }
        if (isTRUE(input$log_sfs_x)) {
            pp <- pp + scale_x_log10()
        }
        if (isTRUE(input$log_sfs_y)) {
            pp <- pp + scale_y_log10()
        }
        pp <- pp + labs(x = "Allele Count", y = "Number of Alleles", fill = "Model", color = "Model")
        ggplotly(pp)
    })

    #######################################################################
    # Calculate MAR only when "Calculate MAR/GDAR" is clicked.
    seed_mar <- eventReactive(input$go3, resolve_seed(input$seed))

    mardf <- eventReactive(input$go3, {
        req(gm())
        withProgress(message = "Calculating MAR ...", {
            MARsampling(
                gm = gm(),
                scheme = input$scheme,
                nrep = input$nrep,
                xfrac = 0.01,
                quorum = TRUE,
                myseed = seed_mar()
            )
        })
    })

    output$seed_used_mar <- renderText({
        req(mardf())
        seed_note(seed_mar())
    })

    # snapshot the settings used for this calculation, so downstream sliders
    # and animations stay consistent even if the user changes the inputs
    nrep_mar <- eventReactive(input$go3, input$nrep)
    mtype_mar <- eventReactive(input$go3, input$Mtype)
    atype_mar <- eventReactive(input$go3, input$Atype)

    # Build MAR
    marres <- eventReactive(input$go3, {
        mar_fit_table(mardf(), mtype_mar(), atype_mar())
    })

    # Print MAR results with 95% CIs (row names and column labels match the plotting panel)
    output$print_marres <- DT::renderDataTable({
        mar_fit_datatable(marres())
    })

    # Option to download mardf()
    output$download_mardf <- downloadHandler(
        filename = function() {
            paste0("mardf_seed", seed_mar(), "_", Sys.Date(), ".csv")
        },
        content = function(file) {
            write.csv(mardf(), file)
        }
    )

    # Restrict the plotting metric to whichever metrics were actually calculated
    output$mtype_plot_mar_ui <- renderUI({
        req(marres())
        choices <- Mchoices[Mchoices %in% mtype_mar()]
        selectInput("Mtype_plot_mar", "Select the genetic diversity metric to plot:",
            choices = choices, selected = choices[1]
        )
    })

    # Plot MAR results
    output$plot_mardf <- renderPlotly({
        req(input$Mtype_plot_mar)
        forplot <- mardf()[, c(atype_mar(), input$Mtype_plot_mar)]
        forplot <- forplot[(forplot[, 2] > 0 & !is.na(forplot[, 2])), ]
        c <- marres()[rownames(marres()) == input$Mtype_plot_mar, "c"]
        z <- marres()[rownames(marres()) == input$Mtype_plot_mar, "z"]
        # make predictions table (since stat_function does not work)
        preddf <- data.frame(x = sort(unique(forplot[, 1]))) %>%
            dplyr::mutate(y = c * x^z)
        colnames(preddf) <- colnames(forplot)
        pp <- ggplot(
            data = forplot,
            mapping = aes(x = .data[[atype_mar()]], y = .data[[input$Mtype_plot_mar]])
        ) +
            geom_point(size = 1, color = "darkgreen") +
            geom_line(data = preddf, color = "darkgray") +
            labs(x = get_name(Achoices, atype_mar()), y = get_name(Mchoices, input$Mtype_plot_mar))

        if (input$log_mar) {
            pp <- pp +
                scale_x_log10() +
                scale_y_log10()
        }
        ggplotly(pp)
    })

    # generate an output slider
    output$slider_mar <- renderUI({
        req(mardf())
        sliderInput(
            inputId = "a_mar",
            label = "Sampling box size:",
            min = 1,
            max = nrow(mardf()) / nrep_mar(),
            step = 1,
            value = 1,
            animate = animationOptions(interval = 1000, loop = FALSE)
        )
    })

    output$anim_mardf <- renderPlot({
        req(input$a_mar, mardf())
        # get the extdf
        bboxlist <- lapply(strsplit(mardf()$extent, ";"), as.integer)
        idx <- (input$a_mar - 1) * nrep_mar() + 1
        sm <- mar:::.get_samplemap(gm()$maps)
        par(mar = c(5.1, 4.1, 4.1, 4.1))
        terra::plot(terra::ext(sm), xlab = "lon", ylab = "lat")
        terra::plot(sm, add = T, col = SAMPLEMAP_PAL, legend.mar = 3, legend.args = list(text = "# of genomes", side = 2))
        for (ii in idx:(idx + nrep_mar() - 1)) {
            terra::plot(rowcol_to_extent(sm, bboxlist[[ii]]),
                col = NA, border = "black", lwd = 2, add = T
            )
        }
    })

    #######################################################################
    # Run extinction simulations only when "Simulate extinction" is clicked.
    seed_ext <- eventReactive(input$go4, resolve_seed(input$seed_ext))

    extdf <- eventReactive(input$go4, {
        req(gm())
        withProgress(message = "Simulating extinction ...", {
            MARextinction(
                gm = gm(),
                scheme = input$scheme_ext,
                nrep = input$nrep_ext,
                xfrac = 0.01,
                myseed = seed_ext()
            )
        })
    })

    output$seed_used_ext <- renderText({
        req(extdf())
        seed_note(seed_ext())
    })

    nrep_ext <- eventReactive(input$go4, input$nrep_ext)
    mtype_ext <- eventReactive(input$go4, input$Mtype_ext)
    atype_ext <- eventReactive(input$go4, input$Atype_ext)

    # Build MAR based on EXT data
    extres <- eventReactive(input$go4, {
        mar_fit_table(extdf(), mtype_ext(), atype_ext())
    })

    # Print MAR results based on EXT data, with 95% CIs
    output$print_extres <- DT::renderDataTable({
        mar_fit_datatable(extres())
    })

    # Option to download extdf()
    output$download_extdf <- downloadHandler(
        filename = function() {
            paste0("extdf_seed", seed_ext(), "_", Sys.Date(), ".csv")
        },
        content = function(file) {
            write.csv(extdf(), file)
        }
    )


    output$mtype_plot_ext_ui <- renderUI({
        req(extres())
        choices <- Mchoices[Mchoices %in% mtype_ext()]
        selectInput("Mtype_plot_ext", "Select the genetic diversity metric to plot:",
            choices = choices, selected = choices[1]
        )
    })

    # Plot EXT results (% lost vs % lost)
    output$plot_extdf <- renderPlotly({
        req(input$Mtype_plot_ext)
        forplot <- extdf()[, c(atype_ext(), input$Mtype_plot_ext, "repid")] %>% na.omit()
        # generate percentages data
        forplot[, 1] <- 1 - forplot[, 1] / max(forplot[, 1])
        forplot[, 2] <- forplot[, 2] / max(forplot[, 2])
        # get c and z again from EXT output
        z <- extres()[rownames(extres()) == input$Mtype_plot_ext, "z"]
        # make predictions table (since stat_function does not work)
        preddf <- data.frame(x = sort(unique(forplot[, 1]))) %>%
            dplyr::mutate(y = (1 - x)^z)
        colnames(preddf) <- colnames(forplot)[1:2]
        pp <- ggplot(
            data = forplot,
            mapping = aes(x = .data[[atype_ext()]], y = .data[[input$Mtype_plot_ext]], color = .data[["repid"]])
        ) +
            geom_point(size = 1) +
            geom_line(data = preddf, color = "darkgray") +
            scale_color_gradient(low = "lightgreen", high = "darkgreen") +
            scale_x_continuous(labels = scales::percent) +
            scale_y_continuous(labels = scales::percent) +
            labs(
                x = paste0("% of ", get_name(Achoices, atype_ext()), " lost"),
                y = paste0("% of ", get_name(Mchoices, input$Mtype_plot_ext), " remained")
            ) +
            theme(legend.position = "none")
        ggplotly(pp)
    })

    # generate an output slider for extinction visualizations
    output$select_ext <- renderUI({
        req(extdf())
        numericInput(
            inputId = "repid_ext", label = "Select which simulations to animate:",
            value = 1, min = 1, max = nrep_ext(), step = 1
        )
    })
    output$slider_ext <- renderUI({
        req(extdf())
        sliderInput(
            inputId = "a_ext",
            label = "Extinction step:",
            min = 1,
            max = nrow(extdf()) / nrep_ext(),
            step = 1,
            value = 1,
            animate = animationOptions(interval = 1000, loop = FALSE)
        )
    })

    output$anim_extdf <- renderPlot({
        req(input$a_ext, input$repid_ext, extdf())
        # get the given extdf()
        extdf0 <- extdf()[extdf()$repid == input$repid_ext, ]
        extl <- lapply(strsplit(extdf0$extl, ";"), as.integer)
        sm <- mar:::.get_samplemap(gm()$maps)
        rr <- sm
        terra::values(rr) <- NA
        par(mar = c(5.1, 4.1, 4.1, 4.1))
        terra::plot(terra::ext(sm), xlab = "lon", ylab = "lat")
        terra::plot(sm, add = T, col = SAMPLEMAP_PAL, legend.mar = 3, legend.args = list(text = "# of genomes", side = 2))
        rr[setdiff(gm()$maps$cellid, extl[[input$a_ext]])] <- 1
        terra::plot(rr, add = T, col = "black", legend = FALSE)
    })
}
