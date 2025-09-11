# app.R
# CyTOF Explorer — upload-enabled, polygon-gating Shiny app

# ---- Packages ----
suppressPackageStartupMessages({
  suppressWarnings({
    library(shiny)
    library(shinyWidgets)
    library(colourpicker)
    library(plotly)
    library(ggplot2)
    library(ggrastr)
    library(ggrepel)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(stringr)
    library(data.table)
    library(ComplexHeatmap)
    library(circlize)
    library(viridis)
    library(scales)
    library(sf)         # point-in-polygon
    library(jsonlite)   # gate import/export
  })
})

options(shiny.maxRequestSize = 100000 * 1024^2)

# ---- Helpers & validation ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

pct_clip <- function(x, p = c(0.01, 0.99)) {
  q <- quantile(x, probs = p, na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

spearman_test <- function(df, freq_col = "freq", cont_var) {
  ok <- complete.cases(df[[freq_col]], df[[cont_var]])
  if (!any(ok)) return(data.frame(rho = NA_real_, p = NA_real_, n = 0))
  ct <- suppressWarnings(cor.test(df[[freq_col]][ok], df[[cont_var]][ok], method = "spearman"))
  data.frame(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

freq_by <- function(cell_ids, meta_cell, group_var, unit_var) {
  dt <- as.data.table(meta_cell)
  dt[, .idx := .I]
  dt[, .in := .idx %in% cell_ids]
  
  tot <- dt[, .(total = .N), by = unit_var]
  inc <- dt[.in == TRUE, .(count = .N), by = unit_var]
  res <- merge(tot, inc, by = unit_var, all.x = TRUE)
  res$count[is.na(res$count)] <- 0
  
  gv <- unique(dt[, c(unit_var, group_var), with = FALSE])
  res <- merge(res, gv, by = unit_var, all.x = TRUE)
  res$freq <- res$count / res$total
  as.data.frame(res)
}

# Parse SVG path "M x,y L x,y ... Z" into numeric coords
parse_svg_path <- function(path) {
  path <- gsub("[MmZz]", "", path)
  coords <- unlist(strsplit(path, "[Ll]"))
  coords <- trimws(coords)
  coords <- coords[nzchar(coords)]
  xy <- do.call(rbind, strsplit(coords, ","))
  data.frame(x = as.numeric(xy[,1]), y = as.numeric(xy[,2]))
}

# Validate minimal structure
validateInput <- function(obj, id_col = NULL) {
  required <- c("data", "source", "metadata", "leiden")
  miss <- setdiff(required, names(obj))
  if (length(miss)) stop("Missing required elements: ", paste(miss, collapse = ", "))
  
  if (!is.matrix(obj$data)) stop("data must be a matrix (cells × features).")
  n <- nrow(obj$data)
  
  if (length(obj$source) != n) stop("source length must equal nrow(data).")
  if (!is.data.frame(obj$metadata)) stop("metadata must be a data.frame.")
  if (!is.list(obj$leiden) || is.null(obj$leiden$clusters)) stop("leiden must be a list with 'clusters'.")
  if (length(obj$leiden$clusters) != n) stop("leiden$clusters length must equal nrow(data).")
  
  # Optional elements: if present, check shape
  if ("umap" %in% names(obj)) {
    if (!is.data.frame(obj$umap$coordinates) || nrow(obj$umap$coordinates) != n)
      stop("umap$coordinates must align with cells.")
  }
  if ("tsne" %in% names(obj)) {
    if (!is.data.frame(obj$tsne$coordinates) || nrow(obj$tsne$coordinates) != n)
      stop("tsne$coordinates must align with cells.")
  }
  if ("run_date" %in% names(obj)) {
    if (length(obj$run_date) != n) stop("run_date must be length nrow(data).")
  }
  
  # id_col presence check (if provided)
  if (!is.null(id_col) && !(id_col %in% colnames(obj$metadata))) {
    stop("Selected ID column '", id_col, "' not found in metadata.")
  }
  invisible(TRUE)
}

# Presence checks
hasUMAP <- function(obj) "umap" %in% names(obj) && !is.null(obj$umap$coordinates)
hasTSNE <- function(obj) "tsne" %in% names(obj) && !is.null(obj$tsne$coordinates)
hasHeatmap <- function(obj) "leiden_heatmap" %in% names(obj) && !is.null(obj$leiden_heatmap$heatmap_tile_data)
hasRunDate <- function(obj) "run_date" %in% names(obj) && length(obj$run_date) == nrow(obj$data)
hasClusterMapping <- function(obj) "cluster_mapping" %in% names(obj)

# Auto-detect metadata ID column that matches source
guess_id_col <- function(metadata, source_vec) {
  # Prioritized candidates by name
  candidates <- intersect(tolower(colnames(metadata)),
                          c("patientid","patient_id","patient","source","subject","sample","id"))
  if (length(candidates)) {
    hit <- candidates[1]
    return(colnames(metadata)[tolower(colnames(metadata)) == hit][1])
  }
  # Fallback: choose column with max overlap
  overlaps <- sapply(metadata, function(col) {
    if (!is.atomic(col)) return(0)
    length(intersect(as.character(col), as.character(source_vec)))
  })
  if (all(overlaps == 0)) return(NULL)
  colnames(metadata)[which.max(overlaps)]
}


# ---- Gate store ----
GateStore <- function() {
  rv <- reactiveValues(.list = list())
  list(
    add = function(g) { rv$.list[[g$name]] <- g },
    list = reactive(rv$.list),
    remove = function(name) { rv$.list[[name]] <- NULL },
    clear = function() { rv$.list <- list() }
  )
}

new_gate <- function(name, cells, embedding, polygon = NULL, color = "#E45756") {
  list(
    name = name,
    cells = sort(unique(as.integer(cells))),
    embedding = embedding, # "UMAP" or "tSNE"
    polygon = polygon,     # list(x=..., y=...) if drawn
    color = color,
    created = as.character(Sys.time())
  )
}


# ---- Embedding module UI ----
EmbeddingUI <- function(id, title = "UMAP") {
  ns <- NS(id)
  tagList(
    h3(title),
    fluidRow(
      column(
        3,
        pickerInput(ns("color_by"), "Color by", choices = NULL),
        checkboxInput(ns("show_labels"), "Show cluster labels", value = FALSE),
        pickerInput(ns("split_by"), "Split by (optional) — metadata", choices = NULL,
                    options = list(`none-selected-text`="None")),
        uiOutput(ns("split_levels_ui")),
        actionButton(ns("plot_facets"), "Plot facets"), 
        sliderInput(ns("max_facets"), "Max facets", min = 2, max = 3, value = 2, step = 1),
        hr(),
        # --- Gating controls (commented out for now) ---
        radioButtons(ns("gate_mode"), "Gating mode",
                     choices = c("Lasso select", "Draw polygon"),
                     selected = "Draw polygon"),
        # textInput(ns("gate_name"), "Gate name", value = ""),
        colourInput(ns("gate_color"), "Gate color", value = "#E45756"),
        # actionButton(ns("save_gate"), "Save gate"),
        # actionButton(ns("clear_selection"), "Clear current selection"),
        # hr(),
        # pickerInput(ns("overlay_gate"), "Overlay gate(s)", choices = NULL, multiple = TRUE)
      ),
      column(
        9,
        plotlyOutput(ns("embed_plot"), height = "650px")
      )
    ),
    # hr(),
    # --- Gate phenotype and abundance section (commented out for now) ---
    # h4("Gate phenotype and abundance"),
    # fluidRow(
    #   column(
    #     4,
    #     pickerInput(ns("phenotype_gate"), "Gate for phenotype", choices = NULL),
    #     pickerInput(ns("phenotype_fun"), "Summary", choices = c("median", "p90"), selected = "median")
    #   ),
    #   column(
    #     8,
    #     plotOutput(ns("inout_heatmap"), height = "400px")
    #   )
    # ),
    # hr(),
    # h4("Abundance testing"),
    # fluidRow(
    #   column(
    #     4,
    #     pickerInput(ns("test_entity"), "Entity",
    #                 choices = c("Selected gate(s)", "Clusters", "Celltypes")),
    #     pickerInput(ns("test_gate"), "Gate(s)", choices = NULL, multiple = TRUE),
    #     pickerInput(ns("group_var"), "Grouping factor (metadata)", choices = NULL),
    #     pickerInput(ns("cont_var"), "Continuous metadata", choices = NULL),
    #     radioButtons(ns("test_type"), "Test",
    #                  choices = c("Wilcoxon (2-group)","Kruskal–Wallis (multi-group)","Spearman (continuous)")),
    #     pickerInput(ns("unit_var"), "Aggregation unit", choices = NULL),
    #     checkboxInput(ns("apply_bh"), "BH adjust across entities", TRUE),
    #     actionButton(ns("run_test"), "Run tests")
    #   ),
    #   column(
    #     8,
    #     plotOutput(ns("abund_plot"), height = "300px"),
    #     tableOutput(ns("test_table"))
    #   )
    # )
  )
}

# ---- Embedding module server ----
EmbeddingServer <- function(id, embedding_name, coords, expr, meta_cell, clusters, cluster_map,
                            gate_store, active_tab, rv) {
  moduleServer(id, function(input, output, session) {
    message(sprintf("EmbeddingServer %s started", embedding_name))
    message(sprintf("coords NULL? %s | expr NULL? %s | meta_cell NULL? %s",
                    is.null(coords), is.null(expr), is.null(meta_cell)))
    
    # One-time picker initialization in a reactive context
    initialized <- FALSE
    
    observeEvent(list(expr(), meta_cell(), clusters()), ignoreInit = FALSE, {
      if (initialized) return()
      expr_val     <- expr()
      meta_val     <- meta_cell()
      clusters_val <- clusters()
      
      req(!is.null(expr_val), !is.null(meta_val))
      
      # Add leiden_cluster column as factor
      if (!"leiden_cluster" %in% colnames(meta_val) && !is.null(clusters_val$assignments)) {
        meta_val$leiden_cluster <- factor(clusters_val$assignments)
      }
      
      numeric_markers <- colnames(expr_val)
      meta_cols <- setdiff(colnames(meta_val), c(".cell"))
      cont_choices <- meta_cols[sapply(meta_val[meta_cols], is.numeric)]
      unit_candidates <- intersect(
        c("PatientID", "patient_ID", "patient", "source", "RunDate", "run_date"),
        colnames(meta_val)
      )
      unit_default <- if (length(unit_candidates)) unit_candidates[1] else meta_cols[1]
      
      updatePickerInput(session, "color_by",
                        choices = c(numeric_markers, meta_cols),
                        selected = if (length(numeric_markers)) numeric_markers[1] else meta_cols[1])
      updatePickerInput(session, "split_by", choices = c("", meta_cols), selected = "")
      updatePickerInput(session, "group_var",
                        choices = c("", meta_cols),
                        selected = "")
      updatePickerInput(session, "cont_var",
                        choices = c("", cont_choices),
                        selected = "")
      updatePickerInput(session, "unit_var", choices = meta_cols, selected = unit_default)
      
      message(sprintf("Picker inputs initialized for %s (reactive observer)", embedding_name))
      initialized <<- TRUE
    })
    
    observeEvent(input$split_by, {
      split_var <- input$split_by
      meta_val <- meta_cell()
      if (!is.null(split_var) && nzchar(split_var) && split_var %in% colnames(meta_val)) {
        levs <- sort(unique(as.character(meta_val[[split_var]])))
        updatePickerInput(session, "split_levels", choices = levs, selected = levs)
      }
    })
    
    output$split_levels_ui <- renderUI({
      req(input$split_by, nzchar(input$split_by))
      pickerInput(ns("split_levels"), "Facet categories to show",
                  choices = NULL, multiple = TRUE)
    })
    
    ns <- session$ns
    
    # Build plotting dataframe defensively, with downsampling for rendering
    df <- reactive({
      coords_val       <- coords()
      expr_val         <- expr()
      meta_val         <- meta_cell()
      clusters_val     <- clusters()
      cluster_map_val  <- cluster_map()
      
      req(!is.null(coords_val))
      dd <- as.data.frame(coords_val)
      
      coords_full <- as.data.frame(coords())
      names(coords_full)[1:2] <- c("UMAP1", "UMAP2")
      coords_full$cluster <- factor(clusters()$assignments)  # factor for discrete colors
      
      color_text_add <- coords_full %>%
        group_by(cluster) %>%
        summarise(
          UMAP1 = mean(UMAP1, na.rm = TRUE),
          UMAP2 = mean(UMAP2, na.rm = TRUE),
          .groups = "drop"
        )
      
      if (ncol(dd) < 2 || nrow(dd) == 0) {
        message(sprintf("[%s] df(): coords invalid — ncol=%s nrow=%s",
                        embedding_name, ncol(dd), nrow(dd)))
        return(tibble::tibble(x = numeric(0), y = numeric(0), .cell = integer(0)))
      }
      
      names(dd)[1:2] <- c("x", "y")
      
      # IMPORTANT: .cell = original row number BEFORE downsampling
      dd$.cell <- seq_len(nrow(dd))
      
      # Add cluster and celltype columns
      if (!is.null(clusters_val) && !is.null(clusters_val$assignments)) {
        if (length(clusters_val$assignments) == nrow(dd)) {
          dd$cluster <- clusters_val$assignments
        } else {
          message(sprintf("[%s] df(): cluster length=%s != nrow(coords)=%s",
                          embedding_name, length(clusters_val$assignments), nrow(dd)))
          dd$cluster <- NA_integer_
        }
      } else {
        dd$cluster <- NA_integer_
      }
      
      if (!is.null(cluster_map_val) &&
          all(c("cluster", "celltype") %in% names(cluster_map_val))) {
        dd$celltype <- as.character(cluster_map_val$celltype[
          match(dd$cluster, cluster_map_val$cluster)
        ])
      } else {
        dd$celltype <- as.character(dd$cluster)
      }
      
      # ---- Downsample for plotting only ----
      if (nrow(dd) > 100000) {
        set.seed(123)  # reproducible sampling
        keep_idx <- sample.int(nrow(dd), 100000)
        dd <- dd[keep_idx, , drop = FALSE]
        message(sprintf("[%s] df(): downsampled to %d rows for plotting",
                        embedding_name, nrow(dd)))
      }
      
      tibble::as_tibble(dd)
    })
    
    current_sel <- reactiveVal(integer(0))
    
    # Lasso selection observer (suspended until first plot)
    sel_obs <- observeEvent(
      event_data("plotly_selected", source = ns("embed")),
      suspended = TRUE, ignoreInit = TRUE, {
        if (input$gate_mode != "Lasso select") return()
        sel <- event_data("plotly_selected", source = ns("embed"))
        req(!is.null(sel), nrow(sel) > 0)
        current_sel(unique(as.integer(sel$customdata)))
      }
    )
    
    # Polygon draw observer (suspended until first plot)
    rel_obs <- observeEvent(
      event_data("plotly_relayout", source = ns("embed")),
      suspended = TRUE, ignoreInit = TRUE, {
        if (input$gate_mode != "Draw polygon") return()
        ev <- event_data("plotly_relayout", source = ns("embed"))
        req(!is.null(ev), length(ev) > 0)
        shape_keys <- names(ev)[grepl("^shapes\\[\\d+\\]\\.path$", names(ev))]
        if (!length(shape_keys)) return()
        path <- ev[[shape_keys[length(shape_keys)]]]
        poly_df <- parse_svg_path(path)
        if (!nrow(poly_df)) return()
        
        # Build sf polygon and test points
        pg <- st_sfc(st_polygon(list(as.matrix(poly_df))), crs = NA)
        dd <- df()
        pts <- st_as_sf(dd, coords = c("x","y"))
        inside <- st_within(pts, pg, sparse = FALSE)[,1]
        current_sel(dd$.cell[inside])
      }
    )
    
    plot_cache_base <- reactiveVal(NULL)  # base plot without gates
    plot_cache      <- reactiveVal(NULL)  # final plot with overlays
    
    # Helper: clip values to reference distribution percentiles (1%–99%)
    clip_to_ref <- function(values, ref, probs = c(0.01, 0.99)) {
      qs <- stats::quantile(ref, probs = probs, na.rm = TRUE)
      pmin(pmax(values, qs[1]), qs[2])
    }
    
    # Add split_by to triggers for this observer
    observeEvent(list(df(), expr(), meta_cell(), input$color_by, input$gate_mode, input$max_facets, input$plot_facets), {
      expr_val <- expr()
      meta_val <- meta_cell()
      req(expr_val, meta_val)
      
      numeric_markers <- colnames(expr_val)
      meta_cols       <- colnames(meta_val)
      valid_cols      <- c(numeric_markers, meta_cols)
      
      color_by <- input$color_by
      if (is.null(color_by) || !(color_by %in% valid_cols)) {
        color_by <- if (length(numeric_markers)) numeric_markers[1] else meta_cols[1]
      }
      
      dd <- df()
      if (nrow(dd) == 0 || all(is.na(dd$x)) || all(is.na(dd$y))) {
        plot_cache_base(
          plotly_empty(type = "scatter", mode = "markers", source = ns("embed")) %>%
            layout(
              xaxis = list(title = paste0(embedding_name, " 1")),
              yaxis = list(title = paste0(embedding_name, " 2"))
            )
        )
        return()
      }
      
      # Colour mapping (unchanged)
      if (color_by %in% colnames(expr_val)) {
        vals_full <- expr_val[, color_by]
        vals_plot <- vals_full[dd$.cell]
        vals_plot_clipped <- clip_to_ref(vals_plot, vals_full)
        dom <- range(clip_to_ref(vals_full, vals_full), na.rm = TRUE)
        cols <- col_numeric(viridis(256), domain = dom)(vals_plot_clipped)
      } else {
        vals_full <- meta_val[[color_by]]
        vals_plot <- vals_full[dd$.cell]
        if (is.numeric(vals_full)) {
          vals_plot_clipped <- clip_to_ref(vals_plot, vals_full)
          dom <- range(clip_to_ref(vals_full, vals_full), na.rm = TRUE)
          cols <- col_numeric(viridis(256), domain = dom)(vals_plot_clipped)
        } else {
          levs <- unique(as.character(vals_full))
          pal  <- setNames(viridis(max(2, length(levs))), levs)
          cols <- pal[as.character(vals_plot)]
        }
      }
      
      # handle split_by
      # handle split_by
      split_var <- input$split_by
      if (!is.null(split_var) && nzchar(split_var)) {
        # Attach split var to plotting data
        dd[[split_var]] <- meta_val[[split_var]][dd$.cell]
        
        # Filter to selected categories (if provided)
        shown_levels <- input$split_levels
        if (!is.null(shown_levels) && length(shown_levels) > 0) {
          dd <- dd[dd[[split_var]] %in% shown_levels, , drop = FALSE]
        }
        
        # Drop NAs and empty levels
        dd <- dd[!is.na(dd[[split_var]]), , drop = FALSE]
        dd[[split_var]] <- droplevels(as.factor(dd[[split_var]]))
        
        # Guard: nothing left after filtering
        if (nrow(dd) == 0 || length(levels(dd[[split_var]])) == 0) {
          plot_cache_base(NULL)
          plot_cache(NULL)
          showNotification("No data points available for the selected facet categories.", type = "error")
          return()
        }
        
        # --- Balanced downsampling across shown, nonempty facets ---
        max_total <- 200000
        
        # Count only nonempty facets
        facet_counts <- table(dd[[split_var]])
        empty_levels <- names(facet_counts)[facet_counts == 0]
        
        # If any selected facets are empty, drop them and notify
        if (length(empty_levels) > 0) {
          dd <- dd[!(dd[[split_var]] %in% empty_levels), , drop = FALSE]
          dd[[split_var]] <- droplevels(dd[[split_var]])
          showNotification(
            paste0("Dropped empty facet categories: ", paste(empty_levels, collapse = ", ")),
            type = "warning"
          )
          facet_counts <- table(dd[[split_var]])  # recompute after drop
        }
        
        n_facets <- length(facet_counts)
        
        # Guard: all selected facets ended up empty
        if (n_facets == 0) {
          plot_cache_base(NULL)
          plot_cache(NULL)
          showNotification("No nonempty facet categories to plot after filtering.", type = "error")
          return()
        }
        
        # Target per facet from nonempty categories
        target_per_facet <- floor(max_total / n_facets)
        
        # Do not exceed the sparsest nonempty facet
        min_facet_size <- min(facet_counts)
        target_per_facet <- min(target_per_facet, min_facet_size)
        
        # Guard: target 0 means at least one chosen facet has 0 points after prior downsample
        if (target_per_facet < 1) {
          showNotification(
            "Selected facets have too few points to plot. Try including more categories or reducing filters.",
            type = "error"
          )
          plot_cache_base(NULL)
          plot_cache(NULL)
          return()
        }
        
        # Sample the same number from each remaining facet
        set.seed(123)
        dd <- dd %>%
          group_by(.data[[split_var]]) %>%
          slice_sample(n = target_per_facet) %>%
          ungroup()
        
        # Guard: ensure we still have data
        if (nrow(dd) == 0) {
          plot_cache_base(NULL)
          plot_cache(NULL)
          showNotification("No data to plot after balanced sampling.", type = "error")
          return()
        }
        
        # --- Color variable ---
        color_by <- input$color_by
        if (color_by %in% colnames(expr_val)) {
          dd[[color_by]] <- expr_val[, color_by][dd$.cell]
        } else {
          dd[[color_by]] <- meta_val[[color_by]][dd$.cell]
        }
        
        # Choose scale
        if (is.numeric(dd[[color_by]])) {
          color_scale <- scale_color_viridis_c()
        } else {
          color_scale <- scale_color_viridis_d()
        }
        
        # Build faceted plot
        p_base <- ggplot(dd, aes(x = x, y = y, color = .data[[color_by]])) +
          ggrastr::geom_point_rast(size = 0.25, alpha = 0.25) +
          guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
          facet_wrap(as.formula(paste("~", split_var)), ncol = input$max_facets) +
          color_scale +
          theme_minimal() +
          theme(legend.position = "right") +
          labs(
            x = paste0(embedding_name, " 1"),
            y = paste0(embedding_name, " 2"),
            color = color_by
          )
        
        # Optional: convert to plotly
        p_base <- plotly::ggplotly(p_base)
      } else {
        # Original single‑panel scatter
        p_base <- plot_ly(
          data = dd,
          x = ~x, y = ~y,
          type = "scatter", mode = "markers",
          marker = list(color = cols, size = 3),
          source = ns("embed"),
          customdata = ~.cell
        ) %>%
          layout(
            xaxis = list(title = paste0(embedding_name, " 1")),
            yaxis = list(title = paste0(embedding_name, " 2")),
            dragmode = if (input$gate_mode == "Lasso select") "lasso" else "zoom"
          )
      }
      
      # If labels are toggled on, add them immediately
      if (isTRUE(input$show_labels)) {
        coords_full <- as.data.frame(coords())
        names(coords_full)[1:2] <- c("x", "y")
        coords_full$cluster <- factor(clusters()$assignments)
        
        label_df <- coords_full %>%
          group_by(cluster) %>%
          summarise(
            x = mean(x, na.rm = TRUE),
            y = mean(y, na.rm = TRUE),
            .groups = "drop"
          )
        
        annots <- lapply(seq_len(nrow(label_df)), function(i) {
          list(
            x = label_df$x[i],
            y = label_df$y[i],
            xref = "x", yref = "y",
            text = as.character(label_df$cluster[i]),
            showarrow = FALSE,
            xanchor = "center", yanchor = "middle", align = "center",
            font = list(color = "black", size = 18),
            bgcolor = "rgba(255,255,255,0.85)",
            bordercolor = "rgba(0,0,0,0)",
            borderpad = 2, opacity = 1
          )
        })
        
        p_base <- p_base %>% layout(annotations = annots)
      }
      
      plot_cache_base(p_base)
      plot_cache(p_base)
    })
    
    # Update overlays (gates + labels) without rebuilding points
    observeEvent(list(input$overlay_gate, input$show_labels), {
      p <- plot_cache_base()
      req(p)
      
      # Add labels if toggled on
      if (isTRUE(input$show_labels)) {
        coords_full <- as.data.frame(coords())
        names(coords_full)[1:2] <- c("x", "y")
        coords_full$cluster <- factor(clusters()$assignments)
        
        label_df <- coords_full %>%
          group_by(cluster) %>%
          summarise(
            x = mean(x, na.rm = TRUE),
            y = mean(y, na.rm = TRUE),
            .groups = "drop"
          )
        
        annots <- lapply(seq_len(nrow(label_df)), function(i) {
          list(
            x = label_df$x[i],
            y = label_df$y[i],
            xref = "x",
            yref = "y",
            text = as.character(label_df$cluster[i]),
            showarrow = FALSE,
            xanchor = "center",
            yanchor = "middle",
            align = "center",
            font = list(color = "black", size = 18),
            bgcolor = "rgba(255,255,255,0.85)",
            bordercolor = "rgba(0,0,0,0)",
            borderpad = 2,
            opacity = 1
          )
        })
        
        p <- p %>% layout(annotations = annots)
      }
      
      # Add gate shapes
      gl <- gate_store$list()
      og <- input$overlay_gate
      if (length(gl) && length(og)) {
        shapes <- list()
        for (nm in og) {
          g <- gl[[nm]]
          if (!is.null(g) && !is.null(g$polygon) && g$embedding == embedding_name) {
            shapes <- c(shapes, list(
              type = "path",
              path = paste0(
                "M ",
                paste(paste(g$polygon$x, g$polygon$y, sep = ","), collapse = " L "),
                " Z"
              ),
              line = list(color = g$color, width = 2)
            ))
          }
        }
        if (length(shapes)) {
          p <- p %>% layout(shapes = shapes)
        }
      }
      
      plot_cache(p)
    })
    
    output$embed_plot <- renderPlotly({
      req(plot_cache())
      plot_cache()
    })
    
    # Keep the plot alive when hidden
    outputOptions(output, "embed_plot", suspendWhenHidden = FALSE)
    
    # Ensure the plot renders even when the tab is hidden (prevents some init races)
    outputOptions(output, "embed_plot", suspendWhenHidden = FALSE)
    
    observeEvent(input$clear_selection, current_sel(integer(0)))
    
    # Save gate
    observeEvent(input$save_gate, {
      cells <- current_sel()
      req(length(cells) > 0)
      name <- input$gate_name
      if (!nzchar(name)) name <- paste0(embedding_name, "_gate_", format(Sys.time(), "%H%M%S"))
      col <- input$gate_color
      polygon <- NULL
      gate <- new_gate(name, cells, embedding_name, polygon = polygon, color = col)
      
      gate_store$add(gate)
      updatePickerInput(session, "overlay_gate", choices = names(gate_store$list()),
                        selected = unique(c(input$overlay_gate, name)))
      updatePickerInput(session, "phenotype_gate", choices = names(gate_store$list()), selected = name)
      updatePickerInput(session, "test_gate", choices = names(gate_store$list()),
                        selected = unique(c(input$test_gate, name)))
      
      showNotification(paste("Saved gate:", name), type = "message")
    })
    
    # Phenotype heatmap
    output$inout_heatmap <- renderPlot({
      gname <- req(input$phenotype_gate)
      glist <- gate_store$list()
      req(gname %in% names(glist))
      gate <- glist[[gname]]
      idx <- gate$cells
      
      fun <- switch(input$phenotype_fun,
                    median = stats::median,
                    p90 = function(z) quantile(z, 0.9, na.rm = TRUE))
      
      in_vals <- apply(expr[idx, , drop = FALSE], 2, fun, na.rm = TRUE)
      out_vals <- apply(expr[-idx, , drop = FALSE], 2, fun, na.rm = TRUE)
      
      M <- rbind(In = in_vals, Out = out_vals)
      Mz <- t(scale(t(M)))
      
      Heatmap(Mz, name = "z",
              cluster_rows = FALSE, cluster_columns = TRUE,
              row_names_side = "left",
              col = colorRamp2(c(-2,0,2), viridis(3)))
    })
    
    # Keep gate pickers synced
    observe({
      updatePickerInput(session, "test_gate", choices = names(gate_store$list()))
      updatePickerInput(session, "overlay_gate", choices = names(gate_store$list()))
      updatePickerInput(session, "phenotype_gate", choices = names(gate_store$list()))
    })
  })
}


# ---- UI ----
ui <- navbarPage(
  "FCView",
  windowTitle = "FCView",
  id = "main_tab",
  
  tabPanel(
    "Home",
    sidebarLayout(
      sidebarPanel(
        h4("Upload inputs"),
        fileInput("rdata_upload", "Upload .RData (contains shinyAppInput)", accept = ".RData"),
        numericInput("max_cells_upload", "Max cells to read in", value = 300000, min = 1000, step = 1000),
        helpText("If the uploaded dataset has more cells than this number, it will be randomly downsampled at load time."), 
        fileInput("json_upload", "Upload gates (.json)", accept = ".json"),
        hr(),
        h4("Gate export"),
        uiOutput("gate_export_ui"),
        downloadButton("export_gates", "Download selected gates (.json)"),
        actionButton("clear_gates", "Clear all gates"),
        width = 3
      ),
      mainPanel(
        h3("Dataset overview"),
        verbatimTextOutput("ds_summary"),
        h4("Available metadata"),
        tableOutput("meta_overview"),
        h4("App capabilities"),
        tags$ul(
          tags$li(tags$b("Embeddings:"), " UMAP and tSNE with marker overlays and metadata faceting (shown if uploaded)"),
          tags$li(tags$b("Gating:"), " Lasso and drawn-polygon gating, multiple gates, gate import/export"),
          tags$li(tags$b("Phenotypes:"), " In-vs-out gate marker heatmaps"),
          tags$li(tags$b("Abundance testing:"), " Wilcoxon, Kruskal–Wallis, Spearman with BH correction"),
          tags$li(tags$b("Clusters:"), " Heatmap of cluster phenotypes when provided")
        )
      )
    )
  ),
  
  tabPanel("UMAP", EmbeddingUI("umap", title = "UMAP")),
  tabPanel("tSNE", EmbeddingUI("tsne", title = "tSNE")),
  
  tabPanel(
    "Clusters",
    fluidRow(
      column(
        3,
        pickerInput("cluster_order_by", "Order markers by",
                    choices = c("variance","none"), selected = "variance")
      ),
      column(
        9,
        plotOutput("cluster_heatmap", height = "650px")
      )
    ),
    hr(),
    h4("Abundance testing"),
    fluidRow(
      column(
        4,
        pickerInput("test_entity", "Entity",
                    choices = c("Clusters", "Celltypes")),  # no "Selected gate(s)" unless you want it here
        pickerInput("test_gate", "Gate(s)", choices = NULL, multiple = TRUE),
        pickerInput("group_var", "Grouping factor (metadata)",
                    choices = NULL,
                    options = list(`none-selected-text` = "None")), 
        pickerInput("cont_var", "Continuous metadata",
                    choices = NULL,
                    options = list(`none-selected-text` = "None")), 
        radioButtons("test_type", "Test",
                     choices = c("Wilcoxon (2-group)",
                                 "Kruskal–Wallis (multi-group)",
                                 "Spearman (continuous)")),
        pickerInput("unit_var", "Aggregation unit", choices = NULL),
        checkboxInput("apply_bh", "BH adjust across entities", TRUE),
        actionButton("run_test", "Run tests")
      ),
      column(
        8,
        plotOutput("abund_plot", height = "300px"),
        tableOutput("test_table")
      )
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  
  # App-wide stores
  rv <- reactiveValues(
    obj = NULL,
    expr = NULL,
    meta_cell = NULL,
    clusters = NULL,
    cluster_map = NULL,
    UMAP = NULL,
    tSNE = NULL,
    cluster_heat = NULL,
    pop_size = NULL,
    rep_used = NULL
  )
  
  gate_store <- GateStore()
  
  observeEvent(input$main_tab, {
    message("Tab changed to: ", input$main_tab)
  })
  
  # Upload RData and initialize datasets immediately (no mapping button)
  observeEvent(input$rdata_upload, {
    req(input$rdata_upload)
    
    # Load all objects from the uploaded .RData
    e <- new.env(parent = emptyenv())
    load(input$rdata_upload$datapath, envir = e)
    
    # Find all objects with the required structure
    candidates <- ls(e)
    matches <- vapply(candidates, function(nm) {
      cand <- e[[nm]]
      is.list(cand) && all(c("data", "source", "metadata", "leiden") %in% names(cand))
    }, logical(1))
    
    if (sum(matches) == 1) {
      obj <- e[[candidates[matches][1]]]
    } else if (sum(matches) > 1) {
      showNotification(
        paste("Multiple valid objects found in .RData:",
              paste(candidates[matches], collapse = ", ")),
        type = "error"
      )
      return()
    } else {
      showNotification("No valid input object found in .RData", type = "error")
      return()
    }
    
    # Ensure embeddings are data.frames
    if ("tsne" %in% names(obj) && is.matrix(obj$tsne$coordinates)) {
      obj$tsne$coordinates <- as.data.frame(obj$tsne$coordinates)
    }
    if ("umap" %in% names(obj) && is.matrix(obj$umap$coordinates)) {
      obj$umap$coordinates <- as.data.frame(obj$umap$coordinates)
    }
    
    # --- Downsample cells if needed (NEVER downsample sample-level metadata) ---
    n_cells <- nrow(obj$data)
    max_cells <- input$max_cells_upload %||% 300000
    
    if (n_cells > max_cells) {
      set.seed(123)
      keep_idx <- sort(sample.int(n_cells, max_cells))
      
      obj$data   <- obj$data[keep_idx, , drop = FALSE]
      obj$source <- obj$source[keep_idx]
      
      if (!is.null(obj$leiden$clusters) && length(obj$leiden$clusters) == n_cells) {
        obj$leiden$clusters <- obj$leiden$clusters[keep_idx]
      }
      if (hasUMAP(obj) && nrow(obj$umap$coordinates) == n_cells) {
        obj$umap$coordinates <- obj$umap$coordinates[keep_idx, , drop = FALSE]
      }
      if (hasTSNE(obj) && nrow(obj$tsne$coordinates) == n_cells) {
        obj$tsne$coordinates <- obj$tsne$coordinates[keep_idx, , drop = FALSE]
      }
      if (!is.null(obj$run_date) && length(obj$run_date) == n_cells) {
        obj$run_date <- obj$run_date[keep_idx]
      }
      
      message(sprintf("Downsampled from %d to %d cells at upload stage", n_cells, max_cells))
      showNotification(sprintf("Downsampled from %d to %d cells", n_cells, max_cells), type = "message")
    }
    
    # Validate the (possibly downsampled) object
    validateInput(obj, id_col = NULL)
    
    expr <- obj$data
    run_date <- obj$run_date %||% NULL
    
    # --- Build meta_cell with robust patient_ID mapping ---
    meta_cell <- data.frame(
      source  = as.character(obj$source),
      RunDate = if (!is.null(run_date)) run_date else NA
    )
    
    # Prepare metadata IDs (character, unique, deduplicated rows)
    if (!("patient_ID" %in% colnames(obj$metadata))) {
      showNotification("metadata does not contain 'patient_ID' column.", type = "error")
      return()
    }
    obj$metadata$patient_ID <- as.character(obj$metadata$patient_ID)
    
    metadata_unique <- obj$metadata %>%
      dplyr::distinct(patient_ID, .keep_all = TRUE)
    
    ids <- unique(metadata_unique$patient_ID)
    if (length(ids) == 0) {
      showNotification("No patient_ID values found in metadata.", type = "error")
      return()
    }
    
    # Escape any regex metacharacters in IDs, then sort by length (longest first)
    ids_escaped <- stringr::str_replace_all(ids, "([\\^$.|?*+()\\[\\]{}\\\\])", "\\\\\\1")
    ids_escaped <- ids_escaped[order(nchar(ids_escaped), decreasing = TRUE)]
    pattern <- paste0("(", paste0(ids_escaped, collapse = "|"), ")")
    
    # Extract patient_ID from source via regex built from metadata IDs
    meta_cell$patient_ID <- stringr::str_extract(meta_cell$source, pattern)
    
    # Post-extraction sanity check (IDs extracted but not present in metadata)
    not_in_meta <- setdiff(unique(na.omit(meta_cell$patient_ID)), metadata_unique$patient_ID)
    if (length(not_in_meta) > 0) {
      showNotification(
        paste0("Extracted patient_IDs not found in metadata: ", paste(not_in_meta, collapse = ", ")),
        type = "warning",
        duration = NULL
      )
    }
    
    # Join metadata by patient_ID (safe: metadata_unique has 1 row per ID)
    meta_cell <- meta_cell %>%
      dplyr::left_join(metadata_unique, by = "patient_ID")
    
    # --- Post-join check for unmatched IDs ---
    unmatched <- sum(is.na(meta_cell$patient_ID))
    if (unmatched > 0) {
      showNotification(
        paste0("Warning: ", unmatched, " cells could not be matched to metadata (patient_ID missing)."),
        type = "warning",
        duration = NULL
      )
      message("First few unmatched source entries:")
      print(utils::head(meta_cell$source[is.na(meta_cell$patient_ID)], 10))
    }
    
    # Ensure PatientID column exists for downstream code
    if (!("PatientID" %in% names(meta_cell))) {
      meta_cell$PatientID <- meta_cell$patient_ID
    }
    
    # Factor RunDate if present
    if ("RunDate" %in% names(meta_cell)) {
      meta_cell$RunDate <- as.factor(meta_cell$RunDate)
    }
    
    # Add leiden_cluster factor if available
    if (!"leiden_cluster" %in% names(meta_cell) &&
        !is.null(obj$leiden$clusters) &&
        length(obj$leiden$clusters) == nrow(obj$data)) {
      meta_cell$leiden_cluster <- factor(obj$leiden$clusters,
                                         levels = sort(unique(obj$leiden$clusters)))
    }
    
    clusters <- list(
      assignments = obj$leiden$clusters,
      settings    = obj$leiden$settings %||% list()
    )
    cluster_map <- if (hasClusterMapping(obj)) obj$cluster_mapping else NULL
    
    UMAP <- if (hasUMAP(obj)) list(coords = obj$umap$coordinates,
                                   settings = obj$umap$settings) else NULL
    tSNE <- if (hasTSNE(obj)) list(coords = obj$tsne$coordinates,
                                   settings = obj$tsne$settings) else NULL
    
    cluster_heat <- if (hasHeatmap(obj)) obj$leiden_heatmap$heatmap_tile_data else NULL
    pop_size     <- if (hasHeatmap(obj)) obj$leiden_heatmap$population_size else NULL
    rep_used     <- if (hasHeatmap(obj)) obj$leiden_heatmap$rep_used else NA
    
    # Store in rv
    rv$expr         <- expr
    rv$meta_cell    <- meta_cell
    rv$clusters     <- clusters
    rv$cluster_map  <- cluster_map
    rv$UMAP         <- UMAP
    rv$tSNE         <- tSNE
    rv$cluster_heat <- cluster_heat
    rv$pop_size     <- pop_size
    rv$rep_used     <- rep_used
    
    message("Upload complete: expr rows=", nrow(rv$expr),
            " meta_cell rows=", nrow(rv$meta_cell),
            " UMAP coords=", if (!is.null(rv$UMAP)) nrow(rv$UMAP$coords) else "NULL",
            " tSNE coords=", if (!is.null(rv$tSNE)) nrow(rv$tSNE$coords) else "NULL")
    
    showNotification("Data loaded and initialized.", type = "message")
  })
  
  observe({
    req(rv$meta_cell)
    
    meta_cols <- colnames(rv$meta_cell)
    cont_choices <- meta_cols[sapply(rv$meta_cell, is.numeric)]
    
    unit_candidates <- intersect(
      c("PatientID", "patient_ID", "patient", "source", "RunDate", "run_date"),
      meta_cols
    )
    unit_default <- if (length(unit_candidates)) unit_candidates[1] else meta_cols[1]
    
    updatePickerInput(session, "group_var",
                      choices = c("", meta_cols),
                      selected = "")
    updatePickerInput(session, "cont_var",
                      choices = c("", cont_choices),
                      selected = "")
    updatePickerInput(session, "unit_var",
                      choices = meta_cols,
                      selected = unit_default)
  })
  
  # Launch embedding modules as soon as data is ready (Option 1)
  observeEvent(rv$UMAP, {
    req(rv$UMAP, rv$expr, rv$meta_cell, rv$clusters)
    message("Launching UMAP module")
    EmbeddingServer(
      "umap", "UMAP",
      reactive(rv$UMAP$coords),
      reactive(rv$expr),
      reactive(rv$meta_cell),
      reactive(rv$clusters),
      reactive(rv$cluster_map),
      gate_store,
      reactive(input$main_tab),
      rv
    )
  })
  
  observeEvent(rv$tSNE, {
    req(rv$tSNE, rv$expr, rv$meta_cell, rv$clusters)
    message("Launching tSNE module")
    EmbeddingServer(
      "tsne", "tSNE",
      reactive(rv$tSNE$coords),
      reactive(rv$expr),
      reactive(rv$meta_cell),
      reactive(rv$clusters),
      reactive(rv$cluster_map),
      gate_store,
      reactive(input$main_tab),
      rv
    )
  })
  
  # Gate JSON upload
  observeEvent(input$json_upload, {
    req(input$json_upload)
    gates <- tryCatch({
      fromJSON(input$json_upload$datapath, simplifyVector = FALSE)
    }, error = function(e) {
      showNotification(paste("Failed to read JSON:", e$message), type = "error"); NULL
    })
    req(!is.null(gates))
    
    n_before <- length(gate_store$list())
    # gates can be a named list or unnamed list of gate objects
    if (is.list(gates) && !is.null(names(gates))) {
      for (nm in names(gates)) gate_store$add(gates[[nm]])
    } else if (is.list(gates)) {
      for (g in gates) gate_store$add(g)
    }
    n_after <- length(gate_store$list())
    showNotification(paste("Loaded", n_after - n_before, "gate(s) from JSON"), type = "message")
  })
  
  # Gate export UI
  output$gate_export_ui <- renderUI({
    gl <- gate_store$list()
    pickerInput("gates_to_export", "Select gates to export",
                choices = names(gl), multiple = TRUE)
  })
  
  # Gate export download
  output$export_gates <- downloadHandler(
    filename = function() paste0("polygon_gates_", format(Sys.Date(), "%Y%m%d"), ".json"),
    content = function(file) {
      sel <- input$gates_to_export
      gl <- gate_store$list()
      if (!length(sel)) {
        write_json(list(), path = file, auto_unbox = TRUE, pretty = TRUE)
      } else {
        write_json(gl[sel], path = file, auto_unbox = TRUE, pretty = TRUE)
      }
    }
  )
  
  observeEvent(input$clear_gates, {
    gate_store$clear()
    showNotification("Cleared all gates.", type = "message")
  })
  
  # Home summaries
  output$ds_summary <- renderPrint({
    req(rv$expr, rv$clusters)
    list(
      n_cells = nrow(rv$expr),
      n_markers = ncol(rv$expr),
      n_clusters = length(unique(rv$clusters$assignments)),
      embeddings_available = c(UMAP = !is.null(rv$UMAP), tSNE = !is.null(rv$tSNE)),
      has_cluster_heatmap = !is.null(rv$cluster_heat),
      rep_used = rv$rep_used
    )
  })
  
  # Metadata overview
  output$meta_overview <- renderTable({
    req(rv$meta_cell)
    data.frame(
      name = colnames(rv$meta_cell),
      type = sapply(rv$meta_cell, function(x) class(x)[1]),
      example = sapply(rv$meta_cell, function(x) paste(utils::head(unique(x), 3), collapse = ", "))
    )
  }, sanitize.text.function = function(x) x)
  
  # Cluster phenotype heatmap
  output$cluster_heatmap <- renderPlot({
    req(rv$cluster_heat)
    M <- rv$cluster_heat
    if (!is.null(input$cluster_order_by) && input$cluster_order_by == "variance") {
      ord <- order(apply(M, 2, stats::var, na.rm = TRUE), decreasing = TRUE)
      M <- M[, ord, drop = FALSE]
    }
    ranno <- NULL
    if (!is.null(rv$pop_size)) {
      ranno <- rowAnnotation(Size = rv$pop_size[,1])
    }
    Heatmap(
      M, name = "expr",
      cluster_rows = TRUE, cluster_columns = TRUE,
      right_annotation = ranno,
      row_names_side = "left",
      col = colorRamp2(c(min(M), median(M), max(M)), viridis(3))
    )
  })
  
  run_tests <- eventReactive(input$run_test, {
    test_type <- input$test_type
    
    # --- Gate-based testing (unchanged) ---
    if (input$test_entity == "Selected gate(s)") {
      req(length(input$test_gate) > 0)
      unit_var  <- req(input$unit_var)
      group_var <- input$group_var
      
      tests <- lapply(input$test_gate, function(gn) {
        gate <- gate_store$list()[[gn]]
        dfreq <- freq_by(gate$cells, rv$meta_cell, group_var, unit_var)
        
        if (test_type == "Wilcoxon (2-group)") {
          if (!nzchar(group_var)) return(data.frame(entity = gn, test = "wilcox", p = NA, n = nrow(dfreq)))
          g <- droplevels(factor(dfreq[[group_var]]))
          if (length(levels(g)) != 2) return(data.frame(entity = gn, test = "wilcox", p = NA, n = nrow(dfreq)))
          wt <- wilcox.test(freq ~ g, data = dfreq)
          data.frame(entity = gn, test = "wilcox", p = wt$p.value, n = nrow(dfreq))
          
        } else if (test_type == "Kruskal–Wallis (multi-group)") {
          if (!nzchar(group_var)) return(data.frame(entity = gn, test = "kruskal", p = NA, n = nrow(dfreq)))
          kw <- kruskal.test(freq ~ dfreq[[group_var]], data = dfreq)
          data.frame(entity = gn, test = "kruskal", p = kw$p.value, n = nrow(dfreq))
          
        } else {
          cont <- input$cont_var
          if (!nzchar(cont)) return(data.frame(entity = gn, test = "spearman", rho = NA, p = NA, n = nrow(dfreq)))
          ct <- spearman_test(dfreq, cont_var = cont)
          cbind(data.frame(entity = gn, test = "spearman"), ct)
        }
      })
      
      out <- do.call(rbind, tests)
      if (nrow(out) && isTRUE(input$apply_bh) && "p" %in% colnames(out))
        out$padj <- p.adjust(out$p, method = "BH")
      out
      
    } else {
      # --- Cluster or celltype testing using abundance matrix ---
      abund <- rv$clusters$abundance
      req(!is.null(abund), nrow(abund) > 0)
      
      # If both group_var and cont_var are blank, return NA row
      if (!nzchar(input$group_var) && !nzchar(input$cont_var)) {
        return(data.frame(entity = NA, test = NA, p = NA, n = NA))
      }
      
      sources <- rownames(abund)
      ids <- unique(rv$meta_cell$patient_ID)
      ids <- ids[order(nchar(ids), decreasing = TRUE)]
      pattern <- paste0("(", paste0(ids, collapse = "|"), ")")
      pid <- stringr::str_extract(sources, pattern)
      
      abund_df <- as.data.frame(abund)
      abund_df$patient_ID <- pid
      meta_unique <- rv$meta_cell %>%
        dplyr::distinct(patient_ID, .keep_all = TRUE)
      abund_df <- abund_df %>%
        dplyr::left_join(meta_unique, by = "patient_ID")
      
      abund_long <- abund_df %>%
        tidyr::pivot_longer(cols = colnames(abund),
                            names_to = "entity",
                            values_to = "freq")
      
      res <- abund_long %>%
        dplyr::group_by(entity) %>%
        dplyr::group_modify(~ {
          if (test_type == "Wilcoxon (2-group)") {
            if (!nzchar(input$group_var)) return(data.frame(test = "wilcox", p = NA, n = nrow(.x)))
            g <- droplevels(factor(.x[[input$group_var]]))
            if (length(levels(g)) != 2) return(data.frame(test = "wilcox", p = NA, n = nrow(.x)))
            wt <- wilcox.test(freq ~ g, data = .x)
            data.frame(test = "wilcox", p = wt$p.value, n = nrow(.x))
            
          } else if (test_type == "Kruskal–Wallis (multi-group)") {
            if (!nzchar(input$group_var)) return(data.frame(test = "kruskal", p = NA, n = nrow(.x)))
            kw <- kruskal.test(freq ~ .x[[input$group_var]], data = .x)
            data.frame(test = "kruskal", p = kw$p.value, n = nrow(.x))
            
          } else {
            if (!nzchar(input$cont_var)) return(data.frame(test = "spearman", rho = NA, p = NA, n = nrow(.x)))
            cont <- input$cont_var
            ct <- spearman_test(.x, cont_var = cont)
            cbind(data.frame(test = "spearman"), ct)
          }
        }) %>%
        dplyr::ungroup()
      
      if (nrow(res) && isTRUE(input$apply_bh) && "p" %in% colnames(res))
        res$padj <- p.adjust(res$p, method = "BH")
      
      res
    }
  })
  
  # --- Plot with message if no variable selected ---
  output$abund_plot <- renderPlot({
    res <- req(run_tests())
    
    # If both vars are None, show message instead of plot
    if (!nzchar(input$group_var) && !nzchar(input$cont_var)) {
      grid::grid.newpage()
      grid::grid.text("Select a grouping or continuous variable to run tests",
                      gp = grid::gpar(fontsize = 14))
      return()
    }
    
    if (!nrow(res)) return(NULL)
    
    if ("rho" %in% names(res)) {
      ggplot(res, aes(x = n, y = rho, color = p)) +
        geom_point() +
        scale_color_viridis_c() +
        theme_minimal() +
        labs(title = "Spearman results", x = "N", y = "rho")
    } else {
      ggplot(res, aes(x = entity %||% "", y = -log10(p))) +
        geom_col(fill = "#2c7fb8") +
        theme_minimal() +
        labs(title = "Abundance test -log10(p)", x = "", y = "-log10(p)")
    }
  })
  
  output$test_table <- renderTable({
    req(run_tests())
  })
}

shinyApp(ui, server)
