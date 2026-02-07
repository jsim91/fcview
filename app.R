# App used to explore FCSimple results
# run by: shiny::runApp(appDir = 'app.R_directory')
# where 'app.R_directory' is the location of this script, app.R

# ---- Packages ----
suppressPackageStartupMessages({
  suppressWarnings({
    library(shiny)
    library(shinyWidgets)
    library(plotly)
    library(ggplot2)
    library(ggrastr)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(ComplexHeatmap)
    library(circlize)
    library(viridis)
    library(scales)
    library(glmnet)
    library(Boruta)
    library(colourpicker)
    library(caret) # train(), trainControl, resampling
    library(nnet) # multinom() for multi-class logistic regression
    library(pROC)
    library(yardstick)
    library(rsample) # vfold_cv(), initial_split - tidy resampling
    library(boot) # bootstrapping CIs for AUC
    library(rlang) # tidy evaluation for !!!syms
    library(randomForest)
    library(multiROC) # for multiclass ROC curves
    library(zip) # for ZIP download handlers
    library(parallel)
    library(patchwork)
    library(sccomp)
    library(survival) # Cox proportional hazards, survfit, cox.zph
    library(survminer) # ggsurvplot for better survival curve visualization
  })
})

options(shiny.maxRequestSize = 5000 * 1024^2)

# ---- Helpers & validation ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

pct_clip <- function(x, p = c(0.01, 0.99)) {
  q <- quantile(x, probs = p, na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

# appendLog <- function(msg) {
#   old <- rv$log()
#   new <- c(old, paste0(format(Sys.time(), "%H:%M:%S"), " | ", msg))
#   if (length(new) > 10) new <- tail(new, 10)  # keep last 10
#   rv$log(new)
# }

resetResults <- function(resultReactive, cacheReactive = NULL) {
  resultReactive(NULL)
  if (!is.null(cacheReactive)) cacheReactive(NULL)
}

align_metadata_abundance <- function(metadata, abundance, notify = NULL) {
  # Extract patient_ID from abundance rownames
  patient_ids <- gsub(
    pattern = "_[0-9]+\\-[A-Za-z]+\\-[0-9]+.*$",
    replacement = "",
    x = rownames(abundance)
  )

  abund_df <- as.data.frame(abundance)
  abund_df$patient_ID <- patient_ids

  # Duplicate checks
  if (anyDuplicated(metadata$patient_ID)) {
    msg <- "Duplicate patient_IDs found in metadata. Consider deduplicating with distinct()."
    warning(msg)
    if (!is.null(notify)) notify(msg, type = "warning")
  }
  if (anyDuplicated(abund_df$patient_ID)) {
    msg <- "Duplicate patient_IDs found in abundance rownames. Check your input files."
    warning(msg)
    if (!is.null(notify)) notify(msg, type = "warning")
  }

  # Mismatch checks
  meta_ids <- unique(metadata$patient_ID)
  abund_ids <- unique(abund_df$patient_ID)

  missing_in_abund <- setdiff(meta_ids, abund_ids)
  missing_in_meta <- setdiff(abund_ids, meta_ids)

  if (length(missing_in_abund) > 0) {
    msg <- paste(
      "These patient_IDs are in metadata but not in abundance:",
      paste(missing_in_abund, collapse = ", ")
    )
    warning(msg)
    if (!is.null(notify)) notify(msg, type = "warning")
  }
  if (length(missing_in_meta) > 0) {
    msg <- paste(
      "These patient_IDs are in abundance but not in metadata:",
      paste(missing_in_meta, collapse = ", ")
    )
    warning(msg)
    if (!is.null(notify)) notify(msg, type = "warning")
  }

  # Merge with metadata (metadata is the anchor)
  merged <- dplyr::left_join(metadata, abund_df, by = "patient_ID")

  return(merged)
}

clean_dummy_names <- function(nms) {
  gsub(pattern = "cluster", replacement = "cluster:", x = nms)
}

spearman_test <- function(df, freq_col = "freq", cont_var) {
  ok <- complete.cases(df[[freq_col]], df[[cont_var]])
  if (!any(ok)) {
    return(data.frame(rho = NA_real_, p = NA_real_, n = 0))
  }
  ct <- suppressWarnings(cor.test(df[[freq_col]][ok], df[[cont_var]][ok], method = "spearman"))
  data.frame(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
}

compute_multiROC <- function(preds) {
  if (!is.factor(preds$obs)) preds$obs <- factor(preds$obs)
  class_levels <- levels(preds$obs)
  if (length(class_levels) < 2) stop("Need ≥2 classes for ROC")

  # Build truth and prob data.frames
  truth_df <- as.data.frame(
    sapply(class_levels, function(cls) as.integer(preds$obs == cls), simplify = "matrix")
  )
  names(truth_df) <- paste0(class_levels, "_true")

  prob_df <- preds[, paste0(".pred_", class_levels), drop = FALSE]
  names(prob_df) <- paste0(class_levels, "_pred")

  # Combine
  multi_df <- cbind(truth_df, prob_df)

  # Force to plain numeric matrix
  multi_mat <- as.matrix(multi_df)
  storage.mode(multi_mat) <- "numeric"

  message("compute_multiROC: multi_mat head")
  print(utils::head(multi_mat))

  # Run multiROC on the matrix
  roc_res <- multiROC::multi_roc(multi_mat, force_diag = TRUE)
  attr(roc_res, "multi_df") <- multi_df
  roc_res
}

# Validate minimal structure
validateInput <- function(obj, id_col = NULL) {
  required <- c("data", "source", "metadata", "cluster")
  miss <- setdiff(required, names(obj))
  if (length(miss)) stop("Missing required elements: ", paste(miss, collapse = ", "))

  if (!is.matrix(obj$data)) stop("data must be a matrix (cells × features).")
  n <- nrow(obj$data)

  if (length(obj$source) != n) stop("source length must equal nrow(data).")
  if (!is.data.frame(obj$metadata)) stop("metadata must be a data.frame.")
  if (!is.list(obj$cluster) || is.null(obj$cluster$clusters)) stop("cluster must be a list with 'clusters'.")
  if (length(obj$cluster$clusters) != n) stop("cluster$clusters length must equal nrow(data).")

  # Optional elements: if present, check shape
  if ("umap" %in% names(obj)) {
    if (!is.data.frame(obj$umap$coordinates) || nrow(obj$umap$coordinates) != n) {
      stop("umap$coordinates must align with cells.")
    }
  }
  if ("tsne" %in% names(obj)) {
    if (!is.data.frame(obj$tsne$coordinates) || nrow(obj$tsne$coordinates) != n) {
      stop("tsne$coordinates must align with cells.")
    }
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
hasHeatmap <- function(obj) "cluster_heatmap" %in% names(obj) && !is.null(obj$cluster_heatmap$heatmap_tile_data)
hasRunDate <- function(obj) "run_date" %in% names(obj) && length(obj$run_date) == nrow(obj$data)
hasClusterMapping <- function(obj) "cluster_mapping" %in% names(obj)

# Auto-detect metadata ID column that matches source
guess_id_col <- function(metadata, source_vec) {
  # Prioritized candidates by name
  candidates <- intersect(
    tolower(colnames(metadata)),
    c("patientid", "patient_id", "patient", "source", "subject", "sample", "id")
  )
  if (length(candidates)) {
    hit <- candidates[1]
    return(colnames(metadata)[tolower(colnames(metadata)) == hit][1])
  }
  # Fallback: choose column with max overlap
  overlaps <- sapply(metadata, function(col) {
    if (!is.atomic(col)) {
      return(0)
    }
    length(intersect(as.character(col), as.character(source_vec)))
  })
  if (all(overlaps == 0)) {
    return(NULL)
  }
  colnames(metadata)[which.max(overlaps)]
}

build_design_matrix <- function(df) {
  # model.matrix will one-hot encode factors/characters; drop intercept
  mm <- model.matrix(~ . - 1, data = df)
  # Ensure a plain numeric matrix
  storage.mode(mm) <- "double"
  mm
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
        pickerInput(ns("split_by"), "Facet by",
          choices = NULL,
          options = list(`none-selected-text` = "None")
        ),
        uiOutput(ns("split_levels_ui")),
        selectInput(
          ns("max_facets"), "Facet columns",
          choices = c(1, 2, 3, 4), selected = 2
        ),
        actionButton(ns("plot_facets"), "Plot facets"),
        hr(),
        # Moved export button here:
        downloadButton(ns("export_embed_pdf"), "Export embedding as PDF"),
        hr()
      ),
      column(
        9,
        plotlyOutput(ns("embed_plot"), height = "650px")
      )
    )
  )
}

# ---- Embedding module server ----
EmbeddingServer <- function(id, embedding_name, coords, expr, meta_cell, clusters, cluster_map,
                            active_tab) {
  moduleServer(id, function(input, output, session) {
    message(sprintf("EmbeddingServer %s started", embedding_name))
    message(sprintf(
      "coords NULL? %s | expr NULL? %s | meta_cell NULL? %s",
      is.null(coords), is.null(expr), is.null(meta_cell)
    ))

    # One-time picker initialization in a reactive context
    initialized <- FALSE
    plot_cache_gg <- reactiveVal(NULL)
    facet_cols_saved <- reactiveVal(2)

    observeEvent(list(expr(), meta_cell(), clusters()), ignoreInit = FALSE, {
      if (initialized) {
        return()
      }
      expr_val <- expr()
      meta_val <- meta_cell()
      clusters_val <- clusters()

      req(!is.null(expr_val), !is.null(meta_val))

      # Add cluster column as factor if missing
      if (!"cluster" %in% colnames(meta_val) && !is.null(clusters_val$assignments)) {
        meta_val$cluster <- factor(clusters_val$assignments)
      }

      numeric_markers <- colnames(expr_val)
      meta_cols <- setdiff(colnames(meta_val), c(".cell"))

      # Ensure cluster is present and comes right after markers
      if (!"cluster" %in% meta_cols && "cluster" %in% colnames(meta_val)) {
        meta_cols <- c("cluster", setdiff(meta_cols, "cluster"))
      }

      # Sort metadata portion alphabetically (excluding cluster)
      meta_cols_sorted <- sort(setdiff(meta_cols, "cluster"))

      # Final order: markers → cluster → sorted metadata
      ordered_choices <- c(numeric_markers, "cluster", meta_cols_sorted)

      # Continuous and categorical choices from metadata
      cont_choices <- sort(meta_cols[sapply(meta_val[meta_cols], is.numeric)])
      factor_cols <- meta_cols[sapply(meta_val[meta_cols], is.factor)]
      char_cols <- meta_cols[sapply(meta_val[meta_cols], is.character)]
      categorical_choices <- sort(c(factor_cols, char_cols))

      # Default unit_var
      unit_default <- if ("patient_ID" %in% meta_cols) "patient_ID" else meta_cols[1]

      updatePickerInput(session, "color_by",
        choices = ordered_choices,
        selected = if (length(numeric_markers)) numeric_markers[1] else ordered_choices[1]
      )
      updatePickerInput(session, "split_by", choices = c("", categorical_choices), selected = "")
      updatePickerInput(session, "group_var", choices = c("", sort(meta_cols)), selected = "")
      updatePickerInput(session, "cont_var", choices = c("", cont_choices), selected = "")
      updatePickerInput(session, "unit_var", choices = sort(meta_cols), selected = unit_default)

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
        choices = NULL, multiple = TRUE
      )
    })

    ns <- session$ns

    # Build plotting dataframe defensively, with downsampling for rendering
    df <- reactive({
      coords_val <- coords()
      expr_val <- expr()
      meta_val <- meta_cell()
      clusters_val <- clusters()
      cluster_map_val <- cluster_map()

      req(!is.null(coords_val))
      dd <- as.data.frame(coords_val)

      coords_full <- as.data.frame(coords())
      names(coords_full)[1:2] <- c("UMAP1", "UMAP2")
      coords_full$cluster <- factor(clusters()$assignments) # factor for discrete colors

      color_text_add <- coords_full %>%
        group_by(cluster) %>%
        summarise(
          UMAP1 = mean(UMAP1, na.rm = TRUE),
          UMAP2 = mean(UMAP2, na.rm = TRUE),
          .groups = "drop"
        )

      if (ncol(dd) < 2 || nrow(dd) == 0) {
        message(sprintf(
          "[%s] df(): coords invalid — ncol=%s nrow=%s",
          embedding_name, ncol(dd), nrow(dd)
        ))
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
          message(sprintf(
            "[%s] df(): cluster length=%s != nrow(coords)=%s",
            embedding_name, length(clusters_val$assignments), nrow(dd)
          ))
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
        set.seed(123) # reproducible sampling
        keep_idx <- sample.int(nrow(dd), 100000)
        dd <- dd[keep_idx, , drop = FALSE]
        message(sprintf(
          "[%s] df(): downsampled to %d rows for plotting",
          embedding_name, nrow(dd)
        ))
      }

      tibble::as_tibble(dd)
    })

    current_sel <- reactiveVal(integer(0))

    plot_cache_base <- reactiveVal(NULL) # base plot
    plot_cache <- reactiveVal(NULL) # final plot with overlays

    # Helper: clip values to reference distribution percentiles (1%–99%)
    clip_to_ref <- function(values, ref, probs = c(0.01, 0.99)) {
      qs <- stats::quantile(ref, probs = probs, na.rm = TRUE)
      pmin(pmax(values, qs[1]), qs[2])
    }

    # Add split_by to triggers for this observer
    observeEvent(
      list(df(), expr(), meta_cell(), input$color_by, input$plot_facets),
      {
        expr_val <- expr()
        meta_val <- meta_cell()
        req(expr_val, meta_val)

        numeric_markers <- colnames(expr_val)
        meta_cols <- colnames(meta_val)
        valid_cols <- c(numeric_markers, meta_cols)

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

        # --- Always add the colour column to dd ---
        if (color_by %in% numeric_markers) {
          dd$.color_val <- as.numeric(expr_val[dd$.cell, color_by])
        } else if (color_by %in% meta_cols) {
          dd$.color_val <- meta_val[[color_by]][dd$.cell]
        } else {
          dd$.color_val <- NA
        }

        # Handle faceting
        split_var <- input$split_by
        if (!is.null(split_var) && nzchar(split_var)) {
          dd[[split_var]] <- meta_val[[split_var]][dd$.cell]
          shown_levels <- input$split_levels
          if (!is.null(shown_levels) && length(shown_levels) > 0) {
            dd <- dd[dd[[split_var]] %in% shown_levels, , drop = FALSE]
          }
          dd <- dd[!is.na(dd[[split_var]]), , drop = FALSE]
          dd[[split_var]] <- droplevels(as.factor(dd[[split_var]]))
          if (nrow(dd) == 0) {
            plot_cache_base(NULL)
            plot_cache(NULL)
            showNotification("No data points available for the selected facet categories.", type = "error")
            return()
          }

          # Balanced downsampling
          max_total <- input$max_cells_upload %||% 100000
          facet_counts <- table(dd[[split_var]])
          n_facets <- length(facet_counts)
          target_per_facet <- floor(max_total / n_facets)
          target_per_facet <- min(target_per_facet, min(facet_counts))
          if (target_per_facet < 1) {
            showNotification("Selected facets have too few points to plot.", type = "error")
            plot_cache_base(NULL)
            plot_cache(NULL)
            return()
          }
          set.seed(123)
          dd <- dd %>%
            group_by(.data[[split_var]]) %>%
            slice_sample(n = target_per_facet) %>%
            ungroup()

          # Build faceted ggplot
          gg <- ggplot(dd, aes(x = x, y = y, color = .data[[".color_val"]])) +
            ggrastr::geom_point_rast(size = 0.25, alpha = 0.5) +
            guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
            facet_wrap(as.formula(paste("~", split_var)), ncol = as.numeric(input$max_facets)) +
            (if (is.numeric(dd$.color_val)) scale_color_viridis_c() else scale_color_viridis_d()) +
            theme_minimal() +
            guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
            theme(legend.position = "right") +
            labs(
              x = paste0(embedding_name, " 1"),
              y = paste0(embedding_name, " 2"),
              color = color_by
            )

          plot_cache_gg(gg)
          p_base <- ggplotly(gg, tooltip = "text") %>% layout(dragmode = "pan")
        } else {
          # Single-panel ggplot for export
          gg <- ggplot(dd, aes(x = x, y = y, color = .data[[".color_val"]])) +
            geom_point(size = 0.25, alpha = 0.25) +
            (if (is.numeric(dd$.color_val)) scale_color_viridis_c() else scale_color_viridis_d()) +
            theme_minimal() +
            guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
            labs(
              x = paste0(embedding_name, " 1"),
              y = paste0(embedding_name, " 2"),
              color = color_by
            )
          plot_cache_gg(gg)

          # Compute colours for plotly
          if (is.numeric(dd$.color_val)) {
            rng <- range(dd$.color_val, na.rm = TRUE)
            if (!is.finite(rng[1]) || !is.finite(rng[2])) rng <- c(0, 1)
            cols <- scales::col_numeric(
              viridis::viridis(256),
              domain = rng
            )(dd$.color_val)
          } else {
            vals <- as.character(dd$.color_val)
            levs <- sort(unique(vals))
            pal <- viridis::viridis(length(levs))
            cols <- setNames(pal, levs)[vals]
          }

          p_base <- plot_ly(
            data = dd,
            x = ~x,
            y = ~y,
            type = "scatter",
            mode = "markers",
            marker = list(color = cols, size = 3),
            source = ns("embed"),
            customdata = ~.cell
          ) %>% layout(
            xaxis = list(title = paste0(embedding_name, " 1")),
            yaxis = list(title = paste0(embedding_name, " 2")),
            dragmode = "zoom"
          )
        }

        plot_cache_base(p_base)
        plot_cache(p_base)
      }
    )

    # Update overlays without rebuilding points
    observeEvent(input$show_labels, {
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

      plot_cache(p)
    })

    output$embed_plot <- renderPlotly({
      req(plot_cache())
      plot_cache()
    })

    output$export_embed_pdf <- downloadHandler(
      filename = function() {
        embed_name <- tolower(embedding_name) # "umap" or "tsne"
        color_var <- tolower(input$color_by %||% "color")

        if (!is.null(input$split_by) && nzchar(input$split_by)) {
          facet_var <- tolower(input$split_by)
          paste0(embed_name, "_", facet_var, "_", color_var, ".pdf")
        } else {
          paste0(embed_name, "_", color_var, ".pdf")
        }
      },
      content = function(file) {
        gg <- plot_cache_gg()
        if (is.null(gg) || !inherits(gg, "ggplot")) {
          showNotification("No ggplot available to export. Try re‑plotting before exporting.", type = "error")
          return()
        }

        if (!is.null(input$split_by) && nzchar(input$split_by)) {
          # Faceted → dynamic sizing
          n_facets <- length(unique(df()[[input$split_by]]))
          if (n_facets < 1) n_facets <- 1

          ncol_facets <- suppressWarnings(as.numeric(input$max_facets))
          if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1

          nrow_facets <- ceiling(n_facets / ncol_facets)
          if (!is.finite(nrow_facets) || nrow_facets < 1) nrow_facets <- 1

          pdf_width <- 6 * ncol_facets
          pdf_height <- 6 * nrow_facets
        } else {
          # Single panel → fixed size
          pdf_width <- 8
          pdf_height <- 8
        }

        ggsave(file,
          plot = gg, device = if (capabilities("cairo")) cairo_pdf else pdf,
          width = pdf_width, height = pdf_height, units = "in"
        )
      },
      contentType = "application/pdf"
    )

    # Ensure the plot renders even when the tab is hidden (prevents some init races)
    outputOptions(output, "embed_plot", suspendWhenHidden = FALSE)

    observeEvent(input$clear_selection, current_sel(integer(0)))
  })
}

# --- Main UI ---
ui <- navbarPage(
  "FCView",
  id = "main_tab",
  header = tags$head(
    tags$style(HTML("
      /* Grey out disabled tabs */
      #main_tab li.disabled > a,
      #main_tab li.disabled > a:hover {
        color: #aaa !important;
        cursor: not-allowed;
        text-decoration: none;
      }
    ")),
    tags$script(HTML("
      (function() {
        var tabsLocked = true;
        function disableTabs() {
          tabsLocked = true;
          $('#main_tab li').addClass('disabled');
        }
        function enableTabs() {
          tabsLocked = false;
          $('#main_tab li').removeClass('disabled');
        }
        $(document).on('show.bs.tab click', '#main_tab a[data-toggle=\"tab\"]', function(e) {
          if (tabsLocked) {
            e.preventDefault();
            e.stopImmediatePropagation();
            return false;
          }
        });
        Shiny.addCustomMessageHandler('enableTabs', function(enable) {
          if (enable) {
            enableTabs();
          } else {
            disableTabs();
            var $home = $('#main_tab a[data-toggle=\"tab\"]').first();
            if ($home.length) $home.tab('show');
          }
        });
        $(document).ready(disableTabs);
        $(document).on('shiny:connected', disableTabs);
      })();
    "))
  ),
  tabPanel(
    "Home",
    sidebarLayout(
      sidebarPanel(
        fileInput("rdata_upload", "Upload .RData (FCSimple analysis object)", accept = ".RData"),
        numericInput("max_cells_upload", "Max cells to read in", value = 300000, min = 1000, step = 1000),
        helpText("If the uploaded dataset has more cells than this number, it will be randomly downsampled after upload. This is done to speed up UMAP and tSNE facet plotting."),
        width = 3
      ),
      mainPanel(
        h3("Dataset overview"),
        verbatimTextOutput("ds_summary"),
        h3("Available metadata"),
        tableOutput("meta_overview")
      )
    )
  ),
  tabPanel(
    "Global Settings",
    h3("Available Metadata Features & Type Coercion"),
    helpText("Select metadata features to update. Categorical features can be converted to continuous where possible and vice versa for downstream testing. You may also choose to hide features from the app by selecting the 'Hide' checkbox. Changes take effect immediately."),
    fluidRow(
      column(
        6,
        pickerInput("features_dropdown", "Select Feature(s) to Update", choices = NULL, multiple = TRUE),
        uiOutput("features_mini_ui")
      ),
      column(
        6,
        verbatimTextOutput("global_settings_summary")
      )
    ),
    hr(),
    h3("Paired Testing Settings"),
    helpText("Select a metadata column to enable paired statistical tests. This column should identify matching samples/patients across conditions."),
    fluidRow(
      column(
        6,
        pickerInput("pairing_var", "Pairing variable",
          choices = NULL,
          options = list(`none-selected-text` = "None (unpaired tests)")
        )
      ),
      column(
        6,
        uiOutput("pairing_summary_ui")
      )
    ),
    hr(),
    h3("Metadata Subsetting"),
    helpText("Subset samples for testing/classification/regression based on metadata values. UMAP and tSNE tabs are not affected."),
    checkboxInput("enable_subsetting", "Enable subsetting", value = FALSE),
    conditionalPanel(
      condition = "input.enable_subsetting",
      fluidRow(
        column(
          6,
          h4("Subsetting Rules"),
          uiOutput("subsetting_rules_ui"),
          br(),
          actionButton("add_subset_rule", "Add Rule", icon = icon("plus")),
          br(), br(),
          radioButtons("subset_logic", "Combine multiple rules using:",
            choices = c("Intersection (AND)" = "intersection", "Union (OR)" = "union"),
            selected = "intersection"
          ),
          br(),
          pickerInput("trim_incomplete_pairs",
            "Trim incomplete pairs by feature(s):",
            choices = NULL,
            multiple = TRUE,
            options = list(`none-selected-text` = "(no trim)")
          ),
          helpText("Optional: Remove samples with incomplete pairing based on selected categorical features. Applied after subsetting rules."),
          br(),
          actionButton("apply_subsetting", "Apply Subsetting", icon = icon("filter"), class = "btn-primary"),
          br(), br()
        ),
        column(
          6,
          h4("Preview"),
          uiOutput("subset_preview_ui")
        )
      )
    ),
    hr(),
    h3("Data Export"),
    fluidRow(
      column(
        12,
        uiOutput("export_buttons_ui"),
        br()
      )
    )
  ),
  tabPanel("UMAP", EmbeddingUI("umap", title = "UMAP")),
  tabPanel("tSNE", EmbeddingUI("tsne", title = "tSNE")),
  tabPanel(
    "Annotation",
    fluidRow(
      column(
        3,
        checkboxInput("cluster_rows", "Cluster rows", value = TRUE),
        checkboxInput("cluster_columns", "Cluster columns", value = TRUE),
        selectInput("heatmap_theme", "Heatmap color theme",
          choices = c("viridis", "heat", "greyscale"),
          selected = "viridis"
        ),
        br(),
        downloadButton("export_heatmap_pdf", "Export heatmap as PDF")
      ),
      column(
        9, 
        plotOutput("cluster_heatmap", height = "700px"),
        br(),
        hr(),
        h4("Cluster Annotation Engine"),
        helpText("Define cell types and assign clusters to them. Annotations are used throughout the app when 'Celltypes' is selected."),
        uiOutput("celltype_annotations_ui"),
        br(),
        actionButton("add_celltype", "Add Cell Type", icon = icon("plus")),
        br(), br(),
        conditionalPanel(
          condition = "output.hasAnnotations",
          actionButton("apply_annotations", "Apply Annotations", icon = icon("check"), class = "btn-success")
        )
      )
    )
  ),
  tabPanel(
    "Testing",
    h4("Test Settings"),
    fluidRow(
      column(
        3,
        conditionalPanel(
          condition = "output.hasClusterMap",
          pickerInput("test_entity", "Entity", choices = c("Clusters", "Celltypes"), selected = "Clusters")
        ),
        pickerInput("group_var", "Categorical metadata", choices = NULL, options = list(`none-selected-text` = "None")),
        pickerInput("cont_var", "Continuous metadata", choices = NULL, options = list(`none-selected-text` = "None")),
        uiOutput("test_type_ui"),
        uiOutput("paired_test_indicator_ui"),
        selectInput("p_adj_method", "P‑value adjustment method",
          choices = c("BH", "bonferroni", "BY", "fdr"), selected = "BH"
        ),
        actionButton("run_test", "Run tests"),
        br(), br(),
        conditionalPanel(
          condition = "output.hasValidResults",
          downloadButton("export_results", "Export results as CSV"),
          br(), br(),
          actionButton("reset_test", "Clear Results")
        )
      ),
      column(
        9,
        conditionalPanel(
          condition = "output.hasResults",
          h4("Test Results"),
          uiOutput("test_error_msg"),
          tableOutput("test_table")
        ),
        textOutput("test_cleared_msg")
      )
    )
  ),
  tabPanel(
    "Categorical",
    h4("Plot Settings"),
    fluidRow(
      column(
        3,
        pickerInput("cat_entity", "Entity", choices = c("Clusters", "Celltypes"), selected = "Clusters"),
        pickerInput("cat_group_var", "Categorical metadata", choices = NULL, options = list(`none-selected-text` = "None")),
        uiOutput("cat_test_type_ui"),
        uiOutput("paired_cat_indicator_ui"),
        uiOutput("cat_pairing_check_ui"),
        selectInput("cat_p_adj_method", "P‑value adjustment method",
          choices = c("BH", "bonferroni", "BY", "fdr"), selected = "BH"
        ),
        checkboxInput("cat_use_adj_p", "Plot adjusted pvalues", value = TRUE),
        selectInput("cat_max_facets", "Facet columns", choices = 2:6, selected = 4),
        selectInput("cat_plot_type", "Plot type", choices = c("Boxplot" = "box", "Violin" = "violin"), selected = "box"),
        uiOutput("cat_points_ui"),
        br(),
        actionButton("cat_populate_colors", "Populate colors for selected group variable"),
        br(),
        uiOutput("cat_color_pickers_ui"),
        br(),
        actionButton("generate_cat_plots", "Generate plots"),
        br(), br(),
        conditionalPanel(
          condition = "output.hasValidCatResults",
          downloadButton("export_cat_pdf", "Export boxplots as PDF"),
          br(), br(),
          actionButton("reset_cat", "Clear Results")
        )
      ),
      column(
        9,
        conditionalPanel(
          condition = "output.hasCatResults",
          h4("Categorical Plots"),
          uiOutput("cat_error_msg"),
          plotOutput("categorical_plot")
        ),
        textOutput("cat_cleared_msg")
      )
    )
  ),
  tabPanel(
    "Continuous",
    h4("Plot Settings"),
    fluidRow(
      column(
        3,
        pickerInput("cont_entity", "Entity", choices = c("Clusters", "Celltypes"), selected = "Clusters"),
        pickerInput("cont_group_var", "Continuous metadata", choices = NULL, options = list(`none-selected-text` = "None")),
        selectInput("cont_p_adj_method", "P‑value adjustment method",
          choices = c("BH", "bonferroni", "BY", "fdr"), selected = "BH"
        ),
        checkboxInput("cont_use_adj_p", "Plot adjusted pvalues", value = TRUE),
        checkboxInput("cont_transpose", "Transpose axes", value = FALSE),
        selectInput("cont_max_facets", "Facet columns", choices = 2:6, selected = 4),
        actionButton("generate_cont_plots", "Generate plots"),
        br(), br(),
        conditionalPanel(
          condition = "output.hasContResults",
          downloadButton("export_cont_pdf", "Export scatter plots as PDF"),
          br(), br(),
          actionButton("reset_cont", "Clear Results")
        )
      ),
      column(
        9,
        conditionalPanel(
          condition = "output.hasContResults",
          h4("Continuous Plots"),
          plotOutput("continuous_plot")
        ),
        textOutput("cont_cleared_msg")
      )
    )
  ),
  tabPanel(
    "Feature Selection",
    h4("Model Settings"),
    fluidRow(
      column(
        2,
        pickerInput("fs_method", "Method",
          choices = c("Ridge Regression", "Elastic Net", "Random Forest (Boruta)")
        ),
        pickerInput("fs_outcome", "Outcome variable", choices = NULL),
        pickerInput("fs_predictors", "Predictor(s)", choices = NULL, multiple = TRUE),
        conditionalPanel(
          condition = "input.fs_predictors.includes('cluster')",
          pickerInput("fs_cluster_subset", "Select clusters to include",
            choices = NULL, multiple = TRUE
          )
        ),
        conditionalPanel(
          condition = "input.fs_method == 'Elastic Net'",
          sliderInput("fs_alpha", "Alpha (0 = Ridge, 1 = Lasso)",
            min = 0, max = 1, value = 0.5, step = 0.05
          )
        ),
        actionButton("run_fs", "Run Feature Selection"),
        br(), br(),
        conditionalPanel(
          condition = "output.hasFSResults",
          downloadButton("export_fs_results", "Download All Results (ZIP)"),
          br(), br(),
          actionButton("reset_fs", "Clear Results")
        )
      ),
      column(
        7,
        conditionalPanel(
          condition = "output.hasFSResults",
          h4("Summary Plot"),
          uiOutput("fs_plot_ui"),
          h4("Details"),
          verbatimTextOutput("fs_summary")
        ),
        textOutput("fs_cleared_msg")
      ),
      column(
        3,
        conditionalPanel(
          condition = "output.hasFSResults",
          h4("Selected Features"),
          tableOutput("fs_results")
        )
      )
    )
  ),
  tabPanel(
    "Classification",
    h4("Model Settings"),
    fluidRow(
      column(
        2,
        pickerInput("lm_outcome", "Outcome variable", choices = NULL),
        pickerInput("lm_predictors", "Predictor(s)", choices = NULL, multiple = TRUE),
        conditionalPanel(
          condition = "input.lm_predictors.includes('cluster')",
          pickerInput("lm_cluster_subset", "Select clusters to include", choices = NULL, multiple = TRUE)
        ),
        radioButtons("lm_model_type", "Model type",
          choices = c("Logistic Regression", "Elastic Net", "Random Forest")
        ),
        conditionalPanel(
          condition = "input.lm_model_type == 'Elastic Net'",
          sliderInput("lm_alpha", "Elastic Net alpha (0 = Ridge, 1 = Lasso)",
            min = 0, max = 1, value = 0.5, step = 0.05
          )
        ),
        radioButtons("lm_validation", "Validation strategy",
          choices = c("Train/Test split", "k-fold CV", "Leave-One-Out")
        ),
        conditionalPanel(
          condition = "input.lm_validation == 'Train/Test split'",
          sliderInput("lm_train_frac", "Train fraction",
            min = 0.5, max = 0.95,
            value = 0.7, step = 0.05
          )
        ),
        conditionalPanel(
          condition = "input.lm_validation == 'k-fold CV'",
          numericInput("lm_k", "Number of folds", value = 5, min = 2, max = 20)
        ),
        actionButton("run_lm", "Run Model"),
        br(), br(),
        conditionalPanel(
          condition = "output.hasLMResults",
          downloadButton("export_lm_zip", "Download All Results (ZIP)"),
          br(), br(),
          actionButton("reset_lm", "Clear Results")
        )
      ),
      column(
        7,
        conditionalPanel(
          condition = "output.hasLMResults",
          h4("ROC Curve"),
          plotOutput("lm_roc_plot", height = "500px"),
          fluidRow(
            column(
              4,
              h4("Model Summary"),
              verbatimTextOutput("lm_summary")
            ),
            column(
              8,
              h4("Performance Metrics"),
              tableOutput("lm_perf_table")
            )
          )
        ),
        textOutput("lm_cleared_msg")
      ),
      column(
        3,
        conditionalPanel(
          condition = "output.hasLMResults",
          h4("Model Features"),
          tableOutput("lm_features")
        )
      )
    )
  ),
  tabPanel(
    "Regression",
    h4("Model Settings"),
    fluidRow(
      column(
        3,
        pickerInput("reg_outcome", "Outcome variable (continuous)", choices = NULL),
        pickerInput("reg_predictors", "Predictor(s)", choices = NULL, multiple = TRUE),
        conditionalPanel(
          condition = "input.reg_predictors.includes('cluster')",
          pickerInput("reg_cluster_subset", "Select clusters to include", choices = NULL, multiple = TRUE)
        ),
        radioButtons("reg_model_type", "Model type",
          choices = c("Linear Regression", "Ridge Regression", "Elastic Net", "Random Forest")
        ),
        conditionalPanel(
          condition = "input.reg_model_type == 'Elastic Net'",
          sliderInput("reg_alpha", "Elastic Net alpha (0 = Ridge, 1 = Lasso)",
            min = 0, max = 1, value = 0.5, step = 0.05
          )
        ),
        radioButtons("reg_validation", "Validation strategy",
          choices = c("Train/Test split", "k-fold CV", "Leave-One-Out")
        ),
        conditionalPanel(
          condition = "input.reg_validation == 'Train/Test split'",
          sliderInput("reg_train_frac", "Train fraction",
            min = 0.5, max = 0.95,
            value = 0.7, step = 0.05
          )
        ),
        conditionalPanel(
          condition = "input.reg_validation == 'k-fold CV'",
          numericInput("reg_k", "Number of folds", value = 5, min = 2, max = 20)
        ),
        actionButton("run_reg", "Run Model"),
        br(), br(),
        conditionalPanel(
          condition = "output.hasRegResults",
          downloadButton("export_reg_zip", "Download All Results (ZIP)"),
          br(), br(),
          actionButton("reset_reg", "Clear Results")
        )
      ),
      column(
        3,
        conditionalPanel(
          condition = "output.hasRegResults",
          h4("Model Summary"),
          verbatimTextOutput("reg_summary"),
          h4("Performance Metrics"),
          tableOutput("reg_perf_table")
        ),
        textOutput("reg_cleared_msg")
      ),
      column(
        3,
        conditionalPanel(
          condition = "output.hasRegResults",
          h4("Observed vs Predicted"),
          plotOutput("reg_obs_pred_plot", height = "400px"),
          h4("Residuals vs Fitted"),
          plotOutput("reg_residual_plot", height = "400px")
        )
      ),
      column(
        3,
        conditionalPanel(
          condition = "output.hasRegResults",
          h4("Model Features"),
          tableOutput("reg_features")
        )
      )
    )
  ),
  tabPanel(
    "sccomp",
    h4("Differential Composition Analysis with sccomp"),
    fluidRow(
      column(
        3,
        helpText("sccomp tests for differential composition of cell types/clusters across sample groups."),
        radioButtons("sccomp_formula_mode", "Formula mode",
          choices = c("Simple" = "simple", "Custom" = "custom"),
          selected = "simple"
        ),
        conditionalPanel(
          condition = "input.sccomp_formula_mode == 'simple'",
          pickerInput("sccomp_group_var", "Grouping variable", choices = NULL),
          conditionalPanel(
            condition = "input.sccomp_group_var != '' && input.sccomp_group_var != null",
            pickerInput("sccomp_reference_level", "Reference level (first level is default)",
              choices = NULL, options = list(`live-search` = TRUE)
            )
          ),
          pickerInput("sccomp_formula_vars", "Additional covariates (optional)", choices = NULL, multiple = TRUE),
          checkboxInput("sccomp_interactions", "Include interactions", FALSE)
        ),
        conditionalPanel(
          condition = "input.sccomp_formula_mode == 'custom'",
          helpText("Formula syntax examples:"),
          helpText("Fixed effects: ~ condition + age"),
          helpText("Interaction: ~ condition * age"),
          helpText("Random effects: ~ condition + (1|patient_ID)"),
          helpText("No intercept: ~ 0 + condition"),
          textInput("sccomp_custom_formula", "Custom formula", value = "~ condition", placeholder = "~ condition + age"),
          textInput("sccomp_custom_reference_levels", "Reference levels (optional)",
            value = "",
            placeholder = "variable1=level1; variable2=level2"
          ),
          helpText("Specify reference levels as: variable=level; separated by semicolons. Example: condition=Control; treatment=Placebo")
        ),
        sliderInput("sccomp_cores", "Number of cores",
          min = 1, max = parallel::detectCores(),
          value = max(1, parallel::detectCores() - 2), step = 1
        ),
        actionButton("run_sccomp", "Run sccomp_estimate", class = "btn-primary"),
        br(), br(),
        conditionalPanel(
          condition = "output.hasSccompResults",
          h5("Post-hoc Contrasts"),
          helpText("After running sccomp_estimate, you can test specific contrasts below."),
          helpText(strong("Available parameters for contrasts:")),
          helpText("Use backticks (`) around parameter names. See options below."),
          verbatimTextOutput("sccomp_available_params"),
          textInput("sccomp_contrast", "Contrast specification",
            placeholder = "e.g., `conditionTreated` - `conditionControl`"
          ),
          actionButton("run_sccomp_test", "Run sccomp_test", class = "btn-secondary"),
          br(), br(),
          downloadButton("export_sccomp_results", "Export results as CSV"),
          br(), br(),
          actionButton("reset_sccomp", "Clear Results")
        )
      ),
      column(
        9,
        conditionalPanel(
          condition = "output.hasSccompResults",
          h4("Results Summary"),
          verbatimTextOutput("sccomp_summary"),
          uiOutput("sccomp_intercept_warning"),
          fluidRow(
            column(
              6,
              h4("Credible Intervals Plot"),
              plotOutput("sccomp_interval_plot", height = "auto", width = "100%"),
              downloadButton("download_sccomp_plot", "Download Interval Plot"),
              br(), br()
            ),
            column(
              6,
              conditionalPanel(
                condition = "output.hasSccompTestResults",
                h4("Contrast Credible Intervals Plot"),
                plotOutput("sccomp_contrast_plot", height = "auto", width = "100%"),
                downloadButton("download_sccomp_contrast_plot", "Download Contrast Interval Plot"),
                br(), br()
              )
            )
          ),
          fluidRow(
            column(
              6,
              uiOutput("sccomp_table_ui")
            ),
            column(
              6,
              conditionalPanel(
                condition = "output.hasSccompTestResults",
                h4("Contrast Test Results"),
                uiOutput("sccomp_test_table_ui")
              )
            )
          )
        ),
        textOutput("sccomp_cleared_msg")
      )
    )
  ),
  tabPanel(
    "Time to Event",
    icon = icon("chart-line"),
    h4("Model Settings"),
    fluidRow(
      column(
        3,
        pickerInput("surv_outcome", "Time-to-event variable (continuous)", choices = NULL),
        hr(),
        selectInput("surv_analysis_mode", "Mode:",
          choices = c("Multivariate" = "multivariate", "Univariate" = "univariate"),
          selected = "multivariate"
        ),
        helpText("Multivariate: all predictors in one model."),
        helpText("Univariate: test each predictor separately."),
        hr(),
        pickerInput("surv_predictors", "Predictor(s)", choices = NULL, multiple = TRUE),
        conditionalPanel(
          condition = "input.surv_predictors.includes('cluster')",
          pickerInput("surv_cluster_subset", "Select clusters to include", choices = NULL, multiple = TRUE)
        ),
        hr(),
        selectInput("surv_split_method", "Split risk score group by:",
          choices = c("Median" = "median", "Mean" = "mean"),
          selected = "median"
        ),
        checkboxInput("surv_show_ci", "Show confidence intervals", value = TRUE),
        helpText("Use for visualization only."),
        br(),
        actionButton("run_surv", "Run Model", class = "btn-primary"),
        br(), br(),
        conditionalPanel(
          condition = "output.hasSurvResults && input.surv_analysis_mode == 'univariate'",
          selectInput("surv_univar_predictor_display", "Display results for:",
            choices = NULL
          ),
          helpText("Select which predictor's univariate results to display.")
        ),
        conditionalPanel(
          condition = "output.hasSurvResults",
          downloadButton("export_surv_zip", "Download All Results (ZIP)"),
          br(), br(),
          actionButton("reset_surv", "Clear Results")
        )
      ),
      column(
        5,
        conditionalPanel(
          condition = "output.hasSurvResults",
          h4("Time to Event Curve"),
          plotOutput("surv_curve_plot", height = "500px"),
          fluidRow(
            column(
              4,
              h4("Model Summary"),
              verbatimTextOutput("surv_summary")
            ),
            column(
              8,
              h4("Performance Metrics"),
              tableOutput("surv_perf_table")
            )
          )
        ),
        uiOutput("surv_error_ui"),
        textOutput("surv_cleared_msg")
      ),
      column(
        4,
        conditionalPanel(
          condition = "output.hasSurvCoefs",
          uiOutput("surv_coef_header"),
          tableOutput("surv_coefs")
        )
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
    rep_used = NULL,
    data_ready = FALSE,
    # Global settings for available features
    available_features = character(0),
    all_meta_cols = character(0),
    # Type coercion
    meta_cached = NULL, # Cached original metadata with original types
    type_coercions = list(), # Named list: col_name -> coercion_type
    # Mini UI persisted states: hide checkbox per feature
    mini_hide_states = list(),
    # Subsetting
    meta_sample_original = NULL, # Original unfiltered metadata
    subsetting_enabled = FALSE,
    subset_rules = list(), # List of subsetting rules
    subset_summary = NULL, # Summary of subsetting results
    subset_id = "000000000", # Unique ID for current subset - "000000000" means no subsetting applied
    # Cluster annotation engine
    celltype_annotations = list() # List of celltype -> cluster assignments
  )
  cat_plot_cache <- reactiveVal(NULL)
  cont_plot_cache <- reactiveVal(NULL)
  test_results_rv <- reactiveVal(NULL)
  cat_state <- reactiveVal(NULL)
  cont_state <- reactiveVal(NULL)
  fs_state <- reactiveVal(NULL)
  lm_state <- reactiveVal(NULL)
  reg_state <- reactiveVal(NULL)
  surv_state <- reactiveVal(NULL)
  celltype_id_counter <- reactiveVal(0) # Counter for unique celltype IDs
  # rv$log <- reactiveVal(character())

  # Disable tabs at startup
  # observe({
  #   if (!isTRUE(rv$data_ready)) {
  #     session$sendCustomMessage("enableTabs", FALSE)
  #   }
  # })

  observeEvent(input$main_tab, {
    message("Tab changed to: ", input$main_tab)
  })

  # Upload RData and initialize datasets immediately (no mapping button)
  observeEvent(input$rdata_upload, {
    session$sendCustomMessage("enableTabs", FALSE) # lock tabs during load
    req(input$rdata_upload)

    # Load all objects from the uploaded .RData
    e <- new.env(parent = emptyenv())
    load(input$rdata_upload$datapath, envir = e)

    # Find all objects with the required structure
    candidates <- ls(e)
    matches <- vapply(candidates, function(nm) {
      cand <- e[[nm]]
      is.list(cand) && all(c("data", "source", "metadata", "cluster") %in% names(cand))
    }, logical(1))

    if (sum(matches) == 1) {
      obj <- e[[candidates[matches][1]]]
    } else if (sum(matches) > 1) {
      showNotification(
        paste(
          "Multiple valid objects found in .RData:",
          paste(candidates[matches], collapse = ", ")
        ),
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

      obj$data <- obj$data[keep_idx, , drop = FALSE]
      obj$source <- obj$source[keep_idx]

      if (!is.null(obj$cluster$clusters) && length(obj$cluster$clusters) == n_cells) {
        obj$cluster$clusters <- obj$cluster$clusters[keep_idx]
      }
      if (!is.null(obj$umap$coordinates) && nrow(obj$umap$coordinates) == n_cells) {
        obj$umap$coordinates <- obj$umap$coordinates[keep_idx, , drop = FALSE]
      }
      if (!is.null(obj$tsne$coordinates) && nrow(obj$tsne$coordinates) == n_cells) {
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

    # --- Build meta_cell (per-cell) with robust patient_ID mapping for UMAP/tSNE tabs ---
    meta_cell <- data.frame(
      source = as.character(obj$source),
      run_date = if (!is.null(run_date)) run_date else NA
    )

    # Prepare metadata IDs (character, unique, deduplicated rows)
    if (!("patient_ID" %in% colnames(obj$metadata))) {
      showNotification("metadata does not contain 'patient_ID' column.", type = "error")
      return()
    }
    obj$metadata$patient_ID <- as.character(obj$metadata$patient_ID)

    metadata_unique <- obj$metadata %>%
      dplyr::distinct(patient_ID, .keep_all = TRUE)

    # Remove run_date from metadata if it exists to avoid .x/.y suffixes during join
    if ("run_date" %in% colnames(metadata_unique)) {
      metadata_unique <- metadata_unique %>% dplyr::select(-run_date)
    }

    ids <- unique(metadata_unique$patient_ID)
    if (length(ids) == 0) {
      showNotification("No patient_ID values found in metadata.", type = "error")
      return()
    }

    # Escape regex metacharacters in IDs, then sort by length (longest first)
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
    meta_cell[] <- lapply(meta_cell, function(col) {
      if (is.character(col) && all(grepl("^\\s*-?\\d*(\\.\\d+)?\\s*$", col[!is.na(col)]))) {
        as.numeric(col)
      } else {
        col
      }
    })

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

    # Convert run_date to factor if present
    if ("run_date" %in% names(meta_cell)) {
      meta_cell$run_date <- as.factor(meta_cell$run_date)
    }

    # Add cluster factor if available
    if (!"cluster" %in% names(meta_cell) &&
      !is.null(obj$cluster$clusters) &&
      length(obj$cluster$clusters) == nrow(obj$data)) {
      meta_cell$cluster <- factor(obj$cluster$clusters,
        levels = sort(unique(obj$cluster$clusters))
      )
    }

    clusters <- list(
      assignments = obj$cluster$clusters,
      settings    = obj$cluster$settings %||% list()
    )
    cluster_map <- if (!is.null(obj$cluster_mapping)) obj$cluster_mapping else NULL

    UMAP <- if (!is.null(obj$umap)) {
      list(
        coords = obj$umap$coordinates,
        settings = obj$umap$settings
      )
    } else {
      NULL
    }
    tSNE <- if (!is.null(obj$tsne)) {
      list(
        coords = obj$tsne$coordinates,
        settings = obj$tsne$settings
      )
    } else {
      NULL
    }

    cluster_heat <- if (!is.null(obj$cluster_heatmap)) obj$cluster_heatmap$heatmap_tile_data else NULL
    pop_size <- if (!is.null(obj$cluster_heatmap)) obj$cluster_heatmap$population_size else NULL
    rep_used <- if (!is.null(obj$cluster_heatmap)) obj$cluster_heatmap$rep_used else NA

    # --- Add per-sample abundance matrix and canonical per-sample metadata for downstream tabs ---
    if (!is.null(obj$cluster$abundance) && is.matrix(obj$cluster$abundance)) {
      clusters$abundance <- obj$cluster$abundance
      rv$abundance_sample <- clusters$abundance

      if (is.null(rownames(rv$abundance_sample)) || any(!nzchar(rownames(rv$abundance_sample)))) {
        showNotification(
          "cluster$abundance has missing rownames; cannot map to metadata. Please set rownames to source strings.",
          type = "error",
          duration = NULL
        )
        message("Abundance matrix rownames are missing or empty; mapping to metadata will fail.")
      } else {
        message(sprintf(
          "Abundance matrix loaded: %d sources × %d entities",
          nrow(rv$abundance_sample), ncol(rv$abundance_sample)
        ))
      }
    } else {
      clusters$abundance <- NULL
      rv$abundance_sample <- NULL
      showNotification("No cluster$abundance matrix found in upload; abundance-based tabs will be disabled.", type = "warning")
    }

    # Load counts data for sccomp
    if (!is.null(obj$cluster$counts) && is.matrix(obj$cluster$counts)) {
      rv$counts_sample <- obj$cluster$counts
      message(sprintf(
        "Counts matrix loaded: %d sources × %d entities",
        nrow(rv$counts_sample), ncol(rv$counts_sample)
      ))
    } else {
      rv$counts_sample <- NULL
      message("No cluster$counts matrix found in upload; sccomp tab will be disabled.")
    }

    # Canonical per-sample metadata (prefer direct sample-level object if provided)
    # If obj already includes a per-sample metadata frame, use it; otherwise derive from metadata_unique.
    if (!is.null(obj$metadata_sample) && is.data.frame(obj$metadata_sample) &&
      "patient_ID" %in% colnames(obj$metadata_sample)) {
      rv$meta_sample <- obj$metadata_sample %>%
        dplyr::distinct(patient_ID, .keep_all = TRUE)
    } else {
      rv$meta_sample <- metadata_unique # one row per patient_ID
    }

    # Store in rv (cell-level objects preserved for UMAP/tSNE/Heatmap)
    rv$expr <- expr
    rv$meta_cell <- meta_cell
    rv$clusters <- clusters
    rv$cluster_map <- cluster_map
    rv$UMAP <- UMAP
    rv$tSNE <- tSNE
    rv$cluster_heat <- cluster_heat
    rv$pop_size <- pop_size
    rv$rep_used <- rep_used
    
    # Initialize celltype annotations from cluster_map if available
    if (!is.null(cluster_map) && all(c("cluster", "celltype") %in% names(cluster_map))) {
      # Extract unique celltypes and their cluster mappings
      annotations <- list()
      celltype_names <- unique(cluster_map$celltype)
      for (i in seq_along(celltype_names)) {
        ct <- celltype_names[i]
        assigned_clusters <- cluster_map$cluster[cluster_map$celltype == ct]
        annotations[[as.character(i)]] <- list(
          id = i,
          name = ct,
          clusters = as.character(assigned_clusters)
        )
      }
      rv$celltype_annotations <- annotations
      celltype_id_counter(length(celltype_names))
      message("Initialized celltype annotations from cluster_map: ", length(celltype_names), " cell types")
    } else {
      rv$celltype_annotations <- list()
      celltype_id_counter(0)
      message("No cluster_map found; celltype annotations initialized empty")
    }

    message(
      "Upload complete: expr rows=", nrow(rv$expr),
      " meta_cell rows=", nrow(rv$meta_cell),
      " meta_sample rows=", nrow(rv$meta_sample),
      " abundance_sample rows=", if (!is.null(rv$abundance_sample)) nrow(rv$abundance_sample) else 0,
      " UMAP coords=", if (!is.null(rv$UMAP)) nrow(rv$UMAP$coords) else 0,
      " tSNE coords=", if (!is.null(rv$tSNE)) nrow(rv$tSNE$coords) else 0
    )

    showNotification("Data loaded and initialized.", type = "message")

    # Initialize FS/LM pickers from per-sample objects
    if (!is.null(rv$meta_sample)) {
      meta_cols <- colnames(rv$meta_sample)

      # For Feature Selection: allow both categorical and continuous outcomes
      fs_outcome_choices <- sort(meta_cols)
      fs_outcome_choices <- setdiff(fs_outcome_choices, c("cluster", "patient_ID", "run_date", "source"))

      # For Classification (LM): only categorical outcomes
      categorical_choices <- sort(meta_cols[sapply(rv$meta_sample, function(x) is.factor(x) || is.character(x))])
      # Exclude cluster, patient_ID, run_date, and source from outcome choices
      categorical_choices <- setdiff(categorical_choices, c("cluster", "patient_ID", "run_date", "source"))

      predictor_choices <- sort(meta_cols)
      # Exclude patient_ID, source, and run_date from predictor choices
      predictor_choices <- setdiff(predictor_choices, c("patient_ID", "source", "run_date"))

      updatePickerInput(session, "fs_outcome", choices = fs_outcome_choices, selected = NULL)
      updatePickerInput(session, "fs_predictors", choices = c(predictor_choices, "cluster"), selected = NULL)
      updatePickerInput(session, "fs_cluster_subset",
        choices = if (!is.null(rv$abundance_sample)) colnames(rv$abundance_sample) else character(0),
        selected = character(0)
      )

      # categorical_choices already excludes cluster and patient_ID
      updatePickerInput(session, "lm_outcome", choices = categorical_choices, selected = NULL)
      updatePickerInput(session, "lm_predictors", choices = c(predictor_choices, "cluster"), selected = NULL)
      updatePickerInput(session, "lm_cluster_subset",
        choices = if (!is.null(rv$abundance_sample)) colnames(rv$abundance_sample) else character(0),
        selected = character(0)
      )
    }

    rv$data_ready <- TRUE

    # Initialize global settings with all metadata columns
    if (!is.null(rv$meta_sample)) {
      all_cols <- setdiff(colnames(rv$meta_sample), c("patient_ID", "source", "run_date"))
      rv$all_meta_cols <- all_cols
      rv$available_features <- all_cols # All selected by default
      # Store original metadata for subsetting
      rv$meta_sample_original <- rv$meta_sample
      # Cache original metadata with original types for type coercion
      rv$meta_cached <- rv$meta_sample
      # Initialize type coercions as empty (no coercions applied)
      rv$type_coercions <- setNames(rep(list(""), length(all_cols)), all_cols)
    }

    session$sendCustomMessage("enableTabs", TRUE)
  })

  # ========== TYPE COERCION HELPERS ==========

  # Function to validate which coercion paths are valid for a given column
  validate_coercions <- function(col_data) {
    valid_types <- character(0)

    # Always allow character
    valid_types <- c(valid_types, "character")

    # Always allow factor
    valid_types <- c(valid_types, "factor")

    # Test numeric by level (factor -> numeric using levels)
    if (is.factor(col_data)) {
      tryCatch(
        {
          test <- as.numeric(col_data)
          if (!all(is.na(test))) valid_types <- c(valid_types, "numeric by level")
        },
        error = function(e) {},
        warning = function(w) {}
      )
    }

    # Test numeric by character (convert to char first, then numeric)
    tryCatch(
      {
        test <- suppressWarnings(as.numeric(as.character(col_data)))
        if (!all(is.na(test))) valid_types <- c(valid_types, "numeric by character")
      },
      error = function(e) {},
      warning = function(w) {}
    )

    # Test integer by level (factor -> integer using levels)
    if (is.factor(col_data)) {
      tryCatch(
        {
          test <- as.integer(col_data)
          if (!all(is.na(test))) valid_types <- c(valid_types, "integer by level")
        },
        error = function(e) {},
        warning = function(w) {}
      )
    }

    # Test integer by character (convert to char first, then integer)
    tryCatch(
      {
        test <- suppressWarnings(as.integer(as.character(col_data)))
        if (!all(is.na(test))) valid_types <- c(valid_types, "integer by character")
      },
      error = function(e) {},
      warning = function(w) {}
    )

    return(valid_types)
  }

  # Function to apply coercion based on type string
  apply_coercion <- function(col_data, coercion_type) {
    if (coercion_type == "" || is.null(coercion_type)) {
      return(col_data)
    }

    result <- switch(coercion_type,
      "character" = as.character(col_data),
      "factor" = as.factor(col_data),
      "numeric by level" = as.numeric(col_data), # Works for factors
      "numeric by character" = suppressWarnings(as.numeric(as.character(col_data))),
      "integer by level" = as.integer(col_data), # Works for factors
      "integer by character" = suppressWarnings(as.integer(as.character(col_data))),
      col_data # Default: no change
    )

    return(result)
  }

  # ========== GLOBAL SETTINGS TAB ==========

  # Render feature edit section as multi-selection dropdown
  output$features_dropdown <- renderUI({
    # req(rv$all_meta_cols, rv$meta_cached)
    pickerInput(
      "features_dropdown",
      "Select Features",
      # choices = c("no selection" = '', rv$all_meta_cols),
      choices = rv$all_meta_cols,
      selected = NULL,
      multiple = TRUE
    )
  })

  # Render per-selected-feature mini UI: hide checkbox + coercion dropdown
  output$features_mini_ui <- renderUI({
    req(rv$all_meta_cols)
    sel <- input$features_dropdown %||% character(0)
    if (length(sel) == 0) {
      return(NULL)
    }

    mini_uis <- lapply(sel, function(col) {
      safe_id <- gsub("[^A-Za-z0-9]", "_", col)
      # Coercion choices derived from cached metadata
      valid_types <- validate_coercions(rv$meta_cached[[col]])
      type_choices <- c("(no coercion)" = "", valid_types)
      current_type <- rv$type_coercions[[col]] %||% ""
      # persisted hide state
      hide_state <- rv$mini_hide_states[[col]] %||% FALSE

      fluidRow(
        column(4, tags$strong(col)),
        column(2, checkboxInput(paste0("hide_mini_", safe_id), "Hide", value = hide_state)),
        column(6, selectInput(paste0("type_mini_", safe_id), NULL, choices = type_choices, selected = current_type, width = "100%"))
      )
    })

    tagList(tags$h5("Feature quick actions (hide / change data type)"), mini_uis)
  })

  # `features_checkboxes` removed — mini UI handles per-feature controls

  # Update available features when checkboxes change
  observe({
    req(rv$all_meta_cols)
    cols <- rv$all_meta_cols

    # Derive selected features from mini hide states (mini UI is the source of truth).
    # If a feature is mini-hidden, it is excluded; otherwise included.
    selected_features <- cols[!vapply(cols, function(col) isTRUE(rv$mini_hide_states[[col]]), logical(1))]
    rv$available_features <- selected_features
  })

  # Wire mini UI inputs directly into reactive state (mini UI is now the single source of truth)
  observeEvent(rv$all_meta_cols, {
    req(rv$all_meta_cols)
    cols <- rv$all_meta_cols
    lapply(cols, function(col) {
      local({
        mycol <- col
        safe_id <- gsub("[^A-Za-z0-9]", "_", mycol)
        mini_hide_id <- paste0("hide_mini_", safe_id)
        mini_type_id <- paste0("type_mini_", safe_id)

        # Persist mini hide checkbox changes
        observeEvent(input[[mini_hide_id]], {
          val <- input[[mini_hide_id]]
          if (is.null(val)) return()
          rv$mini_hide_states[[mycol]] <- isTRUE(val)
          # Reset pairing and subsetting to avoid unexpected behavior when visibility changes
          tryCatch({ updatePickerInput(session, "pairing_var", selected = "") }, error = function(e) {})
          rv$subsetting_enabled <- FALSE
          rv$subset_rules <- list()
          rv$subset_summary <- NULL
          rv$subset_id <- "000000000"
          tryCatch({ subset_rule_ids(integer(0)) }, error = function(e) {})
          tryCatch({ subset_next_id(1) }, error = function(e) {})
          tryCatch({ updateCheckboxInput(session, "enable_subsetting", value = FALSE) }, error = function(e) {})
          showNotification("Pairing variable and subsetting reset due to metadata visibility change.", type = "warning", duration = 4)
        }, ignoreInit = TRUE)

        # Persist mini type selector changes
        observeEvent(input[[mini_type_id]], {
          val <- input[[mini_type_id]]
          if (is.null(val)) return()
          rv$type_coercions[[mycol]] <- val
          # Reset pairing and subsetting to avoid unexpected behavior when types change
          tryCatch({ updatePickerInput(session, "pairing_var", selected = "") }, error = function(e) {})
          rv$subsetting_enabled <- FALSE
          rv$subset_rules <- list()
          rv$subset_summary <- NULL
          rv$subset_id <- "000000000"
          tryCatch({ subset_rule_ids(integer(0)) }, error = function(e) {})
          tryCatch({ subset_next_id(1) }, error = function(e) {})
          tryCatch({ updateCheckboxInput(session, "enable_subsetting", value = FALSE) }, error = function(e) {})
          showNotification("Pairing variable and subsetting reset due to metadata type change.", type = "warning", duration = 4)
        }, ignoreInit = TRUE)
      })
    })
  }, ignoreInit = TRUE)

  # Track if any type coercion has changed (to trigger resets)
  type_coercion_changed <- reactiveVal(FALSE)

  # Update type coercions when dropdowns change
  observe({
    req(rv$all_meta_cols, rv$meta_cached)
    cols <- rv$all_meta_cols

    coercion_changed <- FALSE
    for (col in cols) {
      input_id <- paste0("type_mini_", gsub("[^A-Za-z0-9]", "_", col))
      coercion_type <- input[[input_id]]

      if (!is.null(coercion_type)) {
        # Check if this is a real change (not initialization)
        old_coercion <- rv$type_coercions[[col]]
        if (!is.null(old_coercion) && old_coercion != coercion_type) {
          coercion_changed <- TRUE
        }
        rv$type_coercions[[col]] <- coercion_type
      }
    }

    # Signal that a coercion changed (only if data is ready and it's a real change)
    if (coercion_changed && isTRUE(rv$data_ready)) {
      type_coercion_changed(TRUE)
    }
  })

  # Legacy global select/reset controls removed; use the per-feature mini UI instead.

  # Apply type coercions to metadata reactively
  observe({
    req(rv$meta_cached, rv$type_coercions, rv$all_meta_cols)

    # Start with cached original metadata
    meta_coerced <- rv$meta_cached

    # Apply coercions to each column
    for (col in rv$all_meta_cols) {
      # Never coerce or hide the canonical patient ID
      if (identical(col, "patient_ID")) next
      if (col %in% names(meta_coerced)) {
        coercion_type <- rv$type_coercions[[col]]
        if (!is.null(coercion_type) && nzchar(coercion_type)) {
          # Skip applying coercion for features that are hidden (either via main checkbox
          # not being selected or via the mini hide control). This avoids unexpected
          # type changes for hidden columns and prevents reactive loops.
          is_visible <- col %in% rv$available_features
          is_mini_hidden <- isTRUE(rv$mini_hide_states[[col]])
          if (is_visible && !is_mini_hidden) {
            meta_coerced[[col]] <- apply_coercion(rv$meta_cached[[col]], coercion_type)
          }
        }
      }
    }

    # Update meta_sample_original with coerced types (before subsetting)
    rv$meta_sample_original <- meta_coerced

    # Always update meta_sample with coerced data
    rv$meta_sample <- meta_coerced

    # If a coercion changed, reset pairing and subsetting for safety
    if (isTRUE(type_coercion_changed())) {
      # Reset pairing variable
      updatePickerInput(session, "pairing_var", selected = "")

      # Reset subsetting
      rv$subsetting_enabled <- FALSE
      rv$subset_rules <- list()
      rv$subset_summary <- NULL
      rv$subset_id <- "000000000"
      subset_rule_ids(integer(0))
      subset_next_id(1)
      updateCheckboxInput(session, "enable_subsetting", value = FALSE)

      # Show notification
      showNotification(
        "Type coercion applied. Pairing variable and subsetting have been reset for safety.",
        type = "warning",
        duration = 6
      )

      # Reset the flag
      type_coercion_changed(FALSE)
    }
  })

  # Summary output
  output$global_settings_summary <- renderPrint({
    # React to selected features and mini UI state
    sel <- input$features_dropdown %||% character(0)
    cat("Selected features:", length(sel), "/", length(rv$all_meta_cols), "\n\n")

    if (length(sel) == 0) {
      cat("No features selected. Use the dropdown above to pick features.\n")
      return()
    }

    cat("Selected feature details:\n")
    for (col in sel) {
      hide_state <- rv$mini_hide_states[[col]] %||% FALSE
      coercion <- rv$type_coercions[[col]] %||% ""
      coercion_text <- if (nzchar(coercion)) coercion else "(no coercion)"
      cat(sprintf("- %s : Hide=%s ; Type=%s\n", col, ifelse(isTRUE(hide_state), "TRUE", "FALSE"), coercion_text))
    }
  })

  # Update type dropdowns when metadata loads (preserve current selection)
  observeEvent(rv$meta_sample, {
    req(rv$meta_sample)
    # Preserve current selection where possible when refreshing choices
    # Exclude patient_ID (and other internal columns) from selectable features
    new_choices <- setdiff(colnames(rv$meta_sample), c("patient_ID", "source", "run_date"))
    current_sel <- isolate(input$features_dropdown) %||% character(0)
    # Keep only selections that still exist in new choices
    sel_to_set <- intersect(current_sel, new_choices)
    if (length(sel_to_set) == 0) sel_to_set <- NULL
    updatePickerInput(session, "features_dropdown", choices = c("", new_choices), selected = sel_to_set)
  })

  # Track previous features_dropdown selection to detect deselections
  prev_features_dropdown <- reactiveVal(character(0))

  # When features are deselected from the dropdown, reset their coercion to original
  observeEvent(input$features_dropdown,
    {
      current <- input$features_dropdown %||% character(0)
      prev <- prev_features_dropdown()
      removed <- setdiff(prev, current)
      if (length(removed) > 0) {
        for (col in removed) {
          # Reset coercion to blank (original) if it was set
          if (!is.null(rv$type_coercions[[col]]) && nzchar(rv$type_coercions[[col]])) {
            rv$type_coercions[[col]] <- ""
            mini_type_id <- paste0("type_mini_", gsub("[^A-Za-z0-9]", "_", col))
            tryCatch({ updateSelectInput(session, mini_type_id, selected = "") }, error = function(e) {})
          }
          # Clear any persisted mini-hide flag so feature is not hidden after deselect
          rv$mini_hide_states[[col]] <- FALSE
        }
      }
      prev_features_dropdown(current)
    },
    ignoreNULL = FALSE
  )

  # Update pairing_var choices when metadata loads (preserve current selection)
  observeEvent(rv$meta_sample,
    {
      req(rv$meta_sample)
      all_cols <- colnames(rv$meta_sample)
      current_selection <- input$pairing_var

      # Preserve selection if it's still valid, otherwise reset
      if (!is.null(current_selection) && current_selection %in% all_cols) {
        updatePickerInput(session, "pairing_var", choices = c("", all_cols), selected = current_selection)
      } else {
        updatePickerInput(session, "pairing_var", choices = c("", all_cols), selected = "")
      }

      # Update trim_incomplete_pairs with categorical features only
      categorical_cols <- names(rv$meta_sample)[sapply(rv$meta_sample, function(col) {
        is.character(col) || is.factor(col)
      })]
      updatePickerInput(session, "trim_incomplete_pairs", choices = categorical_cols)
    },
    ignoreInit = TRUE
  )

  # Display pairing summary
  output$pairing_summary_ui <- renderUI({
    pairing_col <- input$pairing_var
    if (is.null(pairing_col) || !nzchar(pairing_col)) {
      return(tags$div(
        style = "margin-top: 20px;",
        tags$strong("Status: "),
        tags$span(style = "color: #999;", "Unpaired tests will be used")
      ))
    }

    req(rv$meta_sample)
    if (!(pairing_col %in% colnames(rv$meta_sample))) {
      return(tags$div(
        style = "margin-top: 20px;",
        tags$strong("Status: "),
        tags$span(style = "color: #d9534f;", "Invalid column selected")
      ))
    }

    vals <- rv$meta_sample[[pairing_col]]
    n_unique <- length(unique(vals[!is.na(vals)]))
    tbl <- table(vals, useNA = "ifany")

    tags$div(
      style = "margin-top: 20px;",
      tags$strong(style = "color: #5bc0de;", "Status: Paired tests enabled"),
      tags$p(
        style = "margin-top: 10px;",
        sprintf("Number of unique values: %d", n_unique)
      ),
      tags$pre(
        style = "background-color: #f5f5f5; padding: 10px; border-radius: 4px; max-height: 200px; overflow-y: auto;",
        paste(capture.output(print(tbl)), collapse = "\n")
      )
    )
  })

  # Paired testing indicators for Testing tab
  output$paired_test_indicator_ui <- renderUI({
    pairing_col <- input$pairing_var
    if (!is.null(pairing_col) && nzchar(pairing_col)) {
      tags$div(
        style = "background-color: #d9edf7; border: 1px solid #bce8f1; border-radius: 4px; padding: 10px; margin: 10px 0;",
        tags$strong(style = "color: #31708f;", "\u2713 Paired testing enabled"),
        tags$br(),
        tags$small(sprintf("Pairing by: %s", pairing_col))
      )
    }
  })

  # ========== METADATA SUBSETTING ==========

  # Track active rule IDs (use unique IDs to preserve input values)
  subset_rule_ids <- reactiveVal(integer(0))
  subset_next_id <- reactiveVal(1)

  # Render subsetting rules UI
  output$subsetting_rules_ui <- renderUI({
    req(rv$meta_sample_original)

    rule_ids <- subset_rule_ids()
    if (length(rule_ids) == 0) {
      return(tags$p(style = "color: #999; font-style: italic;", "Click 'Add Rule' to create a subsetting rule."))
    }

    meta_cols <- colnames(rv$meta_sample_original)

    rule_uis <- lapply(rule_ids, function(rule_id) {
      ns_id <- paste0("rule_", rule_id)

      # Preserve existing column selection if it exists
      current_col <- input[[paste0(ns_id, "_col")]]
      if (is.null(current_col)) current_col <- meta_cols[1]

      tags$div(
        id = paste0("rule_container_", rule_id),
        style = "border: 1px solid #ddd; padding: 10px; margin-bottom: 10px; border-radius: 4px; background-color: #f9f9f9;",
        fluidRow(
          column(
            5,
            selectInput(paste0(ns_id, "_col"), "Column", choices = meta_cols, selected = current_col)
          ),
          column(
            6,
            uiOutput(paste0(ns_id, "_values_ui"))
          ),
          column(
            1,
            br(),
            actionButton(paste0("remove_rule_", rule_id), "",
              icon = icon("trash"), class = "btn-danger btn-sm",
              style = "margin-top: 5px;"
            )
          )
        )
      )
    })

    tagList(rule_uis)
  })

  # Dynamically create value pickers for each rule
  observe({
    req(rv$meta_sample_original)
    rule_ids <- subset_rule_ids()
    if (length(rule_ids) == 0) {
      return()
    }

    lapply(rule_ids, function(rule_id) {
      ns_id <- paste0("rule_", rule_id)
      col_input_id <- paste0(ns_id, "_col")
      values_ui_id <- paste0(ns_id, "_values_ui")
      values_input_id <- paste0(ns_id, "_values")

      # Only re-render when column changes (not when values change)
      observeEvent(input[[col_input_id]],
        {
          selected_col <- input[[col_input_id]]
          if (is.null(selected_col) || !nzchar(selected_col)) {
            output[[values_ui_id]] <- renderUI(NULL)
            return()
          }

          req(selected_col %in% colnames(rv$meta_sample_original))
          col_data <- rv$meta_sample_original[[selected_col]]

          # Check if column is numeric
          if (is.numeric(col_data)) {
            # Numeric column: show operator + slider
            col_range <- range(col_data, na.rm = TRUE)
            col_mean <- mean(col_data, na.rm = TRUE)
            col_median <- median(col_data, na.rm = TRUE)
            step_size <- (col_range[2] - col_range[1]) / 100

            # Preserve existing selections if available
            current_operator <- isolate(input[[paste0(ns_id, "_operator")]])
            current_slider <- isolate(input[[paste0(ns_id, "_slider")]])
            if (is.null(current_operator)) current_operator <- "above"
            if (is.null(current_slider)) current_slider <- col_median

            operator_input_id <- paste0(ns_id, "_operator")
            slider_input_id <- paste0(ns_id, "_slider")
            mean_btn_id <- paste0(ns_id, "_mean_btn")
            median_btn_id <- paste0(ns_id, "_median_btn")

            output[[values_ui_id]] <- renderUI({
              tagList(
                radioButtons(
                  operator_input_id,
                  "Filter type",
                  choices = c(
                    "Above or equal (>=)" = "above",
                    "Below or equal (<=)" = "below",
                    "Range (internal)" = "range_internal",
                    "Range (external)" = "range_external"
                  ),
                  selected = current_operator,
                  inline = FALSE
                ),
                fluidRow(
                  column(6, actionButton(mean_btn_id, "Mean", class = "btn-sm btn-default", style = "width: 100%;")),
                  column(6, actionButton(median_btn_id, "Median", class = "btn-sm btn-default", style = "width: 100%;"))
                ),
                uiOutput(paste0(ns_id, "_slider_ui")),
                checkboxInput(paste0(ns_id, "_exclude_missing"), "Exclude NA values", value = isTRUE(isolate(input[[paste0(ns_id, "_exclude_missing")]])))
              )
            })

            # Render appropriate slider based on operator selection
            output[[paste0(ns_id, "_slider_ui")]] <- renderUI({
              req(input[[operator_input_id]])
              operator <- input[[operator_input_id]]

              if (operator %in% c("range_internal", "range_external")) {
                # Range slider with two handles
                current_range <- isolate(input[[slider_input_id]])
                if (is.null(current_range) || length(current_range) != 2) {
                  # Default to quartiles for initial range
                  current_range <- quantile(col_data, probs = c(0.25, 0.75), na.rm = TRUE)
                }
                sliderInput(
                  slider_input_id,
                  "Value range",
                  min = col_range[1],
                  max = col_range[2],
                  value = current_range,
                  step = step_size
                )
              } else {
                # Single-value slider
                current_val <- isolate(input[[slider_input_id]])
                if (is.null(current_val) || length(current_val) != 1) {
                  current_val <- col_median
                } else if (length(current_val) > 1) {
                  # If switching from range to single, use the midpoint
                  current_val <- mean(current_val)
                }
                sliderInput(
                  slider_input_id,
                  "Threshold value",
                  min = col_range[1],
                  max = col_range[2],
                  value = current_val,
                  step = step_size
                )
              }
            })

            # Add observers for mean/median buttons
            observeEvent(input[[mean_btn_id]],
              {
                operator <- input[[operator_input_id]]
                if (!is.null(operator) && operator %in% c("range_internal", "range_external")) {
                  # For range, center around mean
                  margin <- (col_range[2] - col_range[1]) * 0.25
                  updateSliderInput(session, slider_input_id, value = c(col_mean - margin, col_mean + margin))
                } else {
                  updateSliderInput(session, slider_input_id, value = col_mean)
                }
              },
              ignoreInit = TRUE
            )

            observeEvent(input[[median_btn_id]],
              {
                operator <- input[[operator_input_id]]
                if (!is.null(operator) && operator %in% c("range_internal", "range_external")) {
                  # For range, center around median
                  margin <- (col_range[2] - col_range[1]) * 0.25
                  updateSliderInput(session, slider_input_id, value = c(col_median - margin, col_median + margin))
                } else {
                  updateSliderInput(session, slider_input_id, value = col_median)
                }
              },
              ignoreInit = TRUE
            )
          } else {
            # Categorical column: show pickerInput + option to exclude missing/blank
            unique_vals <- sort(unique(as.character(col_data)))
            unique_vals <- unique_vals[!is.na(unique_vals)]

            # Use isolate to read current values without creating dependency
            current_values <- isolate(input[[values_input_id]])
            if (is.null(current_values)) {
              # Default to all values if no prior selection
              selected_vals <- unique_vals
            } else {
              # Keep existing selection (intersect to handle column changes)
              selected_vals <- intersect(current_values, unique_vals)
              if (length(selected_vals) == 0) selected_vals <- unique_vals
            }

            output[[values_ui_id]] <- renderUI({
              tagList(
                pickerInput(
                  values_input_id,
                  "Values to keep",
                  choices = unique_vals,
                  selected = selected_vals,
                  multiple = TRUE,
                  options = list(
                    `actions-box` = TRUE,
                    `selected-text-format` = "count > 3",
                    `deselect-all-text` = "Deselect All",
                    `select-all-text` = "Select All",
                    `none-selected-text` = "No values selected",
                    `live-search` = TRUE
                  )
                ),
                # (No exclude checkbox for categorical columns)
              )
            })
          }
        },
        ignoreNULL = TRUE,
        ignoreInit = FALSE
      )
    })
  })

  # Add rule button
  observeEvent(input$add_subset_rule, {
    new_id <- subset_next_id()
    subset_rule_ids(c(subset_rule_ids(), new_id))
    subset_next_id(new_id + 1)
  })

  # Remove rule buttons (dynamic observers for each rule ID)
  observe({
    rule_ids <- subset_rule_ids()
    if (length(rule_ids) == 0) {
      return()
    }

    lapply(rule_ids, function(rule_id) {
      observeEvent(input[[paste0("remove_rule_", rule_id)]],
        {
          # Remove this specific rule ID from the list
          current_ids <- subset_rule_ids()
          subset_rule_ids(setdiff(current_ids, rule_id))
        },
        ignoreInit = TRUE,
        ignoreNULL = TRUE,
        once = TRUE
      )
    })
  })

  # Apply subsetting
  observeEvent(input$apply_subsetting, {
    req(rv$meta_sample_original)

    if (!input$enable_subsetting) {
      rv$meta_sample <- rv$meta_sample_original
      rv$subsetting_enabled <- FALSE
      rv$subset_summary <- NULL
      showNotification("Subsetting disabled. All samples are now available.", type = "info")
      return()
    }

    rule_ids <- subset_rule_ids()
    if (length(rule_ids) == 0) {
      showNotification("No subsetting rules defined. Add at least one rule.", type = "warning")
      return()
    }

    # Get original metadata early for rule collection
    meta_orig <- rv$meta_sample_original

    # Collect all rules using rule IDs
    rules_list <- lapply(rule_ids, function(rule_id) {
      ns_id <- paste0("rule_", rule_id)
      col <- input[[paste0(ns_id, "_col")]]

      # Check if column is numeric or categorical
      if (!is.null(col) && col %in% colnames(meta_orig)) {
        col_data <- meta_orig[[col]]

        if (is.numeric(col_data)) {
          # Numeric rule
          operator <- input[[paste0(ns_id, "_operator")]]
          threshold <- input[[paste0(ns_id, "_slider")]]

          # Handle range operators
          if (!is.null(operator) && operator %in% c("range_internal", "range_external")) {
            # threshold should be a vector of length 2 for ranges
            if (!is.null(threshold) && length(threshold) == 2) {
              exclude_missing <- isTRUE(input[[paste0(ns_id, "_exclude_missing")]])
              list(
                column = col, type = "numeric", operator = operator,
                threshold_min = threshold[1], threshold_max = threshold[2],
                exclude_missing = exclude_missing
              )
            } else {
              NULL # Invalid range rule
            }
          } else {
            # Single threshold for above/below
            exclude_missing <- isTRUE(input[[paste0(ns_id, "_exclude_missing")]])
            list(
              column = col, type = "numeric", operator = operator, threshold = threshold,
              exclude_missing = exclude_missing
            )
          }
        } else {
          # Categorical rule
          vals <- input[[paste0(ns_id, "_values")]]
          list(column = col, type = "categorical", values = vals)
        }
      } else {
        NULL
      }
    })

    # Remove any invalid rules
    rules_list <- Filter(function(r) {
      if (is.null(r) || is.null(r$column)) {
        return(FALSE)
      }
      if (r$type == "categorical") {
        # Require at least one selected value for categorical rules
        (!is.null(r$values) && length(r$values) > 0)
      } else if (r$type == "numeric") {
        if (!is.null(r$operator) && r$operator %in% c("range_internal", "range_external")) {
          # Range rule validation
          !is.null(r$threshold_min) && !is.null(r$threshold_max)
        } else {
          # Single threshold rule validation
          !is.null(r$operator) && !is.null(r$threshold)
        }
      } else {
        FALSE
      }
    }, rules_list)

    if (length(rules_list) == 0) {
      showNotification("No valid subsetting rules. Please configure at least one rule.", type = "warning")
      return()
    }

    # Apply subsetting logic
    n_before <- nrow(meta_orig)

    if (input$subset_logic == "intersection") {
      # Intersection: keep rows that match ALL rules
      keep_mask <- rep(TRUE, nrow(meta_orig))
      for (rule in rules_list) {
        if (rule$type == "categorical") {
          col_vals <- as.character(meta_orig[[rule$column]])
          # Build condition: if values specified, match them; otherwise default to TRUE
          if (!is.null(rule$values) && length(rule$values) > 0) {
            cond <- col_vals %in% rule$values
          } else {
            cond <- rep(TRUE, length(col_vals))
          }
          keep_mask <- keep_mask & cond
        } else if (rule$type == "numeric") {
          col_vals <- meta_orig[[rule$column]]
          excl_missing <- !is.null(rule$exclude_missing) && isTRUE(rule$exclude_missing)
          if (rule$operator == "above") {
            if (excl_missing) {
              cond <- (col_vals >= rule$threshold) & !is.na(col_vals)
            } else {
              cond <- (col_vals >= rule$threshold) | is.na(col_vals)
            }
            keep_mask <- keep_mask & cond
          } else if (rule$operator == "below") {
            if (excl_missing) {
              cond <- (col_vals <= rule$threshold) & !is.na(col_vals)
            } else {
              cond <- (col_vals <= rule$threshold) | is.na(col_vals)
            }
            keep_mask <- keep_mask & cond
          } else if (rule$operator == "range_internal") {
            if (excl_missing) {
              cond <- (col_vals >= rule$threshold_min & col_vals <= rule$threshold_max) & !is.na(col_vals)
            } else {
              cond <- (col_vals >= rule$threshold_min & col_vals <= rule$threshold_max) | is.na(col_vals)
            }
            keep_mask <- keep_mask & cond
          } else if (rule$operator == "range_external") {
            if (excl_missing) {
              cond <- (col_vals <= rule$threshold_min | col_vals >= rule$threshold_max) & !is.na(col_vals)
            } else {
              cond <- (col_vals <= rule$threshold_min | col_vals >= rule$threshold_max) | is.na(col_vals)
            }
            keep_mask <- keep_mask & cond
          }
        }
      }
      meta_filtered <- meta_orig[keep_mask, , drop = FALSE]
    } else {
      # Union: keep rows that match ANY rule
      keep_mask <- rep(FALSE, nrow(meta_orig))
      for (rule in rules_list) {
        if (rule$type == "categorical") {
          col_vals <- as.character(meta_orig[[rule$column]])
          if (!is.null(rule$values) && length(rule$values) > 0) {
            cond <- col_vals %in% rule$values
          } else {
            cond <- rep(FALSE, length(col_vals))
          }
          keep_mask <- keep_mask | cond
        } else if (rule$type == "numeric") {
          col_vals <- meta_orig[[rule$column]]
          if (rule$operator == "above") {
            keep_mask <- keep_mask | (col_vals >= rule$threshold & !is.na(col_vals))
          } else if (rule$operator == "below") {
            keep_mask <- keep_mask | (col_vals <= rule$threshold & !is.na(col_vals))
          } else if (rule$operator == "range_internal") {
            # Keep values within range (inclusive)
            keep_mask <- keep_mask | ((col_vals >= rule$threshold_min & col_vals <= rule$threshold_max) & !is.na(col_vals))
          } else if (rule$operator == "range_external") {
            # Keep values outside range (inclusive of boundaries)
            keep_mask <- keep_mask | ((col_vals <= rule$threshold_min | col_vals >= rule$threshold_max) & !is.na(col_vals))
          }
        }
      }
      meta_filtered <- meta_orig[keep_mask, , drop = FALSE]
    }

    n_after <- nrow(meta_filtered)
    n_excluded <- n_before - n_after

    if (n_after == 0) {
      showNotification("Subsetting removed all samples. Please adjust your rules.", type = "error")
      return()
    }

    # Trim incomplete pairs if requested
    n_before_trim <- n_after
    n_trimmed <- 0
    trim_features <- input$trim_incomplete_pairs
    pairing_var <- input$pairing_var

    if (!is.null(trim_features) && length(trim_features) > 0 &&
      !is.null(pairing_var) && nzchar(pairing_var) &&
      pairing_var %in% colnames(meta_filtered)) {
      # For each trim feature, identify complete pairs
      for (trim_feat in trim_features) {
        if (trim_feat %in% colnames(meta_filtered)) {
          # Get all pairing IDs
          pairing_ids <- unique(meta_filtered[[pairing_var]])

          # For each pairing ID, check if all levels of trim_feat are present
          trim_feat_levels <- unique(meta_filtered[[trim_feat]])

          complete_pair_ids <- character(0)
          for (pair_id in pairing_ids) {
            pair_rows <- meta_filtered[meta_filtered[[pairing_var]] == pair_id, ]
            pair_trim_levels <- unique(pair_rows[[trim_feat]])

            # Check if this pair has all levels of the trim feature
            if (length(pair_trim_levels) == length(trim_feat_levels)) {
              complete_pair_ids <- c(complete_pair_ids, pair_id)
            }
          }

          # Keep only complete pairs
          meta_filtered <- meta_filtered[meta_filtered[[pairing_var]] %in% complete_pair_ids, , drop = FALSE]
        }
      }

      n_after_trim <- nrow(meta_filtered)
      n_trimmed <- n_before_trim - n_after_trim
      n_after <- n_after_trim
      n_excluded <- n_before - n_after

      if (n_after == 0) {
        showNotification("Trimming incomplete pairs removed all samples. Please adjust your settings.", type = "error")
        return()
      }
    }

    # Generate unique subset ID based on system time
    # Use default ID if no samples were dropped (equivalent to no subsetting)
    if (n_excluded == 0) {
      subset_id <- "000000000"
    } else {
      time_seed <- as.integer(Sys.time())
      set.seed(time_seed)
      subset_id <- paste0(sample(c(LETTERS, letters, rep(0:9, 2)), size = 9, replace = TRUE), collapse = "")
    }

    # Update reactive values
    rv$meta_sample <- meta_filtered
    rv$subsetting_enabled <- TRUE
    rv$subset_rules <- rules_list
    rv$subset_id <- subset_id
    rv$subset_summary <- list(
      n_before = n_before,
      n_after = n_after,
      n_excluded = n_excluded,
      n_trimmed = n_trimmed,
      subset_id = subset_id
    )

    # Build notification message
    notif_msg <- sprintf("Subsetting applied: %d samples included, %d excluded", n_after, n_excluded)
    if (n_trimmed > 0) {
      notif_msg <- sprintf("%s (%d trimmed for incomplete pairs)", notif_msg, n_trimmed)
    }
    notif_msg <- sprintf("%s. [Subset ID: %s]", notif_msg, subset_id)

    showNotification(notif_msg, type = "message", duration = 5)
  })

  # Reset subsetting when disabled
  observeEvent(input$enable_subsetting, {
    if (!input$enable_subsetting && !is.null(rv$meta_sample_original)) {
      rv$meta_sample <- rv$meta_sample_original
      rv$subsetting_enabled <- FALSE
      rv$subset_summary <- NULL
      # Reset subset_id to default when subsetting is disabled
      rv$subset_id <- "000000000"
      subset_rule_ids(integer(0))
      subset_next_id(1)
    }
  })

  # Render dynamic export buttons
  output$export_buttons_ui <- renderUI({
    is_subsetted <- isTRUE(rv$subsetting_enabled)

    if (is_subsetted) {
      tagList(
        downloadButton("export_subset_meta", "Export Subset Metadata"),
        downloadButton("export_subset_frequencies", "Export Subset Cluster Frequencies"),
        downloadButton("export_subset_counts", "Export Subset Cluster Counts")
      )
    } else {
      tagList(
        downloadButton("export_subset_meta", "Export Metadata"),
        downloadButton("export_subset_frequencies", "Export Cluster Frequencies"),
        downloadButton("export_subset_counts", "Export Cluster Counts")
      )
    }
  })

  # Render subsetting preview
  output$subset_preview_ui <- renderUI({
    if (is.null(rv$subset_summary)) {
      return(tags$div(
        style = "color: #999; font-style: italic;",
        "Configure rules and click 'Apply Subsetting' to see preview."
      ))
    }

    summ <- rv$subset_summary

    # Build preview table showing sample counts per column
    preview_cols <- unique(sapply(rv$subset_rules, function(r) r$column))
    preview_tables <- lapply(preview_cols, function(col) {
      if (!col %in% colnames(rv$meta_sample)) {
        return(NULL)
      }
      tbl <- table(rv$meta_sample[[col]], useNA = "ifany")
      tags$div(
        tags$strong(paste0(col, ":")),
        tags$pre(
          style = "background-color: #f5f5f5; padding: 5px; border-radius: 4px; font-size: 11px; max-height: 150px; overflow-y: auto;",
          paste(capture.output(print(tbl)), collapse = "\n")
        )
      )
    })

    tags$div(
      tags$div(
        style = "background-color: #dff0d8; border: 1px solid #d6e9c6; border-radius: 4px; padding: 10px; margin-bottom: 10px;",
        tags$strong(style = "color: #3c763d;", "Subsetting Applied"),
        tags$p(style = "margin: 5px 0;", sprintf("Samples before: %d", summ$n_before)),
        tags$p(style = "margin: 5px 0;", sprintf("Samples after: %d", summ$n_after)),
        tags$p(style = "margin: 5px 0;", sprintf("Samples excluded: %d", summ$n_excluded)),
        if (!is.null(summ$n_trimmed) && summ$n_trimmed > 0) {
          tags$p(style = "margin: 5px 0; color: #8a6d3b;", sprintf("(includes %d trimmed for incomplete pairs)", summ$n_trimmed))
        } else {
          NULL
        },
        tags$p(style = "margin: 5px 0; font-family: monospace; color: #31708f;", sprintf("Subset ID: %s", summ$subset_id))
      ),
      tags$h5("Sample Distribution"),
      tagList(preview_tables)
    )
  })

  # Reactive to check if subset is active
  output$hasSubset <- reactive({
    # Always show button - even for "000000000" (full dataset)
    !is.null(rv$subset_id)
  })
  outputOptions(output, "hasSubset", suspendWhenHidden = FALSE)

  # ==== Cluster Annotation Engine ====
  
  # Reactive for available clusters
  available_clusters <- reactive({
    if (!is.null(rv$clusters) && !is.null(rv$clusters$assignments)) {
      sort(unique(rv$clusters$assignments))
    } else {
      character(0)
    }
  })
  
  # Check if annotations exist
  output$hasAnnotations <- reactive({
    length(rv$celltype_annotations) > 0
  })
  outputOptions(output, "hasAnnotations", suspendWhenHidden = FALSE)
  
  # Render annotation UI
  output$celltype_annotations_ui <- renderUI({
    if (length(rv$celltype_annotations) == 0) {
      return(tags$p(style = "color: #999; font-style: italic;", "Click 'Add Cell Type' to create annotations."))
    }
    
    clusters <- available_clusters()
    if (length(clusters) == 0) {
      return(tags$p(style = "color: #dc3545;", "No clusters available. Load data with cluster assignments."))
    }
    
    ann_ids <- names(rv$celltype_annotations)
    n_annotations <- length(ann_ids)
    
    # Create grid layout with 2 columns
    n_cols <- 2
    n_rows <- ceiling(n_annotations / n_cols)
    
    rows <- lapply(seq_len(n_rows), function(row_idx) {
      start_idx <- (row_idx - 1) * n_cols + 1
      end_idx <- min(start_idx + n_cols - 1, n_annotations)
      
      cols <- lapply(start_idx:end_idx, function(idx) {
        ann_id <- ann_ids[idx]
        ann <- rv$celltype_annotations[[ann_id]]
        
        column(
          6,  # 12/2 = 6 for 2 columns
          tags$div(
            id = paste0("annotation_container_", ann_id),
            style = "border: 1px solid #ddd; padding: 8px; margin-bottom: 8px; border-radius: 4px; background-color: #f9f9f9;",
            fluidRow(
              column(
                5,
                textInput(paste0("celltype_name_", ann_id), "Cell Type", value = ann$name, placeholder = "e.g., T cells")
              ),
              column(
                6,
                pickerInput(
                  paste0("celltype_clusters_", ann_id),
                  "Clusters",
                  choices = as.character(clusters),
                  selected = ann$clusters,
                  multiple = TRUE,
                  options = list(
                    `actions-box` = TRUE,
                    `selected-text-format` = "count > 3",
                    `count-selected-text` = "{0} clusters"
                  )
                )
              ),
              column(
                1,
                br(),
                actionButton(paste0("remove_celltype_", ann_id), "",
                  icon = icon("trash"), class = "btn-danger btn-sm",
                  style = "margin-top: 5px;"
                )
              )
            )
          )
        )
      })
      
      fluidRow(cols)
    })
    
    tagList(rows)
  })
  
  # Add celltype button
  observeEvent(input$add_celltype, {
    clusters <- available_clusters()
    if (length(clusters) == 0) {
      showNotification("No clusters available. Load data with cluster assignments first.", type = "error")
      return()
    }
    
    # Initialize with cluster->cluster mapping if this is the first celltype
    default_clusters <- if (length(rv$celltype_annotations) == 0) {
      as.character(clusters)
    } else {
      character(0)
    }
    
    new_id <- celltype_id_counter() + 1
    celltype_id_counter(new_id)
    
    rv$celltype_annotations[[as.character(new_id)]] <- list(
      id = new_id,
      name = paste0("Cell Type ", new_id),
      clusters = default_clusters
    )
  })
  
  # Remove celltype buttons (dynamic observers)
  observe({
    if (length(rv$celltype_annotations) == 0) return()
    
    lapply(names(rv$celltype_annotations), function(ann_id) {
      observeEvent(input[[paste0("remove_celltype_", ann_id)]],
        {
          # Remove this annotation
          rv$celltype_annotations[[ann_id]] <- NULL
        },
        ignoreInit = TRUE,
        ignoreNULL = TRUE,
        once = TRUE
      )
    })
  })
  
  # Update annotations in real-time as user makes changes
  observe({
    if (length(rv$celltype_annotations) == 0) return()
    
    lapply(names(rv$celltype_annotations), function(ann_id) {
      # Update name
      observeEvent(input[[paste0("celltype_name_", ann_id)]], {
        if (!is.null(rv$celltype_annotations[[ann_id]])) {
          rv$celltype_annotations[[ann_id]]$name <- input[[paste0("celltype_name_", ann_id)]]
        }
      }, ignoreInit = TRUE)
      
      # Update cluster assignments
      observeEvent(input[[paste0("celltype_clusters_", ann_id)]], {
        if (!is.null(rv$celltype_annotations[[ann_id]])) {
          rv$celltype_annotations[[ann_id]]$clusters <- input[[paste0("celltype_clusters_", ann_id)]]
        }
      }, ignoreInit = TRUE)
    })
  })
  
  # Update dropdown choices to prevent duplicate cluster assignments
  observe({
    if (length(rv$celltype_annotations) == 0) return()
    
    all_clusters <- as.character(available_clusters())
    ann_ids <- names(rv$celltype_annotations)
    
    # For each annotation, update its dropdown choices
    lapply(ann_ids, function(current_ann_id) {
      # Get clusters selected in OTHER annotations
      clusters_selected_elsewhere <- character(0)
      for (other_ann_id in setdiff(ann_ids, current_ann_id)) {
        selected <- input[[paste0("celltype_clusters_", other_ann_id)]]
        if (!is.null(selected) && length(selected) > 0) {
          clusters_selected_elsewhere <- c(clusters_selected_elsewhere, selected)
        }
      }
      
      # Get currently selected clusters for THIS annotation
      current_selected <- input[[paste0("celltype_clusters_", current_ann_id)]]
      
      # Available choices = all clusters - clusters selected elsewhere
      available_choices <- setdiff(all_clusters, clusters_selected_elsewhere)
      
      # Update the picker with new choices, preserving current selection
      updatePickerInput(
        session = session,
        inputId = paste0("celltype_clusters_", current_ann_id),
        choices = available_choices,
        selected = current_selected
      )
    })
  })
  
  # Apply annotations button
  observeEvent(input$apply_annotations, {
    if (length(rv$celltype_annotations) == 0) {
      showNotification("No cell type annotations defined.", type = "warning")
      return()
    }
    
    # Validate: check for empty names
    invalid_names <- sapply(rv$celltype_annotations, function(ann) {
      is.null(ann$name) || nchar(trimws(ann$name)) == 0
    })
    if (any(invalid_names)) {
      showNotification("All cell types must have non-empty names.", type = "error")
      return()
    }
    
    # Validate: check for duplicate names
    all_names <- sapply(rv$celltype_annotations, function(ann) trimws(ann$name))
    if (any(duplicated(all_names))) {
      showNotification("Cell type names must be unique.", type = "error")
      return()
    }
    
    # Build cluster_map from annotations
    cluster_map_list <- list()
    for (ann_id in names(rv$celltype_annotations)) {
      ann <- rv$celltype_annotations[[ann_id]]
      celltype_name <- trimws(ann$name)
      assigned_clusters <- ann$clusters
      
      if (length(assigned_clusters) > 0) {
        for (cl in assigned_clusters) {
          cluster_map_list[[cl]] <- celltype_name
        }
      }
    }
    
    # Handle unassigned clusters: map to themselves
    all_clusters <- available_clusters()
    unassigned_clusters <- setdiff(as.character(all_clusters), names(cluster_map_list))
    for (cl in unassigned_clusters) {
      cluster_map_list[[cl]] <- as.character(cl)
    }
    
    # Convert to data frame
    if (length(cluster_map_list) > 0) {
      rv$cluster_map <- data.frame(
        cluster = as.integer(names(cluster_map_list)),
        celltype = unlist(cluster_map_list),
        stringsAsFactors = FALSE
      )
      
      showNotification(
        sprintf("Applied annotations: %d cell types, %d clusters assigned", 
                length(rv$celltype_annotations), length(all_clusters)),
        type = "message"
      )
    } else {
      rv$cluster_map <- NULL
      showNotification("No cluster assignments made.", type = "warning")
    }
  })

  # Export subset metadata
  output$export_subset_meta <- downloadHandler(
    filename = function() {
      subset_id <- rv$subset_id %||% "000000000"
      paste0("subset_metadata_", subset_id, ".csv")
    },
    content = function(file) {
      req(rv$meta_sample)
      write.csv(rv$meta_sample, file, row.names = FALSE)
    }
  )

  # Export subset cluster frequencies (abundance data)
  output$export_subset_frequencies <- downloadHandler(
    filename = function() {
      subset_id <- rv$subset_id %||% "000000000"
      paste0("subset_cluster_frequencies_", subset_id, ".csv")
    },
    content = function(file) {
      req(rv$abundance_sample, rv$meta_sample)

      # Get patient_IDs from current subsetted metadata in exact order
      subset_ids <- rv$meta_sample$patient_ID

      # Filter and reorder abundance data to match metadata row order exactly
      abundance_subset <- rv$abundance_sample[match(subset_ids, rownames(rv$abundance_sample)), , drop = FALSE]

      # Add patient_ID as first column (in same order as metadata)
      abundance_export <- data.frame(patient_ID = subset_ids, abundance_subset, check.names = FALSE)

      write.csv(abundance_export, file, row.names = FALSE)
    }
  )

  # Export subset cluster counts
  output$export_subset_counts <- downloadHandler(
    filename = function() {
      subset_id <- rv$subset_id %||% "000000000"
      paste0("subset_cluster_counts_", subset_id, ".csv")
    },
    content = function(file) {
      req(rv$counts_sample, rv$meta_sample)

      # Get patient_IDs from current subsetted metadata in exact order
      subset_ids <- rv$meta_sample$patient_ID

      # Filter and reorder counts data to match metadata row order exactly
      counts_subset <- rv$counts_sample[match(subset_ids, rownames(rv$counts_sample)), , drop = FALSE]

      # Add patient_ID as first column (in same order as metadata)
      counts_export <- data.frame(patient_ID = subset_ids, counts_subset, check.names = FALSE)

      write.csv(counts_export, file, row.names = FALSE)
    }
  )

  # Paired testing indicators for Categorical tab
  output$paired_cat_indicator_ui <- renderUI({
    pairing_col <- input$pairing_var
    if (!is.null(pairing_col) && nzchar(pairing_col)) {
      tags$div(
        style = "background-color: #d9edf7; border: 1px solid #bce8f1; border-radius: 4px; padding: 10px; margin: 10px 0;",
        tags$strong(style = "color: #31708f;", "\u2713 Paired testing enabled"),
        tags$br(),
        tags$small(sprintf("Pairing by: %s", pairing_col))
      )
    }
  })

  # Check if selected categorical variable in Testing tab can be paired
  test_can_pair <- reactive({
    pairing_var <- input$pairing_var
    group_var <- input$group_var

    # If no pairing variable or no group variable, can't pair
    if (is.null(pairing_var) || !nzchar(pairing_var) || is.null(group_var) || !nzchar(group_var)) {
      return(FALSE)
    }

    # Check if variables exist in data
    if (!exists("rv") || is.null(rv$meta_sample)) {
      return(FALSE)
    }
    if (!(pairing_var %in% colnames(rv$meta_sample)) || !(group_var %in% colnames(rv$meta_sample))) {
      return(FALSE)
    }

    # Check if any pairing values appear in multiple groups
    df <- rv$meta_sample[, c(pairing_var, group_var)]
    df <- df[complete.cases(df), ]

    if (nrow(df) < 2) {
      return(FALSE)
    }

    # For each unique pairing value, count how many different groups it appears in
    pair_group_counts <- df %>%
      dplyr::group_by(!!rlang::sym(pairing_var)) %>%
      dplyr::summarise(n_groups = dplyr::n_distinct(!!rlang::sym(group_var)), .groups = "drop")

    # If any pairing value appears in 2+ groups, pairing is possible
    any(pair_group_counts$n_groups >= 2)
  })

  # Render test type options for Testing tab based on pairing feasibility
  output$test_type_ui <- renderUI({
    pairing_var <- input$pairing_var
    can_pair <- test_can_pair()

    # If pairing is enabled AND possible for selected variable, show paired tests
    if (!is.null(pairing_var) && nzchar(pairing_var) && can_pair) {
      radioButtons("test_type", "Test",
        choices = c("Pairwise Wilcoxon", "Friedman (multi-group paired)", "Spearman (continuous)"),
        selected = "Pairwise Wilcoxon"
      )
    } else {
      # Otherwise show unpaired tests
      radioButtons("test_type", "Test",
        choices = c("Pairwise Wilcoxon", "Kruskal–Wallis (multi-group unpaired)", "Spearman (continuous)"),
        selected = "Pairwise Wilcoxon"
      )
    }
  })

  # Render test type options for Categorical tab based on pairing feasibility
  output$cat_test_type_ui <- renderUI({
    pairing_var <- input$pairing_var
    can_pair <- cat_can_pair()

    # If pairing is enabled AND possible for selected variable, show paired tests
    if (!is.null(pairing_var) && nzchar(pairing_var) && can_pair) {
      radioButtons("cat_test_type", "Test",
        choices = c("Pairwise Wilcoxon", "Friedman (multi-group paired)"),
        selected = "Pairwise Wilcoxon"
      )
    } else {
      # Otherwise show unpaired tests
      radioButtons("cat_test_type", "Test",
        choices = c("Pairwise Wilcoxon", "Kruskal–Wallis (multi-group unpaired)"),
        selected = "Pairwise Wilcoxon"
      )
    }
  })

  # Check if selected categorical variable can be paired
  cat_can_pair <- reactive({
    pairing_var <- input$pairing_var
    group_var <- input$cat_group_var

    # If no pairing variable or no group variable, can't pair
    if (is.null(pairing_var) || !nzchar(pairing_var) || is.null(group_var) || !nzchar(group_var)) {
      return(FALSE)
    }

    # Check if variables exist in data
    if (!exists("rv") || is.null(rv$meta_sample)) {
      return(FALSE)
    }
    if (!(pairing_var %in% colnames(rv$meta_sample)) || !(group_var %in% colnames(rv$meta_sample))) {
      return(FALSE)
    }

    # Check if any pairing values appear in multiple groups
    df <- rv$meta_sample[, c(pairing_var, group_var)]
    df <- df[complete.cases(df), ]

    if (nrow(df) < 2) {
      return(FALSE)
    }

    # For each unique pairing value, count how many different groups it appears in
    pair_group_counts <- df %>%
      dplyr::group_by(!!rlang::sym(pairing_var)) %>%
      dplyr::summarise(n_groups = dplyr::n_distinct(!!rlang::sym(group_var)), .groups = "drop")

    # If any pairing value appears in 2+ groups, pairing is possible
    any(pair_group_counts$n_groups >= 2)
  })

  # Display warning if pairing not possible for selected variable
  output$cat_pairing_check_ui <- renderUI({
    pairing_var <- input$pairing_var
    group_var <- input$cat_group_var

    # Only show if pairing is enabled but not possible for this variable
    if (!is.null(pairing_var) && nzchar(pairing_var) && !is.null(group_var) && nzchar(group_var)) {
      can_pair <- cat_can_pair()
      if (!can_pair) {
        tags$div(
          style = "background-color: #fcf8e3; border: 1px solid #faebcc; border-radius: 4px; padding: 8px; margin: 8px 0;",
          tags$small(
            style = "color: #8a6d3b;",
            "\u26A0 Selected variable cannot be paired - unpaired test will be used"
          )
        )
      }
    }
  })

  # Render cat_points UI based on whether pairing is possible
  output$cat_points_ui <- renderUI({
    pairing_var <- input$pairing_var
    can_pair <- cat_can_pair()

    # Show "Draw with connection" only if pairing is enabled AND possible for selected variable
    if (!is.null(pairing_var) && nzchar(pairing_var) && can_pair) {
      radioButtons("cat_points", "Show data points",
        choices = c("Draw" = "draw", "Draw with jitter" = "jitter", "Draw with connection" = "connected", "Do not draw" = "none"),
        selected = "draw"
      )
    } else {
      radioButtons("cat_points", "Show data points",
        choices = c("Draw" = "draw", "Draw with jitter" = "jitter", "Do not draw" = "none"),
        selected = "draw"
      )
    }
  })

  # Update all dropdowns when global settings change
  observeEvent(list(rv$available_features, rv$meta_sample),
    {
      req(rv$meta_sample, rv$data_ready)

      meta_cols <- colnames(rv$meta_sample)
      categorical_choices <- sort(meta_cols[sapply(rv$meta_sample, function(x) is.factor(x) || is.character(x))])
      continuous_choices <- sort(meta_cols[sapply(rv$meta_sample, is.numeric)])

      # Base exclusions
      base_exclude_outcomes <- c("cluster", "patient_ID", "run_date", "source")
      base_exclude_predictors <- c("patient_ID", "source", "run_date")

      # Apply base exclusions
      categorical_choices <- setdiff(categorical_choices, base_exclude_outcomes)
      continuous_choices <- setdiff(continuous_choices, base_exclude_outcomes)
      predictor_choices <- setdiff(meta_cols, base_exclude_predictors)

      # Apply global settings filters
      categorical_outcomes <- filter_by_global_settings(categorical_choices)
      continuous_outcomes <- filter_by_global_settings(continuous_choices)
      predictors <- filter_by_global_settings(predictor_choices)

      # Add cluster back to predictors if not already there
      if (!("cluster" %in% predictors)) {
        predictors <- c(predictors, "cluster")
      }

      # Update Feature Selection - allow both categorical and continuous outcomes
      fs_outcomes <- sort(c(categorical_outcomes, continuous_outcomes))
      updatePickerInput(session, "fs_outcome", choices = fs_outcomes)
      updatePickerInput(session, "fs_predictors", choices = predictors)

      # Update Classification
      updatePickerInput(session, "lm_outcome", choices = categorical_outcomes)
      updatePickerInput(session, "lm_predictors", choices = predictors)

      # Update Regression
      updatePickerInput(session, "reg_outcome", choices = continuous_outcomes)
      updatePickerInput(session, "reg_predictors", choices = predictors)

      # Update Testing tab
      updatePickerInput(session, "group_var", choices = c("", filter_by_global_settings(categorical_choices)))
      updatePickerInput(session, "cont_var", choices = c("", filter_by_global_settings(continuous_choices)))

      # Update Categorical tab
      updatePickerInput(session, "cat_group_var", choices = c("", filter_by_global_settings(categorical_choices)))

      # Update Continuous tab
      updatePickerInput(session, "cont_group_var", choices = c("", filter_by_global_settings(continuous_choices)))
    },
    ignoreInit = TRUE,
    ignoreNULL = FALSE
  )

  observeEvent(input$main_tab, {
    if (!isTRUE(rv$data_ready) && !identical(input$main_tab, "Home")) {
      updateNavbarPage(session, "main_tab", selected = "Home")
    }
  })

  get_cluster_clusters <- function() colnames(rv$clusters$abundance)

  # Helper function to filter choices based on global settings
  filter_by_global_settings <- function(choices) {
    if (length(rv$available_features) == 0) {
      return(character(0))
    }
    intersect(choices, rv$available_features)
  }

  # keep choices in sync with data
  observe({
    updatePickerInput(
      session, "fs_cluster_subset",
      choices = get_cluster_clusters(),
      selected = NULL
    )
  })

  # UI-facing flag for conditionalPanel (no nested reactive)
  output$hasClusterMap <- reactive({
    !is.null(rv$cluster_map) && all(c("cluster", "celltype") %in% names(rv$cluster_map))
  })
  outputOptions(output, "hasClusterMap", suspendWhenHidden = FALSE)

  observeEvent(rv$cluster_map,
    {
      cm <- rv$cluster_map
      if (!is.null(cm) && all(c("cluster", "celltype") %in% names(cm))) {
        # nothing to do
      } else {
        updatePickerInput(session, "test_entity", selected = "Clusters")
      }
    },
    ignoreInit = TRUE
  )

  observeEvent(rv$cluster_map,
    {
      cm <- rv$cluster_map
      if (!is.null(cm) && all(c("cluster", "celltype") %in% names(cm))) {
        # nothing to do
      } else {
        updatePickerInput(session, "cat_entity", selected = "Clusters")
      }
    },
    ignoreInit = TRUE
  )

  # Auto-detect categorical vs continuous metadata
  observeEvent(rv$meta_cell,
    {
      meta_cols <- colnames(rv$meta_cell)

      categorical_choices <- sort(meta_cols[sapply(rv$meta_cell, function(x) {
        is.character(x) || is.factor(x)
      })])
      # Exclude cluster, patient_ID, run_date, and source from categorical outcome choices
      categorical_choices <- setdiff(categorical_choices, c("cluster", "patient_ID", "run_date", "source"))
      # Apply global settings filter
      categorical_choices <- filter_by_global_settings(categorical_choices)

      continuous_choices <- sort(meta_cols[sapply(rv$meta_cell, function(x) {
        is.numeric(x) || is.integer(x)
      })])
      # Apply global settings filter
      continuous_choices <- filter_by_global_settings(continuous_choices)

      updatePickerInput(session, "group_var",
        choices = c("", categorical_choices), selected = ""
      )
      updatePickerInput(session, "cont_var",
        choices = c("", continuous_choices), selected = ""
      )
    },
    ignoreInit = TRUE
  )

  observeEvent(rv$meta_cell,
    {
      meta_cols <- colnames(rv$meta_cell)
      continuous_choices <- sort(meta_cols[sapply(rv$meta_cell, is.numeric)])
      # Apply global settings filter
      continuous_choices <- filter_by_global_settings(continuous_choices)
      updatePickerInput(session, "cont_group_var",
        choices = c("", continuous_choices),
        selected = ""
      )
    },
    ignoreInit = TRUE
  )

  observeEvent(rv$meta_cell,
    {
      meta_cols <- sort(colnames(rv$meta_cell))
      updatePickerInput(session, "model_outcome", choices = meta_cols)
      updatePickerInput(session, "model_predictors", choices = meta_cols)
      updatePickerInput(session, "model_covariates", choices = meta_cols)
      updatePickerInput(session, "model_random", choices = meta_cols)
    },
    ignoreInit = TRUE
  )

  # Launch embedding modules as soon as data is ready
  observeEvent(list(rv$UMAP, rv$data_ready),
    {
      req(rv$data_ready, rv$UMAP, rv$expr, rv$meta_cell, rv$clusters)
      message("Launching UMAP module")
      EmbeddingServer(
        "umap", "UMAP",
        reactive(rv$UMAP$coords),
        reactive(rv$expr),
        reactive(rv$meta_cell),
        reactive(rv$clusters),
        reactive(rv$cluster_map),
        reactive(input$main_tab)
      )
    },
    ignoreInit = TRUE
  )

  observeEvent(list(rv$tSNE, rv$data_ready),
    {
      req(rv$data_ready, rv$tSNE, rv$expr, rv$meta_cell, rv$clusters)
      message("Launching tSNE module")
      EmbeddingServer(
        "tsne", "tSNE",
        reactive(rv$tSNE$coords),
        reactive(rv$expr),
        reactive(rv$meta_cell),
        reactive(rv$clusters),
        reactive(rv$cluster_map),
        reactive(input$main_tab)
      )
    },
    ignoreInit = TRUE
  )

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
  output$meta_overview <- renderTable(
    {
      req(rv$meta_cell)
      data.frame(
        name = colnames(rv$meta_cell),
        type = sapply(rv$meta_cell, function(x) class(x)[1]),
        example = sapply(rv$meta_cell, function(x) paste(utils::head(unique(x), 6), collapse = ", "))
      )
    },
    sanitize.text.function = function(x) x
  )

  heatmap_obj <- reactive({
    req(rv$cluster_heat)
    M <- rv$cluster_heat

    # Clean row/col names
    if (!is.null(rownames(M))) rownames(M) <- gsub("\\n", " ", rownames(M))
    if (!is.null(colnames(M))) colnames(M) <- gsub("\\n", " ", colnames(M))

    # Annotation
    ranno <- NULL
    if (!is.null(rv$pop_size)) {
      size_vals <- rv$pop_size[, 1]
      size_col_fun <- circlize::colorRamp2(
        c(min(size_vals, na.rm = TRUE), max(size_vals, na.rm = TRUE)),
        c("white", "red")
      )
      ranno <- rowAnnotation(
        Size = size_vals,
        col = list(Size = size_col_fun),
        gp = grid::gpar(col = "black", lwd = 0.5)
      )
    }

    # Palette
    palette_choice <- switch(input$heatmap_theme,
      "viridis" = {
        col_seq <- seq(min(M), max(M), length.out = n <- 256)
        circlize::colorRamp2(col_seq, viridis(n))
      },
      "heat" = {
        col_seq <- seq(min(M), max(M), length.out = n <- 256)
        circlize::colorRamp2(
          col_seq,
          colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdYlBu")))(n)
        )
      },
      "greyscale" = {
        col_seq <- seq(min(M), max(M), length.out = n <- 256)
        circlize::colorRamp2(
          col_seq,
          colorRampPalette(c("black", "grey50", "white"))(n)
        )
      }
    )

    # Return the Heatmap object
    Heatmap(
      M,
      name = "expr",
      cluster_rows = isTRUE(input$cluster_rows),
      cluster_columns = isTRUE(input$cluster_columns),
      right_annotation = ranno,
      row_names_side = "left",
      rect_gp = gpar(lwd = 0.33, col = "black"),
      col = palette_choice
    )
  })

  output$cluster_heatmap <- renderPlot({
    draw(heatmap_obj())
  })

  output$export_heatmap_pdf <- downloadHandler(
    filename = function() {
      paste0("cluster_heatmap.pdf")
    },
    content = function(file) {
      # Get the underlying matrix from your reactive
      M <- rv$cluster_heat # or heatmap_obj()$matrix if you store it there

      # Calculate dimensions using the same factors as fcs_plot_heatmap
      pdf_width <- (ncol(M) * 0.33) + 3.25
      pdf_height <- (nrow(M) * 0.3) + 2.25

      pdf(file, width = pdf_width, height = pdf_height)
      draw(heatmap_obj())
      dev.off()
    },
    contentType = "application/pdf"
  )

  # Store results + adj_col name from the run
  run_tests <- eventReactive(input$run_test, {
    req(rv$meta_sample, rv$abundance_sample)

    test_type_run <- input$test_type
    p_adj_method_run <- input$p_adj_method
    group_var_run <- input$group_var
    cont_var_run <- input$cont_var
    test_entity_run <- input$test_entity
    pairing_var_run <- input$pairing_var

    abund0 <- rv$abundance_sample
    if (is.null(abund0)) {
      showNotification("No abundance matrix available.", type = "error")
      return(list(df = NULL, adj_col = NULL))
    }

    # Validate pairing completeness if pairing is enabled and a group variable is selected
    if (!is.null(pairing_var_run) && nzchar(pairing_var_run) &&
      !is.null(group_var_run) && nzchar(group_var_run)) {
      if (!pairing_var_run %in% colnames(rv$meta_sample)) {
        showNotification("Pairing variable not found in metadata.", type = "error")
        return(list(df = NULL, adj_col = NULL))
      }

      if (!group_var_run %in% colnames(rv$meta_sample)) {
        showNotification("Group variable not found in metadata.", type = "error")
        return(list(df = NULL, adj_col = NULL))
      }

      # First check if pairing is even possible (any pair values appear in multiple groups)
      pair_group_check <- rv$meta_sample %>%
        dplyr::group_by(!!rlang::sym(pairing_var_run)) %>%
        dplyr::summarise(
          n_groups = dplyr::n_distinct(!!rlang::sym(group_var_run)),
          .groups = "drop"
        )

      # Only validate completeness if pairing is actually possible
      if (any(pair_group_check$n_groups >= 2)) {
        # Pairing is possible - check for complete pairs
        pair_check <- rv$meta_sample %>%
          dplyr::group_by(!!rlang::sym(pairing_var_run)) %>%
          dplyr::summarise(
            n_groups = dplyr::n_distinct(!!rlang::sym(group_var_run)),
            groups = paste(sort(unique(as.character(!!rlang::sym(group_var_run)))), collapse = ", "),
            .groups = "drop"
          )

        total_groups <- dplyr::n_distinct(rv$meta_sample[[group_var_run]])
        incomplete_pairs <- pair_check %>%
          dplyr::filter(n_groups < total_groups)

        if (nrow(incomplete_pairs) > 0) {
          # Build detailed error message
          incomplete_details <- head(incomplete_pairs, 10)
          details_text <- paste(
            sprintf(
              "    %s: %d group(s) present [%s]",
              incomplete_details[[pairing_var_run]],
              incomplete_details$n_groups,
              incomplete_details$groups
            ),
            collapse = "\n"
          )

          error_msg <- sprintf(
            "ERROR: Pairing Validation Failed\n\nCannot run paired test with incomplete pairs.\n\nProblem:\n  - %d subject(s) have incomplete pairs (missing samples in some groups)\n  - Each paired subject must have samples in ALL %d groups\n\nIncomplete pairs detected:\n%s%s\n\nSolution:\n  - Update subsetting rules to ensure complete pairs, OR\n  - Disable pairing (set Pairing Variable to 'None'), OR\n  - Change Group Variable to match available data",
            nrow(incomplete_pairs),
            total_groups,
            details_text,
            if (nrow(incomplete_pairs) > 10) sprintf("\n    ... and %d more", nrow(incomplete_pairs) - 10) else ""
          )

          return(list(df = NULL, adj_col = NULL, error = error_msg))
        }
      }
      # If no pairing is possible (n_groups always 1), fall back to unpaired tests silently
    }

    # Aggregate to celltypes if requested
    abund <- abund0
    if (test_entity_run == "Celltypes" && !is.null(rv$cluster_map)) {
      cm <- rv$cluster_map
      keep <- cm$cluster %in% colnames(abund)
      cm <- cm[keep, , drop = FALSE]
      if (!nrow(cm)) {
        showNotification("No overlapping clusters to aggregate into celltypes.", type = "error")
        return(list(df = NULL, adj_col = NULL))
      }
      split_idx <- split(cm$cluster, cm$celltype)
      abund <- sapply(split_idx, function(cols) rowSums(abund0[, cols, drop = FALSE]))
      abund <- as.matrix(abund)
    }

    # Merge with per-sample metadata
    meta_sub <- rv$meta_sample %>% dplyr::select(patient_ID, dplyr::everything())
    abund_df <- as.data.frame(abund)
    abund_df$patient_ID <- stringr::str_extract(string = rownames(abund_df), pattern = paste0("(", paste0(meta_sub$patient_ID, collapse = "|"), ")"))
    merged <- merge(x = meta_sub, y = abund_df, by = "patient_ID")
    # Long format
    abund_long <- merged %>%
      tidyr::pivot_longer(cols = colnames(abund), names_to = "entity", values_to = "freq")

    # Determine if we're using paired tests - check both if pairing is enabled AND feasible
    use_pairing <- FALSE
    if (!is.null(pairing_var_run) && nzchar(pairing_var_run) && pairing_var_run %in% colnames(abund_long) &&
      !is.null(group_var_run) && nzchar(group_var_run) && group_var_run %in% colnames(abund_long)) {
      # Check if any pairing values appear in multiple groups (indicates pairing is possible)
      pair_group_counts <- abund_long %>%
        dplyr::group_by(!!rlang::sym(pairing_var_run)) %>%
        dplyr::summarise(n_groups = dplyr::n_distinct(!!rlang::sym(group_var_run)), .groups = "drop")
      use_pairing <- any(pair_group_counts$n_groups >= 2)
    }

    # Run tests per entity
    res <- abund_long %>%
      dplyr::group_by(entity) %>%
      dplyr::group_modify(~ {
        if (test_type_run == "Wilcoxon (2-group)") {
          if (!nzchar(group_var_run)) {
            return(data.frame(test = "wilcox", p = NA, n = nrow(.x)))
          }

          g_raw <- .x[[group_var_run]]
          ok <- !is.na(g_raw)
          g <- droplevels(factor(g_raw[ok]))
          freq_ok <- .x$freq[ok]

          if (length(levels(g)) != 2) {
            return(data.frame(test = "wilcox", p = NA, n = sum(ok)))
          }

          if (use_pairing) {
            # Paired Wilcoxon test
            pair_var <- .x[[pairing_var_run]][ok]
            df_test <- data.frame(freq = freq_ok, group = g, pair = pair_var)
            df_test <- df_test[complete.cases(df_test), ]

            # Ensure each pair has exactly 2 observations
            pair_counts <- table(df_test$pair)
            valid_pairs <- names(pair_counts)[pair_counts == 2]
            df_test <- df_test[df_test$pair %in% valid_pairs, ]

            if (nrow(df_test) < 4) {
              return(data.frame(test = "wilcox_paired", p = NA, n = nrow(df_test)))
            }

            # Reshape to wide format for paired test
            df_wide <- tidyr::pivot_wider(df_test, names_from = group, values_from = freq, id_cols = pair)
            if (ncol(df_wide) != 3) {
              return(data.frame(test = "wilcox_paired", p = NA, n = nrow(df_test)))
            }

            x1 <- df_wide[[2]]
            x2 <- df_wide[[3]]
            wt <- suppressWarnings(wilcox.test(x1, x2, paired = TRUE))
            data.frame(test = "wilcox_paired", n = length(valid_pairs), p = wt$p.value)
          } else {
            # Unpaired Wilcoxon test
            wt <- suppressWarnings(wilcox.test(freq_ok ~ g))
            data.frame(test = "wilcox", n = sum(ok), p = wt$p.value)
          }
        } else if (test_type_run == "Pairwise Wilcoxon") {
          if (!nzchar(group_var_run)) {
            return(data.frame(test = "pairwise_wilcox", comparison = NA, p = NA, n = NA))
          }

          g_raw <- .x[[group_var_run]]
          ok <- !is.na(g_raw)
          g <- droplevels(factor(g_raw[ok]))
          freq_ok <- .x$freq[ok]

          if (length(levels(g)) < 2) {
            return(data.frame(test = "pairwise_wilcox", comparison = NA, p = NA, n = sum(ok)))
          }

          grp_levels <- levels(g)
          comparisons <- combn(grp_levels, 2, simplify = FALSE)

          pairwise_results <- lapply(comparisons, function(pair) {
            grp1 <- pair[1]
            grp2 <- pair[2]
            idx1 <- which(g == grp1)
            idx2 <- which(g == grp2)

            if (use_pairing) {
              # Paired pairwise test
              pair_var <- .x[[pairing_var_run]][ok]
              df_test <- data.frame(
                freq = freq_ok[c(idx1, idx2)],
                group = g[c(idx1, idx2)],
                pair = pair_var[c(idx1, idx2)]
              )
              df_test <- df_test[complete.cases(df_test), ]

              pair_counts <- table(df_test$pair)
              valid_pairs <- names(pair_counts)[pair_counts == 2]
              df_test <- df_test[df_test$pair %in% valid_pairs, ]

              if (nrow(df_test) < 4) {
                return(data.frame(
                  test = "pairwise_wilcox_paired",
                  comparison = paste(grp1, "vs", grp2),
                  p = NA,
                  n = nrow(df_test)
                ))
              }

              df_wide <- tidyr::pivot_wider(df_test, names_from = group, values_from = freq, id_cols = pair)
              if (ncol(df_wide) != 3) {
                return(data.frame(
                  test = "pairwise_wilcox_paired",
                  comparison = paste(grp1, "vs", grp2),
                  p = NA,
                  n = nrow(df_test)
                ))
              }

              x1 <- df_wide[[2]]
              x2 <- df_wide[[3]]
              wt <- suppressWarnings(wilcox.test(x1, x2, paired = TRUE))
              data.frame(
                test = "pairwise_wilcox_paired",
                comparison = paste(grp1, "vs", grp2),
                p = wt$p.value,
                n = length(valid_pairs)
              )
            } else {
              # Unpaired pairwise test
              freq1 <- freq_ok[idx1]
              freq2 <- freq_ok[idx2]
              wt <- suppressWarnings(wilcox.test(freq1, freq2))
              data.frame(
                test = "pairwise_wilcox",
                comparison = paste(grp1, "vs", grp2),
                p = wt$p.value,
                n = length(idx1) + length(idx2)
              )
            }
          })

          do.call(rbind, pairwise_results)
        } else if (test_type_run == "Friedman (multi-group paired)") {
          # Friedman test for paired multi-group comparisons
          if (!use_pairing || !nzchar(group_var_run)) {
            return(data.frame(test = "friedman", p = NA, n = NA))
          }

          g_raw <- .x[[group_var_run]]
          ok <- !is.na(g_raw)
          g <- droplevels(factor(g_raw[ok]))
          freq_ok <- .x$freq[ok]

          if (length(levels(g)) < 3) {
            return(data.frame(test = "friedman", p = NA, n = sum(ok)))
          }

          pair_var <- .x[[pairing_var_run]][ok]
          df_test <- data.frame(freq = freq_ok, group = g, pair = pair_var)
          df_test <- df_test[complete.cases(df_test), ]

          # Check that each pair has observations in all groups
          pair_group_counts <- table(df_test$pair, df_test$group)
          # Keep only pairs that have exactly 1 observation per group
          valid_pairs <- rownames(pair_group_counts)[apply(pair_group_counts, 1, function(x) all(x == 1))]

          if (length(valid_pairs) < 3) {
            return(data.frame(test = "friedman", p = NA, n = length(valid_pairs)))
          }

          df_test <- df_test[df_test$pair %in% valid_pairs, ]

          # Reshape to wide format for Friedman test
          df_wide <- tidyr::pivot_wider(df_test, names_from = group, values_from = freq, id_cols = pair)
          if (ncol(df_wide) < 4) {
            return(data.frame(test = "friedman", p = NA, n = length(valid_pairs)))
          }

          # Friedman test requires a matrix with subjects as rows and treatments as columns
          mat <- as.matrix(df_wide[, -1]) # Remove pair column

          ft <- suppressWarnings(friedman.test(mat))
          data.frame(test = "friedman", n = length(valid_pairs), p = ft$p.value)
        } else if (test_type_run == "Kruskal–Wallis (multi-group unpaired)") {
          if (!nzchar(group_var_run)) {
            return(data.frame(test = "kruskal", p = NA, n = nrow(.x)))
          }

          g_raw <- .x[[group_var_run]]
          ok <- !is.na(g_raw)
          g <- droplevels(factor(g_raw[ok]))
          freq_ok <- .x$freq[ok]

          if (length(levels(g)) < 2) {
            return(data.frame(test = "kruskal", p = NA, n = sum(ok)))
          }

          kw <- suppressWarnings(kruskal.test(freq_ok ~ g))
          data.frame(test = "kruskal", n = sum(ok), p = kw$p.value)
        } else if (test_type_run == "Spearman (continuous)") {
          if (!nzchar(cont_var_run)) {
            return(data.frame(test = "spearman", rho = NA, p = NA, n = nrow(.x)))
          }

          x <- .x$freq
          y <- suppressWarnings(as.numeric(.x[[cont_var_run]]))

          ok <- complete.cases(x, y)
          if (sum(ok) < 3) {
            return(data.frame(test = "spearman", rho = NA, p = NA, n = sum(ok)))
          }

          ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
          data.frame(test = "spearman", rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
        }
      }) %>%
      dplyr::ungroup()

    # Adjust p-values
    adj_col <- NULL
    if (nrow(res) && "p" %in% names(res) && nzchar(p_adj_method_run)) {
      adj_col <- paste0(tolower(p_adj_method_run), "_padj")
      res[[adj_col]] <- p.adjust(res$p, method = p_adj_method_run)
    }

    list(df = res, adj_col = adj_col)
  })

  observeEvent(input$run_test, {
    res <- run_tests()
    test_results_rv(res)

    # Store test info for export filename
    rv$last_test_info <- list(
      entity = input$test_entity %||% "clusters",
      test = input$test_type %||% "unknown",
      metadata = if (grepl("spearman|pearson", tolower(input$test_type %||% ""))) {
        input$cont_var %||% "variable"
      } else {
        input$group_var %||% "variable"
      }
    )

    output$test_cleared_msg <- renderText(NULL)
  })

  observeEvent(input$reset_test, {
    test_results_rv(NULL) # clear
    showNotification("Testing results cleared.", type = "message", duration = 5)
    output$test_cleared_msg <- renderText("Results cleared. Run a new test to see results here.")
  })

  output$hasResults <- reactive({
    run <- test_results_rv()
    !is.null(run) && (!is.null(run$error) || (!is.null(run$df) && nrow(run$df) > 0))
  })
  outputOptions(output, "hasResults", suspendWhenHidden = FALSE)

  output$hasValidResults <- reactive({
    run <- test_results_rv()
    !is.null(run) && is.null(run$error) && !is.null(run$df) && nrow(run$df) > 0
  })
  outputOptions(output, "hasValidResults", suspendWhenHidden = FALSE)

  output$test_error_msg <- renderUI({
    run <- test_results_rv()
    if (!is.null(run$error)) {
      tags$div(
        style = "background-color: #2b2b2b; color: #ff6b6b; border: 2px solid #ff6b6b; border-radius: 4px; padding: 15px; margin-bottom: 20px; font-family: 'Courier New', monospace; white-space: pre-wrap; font-size: 13px;",
        tags$pre(
          style = "margin: 0; color: #ff6b6b;",
          run$error
        )
      )
    } else {
      NULL
    }
  })

  output$test_table <- renderTable(
    {
      run <- test_results_rv()
      if (!is.null(run$error)) {
        return(NULL)
      }
      df <- req(run$df)
      adj_col <- run$adj_col

      # sort while numeric
      if (!is.null(adj_col) && adj_col %in% names(df) && is.numeric(df[[adj_col]])) {
        df <- df[order(df[[adj_col]], na.last = TRUE), ]
      }

      # format for display
      df_display <- df
      num_cols <- intersect(c("p", adj_col), names(df_display))
      for (col in num_cols) {
        if (is.numeric(df_display[[col]])) {
          # Use scientific notation for very small p-values
          df_display[[col]] <- sapply(df_display[[col]], function(x) {
            if (is.na(x)) {
              return(NA_character_)
            } else if (x < 0.001) {
              formatC(x, format = "e", digits = 2)
            } else {
              formatC(x, format = "f", digits = 3)
            }
          })
        }
      }
      if ("rho" %in% names(df_display) && is.numeric(df_display$rho)) {
        df_display$rho <- formatC(df_display$rho, format = "f", digits = 2)
      }

      df_display
    },
    sanitize.text.function = function(x) x
  )

  output$export_results <- downloadHandler(
    filename = function() {
      info <- rv$last_test_info
      subset_id <- rv$subset_id %||% "000000000"
      if (is.null(info)) {
        return(paste0("results_", subset_id, ".csv"))
      }

      # Determine if continuous or categorical based on test type
      test_type <- if (grepl("spearman|pearson", tolower(info$test))) {
        "continuous"
      } else {
        "categorical"
      }

      # Format test name (lowercase, underscores, no parentheses)
      test <- tolower(info$test)
      test <- gsub("\\s+", "_", test)
      test <- gsub("_\\(.*\\)", "", test)
      test <- gsub("kruskal_wallis", "kruskal-wallis", test)

      # Format entity and metadata (clusters/outcome)
      entity <- gsub("\\s+", "_", tolower(info$entity))
      metadata <- gsub("\\s+", "_", tolower(info$metadata))

      paste0(test_type, "_", entity, "_", metadata, "_", test, "_", subset_id, ".csv")
    },
    content = function(file) {
      run <- test_results_rv()
      req(!is.null(run), is.null(run$error))
      df <- run$df
      req(!is.null(df), nrow(df) > 0)
      write.csv(df, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )

  # Update cat_group_var choices when metadata arrives
  observeEvent(rv$meta_cell,
    {
      meta_cols <- colnames(rv$meta_cell)
      categorical_choices <- sort(meta_cols[sapply(rv$meta_cell, function(x) is.character(x) || is.factor(x))])
      # Exclude cluster, patient_ID, run_date, and source from categorical outcome choices
      categorical_choices <- setdiff(categorical_choices, c("cluster", "patient_ID", "run_date", "source"))
      # Apply global settings filter
      categorical_choices <- filter_by_global_settings(categorical_choices)
      updatePickerInput(session, "cat_group_var",
        choices = c("", categorical_choices),
        selected = ""
      )
    },
    ignoreInit = TRUE
  )

  # If you also need group_var/cont_var (Testing tab) — keep the same pattern:
  observeEvent(rv$meta_cell,
    {
      meta_cols <- colnames(rv$meta_cell)
      categorical_choices <- sort(meta_cols[sapply(rv$meta_cell, function(x) is.character(x) || is.factor(x))])
      # Exclude cluster, patient_ID, run_date, and source from categorical outcome choices
      categorical_choices <- setdiff(categorical_choices, c("cluster", "patient_ID", "run_date", "source"))
      # Apply global settings filter
      categorical_choices <- filter_by_global_settings(categorical_choices)
      continuous_choices <- sort(meta_cols[sapply(rv$meta_cell, function(x) is.numeric(x) || is.integer(x))])
      # Apply global settings filter
      continuous_choices <- filter_by_global_settings(continuous_choices)
      updatePickerInput(session, "group_var", choices = c("", categorical_choices), selected = "")
      updatePickerInput(session, "cont_var", choices = c("", continuous_choices), selected = "")
    },
    ignoreInit = TRUE
  )

  # Continuous metadata picker for “Continuous” tab
  observeEvent(rv$meta_cell,
    {
      meta_cols <- colnames(rv$meta_cell)
      continuous_choices <- sort(meta_cols[sapply(rv$meta_cell, is.numeric)])
      # Apply global settings filter
      continuous_choices <- filter_by_global_settings(continuous_choices)
      updatePickerInput(session, "cont_group_var",
        choices = c("", continuous_choices),
        selected = ""
      )
    },
    ignoreInit = TRUE
  )

  # Helper to make safe input IDs from arbitrary group values
  sanitize_id <- function(x) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9_\\-]", "_", x)
    x
  }

  # Populate color pickers for the currently-selected categorical grouping
  # This version will use the currently-selected metadata column (input$cat_group_var)
  # even if the user has not yet clicked Generate plots. NA is excluded from the choices.
  observeEvent(input$cat_populate_colors, {
    req(rv$meta_cell)

    gv <- input$cat_group_var
    if (is.null(gv) || !nzchar(gv) || !(gv %in% colnames(rv$meta_cell))) {
      showNotification("No valid grouping variable selected for colors.", type = "error")
      output$cat_color_pickers_ui <- renderUI(NULL)
      rv$cat_colors <- NULL
      return()
    }

    # Always pull levels from the currently selected metadata column
    levels_vec <- sort(unique(as.character(rv$meta_cell[[gv]])))
    levels_vec <- levels_vec[!is.na(levels_vec)]

    if (length(levels_vec) == 0) {
      showNotification("No non-missing group levels found to populate colors for.", type = "warning")
      output$cat_color_pickers_ui <- renderUI(NULL)
      rv$cat_colors <- NULL
      return()
    }

    # Default palette
    default_pal <- viridis::viridis(length(levels_vec))
    names(default_pal) <- levels_vec

    # Build UI inputs
    ui_list <- lapply(seq_along(levels_vec), function(i) {
      lv <- levels_vec[i]
      input_id <- paste0("cat_color_", sanitize_id(lv))
      colourpicker::colourInput(
        inputId = input_id,
        label = paste0("Color for ", lv),
        value = default_pal[i],
        showColour = "both"
      )
    })

    output$cat_color_pickers_ui <- renderUI({
      tagList(
        tags$div(style = "max-height: 300px; overflow-y: auto; padding-right: 6px;", ui_list)
      )
    })

    # Initialize rv$cat_colors
    rv$cat_colors <- setNames(as.character(default_pal), levels_vec)
  })

  # Observe any manual color changes and persist into rv$cat_colors
  observe({
    # If rv$cat_colors exists, watch its keys and update from inputs
    req(!is.null(rv$cat_colors))
    new_colors <- rv$cat_colors
    for (grp in names(rv$cat_colors)) {
      input_name <- paste0("cat_color_", sanitize_id(grp))
      if (!is.null(input[[input_name]])) {
        new_colors[[grp]] <- input[[input_name]]
      }
    }
    rv$cat_colors <- new_colors
  })

  cat_plot_data <- eventReactive(input$generate_cat_plots, {
    req(rv$meta_sample, rv$abundance_sample)

    abund0 <- rv$abundance_sample
    if (is.null(abund0)) {
      showNotification("No abundance matrix available.", type = "error")
      return(NULL)
    }

    # Get plotting parameters early for validation
    group_var <- input$cat_group_var
    pairing_var <- input$pairing_var

    # Validate pairing completeness if pairing is enabled and a group variable is selected
    if (!is.null(pairing_var) && nzchar(pairing_var) &&
      !is.null(group_var) && nzchar(group_var)) {
      if (!pairing_var %in% colnames(rv$meta_sample)) {
        error_msg <- "ERROR: Pairing Validation Failed\n\nPairing variable not found in metadata."
        return(list(error = error_msg))
      }

      if (!group_var %in% colnames(rv$meta_sample)) {
        error_msg <- "ERROR: Group Validation Failed\n\nGroup variable not found in metadata."
        return(list(error = error_msg))
      }

      # First check if pairing is even possible (any pair values appear in multiple groups)
      pair_group_check <- rv$meta_sample %>%
        dplyr::group_by(!!rlang::sym(pairing_var)) %>%
        dplyr::summarise(
          n_groups = dplyr::n_distinct(!!rlang::sym(group_var)),
          .groups = "drop"
        )

      # Only validate completeness if pairing is actually possible
      if (any(pair_group_check$n_groups >= 2)) {
        # Pairing is possible - check for complete pairs
        pair_check <- rv$meta_sample %>%
          dplyr::group_by(!!rlang::sym(pairing_var)) %>%
          dplyr::summarise(
            n_groups = dplyr::n_distinct(!!rlang::sym(group_var)),
            groups = paste(sort(unique(as.character(!!rlang::sym(group_var)))), collapse = ", "),
            .groups = "drop"
          )

        total_groups <- dplyr::n_distinct(rv$meta_sample[[group_var]])
        incomplete_pairs <- pair_check %>%
          dplyr::filter(n_groups < total_groups)

        if (nrow(incomplete_pairs) > 0) {
          # Build detailed error message with proper line breaks for HTML
          incomplete_details <- head(incomplete_pairs, 10)
          details_text <- paste(
            sprintf(
              "    %s: %d group(s) present [%s]",
              incomplete_details[[pairing_var]],
              incomplete_details$n_groups,
              incomplete_details$groups
            ),
            collapse = "\n"
          )

          error_msg <- sprintf(
            "ERROR: Pairing Validation Failed\n\nCannot generate paired plots with incomplete pairs.\n\nProblem:\n  - %d subject(s) have incomplete pairs (missing samples in some groups)\n  - Each paired subject must have samples in ALL %d groups\n\nIncomplete pairs detected:\n%s%s\n\nSolution:\n  - Update subsetting rules to ensure complete pairs, OR\n  - Disable pairing (set Pairing Variable to 'None'), OR\n  - Change Group Variable to match available data",
            nrow(incomplete_pairs),
            total_groups,
            details_text,
            if (nrow(incomplete_pairs) > 10) sprintf("\n    ... and %d more", nrow(incomplete_pairs) - 10) else ""
          )

          return(list(error = error_msg))
        }
      }
      # If no pairing is possible (n_groups always 1), fall back to unpaired plots silently
    }

    # Aggregate to celltypes if needed
    abund <- abund0
    if (input$cat_entity == "Celltypes" && !is.null(rv$cluster_map)) {
      cm <- rv$cluster_map
      keep <- cm$cluster %in% colnames(abund)
      cm <- cm[keep, , drop = FALSE]
      if (!nrow(cm)) {
        return(NULL)
      }
      split_idx <- split(cm$cluster, cm$celltype)
      abund <- sapply(split_idx, function(cols) rowSums(abund0[, cols, drop = FALSE]))
      abund <- as.matrix(abund)
    }

    # Merge with per-sample metadata
    meta_sub <- rv$meta_sample %>% dplyr::select(patient_ID, dplyr::everything())
    abund_df <- as.data.frame(abund, check.names = FALSE, stringsAsFactors = FALSE)
    abund_df$patient_ID <- stringr::str_extract(string = rownames(abund_df), pattern = paste0("(", paste0(meta_sub$patient_ID, collapse = "|"), ")"))
    # merged <- dplyr::left_join(meta_sub, abund_df, by = "patient_ID")
    merged <- merge(x = meta_sub, y = abund_df, by = "patient_ID")

    # Long format and clean
    abund_long <- merged %>%
      tidyr::pivot_longer(cols = colnames(abund), names_to = "entity", values_to = "freq")

    # Filter NA frequencies first (before any other operations)
    abund_long <- abund_long[!is.na(abund_long$freq), ]

    # Clean entity names once (more efficient than gsub on every row)
    if (any(grepl("\n", abund_long$entity, fixed = TRUE))) {
      entity_levels <- unique(abund_long$entity)
      entity_clean <- gsub("\n", " ", entity_levels, fixed = TRUE)
      abund_long$entity <- factor(abund_long$entity, levels = entity_levels, labels = entity_clean)
    }

    # Capture plot-type and point-mode at Generate time
    plot_type_selected <- input$cat_plot_type %||% "box"
    point_mode_selected <- input$cat_points %||% "draw"

    # Run test per entity
    test_type <- input$cat_test_type
    group_var <- input$cat_group_var
    pairing_var <- input$pairing_var

    # Check if pairing is both enabled AND feasible for this categorical variable
    use_pairing <- FALSE
    if (!is.null(pairing_var) && nzchar(pairing_var) && pairing_var %in% colnames(abund_long) &&
      !is.null(group_var) && nzchar(group_var) && group_var %in% colnames(abund_long)) {
      # Check if any pairing values appear in multiple groups (indicates pairing is possible)
      pair_group_counts <- abund_long %>%
        dplyr::group_by(!!rlang::sym(pairing_var)) %>%
        dplyr::summarise(n_groups = dplyr::n_distinct(!!rlang::sym(group_var)), .groups = "drop")
      use_pairing <- any(pair_group_counts$n_groups >= 2)
    }

    if (is.null(group_var) || !nzchar(group_var)) {
      res <- abund_long %>%
        dplyr::group_by(entity) %>%
        dplyr::summarise(p = NA_real_, .groups = "drop")
    } else {
      res <- abund_long %>%
        dplyr::group_by(entity) %>%
        dplyr::group_modify(~ {
          g_raw <- .x[[group_var]]
          ok_rows <- !is.na(g_raw)
          if (!any(ok_rows)) {
            return(data.frame(p = NA_real_))
          }
          g <- droplevels(factor(g_raw[ok_rows]))
          freq_ok <- .x$freq[ok_rows]
          if (length(unique(g)) < 2) {
            return(data.frame(p = NA_real_))
          }

          if (test_type == "Wilcoxon (2-group)" && length(unique(g)) == 2) {
            if (use_pairing) {
              # Paired Wilcoxon
              pair_var <- .x[[pairing_var]][ok_rows]
              df_test <- data.frame(freq = freq_ok, group = g, pair = pair_var)
              df_test <- df_test[complete.cases(df_test), ]

              pair_counts <- table(df_test$pair)
              valid_pairs <- names(pair_counts)[pair_counts == 2]
              df_test <- df_test[df_test$pair %in% valid_pairs, ]

              if (nrow(df_test) < 4) {
                return(data.frame(p = NA_real_))
              }

              df_wide <- tidyr::pivot_wider(df_test, names_from = group, values_from = freq, id_cols = pair)
              if (ncol(df_wide) != 3) {
                return(data.frame(p = NA_real_))
              }

              x1 <- df_wide[[2]]
              x2 <- df_wide[[3]]
              wt <- suppressWarnings(wilcox.test(x1, x2, paired = TRUE))
              data.frame(p = wt$p.value)
            } else {
              # Unpaired Wilcoxon
              wt <- suppressWarnings(wilcox.test(freq_ok ~ g))
              data.frame(p = wt$p.value)
            }
          } else if (test_type == "Pairwise Wilcoxon") {
            # Pairwise comparisons for all group combinations
            grp_levels <- levels(g)
            if (length(grp_levels) < 2) {
              return(data.frame(comparison = NA_character_, p = NA_real_))
            }

            comparisons <- combn(grp_levels, 2, simplify = FALSE)

            pairwise_results <- lapply(comparisons, function(pair) {
              grp1 <- pair[1]
              grp2 <- pair[2]
              idx1 <- which(g == grp1)
              idx2 <- which(g == grp2)

              if (use_pairing) {
                # Paired pairwise test
                pair_var <- .x[[pairing_var]][ok_rows]
                df_test <- data.frame(
                  freq = freq_ok[c(idx1, idx2)],
                  group = g[c(idx1, idx2)],
                  pair = pair_var[c(idx1, idx2)]
                )
                df_test <- df_test[complete.cases(df_test), ]

                pair_counts <- table(df_test$pair)
                valid_pairs <- names(pair_counts)[pair_counts == 2]
                df_test <- df_test[df_test$pair %in% valid_pairs, ]

                if (nrow(df_test) < 4) {
                  return(data.frame(
                    comparison = paste(grp1, "vs", grp2),
                    p = NA_real_
                  ))
                }

                df_wide <- tidyr::pivot_wider(df_test, names_from = group, values_from = freq, id_cols = pair)
                if (ncol(df_wide) != 3) {
                  return(data.frame(
                    comparison = paste(grp1, "vs", grp2),
                    p = NA_real_
                  ))
                }

                x1 <- df_wide[[2]]
                x2 <- df_wide[[3]]
                wt <- suppressWarnings(wilcox.test(x1, x2, paired = TRUE))
                data.frame(
                  comparison = paste(grp1, "vs", grp2),
                  p = wt$p.value
                )
              } else {
                # Unpaired pairwise test
                freq1 <- freq_ok[idx1]
                freq2 <- freq_ok[idx2]
                wt <- suppressWarnings(wilcox.test(freq1, freq2))
                data.frame(
                  comparison = paste(grp1, "vs", grp2),
                  p = wt$p.value
                )
              }
            })

            # Return all pairwise comparisons
            do.call(rbind, pairwise_results)
          } else if (test_type == "Friedman (multi-group paired)") {
            # Friedman test for paired multi-group comparisons
            if (!use_pairing) {
              return(data.frame(p = NA_real_))
            }

            grp_levels <- levels(g)
            if (length(grp_levels) < 3) {
              return(data.frame(p = NA_real_))
            }

            pair_var <- .x[[pairing_var]][ok_rows]
            df_test <- data.frame(freq = freq_ok, group = g, pair = pair_var)
            df_test <- df_test[complete.cases(df_test), ]

            # Check that each pair has observations in all groups
            pair_group_counts <- table(df_test$pair, df_test$group)
            # Keep only pairs that have exactly 1 observation per group
            valid_pairs <- rownames(pair_group_counts)[apply(pair_group_counts, 1, function(x) all(x == 1))]

            if (length(valid_pairs) < 3) {
              return(data.frame(p = NA_real_))
            }

            df_test <- df_test[df_test$pair %in% valid_pairs, ]

            # Reshape to wide format for Friedman test
            df_wide <- tidyr::pivot_wider(df_test, names_from = group, values_from = freq, id_cols = pair)
            if (ncol(df_wide) < 4) {
              return(data.frame(p = NA_real_))
            } # Need pair column + at least 3 group columns

            # Friedman test requires a matrix with subjects as rows and treatments as columns
            mat <- as.matrix(df_wide[, -1]) # Remove pair column

            ft <- suppressWarnings(friedman.test(mat))
            data.frame(p = ft$p.value)
          } else if (test_type == "Kruskal–Wallis (multi-group unpaired)") {
            kw <- kruskal.test(freq_ok ~ as.factor(g))
            data.frame(p = kw$p.value)
          } else {
            data.frame(p = NA_real_)
          }
        }) %>%
        dplyr::ungroup()
    }

    # Adjust p-values if requested
    if (nrow(res) && nzchar(input$cat_p_adj_method) && "p" %in% names(res)) {
      res$padj <- p.adjust(res$p, method = input$cat_p_adj_method)
    }

    # Save info for export
    rv$last_cat_info <- list(
      entity = tolower(input$cat_entity %||% "clusters"),
      group = tolower(input$cat_group_var %||% "group"),
      test_raw = input$cat_test_type %||% "test"
    )

    # Return everything needed for plotting
    list(
      data = abund_long,
      results = res,
      group_var = input$cat_group_var,
      test_type = test_type,
      pairing_var = pairing_var,
      use_pairing = use_pairing, # Whether pairing is feasible for this categorical variable
      use_adj_p = input$cat_use_adj_p,
      facet_cols = as.numeric(input$cat_max_facets),
      plot_type = plot_type_selected,
      point_mode = point_mode_selected,
      colors = if (!is.null(rv$cat_colors)) rv$cat_colors else NULL
    )
  })

  observeEvent(input$generate_cat_plots, {
    cp <- cat_plot_data() # your existing eventReactive
    # Capture current colors at the moment Generate is clicked
    if (!is.null(rv$cat_colors)) {
      cp$colors <- rv$cat_colors
    }
    cat_state(cp)
    output$cat_cleared_msg <- renderText(NULL)
  })

  observeEvent(input$reset_cat, {
    cat_state(NULL) # clear the state
    cat_plot_cache(NULL) # also clear cached ggplot
    showNotification("Categorical plots cleared.", type = "message", duration = 5)
    output$cat_cleared_msg <- renderText("Results cleared. Generate new plots to see them here.")
  })

  output$hasCatResults <- reactive({
    cp <- cat_state()
    !is.null(cp) && (!is.null(cp$error) || (!is.null(cp$data) && nrow(cp$data) > 0))
  })
  outputOptions(output, "hasCatResults", suspendWhenHidden = FALSE)

  output$hasValidCatResults <- reactive({
    cp <- cat_state()
    !is.null(cp) && is.null(cp$error) && !is.null(cp$data) && nrow(cp$data) > 0
  })
  outputOptions(output, "hasValidCatResults", suspendWhenHidden = FALSE)

  # Populate color pickers when button is clicked
  observeEvent(input$cat_populate_colors, {
    req(rv$meta_sample)

    levels_vec <- NULL
    cp <- NULL

    # Try to use the last generated cat_plot_data (ensures exact plotting groups)
    try(
      {
        cp <- cat_plot_data()
      },
      silent = TRUE
    )
    if (!is.null(cp) && is.list(cp) && !is.null(cp$data) && nrow(cp$data) > 0) {
      group_col <- cp$group_var
      if (!is.null(group_col) && nzchar(group_col) && group_col %in% colnames(cp$data)) {
        levels_vec <- unique(as.character(cp$data[[group_col]]))
      }
    }

    # Fallback: use metadata column directly
    if (is.null(levels_vec) || length(levels_vec) == 0) {
      gv <- input$cat_group_var
      if (!is.null(gv) && nzchar(gv) && gv %in% colnames(rv$meta_sample)) {
        levels_vec <- sort(unique(as.character(rv$meta_sample[[gv]])))
      }
    }

    levels_vec <- levels_vec[!is.na(levels_vec)]
    if (length(levels_vec) == 0) {
      showNotification("No non-missing group levels found to populate colors for.", type = "warning")
      output$cat_color_pickers_ui <- renderUI(NULL)
      rv$cat_colors <- NULL
      return()
    }

    # Default palette
    default_pal <- viridis::viridis(length(levels_vec))
    names(default_pal) <- levels_vec

    # Build UI inputs
    ui_list <- lapply(seq_along(levels_vec), function(i) {
      lv <- levels_vec[i]
      input_id <- paste0("cat_color_", sanitize_id(lv))
      colourpicker::colourInput(
        inputId = input_id,
        label = paste0("Color for ", lv),
        value = default_pal[i],
        showColour = "both"
      )
    })

    output$cat_color_pickers_ui <- renderUI({
      tagList(
        tags$div(style = "max-height: 300px; overflow-y: auto; padding-right: 6px;", ui_list)
      )
    })

    # Initialize rv$cat_colors
    rv$cat_colors <- setNames(as.character(default_pal), levels_vec)
  })

  # Observe manual color changes
  observe({
    req(!is.null(rv$cat_colors))
    new_colors <- rv$cat_colors
    for (grp in names(rv$cat_colors)) {
      input_name <- paste0("cat_color_", sanitize_id(grp))
      if (!is.null(input[[input_name]])) {
        new_colors[[grp]] <- input[[input_name]]
      }
    }
    rv$cat_colors <- new_colors
  })

  # Reset colors when group_var changes
  observeEvent(input$cat_group_var, {
    output$cat_color_pickers_ui <- renderUI(NULL)
    rv$cat_colors <- NULL
  })


  output$cat_error_msg <- renderUI({
    cp <- cat_state()
    if (!is.null(cp$error)) {
      tags$div(
        style = "background-color: #2b2b2b; color: #ff6b6b; border: 2px solid #ff6b6b; border-radius: 4px; padding: 15px; margin-bottom: 20px; font-family: 'Courier New', monospace; white-space: pre-wrap; font-size: 13px;",
        tags$pre(
          style = "margin: 0; color: #ff6b6b;",
          cp$error
        )
      )
    } else {
      NULL
    }
  })

  output$categorical_plot <- renderPlot(
    {
      # cp <- cat_plot_data(); req(cp)
      cp <- cat_state()
      req(cp)
      if (!is.null(cp$error)) {
        return(NULL)
      }
      abund_long <- cp$data
      res <- cp$results
      group_var <- cp$group_var
      test_type <- cp$test_type
      use_adj_p <- cp$use_adj_p
      facet_cols <- cp$facet_cols
      plot_type <- cp$plot_type %||% input$cat_plot_type %||% "box"
      point_mode <- cp$point_mode %||% input$cat_points %||% "draw"

      # Validate grouping variable
      if (is.null(group_var) || !nzchar(group_var) || !(group_var %in% colnames(abund_long))) {
        showNotification("No valid grouping variable selected for categorical plotting.", type = "error")
        return(invisible(NULL))
      }

      # Remove rows with NA grouping values for plotting
      abund_long_plot <- abund_long %>% dplyr::filter(!is.na(.data[[group_var]]))

      if (nrow(abund_long_plot) == 0) {
        showNotification("No datapoints after removing NA frequencies or NA group values; nothing to plot.", type = "warning")
        return(invisible(NULL))
      }

      # Determine if we're doing pairwise comparisons
      is_pairwise <- !is.null(test_type) && test_type == "Pairwise Wilcoxon"

      # Prepare p-value annotation dataframe
      if (is_pairwise && "comparison" %in% names(res)) {
        # For pairwise comparisons, adjust p-values across all comparisons
        if (use_adj_p && nrow(res) > 0) {
          res$padj <- p.adjust(res$p, method = input$cat_p_adj_method)
        }
        # We'll use ggpubr for pairwise annotations
        p_df <- NULL
      } else {
        # Single p-value per entity
        p_df <- res %>%
          dplyr::mutate(
            p_to_show = if (isTRUE(use_adj_p) && "padj" %in% names(res)) padj else p,
            label = paste0("p = ", sapply(p_to_show, function(x) {
              if (is.na(x)) {
                NA_character_
              } else if (x < 0.001) {
                formatC(x, format = "e", digits = 2)
              } else {
                as.character(signif(x, 3))
              }
            })),
            x = length(unique(abund_long_plot[[group_var]])) / 2 + 0.5,
            y = tapply(abund_long_plot$freq, abund_long_plot$entity, max, na.rm = TRUE)[entity] * 1.15
          )
      }

      # Determine group levels present in plotting data
      grp_levels <- unique(as.character(abund_long_plot[[group_var]]))
      if (length(grp_levels) == 0) {
        showNotification("No valid group levels found for plotting.", type = "error")
        return(invisible(NULL))
      }

      # Use stored colors from cat_state if available, otherwise generate defaults
      colors_named <- if (!is.null(cp$colors) && length(cp$colors) > 0) {
        matched <- cp$colors[names(cp$colors) %in% grp_levels]
        missing_lvls <- setdiff(grp_levels, names(matched))
        if (length(missing_lvls) > 0) {
          filler <- viridis::viridis(length(missing_lvls))
          names(filler) <- missing_lvls
          matched <- c(matched, as.character(filler))
        }
        cols_ord <- as.character(matched[grp_levels])
        names(cols_ord) <- grp_levels
        cols_ord
      } else {
        pal <- viridis::viridis(length(grp_levels))
        setNames(as.character(pal), grp_levels)
      }

      # Build ggplot
      gg <- ggplot2::ggplot(abund_long_plot, ggplot2::aes(x = .data[[group_var]], y = freq))

      if (identical(plot_type, "violin")) {
        gg <- gg + ggplot2::geom_violin(ggplot2::aes(fill = .data[[group_var]]), alpha = 0.7, scale = "width", trim = TRUE)
      } else {
        gg <- gg + ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[group_var]]), alpha = 0.7, outlier.shape = NA)
      }

      # Add points per captured point_mode (points use pch=21 with black border and fill mapped to group)
      if (identical(point_mode, "draw")) {
        gg <- gg + ggplot2::geom_point(
          ggplot2::aes(fill = .data[[group_var]]),
          position = ggplot2::position_jitter(width = 0.04, height = 0),
          pch = 21, size = 2.6, stroke = 0.3, color = "black", alpha = 0.9, inherit.aes = TRUE
        )
      } else if (identical(point_mode, "jitter")) {
        gg <- gg + ggplot2::geom_jitter(
          ggplot2::aes(fill = .data[[group_var]]),
          width = 0.15, height = 0, pch = 21, size = 2.6, stroke = 0.3, color = "black", alpha = 0.9
        )
      } else if (identical(point_mode, "connected")) {
        # Draw points without jitter and add lines connecting paired samples
        # Only draw connections if pairing is actually feasible for this categorical variable
        use_pairing <- cp$use_pairing %||% FALSE
        pairing_var <- cp$pairing_var

        if (use_pairing && !is.null(pairing_var) && nzchar(pairing_var) && pairing_var %in% colnames(abund_long_plot)) {
          # Add connecting lines first (so they appear behind points)
          gg <- gg + ggplot2::geom_line(
            ggplot2::aes(group = .data[[pairing_var]]),
            color = "black",
            linewidth = 0.3,
            alpha = 0.6
          )
        }
        # Add points without jitter
        gg <- gg + ggplot2::geom_point(
          ggplot2::aes(fill = .data[[group_var]]),
          pch = 21, size = 2.6, stroke = 0.3, color = "black", alpha = 0.9
        )
      } else {
        # none -> no point layer
      }

      gg <- gg +
        ggplot2::facet_wrap(~entity, ncol = facet_cols, scales = "free_y") +
        ggplot2::scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::labs(x = group_var, y = "Frequency") +
        ggplot2::theme(
          strip.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 1.1, b = 1.1)),
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1) # rotate facet x-axis labels 45 degrees
        )

      # Apply manual fill scale; points use black border so no color scale required for border
      gg <- gg + ggplot2::scale_fill_manual(values = colors_named)

      # Add p-value annotations
      if (is_pairwise && "comparison" %in% names(res)) {
        # Build ALL bracket and label data at once (much faster than adding layers in loops)
        all_brackets <- list()
        all_labels <- list()

        for (ent in unique(abund_long_plot$entity)) {
          ent_res <- res %>% dplyr::filter(entity == ent)
          if (nrow(ent_res) == 0) next

          # Build comparison list for this entity
          comparisons_list <- lapply(1:nrow(ent_res), function(i) {
            comp_str <- ent_res$comparison[i]
            groups <- strsplit(comp_str, " vs ")[[1]]
            if (length(groups) == 2) groups else NULL
          })
          comparisons_list <- Filter(Negate(is.null), comparisons_list)

          if (length(comparisons_list) > 0) {
            # Get p-values to display
            p_vals <- if (use_adj_p && "padj" %in% names(ent_res)) {
              ent_res$padj
            } else {
              ent_res$p
            }

            # Get entity data
            ent_data <- abund_long_plot %>% dplyr::filter(entity == ent)
            y_max <- max(ent_data$freq, na.rm = TRUE)
            grp_levels_ordered <- sort(unique(as.character(ent_data[[group_var]])))

            # Determine number of unique categories for this entity
            n_categories <- length(unique(ent_data[[group_var]]))

            # Scale spacing and offset based on number of categories
            # More categories = need more vertical space between brackets
            # Use non-linear scaling: 2 and 4+ categories need more space
            if (n_categories == 2) {
              spacing_multiplier <- 0.12
            } else if (n_categories == 3) {
              spacing_multiplier <- 0.13
            } else {
              # 4+ categories: scale more aggressively
              spacing_multiplier <- 0.13 + 0.07 * (n_categories - 3)
            }

            # Text offset is consistently less than bracket spacing to avoid overlap
            text_offset_multiplier <- spacing_multiplier * 0.55

            # Format p-values: use actual p-values for <3 categories, stars for >=3
            if (n_categories < 3) {
              p_labels <- sapply(p_vals, function(x) {
                if (is.na(x)) {
                  "NA"
                } else if (x < 0.001) {
                  formatC(x, format = "e", digits = 2)
                } else {
                  as.character(signif(x, 3))
                }
              })
            } else {
              p_labels <- sapply(p_vals, function(x) {
                if (is.na(x)) {
                  "ns"
                } else if (x < 0.001) {
                  "***"
                } else if (x < 0.01) {
                  "**"
                } else if (x < 0.05) {
                  "*"
                } else {
                  "ns"
                }
              })
            }

            # Build brackets and labels for all comparisons in this entity
            for (j in seq_along(comparisons_list)) {
              grp1 <- comparisons_list[[j]][1]
              grp2 <- comparisons_list[[j]][2]
              y_pos <- y_max * (1.05 + spacing_multiplier * (j - 1))
              label <- p_labels[j]

              x1 <- which(grp_levels_ordered == grp1)
              x2 <- which(grp_levels_ordered == grp2)

              if (length(x1) > 0 && length(x2) > 0) {
                # Collect bracket segments
                all_brackets[[length(all_brackets) + 1]] <- data.frame(
                  entity = ent,
                  x = c(x1, x1, x2, x2),
                  xend = c(x1, x2, x2, x2),
                  y = c(y_pos - 0.02 * y_max, y_pos, y_pos, y_pos - 0.02 * y_max),
                  yend = c(y_pos, y_pos, y_pos, y_pos - 0.02 * y_max),
                  stringsAsFactors = FALSE
                )

                # Collect labels
                all_labels[[length(all_labels) + 1]] <- data.frame(
                  entity = ent,
                  x = (x1 + x2) / 2,
                  y = y_pos + text_offset_multiplier * y_max,
                  label = label,
                  stringsAsFactors = FALSE
                )
              }
            }
          }
        }

        # Add all brackets as a single layer
        if (length(all_brackets) > 0) {
          bracket_data_all <- do.call(rbind, all_brackets)
          gg <- gg + ggplot2::geom_segment(
            data = bracket_data_all,
            ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
            inherit.aes = FALSE,
            color = "black",
            linewidth = 0.5
          )
        }

        # Add all labels as a single layer
        if (length(all_labels) > 0) {
          label_data_all <- do.call(rbind, all_labels)
          gg <- gg + ggplot2::geom_text(
            data = label_data_all,
            ggplot2::aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 4.5,
            fontface = "bold"
          )
        }
      } else if (!is.null(p_df) && nrow(p_df) > 0) {
        # Single p-value text annotations (filter to plotted entities)
        p_df_plot <- p_df %>% dplyr::filter(entity %in% unique(abund_long_plot$entity))
        if (nrow(p_df_plot) > 0) {
          gg <- gg + ggplot2::geom_text(
            data = p_df_plot,
            ggplot2::aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 5
          )
        }
      }

      # Add caption explaining significance notation for pairwise tests (only if using stars)
      if (is_pairwise && "comparison" %in% names(res)) {
        # Check if any entity has >=3 categories (meaning stars are used)
        use_stars <- any(sapply(unique(abund_long_plot$entity), function(ent) {
          ent_data <- abund_long_plot %>% dplyr::filter(entity == ent)
          length(unique(ent_data[[group_var]])) >= 3
        }))

        if (use_stars) {
          gg <- gg + ggplot2::labs(
            caption = "Significance levels: *** p < 0.001, ** p < 0.01, * p < 0.05, ns = not significant"
          ) + ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0.5))
        }
      }

      # Cache for export
      cat_plot_cache(gg)

      gg
    },
    height = function() {
      gg <- cat_plot_cache()
      cp <- cat_state()
      if (is.null(gg) || is.null(cp)) {
        return(400)
      }
      n_facets <- length(unique(gg$data$entity))
      ncol_facets <- cp$facet_cols
      if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
      nrow_facets <- ceiling(n_facets / ncol_facets)

      # Increase height for pairwise comparisons to accommodate brackets
      is_pairwise <- !is.null(cp$test_type) && cp$test_type == "Pairwise Wilcoxon"
      base_height <- if (is_pairwise) 250 else 200
      base_height * nrow_facets
    },
    width = function() {
      gg <- cat_plot_cache()
      cp <- cat_state()
      if (is.null(gg) || is.null(cp)) {
        return(400)
      }
      ncol_facets <- cp$facet_cols
      if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1

      # Increase width for pairwise comparisons to accommodate annotations
      is_pairwise <- !is.null(cp$test_type) && cp$test_type == "Pairwise Wilcoxon"
      base_width <- if (is_pairwise) 280 else 225
      base_width * ncol_facets
    }
  )

  observeEvent(input$cat_group_var, {
    output$cat_color_pickers_ui <- renderUI(NULL)
    rv$cat_colors <- NULL
  })

  output$export_cat_pdf <- downloadHandler(
    filename = function() {
      subset_id <- rv$subset_id %||% "000000000"
      info <- rv$last_cat_info
      if (is.null(info)) {
        return(paste0("categorical_plots_", subset_id, ".pdf"))
      }
      test <- tolower(info$test_raw)
      test <- gsub("\\s+", "_", test)
      test <- gsub("_\\(.*\\)", "", test)
      test <- gsub("kruskal_wallis", "kruskal-wallis", test)
      paste0("categorical_", info$entity, "_", info$group, "_", test, "_", subset_id, ".pdf")
    },
    content = function(file) {
      gg <- cat_plot_cache()
      cp <- cat_plot_data()
      if (is.null(gg) || is.null(cp)) {
        showNotification("No plot available to export. Please generate the plots first.", type = "error")
        return()
      }
      n_facets <- length(unique(gg$data$entity))
      ncol_facets <- cp$facet_cols
      if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
      nrow_facets <- ceiling(n_facets / ncol_facets)
      if (!is.finite(nrow_facets) || nrow_facets < 1) nrow_facets <- 1
      pdf_width <- 4 * ncol_facets
      pdf_height <- 3 * nrow_facets
      ggsave(file,
        plot = gg, device = if (capabilities("cairo")) cairo_pdf else pdf,
        width = pdf_width, height = pdf_height, units = "in"
      )
    },
    contentType = "application/pdf"
  )

  cont_plot_data <- eventReactive(input$generate_cont_plots, {
    req(rv$meta_sample, rv$abundance_sample)

    abund0 <- rv$abundance_sample
    if (is.null(abund0)) {
      showNotification("No abundance matrix available.", type = "error")
      return(NULL)
    }

    # Aggregate to celltypes if requested
    abund <- abund0
    if (input$cont_entity == "Celltypes" && !is.null(rv$cluster_map)) {
      cm <- rv$cluster_map
      keep <- cm$cluster %in% colnames(abund)
      cm <- cm[keep, , drop = FALSE]
      if (!nrow(cm)) {
        return(NULL)
      }
      split_idx <- split(cm$cluster, cm$celltype)
      abund <- sapply(split_idx, function(cols) rowSums(abund0[, cols, drop = FALSE]))
      abund <- as.matrix(abund)
    }

    # Merge with per-sample metadata
    meta_sub <- rv$meta_sample %>% dplyr::select(patient_ID, dplyr::everything())
    abund_df <- as.data.frame(abund, check.names = FALSE, stringsAsFactors = FALSE)
    abund_df$patient_ID <- stringr::str_extract(string = rownames(abund_df), pattern = paste0("(", paste0(meta_sub$patient_ID, collapse = "|"), ")"))
    # merged <- dplyr::left_join(meta_sub, abund_df, by = "patient_ID")
    merged <- merge(x = meta_sub, y = abund_df, by = "patient_ID")

    # Long format; remove NA freqs immediately
    abund_long <- merged %>%
      tidyr::pivot_longer(cols = colnames(abund), names_to = "entity", values_to = "freq") %>%
      dplyr::mutate(entity = gsub("\\n", " ", entity)) %>%
      dplyr::filter(!is.na(freq))

    # Capture inputs at Generate time
    cont_var <- input$cont_group_var
    transpose_flag <- isTRUE(input$cont_transpose)
    test_padj_method <- input$cont_p_adj_method

    # Run Spearman per entity
    res <- abund_long %>%
      dplyr::group_by(entity) %>%
      dplyr::group_modify(~ {
        if (is.null(cont_var) || !nzchar(cont_var) || !(cont_var %in% colnames(.x))) {
          return(data.frame(rho = NA_real_, p = NA_real_, n = 0))
        }
        ok <- complete.cases(.x$freq, .x[[cont_var]])
        if (!any(ok)) {
          return(data.frame(rho = NA_real_, p = NA_real_, n = 0))
        }
        ct <- suppressWarnings(cor.test(.x$freq[ok], .x[[cont_var]][ok], method = "spearman"))
        data.frame(rho = unname(ct$estimate), p = ct$p.value, n = sum(ok))
      }) %>%
      dplyr::ungroup()

    # Adjust p-values if requested
    if (nrow(res) && nzchar(test_padj_method) && "p" %in% names(res)) {
      res$padj <- p.adjust(res$p, method = test_padj_method)
    }

    # Save info for export
    rv$last_cont_info <- list(
      entity = tolower(input$cont_entity %||% "clusters"),
      group = tolower(cont_var %||% "group"),
      test_raw = "spearman",
      transpose = transpose_flag
    )

    # Return everything needed for plotting
    list(
      data = abund_long,
      results = res,
      cont_var = cont_var,
      use_adj_p = input$cont_use_adj_p,
      facet_cols = as.numeric(input$cont_max_facets),
      transpose = transpose_flag
    )
  })

  observeEvent(input$generate_cont_plots, {
    cp <- cont_plot_data() # your existing eventReactive
    cont_state(cp)
    output$cont_cleared_msg <- renderText(NULL)
  })

  observeEvent(input$reset_cont, {
    cont_state(NULL) # clear the state
    cont_plot_cache(NULL) # clear cached ggplot
    showNotification("Continuous plots cleared.", type = "message", duration = 5)
    output$cont_cleared_msg <- renderText("Results cleared. Generate new plots to see them here.")
  })

  output$hasContResults <- reactive({
    cp <- cont_state()
    !is.null(cp) && !is.null(cp$data) && nrow(cp$data) > 0
  })
  outputOptions(output, "hasContResults", suspendWhenHidden = FALSE)

  output$continuous_plot <- renderPlot(
    {
      # cp <- cont_plot_data(); req(cp)
      cp <- cont_state()
      req(cp)
      abund_long <- cp$data
      res <- cp$results
      cont_var <- cp$cont_var
      use_adj_p <- cp$use_adj_p
      facet_cols <- cp$facet_cols
      transpose_flag <- cp$transpose %||% FALSE

      # Validate continuous variable
      if (is.null(cont_var) || !nzchar(cont_var) || !(cont_var %in% colnames(abund_long))) {
        showNotification("No valid continuous metadata variable selected for continuous plotting.", type = "error")
        return(invisible(NULL))
      }

      # Remove rows with NA in either axis
      plot_df <- abund_long %>% dplyr::filter(!is.na(freq) & !is.na(.data[[cont_var]]))
      if (nrow(plot_df) == 0) {
        showNotification("No datapoints available for plotting after removing missing values.", type = "warning")
        return(invisible(NULL))
      }

      # Prepare p-value annotation dataframe
      p_df <- res %>%
        dplyr::mutate(
          p_to_show = if (isTRUE(use_adj_p) && "padj" %in% names(res)) padj else p,
          label = paste0("p = ", signif(p_to_show, 3), "\n", "rho = ", signif(rho, 3))
        )

      # Choose aesthetics and label positions
      if (!transpose_flag) {
        aes_pt <- ggplot2::aes(x = freq, y = .data[[cont_var]])
        smooth_aes <- ggplot2::aes(x = freq, y = .data[[cont_var]])
        x_lab <- "Abundance"
        y_lab <- cont_var
        p_df <- p_df %>%
          dplyr::mutate(
            x = tapply(plot_df$freq, plot_df$entity, function(v) mean(range(v, na.rm = TRUE)))[entity],
            y = tapply(plot_df[[cont_var]], plot_df$entity, max, na.rm = TRUE)[entity] * 1.05
          )
      } else {
        aes_pt <- ggplot2::aes(x = .data[[cont_var]], y = freq)
        smooth_aes <- ggplot2::aes(x = .data[[cont_var]], y = freq)
        x_lab <- cont_var
        y_lab <- "Abundance"
        p_df <- p_df %>%
          dplyr::mutate(
            x = tapply(plot_df[[cont_var]], plot_df$entity, function(v) mean(range(v, na.rm = TRUE)))[entity],
            y = tapply(plot_df$freq, plot_df$entity, max, na.rm = TRUE)[entity] * 1.05
          )
      }

      # Build ggplot with your preferred aesthetics
      gg <- ggplot2::ggplot(plot_df, mapping = aes_pt) +
        ggplot2::geom_point(
          alpha = 0.75, pch = 21, color = "black", fill = "grey40",
          stroke = 0.1, size = 3
        ) +
        ggplot2::geom_smooth(mapping = smooth_aes, method = "lm", se = FALSE, color = "red2") +
        ggplot2::facet_wrap(~entity, ncol = facet_cols, scales = "free") +
        ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
        ggplot2::scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
        ggplot2::theme_bw(base_size = 18) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(x = x_lab, y = y_lab) +
        ggplot2::theme(strip.text.x = ggplot2::element_text(
          margin = ggplot2::margin(t = 1.1, b = 1.1)
        ))

      # Add p/rho annotations
      if (!is.null(p_df) && nrow(p_df) > 0) {
        p_df_plot <- p_df %>% dplyr::filter(entity %in% unique(plot_df$entity))
        if (nrow(p_df_plot) > 0) {
          gg <- gg + ggplot2::geom_text(
            data = p_df_plot,
            ggplot2::aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            size = 5,
            lineheight = 0.85
          )
        }
      }

      cont_plot_cache(gg)
      gg
    },
    height = function() {
      gg <- cont_plot_cache()
      cp <- cont_plot_data()
      if (is.null(gg) || is.null(cp)) {
        return(400)
      }
      n_facets <- length(unique(gg$data$entity))
      ncol_facets <- cp$facet_cols
      if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
      nrow_facets <- ceiling(n_facets / ncol_facets)
      200 * nrow_facets
    },
    width = function() {
      gg <- cont_plot_cache()
      cp <- cont_plot_data()
      if (is.null(gg) || is.null(cp)) {
        return(400)
      }
      ncol_facets <- cp$facet_cols
      if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
      225 * ncol_facets
    }
  )

  output$export_cont_pdf <- downloadHandler(
    filename = function() {
      subset_id <- rv$subset_id %||% "000000000"
      info <- rv$last_cont_info
      if (is.null(info)) {
        return(paste0("continuous_plots_", subset_id, ".pdf"))
      }
      paste0("continuous_", info$entity, "_", info$group, "_spearman_", subset_id, ".pdf")
    },
    content = function(file) {
      gg <- cont_plot_cache()
      cp <- cont_plot_data()
      if (is.null(gg) || is.null(cp)) {
        showNotification("No plot available to export. Please generate the plots first.", type = "error")
        return()
      }
      n_facets <- length(unique(gg$data$entity))
      ncol_facets <- cp$facet_cols
      if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
      nrow_facets <- ceiling(n_facets / ncol_facets)
      if (!is.finite(nrow_facets) || nrow_facets < 1) nrow_facets <- 1
      pdf_width <- 4 * ncol_facets
      pdf_height <- 3 * nrow_facets
      ggsave(file, plot = gg, device = if (capabilities("cairo")) cairo_pdf else pdf, width = pdf_width, height = pdf_height, units = "in")
    },
    contentType = "application/pdf"
  )

  # Populate Feature Selection dropdowns from same metadata source as other tabs
  observeEvent(rv$meta_cell,
    {
      meta_cols <- sort(colnames(rv$meta_cell))
      # Exclude cluster, patient_ID, run_date, and source from outcome choices
      outcome_choices <- setdiff(meta_cols, c("cluster", "patient_ID", "run_date", "source"))
      # Exclude patient_ID, source, and run_date from predictor choices
      predictor_choices <- setdiff(meta_cols, c("patient_ID", "source", "run_date"))
      updatePickerInput(session, "fs_outcome", choices = outcome_choices)
      updatePickerInput(session, "fs_predictors", choices = predictor_choices)
    },
    ignoreInit = TRUE
  )

  observeEvent(rv$clusters$abundance, {
    cluster_names <- colnames(rv$clusters$abundance)
    updatePickerInput(session, "fs_cluster_subset", choices = cluster_names)
  })

  run_fs <- function() {
    req(rv$meta_sample, rv$abundance_sample, input$fs_outcome, input$fs_predictors)

    # Metadata predictors (exclude placeholder)
    pred_meta <- setdiff(intersect(input$fs_predictors, colnames(rv$meta_sample)), "cluster")

    # Cluster predictors
    all_clusters <- colnames(rv$abundance_sample)
    if ("cluster" %in% input$fs_predictors) {
      if (!is.null(input$fs_cluster_subset) && length(input$fs_cluster_subset) > 0) {
        cluster_predictors <- intersect(input$fs_cluster_subset, all_clusters)
      } else {
        cluster_predictors <- all_clusters
      }
    } else {
      cluster_predictors <- character(0)
    }

    predictors_final <- c(pred_meta, cluster_predictors)
    if (length(predictors_final) == 0) {
      showNotification("Select at least one predictor.", type = "error")
      return(NULL)
    }

    # Subset before merge
    meta_sub <- rv$meta_sample %>%
      dplyr::select(patient_ID, !!input$fs_outcome, dplyr::all_of(pred_meta))
    abund_sub <- rv$abundance_sample[, cluster_predictors, drop = FALSE]

    merged <- align_metadata_abundance(meta_sub, abund_sub, notify = showNotification)

    # Filter missingness
    merged <- merged[!is.na(merged[[input$fs_outcome]]), ]
    y_raw <- merged[[input$fs_outcome]]
    X_raw <- merged[, predictors_final, drop = FALSE]

    n_before <- nrow(merged)
    complete_rows <- complete.cases(data.frame(X_raw, .y = y_raw))
    dropped_ids <- merged$patient_ID[!complete_rows]

    X_raw <- X_raw[complete_rows, , drop = FALSE]
    y_raw <- y_raw[complete_rows]
    merged <- merged[complete_rows, , drop = FALSE]

    n_after <- sum(complete_rows)
    n_dropped <- n_before - n_after
    if (n_after < 3) {
      showNotification("Too few samples after filtering for feature selection.", type = "error")
      return(NULL)
    }

    y <- if (is.character(y_raw)) factor(y_raw) else y_raw

    # User controls
    method <- input$fs_method %||% "Elastic Net"
    tolerance <- 1e-4
    seed_val <- input$fs_seed %||% 123
    reps <- input$fs_reps %||% 1
    boruta_maxruns <- input$fs_maxruns %||% 500
    set.seed(seed_val)

    clean_dummy_names <- function(nms) {
      gsub(pattern = "cluster", replacement = "cluster:", x = nms)
    }

    # --- Boruta branch ---
    if (method == "Random Forest (Boruta)") {
      all_imp <- list()
      for (r in seq_len(reps)) {
        set.seed(seed_val + r - 1)
        X_boruta <- model.matrix(~ . - 1, data = X_raw)
        colnames(X_boruta) <- clean_dummy_names(colnames(X_boruta))
        X_boruta <- as.data.frame(X_boruta)

        if (is.factor(y)) {
          bor <- Boruta::Boruta(x = X_boruta, y = y, doTrace = 2, maxRuns = boruta_maxruns)
        } else {
          bor <- Boruta::Boruta(x = X_boruta, y = as.numeric(y), doTrace = 2, maxRuns = boruta_maxruns)
        }
        imp <- Boruta::attStats(bor)
        imp$meanImp[!is.finite(imp$meanImp)] <- 0
        all_imp[[r]] <- imp
      }
      merged_imp <- Reduce(function(a, b) {
        common <- intersect(rownames(a), rownames(b))
        a[common, "meanImp"] <- (a[common, "meanImp"] + b[common, "meanImp"]) / 2
        a
      }, all_imp)
      merged_imp <- merged_imp[abs(merged_imp$meanImp) >= tolerance, , drop = FALSE]

      sel <- rownames(merged_imp)[merged_imp$decision %in% c("Confirmed", "Tentative")]
      res_df <- data.frame(
        Feature = rownames(merged_imp),
        ImportanceMean = merged_imp$meanImp,
        Decision = merged_imp$decision,
        stringsAsFactors = FALSE
      )

      return(list(
        method = "Boruta",
        results = res_df,
        selected = sel,
        outcome = y,
        predictors = colnames(X_raw),
        tolerance = tolerance,
        details = list(
          samples_before = n_before,
          samples_after = n_after,
          samples_dropped = n_dropped,
          dropped_ids = dropped_ids
        ),
        merged = merged
      ))
    }

    # --- glmnet branch ---
    family <- if (is.factor(y)) {
      if (nlevels(y) == 2) "binomial" else "multinomial"
    } else {
      "gaussian"
    }

    # Create model matrix first
    Xmat <- model.matrix(~ . - 1, data = X_raw)
    colnames(Xmat) <- clean_dummy_names(colnames(Xmat))
    storage.mode(Xmat) <- "double"

    added_dummy <- FALSE
    if (ncol(Xmat) == 1) {
      Xmat <- cbind(Xmat, `__DUMMY__` = 0)
      added_dummy <- TRUE
    }

    # Explicitly standardize all continuous columns in Xmat (mean=0, sd=1)
    # Identify continuous columns (non-dummy variables with range > 0)
    # Dummy variables from factors typically have only 0/1 values
    X_scale_params <- list()
    for (j in seq_len(ncol(Xmat))) {
      col <- Xmat[, j]
      unique_vals <- unique(col)
      # If column has more than 2 unique values, treat as continuous and standardize
      # This preserves dummy variables (0/1) while scaling continuous predictors
      if (length(unique_vals) > 2) {
        col_mean <- mean(col, na.rm = TRUE)
        col_sd <- sd(col, na.rm = TRUE)
        if (col_sd > 0) {
          Xmat[, j] <- (col - col_mean) / col_sd
          X_scale_params[[colnames(Xmat)[j]]] <- list(mean = col_mean, sd = col_sd)
        }
      }
    }

    # Standardize continuous outcomes for better penalized regression performance
    y_scaled <- y
    y_scale_params <- NULL
    if (family == "gaussian" && is.numeric(y)) {
      y_mean <- mean(y, na.rm = TRUE)
      y_sd <- sd(y, na.rm = TRUE)
      if (y_sd > 0) {
        y_scaled <- (y - y_mean) / y_sd
        y_scale_params <- list(mean = y_mean, sd = y_sd)
      }
    }

    alpha_val <- if (method == "Ridge Regression") 0 else (input$fs_alpha %||% 0.5)
    nfolds_val <- input$fs_nfolds %||% 5
    if (nfolds_val > n_after) nfolds_val <- max(3, floor(n_after / 2))

    coef_list <- list()
    for (r in seq_len(reps)) {
      set.seed(seed_val + r - 1)

      # Wrap model fitting in tryCatch to handle errors gracefully
      cvfit <- tryCatch(
        {
          glmnet::cv.glmnet(
            x = Xmat, y = y_scaled, family = family,
            alpha = alpha_val, nfolds = nfolds_val,
            standardize = FALSE # We already standardized manually
          )
        },
        error = function(e) {
          err_msg <- conditionMessage(e)
          if (grepl("one multinomial or binomial class has 1 or 0 observations", err_msg, ignore.case = TRUE)) {
            stop("Insufficient observations in one or more outcome classes. Ensure each class has at least 8 samples, or consider using a different outcome variable or predictor set.", call. = FALSE)
          } else if (grepl("fewer than 8", err_msg, ignore.case = TRUE)) {
            stop("One or more outcome classes have fewer than 8 observations, which is insufficient for reliable model fitting. Consider grouping rare classes or using a different outcome variable.", call. = FALSE)
          } else {
            stop(paste0("Model fitting failed: ", err_msg), call. = FALSE)
          }
        }
      )

      coef_obj <- coef(cvfit, s = "lambda.min")
      if (family == "multinomial") {
        coef_df <- do.call(rbind, lapply(names(coef_obj), function(cls) {
          cm <- as.matrix(coef_obj[[cls]])
          data.frame(Feature = rownames(cm), Class = cls, Coef = as.numeric(cm[, 1]), stringsAsFactors = FALSE)
        }))
        coef_df <- coef_df[coef_df$Feature != "(Intercept)", , drop = FALSE]
        if (added_dummy) coef_df <- coef_df[coef_df$Feature != "__DUMMY__", , drop = FALSE]
        coef_list[[r]] <- coef_df
      } else {
        cm <- as.matrix(coef_obj)
        coef_raw <- as.numeric(cm[, 1])
        feat_names <- rownames(cm)

        # For Feature Selection, keep coefficients in STANDARDIZED space
        # This makes them comparable across features regardless of original units
        # A coefficient of -1.0 means "1 SD increase in predictor → 1 SD decrease in outcome"
        coef_df <- data.frame(Feature = feat_names, Coef = coef_raw, stringsAsFactors = FALSE)
        coef_df <- coef_df[coef_df$Feature != "(Intercept)", , drop = FALSE]
        if (added_dummy) coef_df <- coef_df[coef_df$Feature != "__DUMMY__", , drop = FALSE]
        coef_list[[r]] <- coef_df
      }
    }

    if (family == "multinomial") {
      all_df <- do.call(rbind, coef_list)
      agg <- all_df %>%
        dplyr::group_by(Feature, Class) %>%
        dplyr::summarise(Coef = median(Coef, na.rm = TRUE), .groups = "drop")
      agg <- agg[abs(agg$Coef) >= tolerance, , drop = FALSE]
      res_df <- agg[order(agg$Coef, decreasing = TRUE), ]
      selected <- unique(head(res_df$Feature, 50))
    } else {
      all_df <- do.call(rbind, coef_list)
      agg <- all_df %>%
        dplyr::group_by(Feature) %>%
        dplyr::summarise(Coef = median(Coef, na.rm = TRUE), .groups = "drop")
      agg <- agg[abs(agg$Coef) >= tolerance, , drop = FALSE]
      res_df <- agg[order(agg$Coef, decreasing = TRUE), ]
      selected <- head(res_df$Feature, 50)
    }

    return(list(
      method = if (alpha_val == 0) "Ridge" else "Elastic Net",
      family = family,
      results = res_df,
      selected = selected,
      outcome = y,
      predictors = colnames(Xmat),
      tolerance = tolerance,
      details = list(
        samples_before = n_before,
        samples_after = n_after,
        samples_dropped = n_dropped,
        dropped_ids = dropped_ids
      ),
      merged = merged
    ))
  }

  output$fs_data_head <- renderTable(
    {
      # res <- run_fs(); req(res)
      res <- fs_state()
      req(res)
      head(res$merged, 5)
    },
    sanitize.text.function = function(x) x
  )

  # Dynamic height for Feature Selection plot based on facet rows and features
  fs_plot_height <- reactive({
    res <- fs_state()
    if (is.null(res)) {
      return(550)
    }

    tol <- res$tolerance
    df <- res$results

    # Filter by tolerance to count features that will actually be plotted
    if (identical(res$method, "Boruta")) {
      df <- df[abs(df$ImportanceMean) >= tol, , drop = FALSE]
    } else {
      df <- df[abs(df$Coef) >= tol, , drop = FALSE]
    }

    if (nrow(df) == 0) {
      return(550)
    }

    # For multinomial outcomes, calculate height based on facet rows and features
    if ("Class" %in% names(df)) {
      n_classes <- length(unique(df$Class))
      # Count max features in any class (this is the height per facet)
      max_features_per_class <- max(table(df$Class))

      # ggplot2 facet_wrap default behavior: uses roughly sqrt(n) columns
      # We calculate number of rows based on default wrapping
      ncol_facets <- ceiling(sqrt(n_classes))
      n_rows <- ceiling(n_classes / ncol_facets)

      # Height calculation: base + (features per facet * pixels per feature * number of rows)
      # Reduced values for more compact display
      base_height <- 150
      pixels_per_feature <- 20
      height <- base_height + (max_features_per_class * pixels_per_feature * n_rows)

      return(max(400, min(height, 1500))) # Lower minimum and cap
    } else {
      # For binary/continuous, single plot - base on total features
      n_features <- nrow(df)
      height <- 200 + (n_features * 20)
      return(max(400, min(height, 1200))) # Cap at 1200px
    }
  })

  output$fs_plot_ui <- renderUI({
    height_val <- fs_plot_height()
    plotOutput("fs_plot", height = paste0(height_val, "px"))
  })

  output$fs_plot <- renderPlot({
    # res <- run_fs(); req(res)
    res <- fs_state()
    req(res)

    # Check for error
    if (!is.null(res$error) && isTRUE(res$error)) {
      plot.new()
      text(0.5, 0.5, paste("Error:", res$message), col = "red", cex = 1.2)
      return()
    }

    tol <- res$tolerance

    if (identical(res$method, "Boruta")) {
      df <- res$results
      df <- df[abs(df$ImportanceMean) >= tol, , drop = FALSE]
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No Boruta features pass the tolerance threshold.\nNothing to plot.", cex = 1.1)
        return()
      }
      df <- df[order(df$ImportanceMean, decreasing = TRUE), ]
      # Remove backticks from feature names for cleaner plot labels
      df$Feature <- gsub("`", "", df$Feature, fixed = TRUE)
      # Plot all features that pass tolerance, not just top 30
      ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = reorder(Feature, ImportanceMean),
          y = ImportanceMean,
          fill = Decision
        )
      ) +
        ggplot2::geom_col(color = "black", lwd = 0.4) +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Boruta feature importance",
          x = "Feature", y = "Mean importance"
        ) +
        ggplot2::scale_fill_manual(values = c(
          "Confirmed" = "green4",
          "Rejected" = "red4",
          "Tentative" = "grey40"
        )) +
        ggplot2::theme_bw(base_size = 14)
    } else {
      df <- res$results
      df <- df[abs(df$Coef) >= tol, , drop = FALSE]
      if (nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "All coefficients are zero after regularization.\nTry lowering alpha or adding predictors.", cex = 1.1)
        return()
      }
      # Remove backticks from feature names for cleaner plot labels
      df$Feature <- gsub("`", "", df$Feature, fixed = TRUE)
      if ("Class" %in% names(df)) {
        # Plot all features that pass tolerance per class, not just top 20
        top_n <- df %>%
          dplyr::group_by(Class) %>%
          dplyr::arrange(dplyr::desc(Coef), .by_group = TRUE)
        ggplot2::ggplot(
          top_n,
          ggplot2::aes(
            x = reorder(Feature, Coef),
            y = Coef,
            fill = Coef
          )
        ) +
          ggplot2::geom_col(color = "black", lwd = 0.4) +
          ggplot2::facet_wrap(~Class, scales = "free_y") +
          ggplot2::coord_flip() +
          ggplot2::labs(
            title = paste(res$method, "coefficients (lambda.min)"),
            x = "Feature", y = "Coefficient"
          ) +
          ggplot2::scale_fill_gradient(low = "blue3", high = "red3") +
          ggplot2::theme_bw(base_size = 14)
      } else {
        df <- df[order(df$Coef, decreasing = TRUE), ]
        # Plot all features that pass tolerance, not just top 30
        ggplot2::ggplot(
          df,
          ggplot2::aes(
            x = reorder(Feature, Coef),
            y = Coef,
            fill = Coef
          )
        ) +
          ggplot2::geom_col(color = "black", lwd = 0.4) +
          ggplot2::coord_flip() +
          ggplot2::labs(
            title = paste(res$method, "coefficients (lambda.min)"),
            x = "Feature", y = "Coefficient"
          ) +
          ggplot2::scale_fill_gradient(low = "blue3", high = "red3") +
          ggplot2::theme_bw(base_size = 14)
      }
    }
  })

  output$fs_results <- renderTable(
    {
      # res <- run_fs(); req(res)
      res <- fs_state()
      req(res)

      # Check for error
      if (!is.null(res$error) && isTRUE(res$error)) {
        return(data.frame(Error = res$message))
      }

      df <- res$results
      tol <- res$tolerance

      if (identical(res$method, "Boruta")) {
        df <- df[abs(df$ImportanceMean) >= tol, , drop = FALSE]
        df <- df[order(df$ImportanceMean, decreasing = TRUE), , drop = FALSE]

        if (nrow(df) == 0) {
          return(data.frame(Message = "No Boruta features passed the tolerance threshold."))
        }
      } else {
        if ("Coef" %in% names(df)) {
          df <- df[abs(df$Coef) >= tol, , drop = FALSE]
          df <- df[order(df$Coef, decreasing = TRUE), , drop = FALSE]

          if (nrow(df) == 0) {
            return(data.frame(Message = "All coefficients shrank to zero after regularization."))
          }
        }
      }

      # Remove backticks and backslashes from feature names (added by model.matrix for numeric column names)
      if ("Feature" %in% names(df)) {
        df$Feature <- gsub("`", "", df$Feature, fixed = TRUE)
        df$Feature <- gsub("\\\\", "", df$Feature)
      }

      df
    },
    sanitize.text.function = function(x) x
  )

  output$fs_summary <- renderPrint({
    # res <- run_fs(); req(res)
    res <- fs_state()
    req(res)

    # Check for error
    if (!is.null(res$error) && isTRUE(res$error)) {
      cat("Feature Selection Error:\n\n")
      cat(res$message, "\n\n")
      cat("Suggestions:\n")
      cat("- Ensure each outcome class has at least 8 samples\n")
      cat("- Try a different outcome variable\n")
      cat("- Consider combining rare outcome classes\n")
      cat("- Use fewer predictors or select different predictors\n")
      return()
    }

    det <- res$details %||% list(
      samples_before = NA,
      samples_after = NA,
      samples_dropped = NA,
      dropped_ids = character(0)
    )

    cat("Samples before filtering:", det$samples_before, "\n")
    cat("Samples after filtering:", det$samples_after, "\n")
    cat("Samples dropped:", det$samples_dropped, "\n\n")

    if (length(det$dropped_ids) > 0) {
      cat("Dropped patient_IDs:\n")
      print(det$dropped_ids)
      cat("\n")
    }

    cat("Selected features (top):\n")
    if (identical(res$method, "Boruta")) {
      df <- res$results
      df <- df[order(df$ImportanceMean, decreasing = TRUE), , drop = FALSE]
      if (nrow(df) == 0) {
        cat("No Boruta features passed the tolerance threshold.\n")
      } else {
        print(utils::head(df$Feature, 20))
      }
    } else {
      df <- res$results
      if ("Coef" %in% names(df)) {
        df <- df[order(df$Coef, decreasing = TRUE), , drop = FALSE]
        if (nrow(df) == 0) {
          cat("All coefficients shrank to zero after regularization.\n")
        } else {
          print(utils::head(df$Feature, 20))
        }
      }
    }
  })

  output$hasFSResults <- reactive({
    res <- fs_state()
    !is.null(res) && !is.null(res$results) && nrow(res$results) > 0
  })
  outputOptions(output, "hasFSResults", suspendWhenHidden = FALSE)

  output$export_fs_results <- downloadHandler(
    filename = function() {
      subset_id <- rv$subset_id %||% "000000000"
      res <- fs_state()
      req(res)
      method <- gsub("\\s+", "_", tolower(res$method %||% "unknown"))
      outcome <- gsub("\\s+", "_", tolower(input$fs_outcome %||% "outcome"))

      # Add alpha to filename if Elastic Net
      if (identical(res$method, "Elastic Net")) {
        alpha_val <- input$fs_alpha %||% 0.5
        paste0("feature_selection_", method, "_", outcome, "_alpha", alpha_val, "_", subset_id, ".zip")
      } else {
        paste0("feature_selection_", method, "_", outcome, "_", subset_id, ".zip")
      }
    },
    content = function(file) {
      res <- fs_state()
      req(res)
      subset_id <- rv$subset_id %||% "000000000"

      tmpdir <- tempdir()
      files <- c()

      # Build filename prefix
      alpha_suffix <- if (identical(res$method, "Elastic Net")) {
        alpha_val <- input$fs_alpha %||% 0.5
        paste0("_alpha", alpha_val)
      } else {
        ""
      }

      # 1. Selected Features CSV
      feat_file <- file.path(tmpdir, paste0("feature_selection_selected_features", alpha_suffix, "_", subset_id, ".csv"))
      df <- res$results
      tol <- res$tolerance

      if (identical(res$method, "Boruta")) {
        df <- df[abs(df$ImportanceMean) >= tol, , drop = FALSE]
        df <- df[order(df$ImportanceMean, decreasing = TRUE), , drop = FALSE]
        if (nrow(df) == 0) {
          write.csv(data.frame(Message = "No Boruta features passed the tolerance threshold."),
            feat_file,
            row.names = FALSE
          )
        } else {
          write.csv(df, feat_file, row.names = FALSE)
        }
      } else {
        if ("Coef" %in% names(df)) {
          df <- df[abs(df$Coef) >= tol, , drop = FALSE]
          df <- df[order(df$Coef, decreasing = TRUE), , drop = FALSE]
          if (nrow(df) == 0) {
            write.csv(data.frame(Message = "All coefficients shrank to zero after regularization."),
              feat_file,
              row.names = FALSE
            )
          } else {
            write.csv(df, feat_file, row.names = FALSE)
          }
        }
      }
      files <- c(files, feat_file)

      # 2. Details CSV (key-value format)
      details_file <- file.path(tmpdir, paste0("feature_selection_details", alpha_suffix, "_", subset_id, ".csv"))
      det <- res$details %||% list(samples_before = NA, samples_after = NA, samples_dropped = NA)

      details_kv <- data.frame(
        Key = c(
          "Method", "Outcome variable", "Samples before filtering",
          "Samples after filtering", "Samples dropped"
        ),
        Value = c(
          res$method %||% NA,
          input$fs_outcome %||% NA,
          det$samples_before %||% NA,
          det$samples_after %||% NA,
          det$samples_dropped %||% NA
        ),
        stringsAsFactors = FALSE
      )

      # Add method-specific details
      if (identical(res$method, "Elastic Net") || identical(res$method, "Ridge Regression")) {
        alpha_val <- if (identical(res$method, "Ridge Regression")) 0 else (input$fs_alpha %||% 0.5)
        nfolds_val <- input$fs_nfolds %||% 5
        details_kv <- rbind(
          details_kv,
          data.frame(Key = "Alpha", Value = alpha_val, stringsAsFactors = FALSE),
          data.frame(Key = "CV folds", Value = nfolds_val, stringsAsFactors = FALSE)
        )
      } else if (identical(res$method, "Boruta")) {
        maxruns <- input$fs_maxruns %||% 500
        details_kv <- rbind(
          details_kv,
          data.frame(Key = "Max runs", Value = maxruns, stringsAsFactors = FALSE)
        )
      }

      write.csv(details_kv, details_file, row.names = FALSE)
      files <- c(files, details_file)

      # 3. Summary Plot PDF
      plot_file <- file.path(tmpdir, paste0("feature_selection_summary_plot", alpha_suffix, "_", subset_id, ".pdf"))

      # Regenerate the plot
      tol <- res$tolerance
      if (identical(res$method, "Boruta")) {
        df <- res$results
        df <- df[abs(df$ImportanceMean) >= tol, , drop = FALSE]
        if (nrow(df) > 0) {
          df <- df[order(df$ImportanceMean, decreasing = TRUE), ]
          df$Feature <- gsub("`", "", df$Feature, fixed = TRUE)

          gg <- ggplot2::ggplot(
            df,
            ggplot2::aes(
              x = reorder(Feature, ImportanceMean),
              y = ImportanceMean,
              fill = Decision
            )
          ) +
            ggplot2::geom_col(color = "black", lwd = 0.4) +
            ggplot2::coord_flip() +
            ggplot2::labs(
              title = "Boruta feature importance",
              x = "Feature", y = "Mean importance"
            ) +
            ggplot2::scale_fill_manual(values = c(
              "Confirmed" = "green4",
              "Rejected" = "red4",
              "Tentative" = "grey40"
            )) +
            ggplot2::theme_bw(base_size = 14)

          n_features <- nrow(df)
          plot_height <- max(6, min(n_features * 0.3, 20))
          ggplot2::ggsave(plot_file, gg, width = 10, height = plot_height)
          files <- c(files, plot_file)
        }
      } else {
        df <- res$results
        df <- df[abs(df$Coef) >= tol, , drop = FALSE]
        if (nrow(df) > 0) {
          df$Feature <- gsub("`", "", df$Feature, fixed = TRUE)

          if ("Class" %in% names(df)) {
            # Multinomial
            top_n <- df %>%
              dplyr::group_by(Class) %>%
              dplyr::arrange(dplyr::desc(Coef), .by_group = TRUE)

            gg <- ggplot2::ggplot(
              top_n,
              ggplot2::aes(
                x = reorder(Feature, Coef),
                y = Coef,
                fill = Coef
              )
            ) +
              ggplot2::geom_col(color = "black", lwd = 0.4) +
              ggplot2::facet_wrap(~Class, scales = "free_y") +
              ggplot2::coord_flip() +
              ggplot2::labs(
                title = paste(res$method, "coefficients (lambda.min)"),
                x = "Feature", y = "Coefficient"
              ) +
              ggplot2::scale_fill_gradient(low = "blue3", high = "red3") +
              ggplot2::theme_bw(base_size = 14)

            n_classes <- length(unique(df$Class))
            max_features_per_class <- max(table(df$Class))
            plot_height <- max(8, min((max_features_per_class * 0.3 * n_classes), 20))
            ggplot2::ggsave(plot_file, gg, width = 12, height = plot_height)
            files <- c(files, plot_file)
          } else {
            # Binary/continuous
            df <- df[order(df$Coef, decreasing = TRUE), ]

            gg <- ggplot2::ggplot(
              df,
              ggplot2::aes(
                x = reorder(Feature, Coef),
                y = Coef,
                fill = Coef
              )
            ) +
              ggplot2::geom_col(color = "black", lwd = 0.4) +
              ggplot2::coord_flip() +
              ggplot2::labs(
                title = paste(res$method, "coefficients (lambda.min)"),
                x = "Feature", y = "Coefficient"
              ) +
              ggplot2::scale_fill_gradient(low = "blue3", high = "red3") +
              ggplot2::theme_bw(base_size = 14)

            n_features <- nrow(df)
            plot_height <- max(6, min(n_features * 0.3, 20))
            ggplot2::ggsave(plot_file, gg, width = 10, height = plot_height)
            files <- c(files, plot_file)
          }
        }
      }

      # Create zip file
      zip::zip(file, files = basename(files), root = tmpdir, mode = "cherry-pick")
    },
    contentType = "application/zip"
  )

  observeEvent(list(rv$meta_sample, rv$abundance_sample),
    {
      req(rv$meta_sample, rv$abundance_sample)

      meta_cols <- colnames(rv$meta_sample)
      categorical_choices <- sort(meta_cols[sapply(rv$meta_sample, function(x) is.factor(x) || is.character(x))])
      # Exclude cluster, patient_ID, run_date, and source from outcome choices
      categorical_choices <- setdiff(categorical_choices, c("cluster", "patient_ID", "run_date", "source"))
      predictor_choices <- sort(meta_cols)
      # Exclude patient_ID, source, and run_date from predictor choices
      predictor_choices <- setdiff(predictor_choices, c("patient_ID", "source", "run_date"))
      if (!("cluster" %in% predictor_choices)) {
        predictor_choices <- c(predictor_choices, "cluster")
      }

      updatePickerInput(session, "lm_outcome", choices = categorical_choices, selected = NULL)
      updatePickerInput(session, "lm_predictors", choices = predictor_choices, selected = NULL)
      updatePickerInput(session, "lm_cluster_subset",
        choices = colnames(rv$abundance_sample),
        selected = character(0)
      )
    },
    ignoreInit = TRUE
  )

  # Initialize Regression dropdowns
  observeEvent(list(rv$meta_sample, rv$abundance_sample),
    {
      req(rv$meta_sample, rv$abundance_sample)

      meta_cols <- colnames(rv$meta_sample)
      continuous_choices <- sort(meta_cols[sapply(rv$meta_sample, is.numeric)])
      # Exclude patient_ID, source, run_date, and cluster from continuous outcome choices
      continuous_choices <- setdiff(continuous_choices, c("patient_ID", "source", "run_date", "cluster"))
      predictor_choices <- sort(meta_cols)
      # Exclude patient_ID, source, and run_date from predictor choices
      predictor_choices <- setdiff(predictor_choices, c("patient_ID", "source", "run_date"))
      if (!("cluster" %in% predictor_choices)) {
        predictor_choices <- c(predictor_choices, "cluster")
      }

      updatePickerInput(session, "reg_outcome", choices = continuous_choices, selected = NULL)
      updatePickerInput(session, "reg_predictors", choices = predictor_choices, selected = NULL)
      updatePickerInput(session, "reg_cluster_subset",
        choices = colnames(rv$abundance_sample),
        selected = character(0)
      )
    },
    ignoreInit = TRUE
  )

  observeEvent(input$reset_fs, {
    fs_state(NULL)
    showNotification("Feature Selection results cleared.", type = "message", duration = 5)
    output$fs_cleared_msg <- renderText("Results cleared. Run Feature Selection again to see results here.")
  })

  observeEvent(input$run_fs, {
    result <- tryCatch(
      {
        run_fs()
      },
      error = function(e) {
        # Return error object instead of crashing
        list(error = TRUE, message = conditionMessage(e))
      }
    )

    if (!is.null(result) && isTRUE(result$error)) {
      # Display error in the output area
      fs_state(result)
      showNotification(paste("Feature Selection failed:", result$message), type = "error", duration = 10)
    } else {
      fs_state(result)
    }
    output$fs_cleared_msg <- renderText(NULL)
  })

  run_lm <- function() {
    req(rv$meta_sample, rv$abundance_sample, input$lm_outcome, input$lm_predictors)

    meta_patient <- rv$meta_sample

    # Outcome
    if (!(input$lm_outcome %in% colnames(meta_patient))) {
      showNotification(paste0("Outcome '", input$lm_outcome, "' not found in metadata."), type = "error")
      return(NULL)
    }
    outcome_raw <- meta_patient[[input$lm_outcome]]
    outcome <- factor(outcome_raw)

    # Sanitize levels
    orig_levels <- levels(outcome)
    safe_levels <- make.names(orig_levels)
    levels(outcome) <- safe_levels
    rv$lm_label_map <- setNames(orig_levels, safe_levels)

    if (nlevels(outcome) < 2) {
      showNotification("Outcome must have ≥2 classes.", type = "error")
      return(NULL)
    }

    # Predictors
    pred_meta <- setdiff(intersect(input$lm_predictors, colnames(rv$meta_sample)), "cluster")
    all_clusters <- colnames(rv$abundance_sample)
    if ("cluster" %in% input$lm_predictors) {
      if (!is.null(input$lm_cluster_subset) && length(input$lm_cluster_subset) > 0) {
        cluster_predictors <- intersect(input$lm_cluster_subset, all_clusters)
      } else {
        cluster_predictors <- all_clusters
      }
    } else {
      cluster_predictors <- character(0)
    }

    predictors_final <- c(pred_meta, cluster_predictors)
    if (length(predictors_final) == 0) {
      showNotification("Select at least one predictor.", type = "error")
      return(NULL)
    }

    # Subset before merge
    meta_sub <- meta_patient %>%
      dplyr::select(patient_ID, !!input$lm_outcome, dplyr::all_of(pred_meta))
    abund_sub <- rv$abundance_sample[, cluster_predictors, drop = FALSE]
    merged <- align_metadata_abundance(meta_sub, abund_sub, notify = showNotification)

    # Filter missingness
    merged <- merged[!is.na(merged[[input$lm_outcome]]), ]
    X <- merged[, predictors_final, drop = FALSE]
    y <- factor(merged[[input$lm_outcome]], levels = levels(outcome))

    n_before <- nrow(merged)
    complete_rows <- complete.cases(data.frame(X, .y = y))
    dropped_ids <- merged$patient_ID[!complete_rows]

    X <- X[complete_rows, , drop = FALSE]
    y <- y[complete_rows]
    merged <- merged[complete_rows, , drop = FALSE]

    n_after <- sum(complete_rows)
    n_dropped <- n_before - n_after
    if (n_after < 10) {
      showNotification("Too few samples after filtering. Model not run.", type = "error")
      return(NULL)
    }

    # Null model baseline
    majority_class <- names(which.max(table(y)))
    null_acc <- mean(y == majority_class)

    # Choose model method
    method <- switch(input$lm_model_type,
      "Logistic Regression" = "multinom",
      "Elastic Net"         = "glmnet",
      "Random Forest"       = "rf"
    )
    cat("Running model type:", method, "\n")

    validation <- input$lm_validation
    split_info <- NULL

    # --- Train/Test split ---
    if (validation == "Train/Test split") {
      set.seed(123)
      train_frac <- input$lm_train_frac %||% 0.7
      idx <- caret::createDataPartition(y, p = train_frac, list = FALSE)
      trainX <- X[idx, , drop = FALSE]
      testX <- X[-idx, , drop = FALSE]
      trainY <- y[idx]
      testY <- y[-idx]

      split_info <- list(
        n_train = length(trainY),
        n_test = length(testY),
        train_counts = table(trainY),
        test_counts = table(testY)
      )

      if (method == "glmnet") {
        trainMat <- model.matrix(~ . - 1, data = trainX)
        testMat <- model.matrix(~ . - 1, data = testX)
        storage.mode(trainMat) <- "double"
        storage.mode(testMat) <- "double"
        alpha_val <- input$lm_alpha %||% 0.5
        model <- tryCatch(
          {
            caret::train(
              x = trainMat, y = trainY,
              method = "glmnet",
              trControl = caret::trainControl(classProbs = TRUE, verboseIter = TRUE),
              tuneGrid = expand.grid(
                alpha = alpha_val,
                lambda = 10^seq(-3, 1, length = 20)
              )
            )
          },
          error = function(e) {
            err_msg <- conditionMessage(e)
            if (grepl("one multinomial or binomial class has 1 or 0 observations", err_msg, ignore.case = TRUE)) {
              stop("Insufficient observations in one or more outcome classes. Ensure each class has at least 8 samples.", call. = FALSE)
            } else if (grepl("fewer than 8", err_msg, ignore.case = TRUE)) {
              stop("One or more outcome classes have fewer than 8 observations. Consider grouping rare classes or using a different outcome variable.", call. = FALSE)
            } else {
              stop(paste0("Model training failed: ", err_msg), call. = FALSE)
            }
          }
        )
        probs <- predict(model, newdata = testMat, type = "prob")
        pred_class <- predict(model, newdata = testMat, type = "raw")
      } else {
        model <- tryCatch(
          {
            caret::train(
              x = trainX, y = trainY,
              method = method,
              trControl = caret::trainControl(classProbs = TRUE, verboseIter = TRUE),
              verbose = if (method == "rf") TRUE else FALSE
            )
          },
          error = function(e) {
            stop(paste0("Model training failed: ", conditionMessage(e)), call. = FALSE)
          }
        )
        probs <- predict(model, newdata = testX, type = "prob")
        pred_class <- predict(model, newdata = testX, type = "raw")
      }

      if (is.vector(probs)) {
        pos <- levels(trainY)[2]
        probs <- data.frame(setNames(list(as.numeric(probs)), pos), check.names = FALSE)
      }
      preds_tbl <- dplyr::bind_cols(
        truth = testY,
        dplyr::rename(probs, !!!setNames(names(probs), paste0(".pred_", names(probs)))),
        pred = pred_class
      )
      names(preds_tbl)[1] <- "obs"
      preds_tbl$obs <- factor(preds_tbl$obs, levels = levels(y))

      # Ensure pred column exists
      if (!"pred" %in% names(preds_tbl)) {
        prob_cols <- grep("^\\.pred_", names(preds_tbl), value = TRUE)
        classes <- gsub("^\\.pred_", "", prob_cols)
        max_idx <- apply(as.matrix(preds_tbl[, prob_cols, drop = FALSE]), 1, which.max)
        preds_tbl$pred <- factor(classes[max_idx], levels = classes)
      }

      return(list(
        model = model, preds = preds_tbl, null_acc = null_acc,
        split_info = split_info,
        details = list(
          samples_before = n_before,
          samples_after = n_after,
          samples_dropped = n_dropped,
          dropped_ids = dropped_ids
        )
      ))
    }

    # --- CV or LOO ---
    ctrl <- caret::trainControl(
      method = switch(validation,
        "k-fold CV"     = "cv",
        "Leave-One-Out" = "LOOCV"
      ),
      number = if (validation == "k-fold CV") input$lm_k else 1,
      classProbs = TRUE,
      savePredictions = "final",
      verboseIter = TRUE
    )

    if (method == "glmnet") {
      Xmat <- model.matrix(~ . - 1, data = X)
      storage.mode(Xmat) <- "double"
      alpha_val <- input$lm_alpha %||% 0.5
      model <- caret::train(
        x = Xmat, y = y,
        method = "glmnet",
        trControl = ctrl,
        tuneGrid = expand.grid(
          alpha = alpha_val,
          lambda = 10^seq(-3, 1, length = 20)
        )
      )
    } else {
      model <- caret::train(
        x = X, y = y,
        method = method,
        trControl = ctrl,
        verbose = if (method == "rf") TRUE else FALSE
      )
    }

    preds <- model$pred
    if (is.null(preds) || !nrow(preds)) {
      showNotification("No predictions available from resampling. Try Train/Test split.", type = "error")
      return(NULL)
    }

    level_names <- levels(y)
    prob_cols <- intersect(level_names, names(preds))
    rename_map <- setNames(prob_cols, paste0(".pred_", prob_cols))
    preds <- dplyr::rename(preds, !!!rename_map)
    preds$obs <- factor(preds$obs, levels = levels(y))

    # Ensure pred column exists
    if (!"pred" %in% names(preds)) {
      prob_cols <- grep("^\\.pred_", names(preds), value = TRUE)
      classes <- gsub("^\\.pred_", "", prob_cols)
      max_idx <- apply(as.matrix(preds[, prob_cols, drop = FALSE]), 1, which.max)
      preds$pred <- factor(classes[max_idx], levels = classes)
    }

    return(list(
      model = model, preds = preds, null_acc = null_acc,
      details = list(
        samples_before = n_before,
        samples_after = n_after,
        samples_dropped = n_dropped,
        dropped_ids = dropped_ids
      )
    ))
  }

  observeEvent(input$run_lm, {
    result <- tryCatch(
      {
        run_lm()
      },
      error = function(e) {
        list(error = TRUE, message = conditionMessage(e))
      }
    )

    if (!is.null(result) && isTRUE(result$error)) {
      lm_state(result)
      showNotification(paste("Classification failed:", result$message), type = "error", duration = 10)
      return()
    }

    res <- result
    req(res)

    preds <- res$preds
    req(preds)

    n_classes <- nlevels(preds$obs)
    roc_plot <- NULL

    if (n_classes == 2) {
      # Binary ROC with pROC
      classes <- levels(preds$obs)
      positive <- classes[2]
      roc_obj <- pROC::roc(
        response = preds$obs,
        predictor = preds[[paste0(".pred_", positive)]],
        levels = rev(classes)
      )
      # Create ROC data frame and sort by fpr to ensure smooth curve
      roc_data <- data.frame(
        fpr = 1 - roc_obj$specificities,
        tpr = roc_obj$sensitivities
      )
      roc_data <- roc_data[order(roc_data$fpr, roc_data$tpr), ]

      roc_plot <- ggplot2::ggplot(roc_data, ggplot2::aes(x = fpr, y = tpr)) +
        ggplot2::geom_step(color = "blue", linewidth = 1, direction = "vh") +
        ggplot2::geom_abline(linetype = "dashed", color = "grey50") +
        ggplot2::labs(
          title = sprintf("Binary ROC Curve (AUC = %.3f)", pROC::auc(roc_obj)),
          x = "False Positive Rate", y = "True Positive Rate"
        ) +
        ggplot2::theme_minimal(base_size = 14) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    } else {
      # Multiclass ROC with yardstick
      long_preds <- preds %>%
        dplyr::select(obs, starts_with(".pred_"))
      roc_curves <- yardstick::roc_curve(long_preds, truth = obs, dplyr::starts_with(".pred_"))
      roc_plot <- ggplot2::ggplot(
        roc_curves,
        ggplot2::aes(x = 1 - specificity, y = sensitivity, color = .level)
      ) +
        ggplot2::geom_path(linewidth = 1) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
        ggplot2::labs(title = "Multiclass ROC Curves (yardstick)", color = "Class") +
        ggplot2::theme_bw(base_size = 14) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    }

    # Store everything in lm_state, including cached ROC plot
    lm_state(c(res, list(roc_plot = roc_plot)))

    output$lm_cleared_msg <- renderText(NULL)
  })

  observeEvent(input$reset_lm, {
    lm_state(NULL)
    showNotification("Logistic Modeling results cleared.", type = "message", duration = 5)
    output$lm_cleared_msg <- renderText("Results cleared. Run Logistic Modeling again to see results here.")
  })

  output$lm_roc_plot <- renderPlot({
    s <- lm_state()
    req(s)

    # Check for error
    if (!is.null(s$error) && isTRUE(s$error)) {
      plot.new()
      text(0.5, 0.5, paste("Error:", s$message), col = "red", cex = 1.2)
      return()
    }

    req(s$roc_plot)
    s$roc_plot
  })

  build_perf_table <- function(res, rv) {
    req(res, res$preds)
    preds <- res$preds
    req(preds)

    n_classes <- nlevels(preds$obs)

    obs <- preds$obs
    pred <- preds$pred
    if (is.null(pred)) {
      pred_cols <- grep("^\\.pred_", names(preds), value = TRUE)
      classes <- gsub("^\\.pred_", "", pred_cols)
      max_idx <- apply(as.matrix(preds[, pred_cols, drop = FALSE]), 1, which.max)
      pred <- factor(classes[max_idx], levels = classes)
    }

    model_acc <- mean(pred == obs, na.rm = TRUE)
    successes <- sum(pred == obs, na.rm = TRUE)
    n <- length(obs)

    binom_p <- tryCatch(
      {
        binom.test(successes, n, p = res$null_acc, alternative = "greater")$p.value
      },
      error = function(e) NA
    )

    majority_safe <- names(which.max(table(obs)))
    majority_orig <- rv$lm_label_map[[majority_safe]]

    # Format binomial p-value with scientific notation if < 0.001
    binom_p_formatted <- if (is.na(binom_p)) {
      NA_character_
    } else if (binom_p < 0.001) {
      formatC(binom_p, format = "e", digits = 2)
    } else {
      as.character(signif(binom_p, 3))
    }

    perf_table <- data.frame(
      Metric = c("Null Accuracy (most frequent class)", "Model Accuracy", "Binomial p-value vs Null"),
      Value = c(
        paste0(round(res$null_acc, 3), " (", majority_orig, ")"),
        round(model_acc, 3),
        binom_p_formatted
      ),
      Interpretation = c(
        "Baseline accuracy when always predicting the most frequent class",
        "Proportion of correctly predicted labels by the model (chosen threshold)",
        "One-sided binomial p-value testing model accuracy > null (baseline)"
      ),
      stringsAsFactors = FALSE
    )

    if (n_classes == 2) {
      classes <- levels(preds$obs)
      positive <- classes[2]
      roc_obj <- pROC::roc(
        response = preds$obs,
        predictor = preds[[paste0(".pred_", positive)]],
        levels = rev(classes)
      )
      auc_val <- as.numeric(pROC::auc(roc_obj))
      perf_table <- rbind(
        perf_table,
        data.frame(
          Metric = "AUC (pROC)", Value = round(auc_val, 3),
          Interpretation = "Area under ROC curve; 0.5=random, 1=perfect; measures discrimination",
          stringsAsFactors = FALSE
        )
      )

      # Add 95% CI for AUC using DeLong (fast, nonparametric)
      auc_ci <- tryCatch(
        {
          ci_raw <- pROC::ci.auc(roc_obj, method = "delong")
          # ci.auc often returns c(lower, auc, upper) or c(lower, upper)
          ci_num <- as.numeric(ci_raw)
          if (length(ci_num) >= 3) {
            lower <- ci_num[1]
            upper <- ci_num[3]
          } else if (length(ci_num) == 2) {
            lower <- ci_num[1]
            upper <- ci_num[2]
          } else {
            lower <- NA_real_
            upper <- NA_real_
          }
          c(lower = lower, upper = upper)
        },
        error = function(e) c(lower = NA_real_, upper = NA_real_)
      )

      if (!is.na(auc_ci[1]) && !is.na(auc_ci[2])) {
        ci_text <- sprintf("%.3f - %.3f", auc_ci[1], auc_ci[2])
      } else {
        ci_text <- NA_character_
      }
      perf_table <- rbind(
        perf_table,
        data.frame(
          Metric = "AUC 95% CI (DeLong)", Value = ci_text,
          Interpretation = "95% confidence interval for AUC (DeLong nonparametric)",
          stringsAsFactors = FALSE
        )
      )

      # Efficient AUC p-value vs chance: Wilcoxon rank-sum on positive vs negative scores
      p_val_auc <- tryCatch(
        {
          scores_pos <- preds[[paste0(".pred_", positive)]][preds$obs == positive]
          neg_class <- classes[1]
          scores_neg <- preds[[paste0(".pred_", positive)]][preds$obs == neg_class]
          # require at least one observation in each group
          if (length(scores_pos) < 1 || length(scores_neg) < 1) {
            return(NA_real_)
          }
          wt <- suppressWarnings(wilcox.test(scores_pos, scores_neg, alternative = "greater"))
          wt$p.value
        },
        error = function(e) NA_real_
      )

      pval_formatted <- if (is.na(p_val_auc)) NA_character_ else if (p_val_auc < 0.001) formatC(p_val_auc, format = "e", digits = 2) else as.character(signif(p_val_auc, 3))
      perf_table <- rbind(
        perf_table,
        data.frame(
          Metric = "AUC p-value (Wilcoxon)", Value = pval_formatted,
          Interpretation = "One-sided Wilcoxon rank test p-value testing AUC > 0.5 (uses predicted scores)",
          stringsAsFactors = FALSE
        )
      )
    } else {
      long_preds <- preds %>%
        dplyr::select(obs, starts_with(".pred_"))

      auc_macro <- yardstick::roc_auc(long_preds,
        truth = obs,
        dplyr::starts_with(".pred_"), estimator = "macro"
      )
      auc_macro_wt <- yardstick::roc_auc(long_preds,
        truth = obs,
        dplyr::starts_with(".pred_"), estimator = "macro_weighted"
      )
      auc_ht <- yardstick::roc_auc(long_preds,
        truth = obs,
        dplyr::starts_with(".pred_"), estimator = "hand_till"
      )

      perf_table <- rbind(
        perf_table,
        data.frame(
          Metric = "Macro AUC (yardstick)", Value = round(auc_macro$.estimate, 3),
          Interpretation = "Average AUC across classes (unweighted)", stringsAsFactors = FALSE
        ),
        data.frame(
          Metric = "Macro-weighted AUC (yardstick)", Value = round(auc_macro_wt$.estimate, 3),
          Interpretation = "Average AUC weighted by class prevalence", stringsAsFactors = FALSE
        ),
        data.frame(
          Metric = "Hand-Till AUC (yardstick)", Value = round(auc_ht$.estimate, 3),
          Interpretation = "Hand-Till multiclass generalization of AUC", stringsAsFactors = FALSE
        )
      )
    }

    if (!is.null(res$split_info)) {
      train_summary <- paste(names(res$split_info$train_counts),
        res$split_info$train_counts,
        collapse = "; "
      )
      test_summary <- paste(names(res$split_info$test_counts),
        res$split_info$test_counts,
        collapse = "; "
      )
      extra_rows <- data.frame(
        Metric = c("Train size", "Test size", "Train breakdown", "Test breakdown"),
        Value = c(res$split_info$n_train, res$split_info$n_test, train_summary, test_summary),
        Interpretation = c(
          "Number of samples used for training",
          "Number of samples used for testing",
          "Class counts in training set",
          "Class counts in test set"
        ),
        stringsAsFactors = FALSE
      )
      perf_table <- rbind(perf_table, extra_rows)
    }

    det <- res$details
    extra_rows2 <- data.frame(
      Metric = c("Samples before filtering", "Samples after filtering", "Samples dropped"),
      Value = c(det$samples_before, det$samples_after, det$samples_dropped),
      Interpretation = c(
        "Count of samples before any filtering",
        "Count of samples after filtering applied",
        "Number of samples removed by filtering"
      ),
      stringsAsFactors = FALSE
    )
    perf_table <- rbind(perf_table, extra_rows2)
    if (length(det$dropped_ids) > 0) {
      perf_table <- rbind(
        perf_table,
        data.frame(
          Metric = "Dropped patient_IDs",
          Value = paste(unique(det$dropped_ids), collapse = ", "),
          Interpretation = "List of patient_IDs removed by filtering",
          stringsAsFactors = FALSE
        )
      )
    }

    perf_table
  }

  output$lm_perf_table <- renderTable({
    res <- lm_state()
    req(res)
    build_perf_table(res, rv)
  })

  output$lm_summary <- renderPrint({
    # res <- lm_results(); req(res)
    res <- lm_state()
    req(res)

    # Check for error
    if (!is.null(res$error) && isTRUE(res$error)) {
      cat("Classification Error:\n\n")
      cat(res$message, "\n\n")
      cat("Suggestions:\n")
      cat("- Ensure each outcome class has at least 8 samples\n")
      cat("- Try a different outcome variable\n")
      cat("- Consider combining rare outcome classes\n")
      cat("- Use fewer or different predictors\n")
      return()
    }

    det <- res$details %||% list(
      samples_before = NA,
      samples_after = NA,
      samples_dropped = NA,
      dropped_ids = character(0)
    )

    cat("Samples before filtering:", det$samples_before, "\n")
    cat("Samples after filtering:", det$samples_after, "\n")
    cat("Samples dropped:", det$samples_dropped, "\n\n")

    if (length(det$dropped_ids) > 0) {
      cat("Dropped patient_IDs:\n")
      print(det$dropped_ids)
      cat("\n")
    }

    # Outcome distribution
    preds <- res$preds
    if (!is.null(preds) && "obs" %in% names(preds)) {
      cat("Outcome distribution (observed):\n")
      print(table(preds$obs))
      cat("\n")
    }

    # Train/Test breakdown if available
    if (!is.null(res$split_info)) {
      cat("Train breakdown:\n")
      print(res$split_info$train_counts)
      cat("\nTest breakdown:\n")
      print(res$split_info$test_counts)
      cat("\n")
    }

    cat("Null accuracy (baseline):", round(res$null_acc, 3), "\n")
  })

  order_features <- function(df, raw_col) {
    if (is.null(df) || !raw_col %in% colnames(df)) {
      return(df)
    }

    # Always work with a plain data.frame
    df <- as.data.frame(df)

    # Flag intercept rows
    intercept_row <- grepl("\\(Intercept\\)", df$Feature, fixed = TRUE)
    df_intercept <- df[intercept_row, , drop = FALSE]
    df_notintercept <- df[!intercept_row, , drop = FALSE]

    # Sort non‑intercepts by absolute raw value
    if (nrow(df_notintercept) > 0) {
      raw_vals_notint <- as.numeric(df_notintercept[[raw_col]])
      df_notintercept <- df_notintercept[order(abs(raw_vals_notint), decreasing = TRUE), , drop = FALSE]
    }

    # Sort intercepts (in case there are multiple)
    if (nrow(df_intercept) > 0) {
      raw_vals_int <- as.numeric(df_intercept[[raw_col]])
      df_intercept <- df_intercept[order(abs(raw_vals_int), decreasing = TRUE), , drop = FALSE]
      df <- rbind(df_intercept, df_notintercept)
    } else {
      df <- df_notintercept
    }

    rownames(df) <- NULL
    df
  }

  # --- Extract coefficients for glm / glmnet ---
  # Coefficients for glmnet, glm, multinom with consistent schema and column order
  coef_table <- reactive({
    s <- lm_state()
    req(s, s$model)

    if (s$model$method == "glmnet") {
      coefs <- tryCatch(
        {
          coef(s$model$finalModel, s = s$model$bestTune$lambda)
        },
        error = function(e) NULL
      )
      if (is.null(coefs)) {
        return(NULL)
      }

      if (is.list(coefs)) {
        # Multiclass glmnet: list of matrices, each matrix rows = predictors
        df_list <- lapply(names(coefs), function(cls) {
          mat <- as.matrix(coefs[[cls]])
          data.frame(
            Class = cls, # outcome category
            Feature = rownames(mat), # predictor
            StandardizedCoef = as.numeric(mat),
            stringsAsFactors = FALSE
          )
        })
        df <- do.call(rbind, df_list)
      } else {
        # Binary/continuous glmnet: single matrix, rows = predictors
        mat <- as.matrix(coefs)
        df <- data.frame(
          Feature = rownames(mat),
          StandardizedCoef = as.numeric(mat),
          stringsAsFactors = FALSE
        )
      }

      # Remove backticks from feature names (added by model.matrix for numeric column names)
      df$Feature <- gsub("`", "", df$Feature, fixed = TRUE)

      # Calculate original-scale coefficients for interpretation
      # Formula: beta_original = beta_standardized * (sd_y / sd_x)
      df$OriginalCoef <- df$StandardizedCoef

      if (!is.null(s$model$y_scale_params) || !is.null(s$model$X_scale_params)) {
        y_sd <- if (!is.null(s$model$y_scale_params)) s$model$y_scale_params$sd else 1
        X_scale <- s$model$X_scale_params

        for (i in seq_len(nrow(df))) {
          feat <- df$Feature[i]
          x_sd <- if (!is.null(X_scale[[feat]])) X_scale[[feat]]$sd else 1
          # Rescale to original units
          df$OriginalCoef[i] <- df$StandardizedCoef[i] * (y_sd / x_sd)
        }
      }

      # Enforce column order
      if ("Class" %in% names(df)) {
        df <- df[, c("Class", "Feature", "StandardizedCoef", "OriginalCoef")]
      } else {
        df <- df[, c("Feature", "StandardizedCoef", "OriginalCoef")]
      }
      df <- order_features(df, "StandardizedCoef")
      rownames(df) <- NULL
      df
    } else if (s$model$method %in% c("glm", "multinom")) {
      coefs <- coef(s$model$finalModel)

      if (is.matrix(coefs)) {
        # Multiclass multinom: rows = outcome classes, columns = predictors
        df <- as.data.frame(coefs)
        df$Class <- rownames(coefs) # outcome categories from row names

        df_long <- tidyr::pivot_longer(
          df,
          cols = setdiff(names(df), "Class"),
          names_to = "Feature", # predictors from column names
          values_to = "OriginalCoef"
        )

        # No standardization for glm/multinom, so standardized = original
        df_long$StandardizedCoef <- df_long$OriginalCoef
        df <- df_long[, c("Class", "Feature", "StandardizedCoef", "OriginalCoef")]
      } else {
        # Binary logistic regression: named vector of coefficients (predictors)
        df <- data.frame(
          Feature = names(coefs),
          OriginalCoef = unname(coefs),
          stringsAsFactors = FALSE
        )
        df$StandardizedCoef <- df$OriginalCoef
        df <- df[, c("Feature", "StandardizedCoef", "OriginalCoef")]
      }
      # Remove backticks from feature names
      df$Feature <- gsub("`", "", df$Feature, fixed = TRUE)
      df <- order_features(df, "StandardizedCoef")
      rownames(df) <- NULL
      df
    } else {
      NULL
    }
  })

  # --- Extract variable importance for random forest ---
  rf_importance <- reactive({
    s <- lm_state()
    req(s, s$model)

    if (s$model$method == "rf") {
      imp_scaled <- caret::varImp(s$model, scale = TRUE)$importance
      imp_scaled$Feature <- rownames(imp_scaled)
      colnames(imp_scaled)[1] <- "ScaledImportance"

      rf_model <- s$model$finalModel
      imp_raw <- randomForest::importance(rf_model)
      if ("%IncMSE" %in% colnames(imp_raw)) {
        raw_vals <- imp_raw[, "%IncMSE"]
      } else {
        raw_vals <- imp_raw[, 1]
      }
      imp_raw_df <- data.frame(
        Feature = rownames(imp_raw),
        RawImportance = raw_vals,
        stringsAsFactors = FALSE
      )

      df <- dplyr::left_join(imp_scaled, imp_raw_df, by = "Feature")
      df <- df[, c("Feature", "RawImportance", "ScaledImportance")]
      df <- order_features(df, "RawImportance")
      rownames(df) <- NULL
      df
    } else {
      NULL
    }
  })

  output$lm_features <- renderTable(
    {
      s <- lm_state()
      req(s, s$model)

      if (s$model$method %in% c("glm", "glmnet", "multinom")) {
        coef_table()
      } else if (s$model$method == "rf") {
        rf_importance()
      } else {
        data.frame(Message = "Feature importance not available for this model type.")
      }
    },
    digits = 3
  )

  output$hasLMResults <- reactive({
    res <- lm_state()
    !is.null(res) && !is.null(res$preds) && nrow(res$preds) > 0
  })
  outputOptions(output, "hasLMResults", suspendWhenHidden = FALSE)

  output$export_lm_zip <- downloadHandler(
    filename = function() {
      subset_id <- rv$subset_id %||% "000000000"
      model <- gsub("\\s+", "_", tolower(input$lm_model_type %||% "model"))
      validation <- gsub("\\s+", "_", tolower(input$lm_validation %||% "validation"))
      outcome <- gsub("\\s+", "_", tolower(input$lm_outcome %||% "outcome"))
      paste0(model, "_", validation, "_", outcome, "_", subset_id, ".zip")
    },
    content = function(file) {
      subset_id <- rv$subset_id %||% "000000000"
      s <- lm_state()
      req(s, s$model)

      # Build filename prefix from zip name components
      model <- gsub("\\s+", "_", tolower(input$lm_model_type %||% "model"))
      validation <- gsub("\\s+", "_", tolower(input$lm_validation %||% "validation"))
      outcome <- gsub("\\s+", "_", tolower(input$lm_outcome %||% "outcome"))
      prefix <- paste0(model, "_", validation, "_", outcome)

      tmpdir <- tempdir()
      files <- c()

      # 1. Model Summary (key-value CSV)
      summary_file <- file.path(tmpdir, paste0(prefix, "_model_summary_", subset_id, ".csv"))
      summary_kv <- data.frame(
        Key = c(
          "Model type", "Outcome variable", "Predictors",
          "Validation strategy", "Null accuracy",
          "Samples before filtering", "Samples after filtering", "Samples dropped"
        ),
        Value = c(
          input$lm_model_type,
          input$lm_outcome,
          paste(input$lm_predictors, collapse = "; "),
          input$lm_validation,
          sprintf("%.3f", s$null_acc %||% NA),
          s$details$samples_before %||% NA,
          s$details$samples_after %||% NA,
          s$details$samples_dropped %||% NA
        ),
        stringsAsFactors = FALSE
      )
      write.csv(summary_kv, summary_file, row.names = FALSE)
      files <- c(files, summary_file)

      # 2. Performance Metrics (exactly as in UI)
      perf_file <- file.path(tmpdir, paste0(prefix, "_performance_metrics_", subset_id, ".csv"))
      perf_df <- tryCatch(
        {
          build_perf_table(s, rv)
        },
        error = function(e) data.frame(Message = "Error extracting performance metrics")
      )
      write.csv(perf_df, perf_file, row.names = FALSE)
      files <- c(files, perf_file)

      # 3. Model Features (coefficients / importance with enforced column order)
      feat_file <- file.path(tmpdir, paste0(prefix, "_model_features_", subset_id, ".csv"))
      feat_df <- NULL
      if (s$model$method %in% c("glm", "glmnet", "multinom")) {
        feat_df <- coef_table()
      } else if (s$model$method == "rf") {
        feat_df <- rf_importance()
      }

      if (!is.null(feat_df) && nrow(feat_df) > 0) {
        # Enforce column order
        if (all(c("Class", "Feature", "StandardizedCoef", "OriginalCoef") %in% names(feat_df))) {
          feat_df <- feat_df[, c("Class", "Feature", "StandardizedCoef", "OriginalCoef")]
        } else if (all(c("Feature", "StandardizedCoef", "OriginalCoef") %in% names(feat_df))) {
          feat_df <- feat_df[, c("Feature", "StandardizedCoef", "OriginalCoef")]
        } else if (all(c("Feature", "RawImportance", "ScaledImportance") %in% names(feat_df))) {
          feat_df <- feat_df[, c("Feature", "RawImportance", "ScaledImportance")]
        }
        write.csv(feat_df, feat_file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No feature importance available"),
          feat_file,
          row.names = FALSE
        )
      }
      files <- c(files, feat_file)

      # 4. ROC Curve (PDF) - wider for multi-class
      roc_file <- file.path(tmpdir, paste0("roc_curve_", subset_id, ".pdf"))
      if (!is.null(s$roc_plot)) {
        # Determine if multi-class to adjust width
        n_classes <- nlevels(s$preds$obs)
        plot_width <- if (n_classes > 2) 10.5 else 6 # 1.75x wider for multi-class
        ggsave(roc_file,
          plot = s$roc_plot,
          device = if (capabilities("cairo")) cairo_pdf else pdf,
          width = plot_width, height = 6, units = "in"
        )
      } else {
        pdf(roc_file)
        plot.new()
        text(0.5, 0.5, "ROC plot not available")
        dev.off()
      }
      files <- c(files, roc_file)

      # Bundle into zip
      zip::zip(zipfile = file, files = files, mode = "cherry-pick")
    },
    contentType = "application/zip"
  )

  # ========== REGRESSION TAB LOGIC ==========

  run_reg <- function() {
    req(rv$meta_sample, rv$abundance_sample, input$reg_outcome, input$reg_predictors)

    meta_patient <- rv$meta_sample

    # Outcome - must be numeric
    if (!(input$reg_outcome %in% colnames(meta_patient))) {
      showNotification(paste0("Outcome '", input$reg_outcome, "' not found in metadata."), type = "error")
      return(NULL)
    }
    outcome_raw <- meta_patient[[input$reg_outcome]]
    outcome <- suppressWarnings(as.numeric(outcome_raw))

    if (all(is.na(outcome))) {
      showNotification("Outcome variable could not be converted to numeric.", type = "error")
      return(NULL)
    }

    # Predictors
    pred_meta <- setdiff(intersect(input$reg_predictors, colnames(rv$meta_sample)), "cluster")
    all_clusters <- colnames(rv$abundance_sample)
    if ("cluster" %in% input$reg_predictors) {
      if (!is.null(input$reg_cluster_subset) && length(input$reg_cluster_subset) > 0) {
        cluster_predictors <- intersect(input$reg_cluster_subset, all_clusters)
      } else {
        cluster_predictors <- all_clusters
      }
    } else {
      cluster_predictors <- character(0)
    }

    predictors_final <- c(pred_meta, cluster_predictors)
    if (length(predictors_final) == 0) {
      showNotification("Select at least one predictor.", type = "error")
      return(NULL)
    }

    # Subset before merge
    meta_sub <- meta_patient %>%
      dplyr::select(patient_ID, !!input$reg_outcome, dplyr::all_of(pred_meta))
    abund_sub <- rv$abundance_sample[, cluster_predictors, drop = FALSE]
    merged <- align_metadata_abundance(meta_sub, abund_sub, notify = showNotification)

    # Filter missingness
    merged <- merged[!is.na(merged[[input$reg_outcome]]), ]
    X <- merged[, predictors_final, drop = FALSE]
    y <- as.numeric(merged[[input$reg_outcome]])

    n_before <- nrow(merged)
    complete_rows <- complete.cases(data.frame(X, .y = y))
    dropped_ids <- merged$patient_ID[!complete_rows]

    X <- X[complete_rows, , drop = FALSE]
    y <- y[complete_rows]
    merged <- merged[complete_rows, , drop = FALSE]

    n_after <- sum(complete_rows)
    n_dropped <- n_before - n_after
    if (n_after < 10) {
      showNotification("Too few samples after filtering. Model not run.", type = "error")
      return(NULL)
    }

    # Null model baseline (intercept-only)
    null_rmse <- sqrt(mean((y - mean(y))^2))

    # Choose model method
    method <- switch(input$reg_model_type,
      "Linear Regression"  = "lm",
      "Ridge Regression"   = "glmnet",
      "Elastic Net"        = "glmnet",
      "Random Forest"      = "rf"
    )
    cat("Running regression model type:", method, "\n")

    validation <- input$reg_validation
    split_info <- NULL

    # --- Train/Test split ---
    if (validation == "Train/Test split") {
      set.seed(123)
      train_frac <- input$reg_train_frac %||% 0.7
      idx <- sample(seq_len(length(y)), size = floor(train_frac * length(y)))
      trainX <- X[idx, , drop = FALSE]
      testX <- X[-idx, , drop = FALSE]
      trainY_orig <- y[idx]
      testY_orig <- y[-idx]

      split_info <- list(
        n_train = length(trainY_orig),
        n_test  = length(testY_orig)
      )

      if (method == "glmnet") {
        trainMat <- model.matrix(~ . - 1, data = trainX)
        testMat <- model.matrix(~ . - 1, data = testX)
        storage.mode(trainMat) <- "double"
        storage.mode(testMat) <- "double"

        # Standardize continuous predictors explicitly
        X_scale_params <- list()
        for (j in seq_len(ncol(trainMat))) {
          col <- trainMat[, j]
          unique_vals <- unique(col)
          if (length(unique_vals) > 2) {
            col_mean <- mean(col, na.rm = TRUE)
            col_sd <- sd(col, na.rm = TRUE)
            if (col_sd > 0) {
              trainMat[, j] <- (col - col_mean) / col_sd
              testMat[, j] <- (testMat[, j] - col_mean) / col_sd
              X_scale_params[[colnames(trainMat)[j]]] <- list(mean = col_mean, sd = col_sd)
            }
          }
        }

        # Standardize outcome
        trainY_mean <- mean(trainY_orig, na.rm = TRUE)
        trainY_sd <- sd(trainY_orig, na.rm = TRUE)
        trainY <- if (trainY_sd > 0) (trainY_orig - trainY_mean) / trainY_sd else trainY_orig
        y_scale_params <- list(mean = trainY_mean, sd = trainY_sd)

        alpha_val <- if (input$reg_model_type == "Ridge Regression") 0 else (input$reg_alpha %||% 0.5)
        model <- tryCatch(
          {
            caret::train(
              x = trainMat, y = trainY,
              method = "glmnet",
              trControl = caret::trainControl(verboseIter = TRUE),
              tuneGrid = expand.grid(
                alpha = alpha_val,
                lambda = 10^seq(-3, 1, length = 20)
              )
            )
          },
          error = function(e) {
            stop(paste0("Regression model training failed: ", conditionMessage(e)), call. = FALSE)
          }
        )

        # Store scaling params for coefficient interpretation
        model$X_scale_params <- X_scale_params
        model$y_scale_params <- y_scale_params

        preds_scaled <- predict(model, newdata = testMat)
        # Rescale predictions back to original scale
        preds <- preds_scaled * y_scale_params$sd + y_scale_params$mean
      } else if (method == "rf") {
        model <- tryCatch(
          {
            caret::train(
              x = trainX, y = trainY_orig,
              method = method,
              trControl = caret::trainControl(verboseIter = TRUE),
              verbose = TRUE
            )
          },
          error = function(e) {
            stop(paste0("Regression model training failed: ", conditionMessage(e)), call. = FALSE)
          }
        )
        preds <- predict(model, newdata = testX)
      } else {
        # Linear regression - apply same standardization for consistency
        trainMat <- model.matrix(~ . - 1, data = trainX)
        testMat <- model.matrix(~ . - 1, data = testX)
        storage.mode(trainMat) <- "double"
        storage.mode(testMat) <- "double"

        # Standardize continuous predictors explicitly
        X_scale_params <- list()
        for (j in seq_len(ncol(trainMat))) {
          col <- trainMat[, j]
          unique_vals <- unique(col)
          if (length(unique_vals) > 2) {
            col_mean <- mean(col, na.rm = TRUE)
            col_sd <- sd(col, na.rm = TRUE)
            if (col_sd > 0) {
              trainMat[, j] <- (col - col_mean) / col_sd
              testMat[, j] <- (testMat[, j] - col_mean) / col_sd
              X_scale_params[[colnames(trainMat)[j]]] <- list(mean = col_mean, sd = col_sd)
            }
          }
        }

        # Standardize outcome
        trainY_mean <- mean(trainY_orig, na.rm = TRUE)
        trainY_sd <- sd(trainY_orig, na.rm = TRUE)
        trainY <- if (trainY_sd > 0) (trainY_orig - trainY_mean) / trainY_sd else trainY_orig
        y_scale_params <- list(mean = trainY_mean, sd = trainY_sd)

        model <- tryCatch(
          {
            caret::train(
              x = trainMat, y = trainY,
              method = method,
              trControl = caret::trainControl(verboseIter = TRUE)
            )
          },
          error = function(e) {
            stop(paste0("Regression model training failed: ", conditionMessage(e)), call. = FALSE)
          }
        )

        # Store scaling params for coefficient interpretation
        model$X_scale_params <- X_scale_params
        model$y_scale_params <- y_scale_params

        preds_scaled <- predict(model, newdata = testMat)
        # Rescale predictions back to original scale
        preds <- preds_scaled * y_scale_params$sd + y_scale_params$mean
      }

      preds_tbl <- data.frame(
        obs = testY_orig,
        pred = preds
      )

      return(list(
        model = model, preds = preds_tbl, null_rmse = null_rmse,
        split_info = split_info,
        details = list(
          samples_before = n_before,
          samples_after = n_after,
          samples_dropped = n_dropped,
          dropped_ids = dropped_ids
        )
      ))
    }

    # --- CV or LOO ---
    ctrl <- caret::trainControl(
      method = switch(validation,
        "k-fold CV"     = "cv",
        "Leave-One-Out" = "LOOCV"
      ),
      number = if (validation == "k-fold CV") input$reg_k else 1,
      savePredictions = "final",
      verboseIter = TRUE,
      allowParallel = FALSE
    )

    if (method == "glmnet") {
      Xmat <- model.matrix(~ . - 1, data = X)
      storage.mode(Xmat) <- "double"

      # Standardize continuous predictors explicitly
      X_scale_params <- list()
      for (j in seq_len(ncol(Xmat))) {
        col <- Xmat[, j]
        unique_vals <- unique(col)
        if (length(unique_vals) > 2) {
          col_mean <- mean(col, na.rm = TRUE)
          col_sd <- sd(col, na.rm = TRUE)
          if (col_sd > 0) {
            Xmat[, j] <- (col - col_mean) / col_sd
            X_scale_params[[colnames(Xmat)[j]]] <- list(mean = col_mean, sd = col_sd)
          }
        }
      }

      # Standardize outcome
      y_mean <- mean(y, na.rm = TRUE)
      y_sd <- sd(y, na.rm = TRUE)
      y_scaled <- if (y_sd > 0) (y - y_mean) / y_sd else y
      y_scale_params <- list(mean = y_mean, sd = y_sd)

      alpha_val <- if (input$reg_model_type == "Ridge Regression") 0 else (input$reg_alpha %||% 0.5)
      model <- caret::train(
        x = Xmat, y = y_scaled,
        method = "glmnet",
        trControl = ctrl,
        tuneGrid = expand.grid(
          alpha = alpha_val,
          lambda = 10^seq(-3, 1, length = 20)
        )
      )

      # Store scaling params for coefficient interpretation
      model$X_scale_params <- X_scale_params
      model$y_scale_params <- y_scale_params
    } else if (method == "rf") {
      model <- caret::train(
        x = X, y = y,
        method = method,
        trControl = ctrl,
        verbose = TRUE
      )
    } else {
      # Linear regression - apply same standardization for consistency
      Xmat <- model.matrix(~ . - 1, data = X)
      storage.mode(Xmat) <- "double"

      # Standardize continuous predictors explicitly
      X_scale_params <- list()
      for (j in seq_len(ncol(Xmat))) {
        col <- Xmat[, j]
        unique_vals <- unique(col)
        if (length(unique_vals) > 2) {
          col_mean <- mean(col, na.rm = TRUE)
          col_sd <- sd(col, na.rm = TRUE)
          if (col_sd > 0) {
            Xmat[, j] <- (col - col_mean) / col_sd
            X_scale_params[[colnames(Xmat)[j]]] <- list(mean = col_mean, sd = col_sd)
          }
        }
      }

      # Standardize outcome
      y_mean <- mean(y, na.rm = TRUE)
      y_sd <- sd(y, na.rm = TRUE)
      y_scaled <- if (y_sd > 0) (y - y_mean) / y_sd else y
      y_scale_params <- list(mean = y_mean, sd = y_sd)

      model <- caret::train(
        x = Xmat, y = y_scaled,
        method = method,
        trControl = ctrl
      )

      # Store scaling params for coefficient interpretation
      model$X_scale_params <- X_scale_params
      model$y_scale_params <- y_scale_params
    }

    preds <- model$pred
    if (is.null(preds) || !nrow(preds)) {
      showNotification("No predictions available from resampling. Try Train/Test split.", type = "error")
      return(NULL)
    }

    # Rescale predictions back to original scale for glmnet and lm
    if (method %in% c("glmnet", "lm")) {
      preds$pred <- preds$pred * y_scale_params$sd + y_scale_params$mean
      preds$obs <- preds$obs * y_scale_params$sd + y_scale_params$mean
    }

    preds_tbl <- data.frame(
      obs = preds$obs,
      pred = preds$pred
    )

    return(list(
      model = model, preds = preds_tbl, null_rmse = null_rmse,
      details = list(
        samples_before = n_before,
        samples_after = n_after,
        samples_dropped = n_dropped,
        dropped_ids = dropped_ids
      )
    ))
  }

  observeEvent(input$run_reg, {
    result <- tryCatch(
      {
        run_reg()
      },
      error = function(e) {
        list(error = TRUE, message = conditionMessage(e))
      }
    )

    if (!is.null(result) && isTRUE(result$error)) {
      reg_state(result)
      showNotification(paste("Regression failed:", result$message), type = "error", duration = 10)
      return()
    }

    res <- result
    req(res)

    preds <- res$preds
    req(preds)

    # Create diagnostic plots
    obs_pred_plot <- ggplot2::ggplot(preds, ggplot2::aes(x = obs, y = pred)) +
      ggplot2::geom_point(alpha = 0.6, size = 3) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(
        title = "Observed vs Predicted",
        x = "Observed", y = "Predicted"
      ) +
      ggplot2::theme_minimal(base_size = 14)

    preds$residual <- preds$obs - preds$pred
    residual_plot <- ggplot2::ggplot(preds, ggplot2::aes(x = pred, y = residual)) +
      ggplot2::geom_point(alpha = 0.6, size = 3) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::labs(
        title = "Residuals vs Fitted",
        x = "Fitted Values", y = "Residuals"
      ) +
      ggplot2::theme_minimal(base_size = 14)

    # Store everything in reg_state
    reg_state(c(res, list(obs_pred_plot = obs_pred_plot, residual_plot = residual_plot)))

    output$reg_cleared_msg <- renderText(NULL)
  })

  observeEvent(input$reset_reg, {
    reg_state(NULL)
    showNotification("Regression results cleared.", type = "message", duration = 5)
    output$reg_cleared_msg <- renderText("Results cleared. Run Regression again to see results here.")
  })

  output$reg_obs_pred_plot <- renderPlot({
    s <- reg_state()
    req(s)

    # Check for error
    if (!is.null(s$error) && isTRUE(s$error)) {
      plot.new()
      text(0.5, 0.5, paste("Error:", s$message), col = "red", cex = 1.2)
      return()
    }

    req(s$obs_pred_plot)
    s$obs_pred_plot
  })

  output$reg_residual_plot <- renderPlot({
    s <- reg_state()
    req(s)

    # Check for error
    if (!is.null(s$error) && isTRUE(s$error)) {
      plot.new()
      text(0.5, 0.5, paste("Error:", s$message), col = "red", cex = 1.2)
      return()
    }

    req(s$residual_plot)
    s$residual_plot
  })

  build_reg_perf_table <- function(res) {
    req(res, res$preds)
    preds <- res$preds
    req(preds)

    obs <- preds$obs
    pred <- preds$pred

    # Calculate metrics
    rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
    mae <- mean(abs(obs - pred), na.rm = TRUE)

    # R-squared
    ss_res <- sum((obs - pred)^2, na.rm = TRUE)
    ss_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
    r_squared <- 1 - (ss_res / ss_tot)

    # Correlation
    cor_val <- cor(obs, pred, use = "complete.obs")

    perf_table <- data.frame(
      Metric = c("Null RMSE (intercept-only)", "Model RMSE", "MAE", "R-squared", "Correlation"),
      Value = c(
        sprintf("%.3f", res$null_rmse),
        sprintf("%.3f", rmse),
        sprintf("%.3f", mae),
        sprintf("%.3f", r_squared),
        sprintf("%.3f", cor_val)
      ),
      stringsAsFactors = FALSE
    )

    if (!is.null(res$split_info)) {
      extra_rows <- data.frame(
        Metric = c("Train size", "Test size"),
        Value = c(res$split_info$n_train, res$split_info$n_test),
        stringsAsFactors = FALSE
      )
      perf_table <- rbind(perf_table, extra_rows)
    }

    det <- res$details
    extra_rows2 <- data.frame(
      Metric = c("Samples before filtering", "Samples after filtering", "Samples dropped"),
      Value = c(det$samples_before, det$samples_after, det$samples_dropped),
      stringsAsFactors = FALSE
    )
    perf_table <- rbind(perf_table, extra_rows2)
    if (length(det$dropped_ids) > 0) {
      perf_table <- rbind(
        perf_table,
        data.frame(
          Metric = "Dropped patient_IDs",
          Value = paste(unique(det$dropped_ids), collapse = ", "),
          stringsAsFactors = FALSE
        )
      )
    }

    perf_table
  }

  output$reg_perf_table <- renderTable({
    res <- reg_state()
    req(res)
    build_reg_perf_table(res)
  })

  output$reg_summary <- renderPrint({
    res <- reg_state()
    req(res)

    # Check for error
    if (!is.null(res$error) && isTRUE(res$error)) {
      cat("Regression Error:\n\n")
      cat(res$message, "\n\n")
      cat("Suggestions:\n")
      cat("- Check that your outcome variable has sufficient variation\n")
      cat("- Ensure you have enough samples for model fitting\n")
      cat("- Try a different model type or predictor set\n")
      return()
    }

    det <- res$details %||% list(
      samples_before = NA,
      samples_after = NA,
      samples_dropped = NA,
      dropped_ids = character(0)
    )

    cat("Samples before filtering:", det$samples_before, "\n")
    cat("Samples after filtering:", det$samples_after, "\n")
    cat("Samples dropped:", det$samples_dropped, "\n\n")

    if (length(det$dropped_ids) > 0) {
      cat("Dropped patient_IDs:\n")
      print(det$dropped_ids)
      cat("\n")
    }

    # Outcome distribution
    preds <- res$preds
    if (!is.null(preds) && "obs" %in% names(preds)) {
      cat("Outcome summary (observed):\n")
      print(summary(preds$obs))
      cat("\n")
    }

    # Train/Test breakdown if available
    if (!is.null(res$split_info)) {
      cat("Train size:", res$split_info$n_train, "\n")
      cat("Test size:", res$split_info$n_test, "\n\n")
    }

    cat("Null RMSE (baseline):", round(res$null_rmse, 3), "\n")
  })

  # Coefficients for regression models
  reg_coef_table <- reactive({
    s <- reg_state()
    req(s, s$model)

    if (s$model$method == "glmnet") {
      coefs <- tryCatch(
        {
          coef(s$model$finalModel, s = s$model$bestTune$lambda)
        },
        error = function(e) NULL
      )
      if (is.null(coefs)) {
        return(NULL)
      }

      mat <- as.matrix(coefs)
      df <- data.frame(
        Feature = rownames(mat),
        StandardizedCoef = as.numeric(mat),
        stringsAsFactors = FALSE
      )

      # Remove backticks and backslashes from feature names (added by model.matrix for numeric column names)
      df$Feature <- gsub("`", "", df$Feature, fixed = TRUE)
      df$Feature <- gsub("\\\\", "", df$Feature)

      # Calculate original-scale coefficients for interpretation
      df$OriginalCoef <- df$StandardizedCoef

      if (!is.null(s$model$y_scale_params) || !is.null(s$model$X_scale_params)) {
        y_sd <- if (!is.null(s$model$y_scale_params)) s$model$y_scale_params$sd else 1
        X_scale <- s$model$X_scale_params

        cat("Rescaling coefficients:\n")
        cat("  y_sd:", y_sd, "\n")
        cat("  X_scale params available for:", if (!is.null(X_scale)) names(X_scale) else "NONE", "\n")
        cat("  Features in df:", df$Feature, "\n")

        for (i in seq_len(nrow(df))) {
          feat <- df$Feature[i]
          # Try matching with backticks if exact match fails
          feat_with_backticks <- paste0("`", feat, "`")
          x_sd <- if (!is.null(X_scale[[feat]])) {
            X_scale[[feat]]$sd
          } else if (!is.null(X_scale[[feat_with_backticks]])) {
            X_scale[[feat_with_backticks]]$sd
          } else {
            1
          }
          # Rescale to original units
          df$OriginalCoef[i] <- df$StandardizedCoef[i] * (y_sd / x_sd)
          if (i <= 3) { # Debug first few
            cat("  Feature", feat, ": x_sd =", x_sd, ", rescaled from", df$StandardizedCoef[i], "to", df$OriginalCoef[i], "\n")
          }
        }
      } else {
        cat(
          "No scaling params found - y_scale_params:", !is.null(s$model$y_scale_params),
          "X_scale_params:", !is.null(s$model$X_scale_params), "\n"
        )
      }

      df <- df[, c("Feature", "StandardizedCoef", "OriginalCoef")]
      df <- order_features(df, "StandardizedCoef")
      rownames(df) <- NULL
      df
    } else if (s$model$method == "lm") {
      coefs <- coef(s$model$finalModel)
      df <- data.frame(
        Feature = names(coefs),
        StandardizedCoef = unname(coefs),
        stringsAsFactors = FALSE
      )

      # Remove backticks and backslashes from feature names (added by model.matrix for numeric column names)
      df$Feature <- gsub("`", "", df$Feature, fixed = TRUE)
      df$Feature <- gsub("\\\\", "", df$Feature)

      # Calculate original-scale coefficients for interpretation
      df$OriginalCoef <- df$StandardizedCoef

      if (!is.null(s$model$y_scale_params) || !is.null(s$model$X_scale_params)) {
        y_sd <- if (!is.null(s$model$y_scale_params)) s$model$y_scale_params$sd else 1
        X_scale <- s$model$X_scale_params

        cat("Rescaling lm coefficients:\n")
        cat("  y_sd:", y_sd, "\n")
        cat("  X_scale params available for:", if (!is.null(X_scale)) names(X_scale) else "NONE", "\n")

        for (i in seq_len(nrow(df))) {
          feat <- df$Feature[i]
          # Try matching with backticks if exact match fails
          feat_with_backticks <- paste0("`", feat, "`")
          x_sd <- if (!is.null(X_scale[[feat]])) {
            X_scale[[feat]]$sd
          } else if (!is.null(X_scale[[feat_with_backticks]])) {
            X_scale[[feat_with_backticks]]$sd
          } else {
            1
          }
          # Rescale to original units
          df$OriginalCoef[i] <- df$StandardizedCoef[i] * (y_sd / x_sd)
          if (i <= 3) { # Debug first few
            cat("  Feature", feat, ": x_sd =", x_sd, ", rescaled from", df$StandardizedCoef[i], "to", df$OriginalCoef[i], "\n")
          }
        }
      }

      df <- df[, c("Feature", "StandardizedCoef", "OriginalCoef")]
      df <- order_features(df, "StandardizedCoef")
      rownames(df) <- NULL
      df
    } else {
      NULL
    }
  })

  # Variable importance for random forest regression
  reg_rf_importance <- reactive({
    s <- reg_state()
    req(s, s$model)

    if (s$model$method == "rf") {
      imp_scaled <- caret::varImp(s$model, scale = TRUE)$importance
      imp_scaled$Feature <- rownames(imp_scaled)
      colnames(imp_scaled)[1] <- "ScaledImportance"

      rf_model <- s$model$finalModel
      imp_raw <- randomForest::importance(rf_model)
      if ("%IncMSE" %in% colnames(imp_raw)) {
        raw_vals <- imp_raw[, "%IncMSE"]
      } else {
        raw_vals <- imp_raw[, 1]
      }
      imp_raw_df <- data.frame(
        Feature = rownames(imp_raw),
        RawImportance = raw_vals,
        stringsAsFactors = FALSE
      )

      df <- dplyr::left_join(imp_scaled, imp_raw_df, by = "Feature")
      df <- df[, c("Feature", "RawImportance", "ScaledImportance")]
      df <- order_features(df, "RawImportance")
      rownames(df) <- NULL
      df
    } else {
      NULL
    }
  })

  output$reg_features <- renderTable(
    {
      s <- reg_state()
      req(s, s$model)

      if (s$model$method %in% c("lm", "glmnet")) {
        reg_coef_table()
      } else if (s$model$method == "rf") {
        reg_rf_importance()
      } else {
        data.frame(Message = "Feature importance not available for this model type.")
      }
    },
    digits = 3
  )

  output$hasRegResults <- reactive({
    res <- reg_state()
    !is.null(res) && !is.null(res$preds) && nrow(res$preds) > 0
  })
  outputOptions(output, "hasRegResults", suspendWhenHidden = FALSE)

  output$export_reg_zip <- downloadHandler(
    filename = function() {
      subset_id <- rv$subset_id %||% "000000000"
      model <- gsub("\\s+", "_", tolower(input$reg_model_type %||% "model"))
      validation <- gsub("\\s+", "_", tolower(input$reg_validation %||% "validation"))
      outcome <- gsub("\\s+", "_", tolower(input$reg_outcome %||% "outcome"))
      paste0(model, "_", validation, "_", outcome, "_", subset_id, ".zip")
    },
    content = function(file) {
      subset_id <- rv$subset_id %||% "000000000"
      s <- reg_state()
      req(s, s$model)

      # Build filename prefix from zip name components
      model <- gsub("\\s+", "_", tolower(input$reg_model_type %||% "model"))
      validation <- gsub("\\s+", "_", tolower(input$reg_validation %||% "validation"))
      outcome <- gsub("\\s+", "_", tolower(input$reg_outcome %||% "outcome"))
      prefix <- paste0(model, "_", validation, "_", outcome)

      tmpdir <- tempdir()
      files <- c()

      # 1. Model Summary (key-value CSV)
      summary_file <- file.path(tmpdir, paste0(prefix, "_model_summary_", subset_id, ".csv"))
      summary_kv <- data.frame(
        Key = c(
          "Model type", "Outcome variable", "Predictors",
          "Validation strategy", "Null RMSE",
          "Samples before filtering", "Samples after filtering", "Samples dropped"
        ),
        Value = c(
          input$reg_model_type,
          input$reg_outcome,
          paste(input$reg_predictors, collapse = "; "),
          input$reg_validation,
          sprintf("%.3f", s$null_rmse %||% NA),
          s$details$samples_before %||% NA,
          s$details$samples_after %||% NA,
          s$details$samples_dropped %||% NA
        ),
        stringsAsFactors = FALSE
      )
      write.csv(summary_kv, summary_file, row.names = FALSE)
      files <- c(files, summary_file)

      # 2. Performance Metrics
      perf_file <- file.path(tmpdir, paste0(prefix, "_performance_metrics_", subset_id, ".csv"))
      perf_df <- tryCatch(
        {
          build_reg_perf_table(s)
        },
        error = function(e) data.frame(Message = "Error extracting performance metrics")
      )
      write.csv(perf_df, perf_file, row.names = FALSE)
      files <- c(files, perf_file)

      # 3. Model Features
      feat_file <- file.path(tmpdir, paste0(prefix, "_model_features_", subset_id, ".csv"))
      feat_df <- NULL
      if (s$model$method %in% c("lm", "glmnet")) {
        feat_df <- reg_coef_table()
      } else if (s$model$method == "rf") {
        feat_df <- reg_rf_importance()
      }

      if (!is.null(feat_df) && nrow(feat_df) > 0) {
        write.csv(feat_df, feat_file, row.names = FALSE)
      } else {
        write.csv(data.frame(Message = "No feature importance available"),
          feat_file,
          row.names = FALSE
        )
      }
      files <- c(files, feat_file)

      # 4. Observed vs Predicted plot (PDF)
      obs_pred_file <- file.path(tmpdir, paste0(prefix, "_observed_vs_predicted_", subset_id, ".pdf"))
      if (!is.null(s$obs_pred_plot)) {
        ggsave(obs_pred_file,
          plot = s$obs_pred_plot,
          device = if (capabilities("cairo")) cairo_pdf else pdf,
          width = 6, height = 6, units = "in"
        )
      } else {
        pdf(obs_pred_file)
        plot.new()
        text(0.5, 0.5, "Plot not available")
        dev.off()
      }
      files <- c(files, obs_pred_file)

      # 5. Residual plot (PDF)
      resid_file <- file.path(tmpdir, paste0(prefix, "_residuals_vs_fitted_", subset_id, ".pdf"))
      if (!is.null(s$residual_plot)) {
        ggsave(resid_file,
          plot = s$residual_plot,
          device = if (capabilities("cairo")) cairo_pdf else pdf,
          width = 6, height = 6, units = "in"
        )
      } else {
        pdf(resid_file)
        plot.new()
        text(0.5, 0.5, "Plot not available")
        dev.off()
      }
      files <- c(files, resid_file)

      # Bundle into zip
      zip::zip(zipfile = file, files = files, mode = "cherry-pick")
    },
    contentType = "application/zip"
  )

  # ========== SCCOMP TAB LOGIC ==========

  # Helper functions for sccomp interval plots
  interval_plots <- function(x, subset_id = NULL) {
    # Derive contrasts_string from a 'factor' column if present and non-NA;
    # otherwise leave parameter names untouched. This prevents gsub() being
    # called with NA (which would produce all-NA parameters) when plotting
    # contrast results where 'factor' may be absent or all-NA.
    contrasts_string <- NA_character_
    if ("factor" %in% colnames(x)) {
      non_na <- x$factor[!is.na(x$factor)]
      if (length(non_na) > 0) contrasts_string <- as.character(non_na[1])
    }
    if (!is.na(contrasts_string) && nzchar(contrasts_string)) {
      x$parameter <- gsub(pattern = contrasts_string, replacement = "", x = as.character(x$parameter))
    } else {
      x$parameter <- as.character(x$parameter)
    }

    attr_x <- attributes(x)
    res_x <- as.data.frame(x)

    # Remove Intercept parameter (not meaningful for interpretation)
    res_x <- res_x[!grepl("Intercept", res_x$parameter, ignore.case = TRUE), ]

    if (nrow(res_x) == 0) {
      return(list())
    }

    param_split <- split(x = res_x, f = res_x$parameter)

    iterate_parameters <- function(x2) {
      # Use cell_group column (sccomp standard) or celltype if it exists
      celltype_col <- if ("cell_group" %in% colnames(x2)) "cell_group" else "celltype"

      factor_order <- x2[[celltype_col]][order(x2$c_effect, decreasing = FALSE)]
      x2[[celltype_col]] <- factor(x2[[celltype_col]], levels = factor_order)

      param <- x2$parameter[1]
      # Remove backticks from parameter names (common in contrasts)
      param_clean <- gsub("`", "", param)
      if (F) { # will come back to this later
        # For contrasts like "VariableLevel1 - VariableLevel2", remove prefix from both sides
        if (grepl(" - ", param_clean)) {
          # Split by " - ", remove prefix from each part, then rejoin
          parts <- strsplit(param_clean, " - ")[[1]]
          parts_clean <- sapply(parts, function(p) gsub("^[^_]+_", "", trimws(p)))
          param_clean <- paste(parts_clean, collapse = " - ")
        } else {
          # Remove prefix before first underscore (e.g., "Status" from "StatusDiabetes_Bad_Control")
          param_clean <- gsub("^[^_]+_", "", param_clean)
        }
      }
      param_str <- gsub(pattern = "\\) \\- \\(", replacement = ") -\n(", x = param_clean)
      param_str <- gsub(pattern = "cmv", replacement = "pp65", x = param_str) # hard-coded replacement
      param_str <- gsub(pattern = "HIVp", replacement = "HIV+", x = param_str)
      param_str <- gsub(pattern = "HIVn", replacement = "HIV-", x = param_str)

      # Build title and subtitle
      subtitle_text <- if (!is.null(subset_id)) paste0("Subset ID: ", subset_id) else NULL

      # Title wrapping: only allow breaks at ' - ' separators while respecting
      # the maximum width. If no separators present, fall back to standard wrap.
      title_wrapped <- param_str
      if (nchar(param_str) > 50) {
        if (grepl(" - ", param_str)) {
          parts <- stringr::str_split(param_str, " - ")[[1]]
          parts <- stringr::str_trim(parts)
          lines <- character(0)
          cur <- parts[1]
          if (length(parts) > 1) {
            for (i in seq(2, length(parts))) {
              cand <- paste(cur, parts[i], sep = " - ")
              if (nchar(cand) <= 50) {
                cur <- cand
              } else {
                lines <- c(lines, cur)
                cur <- parts[i]
              }
            }
          }
          lines <- c(lines, cur)
          title_wrapped <- paste(lines, collapse = " -\n")
        } else {
          title_wrapped <- stringr::str_wrap(param_str, width = 50)
        }
      }

      int_pl <- ggplot(data = x2, mapping = aes(x = c_effect, y = .data[[celltype_col]], color = c_FDR < 0.05)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
        geom_errorbar(aes(xmin = c_lower, xmax = c_upper, color = c_FDR < 0.05), linewidth = 0.8, width = 0.5) +
        geom_point(mapping = aes(fill = c_FDR < 0.05), size = 4, pch = 21, stroke = 0.4, color = "black") +
        scale_color_manual(values = c("grey30", "red")) +
        scale_fill_manual(values = c("grey30", "red")) +
        labs(
          x = "Credible interval\n(log-odds scale)",
          y = "Cell Group",
          title = title_wrapped,
          subtitle = subtitle_text
        ) +
        theme_bw(base_size = 14) +
        theme(
          axis.text.y = element_text(size = 16),
          axis.ticks.y = element_blank(), # Remove y-axis tick marks
          axis.text.x = element_text(size = 15),
          axis.title.x = element_text(size = 16),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18, hjust = 0.5),
          plot.subtitle = element_text(size = 14, hjust = 0.5),
          legend.position = "bottom",
          legend.text = element_text(size = 14),
          legend.title.position = "top",
          legend.title = element_text(size = 15, hjust = 0.5)
        )
      return(int_pl)
    }
    iterate_pl <- lapply(X = param_split, FUN = iterate_parameters)
    return(iterate_pl)
  }

  add_caption <- function(plot_obj) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      return(plot_obj) # Return plot without caption if patchwork not available
    }
    combined_plot <- plot_obj +
      patchwork::plot_annotation(
        caption = paste0(
          "Bayesian FDR: Stephens' method (doi: 10.1093/biostatistics/kxw041)\n",
          "FDR-significant populations may cross fold change thresholds because Bayesian FDR considers posterior probabilities rather than p-values.\n",
          "The method sorts null hypothesis probabilities in ascending order and calculates cumulative averages for robust false discovery control."
        ),
        theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0))
      )
    return(combined_plot)
  }

  # Update variable choices
  observeEvent(rv$meta_sample, {
    req(rv$meta_sample)
    meta_cols <- colnames(rv$meta_sample)
    categorical_choices <- sort(meta_cols[sapply(rv$meta_sample, function(x) is.character(x) || is.factor(x))])
    categorical_choices <- setdiff(categorical_choices, c("patient_ID", "run_date", "source"))
    categorical_choices <- filter_by_global_settings(categorical_choices)

    # For simple mode
    updatePickerInput(session, "sccomp_group_var", choices = c("", categorical_choices))
    updatePickerInput(session, "sccomp_formula_vars", choices = c("", categorical_choices))

    # For custom mode - show all variables including continuous
    all_vars <- sort(colnames(rv$meta_sample))
    all_vars <- setdiff(all_vars, c("run_date", "source"))
    all_vars <- filter_by_global_settings(all_vars)
    updatePickerInput(session, "sccomp_available_vars", choices = all_vars)
  })

  # Update reference level choices when grouping variable changes
  observeEvent(input$sccomp_group_var, {
    req(rv$meta_sample, input$sccomp_group_var)

    if (input$sccomp_group_var == "") {
      return()
    }

    # Get unique levels from the selected grouping variable
    if (input$sccomp_group_var %in% colnames(rv$meta_sample)) {
      unique_levels <- sort(unique(as.character(rv$meta_sample[[input$sccomp_group_var]])))
      unique_levels <- unique_levels[!is.na(unique_levels)]

      # Update picker with first level as default
      updatePickerInput(session, "sccomp_reference_level",
        choices = unique_levels,
        selected = unique_levels[1]
      )
    }
  })

  # Reactive to store sccomp results
  sccomp_state <- reactiveVal(NULL)

  # Run sccomp_estimate
  observeEvent(input$run_sccomp, {
    req(rv$meta_sample, rv$counts_sample)

    # Validate counts data exists
    if (is.null(rv$counts_sample)) {
      showNotification("No counts data available. sccomp requires obj$cluster$counts.", type = "error")
      return()
    }

    # Build formula based on mode
    formula_str <- if (input$sccomp_formula_mode == "simple") {
      req(input$sccomp_group_var)
      group_var <- input$sccomp_group_var
      formula_vars <- input$sccomp_formula_vars

      if (length(formula_vars) > 0) {
        if (input$sccomp_interactions) {
          # Use * for interactions
          paste0("~ ", paste(c(group_var, formula_vars), collapse = " * "))
        } else {
          # Use + for additive effects
          paste0("~ ", group_var, " + ", paste(formula_vars, collapse = " + "))
        }
      } else {
        paste0("~ ", group_var)
      }
    } else {
      # Custom mode
      req(input$sccomp_custom_formula)
      input$sccomp_custom_formula
    }

    tryCatch(
      {
        # Log formula information
        message("\n=== sccomp_estimate Configuration ===")
        message("Formula mode: ", input$sccomp_formula_mode)
        message("Formula string: ", formula_str)
        message("Parsed formula: ", deparse(as.formula(formula_str)))
        if (input$sccomp_formula_mode == "simple") {
          message("Primary grouping variable: ", input$sccomp_group_var)
          if (length(input$sccomp_formula_vars) > 0) {
            message("Additional covariates: ", paste(input$sccomp_formula_vars, collapse = ", "))
            message("Interactions enabled: ", input$sccomp_interactions)
          }
        }

        # Prepare data: convert counts matrix to long format with metadata
        message("\n=== Data Preparation ===")
        message("Counts matrix dimensions: ", nrow(rv$counts_sample), " samples × ", ncol(rv$counts_sample), " clusters")
        message("Counts matrix rownames (first 3): ", paste(head(rownames(rv$counts_sample), 3), collapse = ", "))

        counts_long <- as.data.frame(rv$counts_sample) %>%
          tibble::rownames_to_column("sample") %>%
          tidyr::pivot_longer(cols = -sample, names_to = "cell_group", values_to = "count")

        message("Counts long format: ", nrow(counts_long), " rows")
        message("Sample values in counts (first 3): ", paste(head(unique(counts_long$sample), 3), collapse = ", "))

        # Prepare metadata - use patient_ID column if it exists, otherwise use rownames
        meta_df <- rv$meta_sample
        if ("patient_ID" %in% colnames(meta_df)) {
          # patient_ID column exists - use it as the merge key
          message("Using patient_ID column from metadata for merging")
          # Don't add rownames as a column, patient_ID is already there
        } else {
          # No patient_ID column - use rownames
          message("Using metadata rownames for merging")
          meta_df <- meta_df %>% tibble::rownames_to_column("patient_ID")
        }

        message("Metadata dimensions: ", nrow(meta_df), " samples × ", ncol(meta_df), " variables")
        message("Metadata patient_ID values (first 3): ", paste(head(meta_df$patient_ID, 3), collapse = ", "))
        message("Metadata column names: ", paste(colnames(meta_df), collapse = ", "))

        # Check formula variables in metadata BEFORE merge
        formula_vars_check <- all.vars(as.formula(formula_str))
        message("\nChecking formula variables in metadata BEFORE merge:")
        for (var in formula_vars_check) {
          if (var %in% colnames(meta_df)) {
            na_count <- sum(is.na(meta_df[[var]]))
            unique_vals <- length(unique(meta_df[[var]][!is.na(meta_df[[var]])]))
            message(
              "  ", var, ": ", na_count, " NAs out of ", nrow(meta_df), " rows, ",
              unique_vals, " unique non-NA values"
            )
            if (unique_vals > 0 && unique_vals <= 10) {
              message("    Unique values: ", paste(head(unique(meta_df[[var]][!is.na(meta_df[[var]])]), 10), collapse = ", "))
            }
          } else {
            message("  ", var, ": NOT FOUND in metadata")
          }
        }

        sccomp_data <- counts_long %>%
          dplyr::left_join(meta_df, by = c("sample" = "patient_ID"))

        message("Merged data dimensions: ", nrow(sccomp_data), " rows × ", ncol(sccomp_data), " columns")

        # Remove rows with NA in formula variables
        formula_vars_used <- all.vars(as.formula(formula_str))
        message("Formula variables to check: ", paste(formula_vars_used, collapse = ", "))

        rows_before <- nrow(sccomp_data)
        for (var in formula_vars_used) {
          if (var %in% colnames(sccomp_data)) {
            na_count <- sum(is.na(sccomp_data[[var]]))
            if (na_count > 0) {
              message("  Removing ", na_count, " rows with NA in variable: ", var)
            }
            sccomp_data <- sccomp_data %>% dplyr::filter(!is.na(.data[[var]]))
          } else {
            message("  WARNING: Variable '", var, "' not found in data columns")
          }
        }
        message("Rows after NA filtering: ", nrow(sccomp_data), " (removed ", rows_before - nrow(sccomp_data), ")")

        # Check if all rows were filtered out
        if (nrow(sccomp_data) == 0) {
          stop(
            "All rows were filtered out due to NA values in formula variables. ",
            "Check that the selected variables have non-NA values in your metadata. ",
            "Variables used in formula: ", paste(formula_vars_used, collapse = ", ")
          )
        }

        # Clean special characters in formula variables to prevent sccomp errors
        # Only sanitize character/factor variables (do NOT coerce numeric/integer columns)
        # This must happen BEFORE releveling so the reference level matches the cleaned values
        message("\n=== Cleaning Special Characters ===")
        formula_vars_to_clean <- all.vars(as.formula(formula_str))
        for (var in formula_vars_to_clean) {
          if (var %in% colnames(sccomp_data)) {
            original_vals <- unique(sccomp_data[[var]])

            # Only operate on factors or character vectors; leave numeric/integer untouched
            if (is.factor(sccomp_data[[var]])) {
              # Clean factor levels in-place
              levels(sccomp_data[[var]]) <- gsub("\\+", "p", levels(sccomp_data[[var]]))
              levels(sccomp_data[[var]]) <- gsub("\\-", "n", levels(sccomp_data[[var]]))
            } else if (is.character(sccomp_data[[var]])) {
              sccomp_data[[var]] <- gsub("\\+", "p", sccomp_data[[var]])
              sccomp_data[[var]] <- gsub("\\-", "n", sccomp_data[[var]])
            } else {
              # Numeric or other types: skip cleaning to preserve type
              next
            }

            cleaned_vals <- unique(sccomp_data[[var]])
            if (!identical(sort(as.character(original_vals)), sort(as.character(cleaned_vals)))) {
              message("  Cleaned variable '", var, "':")
              message("    Original values: ", paste(head(original_vals, 10), collapse = ", "))
              message("    Cleaned values: ", paste(head(cleaned_vals, 10), collapse = ", "))
            }
          }
        }

        # Relevel variables based on mode
        if (input$sccomp_formula_mode == "simple" && !is.null(input$sccomp_group_var) &&
          input$sccomp_group_var != "" && !is.null(input$sccomp_reference_level)) {
          # Simple mode: relevel the grouping variable
          group_var <- input$sccomp_group_var
          ref_level <- input$sccomp_reference_level

          # Clean the reference level to match the cleaned data
          ref_level_cleaned <- gsub("\\+", "p", ref_level)
          ref_level_cleaned <- gsub("\\-", "n", ref_level_cleaned)

          if (group_var %in% colnames(sccomp_data)) {
            message("\n=== Releveling Grouping Variable ===")
            message("Variable: ", group_var)
            message("Reference level (original): ", ref_level)
            message("Reference level (cleaned): ", ref_level_cleaned)

            sccomp_data[[group_var]] <- factor(sccomp_data[[group_var]])
            sccomp_data[[group_var]] <- relevel(sccomp_data[[group_var]], ref = ref_level_cleaned)

            message("Factor levels after releveling: ", paste(levels(sccomp_data[[group_var]]), collapse = ", "))
          }
        } else if (input$sccomp_formula_mode == "custom" &&
          !is.null(input$sccomp_custom_reference_levels) &&
          input$sccomp_custom_reference_levels != "") {
          # Custom mode: parse and apply reference levels
          message("\n=== Releveling Variables (Custom Mode) ===")
          ref_specs <- strsplit(input$sccomp_custom_reference_levels, ";")[[1]]
          ref_specs <- trimws(ref_specs)

          for (spec in ref_specs) {
            if (grepl("=", spec)) {
              parts <- strsplit(spec, "=")[[1]]
              if (length(parts) == 2) {
                var_name <- trimws(parts[1])
                ref_level <- trimws(parts[2])

                # Clean the reference level to match the cleaned data
                ref_level_cleaned <- gsub("\\+", "p", ref_level)
                ref_level_cleaned <- gsub("\\-", "n", ref_level_cleaned)

                if (var_name %in% colnames(sccomp_data)) {
                  message("Variable: ", var_name)
                  message("Reference level (original): ", ref_level)
                  message("Reference level (cleaned): ", ref_level_cleaned)
                  
                  if(mean(class(sccomp_data[[var_name]]) %in% c('factor','character'))!=0) {
                    sccomp_data[[var_name]] <- factor(sccomp_data[[var_name]])
                    sccomp_data[[var_name]] <- relevel(sccomp_data[[var_name]], ref = ref_level_cleaned)
                    message("Factor levels after releveling: ", paste(levels(sccomp_data[[var_name]]), collapse = ", "))
                  } else {
                    message("Term is not categorical; skipping releveling: ", var_name)
                  }

                } else {
                  message("Warning: Variable '", var_name, "' not found in data")
                }
              }
            }
          }
        }

        # Log data structure
        message("\n=== sccomp_data Structure ===")
        message("Column names: ", paste(colnames(sccomp_data), collapse = ", "))
        message("\nColumn types:")
        for (col in colnames(sccomp_data)) {
          message("  ", col, ": ", class(sccomp_data[[col]])[1])
        }
        message("\nFirst 6 rows of sccomp_data:")
        print(head(sccomp_data))

        # Validate required columns
        if (!"sample" %in% colnames(sccomp_data)) {
          stop("'sample' column missing from sccomp_data")
        }
        if (!"cell_group" %in% colnames(sccomp_data)) {
          stop("'cell_group' column missing from sccomp_data")
        }
        if (!"count" %in% colnames(sccomp_data)) {
          stop("'count' column missing from sccomp_data")
        }

        message("\nUnique samples: ", length(unique(sccomp_data$sample)))
        message("Unique cell groups: ", length(unique(sccomp_data$cell_group)))
        if (nrow(sccomp_data) > 0) {
          message("Count range: ", min(sccomp_data$count), " to ", max(sccomp_data$count))
        }

        # Run sccomp
        if (!requireNamespace("sccomp", quietly = TRUE)) {
          showNotification("sccomp package not installed. Install with: devtools::install_github('MangiolaLaboratory/sccomp')",
            type = "error", duration = NULL
          )
          return()
        }

        message("\n=== Running sccomp_estimate ===")
        message("Using ", input$sccomp_cores, " cores")

        # Run sccomp inside a per-call temporary directory and remove intermediate files
        tmpdir <- tempfile("sccomp_draws_")
        dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
        oldwd <- getwd()
        setwd(tmpdir)
        tmp_result <- NULL
        tryCatch({
          tmp_result <- sccomp::sccomp_estimate(
            .data = sccomp_data,
            formula_composition = as.formula(formula_str),
            sample = "sample",
            cell_group = "cell_group",
            abundance = "count",
            bimodal_mean_variability_association = FALSE,
            cores = input$sccomp_cores
          )
          message("sccomp_estimate completed successfully")

          # Automatically run sccomp_test to get FDR values (required for significance testing)
          message("\n=== Running sccomp_test (automatic) ===")
          tmp_result <- sccomp::sccomp_test(tmp_result)
          message("sccomp_test completed successfully")
        }, error = function(e) {
          stop(e)
        }, finally = {
          # Restore working directory and clean up temporary files
          try(setwd(oldwd), silent = TRUE)
          if (dir.exists(tmpdir)) {
            try(unlink(tmpdir, recursive = TRUE, force = TRUE), silent = TRUE)
          }
        })

        # Use the result from the temporary run
        result <- tmp_result

        # Store results (capture subset_id used for this run so subsequent UI
        # changes to `rv$subset_id` do not retroactively change these plots)
        sccomp_state(list(
          estimate_result = result,
          test_result = NULL,
          formula = formula_str,
          n_samples = length(unique(sccomp_data$sample)),
          n_clusters = length(unique(sccomp_data$cell_group)),
          subset_id = rv$subset_id %||% "data"
        ))

        output$sccomp_cleared_msg <- renderText(NULL)
        showNotification("sccomp analysis completed — intermediate files removed.", type = "message", duration = 5)
      },
      error = function(e) {
        message("\n=== ERROR in sccomp_estimate ===")
        message("Error message: ", e$message)
        message("Error call: ", deparse(e$call))
        message("\nFull traceback:")
        print(traceback())
        showNotification(paste("Error running sccomp_estimate:", e$message), type = "error", duration = NULL)
      }
    )
  })

  # Run sccomp_test for contrasts
  observeEvent(input$run_sccomp_test, {
    s <- sccomp_state()
    req(s, s$estimate_result, input$sccomp_contrast)

    contrast_str <- input$sccomp_contrast

    tryCatch(
      {
        if (!requireNamespace("sccomp", quietly = TRUE)) {
          showNotification("sccomp package not installed.", type = "error")
          return()
        }

        # Run sccomp_test with contrast
        test_result <- sccomp::sccomp_test(
          s$estimate_result,
          contrasts = contrast_str
        )

        # Update state with test results
        s$test_result <- test_result
        sccomp_state(s)

        showNotification("sccomp_test completed!", type = "message", duration = 3)
      },
      error = function(e) {
        showNotification(paste("Error running sccomp_test:", e$message), type = "error", duration = NULL)
      }
    )
  })

  # Clear results
  observeEvent(input$reset_sccomp, {
    sccomp_state(NULL)
    output$sccomp_cleared_msg <- renderText("Results cleared. Configure settings and click 'Run sccomp_estimate'.")
  })

  # Display available parameters for contrasts
  output$sccomp_available_params <- renderPrint({
    s <- sccomp_state()
    req(s, s$estimate_result)

    # Get unique parameters (excluding NA and Intercept)
    params <- unique(s$estimate_result$parameter)
    params <- params[!is.na(params)]
    params <- params[!grepl("Intercept", params, ignore.case = TRUE)]

    if (length(params) > 0) {
      # Show parameters with backticks for clarity
      params_formatted <- paste0("`", params, "`")
      cat("Use these in contrasts:\n")
      cat(paste(params_formatted, collapse = "\n"))
      cat("\n\nExample contrast:\n")
      if (length(params) >= 2) {
        cat("`", params[1], "` - `", params[2], "`", sep = "")
      } else {
        cat("(Need at least 2 parameters for contrasts)")
      }
    } else {
      cat("No parameters available for contrasts.\nNote: Standard model parameterization doesn't support post-hoc contrasts.\nUse '~ 0 + variable' formula instead.")
    }
  })

  # Check if results exist
  output$hasSccompResults <- reactive({
    !is.null(sccomp_state())
  })
  outputOptions(output, "hasSccompResults", suspendWhenHidden = FALSE)

  output$hasSccompTestResults <- reactive({
    s <- sccomp_state()
    !is.null(s) && !is.null(s$test_result)
  })
  outputOptions(output, "hasSccompTestResults", suspendWhenHidden = FALSE)

  # Warning for no-intercept formulas without contrasts
  output$sccomp_intercept_warning <- renderUI({
    s <- sccomp_state()
    req(s)

    # Check if formula has no intercept (contains ~ 0 + or ~ -1 +)
    has_no_intercept <- grepl("~\\s*(0|\\-1)\\s*\\+", s$formula)

    # Show warning whenever no intercept formula is used (regardless of contrast test)
    if (has_no_intercept) {
      tags$div(
        style = "background-color: #fff3cd; border: 1px solid #ffc107; border-radius: 4px; padding: 15px; margin-bottom: 15px;",
        tags$strong(style = "color: #856404;", "⚠ No-Intercept Formula Detected"),
        tags$p(
          style = "color: #856404; margin-top: 10px; margin-bottom: 0;",
          "When using a no-intercept formula (e.g., ", tags$code("~ 0 + condition"), "),",
          " the parameters represent absolute compositions for each group, not differences between groups.",
          " Nearly all clusters will show 'significant' intervals because they test if composition differs from zero,",
          " which is rarely meaningful."
        ),
        tags$p(
          style = "color: #856404; margin-top: 10px; margin-bottom: 0;",
          tags$strong("Recommendation:"), " Specify a custom contrast below (e.g., ",
          tags$code("`conditionTreated` - `conditionControl`"),
          ") to test meaningful differences between groups."
        )
      )
    } else {
      NULL
    }
  })

  # Summary output
  output$sccomp_summary <- renderPrint({
    s <- sccomp_state()
    req(s)

    cat("Formula:", s$formula, "\n")
    cat("Samples analyzed:", s$n_samples, "\n")
    cat("Clusters tested:", s$n_clusters, "\n\n")

    # Extract test results from estimate
    if (!is.null(s$estimate_result)) {
      test_res <- s$estimate_result %>%
        dplyr::filter(!is.na(c_effect))

      # Check if c_FDR column exists (from sccomp_test)
      if ("c_FDR" %in% colnames(test_res)) {
        sig_clusters <- sum(test_res$c_FDR < 0.05, na.rm = TRUE)
        cat("Significant clusters (FDR < 0.05):", sig_clusters, "/", nrow(test_res), "\n")
      } else {
        cat("Note: Run sccomp_test to calculate FDR values\n")
        cat("Total effects estimated:", nrow(test_res), "\n")
      }
    }

    # If contrast test was run
    if (!is.null(s$test_result)) {
      cat("\nContrast test performed\n")
      if ("c_FDR" %in% colnames(s$test_result)) {
        test_sig <- sum(s$test_result$c_FDR < 0.05, na.rm = TRUE)
        cat("Significant in contrast (FDR < 0.05):", test_sig, "\n")
      }
    }
  })

  # Conditional UI for significant clusters table
  output$sccomp_table_ui <- renderUI({
    s <- sccomp_state()
    req(s)

    # Check if formula has no intercept
    has_no_intercept <- grepl("~\\s*(0|\\-1)\\s*\\+", s$formula)

    # If no intercept, always show message instead of table (even after contrast test)
    # The estimate_result table is not meaningful with no-intercept formulas
    # Users should rely on the Contrast Test Results table instead
    if (has_no_intercept) {
      tagList(
        h4("Significant Clusters"),
        tags$div(
          style = "background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 4px; padding: 15px; margin-bottom: 15px;",
          tags$p(
            style = "color: #6c757d; margin-bottom: 0;",
            tags$em(
              "Table hidden: With no-intercept formulas, this table tests if composition differs from zero, ",
              "which is not meaningful. Please specify a custom contrast above to test group differences. ",
              "Significant results from your contrast will appear in the 'Contrast Test Results' section on the right."
            )
          )
        )
      )
    } else {
      # Check if there are any significant results
      has_results <- FALSE
      if ("c_FDR" %in% colnames(s$estimate_result)) {
        sig_df <- s$estimate_result %>%
          dplyr::filter(!is.na(c_effect), c_FDR < 0.05, !grepl("Intercept", parameter, ignore.case = TRUE))
        has_results <- nrow(sig_df) > 0
      }

      if (has_results) {
        tagList(
          h4("Significant Clusters"),
          tableOutput("sccomp_table")
        )
      } else {
        tagList(
          h4("Significant Clusters"),
          tags$p(style = "color: #6c757d; font-style: italic;", "No significant clusters found at FDR < 0.05")
        )
      }
    }
  })

  # Results table
  output$sccomp_table <- renderTable(
    {
      s <- sccomp_state()
      req(s, s$estimate_result)

      # Check if c_FDR column exists
      if (!"c_FDR" %in% colnames(s$estimate_result)) {
        return(NULL)
      }

      # Extract significant results from estimate, excluding Intercept
      df <- s$estimate_result %>%
        dplyr::filter(!is.na(c_effect), c_FDR < 0.05, !grepl("Intercept", parameter, ignore.case = TRUE)) %>%
        dplyr::select(cell_group, parameter, c_effect, c_lower, c_upper, c_FDR) %>%
        dplyr::arrange(c_FDR)

      if (nrow(df) == 0) {
        return(NULL)
      }

      df
    },
    digits = 6,
    sanitize.text.function = function(x) x
  )

  # Conditional UI for contrast test results table
  output$sccomp_test_table_ui <- renderUI({
    s <- sccomp_state()
    req(s, s$test_result)

    # Check if there are any significant results
    has_results <- FALSE
    if ("c_FDR" %in% colnames(s$test_result)) {
      sig_df <- s$test_result %>%
        dplyr::filter(!is.na(c_effect), c_FDR < 0.05)
      has_results <- nrow(sig_df) > 0
    }

    if (has_results) {
      tableOutput("sccomp_test_table")
    } else {
      tags$p(style = "color: #6c757d; font-style: italic;", "No significant clusters found at FDR < 0.05 for this contrast")
    }
  })

  # Contrast test results table
  output$sccomp_test_table <- renderTable(
    {
      s <- sccomp_state()
      req(s, s$test_result)

      # Check if c_FDR column exists
      if (!"c_FDR" %in% colnames(s$test_result)) {
        return(NULL)
      }

      df <- s$test_result %>%
        dplyr::filter(!is.na(c_effect), c_FDR < 0.05) %>%
        dplyr::select(cell_group, parameter, c_effect, c_lower, c_upper, c_FDR) %>%
        dplyr::arrange(c_FDR)

      if (nrow(df) == 0) {
        return(NULL)
      }

      df
    },
    digits = 6,
    sanitize.text.function = function(x) x
  )

  # Generate interval plot
  sccomp_interval_plot <- reactive({
    s <- sccomp_state()
    req(s, s$estimate_result)

    # Check if c_FDR column exists
    if (!"c_FDR" %in% colnames(s$estimate_result)) {
      return(NULL)
    }

    # Filter to rows with effect estimates
    plot_data <- s$estimate_result %>%
      dplyr::filter(!is.na(c_effect))

    if (nrow(plot_data) == 0) {
      return(NULL)
    }

    # Generate interval plots using the subset_id captured when sccomp was run
    s <- sccomp_state()
    subset_id <- s$subset_id %||% "data"
    plots_list <- interval_plots(x = plot_data, subset_id = subset_id)

    if (length(plots_list) == 0) {
      return(NULL)
    }

    # If multiple parameters, combine them in 4 columns
    if (length(plots_list) > 1) {
      # Arrange plots in 4 columns
      combined <- patchwork::wrap_plots(plots_list, ncol = 4)
      plot_with_caption <- add_caption(combined)
    } else {
      plot_with_caption <- add_caption(plots_list[[1]])
    }

    return(plot_with_caption)
  })

  # Calculate dynamic height for interval plot
  sccomp_plot_height <- reactive({
    s <- sccomp_state()
    req(s, s$estimate_result)

    if (!"c_FDR" %in% colnames(s$estimate_result)) {
      return(400)
    }

    plot_data <- s$estimate_result %>%
      dplyr::filter(!is.na(c_effect))

    if (nrow(plot_data) == 0) {
      return(400)
    }

    # Get number of unique parameters (excluding Intercept) and cell groups
    n_params <- length(unique(plot_data$parameter[!grepl("Intercept", plot_data$parameter, ignore.case = TRUE)]))

    if (n_params == 0) {
      return(400)
    }

    n_clusters <- s$n_clusters

    # Calculate more reasonable height based on clusters and rows
    # With 4 columns, height is for rows (ceiling of n_params/4)
    n_rows <- ceiling(n_params / 4)

    # Adjusted formula: base + (clusters * pixels_per_cluster)
    # Reduced from 30px to 18px per cluster for more reasonable sizing
    base_height <- 150
    pixels_per_cluster <- 18
    height_per_row <- base_height + (n_clusters * pixels_per_cluster)
    total_height <- max(300, n_rows * height_per_row)

    return(total_height)
  })

  # Calculate dynamic width for interval plot (for PDF export)
  sccomp_plot_width <- reactive({
    s <- sccomp_state()
    req(s, s$estimate_result)

    if (!"c_FDR" %in% colnames(s$estimate_result)) {
      return(16)
    }

    plot_data <- s$estimate_result %>%
      dplyr::filter(!is.na(c_effect))

    if (nrow(plot_data) == 0) {
      return(16)
    }

    # Get number of unique parameters (excluding Intercept)
    n_params <- length(unique(plot_data$parameter[!grepl("Intercept", plot_data$parameter, ignore.case = TRUE)]))

    if (n_params == 0) {
      return(16)
    }

    # Width based on number of plots
    # Single plot: 4 inches, 2-4 plots: 16 inches for 4-column layout
    if (n_params == 1) {
      return(4) # Single plot width
    } else if (n_params <= 4) {
      return(16) # Full width for up to 4 plots (4 columns)
    } else {
      return(16) # Full width for multiple rows
    }
  })

  output$sccomp_interval_plot <- renderPlot(
    {
      plot <- sccomp_interval_plot()
      req(plot)
      plot
    },
    height = function() sccomp_plot_height(),
    width = function() {
      # Adjust width for plot layout
      s <- sccomp_state()
      if (is.null(s) || is.null(s$estimate_result)) {
        return(1200)
      }

      plot_data <- s$estimate_result %>%
        dplyr::filter(!is.na(c_effect))
      n_params <- length(unique(plot_data$parameter[!grepl("Intercept", plot_data$parameter, ignore.case = TRUE)]))

      if (n_params == 1) {
        return(300) # Single plot width
      } else {
        return(900) # Increased width for 6/6 split
      }
    }
  )

  # Download interval plot
  output$download_sccomp_plot <- downloadHandler(
    filename = function() {
      subset_id <- rv$subset_id %||% "000000000"
      s <- sccomp_state()
      req(s)

      formula_clean <- gsub("[~\\s\\+\\*\\(\\)]", "_", s$formula)
      formula_clean <- gsub("_{2,}", "_", formula_clean)
      formula_clean <- gsub("^_|_$", "", formula_clean)

      paste0("sccomp_intervals_", formula_clean, "_", subset_id, ".pdf")
    },
    content = function(file) {
      plot <- sccomp_interval_plot()
      req(plot)

      s <- sccomp_state()
      req(s, s$estimate_result)

      # Calculate height based on clusters and rows
      plot_data <- s$estimate_result %>%
        dplyr::filter(!is.na(c_effect))

      n_params <- length(unique(plot_data$parameter[!grepl("Intercept", plot_data$parameter, ignore.case = TRUE)]))
      n_clusters <- s$n_clusters
      n_rows <- ceiling(n_params / 4) # 4 columns

      # Height calculation: adjusted for PDF (slightly different from screen)
      base_height <- 4 # inches
      height_per_cluster <- 0.17 # inches per cluster (slightly increased for better spacing)
      height_per_row <- base_height + (n_clusters * height_per_cluster)
      height_val <- max(5, n_rows * height_per_row)

      # Width from reactive function
      width_val <- sccomp_plot_width()

      ggsave(file,
        plot = plot, device = if (capabilities("cairo")) cairo_pdf else pdf,
        width = width_val, height = height_val, units = "in", limitsize = FALSE
      )
    }
  )

  # Export results
  output$export_sccomp_results <- downloadHandler(
    filename = function() {
      subset_id <- rv$subset_id %||% "000000000"
      s <- sccomp_state()
      req(s)

      # Build descriptive name from formula
      formula_clean <- gsub("[~\\s\\+\\*\\(\\)]", "_", s$formula)
      formula_clean <- gsub("_{2,}", "_", formula_clean)
      formula_clean <- gsub("^_|_$", "", formula_clean)
      formula_clean <- tolower(formula_clean)

      paste0("sccomp_", formula_clean, "_", subset_id, ".csv")
    },
    content = function(file) {
      s <- sccomp_state()
      req(s)

      # Export both estimate and test results if available
      if (!is.null(s$test_result)) {
        df <- as.data.frame(s$test_result) %>%
          dplyr::filter(!is.na(c_effect))

        # Select available columns
        available_cols <- c("cell_group", "parameter", "c_effect", "c_lower", "c_upper", "c_p_value", "c_FDR")
        available_cols <- available_cols[available_cols %in% colnames(df)]

        df <- df %>%
          dplyr::select(dplyr::all_of(available_cols)) %>%
          dplyr::arrange(c_FDR)
      } else if (!is.null(s$estimate_result)) {
        df <- as.data.frame(s$estimate_result) %>%
          dplyr::filter(!is.na(c_effect))

        # Select available columns
        available_cols <- c("cell_group", "parameter", "c_effect", "c_lower", "c_upper", "c_p_value", "c_FDR")
        available_cols <- available_cols[available_cols %in% colnames(df)]

        df <- df %>%
          dplyr::select(dplyr::all_of(available_cols)) %>%
          dplyr::arrange(c_FDR)
      } else {
        df <- data.frame(Message = "No results available")
      }

      write.csv(df, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )

  # Generate contrast interval plot
  sccomp_contrast_plot <- reactive({
    s <- sccomp_state()
    req(s, s$test_result)

    # Check if c_FDR column exists
    if (!"c_FDR" %in% colnames(s$test_result)) {
      return(NULL)
    }

    # Filter to rows with effect estimates
    plot_data <- s$test_result %>%
      dplyr::filter(!is.na(c_effect))

    if (nrow(plot_data) == 0) {
      return(NULL)
    }

    # Generate interval plots using the subset_id captured when sccomp was run
    s <- sccomp_state()
    subset_id <- s$subset_id %||% "data"
    plots_list <- interval_plots(x = plot_data, subset_id = subset_id)

    if (length(plots_list) == 0) {
      return(NULL)
    }

    # If multiple parameters, combine them in 4 columns
    if (length(plots_list) > 1) {
      combined <- patchwork::wrap_plots(plots_list, ncol = 4)
      plot_with_caption <- add_caption(combined)
    } else {
      plot_with_caption <- add_caption(plots_list[[1]])
    }

    return(plot_with_caption)
  })

  # Calculate dynamic height for contrast plot - match main plot height
  sccomp_contrast_plot_height <- reactive({
    # Use the same height as the main plot
    return(sccomp_plot_height())
  })

  # Calculate dynamic width for contrast plot
  sccomp_contrast_plot_width <- reactive({
    s <- sccomp_state()
    req(s, s$test_result)

    if (!"c_FDR" %in% colnames(s$test_result)) {
      return(16)
    }

    plot_data <- s$test_result %>%
      dplyr::filter(!is.na(c_effect))

    if (nrow(plot_data) == 0) {
      return(16)
    }

    n_params <- length(unique(plot_data$parameter))

    if (n_params == 0) {
      return(32) # Doubled from 16
    }

    if (n_params == 1) {
      return(8) # Doubled from 4
    } else if (n_params <= 4) {
      return(32) # Doubled from 16
    } else {
      return(32) # Doubled from 16
    }
  })

  output$sccomp_contrast_plot <- renderPlot(
    {
      plot <- sccomp_contrast_plot()
      req(plot)
      plot
    },
    height = function() sccomp_contrast_plot_height(),
    width = function() {
      s <- sccomp_state()
      if (is.null(s) || is.null(s$test_result)) {
        return(2400)
      }

      plot_data <- s$test_result %>%
        dplyr::filter(!is.na(c_effect))
      n_params <- length(unique(plot_data$parameter))

      if (n_params == 1) {
        return(500) # Single plot width
      } else {
        return(1600) # Further increased width for better visibility
      }
    }
  )

  # Download contrast interval plot
  output$download_sccomp_contrast_plot <- downloadHandler(
    filename = function() {
      subset_id <- rv$subset_id %||% "000000000"
      s <- sccomp_state()
      req(s)

      # Get contrast string from test result parameters
      contrast_name <- "contrast"
      if (!is.null(s$test_result) && "parameter" %in% colnames(s$test_result)) {
        params <- unique(s$test_result$parameter)
        if (length(params) > 0) {
          contrast_name <- gsub("[^A-Za-z0-9_]", "_", params[1])
        }
      }

      paste0("sccomp_contrast_intervals_", contrast_name, "_", subset_id, ".pdf")
    },
    content = function(file) {
      plot <- sccomp_contrast_plot()
      req(plot)

      s <- sccomp_state()
      req(s, s$test_result)

      # Calculate height to match main plot - use same formula
      plot_data <- s$estimate_result %>%
        dplyr::filter(!is.na(c_effect))

      n_params <- length(unique(plot_data$parameter[!grepl("Intercept", plot_data$parameter, ignore.case = TRUE)]))
      n_clusters <- s$n_clusters
      n_rows <- ceiling(n_params / 4) # 4 columns to match main plot

      # Height calculation - same as main plot
      base_height <- 4
      height_per_cluster <- 0.17
      height_per_row <- base_height + (n_clusters * height_per_cluster)
      height_val <- max(5, n_rows * height_per_row)

      # Width from reactive function
      width_val <- sccomp_contrast_plot_width()

      ggsave(file,
        plot = plot, device = if (capabilities("cairo")) cairo_pdf else pdf,
        width = width_val, height = height_val, units = "in", limitsize = FALSE
      )
    }
  )

  # ========== TIME TO EVENT ANALYSIS TAB LOGIC ==========

  # Initialize Time to Event Analysis dropdowns
  observeEvent(list(rv$meta_sample, rv$abundance_sample),
    {
      req(rv$meta_sample, rv$abundance_sample)

      meta_cols <- colnames(rv$meta_sample)
      continuous_choices <- sort(meta_cols[sapply(rv$meta_sample, is.numeric)])
      # Exclude patient_ID, source, run_date, and cluster from continuous outcome choices
      continuous_choices <- setdiff(continuous_choices, c("patient_ID", "source", "run_date", "cluster"))
      predictor_choices <- sort(meta_cols)
      # Exclude patient_ID, source, and run_date from predictor choices
      predictor_choices <- setdiff(predictor_choices, c("patient_ID", "source", "run_date"))
      if (!("cluster" %in% predictor_choices)) {
        predictor_choices <- c(predictor_choices, "cluster")
      }

      updatePickerInput(session, "surv_outcome", choices = continuous_choices, selected = NULL)
      updatePickerInput(session, "surv_predictors", choices = predictor_choices, selected = NULL)
      updatePickerInput(session, "surv_cluster_subset",
        choices = colnames(rv$abundance_sample),
        selected = character(0)
      )
    },
    ignoreInit = TRUE
  )

  # No threshold slider needed - split method is selected via dropdown

  # Main time to event analysis function
  run_surv <- function() {
    req(rv$meta_sample, rv$abundance_sample, input$surv_outcome, input$surv_predictors)
    
    # Get analysis mode
    analysis_mode <- input$surv_analysis_mode %||% "multivariate"
    
    if (analysis_mode == "univariate") {
      return(run_surv_univariate())
    } else {
      return(run_surv_multivariate())
    }
  }
  
  # Multivariate analysis (original behavior)
  run_surv_multivariate <- function() {
    req(rv$meta_sample, rv$abundance_sample, input$surv_outcome, input$surv_predictors)

    meta_patient <- rv$meta_sample

    # Validate outcome variable
    if (!(input$surv_outcome %in% colnames(meta_patient))) {
      stop(paste0("Outcome '", input$surv_outcome, "' not found in metadata."))
    }

    outcome_raw <- meta_patient[[input$surv_outcome]]
    outcome <- suppressWarnings(as.numeric(outcome_raw))

    if (all(is.na(outcome))) {
      stop(paste0("Outcome '", input$surv_outcome, "' cannot be converted to numeric."))
    }

    # Get split method (will calculate threshold from risk scores after model fitting)
    split_method <- input$surv_split_method %||% "median"

    # Build predictor matrix
    predictor_list <- input$surv_predictors
    design_df <- data.frame(patient_ID = meta_patient$patient_ID)

    # Add cluster predictors if selected
    if ("cluster" %in% predictor_list) {
      abund <- rv$abundance_sample

      # Filter to selected clusters
      if (!is.null(input$surv_cluster_subset) && length(input$surv_cluster_subset) > 0) {
        selected_clusters <- input$surv_cluster_subset
        keep_cols <- colnames(abund)[colnames(abund) %in% selected_clusters]
        if (length(keep_cols) == 0) {
          stop("No valid clusters selected.")
        }
        abund <- abund[, keep_cols, drop = FALSE]
      }

      # Merge abundance with design
      abund_df <- as.data.frame(abund)
      abund_df$patient_ID <- rownames(abund_df)
      
      # Sanitize cluster column names to avoid formula issues
      # Use custom replacements to preserve meaning: spaces->_, +->p, -->n
      # Also create a mapping to restore original names in output tables
      cluster_cols <- setdiff(names(abund_df), "patient_ID")
      name_mapping <- setNames(cluster_cols, cluster_cols)  # Initialize with identity mapping
      
      for (col in cluster_cols) {
        new_name <- col
        new_name <- gsub(" ", "_", new_name)
        new_name <- gsub("\\+", "p", new_name)
        new_name <- gsub("-", "n", new_name)
        # Remove any other special characters that could break formulas
        new_name <- gsub("[^A-Za-z0-9_]", "_", new_name)
        # Ensure it starts with a letter or underscore
        if (!grepl("^[A-Za-z_]", new_name)) {
          new_name <- paste0("X_", new_name)
        }
        if (new_name != col) {
          names(abund_df)[names(abund_df) == col] <- new_name
          name_mapping[new_name] <- col  # Store sanitized -> original mapping
        }
      }
      
      design_df <- merge(design_df, abund_df, by = "patient_ID", all.x = TRUE)
    }

    # Add metadata predictors
    meta_predictors <- setdiff(predictor_list, "cluster")
    if (length(meta_predictors) > 0) {
      meta_subset <- meta_patient[, c("patient_ID", meta_predictors), drop = FALSE]
      design_df <- merge(design_df, meta_subset, by = "patient_ID", all.x = TRUE)
    }

    # Attach outcome
    design_df$time <- outcome

    # Remove rows with missing values
    n_before <- nrow(design_df)
    design_df_complete <- design_df[complete.cases(design_df), ]
    n_after <- nrow(design_df_complete)
    n_dropped <- n_before - n_after
    dropped_ids <- setdiff(design_df$patient_ID, design_df_complete$patient_ID)

    if (n_after < 10) {
      stop("Insufficient complete cases after removing missing values (n < 10).")
    }

    # For survival analysis, we need an event indicator
    # Since we don't have event data, we'll assume all events occurred (status = 1)
    # This is a simplification - in real survival analysis you'd have censoring info
    design_df_complete$status <- 1

    # Check if outcome variable is being used as a predictor
    predictor_cols <- setdiff(names(design_df_complete), c("patient_ID", "time", "status"))
    if (input$surv_outcome %in% predictor_cols) {
      stop(paste0(
        "The outcome variable '", input$surv_outcome, "' cannot also be used as a predictor. ",
        "This creates perfect collinearity. Please select different predictors."
      ))
    }

    if (length(predictor_cols) == 0) {
      stop("No predictors available after filtering.")
    }

    # Check for high correlation between continuous predictors and outcome
    warnings <- character(0)
    for (pred in predictor_cols) {
      pred_vals <- design_df_complete[[pred]]
      
      # For continuous predictors, check correlation with outcome
      if (is.numeric(pred_vals) && length(unique(pred_vals)) > 10) {
        cor_val <- tryCatch(
          cor(pred_vals, design_df_complete$time, use = "complete.obs"),
          error = function(e) NA
        )
        if (!is.na(cor_val) && abs(cor_val) > 0.95) {
          warnings <- c(warnings, sprintf(
            "WARNING: Predictor '%s' has very high correlation with outcome (r = %.3f). This may cause numerical instability.",
            pred, cor_val
          ))
        }
      }
      
      # For categorical predictors, check if groups perfectly separate outcome ranges
      if (is.factor(pred_vals) || is.character(pred_vals)) {
        groups <- split(design_df_complete$time, pred_vals)
        if (length(groups) >= 2) {
          # Check if group ranges don't overlap (perfect separation)
          group_ranges <- lapply(groups, function(x) c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)))
          max_mins <- max(sapply(group_ranges, function(x) x[1]))
          min_maxs <- min(sapply(group_ranges, function(x) x[2]))
          
          if (max_mins > min_maxs) {
            warnings <- c(warnings, sprintf(
              "WARNING: Categorical predictor '%s' creates perfect or near-perfect separation of outcome values. This may cause convergence issues.",
              pred
            ))
          } else {
            # Check for very high eta-squared (proportion of variance explained)
            fit_temp <- tryCatch(
              lm(design_df_complete$time ~ pred_vals),
              error = function(e) NULL
            )
            if (!is.null(fit_temp)) {
              r_squared <- summary(fit_temp)$r.squared
              if (r_squared > 0.95) {
                warnings <- c(warnings, sprintf(
                  "WARNING: Categorical predictor '%s' explains %.1f%% of outcome variance. This creates near-perfect correlation and may cause numerical instability.",
                  pred, r_squared * 100
                ))
              }
            }
          }
        }
      }
    }

    formula_str <- paste("Surv(time, status) ~", paste(predictor_cols, collapse = " + "))
    formula_obj <- as.formula(formula_str)

    # Fit Cox proportional hazards model
    cox_model <- tryCatch(
      {
        survival::coxph(formula_obj, data = design_df_complete)
      },
      warning = function(w) {
        warnings <<- c(warnings, paste("Cox model warning:", conditionMessage(w)))
        survival::coxph(formula_obj, data = design_df_complete)
      },
      error = function(e) {
        stop(paste0(
          "Cox model fitting failed: ", conditionMessage(e), "\n\n",
          "This often occurs when predictors are too highly correlated with the outcome, ",
          "or when there's perfect separation in categorical predictors. ",
          "Try using different predictors that are not derived from the outcome variable."
        ))
      }
    )

    # Test proportional hazards assumption (Schoenfeld residuals)
    ph_test <- tryCatch(
      {
        survival::cox.zph(cox_model)
      },
      error = function(e) NULL
    )

    # Continuous model results
    cox_summary <- summary(cox_model)

    # Calculate risk scores (linear predictor) for each sample
    risk_scores <- predict(cox_model, type = "lp")  # linear predictor
    design_df_complete$risk_score <- risk_scores

    # Calculate threshold based on risk scores using selected split method
    threshold <- if (split_method == "mean") {
      mean(risk_scores, na.rm = TRUE)
    } else {
      median(risk_scores, na.rm = TRUE)
    }

    # Dichotomize based on RISK SCORE
    # Higher risk score = higher predicted hazard = worse prognosis
    design_df_complete$risk_group <- ifelse(design_df_complete$risk_score > threshold, "High Risk", "Low Risk")
    
    # Check that we have at least 2 groups after risk stratification
    n_groups <- length(unique(design_df_complete$risk_group))
    if (n_groups < 2) {
      # This should rarely happen with median/mean split, but handle edge case
      # Try median if mean failed
      threshold <- median(risk_scores, na.rm = TRUE)
      design_df_complete$risk_group <- ifelse(design_df_complete$risk_score >= threshold, "High Risk", "Low Risk")
      
      n_groups <- length(unique(design_df_complete$risk_group))
      if (n_groups < 2) {
        stop("Cannot create two risk groups. All samples have identical risk scores.")
      }
      
      warnings <- c(warnings, sprintf(
        "WARNING: Mean split resulted in only 1 group. Used median split instead (threshold = %.3f).",
        threshold
      ))
    }
    
    # Also keep outcome-based dichotomization for the logistic model
    # Use median of outcome times for this
    outcome_threshold <- median(design_df_complete$time, na.rm = TRUE)
    design_df_complete$outcome_high <- ifelse(design_df_complete$time > outcome_threshold, 1, 0)

    # Prepare dichotomized formula - replace time with outcome_high
    formula_dichot_str <- paste("outcome_high ~", paste(predictor_cols, collapse = " + "))
    formula_dichot <- as.formula(formula_dichot_str)

    # Fit logistic regression for dichotomized outcome
    dichot_model <- glm(formula_dichot, data = design_df_complete, family = binomial)

    return(list(
      analysis_mode = "multivariate",
      cox_model = cox_model,
      cox_summary = cox_summary,
      ph_test = ph_test,
      dichot_model = dichot_model,
      data = design_df_complete,
      threshold = threshold,
      split_method = split_method,
      predictor_cols = predictor_cols,
      name_mapping = if (exists("name_mapping")) name_mapping else NULL,
      warnings = warnings,
      details = list(
        samples_before = n_before,
        samples_after = n_after,
        samples_dropped = n_dropped,
        dropped_ids = dropped_ids
      )
    ))
  }
  
  # Univariate analysis - run separate model for each predictor
  run_surv_univariate <- function() {
    req(rv$meta_sample, rv$abundance_sample, input$surv_outcome, input$surv_predictors)

    meta_patient <- rv$meta_sample

    # Validate outcome variable
    if (!(input$surv_outcome %in% colnames(meta_patient))) {
      stop(paste0("Outcome '", input$surv_outcome, "' not found in metadata."))
    }

    outcome_raw <- meta_patient[[input$surv_outcome]]
    outcome <- suppressWarnings(as.numeric(outcome_raw))

    if (all(is.na(outcome))) {
      stop(paste0("Outcome '", input$surv_outcome, "' cannot be converted to numeric."))
    }

    # Get split method
    split_method <- input$surv_split_method %||% "median"

    # Build full predictor matrix (same as multivariate)
    predictor_list <- input$surv_predictors
    design_df <- data.frame(patient_ID = meta_patient$patient_ID)
    name_mapping <- NULL

    # Add cluster predictors if selected
    if ("cluster" %in% predictor_list) {
      abund <- rv$abundance_sample

      # Filter to selected clusters
      if (!is.null(input$surv_cluster_subset) && length(input$surv_cluster_subset) > 0) {
        selected_clusters <- input$surv_cluster_subset
        keep_cols <- colnames(abund)[colnames(abund) %in% selected_clusters]
        if (length(keep_cols) == 0) {
          stop("No valid clusters selected.")
        }
        abund <- abund[, keep_cols, drop = FALSE]
      }

      abund_df <- as.data.frame(abund)
      abund_df$patient_ID <- rownames(abund_df)
      
      # Sanitize cluster column names
      cluster_cols <- setdiff(names(abund_df), "patient_ID")
      name_mapping <- setNames(cluster_cols, cluster_cols)
      
      for (col in cluster_cols) {
        new_name <- col
        new_name <- gsub(" ", "_", new_name)
        new_name <- gsub("\\+", "p", new_name)
        new_name <- gsub("-", "n", new_name)
        new_name <- gsub("[^A-Za-z0-9_]", "_", new_name)
        if (!grepl("^[A-Za-z_]", new_name)) {
          new_name <- paste0("X_", new_name)
        }
        if (new_name != col) {
          names(abund_df)[names(abund_df) == col] <- new_name
          name_mapping[new_name] <- col
        }
      }
      
      design_df <- merge(design_df, abund_df, by = "patient_ID", all.x = TRUE)
    }

    # Add metadata predictors
    meta_predictors <- setdiff(predictor_list, "cluster")
    if (length(meta_predictors) > 0) {
      meta_subset <- meta_patient[, c("patient_ID", meta_predictors), drop = FALSE]
      design_df <- merge(design_df, meta_subset, by = "patient_ID", all.x = TRUE)
    }

    # Attach outcome
    design_df$time <- outcome

    # Remove rows with missing values
    n_before <- nrow(design_df)
    design_df_complete <- design_df[complete.cases(design_df), ]
    n_after <- nrow(design_df_complete)
    n_dropped <- n_before - n_after
    dropped_ids <- setdiff(design_df$patient_ID, design_df_complete$patient_ID)

    if (n_after < 10) {
      stop("Insufficient complete cases after removing missing values (n < 10).")
    }

    # Add status indicator
    design_df_complete$status <- 1

    # Get all predictor columns
    all_predictor_cols <- setdiff(names(design_df_complete), c("patient_ID", "time", "status"))
    
    if (input$surv_outcome %in% all_predictor_cols) {
      stop(paste0(
        "The outcome variable '", input$surv_outcome, 
        "' is also selected as a predictor. This creates circular logic. ",
        "Please remove it from the predictor list."
      ))
    }

    # Run separate Cox model for each predictor
    univar_results <- list()
    
    for (pred_col in all_predictor_cols) {
      tryCatch({
        # Build formula for this predictor only
        formula_str <- paste("Surv(time, status) ~", pred_col)
        formula_obj <- as.formula(formula_str)
        
        # Fit Cox model
        cox_model <- survival::coxph(formula_obj, data = design_df_complete)
        cox_summary <- summary(cox_model)
        
        # Test proportional hazards
        ph_test <- tryCatch(survival::cox.zph(cox_model), error = function(e) NULL)
        
        # Calculate risk scores
        risk_scores <- predict(cox_model, type = "lp")
        
        # Calculate threshold
        threshold <- if (split_method == "mean") {
          mean(risk_scores, na.rm = TRUE)
        } else {
          median(risk_scores, na.rm = TRUE)
        }
        
        # Create risk groups
        temp_df <- design_df_complete[, c("patient_ID", "time", "status")]
        temp_df$risk_score <- risk_scores
        temp_df$risk_group <- ifelse(risk_scores > threshold, "High Risk", "Low Risk")
        
        # Check for 2 groups
        if (length(unique(temp_df$risk_group)) < 2) {
          threshold <- median(risk_scores, na.rm = TRUE)
          temp_df$risk_group <- ifelse(risk_scores >= threshold, "High Risk", "Low Risk")
        }
        
        # Store results for this predictor
        univar_results[[pred_col]] <- list(
          predictor = pred_col,
          cox_model = cox_model,
          cox_summary = cox_summary,
          ph_test = ph_test,
          data = temp_df,
          threshold = threshold,
          split_method = split_method
        )
      }, error = function(e) {
        # Skip predictors that fail
        NULL
      })
    }
    
    if (length(univar_results) == 0) {
      stop("No predictors could be successfully modeled.")
    }

    return(list(
      analysis_mode = "univariate",
      univar_results = univar_results,
      predictor_names = names(univar_results),
      name_mapping = name_mapping,
      details = list(
        samples_before = n_before,
        samples_after = n_after,
        samples_dropped = n_dropped,
        dropped_ids = dropped_ids
      )
    ))
  }

  observeEvent(input$run_surv, {
    result <- tryCatch(
      {
        run_surv()
      },
      error = function(e) {
        list(error = TRUE, message = conditionMessage(e))
      }
    )

    if (!is.null(result) && isTRUE(result$error)) {
      surv_state(result)
      showNotification(paste("Time to event analysis failed:", result$message), type = "error", duration = 10)
      return()
    }

    surv_state(result)
    output$surv_cleared_msg <- renderText(NULL)
    output$surv_error_ui <- renderUI(NULL)
    
    # If univariate mode, populate the predictor selector dropdown
    if (!is.null(result$analysis_mode) && result$analysis_mode == "univariate") {
      # Get display names for predictors
      display_names <- result$predictor_names
      if (!is.null(result$name_mapping)) {
        display_names <- sapply(result$predictor_names, function(pn) {
          if (pn %in% names(result$name_mapping)) {
            result$name_mapping[[pn]]
          } else {
            pn
          }
        }, USE.NAMES = FALSE)
      }
      
      # Create named vector for dropdown
      predictor_choices <- setNames(result$predictor_names, display_names)
      updateSelectInput(session, "surv_univar_predictor_display", 
                       choices = predictor_choices,
                       selected = result$predictor_names[1])
    }
    
    # Display warnings if any (multivariate mode only)
    if (!is.null(result$warnings) && length(result$warnings) > 0) {
      showNotification(
        HTML(paste(result$warnings, collapse = "<br><br>")),
        type = "warning",
        duration = 5,  # Auto-dismiss after 5 seconds
        closeButton = TRUE
      )
    }
  })

  observeEvent(input$reset_surv, {
    surv_state(NULL)
    showNotification("Time to event analysis results cleared.", type = "message", duration = 5)
    output$surv_cleared_msg <- renderText("Results cleared. Run time to event analysis again to see results here.")
    output$surv_error_ui <- renderUI(NULL)
  })
  
  # Helper function to get the current result to display
  # In multivariate mode: returns the full result
  # In univariate mode: returns the selected predictor's result
  get_current_result <- function() {
    s <- surv_state()
    req(s)
    
    if (!is.null(s$error) && isTRUE(s$error)) {
      return(s)
    }
    
    if (!is.null(s$analysis_mode) && s$analysis_mode == "univariate") {
      # Get selected predictor
      selected_pred <- input$surv_univar_predictor_display
      if (is.null(selected_pred) || !selected_pred %in% names(s$univar_results)) {
        selected_pred <- s$predictor_names[1]  # Default to first
      }
      
      # Extract that predictor's results and add metadata
      pred_result <- s$univar_results[[selected_pred]]
      pred_result$analysis_mode <- "univariate"
      pred_result$selected_predictor <- selected_pred
      pred_result$name_mapping <- s$name_mapping
      pred_result$details <- s$details
      pred_result$predictor_cols <- selected_pred
      return(pred_result)
    } else {
      # Multivariate mode - return as is
      return(s)
    }
  }

  # Render error message
  output$surv_error_ui <- renderUI({
    s <- surv_state()
    req(s)

    if (!is.null(s$error) && isTRUE(s$error)) {
      div(
        style = "color: red; border: 2px solid red; padding: 10px; margin: 10px;",
        h4("Error in Time to Event Analysis:"),
        p(s$message),
        h5("Suggestions:"),
        tags$ul(
          tags$li("Ensure sufficient samples with complete data (n >= 10)"),
          tags$li("Check that manual threshold is within data range"),
          tags$li("Try different predictors or outcome variable"),
          tags$li("Verify that outcome variable is continuous")
        )
      )
    }
  })

  # Time to event curve plot
  output$surv_curve_plot <- renderPlot({
    s <- get_current_result()
    req(s)

    if (!is.null(s$error) && isTRUE(s$error)) {
      plot.new()
      text(0.5, 0.5, paste("Error:", s$message), col = "red", cex = 1.2)
      return()
    }

    req(s$data, s$threshold)

    # Create risk-based groups for time to event curve
    surv_data <- s$data
    surv_data$group <- factor(surv_data$risk_group, levels = c("Low Risk", "High Risk"))
    
    # Check how many groups we actually have
    n_groups <- length(levels(factor(surv_data$group)))
    actual_groups <- unique(as.character(surv_data$group))
    
    # Handle case where we only have 1 group
    if (length(actual_groups) < 2) {
      plot.new()
      text(0.5, 0.5, "Cannot create time-to-event curves: only 1 risk group present.\nTry using different predictors or a different split method.", 
           col = "red", cex = 1.1)
      return()
    }

    # Fit time to event curves by risk group
    fit <- survival::survfit(Surv(time, status) ~ group, data = surv_data)

    # Build plot title
    plot_title <- if (!is.null(s$analysis_mode) && s$analysis_mode == "univariate" && !is.null(s$selected_predictor)) {
      # Get display name for predictor
      display_name <- s$selected_predictor
      if (!is.null(s$name_mapping) && s$selected_predictor %in% names(s$name_mapping)) {
        display_name <- s$name_mapping[[s$selected_predictor]]
      }
      sprintf("Time to Event Curves: %s (%s split)", display_name, tools::toTitleCase(s$split_method %||% "median"))
    } else {
      sprintf("Time to Event Curves by Predicted Risk (%s split)", tools::toTitleCase(s$split_method %||% "median"))
    }
    
    # Only specify legend.labs if we have exactly 2 groups
    plot_args <- list(
      fit = fit,
      data = surv_data,
      pval = TRUE,
      pval.method = TRUE,
      pval.method.coord = c(0.05, 0.15),
      pval.coord = c(0.05, 0.2),
      pval.hjust = 0.5,
      conf.int = isTRUE(input$surv_show_ci),
      risk.table = FALSE,
      xlab = "Time to Event",
      ylab = "Event-Free Probability",
      title = plot_title,
      subtitle = "P-value: Log-rank test (dichotomized groups)",
      legend.title = "Risk Group",
      palette = c("#00BA38", "#F8766D"),
      ggtheme = theme_minimal(base_size = 14)
    )
    
    # Only add legend.labs if we have exactly 2 groups
    if (length(actual_groups) == 2) {
      plot_args$legend.labs <- c("Low Risk", "High Risk")
    }
    
    do.call(ggsurvplot, plot_args)$plot
  })

  # Model summary
  output$surv_summary <- renderPrint({
    s <- get_current_result()
    req(s)

    if (!is.null(s$error) && isTRUE(s$error)) {
      cat("Time to Event Analysis Error:\n\n")
      cat(s$message, "\n")
      return()
    }

    det <- s$details

    cat("=== Time to Event Analysis Summary ===\n\n")
    cat("Model: Cox Proportional Hazards\n")
    
    if (!is.null(s$analysis_mode) && s$analysis_mode == "univariate") {
      cat("Analysis Mode: Univariate\n")
      # Display original name if available
      display_name <- s$selected_predictor
      if (!is.null(s$name_mapping) && s$selected_predictor %in% names(s$name_mapping)) {
        display_name <- s$name_mapping[[s$selected_predictor]]
      }
      cat("Predictor:", display_name, "\n")
    } else {
      cat("Analysis Mode: Multivariate\n")
      cat("Predictors:", paste(input$surv_predictors, collapse = ", "), "\n")
    }
    
    cat("Outcome:", input$surv_outcome, "\n")
    cat("Risk stratification:", tools::toTitleCase(s$split_method %||% "median"), "split of risk scores\n")
    cat("Risk score threshold:", sprintf("%.3f\n", s$threshold))
    
    # Add risk group counts
    if (!is.null(s$data$risk_group)) {
      group_counts <- table(s$data$risk_group)
      cat("Risk groups: Low Risk n=", group_counts["Low Risk"], 
          "; High Risk n=", group_counts["High Risk"], "\n", sep="")
    }
    
    cat("\nSample counts:\n")
    cat("  Before filtering:", det$samples_before, "\n")
    cat("  After filtering:", det$samples_after, "\n")
    cat("  Dropped:", det$samples_dropped, "\n")
    if (length(det$dropped_ids) > 0) {
      cat("  Dropped IDs:", paste(head(det$dropped_ids, 10), collapse = ", "))
      if (length(det$dropped_ids) > 10) cat(", ...")
      cat("\n")
    }
    
    # Display warnings if any
    if (!is.null(s$warnings) && length(s$warnings) > 0) {
      cat("\n=== WARNINGS ===\n")
      for (i in seq_along(s$warnings)) {
        cat(sprintf("%d. %s\n", i, s$warnings[i]))
      }
      cat("\nNote: High correlation between predictors and outcome can cause\n")
      cat("numerical instability. Consider using predictors that are measured\n")
      cat("independently of the outcome variable.\n")
    }
  })

  # Performance metrics table
  build_surv_perf_table <- function(s, include_icons = TRUE) {
    req(s, s$cox_model)

    cox_summary <- s$cox_summary

    # Extract coefficients for continuous model
    coef_matrix <- cox_summary$coefficients
    
    # Build performance table
    perf_rows <- list()

    # Add validation warnings/status as FIRST row(s)
    if (!is.null(s$warnings) && length(s$warnings) > 0) {
      # Add a header row for warnings
      warning_value <- if (include_icons) "⚠️ WARNINGS DETECTED" else "WARNINGS DETECTED"
      perf_rows[[length(perf_rows) + 1]] <- data.frame(
        Metric = "Model Validity Check",
        Value = warning_value,
        P_value = "---",
        Interpretation = "See warnings below - model may have numerical issues",
        stringsAsFactors = FALSE
      )
      
      # Add each warning as a separate row
      for (i in seq_along(s$warnings)) {
        perf_rows[[length(perf_rows) + 1]] <- data.frame(
          Metric = sprintf("  Warning %d", i),
          Value = "---",
          P_value = "---",
          Interpretation = s$warnings[i],
          stringsAsFactors = FALSE
        )
      }
    } else {
      # No warnings - model is valid
      pass_value <- if (include_icons) "✓ PASSED" else "PASSED"
      perf_rows[[length(perf_rows) + 1]] <- data.frame(
        Metric = "Model Validity Check",
        Value = pass_value,
        P_value = "---",
        Interpretation = "No high correlation or separation issues detected between predictors and outcome",
        stringsAsFactors = FALSE
      )
    }

    # Check for NA coefficients in continuous model (will show in Model Coefficients table)
    # Don't add individual predictor metrics to Performance Metrics table

    # Proportional hazards test (Schoenfeld residuals)
    if (!is.null(s$ph_test)) {
      global_p <- s$ph_test$table["GLOBAL", "p"]
      perf_rows[[length(perf_rows) + 1]] <- data.frame(
        Metric = "Schoenfeld Test (Global)",
        Value = "---",
        P_value = if (global_p < 0.001) formatC(global_p, format = "e", digits = 2) else sprintf("%.4f", global_p),
        Interpretation = "Tests proportional hazards assumption; p > 0.05 suggests assumption holds",
        stringsAsFactors = FALSE
      )
      
      # Count predictors failing Schoenfeld test (p < 0.05)
      n_predictors <- nrow(s$ph_test$table) - 1  # exclude GLOBAL row
      failed_predictors <- sum(s$ph_test$table[1:n_predictors, "p"] < 0.05, na.rm = TRUE)
      
      interpretation_text <- if (failed_predictors == 0) {
        "All predictors satisfy proportional hazards assumption"
      } else if (failed_predictors == 1) {
        "1 predictor violates proportional hazards assumption (p < 0.05). Consider removing this predictor or using time-varying coefficients. See Model Feature Performance table for details."
      } else {
        sprintf("%d predictors violate proportional hazards assumption (p < 0.05). Consider removing these predictors or using time-varying coefficients. See Model Feature Performance table for details.", failed_predictors)
      }
      
      status_value <- if (include_icons) {
        if (failed_predictors == 0) "✓ 0 Failed" else sprintf("⚠️ %d Failed", failed_predictors)
      } else {
        if (failed_predictors == 0) "0 Failed" else sprintf("%d Failed", failed_predictors)
      }
      
      perf_rows[[length(perf_rows) + 1]] <- data.frame(
        Metric = "Predictors Failing Schoenfeld Test",
        Value = status_value,
        P_value = "---",
        Interpretation = interpretation_text,
        stringsAsFactors = FALSE
      )
      # Individual Schoenfeld tests are now in Model Coefficients table
    } else {
      perf_rows[[length(perf_rows) + 1]] <- data.frame(
        Metric = "Schoenfeld Test",
        Value = "---",
        P_value = "N/A",
        Interpretation = "Could not compute proportional hazards test",
        stringsAsFactors = FALSE
      )
    }

    # Add Cox model p-value (for continuous predictor effect)
    if (!is.null(s$predictor_cols) && length(s$predictor_cols) > 0) {
      # Get p-value from Cox model for the predictor (univariate has 1 predictor)
      coef_matrix <- s$cox_summary$coefficients
      if (nrow(coef_matrix) == 1) {
        # Univariate mode - single predictor
        cox_p <- coef_matrix[1, "Pr(>|z|)"]
        perf_rows[[length(perf_rows) + 1]] <- data.frame(
          Metric = "Cox Model P-value (Continuous)",
          Value = "---",
          P_value = if (cox_p < 0.001) formatC(cox_p, format = "e", digits = 2) else sprintf("%.4f", cox_p),
          Interpretation = "Tests whether predictor significantly affects hazard as continuous variable in Cox model",
          stringsAsFactors = FALSE
        )
      }
    }

    # Add dichotomized group counts - one row per group
    group_counts <- table(s$data$risk_group)
    
    # Check if we have both groups
    low_count <- if ("Low Risk" %in% names(group_counts)) as.character(group_counts["Low Risk"]) else "0"
    high_count <- if ("High Risk" %in% names(group_counts)) as.character(group_counts["High Risk"]) else "0"
    
    # Low Risk group
    perf_rows[[length(perf_rows) + 1]] <- data.frame(
      Metric = "Risk Group: Low Risk",
      Value = low_count,
      P_value = "---",
      Interpretation = "Number of samples predicted as low risk based on model",
      stringsAsFactors = FALSE
    )
    
    # High Risk group
    perf_rows[[length(perf_rows) + 1]] <- data.frame(
      Metric = "Risk Group: High Risk",
      Value = high_count,
      P_value = "---",
      Interpretation = "Number of samples predicted as high risk based on model",
      stringsAsFactors = FALSE
    )

    # Log-rank test for dichotomized curves - only if we have 2 groups
    logrank_test <- tryCatch(
      {
        if (length(unique(s$data$risk_group)) >= 2) {
          survdiff_result <- survival::survdiff(Surv(time, status) ~ risk_group, data = s$data)
          logrank_p <- 1 - pchisq(survdiff_result$chisq, length(survdiff_result$n) - 1)
          list(success = TRUE, p_value = logrank_p)
        } else {
          list(success = FALSE, message = "Only 1 risk group present")
        }
      },
      error = function(e) {
        list(success = FALSE, message = conditionMessage(e))
      }
    )
    
    if (logrank_test$success) {
      perf_rows[[length(perf_rows) + 1]] <- data.frame(
        Metric = "Log-rank Test P-value (Dichotomized)",
        Value = "---",
        P_value = if (logrank_test$p_value < 0.001) formatC(logrank_test$p_value, format = "e", digits = 2) else sprintf("%.4f", logrank_test$p_value),
        Interpretation = "Tests curve separation after dichotomizing by risk score. May differ from Cox p-value due to information loss from dichotomization. This p-value appears on the plot.",
        stringsAsFactors = FALSE
      )
    } else {
      perf_rows[[length(perf_rows) + 1]] <- data.frame(
        Metric = "Log-rank Test (Risk-based Curves)",
        Value = "---",
        P_value = "N/A",
        Interpretation = paste("Cannot compute:", logrank_test$message),
        stringsAsFactors = FALSE
      )
    }

    # Sample info
    det <- s$details
    perf_rows[[length(perf_rows) + 1]] <- data.frame(
      Metric = "Samples before filtering",
      Value = as.character(det$samples_before),
      P_value = "---",
      Interpretation = "Count of samples before removing missing values",
      stringsAsFactors = FALSE
    )

    perf_rows[[length(perf_rows) + 1]] <- data.frame(
      Metric = "Samples after filtering",
      Value = as.character(det$samples_after),
      P_value = "---",
      Interpretation = "Count of samples used in analysis",
      stringsAsFactors = FALSE
    )

    perf_rows[[length(perf_rows) + 1]] <- data.frame(
      Metric = "Samples dropped",
      Value = as.character(det$samples_dropped),
      P_value = "---",
      Interpretation = "Number of samples removed due to missing values",
      stringsAsFactors = FALSE
    )

    perf_table <- do.call(rbind, perf_rows)
    
    # Wrap interpretation text to make it less wide (only at whitespace)
    perf_table$Interpretation <- sapply(perf_table$Interpretation, function(text) {
      paste(strwrap(text, width = 50), collapse = "\n")
    }, USE.NAMES = FALSE)
    
    perf_table
  }

  output$surv_perf_table <- renderTable({
    s <- get_current_result()
    req(s)

    if (!is.null(s$error) && isTRUE(s$error)) {
      return(data.frame(Message = "Error in analysis - see message above"))
    }

    build_surv_perf_table(s)
  }, rownames = FALSE, striped = TRUE)

  # Model coefficients table
  surv_coef_table <- reactive({
    include_all_cols = TRUE  # default to all columns for export
  }) 
  
  surv_coef_table <- function(include_all_cols = TRUE) {
    s <- isolate(surv_state())  # Use full state, not filtered result
    req(s)

    if (!is.null(s$error) && isTRUE(s$error)) {
      return(NULL)
    }
    
    # Check if univariate mode
    if (!is.null(s$analysis_mode) && s$analysis_mode == "univariate") {
      # Build table with all univariate results
      all_results <- list()
      
      for (pred_name in names(s$univar_results)) {
        pred_result <- s$univar_results[[pred_name]]
        coef_matrix <- pred_result$cox_summary$coefficients
        
        # Get display name
        display_name <- pred_name
        if (!is.null(s$name_mapping) && pred_name %in% names(s$name_mapping)) {
          display_name <- s$name_mapping[[pred_name]]
        }
        
        # Calculate log-rank p-value for this predictor
        logrank_p <- tryCatch({
          if (length(unique(pred_result$data$risk_group)) >= 2) {
            survdiff_result <- survival::survdiff(Surv(time, status) ~ risk_group, data = pred_result$data)
            1 - pchisq(survdiff_result$chisq, length(survdiff_result$n) - 1)
          } else {
            NA_real_
          }
        }, error = function(e) NA_real_)
        
        # Extract coefficients
        all_results[[pred_name]] <- data.frame(
          Feature = display_name,
          Coefficient = coef_matrix[1, "coef"],
          HazardRatio = coef_matrix[1, "exp(coef)"],
          SE = coef_matrix[1, "se(coef)"],
          Z_value = coef_matrix[1, "z"],
          P_value_Cox = coef_matrix[1, "Pr(>|z|)"],
          P_value_Logrank = logrank_p,
          Schoenfeld_P = if (!is.null(pred_result$ph_test)) {
            pred_result$ph_test$table[1, "p"]
          } else {
            NA_real_
          },
          stringsAsFactors = FALSE
        )
      }
      
      coef_df <- do.call(rbind, all_results)
      rownames(coef_df) <- NULL
      
      # Format p-values
      coef_df$P_value_Cox <- sapply(coef_df$P_value_Cox, function(p) {
        if (is.na(p)) return(NA_character_)
        if (p < 0.001) formatC(p, format = "e", digits = 2) else sprintf("%.4f", p)
      })
      
      coef_df$P_value_Logrank <- sapply(coef_df$P_value_Logrank, function(p) {
        if (is.na(p)) return(NA_character_)
        if (p < 0.001) formatC(p, format = "e", digits = 2) else sprintf("%.4f", p)
      })
      
      coef_df$Schoenfeld_P <- sapply(coef_df$Schoenfeld_P, function(p) {
        if (is.na(p)) return(NA_character_)
        if (p < 0.001) formatC(p, format = "e", digits = 2) else sprintf("%.4f", p)
      })
      
      # Order by absolute coefficient
      coef_df <- coef_df[order(-abs(coef_df$Coefficient)), ]
      
      # Remove technical columns if requested
      if (!include_all_cols) {
        coef_df <- coef_df[, !(names(coef_df) %in% c("Coefficient", "SE", "Z_value"))]
      }
      
      return(coef_df)
    }
    
    # Multivariate mode - use get_current_result
    s <- get_current_result()
    req(s, s$cox_model)

    coef_matrix <- s$cox_summary$coefficients

    # Store sanitized names before restoration (these match Schoenfeld test results)
    sanitized_names <- rownames(coef_matrix)
    
    # Restore original names using name_mapping if available
    feature_names <- rownames(coef_matrix)
    if (!is.null(s$name_mapping)) {
      feature_names <- sapply(feature_names, function(fn) {
        if (fn %in% names(s$name_mapping)) {
          s$name_mapping[[fn]]
        } else {
          fn
        }
      }, USE.NAMES = FALSE)
    }
    
    coef_df <- data.frame(
      Feature = feature_names,
      Coefficient = coef_matrix[, "coef"],
      HazardRatio = coef_matrix[, "exp(coef)"],
      SE = coef_matrix[, "se(coef)"],
      Z_value = coef_matrix[, "z"],
      P_value_Cox = coef_matrix[, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )

    # Add Schoenfeld test p-values if available (use sanitized names for lookup)
    if (!is.null(s$ph_test)) {
      schoenfeld_p <- rep(NA_real_, nrow(coef_df))
      for (i in 1:nrow(coef_df)) {
        san_name <- sanitized_names[i]  # Use sanitized name for lookup
        if (san_name %in% rownames(s$ph_test$table)) {
          schoenfeld_p[i] <- s$ph_test$table[san_name, "p"]
        }
      }
      coef_df$Schoenfeld_P <- schoenfeld_p
    } else {
      coef_df$Schoenfeld_P <- NA_real_
    }

    # Format p-values
    coef_df$P_value_Cox <- sapply(coef_df$P_value_Cox, function(p) {
      if (is.na(p)) return(NA_character_)
      if (p < 0.001) formatC(p, format = "e", digits = 2) else sprintf("%.4f", p)
    })
    
    coef_df$Schoenfeld_P <- sapply(coef_df$Schoenfeld_P, function(p) {
      if (is.na(p)) return(NA_character_)
      if (p < 0.001) formatC(p, format = "e", digits = 2) else sprintf("%.4f", p)
    })

    # Order by absolute coefficient
    coef_df <- coef_df[order(-abs(coef_df$Coefficient)), ]

    # Remove technical columns if requested
    if (!include_all_cols) {
      coef_df <- coef_df[, !(names(coef_df) %in% c("Coefficient", "SE", "Z_value"))]
    }
    
    coef_df
  }

  output$surv_coefs <- renderTable({
    surv_coef_table(include_all_cols = FALSE)
  })
  
  # Dynamic header for coefficients table
  output$surv_coef_header <- renderUI({
    s <- surv_state()
    if (!is.null(s) && !is.null(s$analysis_mode) && s$analysis_mode == "univariate") {
      tagList(
        h4("Model Feature Performance (univariate)"),
        helpText("Showing univariate Cox model results for each predictor independently. Each predictor was tested in a separate model.")
      )
    } else {
      tagList(
        h4("Model Feature Performance (multivariate)"),
        helpText("Showing results from a single model with all predictors included together.")
      )
    }
  })

  output$hasSurvResults <- reactive({
    s <- surv_state()
    !is.null(s) && is.null(s$error)
  })
  outputOptions(output, "hasSurvResults", suspendWhenHidden = FALSE)

  output$hasSurvCoefs <- reactive({
    s <- surv_state()
    if (is.null(s) || !is.null(s$error)) {
      return(FALSE)
    }
    # Check for multivariate mode or univariate mode
    if (!is.null(s$analysis_mode) && s$analysis_mode == "univariate") {
      return(!is.null(s$univar_results) && length(s$univar_results) > 0)
    } else {
      return(!is.null(s$cox_model))
    }
  })
  outputOptions(output, "hasSurvCoefs", suspendWhenHidden = FALSE)

  # ZIP download handler
  output$export_surv_zip <- downloadHandler(
    filename = function() {
      subset_id <- rv$subset_id %||% "000000000"
      outcome <- gsub("\\s+", "_", tolower(input$surv_outcome %||% "outcome"))
      s <- surv_state()
      mode_suffix <- if (!is.null(s$analysis_mode) && s$analysis_mode == "univariate") "univar" else "multivar"
      
      # Add predictor name for univariate mode
      if (!is.null(s$analysis_mode) && s$analysis_mode == "univariate" && !is.null(input$surv_univar_predictor_display)) {
        predictor_name <- input$surv_univar_predictor_display
        # Get display name if available
        if (!is.null(s$name_mapping) && predictor_name %in% names(s$name_mapping)) {
          predictor_name <- s$name_mapping[[predictor_name]]
        }
        # Sanitize for filename
        predictor_name <- gsub("[^A-Za-z0-9_-]", "_", predictor_name)
        paste0("time_to_event_cox_", outcome, "_", mode_suffix, "_", predictor_name, "_", subset_id, ".zip")
      } else {
        paste0("time_to_event_cox_", outcome, "_", mode_suffix, "_", subset_id, ".zip")
      }
    },
    content = function(file) {
      subset_id <- rv$subset_id %||% "000000000"
      s <- get_current_result()
      req(s, s$cox_model)

      outcome <- gsub("\\s+", "_", tolower(input$surv_outcome %||% "outcome"))
      
      # Determine mode and predictors
      analysis_mode <- if (!is.null(s$analysis_mode)) s$analysis_mode else "multivariate"
      mode_suffix <- if (analysis_mode == "univariate") "univar" else "multivar"
      
      # Get predictor text and sanitized name for filename
      predictor_text <- if (analysis_mode == "univariate") {
        display_name <- s$selected_predictor
        if (!is.null(s$name_mapping) && s$selected_predictor %in% names(s$name_mapping)) {
          display_name <- s$name_mapping[[s$selected_predictor]]
        }
        display_name
      } else {
        paste(input$surv_predictors, collapse = "; ")
      }
      
      # Build prefix with predictor name for univariate
      if (analysis_mode == "univariate") {
        predictor_filename <- gsub("[^A-Za-z0-9_-]", "_", predictor_text)
        prefix <- paste0("time_to_event_cox_", outcome, "_", mode_suffix, "_", predictor_filename)
      } else {
        prefix <- paste0("time_to_event_cox_", outcome, "_", mode_suffix)
      }

      tmpdir <- tempdir()
      files <- c()

      # 1. Model Summary (key-value CSV)
      summary_file <- file.path(tmpdir, paste0(prefix, "_model_summary_", subset_id, ".csv"))
      summary_kv <- data.frame(
        Key = c(
          "Model type", "Analysis mode", "Outcome variable", "Predictors", "Risk stratification method", 
          "Risk score threshold", "Samples before filtering", "Samples after filtering", "Samples dropped"
        ),
        Value = c(
          "Cox Proportional Hazards",
          tools::toTitleCase(analysis_mode),
          input$surv_outcome,
          predictor_text,
          tools::toTitleCase(s$split_method %||% "median"),
          sprintf("%.3f", s$threshold),
          s$details$samples_before %||% NA,
          s$details$samples_after %||% NA,
          s$details$samples_dropped %||% NA
        ),
        stringsAsFactors = FALSE
      )
      write.csv(summary_kv, summary_file, row.names = FALSE)
      files <- c(files, summary_file)

      # 2. Performance Metrics
      perf_file <- file.path(tmpdir, paste0(prefix, "_performance_metrics_", subset_id, ".csv"))
      perf_df <- tryCatch(
        {
          build_surv_perf_table(s, include_icons = FALSE)
        },
        error = function(e) data.frame(Message = "Error extracting performance metrics")
      )
      write.csv(perf_df, perf_file, row.names = FALSE)
      files <- c(files, perf_file)

      # 3. Model Coefficients
      coef_file <- file.path(tmpdir, paste0(prefix, "_model_coefficients_", subset_id, ".csv"))
      coef_df <- tryCatch(
        {
          surv_coef_table()
        },
        error = function(e) data.frame(Message = "Error extracting coefficients")
      )
      write.csv(coef_df, coef_file, row.names = FALSE)
      files <- c(files, coef_file)

      # 4. Time to Event Curve (PDF)
      curve_file <- file.path(tmpdir, paste0(prefix, "_time_to_event_curve_", subset_id, ".pdf"))
      tryCatch(
        {
          surv_data <- s$data
          surv_data$group <- factor(surv_data$risk_group, levels = c("Low Risk", "High Risk"))
          
          # Check if we have at least 2 groups
          actual_groups <- unique(as.character(surv_data$group))
          
          if (length(actual_groups) < 2) {
            # Only 1 group - create error plot
            pdf(curve_file, width = 8, height = 6)
            plot.new()
            text(0.5, 0.5, "Cannot create time-to-event curves:\nonly 1 risk group present", cex = 1.2, col = "red")
            dev.off()
          } else {
            # We have 2 groups - proceed with plotting
            fit <- survival::survfit(Surv(time, status) ~ group, data = surv_data)
            
            # Build plot title
            plot_title <- if (!is.null(s$analysis_mode) && s$analysis_mode == "univariate" && !is.null(s$selected_predictor)) {
              # Get display name for predictor
              display_name <- s$selected_predictor
              if (!is.null(s$name_mapping) && s$selected_predictor %in% names(s$name_mapping)) {
                display_name <- s$name_mapping[[s$selected_predictor]]
              }
              sprintf("Time to Event Curves: %s (%s split)", display_name, tools::toTitleCase(s$split_method %||% "median"))
            } else {
              sprintf("Time to Event Curves by Predicted Risk (%s split)", tools::toTitleCase(s$split_method %||% "median"))
            }
            
            plot_args <- list(
              fit = fit,
              data = surv_data,
              pval = TRUE,
              pval.method = TRUE,
              pval.method.coord = c(0.05, 0.15),
              pval.coord = c(0.05, 0.2),
              pval.hjust = 0.5,
              conf.int = isTRUE(input$surv_show_ci),
              risk.table = TRUE,
              risk.table.height = 0.25,
              xlab = "Time to Event",
              ylab = "Event-Free Probability",
              title = plot_title,
              subtitle = "P-value: Log-rank test (dichotomized groups)",
              legend.title = "Risk Group",
              palette = c("#00BA38", "#F8766D"),
              ggtheme = theme_minimal(base_size = 14)
            )
            
            # Only add legend.labs if we have exactly 2 groups
            if (length(actual_groups) == 2) {
              plot_args$legend.labs <- c("Low Risk", "High Risk")
            }
            
            surv_plot <- do.call(ggsurvplot, plot_args)
            
            ggsave(curve_file,
              plot = surv_plot$plot,
              device = if (capabilities("cairo")) cairo_pdf else pdf,
              width = 8, height = 6, units = "in"
            )
          }
        },
        error = function(e) {
          pdf(curve_file, width = 8, height = 6)
          plot.new()
          text(0.5, 0.5, paste("Error creating time-to-event curve:", conditionMessage(e)), cex = 0.9, col = "red")
          dev.off()
        }
      )
      files <- c(files, curve_file)

      # Bundle into zip
      zip::zip(zipfile = file, files = files, mode = "cherry-pick")
    },
    contentType = "application/zip"
  )
}

shinyApp(ui, server)
