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
    library(rsample) # vfold_cv(), initial_split (if you want tidy resampling)
    library(boot) # bootstrapping CIs for AUC
    library(rlang) # tidy evaluation for !!!syms
  })
})

options(shiny.maxRequestSize = 100000 * 1024^2)

# ---- Helpers & validation ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

pct_clip <- function(x, p = c(0.01, 0.99)) {
  q <- quantile(x, probs = p, na.rm = TRUE)
  pmin(pmax(x, q[1]), q[2])
}

appendLog <- function(msg) {
  old <- rv$log()
  new <- c(old, paste0(format(Sys.time(), "%H:%M:%S"), " | ", msg))
  if (length(new) > 10) new <- tail(new, 10)  # keep last 10 lines
  rv$log(new)
}

align_metadata_abundance <- function(metadata, abundance) {
  # Extract patient_ID from abundance rownames ("patientID_runDate.fcs" is the expected format; the pattern "_[0-9]+\\-[A-Za-z]+\\-[0-9]+.*$" will always match the _runDate.fcs part)
  patient_ids <- gsub(pattern = "_[0-9]+\\-[A-Za-z]+\\-[0-9]+.*$", replacement = "", x = rownames(abundance))
  abund_df <- as.data.frame(abundance)
  abund_df$patient_ID <- patient_ids
  
  # Merge with metadata
  merged <- dplyr::left_join(metadata, abund_df, by = "patient_ID")
  return(merged)
}

clean_dummy_names <- function(nms) {
  gsub(pattern = 'leiden_cluster', replacement = 'leiden_cluster:', x = nms)
}

spearman_test <- function(df, freq_col = "freq", cont_var) {
  ok <- complete.cases(df[[freq_col]], df[[cont_var]])
  if (!any(ok)) return(data.frame(rho = NA_real_, p = NA_real_, n = 0))
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
        pickerInput(ns("split_by"), "Facet by", choices = NULL,
                    options = list(`none-selected-text` = "None")),
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
    message(sprintf("coords NULL? %s | expr NULL? %s | meta_cell NULL? %s",
                    is.null(coords), is.null(expr), is.null(meta_cell)))
    
    # One-time picker initialization in a reactive context
    initialized <- FALSE
    plot_cache_gg <- reactiveVal(NULL)
    facet_cols_saved <- reactiveVal(2)
    
    observeEvent(list(expr(), meta_cell(), clusters()), ignoreInit = FALSE, {
      if (initialized) return()
      expr_val     <- expr()
      meta_val     <- meta_cell()
      clusters_val <- clusters()
      
      req(!is.null(expr_val), !is.null(meta_val))
      
      # Add leiden_cluster column as factor if missing
      if (!"leiden_cluster" %in% colnames(meta_val) && !is.null(clusters_val$assignments)) {
        meta_val$leiden_cluster <- factor(clusters_val$assignments)
      }
      
      numeric_markers <- colnames(expr_val)
      meta_cols <- setdiff(colnames(meta_val), c(".cell"))
      
      # Ensure leiden_cluster is present and comes right after markers
      if (!"leiden_cluster" %in% meta_cols && "leiden_cluster" %in% colnames(meta_val)) {
        meta_cols <- c("leiden_cluster", setdiff(meta_cols, "leiden_cluster"))
      }
      
      # Sort metadata portion alphabetically (excluding leiden_cluster)
      meta_cols_sorted <- sort(setdiff(meta_cols, "leiden_cluster"))
      
      # Final order: markers → leiden_cluster → sorted metadata
      ordered_choices <- c(numeric_markers, "leiden_cluster", meta_cols_sorted)
      
      # Continuous and categorical choices from metadata
      cont_choices <- sort(meta_cols[sapply(meta_val[meta_cols], is.numeric)])
      factor_cols  <- meta_cols[sapply(meta_val[meta_cols], is.factor)]
      char_cols    <- meta_cols[sapply(meta_val[meta_cols], is.character)]
      categorical_choices <- sort(c(factor_cols, char_cols))
      
      # Default unit_var
      unit_default <- if ("PatientID" %in% meta_cols) "PatientID" else meta_cols[1]
      
      updatePickerInput(session, "color_by",
                        choices = ordered_choices,
                        selected = if (length(numeric_markers)) numeric_markers[1] else ordered_choices[1])
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
    
    plot_cache_base <- reactiveVal(NULL)  # base plot
    plot_cache      <- reactiveVal(NULL)  # final plot with overlays
    
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
            plot_cache_base(NULL); plot_cache(NULL)
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
            plot_cache_base(NULL); plot_cache(NULL)
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
        embed_name <- tolower(embedding_name)  # "umap" or "tsne"
        color_var  <- tolower(input$color_by %||% "color")
        
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
          
          pdf_width  <- 6 * ncol_facets
          pdf_height <- 6 * nrow_facets
        } else {
          # Single panel → fixed size
          pdf_width  <- 8
          pdf_height <- 8
        }
        
        ggsave(file, plot = gg, device = cairo_pdf,
               width = pdf_width, height = pdf_height, units = "in")
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
  tabPanel("Home",
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
  
  tabPanel("UMAP", EmbeddingUI("umap", title = "UMAP")),
  tabPanel("tSNE", EmbeddingUI("tsne", title = "tSNE")),
  
  tabPanel(
    "Heatmap",
    fluidRow(
      column(3,
             checkboxInput("cluster_rows", "Cluster rows", value = TRUE),
             checkboxInput("cluster_columns", "Cluster columns", value = TRUE),
             selectInput("heatmap_theme", "Heatmap color theme",
                         choices = c("viridis", "heat", "greyscale"),
                         selected = "viridis"),
             br(),
             downloadButton("export_heatmap_pdf", "Export heatmap as PDF")
      ),
      column(9, plotOutput("cluster_heatmap", height = "700px"))
    )
  ),
  tabPanel(
    "Testing",
    h4("Abundance testing"),
    fluidRow(
      column(
        3,
        conditionalPanel(
          condition = "output.hasClusterMap",
          pickerInput("test_entity", "Entity",
                      choices = c("Clusters", "Celltypes"),
                      selected = "Clusters")
        ),
        pickerInput("group_var", "Categorical metadata", choices = NULL,
                    options = list(`none-selected-text` = "None")),
        pickerInput("cont_var", "Continuous metadata", choices = NULL,
                    options = list(`none-selected-text` = "None")),
        radioButtons("test_type", "Test",
                     choices = c("Wilcoxon (2-group)",
                                 "Kruskal–Wallis (multi-group)",
                                 "Spearman (continuous)")),
        selectInput("p_adj_method", "P‑value adjustment method",
                    choices = c("BH", "bonferroni", "BY", "fdr"),
                    selected = "BH"),
        actionButton("run_test", "Run tests"),
        br(), br(),
        conditionalPanel(
          condition = "output.hasResults",
          downloadButton("export_results", "Export results as CSV")
        )
      ),
      column(9, tableOutput("test_table"))
    )
  ), 
  tabPanel(
    "Categorical",
    h4("Categorical plotting"),
    fluidRow(
      column(
        3,
        conditionalPanel(
          condition = "output.hasClusterMap",
          pickerInput("cat_entity", "Entity",
                      choices = c("Clusters", "Celltypes"),
                      selected = "Clusters")
        ),
        # If no cluster map, show a disabled picker fixed to Clusters
        conditionalPanel(
          condition = "!output.hasClusterMap",
          pickerInput("cat_entity", "Entity",
                      choices = c("Clusters"),
                      selected = "Clusters")
        ),
        pickerInput("cat_group_var", "Categorical metadata",
                    choices = NULL,
                    options = list(`none-selected-text` = "None")),
        radioButtons("cat_test_type", "Test",
                     choices = c("Wilcoxon (2-group)",
                                 "Kruskal–Wallis (multi-group)")),
        selectInput("cat_p_adj_method", "P‑value adjustment method",
                    choices = c("BH", "bonferroni", "BY", "fdr"),
                    selected = "BH"),
        checkboxInput("cat_use_adj_p", "Plot adjusted pvalues", value = TRUE),
        selectInput("cat_max_facets", "Facet columns",
                    choices = 2:6, selected = 4), 
        selectInput("cat_plot_type", "Plot type",
                    choices = c("Boxplot" = "box", "Violin" = "violin"), selected = "box"),
        radioButtons("cat_points", "Show data points", 
                     choices = c("Draw" = "draw", "Draw with jitter" = "jitter", "Do not draw" = "none"), selected = "draw"), 
        actionButton("cat_populate_colors", "Populate colors for selected group variable"), 
        uiOutput("cat_color_pickers_ui"),  # dynamic UI container for per-group color pickers
        br(), 
        actionButton("generate_cat_plots", "Generate plots"), 
        br(), br(),
        downloadButton("export_cat_pdf", "Export boxplots as PDF")
      ),
      column(
        9,
        column(9, plotOutput("categorical_plot"))
      )
    )
  ), 
  tabPanel(
    "Continuous",
    h4("Continuous metadata scatter plots"),
    fluidRow(
      column(
        3,
        conditionalPanel(
          condition = "output.hasClusterMap",
          pickerInput("cont_entity", "Entity",
                      choices = c("Clusters", "Celltypes"),
                      selected = "Clusters")
        ),
        conditionalPanel(
          condition = "!output.hasClusterMap",
          pickerInput("cont_entity", "Entity",
                      choices = c("Clusters"),
                      selected = "Clusters")
        ),
        pickerInput("cont_group_var", "Continuous metadata",
                    choices = NULL,
                    options = list(`none-selected-text` = "None")),
        selectInput("cont_p_adj_method", "P‑value adjustment method",
                    choices = c("BH", "bonferroni", "BY", "fdr"),
                    selected = "BH"),
        checkboxInput("cont_use_adj_p", "Plot adjusted pvalues", value = TRUE),
        checkboxInput("cont_transpose", "Transpose axes", value = FALSE), 
        selectInput("cont_max_facets", "Facet columns", choices = 2:6, selected = 4),
        actionButton("generate_cont_plots", "Generate plots"),
        br(), br(),
        downloadButton("export_cont_pdf", "Export scatter plots as PDF")
      ),
      column(
        9,
        column(9, plotOutput("continuous_plot"))
      )
    )
  ), 
  tabPanel(
    "Feature Selection",
    sidebarLayout(
      sidebarPanel(
        tabsetPanel(
          tabPanel("Controls",
                   pickerInput("fs_method", "Method", 
                               choices = c("Ridge Regression", "Elastic Net", "Random Forest (Boruta)")),
                   pickerInput("fs_outcome", "Outcome variable", choices = NULL),
                   pickerInput("fs_predictors", "Predictor(s)", choices = NULL, multiple = TRUE),
                   conditionalPanel(
                     condition = "input.fs_predictors.includes('leiden_cluster')",
                     pickerInput("fs_leiden_subset", "Select clusters to include",
                                 choices = NULL, multiple = TRUE)
                   ), 
                   conditionalPanel(
                     condition = "input.fs_method == 'Elastic Net'",
                     sliderInput("fs_alpha", "Alpha (0 = Ridge, 1 = Lasso)", 
                                 min = 0, max = 1, value = 0.5, step = 0.05)
                   ),
                   actionButton("run_fs", "Run Feature Selection"),
                   br(), br(),
                   conditionalPanel(
                     condition = "output.hasFSResults",
                     downloadButton("export_fs_results", "Export results as CSV")
                   )
          ),
          tabPanel("Console Log",
                   verbatimTextOutput("console_log", placeholder = TRUE)
          )
        )
      ),
      mainPanel(
        h4("Summary Plot"),
        plotOutput("fs_plot", height = "550px"),
        
        h4("Selected Features"),
        tableOutput("fs_results"),
        
        h4("Details"),
        verbatimTextOutput("fs_summary")
      )
    )
  ), 
  tabPanel("Logistic Modeling",
           sidebarLayout(
             sidebarPanel(
               pickerInput("lm_outcome", "Outcome variable", choices = NULL),
               pickerInput("lm_predictors", "Predictor(s)", choices = NULL, multiple = TRUE),
               conditionalPanel(
                 condition = "input.lm_predictors.includes('leiden_cluster')",
                 pickerInput("lm_leiden_subset", "Select clusters to include", choices = NULL, multiple = TRUE)
               ),
               radioButtons("lm_model_type", "Model type",
                            choices = c("Logistic Regression", "Elastic Net", "Random Forest")),
               conditionalPanel(
                 condition = "input.lm_model_type == 'Elastic Net'",
                 sliderInput("lm_alpha", "Elastic Net alpha (0 = Ridge, 1 = Lasso)",
                             min = 0, max = 1, value = 0.5, step = 0.05)
               ),
               radioButtons("lm_validation", "Validation strategy",
                            choices = c("Train/Test split", "k-fold CV", "Leave-One-Out")),
               conditionalPanel(
                 condition = "input.lm_validation == 'Train/Test split'",
                 sliderInput("lm_train_frac", "Train fraction",
                             min = 0.5, max = 0.95, value = 0.7, step = 0.05)
               ),
               conditionalPanel(
                 condition = "input.lm_validation == 'k-fold CV'",
                 numericInput("lm_k", "Number of folds", value = 5, min = 2, max = 20)
               ),
               actionButton("run_lm", "Run Model")
             ),
             mainPanel(
               h4("ROC Curves"),
               plotOutput("lm_roc_plot", height = "500px"),
               h4("Performance Summary"),
               tableOutput("lm_perf_table")
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
    data_ready = FALSE
  )
  cat_plot_cache <- reactiveVal(NULL)
  cont_plot_cache <- reactiveVal(NULL)
  rv$log <- reactiveVal(character())
  
  output$console_log <- renderText({
    paste(rv$log(), collapse = "\n")
  })
  
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
    session$sendCustomMessage("enableTabs", FALSE)  # lock tabs during load
    
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
    meta_cell[] <- lapply(meta_cell, function(col) {
      if (is.character(col) && all(grepl("^\\s*-?\\d*(\\.\\d+)?\\s*$", col[!is.na(col)]))) {
        as.numeric(col)
      } else {
        col
      }
    })
    
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
    
    # Add abundance matrix if present and well-formed
    if (!is.null(obj$leiden$abundance) && is.matrix(obj$leiden$abundance)) {
      clusters$abundance <- obj$leiden$abundance
      
      # Strongly recommended: require rownames(sources) for mapping to metadata
      if (is.null(rownames(clusters$abundance)) || any(!nzchar(rownames(clusters$abundance)))) {
        showNotification(
          "leiden$abundance has missing rownames; cannot map to metadata. Please set rownames to source strings.",
          type = "error",
          duration = NULL
        )
        message("Abundance matrix rownames are missing or empty; mapping to metadata will fail.")
      } else {
        message(sprintf("Abundance matrix loaded: %d sources × %d entities",
                        nrow(clusters$abundance), ncol(clusters$abundance)))
      }
    } else {
      showNotification("No leiden$abundance matrix found in upload; abundance testing will be disabled.", type = "warning")
      message("No obj$leiden$abundance; rv$clusters$abundance will remain NULL.")
    }
    
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
    
    rv$cluster_map <- cluster_map
    
    message("Upload complete: expr rows=", nrow(rv$expr),
            " meta_cell rows=", nrow(rv$meta_cell),
            " UMAP coords=", if (!is.null(rv$UMAP)) nrow(rv$UMAP$coords) else "NULL",
            " tSNE coords=", if (!is.null(rv$tSNE)) nrow(rv$tSNE$coords) else "NULL")
    
    showNotification("Data loaded and initialized.", type = "message")
    
    rv$data_ready <- TRUE
    session$sendCustomMessage("enableTabs", TRUE)
  })
  
  observeEvent(input$main_tab, {
    if (!isTRUE(rv$data_ready) && !identical(input$main_tab, "Home")) {
      updateNavbarPage(session, "main_tab", selected = "Home")
    }
  })
  
  # UI-facing flag for conditionalPanel (no nested reactive)
  output$hasClusterMap <- reactive({
    !is.null(rv$cluster_map) && all(c("cluster", "celltype") %in% names(rv$cluster_map))
  })
  outputOptions(output, "hasClusterMap", suspendWhenHidden = FALSE)
  
  observeEvent(rv$cluster_map, {
    cm <- rv$cluster_map
    if (!is.null(cm) && all(c("cluster","celltype") %in% names(cm))) {
      # nothing to do
    } else {
      updatePickerInput(session, "test_entity", selected = "Clusters")
    }
  }, ignoreInit = TRUE)
  
  observeEvent(rv$cluster_map, {
    cm <- rv$cluster_map
    if (!is.null(cm) && all(c("cluster","celltype") %in% names(cm))) {
      # nothing to do
    } else {
      updatePickerInput(session, "cat_entity", selected = "Clusters")
    }
  }, ignoreInit = TRUE)
  
  # Auto-detect categorical vs continuous metadata
  observeEvent(rv$meta_cell, {
    meta_cols <- colnames(rv$meta_cell)
    
    categorical_choices <- sort(meta_cols[sapply(rv$meta_cell, function(x)
      is.character(x) || is.factor(x)
    )])
    continuous_choices <- sort(meta_cols[sapply(rv$meta_cell, function(x)
      is.numeric(x) || is.integer(x)
    )])
    
    updatePickerInput(session, "group_var",
                      choices = c("", categorical_choices), selected = "")
    updatePickerInput(session, "cont_var",
                      choices = c("", continuous_choices), selected = "")
  }, ignoreInit = TRUE)
  
  observeEvent(rv$meta_cell, {
    meta_cols <- colnames(rv$meta_cell)
    continuous_choices <- sort(meta_cols[sapply(rv$meta_cell, is.numeric)])
    updatePickerInput(session, "cont_group_var",
                      choices = c("", continuous_choices),
                      selected = "")
  }, ignoreInit = TRUE)
  
  observeEvent(rv$meta_cell, {
    meta_cols <- sort(colnames(rv$meta_cell))
    updatePickerInput(session, "model_outcome", choices = meta_cols)
    updatePickerInput(session, "model_predictors", choices = meta_cols)
    updatePickerInput(session, "model_covariates", choices = meta_cols)
    updatePickerInput(session, "model_random", choices = meta_cols)
  }, ignoreInit = TRUE)
  
  # Launch embedding modules as soon as data is ready
  observeEvent(list(rv$UMAP, rv$data_ready), {
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
  }, ignoreInit = TRUE)
  
  observeEvent(list(rv$tSNE, rv$data_ready), {
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
  }, ignoreInit = TRUE)
  
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
    palette_choice <- switch(
      input$heatmap_theme,
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
      paste0("cluster_heatmap_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      # Get the underlying matrix from your reactive
      M <- rv$cluster_heat  # or heatmap_obj()$matrix if you store it there
      
      # Calculate dimensions using the same factors as fcs_plot_heatmap
      pdf_width  <- (ncol(M) * 0.33) + 3.25
      pdf_height <- (nrow(M) * 0.3)  + 2.25
      
      pdf(file, width = pdf_width, height = pdf_height)
      draw(heatmap_obj())
      dev.off()
    },
    contentType = "application/pdf"
  )
  
  output$hasResults <- reactive({
    df <- run_tests()
    !is.null(df) && nrow(df) > 0
  })
  outputOptions(output, "hasResults", suspendWhenHidden = FALSE)
  
  # Store results + adj_col name from the run
  run_tests <- eventReactive(input$run_test, {
    req(rv$meta_cell)
    
    # Capture all relevant inputs at run time
    test_type_run     <- input$test_type
    p_adj_method_run  <- input$p_adj_method
    group_var_run     <- input$group_var
    cont_var_run      <- input$cont_var
    test_entity_run   <- input$test_entity
    
    abund0 <- rv$clusters$abundance
    if (is.null(abund0)) {
      showNotification("No abundance matrix available. Ensure obj$leiden$abundance is present in the upload.",
                       type = "error", duration = NULL)
      return(list(df = NULL, adj_col = NULL))
    }
    if (is.null(rownames(abund0)) || any(!nzchar(rownames(abund0)))) {
      showNotification("Abundance matrix has no valid rownames; cannot map to metadata.",
                       type = "error", duration = NULL)
      return(list(df = NULL, adj_col = NULL))
    }
    
    abund <- abund0
    if (test_entity_run == "Celltypes" && !is.null(rv$cluster_map) &&
        all(c("cluster", "celltype") %in% names(rv$cluster_map))) {
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
    
    if (!nzchar(group_var_run) && !nzchar(cont_var_run)) {
      return(list(df = data.frame(entity = NA, test = NA, p = NA, n = NA), adj_col = NULL))
    }
    
    # Map abundance rows to patient_ID
    sources <- rownames(abund)
    ids <- unique(rv$meta_cell$patient_ID)
    ids_esc <- stringr::str_replace_all(ids, "([\\^$.|?*+()\\[\\]{}\\\\])", "\\\\\\1")
    ids_esc <- ids_esc[order(nchar(ids_esc), decreasing = TRUE)]
    pattern <- paste0("(", paste0(ids_esc, collapse = "|"), ")")
    pid <- stringr::str_extract(sources, pattern)
    
    abund_df <- as.data.frame(abund)
    abund_df$patient_ID <- pid
    meta_unique <- rv$meta_cell %>% dplyr::distinct(patient_ID, .keep_all = TRUE)
    abund_df <- abund_df %>% dplyr::left_join(meta_unique, by = "patient_ID")
    
    abund_long <- abund_df %>%
      tidyr::pivot_longer(cols = colnames(abund), names_to = "entity", values_to = "freq")
    
    res <- abund_long %>%
      dplyr::group_by(entity) %>%
      dplyr::group_modify(~ {
        if (test_type_run == "Wilcoxon (2-group)") {
          if (!nzchar(group_var_run)) return(data.frame(test = "wilcox", p = NA, n = nrow(.x)))
          g <- droplevels(factor(.x[[group_var_run]]))
          ok <- !is.na(g)
          g <- g[ok]; freq_ok <- .x$freq[ok]
          if (length(levels(g)) != 2) return(data.frame(test = "wilcox", p = NA, n = sum(ok)))
          summaries <- tapply(freq_ok, g, function(v) {
            med <- median(v, na.rm = TRUE)
            q25 <- quantile(v, 0.25, na.rm = TRUE)
            q75 <- quantile(v, 0.75, na.rm = TRUE)
            n_grp <- sum(!is.na(v))
            sprintf("%.2f (%.2f-%.2f, n=%d)", med, q25, q75, n_grp)
          })
          sum_df <- as.data.frame(as.list(summaries), stringsAsFactors = FALSE)
          names(sum_df) <- paste0(names(sum_df), "_IQR")
          wt <- suppressWarnings(wilcox.test(freq_ok ~ g))
          cbind(data.frame(test = "wilcox", n = sum(ok), p = wt$p.value), sum_df)
          
        } else if (test_type_run == "Kruskal–Wallis (multi-group)") {
          if (!nzchar(group_var_run)) return(data.frame(test = "kruskal", p = NA, n = nrow(.x)))
          g <- .x[[group_var_run]]
          ok <- !is.na(g)
          if (!any(ok)) return(data.frame(test = "kruskal", p = NA, n = 0))
          summaries <- tapply(.x$freq[ok], g[ok], function(v) {
            med <- median(v, na.rm = TRUE)
            q25 <- quantile(v, 0.25, na.rm = TRUE)
            q75 <- quantile(v, 0.75, na.rm = TRUE)
            n_grp <- sum(!is.na(v))
            sprintf("%.2f (%.2f-%.2f, n=%d)", med, q25, q75, n_grp)
          })
          sum_df <- as.data.frame(as.list(summaries), stringsAsFactors = FALSE)
          names(sum_df) <- paste0(names(sum_df), "_IQR")
          kw <- kruskal.test(.x$freq[ok] ~ as.factor(g[ok]))
          cbind(data.frame(test = "kruskal", n = sum(ok), p = kw$p.value), sum_df)
          
        } else { # Spearman
          if (!nzchar(cont_var_run)) return(data.frame(test = "spearman", rho = NA, p = NA, n = nrow(.x)))
          cont <- cont_var_run
          ct <- spearman_test(.x, cont_var = cont)
          cbind(data.frame(test = "spearman"), ct)
        }
      }) %>%
      dplyr::ungroup()
    
    adj_col <- NULL
    if (nrow(res) && "p" %in% colnames(res) && nzchar(p_adj_method_run)) {
      adj_col <- paste0(tolower(p_adj_method_run), "_padj")
      res[[adj_col]] <- p.adjust(res$p, method = p_adj_method_run)
    }
    
    # Add metadata column name being tested
    tested_var <- if (test_type_run == "Spearman (continuous)") {
      if (nzchar(cont_var_run)) cont_var_run else NA_character_
    } else if (test_type_run %in% c("Wilcoxon (2-group)", "Kruskal–Wallis (multi-group)")) {
      if (nzchar(group_var_run)) group_var_run else NA_character_
    } else {
      NA_character_
    }
    
    # Determine entity used in this run
    has_map <- !is.null(rv$cluster_map) && all(c("cluster","celltype") %in% names(rv$cluster_map))
    entity_used <- if (!has_map) "clusters" else tolower(test_entity_run)
    
    # Determine test short name
    test_map <- c(
      "Wilcoxon (2-group)" = "wilcoxon",
      "Kruskal–Wallis (multi-group)" = "kruskal_wallis",
      "Spearman (continuous)" = "spearman"
    )
    test_used <- test_map[[test_type_run]]
    
    # Determine metadata used
    metadata_used <- if (test_type_run == "Spearman (continuous)") {
      if (nzchar(cont_var_run)) cont_var_run else "none"
    } else if (test_type_run %in% c("Wilcoxon (2-group)", "Kruskal–Wallis (multi-group)")) {
      if (nzchar(group_var_run)) group_var_run else "none"
    } else {
      "none"
    }
    
    # Save for later use in downloadHandler
    rv$last_test_info <- list(
      entity = entity_used,
      test = test_used,
      metadata = metadata_used
    )
    
    res$metadata <- tested_var
    
    # Ensure column order: entity, metadata, test, n, p, padj, rho (if present), then all *_IQR
    iqr_cols <- grep("_IQR$", names(res), value = TRUE)
    base_cols <- c("entity", "metadata", "test", "n", "p")
    if (!is.null(adj_col) && adj_col %in% names(res)) base_cols <- c(base_cols, adj_col)
    if ("rho" %in% names(res)) base_cols <- c(base_cols, "rho")
    res <- res[, c(base_cols, iqr_cols), drop = FALSE]
    
    list(df = res, adj_col = adj_col)
  })
  
  output$test_table <- renderTable({
    run <- run_tests()
    df <- req(run$df)
    adj_col <- run$adj_col
    
    # Format p and adjusted p (only if numeric)
    num_cols <- intersect(c("p", adj_col), names(df))
    for (col in num_cols) {
      if (is.numeric(df[[col]])) {
        df[[col]] <- formatC(df[[col]], format = "f", digits = 3)
      }
    }
    
    # Format rho if present and numeric
    if ("rho" %in% names(df) && is.numeric(df$rho)) {
      df$rho <- formatC(df$rho, format = "f", digits = 2)
    }
    
    # Order by adjusted p if present
    if (!is.null(adj_col) && adj_col %in% names(df) && is.numeric(as.numeric(df[[adj_col]]))) {
      df <- df[order(as.numeric(df[[adj_col]]), na.last = TRUE), ]
    }
    
    df
  }, sanitize.text.function = function(x) x)
  
  output$hasResults <- reactive({
    run <- run_tests()
    !is.null(run$df) && nrow(run$df) > 0
  })
  outputOptions(output, "hasResults", suspendWhenHidden = FALSE)
  
  output$export_results <- downloadHandler(
    filename = function() {
      info <- rv$last_test_info
      if (is.null(info)) return("results.csv")
      fname <- paste(info$entity, info$test, info$metadata, sep = "_")
      fname <- gsub(" ", "_", fname)
      fname <- tolower(fname)
      paste0(fname, ".csv")
    },
    content = function(file) {
      run <- run_tests()
      df <- run$df
      req(!is.null(df), nrow(df) > 0)
      write.csv(df, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  # Update cat_group_var choices when metadata arrives
  observeEvent(rv$meta_cell, {
    meta_cols <- colnames(rv$meta_cell)
    categorical_choices <- sort(meta_cols[sapply(rv$meta_cell, function(x) is.character(x) || is.factor(x))])
    updatePickerInput(session, "cat_group_var",
                      choices = c("", categorical_choices),
                      selected = "")
  }, ignoreInit = TRUE)
  
  # If you also need group_var/cont_var (Testing tab) — keep the same pattern:
  observeEvent(rv$meta_cell, {
    meta_cols <- colnames(rv$meta_cell)
    categorical_choices <- sort(meta_cols[sapply(rv$meta_cell, function(x) is.character(x) || is.factor(x))])
    continuous_choices  <- sort(meta_cols[sapply(rv$meta_cell, function(x) is.numeric(x) || is.integer(x))])
    updatePickerInput(session, "group_var", choices = c("", categorical_choices), selected = "")
    updatePickerInput(session, "cont_var",  choices = c("", continuous_choices),  selected = "")
  }, ignoreInit = TRUE)
  
  # Continuous metadata picker for “Continuous” tab
  observeEvent(rv$meta_cell, {
    meta_cols <- colnames(rv$meta_cell)
    continuous_choices <- sort(meta_cols[sapply(rv$meta_cell, is.numeric)])
    updatePickerInput(session, "cont_group_var",
                      choices = c("", continuous_choices),
                      selected = "")
  }, ignoreInit = TRUE)
  
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
    req(rv$meta_cell)   # require metadata to be present
    
    # Preferred source: use the last generated cat_plot_data if available (ensures exact plotting groups)
    levels_vec <- NULL
    cp <- NULL
    # try to use the precomputed eventReactive if it exists and has data
    if (exists("cat_plot_data", mode = "function")) {
      # guard around possible NULL result if not run yet
      try({
        cp <- cat_plot_data()
      }, silent = TRUE)
    }
    if (!is.null(cp) && is.list(cp) && !is.null(cp$data) && nrow(cp$data) > 0) {
      # use group values from the prepared plotting data (guarantees match)
      group_col <- cp$group_var
      if (!is.null(group_col) && nzchar(group_col) && group_col %in% colnames(cp$data)) {
        levels_vec <- unique(as.character(cp$data[[group_col]]))
      }
    }
    
    # If not available, fall back to the metadata column selected by the user (cat_group_var)
    if (is.null(levels_vec) || length(levels_vec) == 0) {
      gv <- input$cat_group_var
      if (!is.null(gv) && nzchar(gv) && gv %in% colnames(rv$meta_cell)) {
        levels_vec <- sort(unique(as.character(rv$meta_cell[[gv]])))
      } else {
        # nothing to show
        levels_vec <- character(0)
      }
    }
    
    # Exclude NA from choices so no color is created for NA
    levels_vec <- levels_vec[!is.na(levels_vec)]
    
    if (length(levels_vec) == 0) {
      showNotification("No non-missing group levels found to populate colors for.", type = "warning")
      output$cat_color_pickers_ui <- renderUI(NULL)
      rv$cat_colors <- NULL
      return()
    }
    
    # Default palette sized to number of levels
    n_levels <- length(levels_vec)
    default_pal <- viridis::viridis(n_levels)
    names(default_pal) <- levels_vec
    
    # Prefer GUI picker from colourpicker package if available; fallback to textInput
    use_colourpicker <- requireNamespace("colourpicker", quietly = TRUE)
    if (!use_colourpicker) {
      showNotification("Package 'colourpicker' not installed; using text inputs for hex colors. Install with install.packages('colourpicker') for a GUI picker.", type = "message", duration = 8)
    }
    
    # Build UI inputs (one per non-missing group)
    ui_list <- lapply(seq_along(levels_vec), function(i) {
      lv <- levels_vec[i]
      input_id <- paste0("cat_color_", sanitize_id(lv))
      label_text <- paste0("Color for ", lv)
      if (use_colourpicker) {
        colourpicker::colourInput(
          inputId = input_id,
          label = label_text,
          value = default_pal[i],
          showColour = "both"
        )
      } else {
        textInput(
          inputId = input_id,
          label = label_text,
          value = default_pal[i],
          placeholder = "#RRGGBB"
        )
      }
    })
    
    output$cat_color_pickers_ui <- renderUI({
      tagList(
        tags$div(
          style = "max-height: 300px; overflow-y: auto; padding-right: 6px;",
          ui_list
        )
      )
    })
    
    # Initialize rv$cat_colors with the defaults (named vector keyed by group label)
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
    req(rv$meta_cell, rv$clusters$abundance)
    
    abund0 <- rv$clusters$abundance
    if (is.null(rownames(abund0)) || any(!nzchar(rownames(abund0)))) {
      showNotification("Abundance matrix has no valid rownames; cannot map to metadata.", type = "error")
      return(NULL)
    }
    
    # Aggregate to celltypes if needed
    abund <- abund0
    if (input$cat_entity == "Celltypes" && !is.null(rv$cluster_map)) {
      cm <- rv$cluster_map
      keep <- cm$cluster %in% colnames(abund)
      cm <- cm[keep, , drop = FALSE]
      if (!nrow(cm)) return(NULL)
      split_idx <- split(cm$cluster, cm$celltype)
      abund <- sapply(split_idx, function(cols) rowSums(abund0[, cols, drop = FALSE]))
      abund <- as.matrix(abund)
    }
    
    # Map abundance rows to patient_ID via robust escaped regex
    sources <- rownames(abund)
    ids <- unique(rv$meta_cell$patient_ID)
    ids_esc <- stringr::str_replace_all(ids, "([\\\\^$.|?*+()\\[\\]{}\\\\])", "\\\\\\1")
    ids_esc <- ids_esc[order(nchar(ids_esc), decreasing = TRUE)]
    pattern <- paste0("(", paste0(ids_esc, collapse = "|"), ")")
    pid <- stringr::str_extract(sources, pattern)
    
    abund_df <- as.data.frame(abund, check.names = FALSE, stringsAsFactors = FALSE)
    abund_df$patient_ID <- pid
    meta_unique <- rv$meta_cell %>% dplyr::distinct(patient_ID, .keep_all = TRUE)
    abund_df <- abund_df %>% dplyr::left_join(meta_unique, by = "patient_ID")
    
    # Long format and remove NA abundance rows immediately (so testing/plotting ignore NA freqs)
    abund_long <- abund_df %>%
      tidyr::pivot_longer(cols = colnames(abund), names_to = "entity", values_to = "freq") %>%
      dplyr::mutate(entity = gsub(pattern = "\\n", replacement = " ", x = entity))
    
    # Remove rows with NA freq up-front and also rows with NA grouping variable later
    abund_long <- abund_long %>% dplyr::filter(!is.na(freq))
    
    # Capture plot-type and point-mode at Generate time so later UI changes do not re-trigger calculations
    plot_type_selected <- input$cat_plot_type %||% "box"    # "box" or "violin"
    point_mode_selected <- input$cat_points %||% "draw"    # "draw", "jitter", "none"
    
    # Run test per entity; ensure we drop NA grouping values for the tested grouping variable
    test_type <- input$cat_test_type
    group_var <- input$cat_group_var
    
    # If group_var is blank, return NA p-values for each entity
    if (is.null(group_var) || !nzchar(group_var)) {
      res <- abund_long %>%
        dplyr::group_by(entity) %>%
        dplyr::summarise(p = NA_real_, .groups = "drop")
    } else {
      res <- abund_long %>%
        dplyr::group_by(entity) %>%
        dplyr::group_modify(~ {
          # drop rows where group_var is NA
          if (!(group_var %in% colnames(.x))) return(data.frame(p = NA_real_))
          g_raw <- .x[[group_var]]
          ok_rows <- !is.na(g_raw)
          if (!any(ok_rows)) return(data.frame(p = NA_real_))
          g <- droplevels(factor(g_raw[ok_rows]))
          freq_ok <- .x$freq[ok_rows]
          
          # guard against too few groups
          if (length(unique(g)) < 2) return(data.frame(p = NA_real_))
          
          if (test_type == "Wilcoxon (2-group)" && length(unique(g)) == 2) {
            wt <- suppressWarnings(wilcox.test(freq_ok ~ g))
            data.frame(p = wt$p.value)
          } else if (test_type == "Kruskal\u2013Wallis (multi-group)") {
            kw <- kruskal.test(freq_ok ~ as.factor(g))
            data.frame(p = kw$p.value)
          } else {
            data.frame(p = NA_real_)
          }
        }) %>%
        dplyr::ungroup()
    }
    
    # Adjust p-values if appropriate
    if (nrow(res) && nzchar(input$cat_p_adj_method) && "p" %in% names(res)) {
      res$padj <- p.adjust(res$p, method = input$cat_p_adj_method)
    }
    
    # Save info for export
    rv$last_cat_info <- list(
      entity = tolower(input$cat_entity %||% "clusters"),
      group = tolower(input$cat_group_var %||% "group"),
      test_raw = input$cat_test_type %||% "test"
    )
    
    # Return everything needed for plotting; include captured plot-type and point-mode
    list(
      data = abund_long,
      results = res,
      group_var = input$cat_group_var,
      use_adj_p = input$cat_use_adj_p,
      facet_cols = as.numeric(input$cat_max_facets),
      plot_type = plot_type_selected,
      point_mode = point_mode_selected
    )
  })
  
  output$categorical_plot <- renderPlot({
    cp <- cat_plot_data()
    req(cp)
    abund_long <- cp$data
    res <- cp$results
    group_var <- cp$group_var
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
    
    # Prepare p-value annotation dataframe
    p_df <- res %>%
      dplyr::mutate(
        p_to_show = if (isTRUE(use_adj_p) && "padj" %in% names(res)) padj else p,
        label = paste0("p = ", signif(p_to_show, 3)),
        x = length(unique(abund_long_plot[[group_var]])) / 2 + 0.5,
        y = tapply(abund_long_plot$freq, abund_long_plot$entity, max, na.rm = TRUE)[entity] * 1.05
      )
    
    # Determine group levels present in plotting data
    grp_levels <- unique(as.character(abund_long_plot[[group_var]]))
    if (length(grp_levels) == 0) {
      showNotification("No valid group levels found for plotting.", type = "error")
      return(invisible(NULL))
    }
    
    # Resolve manual colors non-reactively (isolate => changing pickers won't auto-trigger replot)
    colors_named <- isolate({
      if (!is.null(rv$cat_colors) && length(rv$cat_colors) > 0) {
        matched <- rv$cat_colors[names(rv$cat_colors) %in% grp_levels]
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
    })
    
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
    } else {
      # none -> no point layer
    }
    
    gg <- gg +
      ggplot2::facet_wrap(~entity, ncol = facet_cols, scales = "free_y") +
      ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      ggplot2::theme_bw(base_size = 18) +
      ggplot2::labs(x = group_var, y = "Frequency") +
      ggplot2::theme(
        strip.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 1.1, b = 1.1)),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1) # rotate facet x-axis labels 45 degrees
      )
    
    # Apply manual fill scale; points use black border so no color scale required for border
    gg <- gg + ggplot2::scale_fill_manual(values = colors_named)
    
    # Add p-value text annotations (filter to plotted entities)
    if (!is.null(p_df) && nrow(p_df) > 0) {
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
    
    # Cache for export
    cat_plot_cache(gg)
    
    gg
  },
  height = function() {
    gg <- cat_plot_cache()
    cp <- cat_plot_data()
    if (is.null(gg) || is.null(cp)) return(400)
    n_facets <- length(unique(gg$data$entity))
    ncol_facets <- cp$facet_cols
    if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
    nrow_facets <- ceiling(n_facets / ncol_facets)
    200 * nrow_facets
  },
  width = function() {
    gg <- cat_plot_cache()
    cp <- cat_plot_data()
    if (is.null(gg) || is.null(cp)) return(400)
    ncol_facets <- cp$facet_cols
    if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
    225 * ncol_facets
  })
  
  observeEvent(input$cat_group_var, {
    output$cat_color_pickers_ui <- renderUI(NULL)
    rv$cat_colors <- NULL
  })
  
  output$export_cat_pdf <- downloadHandler(
    filename = function() {
      info <- rv$last_cat_info
      if (is.null(info)) return("categorical_plots.pdf")
      # Normalise test name
      test <- tolower(info$test_raw)
      test <- gsub("\\s+", "_", test)           # spaces → underscores
      test <- gsub("_\\(.*\\)", "", test)       # remove "(...)" parts
      test <- gsub("kruskal_wallis", "kruskal-wallis", test)  # special case
      paste0("categorical_", info$entity, "_", info$group, "_", test, ".pdf")
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
      pdf_width  <- 4 * ncol_facets
      pdf_height <- 3 * nrow_facets
      ggsave(file, plot = gg, device = cairo_pdf,
             width = pdf_width, height = pdf_height, units = "in")
    },
    contentType = "application/pdf"
  )
  
  cont_plot_data <- eventReactive(input$generate_cont_plots, {
    req(rv$meta_cell, rv$clusters$abundance)
    
    abund0 <- rv$clusters$abundance
    if (is.null(rownames(abund0)) || any(!nzchar(rownames(abund0)))) {
      showNotification("Abundance matrix has no valid rownames; cannot map to metadata.", type = "error")
      return(NULL)
    }
    
    # Aggregate to celltypes if needed
    abund <- abund0
    if (input$cont_entity == "Celltypes" && !is.null(rv$cluster_map)) {
      cm <- rv$cluster_map
      keep <- cm$cluster %in% colnames(abund)
      cm <- cm[keep, , drop = FALSE]
      if (!nrow(cm)) return(NULL)
      split_idx <- split(cm$cluster, cm$celltype)
      abund <- sapply(split_idx, function(cols) rowSums(abund0[, cols, drop = FALSE]))
      abund <- as.matrix(abund)
    }
    
    # Map abundance rows to patient_ID
    sources <- rownames(abund)
    ids <- unique(rv$meta_cell$patient_ID)
    ids_esc <- stringr::str_replace_all(ids, "([\\\\^$.|?*+()\\[\\]{}\\\\])", "\\\\\\1")
    ids_esc <- ids_esc[order(nchar(ids_esc), decreasing = TRUE)]
    pattern <- paste0("(", paste0(ids_esc, collapse = "|"), ")")
    pid <- stringr::str_extract(sources, pattern)
    
    abund_df <- as.data.frame(abund, check.names = FALSE, stringsAsFactors = FALSE)
    abund_df$patient_ID <- pid
    meta_unique <- rv$meta_cell %>% dplyr::distinct(patient_ID, .keep_all = TRUE)
    abund_df <- abund_df %>% dplyr::left_join(meta_unique, by = "patient_ID")
    
    # Long format; remove NA freq immediately so tests and plotting ignore NA abundances
    abund_long <- abund_df %>%
      tidyr::pivot_longer(cols = colnames(abund), names_to = "entity", values_to = "freq") %>%
      dplyr::mutate(entity = gsub("\\n", " ", entity)) %>%
      dplyr::filter(!is.na(freq))
    
    # Capture transpose and cont_var at Generate time
    cont_var <- input$cont_group_var
    transpose_flag <- isTRUE(input$cont_transpose)
    test_padj_method <- input$cont_p_adj_method
    
    # Run Spearman per entity: drop rows where freq or cont_var is NA
    res <- abund_long %>%
      dplyr::group_by(entity) %>%
      dplyr::group_modify(~ {
        if (is.null(cont_var) || !nzchar(cont_var) || !(cont_var %in% colnames(.x))) {
          return(data.frame(rho = NA_real_, p = NA_real_, n = 0))
        }
        ok <- complete.cases(.x$freq, .x[[cont_var]])
        if (!any(ok)) return(data.frame(rho = NA_real_, p = NA_real_, n = 0))
        # Spearman is symmetric; orientation does not affect rho/p
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
    
    # Return everything needed for plotting; include the transpose flag and captured cont_var
    list(
      data = abund_long,
      results = res,
      cont_var = cont_var,
      use_adj_p = input$cont_use_adj_p,
      facet_cols = as.numeric(input$cont_max_facets),
      transpose = transpose_flag
    )
  })
  
  output$continuous_plot <- renderPlot({
    cp <- cont_plot_data()
    req(cp)
    
    abund_long <- cp$data
    res <- cp$results
    cont_var <- cp$cont_var
    use_adj_p <- cp$use_adj_p
    facet_cols <- cp$facet_cols
    transpose_flag <- cp$transpose %||% FALSE
    
    # Validate cont_var
    if (is.null(cont_var) || !nzchar(cont_var) || !(cont_var %in% colnames(abund_long))) {
      showNotification("No valid continuous metadata variable selected for continuous plotting.", type = "error")
      return(invisible(NULL))
    }
    
    # Remove rows with NA in either plotted axis (safe: cont_plot_data already removed NA freq but cont_var could have NA)
    if (!transpose_flag) {
      # Default: Abundance = x, continuous metadata = y
      plot_df <- abund_long %>% dplyr::filter(!is.na(freq) & !is.na(.data[[cont_var]]))
    } else {
      # Transposed: continuous metadata = x, Abundance = y
      plot_df <- abund_long %>% dplyr::filter(!is.na(freq) & !is.na(.data[[cont_var]]))
    }
    
    if (nrow(plot_df) == 0) {
      showNotification("No datapoints available for plotting after removing missing values.", type = "warning")
      return(invisible(NULL))
    }
    
    # Prepare p-value annotation dataframe (use same results computed at generate time)
    p_df <- res %>%
      dplyr::mutate(
        p_to_show = if (isTRUE(use_adj_p) && "padj" %in% names(res)) padj else p,
        label = paste0("p = ", signif(p_to_show, 3), "\n", "rho = ", signif(rho, 3))
      )
    
    # Choose aesthetics and trendline according to transpose_flag
    if (!transpose_flag) {
      # Abundance on x, continuous on y
      aes_pt <- ggplot2::aes(x = freq, y = .data[[cont_var]])
      smooth_aes <- ggplot2::aes(x = freq, y = .data[[cont_var]])
      x_lab <- "Abundance"
      y_lab <- cont_var
      # p-value label positions: put near top-right by entity
      p_df <- p_df %>%
        dplyr::mutate(
          x = tapply(plot_df$freq, plot_df$entity, function(v) mean(range(v, na.rm = TRUE)))[entity],
          y = tapply(plot_df[[cont_var]], plot_df$entity, max, na.rm = TRUE)[entity] * 1.05
        )
    } else {
      # Transposed: continuous on x, Abundance on y
      aes_pt <- ggplot2::aes(x = .data[[cont_var]], y = freq)
      smooth_aes <- ggplot2::aes(x = .data[[cont_var]], y = freq)
      x_lab <- cont_var
      y_lab <- "Abundance"
      # p-value label positions: put near top-right by entity
      p_df <- p_df %>%
        dplyr::mutate(
          x = tapply(plot_df[[cont_var]], plot_df$entity, function(v) mean(range(v, na.rm = TRUE)))[entity],
          y = tapply(plot_df$freq, plot_df$entity, max, na.rm = TRUE)[entity] * 1.05
        )
    }
    
    # Build the ggplot
    gg <- ggplot2::ggplot(plot_df, mapping = aes_pt) +
      ggplot2::geom_point(alpha = 0.75, pch = 21, color = 'black', fill = 'grey40', stroke = 0.1, size = 3) +
      ggplot2::geom_smooth(mapping = smooth_aes, method = "lm", se = FALSE, color = "red2") +
      ggplot2::facet_wrap(~entity, ncol = facet_cols, scales = "free") +
      ggplot2::scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      ggplot2::scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
      ggplot2::theme_bw(base_size = 18) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(x = x_lab, y = y_lab) +
      ggplot2::theme(strip.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 1.1, b = 1.1)))
    
    # Add p-value / rho annotations (filter to entities actually plotted)
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
    
    # Cache for export
    cont_plot_cache(gg)
    
    gg
  },
  height = function() {
    gg <- cont_plot_cache()
    cp <- cont_plot_data()
    if (is.null(gg) || is.null(cp)) return(400)
    n_facets <- length(unique(gg$data$entity))
    ncol_facets <- cp$facet_cols
    if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
    nrow_facets <- ceiling(n_facets / ncol_facets)
    200 * nrow_facets
  },
  width = function() {
    gg <- cont_plot_cache()
    cp <- cont_plot_data()
    if (is.null(gg) || is.null(cp)) return(400)
    ncol_facets <- cp$facet_cols
    if (is.na(ncol_facets) || ncol_facets < 1) ncol_facets <- 1
    225 * ncol_facets
  })
  
  output$export_cont_pdf <- downloadHandler(
    filename = function() {
      info <- rv$last_cont_info
      if (is.null(info)) return("continuous_plots.pdf")
      paste0("continuous_", info$entity, "_", info$group, "_", info$test_raw, ".pdf")
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
      pdf_width  <- 4 * ncol_facets
      pdf_height <- 3 * nrow_facets
      ggsave(file, plot = gg, device = cairo_pdf,
             width = pdf_width, height = pdf_height, units = "in")
    },
    contentType = "application/pdf"
  )
  
  # Populate Feature Selection dropdowns from same metadata source as other tabs
  observeEvent(rv$meta_cell, {
    meta_cols <- sort(colnames(rv$meta_cell))
    updatePickerInput(session, "fs_outcome", choices = meta_cols)
    updatePickerInput(session, "fs_predictors", choices = meta_cols)
  }, ignoreInit = TRUE)
  
  observeEvent(rv$clusters$abundance, {
    cluster_names <- colnames(rv$clusters$abundance)
    updatePickerInput(session, "fs_leiden_subset", choices = cluster_names)
  })
  
  run_fs <- eventReactive(input$run_fs, {
    req(rv$meta_cell, rv$clusters$abundance, input$fs_outcome, input$fs_predictors)
    
    # Align metadata and abundance
    merged <- align_metadata_abundance(rv$meta_cell, rv$clusters$abundance)
    
    # Drop NA outcomes
    merged <- merged[!is.na(merged[[input$fs_outcome]]), ]
    
    # Outcome
    y_raw <- merged[[input$fs_outcome]]
    
    # Predictors
    pred_meta <- intersect(input$fs_predictors, colnames(merged))
    if (length(pred_meta) == 0) {
      showNotification("Select at least one predictor.", type = "error")
      return(NULL)
    }
    X_raw <- merged[, pred_meta, drop = FALSE]
    
    # Sample counts
    n_before <- length(y_raw)
    complete_rows <- stats::complete.cases(data.frame(X_raw, .y = y_raw))
    X_raw <- X_raw[complete_rows, , drop = FALSE]
    y_raw <- y_raw[complete_rows]
    n_after <- length(y_raw)
    n_dropped <- n_before - n_after
    if (n_after < 3) {
      showNotification("Too few samples after filtering for feature selection.", type = "error")
      return(NULL)
    }
    
    # Outcome coercion
    y <- y_raw
    if (is.character(y)) y <- factor(y)
    
    # User controls
    method <- input$fs_method %||% "Elastic Net"
    tolerance <- 1e-4
    seed_val <- input$fs_seed %||% 123
    reps <- input$fs_reps %||% 1
    boruta_maxruns <- input$fs_maxruns %||% 500
    
    set.seed(seed_val)
    
    # Helper: clean dummy names
    clean_dummy_names <- function(nms) {
      gsub(pattern = "leiden_cluster", replacement = "leiden_cluster:", x = nms)
    }
    
    if (method == "Random Forest (Boruta)") {
      all_imp <- list()
      bor_models <- list()
      for (r in seq_len(reps)) {
        set.seed(seed_val + r - 1)
        X_boruta <- model.matrix(~ . - 1, data = X_raw)
        colnames(X_boruta) <- clean_dummy_names(colnames(X_boruta))
        X_boruta <- as.data.frame(X_boruta)
        
        if (is.factor(y)) {
          bor <- Boruta::Boruta(x = X_boruta, y = y, doTrace = 0, maxRuns = boruta_maxruns)
        } else {
          bor <- Boruta::Boruta(x = X_boruta, y = as.numeric(y), doTrace = 0, maxRuns = boruta_maxruns)
        }
        imp <- Boruta::attStats(bor)
        imp$meanImp[!is.finite(imp$meanImp)] <- 0
        all_imp[[r]] <- imp
        bor_models[[r]] <- bor
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
        summary = list(
          n_before = n_before,
          n_after = n_after,
          n_dropped = n_dropped,
          model = bor_models[[1]]  # return first Boruta model for inspection
        )
      ))
    }
    
    # glmnet-based methods
    family <- if (is.factor(y)) {
      if (nlevels(y) == 2) "binomial" else "multinomial"
    } else {
      "gaussian"
    }
    
    Xmat <- model.matrix(~ . - 1, data = X_raw)
    colnames(Xmat) <- clean_dummy_names(colnames(Xmat))
    storage.mode(Xmat) <- "double"
    
    alpha_val <- if (method == "Ridge Regression") 0 else (input$fs_alpha %||% 0.5)
    nfolds_val <- input$fs_nfolds %||% 5
    if (nfolds_val > n_after) nfolds_val <- max(3, floor(n_after / 2))
    
    coef_list <- list()
    cv_models <- list()
    for (r in seq_len(reps)) {
      set.seed(seed_val + r - 1)
      cvfit <- glmnet::cv.glmnet(
        x = Xmat, y = y, family = family,
        alpha = alpha_val, nfolds = nfolds_val
      )
      coef_obj <- coef(cvfit, s = "lambda.min")
      if (family == "multinomial") {
        coef_df <- do.call(rbind, lapply(names(coef_obj), function(cls) {
          cm <- as.matrix(coef_obj[[cls]])
          data.frame(Feature = rownames(cm), Class = cls, Coef = as.numeric(cm[, 1]), stringsAsFactors = FALSE)
        }))
        coef_df <- coef_df[coef_df$Feature != "(Intercept)", , drop = FALSE]
        coef_list[[r]] <- coef_df
      } else {
        cm <- as.matrix(coef_obj)
        coef_df <- data.frame(Feature = rownames(cm), Coef = as.numeric(cm[, 1]), stringsAsFactors = FALSE)
        coef_df <- coef_df[coef_df$Feature != "(Intercept)", , drop = FALSE]
        coef_list[[r]] <- coef_df
      }
      cv_models[[r]] <- cvfit
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
      summary = list(
        n_before = n_before,
        n_after = n_after,
        n_dropped = n_dropped,
        model = cv_models[[1]]  # return first cv.glmnet model for inspection
      )
    ))
  })
  
  observeEvent(input$run_fs, {
    rv$log(character())  # clear log
    appendLog("Starting feature selection...")
    
    withCallingHandlers(
      {
        res <- run_fs()
        rv$fs_result <- res
      },
      message = function(m) {
        appendLog(conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        appendLog(paste("Warning:", conditionMessage(w)))
        invokeRestart("muffleWarning")
      },
      error = function(e) {
        appendLog(paste("Error:", conditionMessage(e)))
        stop(e)
      }
    )
    
    appendLog("Finished feature selection.")
  })
  
  output$fs_plot <- renderPlot({
    # res <- run_fs(); req(res)
    res <- rv$fs_result; req(res)
    tol <- res$tolerance
    
    if (identical(res$method, "Boruta")) {
      # --- Boruta importance plot ---
      df <- res$results
      df <- df[abs(df$ImportanceMean) >= tol, , drop = FALSE]
      df <- df[order(df$ImportanceMean, decreasing = TRUE), ]
      top_n <- head(df, 30)
      
      ggplot2::ggplot(top_n,
                      ggplot2::aes(x = reorder(Feature, ImportanceMean),
                                   y = ImportanceMean,
                                   fill = Decision)) +
        ggplot2::geom_col(color = 'black', lwd = 0.4) +
        ggplot2::coord_flip() +
        ggplot2::labs(title = "Boruta feature importance",
                      x = "Feature", y = "Mean importance") +
        scale_fill_manual(values = c('Confirmed' = 'green4', 'Rejected' = 'red4', 'Tentative' = 'grey40')) + 
        ggplot2::theme_bw(base_size = 14)
      
    } else {
      # --- Elastic Net / Ridge ---
      df <- res$results
      df <- df[abs(df$Coef) >= tol, , drop = FALSE]
      
      if ("Class" %in% names(df)) {
        # Multinomial: facet by Class
        top_n <- df %>%
          dplyr::group_by(Class) %>%
          dplyr::arrange(dplyr::desc(Coef), .by_group = TRUE) %>%
          dplyr::slice_head(n = 20)
        
        ggplot2::ggplot(top_n,
                        ggplot2::aes(x = reorder(Feature, Coef),
                                     y = Coef,
                                     fill = Coef)) +
          ggplot2::geom_col(color = 'black', lwd = 0.4) +
          ggplot2::facet_wrap(~Class, scales = "free_y") +
          ggplot2::coord_flip() +
          ggplot2::labs(title = paste(res$method, "coefficients (lambda.min)"),
                        x = "Feature", y = "Coefficient") + 
          scale_fill_gradient(low = 'blue3', high = 'red3') + 
          ggplot2::theme_bw(base_size = 14)
        
      } else {
        # Binary/continuous outcome
        df <- df[order(df$Coef, decreasing = TRUE), ]
        top_n <- head(df, 30)
        
        ggplot2::ggplot(top_n,
                        ggplot2::aes(x = reorder(Feature, Coef),
                                     y = Coef,
                                     fill = Coef)) +
          ggplot2::geom_col(color = 'black', lwd = 0.4) +
          ggplot2::coord_flip() +
          ggplot2::labs(title = paste(res$method, "coefficients (lambda.min)"),
                        x = "Feature", y = "Coefficient") + 
          scale_fill_gradient(low = 'blue3', high = 'red3') + 
          ggplot2::theme_bw(base_size = 14)
      }
    }
  })
  
  output$fs_results <- renderTable({
    # res <- run_fs(); req(res)
    res <- rv$fs_result; req(res)
    df <- res$results
    tol <- res$tolerance
    
    if (identical(res$method, "Boruta")) {
      df <- df[abs(df$ImportanceMean) >= tol, , drop = FALSE]
      df <- df[order(df$ImportanceMean, decreasing = TRUE), ]
    } else {
      if ("Coef" %in% names(df)) {
        df <- df[abs(df$Coef) >= tol, , drop = FALSE]
        df <- df[order(df$Coef, decreasing = TRUE), ]
      }
    }
    df
  }, sanitize.text.function = function(x) x)
  
  # output$fs_summary <- renderPrint({
  #   res <- run_fs(); req(res)
  #   det <- res$details %||% list(samples_before = NA, samples_after = NA, samples_dropped = NA)
  #   tol <- res$tolerance
  #   
  #   cat("Samples before filtering:", det$samples_before, "\n")
  #   cat("Samples after filtering:", det$samples_after, "\n")
  #   cat("Samples dropped:", det$samples_dropped, "\n\n")
  #   
  #   cat("Selected features (top):\n")
  #   if (identical(res$method, "Boruta")) {
  #     df <- res$results
  #     df <- df[abs(df$ImportanceMean) >= tol, , drop = FALSE]
  #     df <- df[order(df$ImportanceMean, decreasing = TRUE), ]
  #     print(utils::head(df$Feature, 20))
  #   } else {
  #     df <- res$results
  #     if ("Coef" %in% names(df)) {
  #       df <- df[abs(df$Coef) >= tol, , drop = FALSE]
  #       df <- df[order(df$Coef, decreasing = TRUE), ]
  #       print(utils::head(df$Feature, 20))
  #     }
  #   }
  # })
  
  output$fs_summary <- renderPrint({
    res <- rv$fs_result; req(res)
    sumry <- res$summary
    tol <- res$tolerance
    
    cat("Samples before filtering:", sumry$n_before, "\n")
    cat("Samples after filtering:", sumry$n_after, "\n")
    cat("Samples dropped:", sumry$n_dropped, "\n\n")
    
    cat("Selected features (top):\n")
    if (identical(res$method, "Boruta")) {
      df <- res$results
      df <- df[abs(df$ImportanceMean) >= tol, , drop = FALSE]
      df <- df[order(df$ImportanceMean, decreasing = TRUE), ]
      print(utils::head(df$Feature, 20))
    } else {
      df <- res$results
      if ("Coef" %in% names(df)) {
        df <- df[abs(df$Coef) >= tol, , drop = FALSE]
        df <- df[order(df$Coef, decreasing = TRUE), ]
        print(utils::head(df$Feature, 20))
      }
    }
    
    cat("\nModel object:\n")
    print(sumry$model)
  })
  
  output$hasFSResults <- reactive({
    # !is.null(run_fs())
    !is.null(rv$fs_result)
  })
  outputOptions(output, "hasFSResults", suspendWhenHidden = FALSE)
  
  output$export_fs_results <- downloadHandler(
    filename = function() {
      info <- rv$last_fs_info
      if (is.null(info)) {
        return(paste0("feature_selection_", Sys.Date(), ".csv"))
      }
      paste0(info$method, "_feature_selection_with_outcome_", info$outcome, ".csv")
    },
    content = function(file) {
      # res <- run_fs(); req(res)
      res <- rv$fs_result; req(res)
      # res$results already has Method, Outcome, Alpha as last columns
      utils::write.csv(res$results, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  observeEvent(rv$meta_cell, {
    req(rv$meta_cell)
    
    meta_cols <- colnames(rv$meta_cell)
    
    # Outcomes = categorical metadata only (factor/character)
    categorical_choices <- sort(meta_cols[sapply(rv$meta_cell, function(x) is.factor(x) || is.character(x))])
    
    # Predictors = metadata only (no expression markers)
    predictor_choices <- sort(meta_cols)
    
    updatePickerInput(session, "lm_outcome",
                      choices = categorical_choices,
                      selected = NULL)
    
    updatePickerInput(session, "lm_predictors",
                      choices = predictor_choices,
                      selected = NULL)
    
    # Populate leiden_cluster subset choices if present; default to none selected
    if ("leiden_cluster" %in% meta_cols) {
      levs <- sort(unique(as.character(rv$meta_cell$leiden_cluster)))
      updatePickerInput(session, "lm_leiden_subset",
                        choices = levs,
                        selected = character(0))
    }
  }, ignoreInit = TRUE)
  
  lm_results <- eventReactive(input$run_lm, {
    req(rv$meta_cell, input$lm_outcome, input$lm_predictors)
    
    # --- Patient-level metadata ---
    meta_patient <- rv$meta_cell %>%
      dplyr::distinct(patient_ID, .keep_all = TRUE)
    
    if (!("patient_ID" %in% colnames(meta_patient))) {
      showNotification("No patient_ID column available in metadata.", type = "error")
      return(NULL)
    }
    
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
    meta_cols <- colnames(meta_patient)
    pred_meta <- intersect(input$lm_predictors, meta_cols)
    if (length(pred_meta) == 0) {
      showNotification("Select at least one predictor.", type = "error")
      return(NULL)
    }
    X <- meta_patient[, pred_meta, drop = FALSE]
    
    # Optional cluster filter
    if ("leiden_cluster" %in% input$lm_predictors &&
        "leiden_cluster" %in% colnames(meta_patient) &&
        !is.null(input$lm_leiden_subset) && length(input$lm_leiden_subset) > 0) {
      keep <- as.character(meta_patient$leiden_cluster) %in% input$lm_leiden_subset
      X <- X[keep, , drop = FALSE]
      outcome <- outcome[keep]
      meta_patient <- meta_patient[keep, , drop = FALSE]
    }
    
    # NA filtering
    model_df <- cbind(X, .outcome = outcome, patient_ID = meta_patient$patient_ID)
    na_rows <- which(!complete.cases(model_df[, setdiff(names(model_df), "patient_ID")]))
    if (length(na_rows) > 0) {
      dropped_ids <- model_df$patient_ID[na_rows]
      showNotification(
        paste("Dropped", length(na_rows), "patients due to missing values. Affected patient_IDs:",
              paste(unique(dropped_ids), collapse = ", ")),
        type = "warning", duration = 10
      )
      model_df <- model_df[-na_rows, , drop = FALSE]
    }
    
    X <- model_df[, setdiff(names(model_df), c(".outcome", "patient_ID")), drop = FALSE]
    outcome <- model_df$.outcome
    
    if (length(outcome) < 10) {
      showNotification("Too few samples after filtering. Model not run.", type = "error")
      return(NULL)
    }
    
    # Null model baseline
    majority_class <- names(which.max(table(outcome)))
    null_acc <- mean(outcome == majority_class)
    
    # Choose model method
    method <- switch(input$lm_model_type,
                     "Logistic Regression" = "multinom",
                     "Elastic Net"         = "glmnet",
                     "Random Forest"       = "rf")
    
    validation <- input$lm_validation
    
    # --- Train/Test split ---
    if (validation == "Train/Test split") {
      set.seed(123)
      train_frac <- input$lm_train_frac %||% 0.7
      idx <- caret::createDataPartition(outcome, p = train_frac, list = FALSE)
      trainX <- X[idx, , drop = FALSE]
      testX  <- X[-idx, , drop = FALSE]
      trainY <- outcome[idx]
      testY  <- outcome[-idx]
      
      train_counts <- table(trainY)
      test_counts  <- table(testY)
      
      split_info <- list(
        n_train = length(trainY),
        n_test  = length(testY),
        train_counts = train_counts,
        test_counts  = test_counts
      )
      
      if (method == "glmnet") {
        trainMat <- model.matrix(~ . - 1, data = trainX)
        testMat  <- model.matrix(~ . - 1, data = testX)
        storage.mode(trainMat) <- "double"
        storage.mode(testMat)  <- "double"
        
        model <- caret::train(
          x = trainMat, y = trainY,
          method = "glmnet",
          trControl = caret::trainControl(classProbs = TRUE),
          tuneGrid = expand.grid(alpha = input$lm_alpha,
                                 lambda = 10^seq(-3, 1, length = 20))
        )
        probs <- predict(model, newdata = testMat, type = "prob")
        pred_class <- predict(model, newdata = testMat, type = "raw")
      } else {
        model <- caret::train(
          x = trainX, y = trainY,
          method = method,
          trControl = caret::trainControl(classProbs = TRUE)
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
      preds_tbl$obs <- factor(preds_tbl$obs, levels = levels(outcome))
      
      # Binary probability completion
      classes <- levels(outcome)
      have_cols <- grep("^\\.pred_", names(preds_tbl), value = TRUE)
      if (length(classes) == 2 && length(have_cols) == 1) {
        missing_cls <- setdiff(classes, gsub("^\\.pred_", "", have_cols))
        preds_tbl[[paste0(".pred_", missing_cls)]] <- 1 - preds_tbl[[have_cols[1]]]
      }
      
      # --- Console logging ---
      message("=== Logistic Model (Train/Test) complete ===")
      message("Outcome levels: ", paste(levels(outcome), collapse = ", "))
      message("Prediction columns: ", paste(grep("^\\.pred_", names(preds_tbl), value = TRUE), collapse = ", "))
      utils::str(utils::head(preds_tbl, 5))
      
      return(list(model = model, preds = preds_tbl, null_acc = null_acc,
                  split_info = split_info))
    }
    
    # --- CV or LOO ---
    ctrl <- caret::trainControl(
      method = switch(validation,
                      "k-fold CV"     = "cv",
                      "Leave-One-Out" = "LOOCV"
      ),
      number = if (validation == "k-fold CV") input$lm_k else 1,
      classProbs = TRUE,
      savePredictions = "final"
    )
    
    if (method == "glmnet") {
      Xmat <- model.matrix(~ . - 1, data = X)
      storage.mode(Xmat) <- "double"
      model <- caret::train(
        x = Xmat, y = outcome,
        method = "glmnet",
        trControl = ctrl,
        tuneGrid = expand.grid(alpha = input$lm_alpha,
                               lambda = 10^seq(-3, 1, length = 20))
      )
    } else {
      model <- caret::train(
        x = X, y = outcome,
        method = method,
        trControl = ctrl
      )
    }
    
    preds <- model$pred
    if (is.null(preds) || !nrow(preds)) {
      showNotification("No predictions available from resampling. Try Train/Test split.", type = "error")
      return(NULL)
    }
    
    level_names <- levels(outcome)
    prob_cols <- intersect(level_names, names(preds))
    rename_map <- setNames(prob_cols, paste0(".pred_", prob_cols))
    preds <- dplyr::rename(preds, !!!rename_map)
    preds$obs <- factor(preds$obs, levels = levels(outcome))
    
    classes <- levels(outcome)
    have_cols <- grep("^\\.pred_", names(preds), value = TRUE)
    if (length(classes) == 2 && length(have_cols) == 1) {
      missing_cls <- setdiff(classes, gsub("^\\.pred_", "", have_cols))
      preds[[paste0(".pred_", missing_cls)]] <- 1 - preds[[have_cols[1]]]
    }
    
    # --- Console logging ---
    message("=== Logistic Model (CV/LOO) complete ===")
    message("Outcome levels: ", paste(levels(outcome), collapse = ", "))
    message("Prediction columns: ", paste(grep("^\\.pred_", names(preds), value = TRUE), collapse = ", "))
    utils::str(utils::head(preds, 5))
    
    return(list(model = model, preds = preds, null_acc = null_acc))
  })
  
  output$lm_roc_plot <- renderPlot({
    res <- lm_results(); req(res)
    preds <- res$preds; req(preds)
    
    n_classes <- nlevels(preds$obs)
    
    if (n_classes == 2) {
      # Binary ROC with pROC
      classes <- levels(preds$obs)
      positive <- classes[2]
      roc_obj <- pROC::roc(response = preds$obs,
                           predictor = preds[[paste0(".pred_", positive)]],
                           levels = rev(classes))
      plot(roc_obj, main = "Binary ROC Curve", col = "blue", lwd = 2, print.auc = TRUE)
      abline(a = 0, b = 1, lty = 2, col = "red")
    } else {
      # Multiclass ROC with yardstick (one-vs-rest curves per class)
      long_preds <- preds %>%
        dplyr::select(obs, starts_with(".pred_"))
      
      roc_curves <- yardstick::roc_curve(long_preds, truth = obs, dplyr::starts_with(".pred_"))
      
      ggplot2::ggplot(roc_curves, ggplot2::aes(x = 1 - specificity, y = sensitivity, color = .level)) +
        ggplot2::geom_path(size = 1) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
        ggplot2::labs(title = "Multiclass ROC Curves (yardstick)", color = "Class") +
        ggplot2::theme_bw(base_size = 14)
    }
  })
  
  output$lm_perf_table <- renderTable({
    res <- lm_results(); req(res)
    preds <- res$preds; req(preds)
    
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
    
    binom_p <- tryCatch({
      binom.test(successes, n, p = res$null_acc, alternative = "greater")$p.value
    }, error = function(e) NA)
    
    majority_safe <- names(which.max(table(obs)))
    majority_orig <- rv$lm_label_map[[majority_safe]]
    
    perf_table <- data.frame(
      Metric = c("Null Accuracy (most frequent class)", "Model Accuracy", "Binomial p-value vs Null"),
      Value  = c(paste0(round(res$null_acc, 3), " (", majority_orig, ")"),
                 round(model_acc, 3),
                 signif(binom_p, 3)),
      stringsAsFactors = FALSE
    )
    
    if (n_classes == 2) {
      # Binary AUC with pROC
      classes <- levels(preds$obs)
      positive <- classes[2]
      roc_obj <- pROC::roc(response = preds$obs,
                           predictor = preds[[paste0(".pred_", positive)]],
                           levels = rev(classes))
      auc_val <- as.numeric(pROC::auc(roc_obj))
      perf_table <- rbind(perf_table,
                          data.frame(Metric = "AUC (pROC)", Value = round(auc_val, 3), stringsAsFactors = FALSE))
    } else {
      # Multiclass AUCs with yardstick
      long_preds <- preds %>%
        dplyr::select(obs, starts_with(".pred_"))
      
      auc_macro <- yardstick::roc_auc(long_preds, truth = obs,
                                      dplyr::starts_with(".pred_"), estimator = "macro")
      auc_macro_wt <- yardstick::roc_auc(long_preds, truth = obs,
                                         dplyr::starts_with(".pred_"), estimator = "macro_weighted")
      auc_ht <- yardstick::roc_auc(long_preds, truth = obs,
                                   dplyr::starts_with(".pred_"), estimator = "hand_till")
      
      perf_table <- rbind(
        perf_table,
        data.frame(Metric = "Macro AUC (yardstick)", Value = round(auc_macro$.estimate, 3), stringsAsFactors = FALSE),
        data.frame(Metric = "Macro-weighted AUC (yardstick)", Value = round(auc_macro_wt$.estimate, 3), stringsAsFactors = FALSE),
        data.frame(Metric = "Hand-Till AUC (yardstick)", Value = round(auc_ht$.estimate, 3), stringsAsFactors = FALSE)
      )
    }
    
    # Train/Test split info
    if (!is.null(res$split_info)) {
      train_summary <- paste(names(res$split_info$train_counts),
                             res$split_info$train_counts, collapse = "; ")
      test_summary  <- paste(names(res$split_info$test_counts),
                             res$split_info$test_counts,  collapse = "; ")
      extra_rows <- data.frame(
        Metric = c("Train size", "Test size", "Train breakdown", "Test breakdown"),
        Value  = c(res$split_info$n_train, res$split_info$n_test, train_summary, test_summary),
        stringsAsFactors = FALSE
      )
      perf_table <- rbind(perf_table, extra_rows)
    }
    
    perf_table
  })
  
}

shinyApp(ui, server)
