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
                            active_tab, rv) {
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
      factor_cols <- meta_cols[sapply(meta_val[meta_cols], is.factor)]
      char_cols   <- meta_cols[sapply(meta_val[meta_cols], is.character)]
      categorical_choices <- c(factor_cols, char_cols)
      updatePickerInput(session, "split_by", choices = c("", categorical_choices), selected = "")
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
        h4("Upload inputs"),
        fileInput("rdata_upload", "Upload .RData (FCSimple analysis object)", accept = ".RData"),
        numericInput("max_cells_upload", "Max cells to read in", value = 300000, min = 1000, step = 1000),
        helpText("If the uploaded dataset has more cells than this number, it will be randomly downsampled at load time. This is done to speed up UMAP and tSNE facet plotting."), 
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
    "Clusters",
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
  rv$data_ready <- reactiveVal(FALSE)
  cat_plot_cache <- reactiveVal(NULL)
  cont_plot_cache <- reactiveVal(NULL)
  
  # Disable tabs at startup
  observe({
    if (!rv$data_ready()) {
      session$sendCustomMessage("enableTabs", FALSE)
    }
  })
  
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
    
    rv$data_ready(TRUE)
    session$sendCustomMessage("enableTabs", TRUE) # unlock tabs when ready
  })
  
  observeEvent(input$main_tab, {
    if (!rv$data_ready() && !identical(input$main_tab, "Home")) {
      updateNavbarPage(session, "main_tab", selected = "Home")
    }
  })
  
  # Server-side reactive flag
  hasClusterMap <- reactive({
    !is.null(rv$cluster_map) &&
      all(c("cluster", "celltype") %in% names(rv$cluster_map))
  })
  
  # UI-facing flag for conditionalPanel
  output$hasClusterMap <- reactive({
    hasClusterMap()
  })
  outputOptions(output, "hasClusterMap", suspendWhenHidden = FALSE)
  
  
  observe({
    if (!hasClusterMap()) {
      updatePickerInput(session, "test_entity", selected = "Clusters")
    }
  })
  observe({
    if (!hasClusterMap()) {
      updatePickerInput(session, "cat_entity", selected = "Clusters")
    }
  })
  
  # Auto-detect categorical vs continuous metadata
  observe({
    req(rv$meta_cell)
    meta_cols <- colnames(rv$meta_cell)
    
    categorical_choices <- meta_cols[sapply(rv$meta_cell, function(x)
      is.character(x) || is.factor(x)
    )]
    continuous_choices <- meta_cols[sapply(rv$meta_cell, function(x)
      is.numeric(x) || is.integer(x)
    )]
    
    updatePickerInput(session, "group_var",
                      choices = c("", categorical_choices), selected = "")
    updatePickerInput(session, "cont_var",
                      choices = c("", continuous_choices), selected = "")
  })
  observe({
    req(rv$meta_cell)
    meta_cols <- colnames(rv$meta_cell)
    continuous_choices <- meta_cols[sapply(rv$meta_cell, is.numeric)]
    updatePickerInput(session, "cont_group_var",
                      choices = c("", continuous_choices),
                      selected = "")
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
      reactive(input$main_tab),
      rv
    )
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
  
  run_tests <- eventReactive(input$run_test, {
    req(rv$meta_cell)
    test_type <- input$test_type
    
    # --- Cluster or celltype testing using abundance matrix ---
    abund0 <- rv$clusters$abundance
    if (is.null(abund0)) {
      showNotification("No abundance matrix available. Ensure obj$leiden$abundance is present in the upload.", type = "error", duration = NULL)
      message("run_tests: rv$clusters$abundance is NULL; aborting.")
      return(data.frame(entity = NA, test = NA, p = NA, n = NA))
    }
    if (is.null(rownames(abund0)) || any(!nzchar(rownames(abund0)))) {
      showNotification("Abundance matrix has no valid rownames; cannot map to metadata.", type = "error", duration = NULL)
      message("run_tests: abundance rownames missing/empty; aborting.")
      return(data.frame(entity = NA, test = NA, p = NA, n = NA))
    }
    
    abund <- abund0
    if (input$test_entity == "Celltypes" &&
        !is.null(rv$cluster_map) &&
        all(c("cluster", "celltype") %in% names(rv$cluster_map))) {
      
      cm <- rv$cluster_map
      keep <- cm$cluster %in% colnames(abund)
      cm <- cm[keep, , drop = FALSE]
      if (!nrow(cm)) {
        showNotification("No overlapping clusters to aggregate into celltypes.", type = "error")
        return(data.frame(entity = NA, test = NA, p = NA, n = NA))
      }
      split_idx <- split(cm$cluster, cm$celltype)
      abund <- sapply(split_idx, function(cols) rowSums(abund0[, cols, drop = FALSE]))
      abund <- as.matrix(abund)
    }
    
    if (!nzchar(input$group_var) && !nzchar(input$cont_var)) {
      return(data.frame(entity = NA, test = NA, p = NA, n = NA))
    }
    
    # Map abundance rows (sources) to patient_ID via escaped regex
    sources <- rownames(abund)
    ids <- unique(rv$meta_cell$patient_ID)
    ids_esc <- stringr::str_replace_all(ids, "([\\^$.|?*+()\\[\\]{}\\\\])", "\\\\\\\\1")
    ids_esc <- ids_esc[order(nchar(ids_esc), decreasing = TRUE)]
    pattern <- paste0("(", paste0(ids_esc, collapse = "|"), ")")
    pid <- stringr::str_extract(sources, pattern)
    
    abund_df <- as.data.frame(abund)
    abund_df$patient_ID <- pid
    meta_unique <- rv$meta_cell %>% dplyr::distinct(patient_ID, .keep_all = TRUE)
    abund_df <- abund_df %>% dplyr::left_join(meta_unique, by = "patient_ID")
    
    # Debug: confirm grouping column presence
    if (nzchar(input$group_var)) {
      gv <- input$group_var
      if (!(gv %in% names(abund_df))) {
        message("run_tests: grouping variable not found after join: ", gv)
        showNotification(paste0("Grouping variable '", gv, "' not found after join."), type = "error")
      } else {
        message("run_tests: grouping var head after join:")
        print(utils::head(abund_df[, c("patient_ID", gv)], 10))
      }
    }
    
    abund_long <- abund_df %>%
      tidyr::pivot_longer(cols = colnames(abund), names_to = "entity", values_to = "freq")
    
    res <- abund_long %>%
      dplyr::group_by(entity) %>%
      dplyr::group_modify(~ {
        if (test_type == "Wilcoxon (2-group)") {
          if (!nzchar(input$group_var)) return(data.frame(test = "wilcox", p = NA, n = nrow(.x)))
          g <- droplevels(factor(.x[[input$group_var]]))
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
        } else if (test_type == "Kruskal–Wallis (multi-group)") {
          if (!nzchar(input$group_var)) return(data.frame(test = "kruskal", p = NA, n = nrow(.x)))
          g <- .x[[input$group_var]]
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
        } else {
          if (!nzchar(input$cont_var)) return(data.frame(test = "spearman", rho = NA, p = NA, n = nrow(.x)))
          cont <- input$cont_var
          ct <- spearman_test(.x, cont_var = cont)
          cbind(data.frame(test = "spearman"), ct)
        }
      }) %>%
      dplyr::ungroup()
    
    if (nrow(res) && "p" %in% colnames(res) && nzchar(input$p_adj_method)) {
      adj_col <- paste0(tolower(input$p_adj_method), "_padj")
      res[[adj_col]] <- p.adjust(res$p, method = input$p_adj_method)
    }
    
    # Add metadata column name being tested based on test type
    tested_var <- if (test_type == "Spearman (continuous)") {
      if (nzchar(input$cont_var)) input$cont_var else NA_character_
    } else if (test_type %in% c("Wilcoxon (2-group)", "Kruskal–Wallis (multi-group)")) {
      if (nzchar(input$group_var)) input$group_var else NA_character_
    } else {
      NA_character_
    }
    
    # Determine entity used in this run
    entity_used <- if (!hasClusterMap()) {
      "clusters"
    } else {
      tolower(input$test_entity)
    }
    
    # Determine test short name
    test_map <- c(
      "Wilcoxon (2-group)" = "wilcoxon",
      "Kruskal–Wallis (multi-group)" = "kruskal_wallis",
      "Spearman (continuous)" = "spearman"
    )
    test_used <- test_map[[input$test_type]]
    
    # Determine metadata used based on test type
    metadata_used <- if (input$test_type == "Spearman (continuous)") {
      if (nzchar(input$cont_var)) input$cont_var else "none"
    } else if (input$test_type %in% c("Wilcoxon (2-group)", "Kruskal–Wallis (multi-group)")) {
      if (nzchar(input$group_var)) input$group_var else "none"
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
    adj_col <- paste0(tolower(input$p_adj_method), "_padj")
    base_cols <- c("entity", "metadata", "test", "n", "p")
    if (adj_col %in% names(res)) {
      base_cols <- c(base_cols, adj_col)
    }
    if ("rho" %in% names(res)) {
      base_cols <- c(base_cols, "rho")
    }
    res <- res[, c(base_cols, iqr_cols), drop = FALSE]
    
    res
  })
  
  output$test_table <- renderTable({
    df <- req(run_tests())
    
    # Detect adjusted p-value column
    adj_col <- paste0(tolower(input$p_adj_method), "_padj")
    
    # Format p and adjusted p
    num_cols <- intersect(c("p", adj_col), names(df))
    df[num_cols] <- lapply(df[num_cols], function(x) formatC(x, format = "f", digits = 3))
    
    # Format rho if present
    if ("rho" %in% names(df)) {
      df$rho <- formatC(df$rho, format = "f", digits = 2)
    }
    
    # Order by adjusted p if present
    if (adj_col %in% names(df)) {
      df <- df[order(as.numeric(df[[adj_col]]), na.last = TRUE), ]
    }
    
    df
  }, sanitize.text.function = function(x) x)
  
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
      df <- run_tests()
      req(!is.null(df), nrow(df) > 0)
      # Remove BOM (this was crashing your connection)
      write.csv(df, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  # Update categorical metadata choices when data is ready
  observe({
    req(rv$meta_cell)
    meta_cols <- colnames(rv$meta_cell)
    categorical_choices <- meta_cols[sapply(rv$meta_cell, function(x) is.character(x) || is.factor(x))]
    updatePickerInput(session, "cat_group_var",
                      choices = c("", categorical_choices),
                      selected = "")
  })
  
  cat_plot_data <- eventReactive(input$generate_cat_plots, {
    req(rv$meta_cell, rv$clusters$abundance)
    abund0 <- rv$clusters$abundance
    if (is.null(rownames(abund0)) || any(!nzchar(rownames(abund0)))) {
      showNotification("Abundance matrix has no valid rownames; cannot map to metadata.",
                       type = "error")
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
    
    # Long format
    abund_long <- abund_df %>%
      tidyr::pivot_longer(cols = colnames(abund), names_to = "entity", values_to = "freq") %>%
      mutate(entity = gsub(pattern = "\\n", replacement = " ", x = entity))
    
    # Run test per entity
    test_type <- input$cat_test_type
    group_var <- input$cat_group_var
    res <- abund_long %>%
      dplyr::group_by(entity) %>%
      dplyr::group_modify(~ {
        g <- droplevels(factor(.x[[group_var]]))
        ok <- !is.na(g)
        g <- g[ok]; freq_ok <- .x$freq[ok]
        if (length(unique(g)) < 2) return(data.frame(p = NA))
        if (test_type == "Wilcoxon (2-group)" && length(unique(g)) == 2) {
          wt <- suppressWarnings(wilcox.test(freq_ok ~ g))
          data.frame(p = wt$p.value)
        } else if (test_type == "Kruskal–Wallis (multi-group)") {
          kw <- kruskal.test(freq_ok ~ g)
          data.frame(p = kw$p.value)
        } else {
          data.frame(p = NA)
        }
      }) %>%
      dplyr::ungroup()
    
    # Adjust p-values
    if (nrow(res) && nzchar(input$cat_p_adj_method)) {
      res$padj <- p.adjust(res$p, method = input$cat_p_adj_method)
    }
    
    # Save info for export
    rv$last_cat_info <- list(
      entity   = tolower(input$cat_entity %||% "clusters"),
      group    = tolower(input$cat_group_var %||% "group"),
      test_raw = input$cat_test_type %||% "test"
    )
    
    # Return everything needed for plotting without reading live inputs later
    list(
      data       = abund_long,
      results    = res,
      group_var  = input$cat_group_var,
      use_adj_p  = input$cat_use_adj_p,
      facet_cols = as.numeric(input$cat_max_facets)
    )
  })
  
  output$categorical_plot <- renderPlot({
    cp <- cat_plot_data()
    req(cp)
    abund_long <- cp$data
    res        <- cp$results
    group_var  <- cp$group_var
    use_adj_p  <- cp$use_adj_p
    facet_cols <- cp$facet_cols
    
    p_df <- res %>%
      mutate(
        p_to_show = if (isTRUE(use_adj_p) && "padj" %in% names(res)) padj else p,
        label = paste0("p = ", signif(p_to_show, 3)),
        x = length(unique(abund_long[[group_var]])) / 2 + 0.5,
        y = tapply(abund_long$freq, abund_long$entity, max, na.rm = TRUE)[entity] * 1.05
      )
    
    gg <- ggplot(abund_long, aes(x = .data[[group_var]], y = freq)) +
      geom_boxplot(aes(fill = .data[[group_var]]), alpha = 0.7) +
      facet_wrap(~entity, ncol = facet_cols, scales = "free_y") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      theme_bw(base_size = 18) +
      theme(legend.position = "none") +
      labs(x = group_var, y = "Frequency") +
      geom_text(data = p_df, aes(x = x, y = y, label = label),
                inherit.aes = FALSE, size = 6) +
      theme(strip.text.x = element_text(margin = margin(t = 1.1, b = 1.1)))
    
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
  
  # Event to generate plots
  cont_plot_data <- eventReactive(input$generate_cont_plots, {
    req(rv$meta_cell, rv$clusters$abundance)
    abund0 <- rv$clusters$abundance
    if (is.null(rownames(abund0)) || any(!nzchar(rownames(abund0)))) {
      showNotification("Abundance matrix has no valid rownames; cannot map to metadata.",
                       type = "error")
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
    ids_esc <- stringr::str_replace_all(ids, "([\\^$.|?*+()\\[\\]{}\\\\])", "\\\\\\1")
    ids_esc <- ids_esc[order(nchar(ids_esc), decreasing = TRUE)]
    pattern <- paste0("(", paste0(ids_esc, collapse = "|"), ")")
    pid <- stringr::str_extract(sources, pattern)
    
    abund_df <- as.data.frame(abund)
    abund_df$patient_ID <- pid
    meta_unique <- rv$meta_cell %>% dplyr::distinct(patient_ID, .keep_all = TRUE)
    abund_df <- abund_df %>% dplyr::left_join(meta_unique, by = "patient_ID")
    
    # Long format
    abund_long <- abund_df %>%
      tidyr::pivot_longer(cols = colnames(abund), names_to = "entity", values_to = "freq") %>%
      mutate(entity = gsub("\\n", " ", entity))
    
    # Run Spearman per entity
    cont_var <- input$cont_group_var
    res <- abund_long %>%
      dplyr::group_by(entity) %>%
      dplyr::group_modify(~ {
        ok <- complete.cases(.x$freq, .x[[cont_var]])
        if (!any(ok)) return(data.frame(rho = NA, p = NA))
        ct <- suppressWarnings(cor.test(.x$freq[ok], .x[[cont_var]][ok], method = "spearman"))
        data.frame(rho = unname(ct$estimate), p = ct$p.value)
      }) %>%
      dplyr::ungroup()
    
    # Adjust p-values
    if (nrow(res) && nzchar(input$cont_p_adj_method)) {
      res$padj <- p.adjust(res$p, method = input$cont_p_adj_method)
    }
    
    # Save info for export
    rv$last_cont_info <- list(
      entity   = tolower(input$cont_entity %||% "clusters"),
      group    = tolower(input$cont_group_var %||% "group"),
      test_raw = "spearman"
    )
    
    # Return everything needed for plotting without reading live inputs later
    list(
      data         = abund_long,
      results      = res,
      cont_var     = input$cont_group_var,
      use_adj_p    = input$cont_use_adj_p,
      facet_cols   = as.numeric(input$cont_max_facets)  # capture facet cols here
    )
  })
  
  output$continuous_plot <- renderPlot({
    cp <- cont_plot_data()
    req(cp)
    abund_long <- cp$data
    res        <- cp$results
    cont_var   <- cp$cont_var
    use_adj_p  <- cp$use_adj_p
    facet_cols <- cp$facet_cols
    
    # Annotation data
    p_df <- res %>%
      mutate(
        p_to_show = if (isTRUE(use_adj_p) && "padj" %in% names(res)) padj else p,
        label = paste0("p = ", signif(p_to_show, 3),
                       "\n",
                       "rho = ", signif(rho, 3)),
        x = tapply(abund_long$freq, abund_long$entity,
                   function(v) mean(range(v, na.rm = TRUE)))[entity],
        y = tapply(abund_long[[cont_var]], abund_long$entity,
                   max, na.rm = TRUE)[entity] * 1.05
      )
    
    gg <- ggplot(abund_long, aes(x = freq, y = .data[[cont_var]])) +
      geom_point(alpha = 0.75, pch = 21, color = 'black',
                 fill = 'grey40', stroke = 0.1, size = 3) +
      geom_smooth(method = "lm", se = FALSE, color = "red2") +
      facet_wrap(~entity, ncol = facet_cols, scales = "free") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
      theme_bw(base_size = 18) +
      theme(legend.position = "none") +
      labs(x = "Abundance", y = cont_var) +
      geom_text(data = p_df,
                aes(x = x, y = y, label = label),
                inherit.aes = FALSE,
                size = 5,
                lineheight = 0.85) +
      theme(strip.text.x = element_text(margin = margin(t = 1.1, b = 1.1)))
    
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
  
}

shinyApp(ui, server)
