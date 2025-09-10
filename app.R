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
    library(ggrepel)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(data.table)
    library(ComplexHeatmap)
    library(circlize)
    library(viridis)
    library(scales)
    library(sf)         # point-in-polygon
    library(jsonlite)   # gate import/export
  })
})

options(shiny.maxRequestSize = .Machine$integer.max) # allow large uploads; adjust as needed

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
    if (!is.data.frame(obj$umap$coordinates) || nrow(obj$umap$coordinates) != n) stop("umap$coordinates must align with cells.")
  }
  if ("tsne" %in% names(obj)) {
    if (!is.data.frame(obj$tsne$coordinates) || nrow(obj$tsne$coordinates) != n) stop("tsne$coordinates must align with cells.")
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
  candidates <- intersect(tolower(colnames(metadata)), c("patientid","patient_id","patient","source","subject","sample","id"))
  if (length(candidates)) {
    hit <- candidates[which.max(nchar(candidates) == nchar(candidates))] # just pick first
    # Map back to original case
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
    embedding = embedding,               # "UMAP" or "tSNE"
    polygon = polygon,                   # list(x=..., y=...) if drawn
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
      column(3,
             pickerInput(ns("color_by"), "Color by", choices = NULL),
             pickerInput(ns("split_by"), "Split by (optional) — metadata", choices = NULL, options = list(`none-selected-text`="None")),
             sliderInput(ns("max_facets"), "Max facets", min = 2, max = 16, value = 6, step = 1),
             hr(),
             radioButtons(ns("gate_mode"), "Gating mode", choices = c("Lasso select", "Draw polygon"), selected = "Draw polygon"),
             textInput(ns("gate_name"), "Gate name", value = ""),
             colourInput(ns("gate_color"), "Gate color", value = "#E45756"), 
             actionButton(ns("save_gate"), "Save gate"),
             actionButton(ns("clear_selection"), "Clear current selection"),
             hr(),
             pickerInput(ns("overlay_gate"), "Overlay gate(s)", choices = NULL, multiple = TRUE)
      ),
      column(9,
             plotlyOutput(ns("embed_plot"), height = "650px")
      )
    ),
    hr(),
    h4("Gate phenotype and abundance"),
    fluidRow(
      column(4,
             pickerInput(ns("phenotype_gate"), "Gate for phenotype", choices = NULL),
             pickerInput(ns("phenotype_fun"), "Summary", choices = c("median", "p90"), selected = "median")
      ),
      column(8, plotOutput(ns("inout_heatmap"), height = "400px"))
    ),
    hr(),
    h4("Abundance testing"),
    fluidRow(
      column(4,
             pickerInput(ns("test_entity"), "Entity", choices = c("Selected gate(s)", "Clusters", "Celltypes")),
             pickerInput(ns("test_gate"), "Gate(s)", choices = NULL, multiple = TRUE),
             pickerInput(ns("group_var"), "Grouping factor (metadata)", choices = NULL),
             pickerInput(ns("cont_var"), "Continuous metadata", choices = NULL),
             radioButtons(ns("test_type"), "Test", choices = c("Wilcoxon (2-group)","Kruskal–Wallis (multi-group)","Spearman (continuous)")),
             pickerInput(ns("unit_var"), "Aggregation unit", choices = NULL),
             checkboxInput(ns("apply_bh"), "BH adjust across entities", TRUE),
             actionButton(ns("run_test"), "Run tests")
      ),
      column(8,
             plotOutput(ns("abund_plot"), height = "300px"),
             tableOutput(ns("test_table"))
      )
    )
  )
}

# ---- Embedding module server ----
EmbeddingServer <- function(id, embedding_name, coords, expr, meta_cell, clusters, cluster_map, gate_store, active_tab, rv) {
  moduleServer(id, function(input, output, session) {
    inputs_initialized <- reactiveVal(FALSE)
    
    observe({
      req(!inputs_initialized())
      req(expr, meta_cell)
      req(active_tab() == embedding_name)  # Ensures tab is visible
      
      isolate({
        numeric_markers <- colnames(expr)
        meta_cols <- setdiff(colnames(meta_cell), c(".cell"))
        cont_choices <- meta_cols[sapply(meta_cell[meta_cols], is.numeric)]
        unit_candidates <- intersect(c("PatientID", "patient_ID", "patient", "source", "RunDate", "run_date"), colnames(meta_cell))
        unit_default <- if (length(unit_candidates)) unit_candidates[1] else meta_cols[1]
        
        updatePickerInput(session, "color_by",
                          choices = c(numeric_markers, meta_cols),
                          selected = numeric_markers[1])
        updatePickerInput(session, "split_by", choices = c("", meta_cols), selected = "")
        updatePickerInput(session, "group_var", choices = meta_cols)
        updatePickerInput(session, "cont_var", choices = cont_choices)
        updatePickerInput(session, "unit_var", choices = meta_cols, selected = unit_default)
        
        message(sprintf("Picker inputs initialized for %s", embedding_name))
        inputs_initialized(TRUE)
      })
    })
    
    
    ui_ready <- reactive({
      input$color_by  # Will be NULL until the UI is rendered
      TRUE
    })
    
    ns <- session$ns
    df <- reactive({
      req(nrow(coords) == nrow(meta_cell), nrow(coords) == nrow(expr))
      out <- tibble::as_tibble(coords)
      names(out)[1:2] <- c("x","y")
      out$.cell <- seq_len(nrow(out))
      out$cluster <- clusters$assignments
      out$celltype <- if (!is.null(cluster_map)) {
        cluster_map$celltype[match(out$cluster, cluster_map$cluster)]
      } else {
        as.character(out$cluster)
      }
      out
    })
    
    # observe({
    #   req(expr, meta_cell)
    #   invalidateLater(500, session)  # Give the UI time to render
    #   
    #   isolate({
    #     numeric_markers <- colnames(expr)
    #     meta_cols <- setdiff(colnames(meta_cell), c(".cell"))
    #     cont_choices <- meta_cols[sapply(meta_cell[meta_cols], is.numeric)]
    #     unit_candidates <- intersect(c("PatientID", "patient_ID", "patient", "source", "RunDate", "run_date"), colnames(meta_cell))
    #     unit_default <- if (length(unit_candidates)) unit_candidates[1] else meta_cols[1]
    #     
    #     updatePickerInput(session, "color_by",
    #                       choices = c(numeric_markers, meta_cols),
    #                       selected = numeric_markers[1])
    #     updatePickerInput(session, "split_by", choices = c("", meta_cols), selected = "")
    #     updatePickerInput(session, "group_var", choices = meta_cols)
    #     updatePickerInput(session, "cont_var", choices = cont_choices)
    #     updatePickerInput(session, "unit_var", choices = meta_cols, selected = unit_default)
    #     
    #     message("Picker inputs updated after UI became ready")
    #   })
    # })
    
    current_sel <- reactiveVal(integer(0))
    
    # Observers start suspended; we resume them after the first real plot is rendered
    sel_observer_resumed  <- reactiveVal(FALSE)
    rel_observer_resumed  <- reactiveVal(FALSE)
    
    sel_obs <- observeEvent(
      event_data("plotly_selected", source = ns("embed")),
      suspended = TRUE,
      ignoreInit = TRUE,
      {
        if (input$gate_mode != "Lasso select") return()
        sel <- event_data("plotly_selected", source = ns("embed"))
        req(!is.null(sel), nrow(sel) > 0)
        current_sel(unique(as.integer(sel$customdata)))
      }
    )
    
    rel_obs <- observeEvent(
      event_data("plotly_relayout", source = ns("embed")),
      suspended = TRUE,
      ignoreInit = TRUE,
      {
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
    
    output$embed_plot <- renderPlotly({
      message(sprintf("Tab: %s | input$color_by: %s", active_tab(), input$color_by %||% "NULL"))
      
      message(sprintf("renderPlotly triggered for %s", embedding_name))
      message(sprintf("color_by: %s", input$color_by))
      message(sprintf("df rows: %d", nrow(df())))
      message(sprintf("df x range: [%f, %f]", min(df()$x, na.rm = TRUE), max(df()$x, na.rm = TRUE)))
      message(sprintf("df y range: [%f, %f]", min(df()$y, na.rm = TRUE), max(df()$y, na.rm = TRUE)))
      message(sprintf("Rendering plot for %s with color_by: %s", embedding_name, input$color_by))
      
      if (is.null(input$color_by)) {
        message("input$color_by is NULL — skipping plot render")
        return(NULL)
      }
      
      # Only render when the corresponding tab is active
      req(active_tab() == embedding_name)
      req(expr, meta_cell, coords, input$color_by)
      
      color_by <- input$color_by
      
      # Validate color_by column
      valid_cols <- c(colnames(expr), colnames(meta_cell))
      message(sprintf("Is color_by valid? %s", color_by %in% valid_cols))
      message(sprintf("Valid columns: %s", paste(valid_cols, collapse = ", ")))
      message(sprintf("input$color_by: %s", color_by))
      if (!(color_by %in% valid_cols)) {
        message("Invalid color_by:", color_by)
        return(plotly_empty(type = "scatter", mode = "markers", source = ns("embed")) %>%
                 layout(xaxis = list(title = paste0(embedding_name, "1")),
                        yaxis = list(title = paste0(embedding_name, "2"))))
      }
      
      dd <- df()
      if (nrow(dd) == 0 || all(is.na(dd$x)) || all(is.na(dd$y))) {
        message("Empty or invalid coordinates in df()")
        return(plotly_empty(type = "scatter", mode = "markers", source = ns("embed")) %>%
                 layout(xaxis = list(title = paste0(embedding_name, "1")),
                        yaxis = list(title = paste0(embedding_name, "2"))))
      }
      
      # Colors
      if (color_by %in% colnames(expr)) {
        vals <- expr[, color_by]
        cols <- col_numeric(viridis(256), domain = range(vals, na.rm = TRUE))(pct_clip(vals))
      } else {
        vals <- meta_cell[[color_by]]
        if (is.numeric(vals)) {
          cols <- col_numeric(viridis(256), domain = range(vals, na.rm = TRUE))(vals)
        } else {
          levs <- unique(as.character(vals))
          pal <- setNames(viridis(max(2, length(levs))), levs)
          cols <- pal[as.character(vals)]
        }
      }
      
      # Create base plot
      p <- plot_ly(source = ns("embed"), type = "scattergl", mode = "markers",
                   x = dd$x, y = dd$y, customdata = dd$.cell,
                   marker = list(size = 4, color = cols, opacity = 0.85)) %>%
        layout(xaxis = list(title = paste0(embedding_name, "1")),
               yaxis = list(title = paste0(embedding_name, "2")),
               dragmode = if (input$gate_mode == "Lasso select") "lasso" else "zoom")
      
      # Overlay gates
      glist <- gate_store$list()
      if (length(glist) > 0) {
        overlay_names <- input$overlay_gate
        active <- if (length(overlay_names)) glist[names(glist) %in% overlay_names] else list()
        for (g in active) {
          idx <- dd$.cell %in% g$cells
          p <- add_trace(p, x = dd$x[idx], y = dd$y[idx], type = "scattergl", mode = "markers",
                         marker = list(size = 5, color = g$color, opacity = 0.95),
                         name = g$name, inherit = FALSE, showlegend = TRUE)
        }
      }
      
      # Add draw tools
      mode_buttons <- if (input$gate_mode == "Draw polygon") {
        c("drawclosedpath", "eraseshape", "lasso2d", "select2d", "toImage")
      } else c("lasso2d", "select2d", "toImage")
      
      p <- config(p, modeBarButtonsToAdd = mode_buttons, edits = list(shapePosition = TRUE))
      
      # Optional: log rendering
      message(sprintf("Rendering plot for %s with color_by: %s", embedding_name, color_by))
      
      return(p)
    })
    
    # Ensure the plot renders even when the tab is hidden
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
      in_vals  <- apply(expr[idx, , drop = FALSE], 2, fun, na.rm = TRUE)
      out_vals <- apply(expr[-idx, , drop = FALSE], 2, fun, na.rm = TRUE)
      M <- rbind(In = in_vals, Out = out_vals)
      Mz <- t(scale(t(M)))
      Heatmap(Mz, name = "z", cluster_rows = FALSE, cluster_columns = TRUE,
              row_names_side = "left", col = colorRamp2(c(-2,0,2), viridis(3)))
    })
    
    # Abundance testing
    observe({
      updatePickerInput(session, "test_gate", choices = names(gate_store$list()))
      updatePickerInput(session, "overlay_gate", choices = names(gate_store$list()))
      updatePickerInput(session, "phenotype_gate", choices = names(gate_store$list()))
    })
    
    run_tests <- eventReactive(input$run_test, {
      unit_var <- req(input$unit_var)
      group_var <- req(input$group_var)
      test_type <- input$test_type
      
      if (input$test_entity == "Selected gate(s)") {
        req(length(input$test_gate) > 0)
        tests <- lapply(input$test_gate, function(gn) {
          gate <- gate_store$list()[[gn]]
          dfreq <- freq_by(gate$cells, meta_cell, group_var, unit_var)
          if (test_type == "Wilcoxon (2-group)") {
            g <- droplevels(factor(dfreq[[group_var]]))
            if (length(levels(g)) != 2) return(data.frame(entity = gn, test = "wilcox", p = NA, n = nrow(dfreq)))
            wt <- wilcox.test(freq ~ g, data = dfreq)
            data.frame(entity = gn, test = "wilcox", p = wt$p.value, n = nrow(dfreq))
          } else if (test_type == "Kruskal–Wallis (multi-group)") {
            kw <- kruskal.test(freq ~ dfreq[[group_var]], data = dfreq)
            data.frame(entity = gn, test = "kruskal", p = kw$p.value, n = nrow(dfreq))
          } else {
            cont <- req(input$cont_var)
            ct <- spearman_test(dfreq %>% rename(!!cont := all_of(cont)), cont_var = cont)
            cbind(data.frame(entity = gn, test = "spearman"), ct)
          }
        })
        out <- do.call(rbind, tests)
        if (nrow(out) && isTRUE(input$apply_bh) && "p" %in% colnames(out))
          out$padj <- p.adjust(out$p, method = "BH")
        out
      } else {
        entity_var <- if (input$test_entity == "Clusters") "cluster" else "celltype"
        dd <- df()
        dd$entity <- if (entity_var == "cluster") dd$cluster else dd$celltype
        # ... continue as in your original
      }
    })
  })
}

# ---- UI ----
ui <- navbarPage("CyTOF Explorer", id = "main_tab",
                 tabPanel("Home",
                          sidebarLayout(
                            sidebarPanel(
                              h4("Upload inputs"),
                              fileInput("rdata_upload", "Upload .RData (contains shinyAppInput)", accept = ".RData"),
                              fileInput("json_upload", "Upload gates (.json)", accept = ".json"),
                              hr(),
                              uiOutput("id_column_ui"),
                              actionButton("apply_id_mapping", "Apply ID mapping"),
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
                 # tabPanel("UMAP", uiOutput("umap_content")),
                 # tabPanel("tSNE", uiOutput("tsne_content")),
                 tabPanel("Clusters",
                          fluidRow(
                            column(3,
                                   pickerInput("cluster_order_by", "Order markers by", choices = c("variance","none"), selected = "variance")
                            ),
                            column(9, plotOutput("cluster_heatmap", height = "650px"))
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
    id_col = NULL
  )
  gate_store <- GateStore()
  
  # Handle RData upload
  observeEvent(input$rdata_upload, {
    req(input$rdata_upload)
    e <- new.env(parent = emptyenv())
    load(input$rdata_upload$datapath, envir = e)
    if (!"shinyAppInput" %in% ls(e)) {
      showNotification("No object named 'shinyAppInput' found in the .RData", type = "error")
      return()
    }
    obj <- e$shinyAppInput
    # Try to guess ID column now; do initial validation without id_col name requirement
    validateInput(obj, id_col = NULL)
    rv$obj <- obj
    
    # Suggest ID column mapping
    guessed <- guess_id_col(obj$metadata, obj$source)
    rv$id_col <- guessed
    
    showNotification(paste0("Loaded shinyAppInput. Guessed ID column: ", guessed %||% "none"), type = "message")
  }, ignoreInit = FALSE)
  
  # UI to choose metadata ID column to match `source`
  output$id_column_ui <- renderUI({
    req(rv$obj)
    selectInput("id_col_select", "Metadata ID column matching 'source'",
                choices = colnames(rv$obj$metadata), selected = rv$id_col %||% colnames(rv$obj$metadata)[1])
  })
  
  # Apply ID mapping and construct app datasets
  observeEvent(input$apply_id_mapping, {
    req(rv$obj)
    id_col <- input$id_col_select
    validateInput(rv$obj, id_col = id_col)
    
    expr <- rv$obj$data
    source_vec <- rv$obj$source
    run_date <- rv$obj$run_date %||% NULL
    meta <- rv$obj$metadata
    
    meta_cell <- data.frame(
      source = source_vec,
      RunDate = if (!is.null(run_date)) run_date else NA
    )
    colnames(meta_cell)[1] <- "source"
    meta_cell <- meta_cell %>%
      left_join(meta %>% mutate(.__id = .data[[id_col]]), by = c("source" = ".__id"))
    
    if (!("PatientID" %in% names(meta_cell))) {
      meta_cell$PatientID <- meta_cell$source
    }
    if ("patient_ID" %in% names(meta_cell) && !"patient_ID" %in% c("source","PatientID")) {
      # keep original too
    }
    if ("RunDate" %in% names(meta_cell)) meta_cell$RunDate <- as.factor(meta_cell$RunDate)
    
    clusters <- list(assignments = rv$obj$leiden$clusters,
                     settings = rv$obj$leiden$settings %||% list())
    cluster_map <- if (hasClusterMapping(rv$obj)) rv$obj$cluster_mapping else NULL
    
    UMAP <- if (hasUMAP(rv$obj)) {
      list(coords = rv$obj$umap$coordinates, settings = rv$obj$umap$settings)
    } else NULL
    
    tSNE <- if (hasTSNE(rv$obj)) {
      list(coords = rv$obj$tsne$coordinates, settings = rv$obj$tsne$settings)
    } else NULL
    
    cluster_heat <- if (hasHeatmap(rv$obj)) rv$obj$leiden_heatmap$heatmap_tile_data else NULL
    pop_size <- if (hasHeatmap(rv$obj)) rv$obj$leiden_heatmap$population_size else NULL
    rep_used <- if (hasHeatmap(rv$obj)) rv$obj$leiden_heatmap$rep_used else NA
    
    # Store
    rv$expr <- expr
    rv$meta_cell <- meta_cell
    rv$clusters <- clusters
    rv$cluster_map <- cluster_map
    rv$UMAP <- UMAP
    rv$tSNE <- tSNE
    rv$cluster_heat <- cluster_heat
    rv$pop_size <- pop_size
    rv$rep_used <- rep_used
    
    EmbeddingServer("umap", "UMAP", rv$UMAP$coords, rv$expr, rv$meta_cell, rv$clusters, rv$cluster_map, gate_store, reactive(input$main_tab), rv)
    EmbeddingServer("tsne", "tSNE", rv$tSNE$coords, rv$expr, rv$meta_cell, rv$clusters, rv$cluster_map, gate_store, reactive(input$main_tab), rv)
    
    # ✅ Reset color_by to a valid default
    valid_cols <- c(colnames(expr), colnames(meta_cell))
    default_color <- valid_cols[1]
    updatePickerInput(session, "color_by", choices = valid_cols, selected = default_color)
    
    showNotification("Metadata mapping applied. App initialized.", type = "message")
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
    pickerInput("gates_to_export", "Select gates to export", choices = names(gl), multiple = TRUE)
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
  
  # Conditional tabs for embeddings
  output$umap_content <- renderUI({
    req(input$main_tab == "UMAP")  # Only render when UMAP tab is active
    req(rv$UMAP, rv$expr, rv$meta_cell, rv$clusters)
    message("Rendering UMAP UI")
    EmbeddingUI("umap", "UMAP")
  })
  output$tsne_content <- renderUI({
    req(input$main_tab == "tSNE")  # Only render when tSNE tab is active
    req(rv$tSNE, rv$expr, rv$meta_cell, rv$clusters)
    message("Rendering tSNE UI")
    EmbeddingUI("tsne", "tSNE")
  })
  
  # Launch embedding modules as soon as data is ready
  # observeEvent(rv$UMAP, {
  #   req(rv$UMAP, rv$expr, rv$meta_cell, rv$clusters)
  #   EmbeddingServer("umap", "UMAP", rv$UMAP$coords, rv$expr, rv$meta_cell, rv$clusters, rv$cluster_map, gate_store, reactive(input$main_tab), rv)
  # })
  # observeEvent(rv$tSNE, {
  #   req(rv$tSNE, rv$expr, rv$meta_cell, rv$clusters)
  #   EmbeddingServer("tsne", "tSNE", rv$tSNE$coords, rv$expr, rv$meta_cell, rv$clusters, rv$cluster_map, gate_store, reactive(input$main_tab), rv)
  # })
  
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
    Heatmap(M, name = "expr",
            cluster_rows = TRUE, cluster_columns = TRUE,
            right_annotation = ranno,
            row_names_side = "left",
            col = colorRamp2(c(min(M), median(M), max(M)), viridis(3)))
  })
}

shinyApp(ui, server)
