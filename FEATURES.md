# FCView App Features

## General Features

- Data Upload: load an `.RData` object containing required fields (`data`, `source`, `metadata`, `cluster`) via `rdata_upload`.
- Automatic Initialization: parse, validate, and initialize expression, metadata, clusters, embeddings, and abundance/count matrices on upload.
- Downsampling: optional per-upload cell downsampling with user-configurable max.
- Metadata Mapping: robust `patient_ID` extraction and join between per-cell `source` and per-sample metadata.
- Feature Selection UI: multi-select `features_dropdown` with compact per-feature mini UI (`features_mini_ui`) for quick actions.
- Per-Feature Quick Actions: hide features and set type coercion from the mini UI (`hide_mini_<id>`, `type_mini_<id>`).
- Type Coercion: validate and apply coercions (character, factor, numeric/integer by level/character).
- Subsetting Rules: build, apply, preview, and export per-sample subsetting rules (numeric and categorical).
- Pairing Support: picker for a pairing variable, paired-test detection and paired/unpaired UI adjustments.
- Testing & Statistics: run abundance tests (with paired/unpaired options) and view/export results.
- sccomp Integration: run `sccomp_estimate` and `sccomp_test` with custom/simple formula support, plotting, and exports.
- Embeddings & Heatmaps: UMAP/tSNE embedding support and cluster heatmap plotting.
- Feature Export: download subset metadata, cluster frequencies, counts, and sccomp/export artifacts.
- Reactive Safety: mini UI is authoritative; changes reset pairing/subsetting to avoid inconsistent state.

## Comprehensive Features

### Input Validation & Loading

- Structured Input Detection: accepts `.RData` containing at least `data`, `source`, `metadata`, and `cluster`; warns on multiple/none matches.
- Environment-safe Load: loads into a temporary env and selects the single matching object.
- Embedding Handling: converts matrix embeddings to `data.frame` automatically (`umap`, `tsne`).
- Downsampling on Upload: configurable `max_cells_upload` threshold; deterministic sampling with seed and notifications.
- Comprehensive Sanity Checks: validates presence/format of `patient_ID`, abundance rownames, counts matrix, and emits user-facing notifications.

### Data Structures & Caching

- Per-Cell and Per-Sample Stores: stores `rv$expr`, `rv$meta_cell`, `rv$meta_sample`, `rv$meta_sample_original`, `rv$meta_cached`.
- Clusters & Mapping: preserves `clusters` (assignments, settings, abundance), `cluster_map`, and `cluster_heat` objects.
- Counts & Abundance: retains `rv$counts_sample` for `sccomp` and `rv$abundance_sample` for frequency-based analyses.
- Type Coercion State: `rv$type_coercions` stores requested coercion per column; `type_coercion_changed` flags changes.
- Mini UI Persistence: `rv$mini_hide_states` persists per-feature hide toggles across UI interactions.

### Global Settings & Feature Controls

- `features_dropdown`: select which metadata columns to surface in downstream UI.
- `features_mini_ui`: per-selected-feature compact rows showing hide checkbox and type selector; `rv$available_features` derives from mini hide states.
- Single Source of Truth: `features_mini_ui` replaces legacy checkboxes—mini inputs persist directly into `rv` state.

### Type Coercion System

- Validation: `validate_coercions()` inspects column data and offers valid conversion paths.
- Application: `apply_coercion()` performs conversions; coercions applied only to visible (non-hidden) features and never to `patient_ID`.
- Coercion Safety: when coercions change, `pairing_var` and all subsetting are reset and user is notified.

### Subsetting Engine

- Rule-based UI: add multiple rules with stable IDs; supports numeric operators and categorical value picks.
- Dynamic Value UI: numeric sliders and categorical value pickers render depending on selected column type.
- Apply & Preview: intersection/union logic, optional trimming of incomplete pairs, preview sample distribution and subset summary, and compute a unique `subset_id`.
- Exports: `export_subset_meta`, `export_subset_frequencies`, and `export_subset_counts` write current subset to CSV with ordered rows aligned to `patient_ID`.

### Pairing & Paired-Testing Support

- Pairing Variable Picker: `pairing_var` updated with metadata columns; preserved when possible.
- Feasibility Checks: `test_can_pair()` and `cat_can_pair()` verify pairing feasibility.
- UI Adaptation: test selection UIs adapt to show paired tests only when pairing is both selected and feasible.

### Testing & Statistical Workflows

- Abundance Testing: pipeline to merge per-sample abundance with metadata, reshape to long, and run statistical comparisons.
- Paired/Unpaired Options: automatically uses pairing when selected and feasible.
- P-value Adjustment & Results: supports `p_adj_method`; outputs test tables, summaries, and exports.

### sccomp Integration

- Flexible Formula Input: simple mode (group ± covariates) or custom formula text.
- Reference Level Selection and factor cleaning for sccomp.
- Execution & Resource Control: runs `sccomp_estimate()` with specified `cores`, stores results, runs `sccomp_test()`, and cleans temporary files.
- Plots & Downloads: interval plots for estimates/contrasts, CSV and PDF exports.

### Embeddings & Visualization

- UMAP/tSNE: accepts and displays embeddings used in embedding tabs.
- Cluster Heatmap: builds heatmap with configurable palettes and annotation; downloadable PDF.
- Interval & Contrast Plots: generated with patchwork and captioning.

### Modeling & Machine Learning

- Feature Selection Module: `fs_outcome`, `fs_predictors`, cluster-subset pickers.
- Classification & Regression: `lm_outcome`, `lm_predictors`, `reg_outcome`, `reg_predictors`; model training and feature importance outputs.
- Export Models & Results: zipped export of model summaries, performance, features, and figures.

### UX & Notifications

- User Notifications: `showNotification` used extensively to relay status, errors, and warnings.
- Tab Locking: tabs disabled while upload/initialization is in progress via `enableTabs` custom message.
- Status Indicators: reactive outputs such as `hasResults`, `hasSubset`, `hasSccompResults` drive conditional panels and downloads.

### Robust Reactivity & Safety

- Single Source of Truth: mini UI authoritative to minimize reactive loops.
- Ignore `patient_ID`: excluded from coercion and from selectable/coercible feature lists.
- Reset Behavior: changing visibility or type via mini UI resets pairing and clears subsetting.
- Try/Catch Guards: defensive `tryCatch` around UI updates to avoid errors when inputs are not yet present.

### Exporting & Reporting

- Multiple Download Handlers: metadata, counts, frequencies, sccomp CSVs, sccomp PDF plots, model ZIP bundles, and heatmap PDFs.
- Human-Readable Summaries: `global_settings_summary`, pairing summaries, and subsetting preview UI provide concise overviews.

### Developer & Debugging Aids

- Verbose Messages: server prints and messages for tracing upload and sccomp steps.
- Caching for Plots: `cat_plot_cache` and `cont_plot_cache` reactiveVals for expensive plotting reuse.
- Internal Flags: `rv$data_ready` used to gate expensive/reactive updates.

---

If you'd like, I can also add a short link in the `README.md` pointing to this file or embed these sections directly into `README.md`.
