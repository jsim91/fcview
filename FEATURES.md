# fcview App Features

## General Features

- Data Upload: load an `.RData` object containing required fields `data`, `source`, `metadata`, and `cluster`.
- Automatic Initialization: parse, validate, and initialize expression, metadata, clusters, embeddings, and abundance/count matrices on upload.
- Dataset Preview: display detected metadata features available for use along with data type and example values.
- Downsampling: optional per-upload cell downsampling.
- Metadata Mapping: coherent mapping of per-sample metadata with per-cell metadata.
- Metadata Feature Quick Actions: hide features or set custom data types for features (continuous -> categorical or categorical -> continuous).
- Subsetting Rules: build, apply, preview, and export per-sample/feature subsetting rules (numeric and categorical).
- Pairing Support: select the metadata variable that should be used for paired testing if repeated measures were taken.
- Paired Testing Autodetection: app backend detects if paired testing can be done and updates UI options in real-time.
- Data Modeling Tasks: mult-method support for feature selection and prediction modeling along continuous and categorical endpoints.
- Testing & Statistics: run non-parameteric abundance tests including Wilcoxon rank sum (Mann-Whitney, unpaired), Wilcoxon signed rank (paired), Kruskal-Wallis (multi-group unpaired), Friedman (multi-group paired).
- sccomp Integration: run sccomp's Bayesian testing framework with custom/simple formula support, plotting, and exports.
- Embeddings & Heatmaps: UMAP/tSNE embedding support with custom faceting and cluster heatmap plotting.
- Feature Export: download metadata matching current subset rules, cluster frequencies, and cluster counts.

## Comprehensive Features & Backend References

### Input Validation & Loading

- Structured Input Detection: accepts `.RData` containing at least `data`, `source`, `metadata`, and `cluster`; warns on multiple/none matches.
- Environment-safe Load: loads into a temporary env and places list elements into either a reactive or static context, as needed.
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
- Application: `apply_coercion()` performs data type conversions; coercions applied only to exposed features (excludes `run_date` and `patient_ID`).
- Real-time Effects: testing and modeling methods are updated when data types change in real-time.
- Coercion Safety: when coercions change, `pairing_var` and all subsetting are reset and user is notified.

### Pairing & Paired-Testing Support

- Pairing Variable Picker: `pairing_var` updated with metadata columns; preserved when possible.
- Feasibility Checks: `test_can_pair()` and `cat_can_pair()` verify pairing feasibility.
- UI Adaptation: test selection UIs adapt to show paired tests only when pairing is both selected and feasible.

### Subsetting Engine

- Rule-based UI: add multiple rules with stable IDs; supports numeric operators and categorical value picks.
- Dynamic Value UI: numeric sliders and categorical value pickers render depending on selected column type.
- Apply & Preview: intersection/union logic, optional trimming of incomplete pairs to prevent paired testing errors, preview sample distribution and subset summary, and compute a unique `subset_id`.
- Exports: `export_subset_meta`, `export_subset_frequencies`, and `export_subset_counts` write current subset to CSV with ordered rows aligned to `patient_ID`.

### Testing & Statistical Workflows

- Abundance Testing: pipeline to merge per-sample abundance with metadata, reshape to long, and run statistical comparisons.
- Paired/Unpaired Options: automatically uses pairing when selected and feasible.
- P-value Adjustment & Results: supports `p_adj_method`; outputs test tables, summaries, and exports.
- Customization: select custom colors and draw options for testing figures.

### sccomp Integration

- Flexible Formula Input: simple mode (group ± covariates) or custom formula text for complex mixed effects modeling or zero-intercept format.
- Reference Control: set custom reference levels for all terms in the model.
- Execution & Resource Control: runs `sccomp_estimate()` with specified `cores`, stores results, and runs `sccomp_test()`.
- Custom Contrasts: run explicit contrasts derived from current `sccomp_estimate` result.
- Plots & Downloads: interval plots for estimates/contrasts, CSV and PDF exports.

### Embeddings & Cluster Phenotyping

- UMAP/tSNE: accepts and displays embeddings used in embedding tabs.
- Faceting: tile embedding plots by metadata feature with equal sampling to avoid overrepresentation bias.
- Cluster Heatmap: builds heatmap with configurable palettes and annotation; downloadable PDF.
- Figure Creation: full control over PDF export of rendered facets and heatmap.

### Modeling & Machine Learning

- Supported Models: Strict Ridge, Elastic net with tunable alpha, Random Forest, and simple linear models depending on current task.
- Feature Selection Module: `fs_outcome`, `fs_predictors`, cluster-subset pickers.
- Classification & Regression: `lm_outcome`, `lm_predictors`, `reg_outcome`, `reg_predictors`; model training and feature importance outputs.
- Export Models & Results: zipped export of model summaries, performance, features, and figures.

### User Experience & Notifications

- User Notifications: `showNotification` used extensively to relay status, errors, and warnings.
- Tab Locking: tabs disabled while upload/initialization is in progress via `enableTabs` custom message.
- Status Indicators: reactive outputs such as `hasResults`, `hasSubset`, `hasSccompResults` drive conditional panels and downloads.

### Robust Reactivity & Safety

- Robust Feature Control: feature quick actions selections take precedence over default settings for statistical testing and modeling tasks.
- Persistent Metadata Description: Home tab available metadata table held in a static context for reference.
- Ignore `patient_ID`: excluded from coercion and from selectable/coercible feature lists.
- Reset Behavior: changing visibility or type via mini UI resets pairing and clears subsetting to prevent unexpected behavior.
- Try/Catch Guards: defensive `tryCatch` around UI updates to avoid errors when inputs are not yet present.

### Exporting & Reporting

- Multiple Download Handlers: metadata, counts, frequencies, sccomp CSVs, sccomp PDF plots, model ZIP bundles, heatmap PDFs, and rendered embeddings PDFs.
- Human-Readable Summaries: `global_settings_summary`, pairing summaries, and subsetting preview UI provide concise overviews.
