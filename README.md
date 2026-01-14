# shiny app for exploring [FCSimple](https://github.com/jsim91/FCSimple) results
The input object, saved in .RData format, should contain the following elements (some are required), as formatted by FCSimple functions:
  - data (FCSimple::fcs_join)
  - source (FCSimple::fcs_join)
  - run_date (FCSimple::fcs_join)
  - metadata (FCSimple::fcs_join or FCSimple::fcs_update)
  - leiden including calculated abundance matrix (other algorithms not supported, FCSimple::fcs_cluster, FCSimple::fcs_calculate_abundance)
  - leiden clustering heatmap (FCSimple::fcs_cluster_heatmap)
  - umap (FCSimple::fcs_reduce_dimensions)
  - tsne (FCSimple::fcs_reduce_dimensions)

Load time scales with object size. Consider subsampling obj$data and all element-wise matching list entries to reduce object weight. The app will use the stored abundance matrix in obj$leiden$abundance for statistical testing. Downsampling the object elements will make loading faster and will not impact the statistical testing results. Recommended downsample size: 300000 if nrow(obj$data)>300000.

**Recommended Object Layout**

This app expects a named list-like object (example: `obj_small`) with the following top-level elements.

Top-level structure

| Element | Type / Length | Description |
| --- | --- | --- |
| `data` | numeric (example length: 2292768) | Main measurement vector / matrix |
| `source` | character (example length: 95532) | Source identifiers |
| `run_date` | character (example length: 95532) | Run timestamps |
| `cluster` | list | Clustering results (see below) |
| `umap` | list | UMAP coordinates + settings (see below) |
| `cluster_heatmap` | list | Heatmap tile data and related objects |
| `versions` | list | R/Python/session info |
| `metadata` | data.frame / list | Metadata mapped to events or samples |
| `collection_instrument` | character | Instrument identifier |
| `object_history` | character | Notes / provenance |

`cluster` sub-structure

| Element | Type / Length | Description |
| --- | --- | --- |
| `clusters` | factor (example length: 95532) | Cluster assignment per event |
| `settings` | list (4) | Flags and `nClus` number |
| `abundance` | numeric (example length: 5490) | Abundance vector |

`cluster$settings` (example)

| Name | Type | Meaning |
| --- | --- | --- |
| `compensate` | logical | Compensation applied? |
| `transform` | logical | Transform applied? |
| `silent` | logical | Suppress messages? |
| `nClus` | numeric | Number of clusters |

`umap` sub-structure

| Element | Type / Length | Description |
| --- | --- | --- |
| `coordinates` | data.frame | X/Y (or multi-dim) coordinates |
| `settings` | list (see below) | UMAP parameters |

`umap$settings` (example keys)

- `use_rep` (character)
- `language` (character)
- `init` (character)
- `low_memory` (character)
- `random_state` (numeric)
- `num_neighbors` (numeric)
- `min_dist` (numeric)
- `transform_seed` (numeric)
- `verbose` (character)

`versions` (top-level view)

| Element | Type | Description |
| --- | --- | --- |
| `R` | list | R system version info |
| `Python` | character | Python version string |
| `Rsession` | list | `sessionInfo()`-like structure |
| `pip_list` | data.frame | pip packages / versions |
