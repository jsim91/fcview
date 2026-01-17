# shiny app for exploring [FCSimple](https://github.com/jsim91/FCSimple) results

<b>TODO:</b></br>
- update interval plots to correctly identify contrast string for removal from titles

The input object, saved in .RData format, should contain the following elements, as formatted by FCSimple functions unless otherwise noted:
  - data `[required]`
  - source `[required]`
  - run_date `[required]`
  - metadata `[required]`
  - umap `[recommended]`
  - tsne `[recommended]`
  - cluster `[required]`
  - cluster_heatmap `[required]`

For compatibility, the cluster element should be renamed to "cluster" and the heatmap element renamed to "cluster_heatmap", no matter what algorithm was used during the calculation phase. The heatmap plot object can be dropped to shrink file size. Only the underlying heatmap matrix data is required. See layout below for more on this.

**Recommended Object Layout**

Top-level structure

| Element | Type | Description |
| --- | --- | --- |
| `data` | numeric matrix | Feature measurements matrix |
| `source` | character | Sample source identifiers |
| `run_date` | character | Run timestamps |
| `metadata` | data.frame | Sample source level metadata |
| `umap` | list | UMAP coordinates + settings (see below) |
| `tsne` | list | tSNE coordinates + settings (see below) |
| `cluster` | list | Clustering result (see below) |
| `cluster_heatmap` | list | Cluster heatmap result (see below) |

`umap` sub-structure

| Element | Type | Description | Function |
| --- | --- | --- | --- |
| `coordinates` | data.frame | Embedding coordinates | FCSimple::fcs_reduce_dimensions |
| `settings` | list | UMAP parameters | FCSimple::fcs_reduce_dimensions |

`tsne` sub-structure

| Element | Type | Description | Function |
| --- | --- | --- | --- |
| `coordinates` | data.frame | Embedding coordinates | FCSimple::fcs_reduce_dimensions |
| `settings` | list | tSNE parameters | FCSimple::fcs_reduce_dimensions |

`cluster` sub-structure

| Element | Type | Description | Function |
| --- | --- | --- | --- |
| `clusters` | factor | Cluster membership | FCSimple::fcs_cluster |
| `settings` | list | Clustering parameters | FCSimple::fcs_cluster |
| `abundance` | numeric matrix | Sample cluster frequency array | FCSimple::fcs_calculate_abundance |
| `counts` | numeric matrix | Sample cluster counts array | FCSimple::fcs_calculate_abundance |

`cluster_heatmap` sub-structure

| Element | Type | Description | Function |
| --- | --- | --- | --- |
| `heatmap_tile_data` | numeric matrix | Heatmap data | FCSimple::fcs_cluster_heatmap |

<b>NOTE:</b> </br>
You may delete the `obj$cluster_heatmap$heatmap` and `obj$cluster_heatmap$rep_used` elements of `obj$cluster_heatmap` if present; `heatmap` and `rep_used` are not needed for fcview and it will cut down on file size if these are set to NULL. Reminder: `obj$cluster_heatmap` must remain type list after dropping `heatmap` and `rep_used` elements.
