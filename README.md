# shiny app for exploring results of generated from the FCSimple package
Object, saved as .RData, should contain the following elements, as formatted by FCSimple functions:
  - data (FCSimple::fcs_join)
  - source (FCSimple::fcs_join)
  - run_date (FCSimple::fcs_join)
  - metadata (FCSimple::fcs_join or FCSimple::fcs_update)
  - leiden including calculated abundance matrix (other algorithms not supported, FCSimple::fcs_cluster, FCSimple::fcs_calculate_abundance)
  - leiden clustering heatmap (FCSimple::fcs_cluster_heatmap)
  - umap (FCSimple::fcs_reduce_dimensions)
  - tsne (FCSimple::fcs_reduce_dimensions)
