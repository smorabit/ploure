#' Simulated single-cell table
#'
#' 3,000 simulated cells for examples and tests (no real data): nine clusters of three lineages on a 2D embedding,
#' a separate embedding of the lymphoid cells, nine patients in three disease groups, QC metrics, a continuous score
#' and eight marker genes. Made by `data-raw/sim_cells.R`.
#'
#' @format A data frame with one row per cell:
#' \describe{
#'   \item{cell}{Cell name.}
#'   \item{umap_1, umap_2}{Embedding of all cells.}
#'   \item{sub_umap_1, sub_umap_2}{Embedding of the lymphoid cells (NA for the others).}
#'   \item{cluster, lineage}{Cluster (factor) and its lineage.}
#'   \item{patient, disease, sample}{Patient, disease group (LC, BC, OC) and sample.}
#'   \item{n_umi, n_genes, pct_mt}{QC metrics.}
#'   \item{score}{A continuous score, high in CD8 T and NK cells.}
#'   \item{CD3E, CD8A, NKG7, MS4A1, LYZ, CD14, COL1A1, PECAM1}{Log-normalized-like expression of marker genes.}
#' }
'sim_cells'

#' Simulated QC flow tables
#'
#' Example inputs of [QCFlowChart()] (no real data): `sim_qc_steps` has a trunk of three shared steps (branch
#' `'all'`) that fans out into three lineages, one of which has a second review stage; `sim_qc_removals` gives the
#' removed cells of some steps by category. Made by `data-raw/sim_qc_flow.R`.
#'
#' @format `sim_qc_steps`: a data frame with columns `branch`, `step`, `n_in`, `n_out`, `label` and `excluded`.
#'   `sim_qc_removals`: a data frame with columns `branch`, `step`, `cluster`, `category` and `n_removed`.
'sim_qc_steps'

#' @rdname sim_qc_steps
'sim_qc_removals'
