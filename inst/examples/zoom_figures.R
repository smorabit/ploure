#---------------------------------------------------------#
# Zoom-in embeddings from simulated data
#
# The four kinds of PLEURA zoom figure (v2/overview/malignant_zoom_umap.Rmd), drawn with simulated cells:
# - the lymphoid subset by cluster
# - by patient, colored by a grouped palette (hue = disease, shade = patient), labeled with each patient's disease,
#   with the disease families as the legend
# - by a continuous score, with percentile limits and low cells in grey
# - the Seurat path, with a gene
#
# usage: Rscript inst/examples/zoom_figures.R [output_dir]
#---------------------------------------------------------#

library(ploure)
library(ggplot2)

out_dir <- commandArgs(TRUE)[1]
if(is.na(out_dir)) out_dir <- tempdir()
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cluster_colors <- setNames(scales::hue_pal()(nlevels(sim_cells$cluster)), levels(sim_cells$cluster))
lymphoid <- sim_cells[sim_cells$lineage == 'Lymphoid', ]

zoom_args <- list(
    full = sim_cells, full_reduction = 'umap', full_color_by = 'cluster', full_colors = cluster_colors,
    subset = lymphoid, subset_reduction = c('sub_umap_1', 'sub_umap_2'),
    titles = c('All cells', sprintf('Lymphoid cells, n = %s', scales::comma(nrow(lymphoid))))
)

#---------------------------------------------------------#
# categorical
#---------------------------------------------------------#

p_cluster <- do.call(ZoomEmbedding, c(zoom_args, list(
    subset_color_by = 'cluster', subset_colors = cluster_colors, legend_titles = c('Cluster', 'Cluster')
)))

cairo_pdf(file.path(out_dir, 'zoom_cluster.pdf'), width = 13, height = 7.5)
print(p_cluster)
invisible(dev.off())

#---------------------------------------------------------#
# by patient: grouped palette, labels and a family legend
#---------------------------------------------------------#

disease_names <- c(LC = 'Lung cancer (LC)', BC = 'Breast cancer (BC)', OC = 'Ovarian cancer (OC)')
patient_colors <- GroupedPalette(
    sim_cells$patient, sim_cells$disease,
    ramps = list(
        LC = c('#c7e9c0', '#41ab5d', '#00441b'),
        BC = c('#ffeda0', '#fd8d3c', '#800026'),
        OC = c('#f1b6da', '#c51b7d', '#8c510a')
    ),
    family_labels = disease_names
)
patients <- unique(sim_cells[c('patient', 'disease')])
patient_labels <- setNames(paste(patients$patient, patients$disease), patients$patient)

p_patient <- do.call(ZoomEmbedding, c(zoom_args, list(
    subset_color_by = 'patient', subset_colors = patient_colors, subset_label_by = 'patient',
    label_map = patient_labels, legend_colors = list(full = NULL, subset = attr(patient_colors, 'families')),
    legend_ncol = c(3, 1), legend_titles = c('Cluster', 'Patient, colored by disease')
)))

cairo_pdf(file.path(out_dir, 'zoom_patient.pdf'), width = 13, height = 7.5)
print(p_patient)
invisible(dev.off())

#---------------------------------------------------------#
# continuous score
#---------------------------------------------------------#

p_score <- do.call(ZoomEmbedding, c(zoom_args, list(
    subset_color_by = 'score', plot_min = 'q1', plot_max = 'q99', legend_titles = c('Cluster', 'Score')
)))

cairo_pdf(file.path(out_dir, 'zoom_score.pdf'), width = 13, height = 7.5)
print(p_score)
invisible(dev.off())

#---------------------------------------------------------#
# Seurat object with a gene
#---------------------------------------------------------#

if(requireNamespace('SeuratObject', quietly = TRUE)){
    genes <- c('CD3E', 'CD8A', 'NKG7', 'MS4A1', 'LYZ', 'CD14', 'COL1A1', 'PECAM1')
    X <- Matrix::Matrix(t(as.matrix(lymphoid[, genes])), sparse = TRUE)
    colnames(X) <- lymphoid$cell
    meta <- lymphoid[, setdiff(colnames(lymphoid), genes)]
    rownames(meta) <- lymphoid$cell
    seurat_obj <- SeuratObject::CreateSeuratObject(counts = round(expm1(X)), meta.data = meta)
    seurat_obj <- SeuratObject::SetAssayData(seurat_obj, layer = 'data', new.data = X)
    emb <- as.matrix(meta[, c('sub_umap_1', 'sub_umap_2')])
    colnames(emb) <- c('umap_1', 'umap_2')
    seurat_obj[['umap']] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = 'umap_', assay = 'RNA')

    p_gene <- ZoomEmbedding(
        full = sim_cells, full_reduction = 'umap', full_color_by = 'cluster', full_colors = cluster_colors,
        subset = seurat_obj, subset_reduction = 'umap', subset_color_by = 'CD8A',
        plot_min = 0, plot_max = 'q99',
        titles = c('All cells', sprintf('Lymphoid cells, n = %s', scales::comma(ncol(seurat_obj)))),
        legend_titles = c('Cluster', 'CD8A')
    )

    cairo_pdf(file.path(out_dir, 'zoom_gene.pdf'), width = 13, height = 7.5)
    print(p_gene)
    invisible(dev.off())
}

cat('figures written to', out_dir, '\n')
