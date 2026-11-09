# small inputs built from the simulated table
genes <- c('CD3E', 'CD8A', 'NKG7', 'MS4A1', 'LYZ', 'CD14', 'COL1A1', 'PECAM1')

sim_small <- function(n = 600, seed = 1){
    set.seed(seed)
    sim_cells[sort(sample(nrow(sim_cells), n)), ]
}

# a Seurat object with the marker genes in the data layer and the umap as a reduction
sim_seurat <- function(df = sim_small()){
    X <- t(as.matrix(df[, genes]))
    colnames(X) <- df$cell
    X <- Matrix::Matrix(X, sparse = TRUE)
    meta <- df[, setdiff(colnames(df), genes)]
    rownames(meta) <- df$cell
    obj <- SeuratObject::CreateSeuratObject(counts = round(expm1(X)), meta.data = meta)
    obj <- SeuratObject::SetAssayData(obj, layer = 'data', new.data = X)
    emb <- as.matrix(df[, c('umap_1', 'umap_2')])
    rownames(emb) <- df$cell
    colnames(emb) <- c('umap_1', 'umap_2')
    obj[['umap']] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = 'umap_', assay = 'RNA')
    obj
}
