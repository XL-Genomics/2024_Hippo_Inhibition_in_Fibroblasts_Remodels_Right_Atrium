####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Part <- 'PART06_scVI'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2019_atrial_latscko_crtsai/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scATACseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('local', Ver, Part, 'mouse', Project_dir = '2019_atrial_latscko_crtsai', Data_drive = 'bree')

library('reticulate')
library('sceasy')
use_condaenv(condaenv = "scvi", conda = "/Users/Felix/Conda/anaconda3/bin/python") ## conda env on Gondor
sc <- import("scanpy")
scvi <- import("scvi", convert = FALSE)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
srt <- readRDS('integrated/PART03.merged.flt.srt.rds')
srt@meta.data <- readRDS('integrated/PART04.merged.dlt.srt_meta.rds')
DefaultAssay(srt) <- 'CBN'

genes_include <- RemoveRiboMito(rownames(srt), human_or_mouse = 'mouse', ribo = T, mito = T)

rna.srt <- DietSeurat(srt, layers = 'counts', assays = 'CBN', features = genes_include)
rna.srt <- RenameAssays(rna.srt, 'CBN', 'RNA')
rna.srt[['RNA']] <- as(object = rna.srt[['RNA']], Class = "Assay") ## switch to Seurat v3 assay
adata <- convertFormat(rna.srt, 
                       from = "seurat",
                       to = "anndata", 
                       main_layer = "counts", 
                       assay = "RNA", 
                       drop_single_values = T)
print(adata) 
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Train scVI Model  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
scvi$model$SCVI$setup_anndata(adata, 
                              batch_key = 'Sample_id',
                              continuous_covariate_keys = c('nCount_CBN', 'nFeature_CBN'))
svi <- scvi$model$SCVI(adata)
PYTORCH_ENABLE_MPS_FALLBACK=1
svi$train(use_gpu = F) ## real training time 4:38:31

# get the latent representation
latent = svi$get_latent_representation()
latent <- as.matrix(latent)
rownames(latent) = colnames(rna.srt)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(latent, 'analysis/PART06.scvi_latent.mat.rds')
svi$save(dir_path = 'analysis/PART06.trained_scvi.model')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Reduce latent dimensions in Seurat  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# put it back in our original Seurat object
ndims <- ncol(latent) ## 10
srt[["scVI"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = "CBN")

# Find clusters, then run UMAP, and visualize
DefaultAssay(srt) <- 'CBN'
srt <- FindNeighbors(srt, reduction = "scVI", dims = 1:ndims)
srt <- FindClusters(srt, resolution = seq(0.1, 1, 0.1))
srt <- RunUMAP(srt, reduction = "scVI", dims = 1:ndims, reduction.name = 'scVI_umap', reduction.key = 'scVIumap_')

PlotPDF('1.umap.cluster_scvi', 20, 20)
DimPlot2(srt, group.by = paste0('CBN_snn_res.', seq(0.1, 1, 0.1)), reduction = 'scVI_umap', pt.size = 1, raster = T,
         label = T, cols = c(mycol_40, 'black')) &
    NoLegend()
DimPlot2(srt, group.by = 'Group2', reduction = 'scVI_umap', label = F, cols = mycol_20, pt.size = 0.1, raster = F)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(srt, 'integrated/PART06.scvi_integrated.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
