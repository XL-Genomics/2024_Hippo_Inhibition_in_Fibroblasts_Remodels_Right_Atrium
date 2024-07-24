####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Part <- 'PART08_RNA_Integration'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2019_atrial_latscko_crtsai/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scATACseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('local', Ver, Part, 'mouse', Project_dir = '2019_atrial_latscko_crtsai', Data_drive = 'bree')

plan("multisession", workers = 4)
options(future.globals.maxSize = 100*1024^3) # for 100 Gb RAM
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
srt <- readRDS('integrated/PART03.merged.flt.srt.rds')
srt@meta.data <- readRDS('integrated/PART04.merged.dlt.srt_meta.rds')
DefaultAssay(srt) <- 'CBN'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Harmony-batch correction RNA integration  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
genes_include <- RemoveRiboMito(rownames(srt), human_or_mouse = 'mouse', ribo = T, mito = T)

tmp.srt <- DietSeurat(srt, assays = 'CBN', features = genes_include)

tmp.srt <- FindVariableFeatures(tmp.srt)
tmp.srt <- ScaleData(tmp.srt, vars.to.regress = c('nCount_CBN', 'nFeature_CBN'))
tmp.srt <- RunPCA(tmp.srt)
tmp.srt <- RunUMAP(tmp.srt, reduction = 'pca', dims = 1:50, 
                   reduction.name = 'RNA_pc_umap', reduction.key = 'pcrnaumap_', 
                   min.dist = 0.5, n.neighbors = 50)
tmp.srt <- RunHarmony(tmp.srt, group.by.vars = 'Group1', assay.use = 'CBN', reduction = 'pca',  
                      kmeans_init_nstart=20, kmeans_init_iter_max=100,
                      project.dim = F, reduction.save = 'RNA_harmony')
tmp.srt <- RunUMAP(tmp.srt, reduction = 'RNA_harmony', dims = 1:50, 
                   reduction.name = 'RNA_harmony_umap', reduction.key = 'hmnrnaumap_', 
                   min.dist = 0.5, n.neighbors = 50)

srt@reductions$pca <- tmp.srt@reductions$pca
srt@reductions$RNA_harmony <- tmp.srt@reductions$RNA_harmony
srt@reductions$RNA_pc_umap <- tmp.srt@reductions$RNA_pc_umap
srt@reductions$RNA_harmony_umap <- tmp.srt@reductions$RNA_harmony_umap

p1 <- DimPlot2(srt, reduction = 'RNA_harmony_umap', pt.size = 0.1, group.by = 'Group1', cols = mycol_20) +
    labs(title = 'Harmony-Corrected RNA UMAP')
p2 <- DimPlot2(srt, reduction = 'umap', pt.size = 0.1, group.by = 'Group1', cols = mycol_20) +
    labs(title = 'Uncorrected Joint UMAP')
p1 + p2

PlotPDF('2.harmony_umap.rna', 12, 6)
print(p1 + p2)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  FastMNN-batch correction RNA integration  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp.srt <- DietSeurat(srt, assays = 'CBN', features = genes_include)
tmp.srt[['CBN']] <- as(tmp.srt[['CBN']], Class = 'Assay') ## Switch to v3 assay

srt.list <- SplitObject(tmp.srt, split.by = 'Group1')
for(i in 1:L(srt.list)){
    srt.list[[i]] <- SCTransform(object = srt.list[[i]], method ="glmGamPoi", assay = 'CBN',
                                 variable.features.n = 3000)
}

tmp.srt <- RunFastMNN(object.list = srt.list, assay = 'SCT', features = 3000)
tmp.srt <- RunUMAP(tmp.srt, reduction = 'mnn', dims = 1:50, 
                   reduction.name = 'RNA_FastMNN_umap', reduction.key = 'mmnrnaumap_',
                   min.dist = 0.5, n.neighbors = 50)
tmp.srt$Group1 <- factor(tmp.srt$Group1, levels = levels(srt$Group1))
p1 <- DimPlot2(srt, reduction = 'RNA_harmony_umap', pt.size = 0.1, group.by = 'Group1', cols = mycol_20) +
    labs(title = 'Harmony-Corrected RNA UMAP')
p2 <- DimPlot2(tmp.srt, reduction = 'RNA_FastMNN_umap', pt.size = 0.1, group.by = 'Group1', cols = mycol_20) +
    labs(title = 'FastMNN-Corrected RNA UMAP')
p1 + p2
PlotPDF('3.mnn_umap.rna', 12, 6)
print(p1 + p2)
dev.off()

srt@assays$SCT <- tmp.srt@assays$SCT
srt@reductions$mnn <- tmp.srt@reductions$mnn
srt@reductions$RNA_FastMNN_umap <- tmp.srt@reductions$RNA_FastMNN_umap
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Pick FastMNN for Embedding and Clustering  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DefaultAssay(srt) <- 'SCT'
srt <- FindNeighbors(srt, reduction = 'mnn', dims = 1:50)
srt <- FindClusters(srt, resolution = seq(0.1, 0.9, 0.1))
p <- DimPlot2(srt, reduction = 'RNA_FastMNN_umap', pt.size = 0.1, label = T,
              group.by = paste0('SCT_snn_res.', seq(0.1, 0.9, 0.1)), cols = mycol_40s, ncol = 3) +
    labs(caption = 'FastMNN-Corrected Joint UMAP')
p

PlotPDF('4.umap.fastmnn_cluster', 15, 13)
p
dev.off()

DefaultAssay(srt) <- 'CBN'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(srt, 'integrated/PART08.rna_integrated.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----