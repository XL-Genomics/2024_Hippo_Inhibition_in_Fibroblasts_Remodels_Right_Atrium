####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Part <- 'PART12_Annotation_CM'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2019_atrial_latscko_crtsai/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scATACseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('local', Ver, Part, 'mouse', Project_dir = '2019_atrial_latscko_crtsai', Data_drive = 'bree')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
srt <- readRDS('integrated/PART08.rna_integrated.srt.rds')
meta <- readRDS('integrated/PART10.annotated.srt_meta.rds')
srt <- AddMetaData(srt, metadata = meta)

sub.srt <- srt[, srt$Cell_type %in% c('CM')]
sub.srt@reductions$mnn@assay.used <- 'CBN'
DimPlot2(sub.srt, group.by = 'Group1', reduction = 'RNA_FastMNN_umap')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  CM cell annotation by WNN  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sub.srt <- RunUMAP(sub.srt, reduction = 'mnn', dims = 1:50, min.dist = 0.5,
                   reduction.name = 'sub_RNA_umap', reduction.key = 'subRNAUMAP_')

## Check previous annotation
sort(Table(sub.srt$Cell_type_prev))
sort(Table(sub.srt$Cell_state_prev))
sub.srt$tmp <- factor(sub.srt$Cell_state_prev, levels = c(levels(sub.srt$Cell_state_prev), 'Other'))
sub.srt$tmp[! sub.srt$Cell_type_prev %in% c('CM', 'New', 'Doublet')] <- 'Other'
sub.srt$tmp <- droplevels(sub.srt$tmp)
DimPlot2(sub.srt[, sub.srt$Cell_type_prev != 'New'], group.by = 'tmp', reduction = 'sub_RNA_umap', cols = mycol_10)

## Cluster by WNN
sub.srt <- FindNeighbors(sub.srt, reduction = 'mnn', dims = 1:50, assay = 'CBN')
sub.srt <- FindClusters(sub.srt, resolution = seq(0.1, 1, 0.1))
p <- list(
        DimPlot2(sub.srt, reduction = 'sub_RNA_umap', group.by = 'Group1', cols = mycol_20s),
        DimPlot2(sub.srt, reduction = 'sub_RNA_umap', group.by = 'Doublet_SC', cols = mycol_10),
        DimPlot2(sub.srt, reduction = 'sub_RNA_umap', group.by = 'tmp', cols = mycol_10),
        DimPlot2(sub.srt, reduction = 'sub_RNA_umap', group.by = 'CBN_snn_res.0.1', label = T, cols = mycol_20s),
        DimPlot2(sub.srt, reduction = 'sub_RNA_umap', group.by = 'CBN_snn_res.0.2', label = T, cols = mycol_20s),
        DimPlot2(sub.srt, reduction = 'sub_RNA_umap', group.by = 'CBN_snn_res.0.3', label = T, cols = mycol_20s),
        DimPlot2(sub.srt, reduction = 'sub_RNA_umap', group.by = 'CBN_snn_res.0.4', label = T, cols = mycol_20s),
        DimPlot2(sub.srt, reduction = 'sub_RNA_umap', group.by = 'CBN_snn_res.0.5', label = T, cols = mycol_20s),
        DimPlot2(sub.srt, reduction = 'sub_RNA_umap', group.by = 'CBN_snn_res.0.6', label = T, cols = mycol_20s),
        DimPlot2(sub.srt, reduction = 'sub_RNA_umap', group.by = 'CBN_snn_res.0.7', label = T, cols = mycol_20s),
        DimPlot2(sub.srt, reduction = 'sub_RNA_umap', group.by = 'CBN_snn_res.0.8', label = T, cols = mycol_20s),
        DimPlot2(sub.srt, reduction = 'sub_RNA_umap', group.by = 'CBN_snn_res.0.9', label = T, cols = mycol_30s),
        DimPlot2(sub.srt, reduction = 'sub_RNA_umap', group.by = 'CBN_snn_res.1', label = T, cols = mycol_30s)
)
p <- wrap_plots(p, nrow = 4)
PlotPDF('1.0.umap.clusters', 20, 20)
p
dev.off()

PlotPDF('1.1.feat.markers_lvl1', 20, 20)
FeaturePlot2(sub.srt, features = markers_lvl1, reduction = 'sub_RNA_umap', raster = F)
dev.off()

PlotPDF('1.2.feat.markers_prev', 10, 10)
FeaturePlot2(sub.srt, 
             features = c(
                     'Angpt1', 'Nav3', 'Gja1',
                     'Nppa', 'Clu', 'Ryr3',
                     'Hcn1', 'Unc5c', 'Nt5dc2',
                     'Cntn2', 'Gabrg3', 'Smoc2',
                     'Acta2', 'Nrg3', 'Myh11'),
             reduction = 'sub_RNA_umap', raster = F)
dev.off()

Idents(sub.srt) <- 'CBN_snn_res.0.4'
mk <- FindAllMarkers(sub.srt, only.pos = T, logfc.threshold = 0.75)
PlotPDF('1.3.heat.cluster_markers', 15, 15)
MarkerHeatmap(srt = sub.srt, marker.df = mk, disp.min = 0.5) 
dev.off()

## Re-annotate
sub.srt$Cell_state <- NA
sub.srt$Cell_state[sub.srt$CBN_snn_res.0.4 %in% c(0, 1, 2, 4)] <- 'CM1' 
sub.srt$Cell_state[sub.srt$CBN_snn_res.0.4 %in% c(6)] <- 'CM2-SAN'
sub.srt$Cell_state[sub.srt$CBN_snn_res.0.4 %in% c(5, 7, 8, 3)] <- 'Doublet'
sub.srt$Cell_state[sub.srt$Doublet_SC] <- 'Doublet'
sub.srt$Cell_state <- factor(sub.srt$Cell_state, levels = c('CM1', 'CM2-SAN', 'Doublet'))

PlotPDF('2.1.dim.sub_cluster_annotated', 5, 5)
DimPlot2(sub.srt, reduction = 'sub_RNA_umap', label = T, cols = mycol_10, group.by = 'Cell_state')
dev.off()

Idents(sub.srt) <- 'Cell_state'
mk <- FindAllMarkers(sub.srt, only.pos = T, logfc.threshold = 1)
PlotPDF('2.2.heat.cluster_markers_annotated', 15, 15)
MarkerHeatmap(srt = sub.srt, marker.df = mk, disp.min = 0.5)
dev.off()

PlotPDF('3.dim.sub_cluster_annotated_wnn', 4, 10)
DimPlot2(sub.srt, reduction = 'sub_RNA_umap', label = T, cols = mycol_10, group.by = 'Cell_state')
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Re-embed without ambiguous  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
tmp.srt <- sub.srt[, sub.srt$Cell_state != 'Doublet']

## Re-embed RNA
tmp.srt <- RunUMAP(tmp.srt, reduction = 'mnn', dims = 1:50, min.dist = 0.5,
                   reduction.name = 'sub_clean_RNA_umap', reduction.key = 'subcleanRNAUMAP_')
DimPlot2(tmp.srt, reduction = 'sub_clean_RNA_umap', group.by = 'Cell_state', cols = mycol_10,  label = T)

## Add dimensions back to full seurat
sub.srt@reductions$sub_clean_RNA_umap <- sub.srt@reductions$sub_RNA_umap
sub.srt@reductions$sub_clean_RNA_umap@cell.embeddings[, c(1,2)] <- NA
sub.srt@reductions$sub_clean_RNA_umap@cell.embeddings[Cells(tmp.srt), c(1,2)] <-
        tmp.srt@reductions$sub_clean_RNA_umap@cell.embeddings
colnames(sub.srt@reductions$sub_clean_RNA_umap@cell.embeddings) <-
        colnames(tmp.srt@reductions$sub_clean_RNA_umap@cell.embeddings)
sub.srt@reductions$sub_clean_RNA_umap@key <- tmp.srt@reductions$sub_clean_RNA_umap@key

PlotPDF('4.umap.non_ambigouse_annotated_reembed', 4, 12)
DimPlot2(sub.srt, reduction = 'sub_clean_RNA_umap', group.by = 'Cell_state', cols = mycol_10,  label = T)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Check Abundance  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p <- CountCellBarPlot(sub.srt, group.by = 'Group2', stack.by = 'Cell_state', cols = mycol_10, 
                      stack.subset = levels(sub.srt$Cell_state)[levels(sub.srt$Cell_state) != 'Doublet'])
p
PlotPDF('5.bar.cell_state_composition', 5, 5)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
saveRDS(sub.srt, 'integrated/PART12.cm_cells.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----

