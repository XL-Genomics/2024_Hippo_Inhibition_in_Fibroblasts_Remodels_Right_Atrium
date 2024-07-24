####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Part <- 'PART10_Annotation_Global'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2019_atrial_latscko_crtsai/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('local', Ver, Part, 'mouse', Project_dir = '2019_atrial_latscko_crtsai', Data_drive = 'bree')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
srt <- readRDS('integrated/PART08.rna_integrated.srt.rds')
DefaultAssay(srt) <- 'CBN'
srt.bkp <- srt
srt.bkp -> srt

prev.srt <- readRDS('../mouse_v1/integrated/PART19.consolidated.srt.rds')
O(colnames(srt), colnames(prev.srt))
srt$Cell_type_prev <- factor('New', levels = c(levels(prev.srt$Cell_type), 'New'))
srt$Cell_type_prev[colnames(prev.srt)] <- prev.srt$Cell_type
srt$Cell_state_prev <- factor('New', levels = c(levels(prev.srt$Cell_state), 'New'))
srt$Cell_state_prev[colnames(prev.srt)] <- prev.srt$Cell_state
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Annotate Cell Type using WNN ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p1 <- FeaturePlot2(srt, features = c(markers_lvl1, 'Doublet_SC_score'), ncol = 6, pt.size = 1, order = F, raster = T,
                   reduction = 'RNA_FastMNN_umap', min.cutoff = 'q2', max.cutoff = 'q98')
PlotPDF('1.1.feat.lvl1_markers', 20, 20)
p1
dev.off()

tmp.srt <- srt[, unlist(DownsampleByMeta(srt, meta_var = 'Group1', down_to_min_group = T, random = T))]
p2 <- DimPlot2(srt, reduction = 'RNA_FastMNN_umap', pt.size = 1, label = T, raster = T,
               group.by = c(paste0('SCT_snn_res.', seq(0.1, 0.9, 0.1)), 'Doublet_SC', 'Cell_type_prev'),
               cols = mycol_40s, ncol = 3) +
    DimPlot2(tmp.srt, reduction = 'RNA_FastMNN_umap', pt.size = 1, label = F, raster = T,
             group.by = 'Group1', cols = mycol_20)
PlotPDF('1.2.dim.global_view', 20, 20)
p2
dev.off()

Idents(srt) <- 'SCT_snn_res.0.2'
mk <- FindAllMarkers(srt, only.pos = T, logfc.threshold = 2)
mk <- mk[mk$p_val_adj < 1e-5, ]


MappingHeatmap(srt, que_var = 'SCT_snn_res.0.2', ref_var = 'Cell_type_prev')

srt$Cell_type <- NA
srt$Cell_type[srt$SCT_snn_res.0.2 %in% c(5)] <- 'CM'
srt$Cell_type[srt$SCT_snn_res.0.2 %in% c(1, 8, 4)] <- 'FB'
srt$Cell_type[srt$SCT_snn_res.0.2 %in% c(0, 6, 13)] <- 'EC'
srt$Cell_type[srt$SCT_snn_res.0.2 %in% c(2, 7, 12)] <- 'EpiC'
srt$Cell_type[srt$SCT_snn_res.0.2 %in% c(11)] <- 'Mural'
srt$Cell_type[srt$SCT_snn_res.0.2 %in% c(3, 9, 10)] <- 'Immune'
srt$Cell_type[srt$SCT_snn_res.0.2 %in% c(14)] <- 'Schwann'
srt$Cell_type[srt$SCT_snn_res.0.3 %in% c(17)] <- 'Adipo'
srt$Cell_type[srt$SCT_snn_res.0.8 %in% c(27)] <- 'FB'

MappingHeatmap(srt, 'Cell_type', 'Cell_type_prev')

srt$Cell_type <- droplevels(factor(srt$Cell_type,
                                   levels = c('CM', 'FB', 'EpiC', 'EC', 'Mural',
                                              'Immune', 'Schwann', 'Adipo')))

p4 <- DimPlot2(srt, group.by = 'Cell_type', cols = c(mycol_10[1:8], 'grey75'), label = T, pt.size = 0.1,
               reduction = 'RNA_FastMNN_umap')
p4
PlotPDF('1.3.umap.main_cell_type_annotated', 6, 5)
p4
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df <- srt@meta.data[, !grepl('^SCT_snn_res.', colnames(srt@meta.data))]
saveRDS(df, 'integrated/PART10.annotated.srt_meta.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## Test:
CountCellBarPlot(srt, group.by = 'Group1', stack.by = 'Cell_type')
