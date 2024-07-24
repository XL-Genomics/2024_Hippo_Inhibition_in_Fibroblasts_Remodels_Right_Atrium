####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Part <- 'PART19_Consolidation'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2019_atrial_latscko_crtsai/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scATACseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('local', Ver, Part, 'mouse', Project_dir = '2019_atrial_latscko_crtsai', Data_drive = 'bree')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
srt <- readRDS('integrated/PART08.rna_integrated.srt.rds')
meta <- readRDS('integrated/PART10.annotated.srt_meta.rds')
srt <- AddMetaData(srt, metadata = meta)

srt$Cell_state <- NA
srt$Cell_type <- as.vector(srt$Cell_type)

Table(srt$Cell_state, srt$Cell_type)

srt$Cell_type <- factor(srt$Cell_type, 
                        levels = c('CM', 'FB', 'EpiC', 'EC', 'Mural', 'Immune', 'Schwann', 'Adipo', 'Doublet'))


## This is where the cell type-specific RNA umaps (with ambiguous cells) are stored
srt@reductions$sub_RNA_umap <- srt@reductions$umap
srt@reductions$sub_RNA_umap@cell.embeddings[,1] <- NA
srt@reductions$sub_RNA_umap@cell.embeddings[,2] <- NA
srt@reductions$sub_RNA_umap@key <- 'subRNAUMAP_'
colnames(srt@reductions$sub_RNA_umap@cell.embeddings) <- c('subRNAUMAP_1', 'subRNAUMAP_2')

## This is where the cell type-specific RNA umaps (without ambiguous cells) are stored
srt@reductions$sub_clean_RNA_umap <- srt@reductions$umap
srt@reductions$sub_clean_RNA_umap@cell.embeddings[,1] <- NA
srt@reductions$sub_clean_RNA_umap@cell.embeddings[,2] <- NA
srt@reductions$sub_clean_RNA_umap@key <- 'subcleanRNAUMAP_'
colnames(srt@reductions$sub_clean_RNA_umap@cell.embeddings) <- c('subcleanRNAUMAP_1', 'subcleanRNAUMAP_2')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Consolidate each cell types  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ct <- c('immune', 'cm', 'fb', 'ec', 'epic', 'other')
for(i in 1:6) {
        message(ct[i])
        
        sub.srt <- readRDS(paste0('integrated/PART1', i, '.', ct[i], '_cells.srt.rds'))
        
        srt$Cell_state[Cells(sub.srt)] <- as.vector(sub.srt$Cell_state)
        print(Table(srt$Cell_state, srt$Cell_type))
        
        srt@reductions$sub_RNA_umap@cell.embeddings[Cells(sub.srt), ] <-
                sub.srt@reductions$sub_RNA_umap@cell.embeddings
        srt@reductions$sub_clean_RNA_umap@cell.embeddings[Cells(sub.srt), ] <-
                sub.srt@reductions$sub_clean_RNA_umap@cell.embeddings
}
srt$Cell_state <- str_replace(as.vector(srt$Cell_state), '-', '_')
srt$Cell_state <- factor(srt$Cell_state, 
                         levels = c(
                                 'CM1', 'CM2_SAN', 
                                 'FB1',  'FB2',  'FB3', 'FB4_YapHi',
                                 'EpiC1_RA', 'EpiC2_LA',
                                 'EC1_Vessel', 'EC2_LAEndoC',  'EC3_RAEndoC', 'EC4_Prol', 'EC5_LEC',
                                 'Peri', 'SMC',  
                                 'MP1_Recruit',  'MP2_Resident',  'MP3_Prol', 'DC',  'TC', 'BC',
                                 'Schwann', 'Adipo', 
                                 'Doublet'))
srt$Cell_state[is.na(srt$Cell_state)] <- srt$Cell_type[is.na(srt$Cell_state)]

Table(srt$Cell_state, srt$Cell_type)

## Visualize all UMAPs
data <- data.frame(X = srt@reductions$sub_RNA_umap@cell.embeddings[, 1],
                   Y = srt@reductions$sub_RNA_umap@cell.embeddings[, 2],
                   Cell_type = srt$Cell_type,
                   Cell_state = srt$Cell_state)
p1 <- ggplot(data, aes(x = X, y = Y, color = Cell_state)) + 
        geom_point(size = 0.1) + 
        scale_color_manual(values = mycol_30s) + 
        facet_wrap(~Cell_type) + 
        theme_classic() +
        NoAxes() +
        labs(caption = 'sub_RNA_umap') +
        theme(aspect.ratio = 1)

data <- data.frame(X = srt@reductions$sub_clean_RNA_umap@cell.embeddings[, 1],
                   Y = srt@reductions$sub_clean_RNA_umap@cell.embeddings[, 2],
                   Cell_type = srt$Cell_type,
                   Cell_state = srt$Cell_state)
p2 <- ggplot(data, aes(x = X, y = Y, color = Cell_state)) + 
        geom_point(size = 0.1) + 
        scale_color_manual(values = mycol_30s) + 
        facet_wrap(~Cell_type) + 
        theme_classic() +
        NoAxes() +
        labs(caption = 'sub_clean_RNA_umap') +
        theme(aspect.ratio = 1)

PlotPDF('0.umap.consolidated_sub_umaps', 10, 10)
p1
p2
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Organize annotation order  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <- as.data.frame(Table(srt$Cell_state, srt$Cell_type))
p <- ggplot(data, aes(x = Var1, y = Var2, fill = log10(Freq+1))) +
        geom_tile() +
        scale_y_discrete(limits = rev) +
        scale_fill_distiller(palette = 'RdYlBu') +
        theme_classic() +
        theme(aspect.ratio = LU(srt$Cell_type)/LU(srt$Cell_state), axis.ticks = element_blank()) +
        labs(x = '', y = '') +
        RotatedAxis()
p
PlotPDF('1.1.heat.celltype_vs_cellstate', 12, 5)
p
dev.off()

p1 <- DimPlot2(srt, group.by = 'Cell_type', cols = c(mycol_10[1:8], 'grey80'), 
               reduction = 'RNA_FastMNN_umap', label = T, pt.size = 0.1)
p2 <- DimPlot2(srt, group.by = 'Cell_state', cols = c(mycol_30s[1:23], 'grey80'), 
               reduction = 'RNA_FastMNN_umap', label = T, pt.size = 0.1)
p1 + p2
PlotPDF('1.2.umap.annotated', 18, 7)
p1 + p2
dev.off()

tmp.srt <- srt[, unlist(DownsampleByMeta(srt, meta_var = 'Group1', down_to_min_group = T, random = T))]
mtx <- as.data.frame(Table(tmp.srt$Group1, tmp.srt$Cell_state))
total <- as.vector(Table(tmp.srt$Cell_state))
names(total) <- names(Table(tmp.srt$Cell_state))
mtx$Freq <- mtx$Freq/total[mtx$Var2]
p <- ggplot(melt(mtx)) +
        geom_tile(aes(x = Var1, y = Var2, fill = value*100)) +
        scale_fill_viridis_c() +
        theme_classic() +
        labs(x = '', y = '', fill = '%') +
        theme(aspect.ratio = 24/14) +
        scale_y_discrete(limits = rev) +
        RotatedAxis()
p
PlotPDF('1.3.heat.cell_state_distribution_by_sample', 10, 10)
p
dev.off()

## No Cell State is likely originated from single sample or technical artifact
srt$Non_ambiguous <- T
srt$Non_ambiguous[srt$Cell_type %in% c('Doublet')] <- F
srt$Non_ambiguous[srt$Cell_state %in% c('Doublet')] <- F
Table(srt$Non_ambiguous, srt$Cell_type)
Table(srt$Non_ambiguous, srt$Cell_state)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Re-embed global RNA UMAP without ambiguous cells  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Color_cell_type <- mycol_10[1:8]
Color_cell_state <- c(mycol_30s[1:23], 'grey75')

clean.srt <- srt[, srt$Non_ambiguous]

clean.srt <- RunUMAP(clean.srt, reduction = 'mnn', dims = 1:50,
                     reduction.name = 'clean_RNA_umap', reduction.key = 'cRNAUMAP_',
                     n.epochs = 200, min.dist = 0.8, n.neighbors = 50)
DimPlot2(clean.srt, reduction = 'clean_RNA_umap', group.by = 'Cell_type',
         cols = Color_cell_type, label = T, raster = T) /
        DimPlot2(clean.srt, reduction = 'clean_RNA_umap', group.by = 'Cell_state',
                 cols = Color_cell_state, label = F, raster = T)

## Clean RNA umap without ambiguous cells
srt@reductions$clean_RNA_umap <- srt@reductions$umap
srt@reductions$clean_RNA_umap@cell.embeddings[,1] <- NA
srt@reductions$clean_RNA_umap@cell.embeddings[,2] <- NA
srt@reductions$clean_RNA_umap@cell.embeddings[Cells(clean.srt), 1:2] <-
        clean.srt@reductions$clean_RNA_umap@cell.embeddings[,1:2]
srt@reductions$clean_RNA_umap@key <- 'cRUMAP_'
colnames(srt@reductions$clean_RNA_umap@cell.embeddings) <- c('cRUMAP_1', 'cRUMAP_2')

## Full RNA umap with doublets cells
srt@reductions$full_RNA_umap <- srt@reductions$RNA_FastMNN_umap
srt@reductions$full_RNA_umap@key <- 'fRUMAP_'
colnames(srt@reductions$sub_clean_RNA_umap@cell.embeddings) <- c('fRUMAP_1', 'fRUMAP_2')

## Sub RNA umap for all cell types
srt@reductions$sub_full_RNA_umap <- srt@reductions$sub_RNA_umap
srt@reductions$sub_full_RNA_umap@key <- 'fSUBRUMAP_'
colnames(srt@reductions$sub_full_RNA_umap@cell.embeddings) <- c('fSUBRUMAP_1', 'fSUBRUMAP_2')

## Sub RNA umap for all non-ambiguous cell types
srt@reductions$sub_clean_RNA_umap <- srt@reductions$sub_clean_RNA_umap
srt@reductions$sub_clean_RNA_umap@key <- 'cSUBRUMAP_'
colnames(srt@reductions$sub_clean_RNA_umap@cell.embeddings) <- c('cSUBRUMAP_1', 'cSUBRUMAP_2')


plist <- list(
        DimPlot2(srt, reduction = 'full_RNA_umap', group.by = 'Cell_type',
                 cols = Color_cell_type, label = T, raster = TRUE, pt.size = 1),
        DimPlot2(srt, reduction = 'full_RNA_umap', group.by = 'Cell_state',
                 cols = Color_cell_state, label = T, raster = TRUE, pt.size = 1),
        DimPlot2(srt, reduction = 'clean_RNA_umap', group.by = 'Cell_type',
                 cols = Color_cell_type, label = T, raster = TRUE, pt.size = 1),
        DimPlot2(srt, reduction = 'clean_RNA_umap', group.by = 'Cell_state',
                 cols = Color_cell_state, label = T, raster = TRUE, pt.size = 1)
)
p <- wrap_plots(plist, ncol = 4)
PlotPDF('2.1.umap.annotation_various_umaps', 26, 6)
p
dev.off()


PlotPDF('2.2.umap.sub_annotation_split', 12, 16)
DimPlot2(srt, reduction = 'sub_full_RNA_umap', group.by = 'Cell_state', split.by = 'Cell_type',
         cols = Color_cell_state, label = T, raster = TRUE, pt.size = 3, ncol = 3)
DimPlot2(srt, reduction = 'sub_clean_RNA_umap', group.by = 'Cell_state', split.by = 'Cell_type',
         cols = Color_cell_state, label = T, raster = TRUE, pt.size = 3, ncol = 3)
dev.off()

## rename mnn dimension
srt@reductions$RNA_mnn <- srt@reductions$mnn
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Keep important dimensions  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
all_reductions <- srt@reductions
srt <- RenameAssays(srt, 'CBN', 'RNA')
srt <- DietSeurat(srt, 
                  layers = c('counts', 'data'),
                  assays = c('RNA', 'ATAC'),
                  dimreducs = c('full_RNA_umap', 'clean_RNA_umap',  'RNA_mnn',
                                'sub_full_RNA_umap', 'sub_clean_RNA_umap'))
srt@reductions <- srt@reductions[c('full_RNA_umap', 'clean_RNA_umap',  'RNA_mnn',
                                   'sub_full_RNA_umap', 'sub_clean_RNA_umap')]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(all_reductions, 'integrated/PART19.all_reductions.srt_dimreducs.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Editing Meta Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
saveRDS(srt@meta.data, 'integrated/PART19.all_meta.srt_meta.rds')

meta <- srt@meta.data
colnames(meta)[colnames(meta) == 'Genotype_l'] <- 'Genotype_long'
colnames(meta)[colnames(meta) == 'Genotype_s'] <- 'Genotype_short'
colnames(meta)[colnames(meta) == 'pct_mito_CBN'] <- 'pct_mito_RNA'
colnames(meta)[colnames(meta) == 'S.Score'] <- 'S_Score'
colnames(meta)[colnames(meta) == 'G2M.Score'] <- 'G2M_Score'
meta <- meta[, c(
        'Study',
        'Library',
        'Cell_barcode',
        'Sample_id',
        'Method',
        'Platform', 
        'Protocol',
        'Tissue', 
        'Enrichment',
        'Preparation',
        'Condition',
        'Genotype_long',
        'Genotype_short',
        'Sex', 
        'Age',
        'Replicate', 
        'nCount_ATAC',
        'nFeature_ATAC', 
        'nCount_RNA',
        'nFeature_RNA',
        'pct_mito_RNA',
        'Doublet_SC',
        'Doublet_SC_score',
        'S_Score',
        'G2M_Score',
        'Cell_type',
        'Cell_state',
        'Non_ambiguous',
        'Group1',
        'Group2',
        'Group3'
)]
srt@meta.data <- meta

## Check cell state distribution
p1 <- MappingHeatmap(srt, que_var = 'Cell_type', ref_var = 'Group1',
                     percentage = T, log10_scale = T, center = F,
                     que_order = levels(srt$Cell_type), ref_order = levels(srt$Group1), ref.disp.min = 0)
p2 <- MappingHeatmap(srt, que_var = 'Cell_state', ref_var = 'Group1',
                     percentage = T, log10_scale = T, center = F,
                     que_order = levels(srt$Cell_state), ref_order = levels(srt$Group1), ref.disp.min = 0)
p1/p2
PlotPDF('3.1.heat.cell_type_vs_donor_distribution', 7, 10)
print(p1/p2)
dev.off()

data1 <- p1$data |> group_by(Query) |> transmute(max(Value)) |> distinct()
data2 <- p2$data |> group_by(Query) |> transmute(max(Value)) |> distinct()
p1 <- ggplot(data1) +
        geom_bar(stat = 'identity', aes(x = Query, y = `max(Value)`, fill = `max(Value)` < 0.5)) +
        geom_hline(yintercept = 0.5)
p2 <- ggplot(data2) +
        geom_bar(stat = 'identity', aes(x = Query, y = `max(Value)`, fill = `max(Value)` < 0.5)) +
        geom_hline(yintercept = 0.5)
p <- p1/p2 &
        theme_classic() &
        RotatedAxis()
p
PlotPDF('3.2.bar.cell_type_vs_donor_distribution', 10, 10)
p
dev.off()


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DefaultAssay(srt) <- 'RNA'
saveRDS(srt, 'integrated/PART19.consolidated.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#### Test  ####
srt1 <- srt[, srt$Group1 %in% c(
        'Control_LA_1',
        'Control_LA_2',
        'LatsCKO_LA_1',
        'LatsCKO_LA_2',
        'Control_RA_1',
        'Control_RA_2',
        'LatsCKO_RA_1',
        'LatsCKO_RA_2'
)]

srt1 <- RunUMAP(srt1[, srt1$Non_ambiguous], dims = 1:50, reduction = 'RNA_mnn')

DimPlot2(srt1, reduction = 'umap', group.by = 'Cell_state', cols = mycol_30s)

srt2 <- srt[, srt$Group1 %in% c(
        'Control_RA_Veh',
        'LatsCKO_RA_Veh',
        'Control_RA_Inh_1',
        'LatsCKO_RA_Inh_1',
        'Control_RA_Inh_2',
        'LatsCKO_RA_Inh_2'
)]
srt2 <- RunUMAP(srt2[, srt2$Non_ambiguous], dims = 1:50, reduction = 'RNA_mnn')
