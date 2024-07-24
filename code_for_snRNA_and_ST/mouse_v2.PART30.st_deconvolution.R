####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Part <- 'PART30_ST_Deconvolution'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2019_atrial_latscko_crtsai/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scATACseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/stRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('local', Ver, Part, 'mouse', Project_dir = '2019_atrial_latscko_crtsai', Data_drive = 'bree')
pt.sizes <- c(3.2, 1.2, 1.3, 1.5, 1.7, 2)
pt.sizes_wh <- c(2, 1.5, 1.5, 1.5, 1.5, 1.8)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Batch 1 -- High TAM ####
st_hitam.srt <- readRDS(paste0('/Volumes/shire/project/2022_st_compendium/rdata/mouse_v2/individual/',
                               '01.2022_LatsCKO.srt.rds'))

#### Batch 2 -- Low TAM ####
st_lotam.srt <- readRDS(paste0('/Volumes/shire/project/2022_st_compendium/rdata/mouse_v2/individual/',
                               '01.2023_Csf1r_CTsai.srt.rds'))

st.srt <- merge(st_hitam.srt, st_lotam.srt, merge.dr = T, merge.data = T)
st.srt$Sample <- factor(st.srt$Sample, levels = c(
    "Control",
    "LatsCKO", 
    "Control_Veh",
    "Control_Inh",
    "LatsCKO_Veh",
    "LatsCKO_Inh"
))
st.srt <- st.srt[, st.srt$Zone != '']
st.srt$Zone <- factor(st.srt$Zone, levels = c('LA', 'RA', 'IAS', 'AoV', 'Ao', 'MV', 'PA', 'Fat', 
                                              'IVS', 'LV', 'RV', 'TV'))
Idents(st.srt) <- 'Zone'
DimPlotST(st.srt, group.by = 'Zone', pt.sizes = pt.sizes_wh, 
          cols = mycol_20s, ncol = 6, legend = 3)

st_atr.srt <- DropMetaLevels(st.srt[, st.srt$Zone %in% c('Ao', 'AoV', 'Fat', 'IAS', 'MV', 'PA', 'RA', 'LA')])
Idents(st_atr.srt) <- 'Zone'
st_atr.srt$Zone <- factor(st_atr.srt$Zone, levels = c('LA', 'RA', 'IAS', 'AoV', 'Ao', 'MV', 'PA', 'Fat'))
DimPlotST(st_atr.srt, group.by = 'Zone', pt.sizes = pt.sizes_wh, 
          cols = mycol_10, ncol = 6, legend = 3)

#### Full snRNA Reference data ####
sn.srt <- readRDS('integrated/PART19.consolidated.srt.rds')
sn.srt <- DietSeurat(sn.srt, assays = 'RNA')
sn.srt <- DropMetaLevels(sn.srt[, sn.srt$Non_ambiguous & 
                                    !sn.srt$Cell_state %in% c('EC4_Prol', 'MP3_Prol', 'DC', 'Schwann')])
Table(sn.srt$Cell_state, sn.srt$Cell_type)

sn.srt$Cell_type_decon <- as.vector(sn.srt$Cell_type)
sn.srt$Cell_type_decon[sn.srt$Cell_state == 'FB4_YapHi'] <- 'FB_YapHi'
sn.srt$Cell_type_decon[sn.srt$Cell_state %in% c('FB1', 'FB2', 'FB3')] <- 'FB_YapLo'
sn.srt$Cell_type_decon[sn.srt$Cell_state %in% c('MP1_Recruit')] <- 'MP_Recruit'
sn.srt$Cell_type_decon[sn.srt$Cell_state %in% c('MP2_Resident')] <- 'MP_Resident'
sn.srt$Cell_type_decon[sn.srt$Cell_type_decon == 'Immune'] <- 'Lymphoid'
sn.srt$Cell_type_decon <- factor(sn.srt$Cell_type_decon, levels = c(
    'CM',
    'FB_YapLo',
    'FB_YapHi',
    'EpiC',
    'EC',
    'Mural',
    'MP_Recruit',
    'MP_Resident',
    'Lymphoid',
    'Adipo'
))
Idents(sn.srt) <- 'Cell_type_decon'

sn_dns.srt <- DropMetaLevels(sn.srt[, unlist(DownsampleByMeta(sn.srt, 
                                                              meta_var = 'Cell_type_decon',
                                                              n = 100,
                                                              down_to_min_group = F, 
                                                              random = T))]) ## Downsample for efficiency
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Deconvolve Atrial Spots with SPOTlight  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
st <- DietSeurat(st_atr.srt, assays = 'SCT')
st <- RenameAssays(st, 'SCT' = 'RNA')
st.sce <- as.SingleCellExperiment(st, assay = 'RNA')

p1 <- DimPlotST(st_atr.srt, group.by = 'Zone', pt.sizes = c(2, 1.5, 1.5, 1.5, 1.5, 1.8), 
          cols = mycol_10, ncol = 6, legend = 3)
p2 <- SpatialDimPlot(st_atr.srt, group.by = 'Zone') &
    scale_fill_manual(values = mycol_10)
PlotPDF('1.0.st_dim.atrial_zone', 18, 6)
print(p1/p2)
dev.off()


## Decon Cell Type ####
Table(sn_dns.srt$Cell_type_decon)
mk <- FindAllMarkers(sn_dns.srt, logfc.threshold = 1, min.pct = 0.3, only.pos = T, return.thresh = 0.001)
mk_ct <- mk[mk$avg_log2FC > 2 & mk$p_val_adj < 1e-10,]
mk_ct <- mk_ct[!duplicated(mk_ct$gene), ]
Table(mk_ct$cluster)
mk_ct <- mk_ct |> group_by(cluster) |> top_n(n = 50, wt = -p_val_adj) ## top 50 genes as marker
mk_ct2 <- mk_ct |> group_by(cluster) |> top_n(n = 10, wt = -p_val_adj) ## top 10 genes as marker
mk_ct$neg_log10_p_adj <- -log10(mk_ct$p_val_adj)
sc.sce <- as.SingleCellExperiment(DietSeurat(sn_dns.srt, assays = 'RNA'), assay = 'RNA')

res_cell_type <- SPOTlight(
    x = sc.sce,
    y = st.sce,
    groups = as.character(sc.sce$Cell_type_decon),
    mgs = as.data.frame(mk_ct),
    hvg = NULL,
    weight_id = "neg_log10_p_adj",
    group_id = "cluster",
    gene_id = "gene",
    min_prop = 0,
    assay_sc = 'RNA',
    assay_sp = 'RNA',
    slot_sc = 'counts',
    slot_sp = 'counts')
## 0.18min

mod_ct <- res_cell_type$NMF
decon_ct <- res_cell_type$mat

p <- plotTopicProfiles(
    x = mod_ct,
    y = sc.sce$Cell_type_decon,
    facet = F,
    min_prop = 0.01,
    ncol = 1) +
    theme(aspect.ratio = 1)
p
PlotPDF('1.1.dot.cell_type_topic_profile_plts', 4, 4)
print(p)
dev.off()

st_atr.srt <- AddModuleScore2(st_atr.srt, features = split(mk_ct2$gene, mk_ct2$cluster),
                              names = paste0('CT_Score_', levels(mk_ct2$cluster)), return_z = T)
p <- FeaturePlotST(st_atr.srt, 
                   features = paste0('CT_Score_', levels(mk_ct$cluster)),
                   minvals = rep(-1, 10),
                   maxvals = rep(3, 10),
                   pt.sizes = pt.sizes*0.3,
                   ncol = 6) &
    scale_fill_distiller(palette = 'Spectral', limits = c(-1, 3)) &
    theme(aspect.ratio = 0.5)
PlotPDF('1.2.spatial.cell_type_marker_score', 15, 20)
print(p)
dev.off()


identical(rownames(decon_ct), Cells(st_atr.srt))
decon_df <- data.frame(decon_ct)
colnames(decon_df) <- paste0('SPL_', colnames(decon_df))
tmp.srt <- AddMetaData(st_atr.srt, metadata = decon_df)
p <- FeaturePlotST(tmp.srt,
                   features = colnames(decon_df),
                   minvals = rep(0, 10),
                   maxvals = rep(0.3, 10),
                   pt.sizes = pt.sizes*0.3,
                   ncol = 6) &
    scale_fill_distiller(palette = 'Spectral', limits = c(0, 0.3)) &
    theme(aspect.ratio = 0.5)
PlotPDF('1.3.spatial.spl_deconv_celltype', 15, 20)
print(p)
dev.off()

pcc <- matrix(NA, 10, 10)
rownames(pcc) <- str_replace(paste0('SPL_', levels(mk_ct$cluster)), ' ', '.')
colnames(pcc) <- paste0('CT_Score_', levels(mk_ct$cluster))
for(i in 1:10){
    name_i <- rownames(pcc)[i]
    for(j in 1:10){
        name_j <- colnames(pcc)[j]
        pcc[i, j] <- cor(tmp.srt@meta.data[, name_i], st_atr.srt@meta.data[, name_j])
        if(pcc[i, j] < 0){pcc[i, j] <- 0}
    }
}
p <- ggplot(as.data.frame(melt(pcc)))+
    geom_tile(aes(x = Var1, y = Var2, fill = value)) +
    scale_fill_viridis_c() +
    theme_classic() +
    theme(aspect.ratio = 1) +
    RotatedAxis()
p
PlotPDF('1.4.heat.marker_score_decon_pcc', 10, 10)
p
dev.off()

## This part is for correcting a bug in SPOTLIGHTv1.0 as cell types were incorrectly named
decon_df2 <- decon_df
decon_df2$SPL_CM <- decon_df$SPL_FB_YapHi
decon_df2$SPL_FB_YapLo <- decon_df$SPL_CM
decon_df2$SPL_FB_YapHi <- decon_df$SPL_FB_YapLo
decon_df2$SPL_EpiC <- decon_df$SPL_EC
decon_df2$SPL_EC <- decon_df$SPL_Mural
decon_df2$SPL_Mural <- decon_df$SPL_EpiC
decon_df2$SPL_Resd <- res_cell_type$res_ss
st_atr.srt <- AddMetaData(st_atr.srt, metadata = decon_df2)

pcc <- matrix(NA, 10, 10)
rownames(pcc) <- str_replace(paste0('SPL_', levels(mk_ct$cluster)), ' ', '.')
colnames(pcc) <- paste0('CT_Score_', levels(mk_ct$cluster))
for(i in 1:10){
    name_i <- rownames(pcc)[i]
    for(j in 1:10){
        name_j <- colnames(pcc)[j]
        pcc[i, j] <- cor(st_atr.srt@meta.data[, name_i], st_atr.srt@meta.data[, name_j])
        if(pcc[i, j] < 0){pcc[i, j] <- 0}
    }
}
p <- ggplot(as.data.frame(melt(pcc)))+
    geom_tile(aes(x = Var1, y = Var2, fill = value)) +
    scale_fill_viridis_c() +
    theme_classic() +
    theme(aspect.ratio = 1) +
    RotatedAxis()
p
PlotPDF('1.5.heat.corrected_marker_score_decon_pcc', 10, 10)
p
dev.off()

p <- FeaturePlotST(st_atr.srt,
                   features = c(paste0('SPL_', levels(sc.sce$Cell_type_decon)), 'SPL_Resd'),
                   minvals = rep(0, 10),
                   maxvals = rep(0.4, 10),
                   pt.sizes = pt.sizes*0.3,
                   ncol = 6) &
    scale_fill_distiller(palette = 'Spectral', limits = c(0, 0.3)) &
    theme(aspect.ratio = 0.5)
PlotPDF('1.6.spatial.corrected_spl_deconv_celltype', 15, 20)
print(p)
dev.off()

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(res_cell_type, 'analysis/PART30.st_decon_by_cell_type.splotlight.rds')
saveRDS(sn_dns.srt, 'analysis/PART30.downsampled_sc_for_decon_by_celltype.srt.rds')
WriteCSV(mk_ct, 'PART30.top_50_markers_for_deconvolution')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Return Results to Whole Heart ST ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meta_a <- st_atr.srt@meta.data[, c('CT_Score_CM',
                                   'CT_Score_FB_YapLo',
                                   'CT_Score_FB_YapHi',
                                   'CT_Score_EpiC',
                                   'CT_Score_EC',
                                   'CT_Score_Mural',
                                   'CT_Score_MP_Recruit',
                                   'CT_Score_MP_Resident',
                                   'CT_Score_Lymphoid',
                                   'CT_Score_Adipo',
                                   'SPL_Adipo',
                                   'SPL_CM',
                                   'SPL_EC',
                                   'SPL_EpiC',
                                   'SPL_FB_YapHi',
                                   'SPL_FB_YapLo',
                                   'SPL_Lymphoid',
                                   'SPL_MP_Recruit',
                                   'SPL_MP_Resident',
                                   'SPL_Mural',
                                   'SPL_Resd')]
st.srt <- AddMetaData(st.srt, metadata = meta_a)

p <- FeaturePlotST(st.srt,
                   features = c('SPL_Adipo',
                                'SPL_CM',
                                'SPL_EC',
                                'SPL_EpiC',
                                'SPL_FB_YapHi',
                                'SPL_FB_YapLo',
                                'SPL_Lymphoid',
                                'SPL_MP_Recruit',
                                'SPL_MP_Resident',
                                'SPL_Mural',
                                'SPL_Resd'),
                   minvals = rep(0, 10),
                   maxvals = rep(0.4, 10),
                   pt.sizes = rep(0.2, 6),
                   ncol = 6, 
                   asp = 1) &
    scale_fill_distiller(palette = 'Spectral', na.value = 'grey90')
PlotPDF('2.1.spatial.corrected_spl_deconv_celltype_whole_heart', 15, 23)
print(p)
dev.off()

st.srt <- JoinLayers(st.srt, assay = 'ST')
st.srt <- NormalizeData(st.srt, assay = 'ST')
st.srt <- FindVariableFeatures(st.srt, assay = 'ST')
st_atr.srt <- JoinLayers(st_atr.srt, assay = 'ST')
st_atr.srt <- NormalizeData(st_atr.srt, assay = 'ST')
st_atr.srt <- FindVariableFeatures(st_atr.srt, assay = 'ST')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(st.srt, 'integrated/PART30.whole_heart_st_post_decon.srt.rds')
saveRDS(st_atr.srt, 'integrated/PART30.atrial_st_post_decon.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
