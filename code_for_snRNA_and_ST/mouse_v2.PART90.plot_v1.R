####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Part <- 'PART90_Plot_V1'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2019_atrial_latscko_crtsai/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scATACseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/stRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('local', Ver, Part, 'mouse', Project_dir = '2019_atrial_latscko_crtsai', Data_drive = 'bree')

Color_cell_type <- c(mycol_10[1:8], 'black')
Color_cell_state <- c(mycol_30s[1:24], 'grey70')
Color_sample <- mycol_50[c(4, 1, 9, 6)]
Color_sample_inh <- c( 'steelblue2', '#3C5488', '#FC8D62FF', 'plum3')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~ Library Info ####
lib_meta.df <- read.csv(paste0(Docu_dir, 'mouse_sample_meta.csv'))


##~~  Full RNA+ATAC data with ambiguous cells ####
full_amb.srt <- readRDS('integrated/PART19.consolidated.srt.rds')
Idents(full_amb.srt) <- 'Cell_type'

Table(full_amb.srt$Cell_state, full_amb.srt$Cell_type)

##~~  Clean RNA+ATAC data without ambiguous cells ####
Table(full_amb.srt$Cell_state, full_amb.srt$Non_ambiguous)
full.srt <- DropMetaLevels(full_amb.srt[, full_amb.srt$Non_ambiguous])
Table(full.srt$Cell_state, full.srt$Cell_type)


##~~  Clean RNA data of LA+RA Cont+LatsCKO ####
# lr.srt <- DropMetaLevels(full.srt[, full.srt$Study == '2023_Latscko_CTsai'])
# lr.srt <- RunUMAP(lr.srt, dims = 1:50, reduction = 'RNA_mnn')
# saveRDS(lr.srt, 'integrated/PART90.lr_snrna_for_plotting.srt.rds')
lr.srt <- readRDS('integrated/PART90.lr_snrna_for_plotting.srt.rds')
lr_fb.srt <- DropMetaLevels(lr.srt[, lr.srt$Cell_type == 'FB'])
lr_mp.srt <- DropMetaLevels(lr.srt[, lr.srt$Cell_type == 'Immune'])


##~~  Clean RNA data of RA Csf1r Inhibitor Cont+LatsCKO ####
# inh.srt <- DropMetaLevels(full.srt[, full.srt$Study == '2023_Csf1r_CTsai'])
# inh.srt <- RunUMAP(inh.srt, dims = 1:50, reduction = 'RNA_mnn')
# saveRDS(inh.srt, 'integrated/PART90.csf1rinh_snrna_for_plotting.srt.rds')
inh.srt <- readRDS('integrated/PART90.csf1rinh_snrna_for_plotting.srt.rds')
inh_fb.srt <- DropMetaLevels(inh.srt[, inh.srt$Cell_type == 'FB'])
inh_mp.srt <- DropMetaLevels(inh.srt[, inh.srt$Cell_type == 'Immune'])


##~~  Clean RNA+ATAC cell type and state markers ####
markers <- readRDS('analysis/PART20.markers.srt_misc.rds')


##~~  Gene sets ####
gl1.mu <- Yap_target$Yap_target_cf_cutrun
gl1.hs <- ConvertGeneSpecies(gl1.mu, 'mouse', 'human')
gl1.hs <- revalue(gl1.hs, replace = c(CCN1 = 'CYR61', CCN5 = 'WISP2', PWWP3B = 'MUM1L1'))

gl2.hs <- read.csv('external/HALLMARK_GLYCOLYSIS.v2023.2.Hs.grp', header = T, comment.char = '#')[, 1]
gl2.mu <- ConvertGeneSpecies(gl2.hs, 'human', 'mouse')

gl3.hs <- read.csv('external/KEGG_GLYCOLYSIS_GLUCONEOGENESIS.v2023.2.Hs.grp', header = T, comment.char = '#')[, 1]
gl3.mu <- ConvertGeneSpecies(gl3.hs, 'human', 'mouse')

gl4.hs <- read.csv('external/REACTOME_GLYCOLYSIS.v2023.2.Hs.grp', header = T, comment.char = '#')[, 1]
gl4.mu <- ConvertGeneSpecies(gl4.hs, 'human', 'mouse')

gl5.mu <- U(read_excel('external/Glycolysis_genes.xlsx')$Symbol) ## Glycolysis GO
gl5.hs <- ConvertGeneSpecies(gl5.mu, 'mouse', 'human')

gl6.mu <- U(read_excel('external/Glycolysis_genes.xlsx', sheet = 2)$Gene) ## Glycolysis GO Overlap with Bulk CKO Up
gl6.hs <- ConvertGeneSpecies(gl6.mu, 'mouse', 'human')

gl7.hs <- U(read_excel('external/Table1_Chondrocyte_Osteoblast_genes.xlsx')$Chondrocyte) ## Chondrocyte markers
gl7.hs <- gl7.hs[!is.na(gl7.hs)]
gl7.mu <- ConvertGeneSpecies(gl7.hs, 'human', 'mouse')

gl8.hs <- U(read_excel('external/Table1_Chondrocyte_Osteoblast_genes.xlsx')$Osteoblast) ## Osteoblast markers
gl8.hs <- gl8.hs[!is.na(gl8.hs)]
gl8.mu <- ConvertGeneSpecies(gl8.hs, 'human', 'mouse')

gl9.hs <- union(gl8.hs, gl7.hs) ## Combined OCP markers
gl9.mu <- ConvertGeneSpecies(gl9.hs, 'human', 'mouse')


##~~  Clean Teichmann Human Heart ####
# human.srt <- readRDS('/Volumes/shire/data/resource_scrna/STeichmann_heart_atlas/2020_v1/raw.srt.rds')
# human_cf.srt <- DropMetaLevels(human.srt[, human.srt$cell_type == 'Fibroblast' &
#                                              human.srt$cell_source %in% c('Harvard-Nuclei', 'Sanger-Nuclei') &
#                                              human.srt$region %in% c('LA', 'RA', 'LV', 'RV') &
#                                              human.srt$donor != 'D1'])  ## D1 doesn't have RA sample
# human_cf.srt$region <- factor(human_cf.srt$region, levels = c('LV', 'RV', 'LA', 'RA'))
# human_cf.srt <- AddModuleScore2(human_cf.srt,
#                                 features = list(gl1.hs,
#                                                 gl2.hs,
#                                                 gl3.hs,
#                                                 gl4.hs,
#                                                 gl5.hs,
#                                                 gl6.hs,
#                                                 U(c(gl2.hs, gl3.hs, gl4.hs)),
#                                                 gl9.hs
#                                                 ),
#                                 names = c('Yap_CF_Target',  ## gl1.hs
#                                           'Hallmark_Glycolysis',  ## gl2.hs
#                                           'KEGG_Glycolysis_Glucogenesis',  ## gl3.hs
#                                           'Reactome_Glycolysis',  ## gl4.hs
#                                           'GO_Glycolysis',  ## gl5.hs
#                                           'GO_Glycolysis_Up_in_LatsCKO',  ## gl6.hs
#                                           'Union_Glycolysis',  ## U(c(gl2.hs, gl3.hs, gl4.hs, gl5.hs)
#                                           'OCP_Marker'  ## gl9.hs
#                                           ),
#                                 return_z = T)
# saveRDS(human_cf.srt, 'integrated/PART90.human_cf.srt.rds')
human_cf.srt <- readRDS('integrated/PART90.human_cf.srt.rds')


##~~  ST Data - Whole Heart ####
st.srt <- readRDS('integrated/PART30.whole_heart_st_post_decon.srt.rds')
Idents(st.srt) <- 'Zone'
pt.size.wh <- c(1.15, 1, 1, 1, 1.1, 1)


##~~  ST Data - Atrial Only ####
st_atr.srt <- readRDS('integrated/PART30.atrial_st_post_decon.srt.rds')
Idents(st_atr.srt) <- 'Zone'
pt.size.atr <- c(3.2, 1.2, 1.3, 1.5, 1.7, 2)*0.3
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  For Publications 2024-02  ####
####  Figure 1  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~ Plot 1 Box plot - Human RA vs LA YAP activity and Glycolysis activity ####
data <- human_cf.srt@meta.data[, c('Yap_CF_Target',
                                   'Hallmark_Glycolysis',
                                   'KEGG_Glycolysis_Glucogenesis',
                                   'Reactome_Glycolysis',
                                   'GO_Glycolysis',
                                   'GO_Glycolysis_Up_in_LatsCKO',
                                   'Union_Glycolysis',
                                   'OCP_Marker',
                                   'region',
                                   'donor')]
plist <- list()
for(x in 1:8) {
    lst <- split(data[, x], data$region)
    sub_data <- NULL
    for(i in 1:L(lst)){
        qtl <- quantile(lst[[i]], probs = c(0.1, 0.9))
        lst[[i]] <- lst[[i]][lst[[i]] < qtl[2] & lst[[i]] > qtl[1]]
        sub_data <- rbind(sub_data, data.frame(Region = rep(names(lst[i]), L(lst[[i]])), Score = lst[[i]]))
    }
    sub_data$Region <- factor(sub_data$Region, levels = c('LV', 'RV', 'LA', 'RA'))
    plist[[x]] <- ggplot(sub_data) +
        geom_boxplot(aes(x = Region, y = Score, fill = Region), outlier.shape = NA) +
        scale_fill_manual(values = c(rep('grey80', 3), mycol_10[2])) +
        theme_classic() +
        labs(title = colnames(data)[x], y = 'Z-score') +
        theme(aspect.ratio = 2) +
        NoLegend()
}
wrap_plots(plist[1:8], ncol = 4)
x <- split(data$Yap_CF_Target, data$region)
pval1 <- t.test(x$LA, x$RA)$p.value
pval1 # < 1e-20
plist[[1]] <- plist[[1]] +
    labs(caption = 'T test p < 1e-20')

y <- split(data$Reactome_Glycolysis, data$region) ## Use Reactome_Glycolysis
pval2 <- t.test(y$LA, y$RA)$p.value
pval2 # < 1e-20
plist[[4]] <- plist[[4]] +
    labs(caption = 'T test p < 1e-20')
p1 <- wrap_plots(plist[c(1,4)], ncol = 2)
p1
PlotPDF('1.01.box.human_ra_yap__glycolysis_score', 5, 5)
p1
dev.off()



##~~ Plot 2 Baton plot - Human RA vs LA YAP activity and Glycolysis activity ####
data <- human_cf.srt@meta.data[, c('Yap_CF_Target',
                                   'Reactome_Glycolysis',
                                   'region',
                                   'donor')]
data2 <- data |> group_by(donor, region) |>  mutate(Mean_Yap_CF_Target = mean(Yap_CF_Target),
                                                    Mean_Glycolysis = mean(Reactome_Glycolysis))
data2 <- U(data2[, c('donor', 'region', 'Mean_Yap_CF_Target', 'Mean_Glycolysis')])
data3 <- data2[data2$region %in% c('RA', 'LA') & data2$donor != 'D1', ]

p2.1 <- ggplot(data3, aes(x = region, y = Mean_Yap_CF_Target)) +
    geom_point(size = 2) +
    geom_line(aes(group = donor), linewidth = 1, alpha = 0.5) +
    scale_color_manual(values = mycol_20) +
    #scale_y_continuous(limits = c(-1.4, 1.4)) +
    labs(x = '', y = 'Mean Z-score', color = 'Donor', title = 'Yap Target in CF',
         caption = paste0('Wilcoxon signed-rank test p < ',
                          round(wilcox.test(split(data3$Mean_Yap_CF_Target, data3$region)[['LA']],
                                            split(data3$Mean_Yap_CF_Target, data3$region)[['RA']],
                                            paired = T)$p.value,
                                5)))+
    theme_classic() +
    theme(aspect.ratio = 2) +
    NoLegend()
p2.2 <- ggplot(data3, aes(x = region, y = Mean_Glycolysis)) +
    geom_point(size = 2) +
    geom_line(aes(group = donor), linewidth = 1, alpha = 0.5) +
    #scale_y_continuous(limits = c(-1.5, 0.5)) +
    labs(x = '', y = 'Mean Z-score', color = 'Donor', title = 'HALLMARK GLYCOLYSIS',
         caption = paste0('Wilcoxon signed-rank test p < ',
                          round(wilcox.test(split(data3$Mean_Glycolysis, data3$region)[['LA']],
                                            split(data3$Mean_Glycolysis, data3$region)[['RA']],
                                            paired = T)$p.value,
                                5)))+
    theme_classic() +
    theme(aspect.ratio = 2)
p2.1 + p2.2
PlotPDF('1.02.baton.human_ra_yap_score_glycolysis', 6, 6)
p2.1 + p2.2
dev.off()



##~~ Plot 3 UMAP LA+RA snRNA Global  ####
p3 <- DimPlot2(lr.srt, group.by = 'Cell_type', cols = Color_cell_type, pt.size = 0.1, alpha = 0.3, reduction = 'umap')
p3
PlotPDF('1.03.umap.la_ra_global_cell_type', 5, 5)
p3
dev.off()



##~~ Plot 4 Heatmap LA+RA snRNA Cell Type Markers  ####
df <- markers$Cell_type_rna_marker
df <- df[df$p_val_adj == 0 & df$pct.1 > 0.3, ]
Table(df$cluster)
p4 <- MarkerHeatmap(lr.srt, marker.df = df, group.cols = Color_cell_type, top = 15)
p4
PlotPDF('1.04.heat.la_ra_cell_type_markers', 6, 6)
p4
dev.off()



##~~ Plot 5 UMAP LA+RA snRNA Sub UMAP  ####
data <- data.frame(X = lr.srt@reductions$sub_clean_RNA_umap@cell.embeddings[, 1],
                   Y = lr.srt@reductions$sub_clean_RNA_umap@cell.embeddings[, 2],
                   Cell_type = lr.srt$Cell_type,
                   Cell_state = lr.srt$Cell_state)
data <- data[! data$Cell_type %in% c('Schwann', 'Adipo'), ]
p5 <- ggplot(data, aes(x = X, y = Y, color = Cell_state)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c(mycol_20s, 'grey')) +
    facet_wrap(~Cell_type, scales = 'free', nrow = 1) +
    theme_classic() +
    NoAxes() +
    theme(aspect.ratio = 1)
p5
PlotPDF('1.05.umap.la_ra_sub_cell_state', 20, 3)
p5
dev.off()



##~~ Plot 6 Heatmap LA+RA snRNA Cell State Markers  ####
data <- markers$Cell_state_rna_marker
plist <- list()
for(i in 1:6){
    mk <- data[[i]]
    mk <- mk[df$p_val_adj == 0 & mk$pct.1 > 0.3, ]
    tmp.srt <- DropMetaLevels(lr.srt[, lr.srt$Cell_type == names(data)[i]])
    Idents(tmp.srt) <- 'Cell_state'
    plist[[i]] <- MarkerHeatmap(tmp.srt, marker.df = mk, top = 30, n_cells = 100, disp.min = 0.5, disp.max = 1.5) +
        theme(aspect.ratio = 1)
}
p6 <- wrap_plots(plist, nrow = 1)
p6
PlotPDF('1.06.heat.la_ra_cell_state_markers', 35, 5)
p6
dev.off()



##~~ Plot 7 Box Plot LA+RA snRNA Abundance  ####
mtx <- as.matrix(Table(lr.srt$Cell_type, lr.srt$Group1) + 1)
for(i in 1:ncol(mtx)){mtx[, i] <- mtx[, i]/sum(mtx[, i])}
data <- data.frame(LA_1 = mtx[, 3]/mtx[, 1],
                   LA_2 = mtx[, 4]/mtx[, 2],
                   RA_1 = mtx[, 7]/mtx[, 5],
                   RA_2 = mtx[, 8]/mtx[, 6])
data$Cell_type <- rownames(data)
data <- melt(data)
data$LR <- str_split(data$variable, '_', simplify = T)[, 1]
data <- data[data$Cell_type %in% c('CM', 'EC', 'EpiC', 'FB', 'Immune'), ]
data$Cell_type <- factor(data$Cell_type, levels = c('CM', 'FB', 'EpiC', 'EC', 'Immune'))
p7 <- ggplot(data, aes(x = Cell_type, y = value, group.by = LR)) +
    geom_point() +
    geom_boxplot(outlier.shape = NA, aes(fill = Cell_type)) +
    scale_fill_manual(values = Color_cell_type) +
    scale_y_continuous(breaks = seq(0, 10, 1))+
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'red3') +
    facet_wrap(~LR, nrow = 1) +
    RotatedAxis() +
    labs(x = '', y = 'Abundance FC vs Controls') +
    theme_classic() +
    theme(aspect.ratio = 2)
p7
PlotPDF('1.07.box.la_ra_cell_type_abundance', 5, 5)
p7
dev.off()



##~~ Plot 8 UMAP LA+RA snRNA 4 Group Global  ####
data <- data.frame(X = lr.srt@reductions$umap@cell.embeddings[, 1],
                   Y = lr.srt@reductions$umap@cell.embeddings[, 2],
                   Group2 = lr.srt$Group2)
p8 <- ggplot(data, aes(x = X, y = Y, color = Group2)) +
    geom_point(alpha = 0.5, size = 0.2) +
    scale_color_manual(values = Color_sample) +
    facet_wrap(~Group2, nrow = 1) +
    theme_classic() +
    NoAxes() +
    NoLegend() +
    theme(aspect.ratio = 1)
p8
PlotPDF('1.08.umap.la_ra_samples', 12, 3)
p8
dev.off()



##~~ Plot 9 Milo LasCKO vs Control in LA+RA snRNA Abundance Milo  ####
Table(lr.srt$Group1, lr.srt$Genotype_short)
# tmp.srt <- DropMetaLevels(lr.srt[, lr.srt$Tissue == 'LA'])
# tmp.srt$Group2 <- tmp.srt$Genotype_short
# milo_la <- RunMilo(tmp.srt, k = 20, d = 50, alpha = 0.05, prop = 0.2,
#                     umap_dim = 'umap', pc_dim = 'RNA_mnn')
# tmp.srt <- DropMetaLevels(lr.srt[, lr.srt$Tissue == 'RA'])
# tmp.srt$Group2 <- tmp.srt$Genotype_short
# milo_ra <- RunMilo(tmp.srt, k = 20, d = 50, alpha = 0.05, prop = 0.2,
#                     umap_dim = 'umap', pc_dim = 'RNA_mnn')
# ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# saveRDS(milo_la, 'analysis/PART90.global_la_sn_lats_vs_cont.milo_result.rds')
# saveRDS(milo_ra, 'analysis/PART90.global_ra_sn_lats_vs_cont.milo_result.rds')
# ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
milo_la <- readRDS('analysis/PART90.global_la_sn_lats_vs_cont.milo_result.rds')
milo_ra <- readRDS('analysis/PART90.global_ra_sn_lats_vs_cont.milo_result.rds')
cutoff <- 8
milo_la[[2]]$logFC[milo_la[[2]]$logFC > cutoff] <- cutoff
milo_la[[2]]$logFC[milo_la[[2]]$logFC < -cutoff] <- -cutoff
milo_ra[[2]]$logFC[milo_ra[[2]]$logFC > cutoff] <- cutoff
milo_ra[[2]]$logFC[milo_ra[[2]]$logFC < -cutoff] <- -cutoff
p9.1 <- plotNhoodGraphDA.my(milo_la[[1]], milo_res = milo_la[[2]],
                            alpha = 0.1, edge_alpha = 0.1, edge_color = 'grey90') +
    labs(title = 'LA')
p9.2 <- plotNhoodGraphDA.my(milo_ra[[1]], milo_res = milo_ra[[2]],
                            alpha = 0.1, edge_alpha = 0.1, edge_color = 'grey90') +
    labs(title = 'RA')
p9 <- p9.1 + p9.2 &
    scale_fill_distiller(palette = 'RdBu', limits = c(-cutoff, cutoff)) &
    theme(aspect.ratio = 1,
          axis.line = element_line(color = 'black'))
p9
PlotPDF('1.09.milo.lats_vs_cko_split_by_la_ra', 10, 4)
p9
dev.off()



##~~ Plot 10 Milo LA vs RA in Control snRNA Abundance Milo  ####
# Table(lr.srt$Group1, lr.srt$Genotype_short)
# tmp.srt <- DropMetaLevels(lr.srt[, lr.srt$Genotype_short == 'Control'])
# tmp.srt$Group2 <- tmp.srt$Tissue
# milo_cont <- RunMilo(tmp.srt, k = 20, d = 50, alpha = 0.05, prop = 0.2,
#                      umap_dim = 'umap', pc_dim = 'RNA_mnn')
# ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# saveRDS(milo_cont, 'analysis/PART90.global_control_sn_ra_vs_la.milo_result.rds')
# ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
milo_cont <- readRDS('analysis/PART90.global_control_sn_ra_vs_la.milo_result.rds')
cutoff <- 6
milo_cont[[2]]$logFC[milo_cont[[2]]$logFC > cutoff] <- cutoff
milo_cont[[2]]$logFC[milo_cont[[2]]$logFC < -cutoff] <- -cutoff
p10.1 <- plotNhoodGraphDA.my(milo_cont[[1]], milo_res = milo_cont[[2]],
                            alpha = 0.1, edge_alpha = 0.1, edge_color = 'grey90') +
    labs(title = 'Control RA vs LA') +
    scale_fill_distiller(palette = 'RdBu', limits = c(-cutoff, cutoff)) +
    theme(aspect.ratio = 1,
          axis.line = element_line(color = 'black'))
p10.2 <- DimPlot2(tmp.srt, group.by = 'Cell_state', cols = Color_cell_state,
                  pt.size = 0.1, alpha = 0.3, reduction = 'umap')
p10.1 + p10.2
PlotPDF('1.10.milo.ra_vs_la_in_control', 14, 6)
p10.1 + p10.2
dev.off()



##~~ Plot 11 UMAP LA+RA snRNA Cell State UMAP  ####
p11 <- DimPlot2(lr.srt, group.by = 'Cell_state', cols = Color_cell_state,
                pt.size = 0.1, alpha = 0.3, reduction = 'umap')
p11
PlotPDF('1.11.umap.la_ra_global_cell_type', 7, 6)
p11
dev.off()


##~~ Plot 12 ST All ST Zone  ####
p12.1 <- DimPlotST(st.srt, group.by = 'Zone', pt.sizes = pt.size*1.2, cols = mycol_20s, ncol = 6, legend = 6) &
    theme(aspect.ratio = 1) &
    RestoreLegend()
p12.2 <-SpatialDimPlot(st.srt, alpha = 0) &
    theme(aspect.ratio = 1)  &
    NoLegend()
p12 <- p12.1/p12.2
PlotPDF('1.12.st.all_st_by_zone', 28, 8)
p12
dev.off()


##~~ Plot 13 ST Cell type deconvolution  ####
gl <- grep('^SPL_', colnames(st_atr.srt@meta.data), value = T)
tmp.srt <- st_atr.srt
tmp.srt$SPL_FB <- rowSums(tmp.srt@meta.data[, c("SPL_FB_YapHi", "SPL_FB_YapLo")])
tmp.srt$SPL_IC <- rowSums(tmp.srt@meta.data[, c("SPL_Lymphoid", "SPL_MP_Recruit", 'SPL_MP_Resident')])
p13 <- FeaturePlotST_Dark(tmp.srt,
                          features = grep('^SPL_', colnames(tmp.srt@meta.data), value = T),
                          minvals = rep(0, 10),
                          maxvals = rep(0.5, 10),
                          pt.sizes = pt.sizes*0.3,
                          ncol = 6) &
    theme(aspect.ratio = 0.5)
p13
PlotPDF('1.13.st.decon_feature', 15, 25)
p13
dev.off()


##~~ Plot 14 ST Glycolysis and Yap Scores  ####
st.srt <- FindVariableFeatures(st.srt, assay = 'ST', nfeatures = 10000)
gl <- U(intersect(VariableFeatures(st.srt, assay = 'ST', nfeatures = 10000), c(gl2.mu, gl3.mu, gl4.mu, gl5.mu))) # 138

tmp.srt <- AddModuleScore2(st.srt, features = list(Yap_target$Yap_target_cf_cutrun, gl), assay = 'SCT',
                           names = c('Yap_Activity', 'Glycolysis'), return_z = T)
p14 <- FeaturePlotST_Dark(tmp.srt, features = c('Yap_Activity', 'Glycolysis'), pt.sizes = pt.size.wh,
                          minvals = c(1, 1, 0, 0, 0, 0),
                          maxvals = c(4, 4, 4, 4, 4, 4),
                          ncol = 6, asp = 1)
p14
PlotPDF('1.14.st.yap_glycolysis_score', 28, 8)
p14
dev.off()

WriteCSV(data.frame(Glycolysis_genes_used = paste0('"', gl, '"')), 'PART90.Glycolysis_ST_score_genes')


##~~ Plot 15 ST Glycolysis and Yap Scores Correlation  ####
st.srt <- FindVariableFeatures(st.srt, assay = 'ST', nfeatures = 10000)
gl <- U(intersect(VariableFeatures(st.srt, assay = 'ST', nfeatures = 10000), c(gl2.mu, gl3.mu, gl4.mu, gl5.mu))) # 138

tmp.srt <- AddModuleScore2(st.srt, features = list(Yap_target$Yap_target_cf_cutrun, gl), assay = 'SCT',
                           names = c('Yap_Activity', 'Glycolysis'), return_z = T)

data <- tmp.srt@meta.data[tmp.srt$Condition_group != 'Csf1r_inhibitor' & Cells(tmp.srt) %in% Cells(st_atr.srt), ]
data$Yap_Activity[data$Yap_Activity < 0] <- 0
data$Glycolysis[data$Glycolysis < 0] <- 0
p15 <- ggplot(data) +
    geom_point(aes(x = Yap_Activity, y = Glycolysis, color = Zone), size = 0.5) +
    scale_color_manual(values = mycol_20s) +
    facet_wrap(~Sample) +
    theme_classic() +
    theme(aspect.ratio = 1)
p15
PlotPDF('1.15.scatter.st_yap_glycolysis_score_correlation', 6, 6)
p15
dev.off()


##~~ Plot 16 ST Atrial cell type fraction quantification  ####
gl <- grep('^SPL_', colnames(st_atr.srt@meta.data), value = T)
tmp.srt <- st_atr.srt
tmp.srt$SPL_FB <- rowSums(tmp.srt@meta.data[, c("SPL_FB_YapHi", "SPL_FB_YapLo")])
tmp.srt$SPL_IC <- rowSums(tmp.srt@meta.data[, c("SPL_Lymphoid", "SPL_MP_Recruit", 'SPL_MP_Resident')])
tmp.srt$tmp <- paste(tmp.srt$Zone, tmp.srt$Genotype)
p16 <- BoxPlot2(tmp.srt[, tmp.srt$Sample == c('Control_Veh', 'LatsCKO_Veh') & tmp.srt$Zone %in% c('RA', 'LA')],
                feature = 'SPL_FB', group.by = 'tmp') +
    BoxPlot2(tmp.srt[, tmp.srt$Sample == c('Control_Veh', 'LatsCKO_Veh') & tmp.srt$Zone %in% c('RA', 'LA')],
             feature = 'SPL_IC', group.by = 'tmp') &
    scale_y_continuous(limits = c(0, 0.6))
p16
PlotPDF('1.16.box.st_atrial_composition', 6, 6)
p16
dev.off()


##~~ Plot 17 Box ST Glycolysis and Yap Scores  ####
st.srt <- FindVariableFeatures(st.srt, assay = 'ST', nfeatures = 10000)
gl <- U(intersect(VariableFeatures(st.srt, assay = 'ST', nfeatures = 10000), c(gl2.mu, gl3.mu, gl4.mu, gl5.mu))) # 138

tmp.srt <- AddModuleScore2(st.srt, features = list(Yap_target$Yap_target_cf_cutrun, gl), assay = 'SCT',
                           names = c('Yap_Activity', 'Glycolysis'), return_z = T)
p17.1 <- BoxPlot(tmp.srt[, tmp.srt$Zone %in% c('LA', 'RA' ,'LV', 'RV') & 
                             tmp.srt$Sample %in% c('Control_Veh', 'LatsCKO_Veh')],
                 feature = 'Yap_Activity', group.by = 'Zone', cols = mycol_20s, split.by = 'Sample')

p17.2 <- BoxPlot(tmp.srt[, tmp.srt$Zone %in% c('LA', 'RA' ,'LV', 'RV') &
                             tmp.srt$Sample %in% c('Control_Veh', 'LatsCKO_Veh')],
        feature = 'Glycolysis', group.by = 'Zone', cols = mycol_20s, split.by = 'Sample')
p17 <- p17.1/p17.2
PlotPDF('1.17.box.st_yap_and_glycolysis_score', 6, 6)
p17
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####  Figure 3  ####
##~~ Plot 1 ST OCP Score  ####
tmp.srt <- AddModuleScore2(st.srt, features = list(gl7.mu, gl8.mu), assay = 'SCT',
                           names = c('Chondro_Score', 'Osteo_Score'), return_z = T)
p1 <- FeaturePlotST_Dark(tmp.srt, features = c('Chondro_Score', 'Osteo_Score'),
                          pt.sizes = pt.size,
                          minvals = c(1, 1, 0, 0, 0, 0),
                          maxvals = rep(2.5, 6),
                          ncol = 6, asp = 1)
p1
PlotPDF('2.01.st_feat.st_chondro_osteo_score', 28, 8)
p1
dev.off()


##~~ Plot 2 UMAP CF Cell State  ####
p2 <- DimPlot2(lr_fb.srt, group.by = 'Cell_state', cols = mycol_10, reduction = 'sub_clean_RNA_umap', alpha = 1) |
    DimPlot2(lr_fb.srt, group.by = 'Group2', cols = Color_sample, reduction = 'sub_clean_RNA_umap', alpha = 1)
p2
PlotPDF('2.02.umap.cf_cell_state', 8, 5)
p2
dev.off()


##~~ Plot 3 Feature CF Yap and OCP Score  ####
tmp.srt <- AddModuleScore2(lr_fb.srt, features = list(Yap_target$Yap_target_cf_cutrun, gl9.mu), return_z = T,
                           names = c('Yap_score', 'OCP_score'))
p3 <- FeaturePlot2(tmp.srt, features = c('Yap_score', 'OCP_score'),
                   reduction = 'sub_clean_RNA_umap', min.cutoff = 0, max.cutoff = 2)
p3
PlotPDF('2.03.feat.cf_yap_ocp_score', 8, 5)
p3
dev.off()


##~~ Plot 4 Scatter CF Yap and OCP Score Correlation  ####
tmp.srt <- AddModuleScore2(lr.srt, features = list(Yap_target$Yap_target_cf_cutrun, gl9.mu), return_z = T,
                           names = c('Yap_score', 'OCP_score'))
data <- tmp.srt@meta.data[tmp.srt$Group2 %in% c('Control_RA', 'LatsCKO_RA'), ]
data$Yap_score[data$Yap_score < 0] <- 0
data$OCP_score[data$OCP_score < 0] <- 0
p4 <- ggplot(data) +
    geom_point(aes(x = Yap_score, y = OCP_score, color = Cell_type), size = 0.5) +
    scale_color_manual(values = mycol_10) +
    facet_wrap(~Group2) +
    theme_classic() +
    theme(aspect.ratio = 1)
p4
PlotPDF('2.04.scatter.cf_yap_ocp_score', 8, 4)
p4
dev.off()


##~~ Plot 5 Scatter ST Yap and OCP Score Correlation  ####
tmp.srt <- AddModuleScore2(st.srt, features = list(Yap_target$Yap_target_cf_cutrun, gl9.mu), assay = 'SCT',
                           names = c('Yap_Activity', 'OCP_score'), return_z = T)
data <- tmp.srt@meta.data[tmp.srt$Condition_group != 'Csf1r_inhibitor' & Cells(tmp.srt) %in% Cells(st_atr.srt), ]
data$Yap_Activity[data$Yap_Activity < 0] <- 0
data$OCP_score[data$OCP_score < 0] <- 0
p5 <- ggplot(data) +
    geom_point(aes(x = Yap_Activity, y = OCP_score, color = Zone), size = 0.5) +
    scale_color_manual(values = mycol_20s) +
    facet_wrap(~Sample) +
    theme_classic() +
    theme(aspect.ratio = 1)
p5
PlotPDF('2.05.scatter.st_yap_ocp_score_correlation', 6, 6)
p5
dev.off()


##~~ Plot 6 ST+snRNA Feature Example  ####
p6.1 <- FeaturePlot3(lr_fb.srt, features = c('Gsn', 'Amotl2'),
                     reduction = 'sub_clean_RNA_umap', adjust = 2, pt.size = 2) &
    scale_color_distiller(palette = 'RdYlBu')
p6.2 <- FeaturePlotST_Dark(st.srt, features = c('Gsn', 'Amotl2'),
                           pt.sizes = pt.size,
                           minvals = rep(2, 6),
                           maxvals = rep(4, 6),
                           ncol = 6, asp = 1)
p6 <- wrap_plots(p6.1[[1]], p6.2[[3]], p6.2[[5]],
                 p6.1[[2]], p6.2[[9]], p6.2[[11]],
                 ncol = 3)
p6
PlotPDF('2.06.feat.sn_st_example_gene', 12, 8)
p6
dev.off()


##~~ Plot 7 ST+snRNA Feature OCP Example  ####
p7.1 <- FeaturePlot3(lr_fb.srt, features = c('Sox9', 'Col2a1', 'Acan', 'Runx2'),
                     reduction = 'sub_clean_RNA_umap', adjust = 2, pt.size = 2) &
    scale_color_distiller(palette = 'RdYlBu')
p7.2 <- FeaturePlotST_Dark(st.srt, features = c('Sox9', 'Col2a1', 'Acan', 'Runx2'),
                           pt.sizes = pt.size,
                           minvals = rep(0.5, 6),
                           maxvals = rep(2.5, 6),
                           ncol = 6, asp = 1)
p7 <- wrap_plots(p7.1[[1]], p7.2[[3]], p7.2[[5]],
                 p7.1[[2]], p7.2[[9]], p7.2[[11]],
                 p7.1[[3]], p7.2[[15]], p7.2[[17]],
                 p7.1[[4]], p7.2[[21]], p7.2[[23]],
                 ncol = 3)
p7
PlotPDF('2.07.feat.sn_st_example_ocp_gene', 12, 16)
p7
dev.off()

##~~ Plot 8 Feature OCP global umap  ####
tmp.srt <- AddModuleScore2(lr.srt, features = list(gl7.mu, gl8.mu), return_z = T,
                           names = c('Chondro_score', 'Osteo_score'))
p8 <- FeaturePlot2(tmp.srt, features = c('Chondro_score', 'Osteo_score'),
                   reduction = 'umap', min.cutoff = 1, max.cutoff = 5, pt.size = 0.1)
p8
PlotPDF('2.08.feat.global_chondro_osteo_score', 16, 8)
p8
dev.off()


##~~ Plot 9 Feature Plot CF4 marker expression  ####
p9 <- FeaturePlot2(inh_fb.srt, features = c('Csf1', 'Amotl2', 'Ccn2',  'Vgll3'),
                   reduction = 'sub_clean_RNA_umap', min.cutoff = 0, max.cutoff = 3.5, split.by = 'Group2', order = F)
p9[[16]] <- p9[[16]] + RestoreLegend()
p9
PlotPDF('2.09.feat.cf4_marker', 12, 12)
p9
dev.off()


##~~ Plot 10 Dot Plot CF4 marker expression  ####
p10 <- DotPlot2(inh_fb.srt, features = c('Csf1', 'Amotl2', 'Ccn2',  'Vgll3'),
                group.by = 'Group2', split.by = 'Cell_state', cols = 'RdYlBu') +
        theme_Publication(aspect.ratio = 1)
p10
PlotPDF('2.10.dot.cf4_marker', 6, 6)
p10
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####  Figure 4  ####
##~~ Plot 1 ST Colocalization of Immune cells and YapHi CF  ####
tmp.srt <- st_atr.srt
tmp.srt$Colocal_FBYapHi_ResMP <- GetColocalProb(tmp.srt, meta_features = c('SPL_FB_YapHi', 'SPL_MP_Resident'))
tmp.srt$Csf1_score <- scale(GetAssayData(tmp.srt, assay = 'SCT', layer = 'data')['Csf1', ])
tmp.srt$Csf1r_score <- scale(GetAssayData(tmp.srt, assay = 'SCT', layer = 'data')['Csf1r', ])
tmp.srt$Colocal_Csf1_Csf1r <- GetColocalProb(tmp.srt, meta_features = c('Csf1_score', 'Csf1r_score'))

p1.1 <- FeaturePlotST_Dark(tmp.srt,
                         features = c('SPL_FB_YapHi', 'SPL_MP_Resident',
                                      'Csf1_score', 'Csf1r_score'),
                         minvals = c(0, 0, 0, 0),
                         maxvals = c(0.5, 0.5, 3, 3),
                         pt.sizes = pt.size.atr,
                         ncol = 6) &
    theme(aspect.ratio = 0.5)
tmp.srt$Colocal_FBYapHi_ResMP <- tmp.srt$Colocal_FBYapHi_ResMP + sample(1:10, size = ncol(tmp.srt), replace = T)*1e-10
p1.2 <- FeaturePlotST_Dark(tmp.srt,
                      features = c('Colocal_FBYapHi_ResMP',
                                   'Colocal_Csf1_Csf1r'),
                      minvals = c(0, 0),
                      maxvals = c(3, 3),
                      pt.sizes = pt.size.atr,
                      ncol = 6) &
    scale_colour_viridis_c(option = 'magma', begin = 0.1) &
    theme(aspect.ratio = 0.5)
p1.1
p1.2
PlotPDF('3.01.st_feat.immune_prop_vs_yaphi_cf_prop_colocal_pval', 15, 10)
p1.1
p1.2
dev.off()


##~~ Plot 2 UMAP MP Cell State  ####
p2 <- DimPlot2(lr_mp.srt, group.by = 'Cell_state', cols = mycol_10, reduction = 'sub_clean_RNA_umap', alpha = 1) |
    DimPlot2(lr_mp.srt, group.by = 'Group2', cols = Color_sample, reduction = 'sub_clean_RNA_umap', alpha = 1)
p2
PlotPDF('3.02.umap.immune_cell_state', 12, 5)
p2
dev.off()


##~~ Plot 3 Box Plot MP state composition   ####
mtx <- as.matrix(Table(lr.srt$Cell_state, lr.srt$Group1)+10)
for(i in 1:ncol(mtx)){mtx[, i] <- mtx[, i]/sum(mtx[, i])}
data <- data.frame(LA_1 = mtx[, 3]/mtx[, 1],
                   LA_2 = mtx[, 4]/mtx[, 2],
                   RA_1 = mtx[, 7]/mtx[, 5],
                   RA_2 = mtx[, 8]/mtx[, 6])
data$Cell_state <- rownames(data)
data <- melt(data)
data$LR <- str_split(data$variable, '_', simplify = T)[, 1]
data <- data[data$Cell_state %in% c('MP1_Recruit', 'MP2_Resident', 'MP3_Prol', 'TC'), ]
data$Cell_state <- factor(data$Cell_state, levels = c('MP1_Recruit', 'MP2_Resident', 'MP3_Prol', 'TC'))
p3 <- ggplot(data, aes(x = Cell_state, y = value, group.by = LR)) +
    geom_point() +
    geom_boxplot(outlier.shape = NA, aes(fill = Cell_state)) +
    scale_fill_manual(values = Color_cell_type) +
    #scale_y_continuous(breaks = seq(0, 10, 1))+
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'red3') +
    facet_wrap(~LR, nrow = 1) +
    RotatedAxis() +
    labs(x = '', y = 'Abundance FC vs Controls') +
    theme_classic() +
    theme(aspect.ratio = 2)
p3
PlotPDF('3.03.box.la_ra_immune_cell_subtype_abundance', 5, 5)
p3
dev.off()


##~~ Plot 4 UMAP Inhibitor data global   ####
p4.1 <- DimPlot2(inh.srt, group.by = 'Cell_type', cols = Color_cell_type,
                 pt.size = 0.1, alpha = 1, reduction = 'umap')
tmp.srt <- inh.srt[, unlist(DownsampleByMeta(inh.srt, meta_var = 'Group2', down_to_min_group = T, random = T))]
p4.2 <- DimPlot2(tmp.srt, group.by = 'Group2', cols = Color_sample_inh, pt.size = 0.3, alpha = 1, reduction = 'umap')
p4.1 + p4.2
PlotPDF('3.04.umap.inhibitor_global_cell_type', 15, 7)
p4.1 + p4.2
dev.off()


##~~ Plot 5 Box Inhibitor vs Veh control Cell Type   ####
mtx <- as.matrix(Table(inh.srt$Cell_state, inh.srt$Group1))
for(i in 1:ncol(mtx)){mtx[, i] <- mtx[, i]/sum(mtx[, i])}
data <- data.frame(Cont_1 = mtx[, 3]/mtx[, 1],
                   Cont_2 = mtx[, 5]/mtx[, 1],
                   Lats_1 = mtx[, 4]/mtx[, 2],
                   Lats_2 = mtx[, 6]/mtx[, 2])
data$Cell_type <- rownames(data)
data <- melt(data)
data$GT <- str_split(data$variable, '_', simplify = T)[, 1]
#data <- data[data$Cell_type %in% c('CM', 'EC', 'EpiC', 'FB', 'Immune'), ]
#data$Cell_type <- factor(data$Cell_type, levels = c('CM', 'FB', 'EpiC', 'EC', 'Immune'))
p5 <- ggplot(data, aes(x = Cell_type, y = log2(value), group.by = GT)) +
    #geom_point() +
    geom_boxplot(outlier.shape = NA, aes(fill = Cell_type)) +
    scale_fill_manual(values = Color_cell_state) +
    scale_y_continuous(breaks = seq(-3, 3, 1), limits = c(-3, 1))+
    geom_hline(yintercept = c(-1, 1), linetype = 'dashed', color = 'red3') +
    facet_wrap(~GT, nrow = 1) +
    RotatedAxis() +
    labs(x = '', y = 'Abundance FC vs Vehicle') +
    theme_classic() +
    theme(aspect.ratio = 2)
p5
PlotPDF('3.05.box.cont_lats_cell_type_abundance_inh_vs_veh', 5, 5)
p5
dev.off()

##~~ Plot 6 Box Inhibitor vs Veh control Myeloid Cell State   ####
mtx <- as.matrix(Table(inh.srt$Cell_state, inh.srt$Group1))
for(i in 1:ncol(mtx)){mtx[, i] <- mtx[, i]/sum(mtx[, i])}
data <- data.frame(Cont_1 = mtx[, 3]/mtx[, 1],
                   Cont_2 = mtx[, 5]/mtx[, 1],
                   Lats_1 = mtx[, 4]/mtx[, 2],
                   Lats_2 = mtx[, 6]/mtx[, 2])
data$Cell_type <- rownames(data)
data <- melt(data)
data$GT <- str_split(data$variable, '_', simplify = T)[, 1]
data <- data[data$Cell_type %in% c('MP1_Recruit', 'MP2_Resident', 'MP3_Prol',  'DC', 'BC', 'TC'), ]
data$Cell_type <- factor(data$Cell_type, levels = c('MP1_Recruit', 'MP2_Resident', 'MP3_Prol',  'DC', 'BC', 'TC'))
p6 <- ggplot(data, aes(x = Cell_type, y = log2(value), group.by = GT)) +
    #geom_point() +
    geom_boxplot(outlier.shape = NA, aes(fill = Cell_type)) +
    scale_fill_manual(values = Color_cell_state) +
    #scale_y_continuous(breaks = seq(-3, 3, 1), limits = c(-3, 1))+
    geom_hline(yintercept = c(-1, 1), linetype = 'dashed', color = 'red3') +
    facet_wrap(~GT, nrow = 1) +
    RotatedAxis() +
    labs(x = '', y = 'Abundance FC vs Vehicle') +
    theme_classic() +
    theme(aspect.ratio = 2)
p6
PlotPDF('3.06.box.cont_lats_immune_cell_state_abundance_inh_vs_veh', 5, 5)
p6
dev.off()


##~~ Plot 7 Bar ST Colocalization quantification    ####
tmp.srt <- st_atr.srt
tmp.srt$Colocal_FBYapHi_ResMP <- GetColocalProb(tmp.srt, meta_features = c('SPL_FB_YapHi', 'SPL_MP_Resident'))
tmp.srt$Csf1_score <- scale(GetAssayData(tmp.srt, assay = 'SCT', layer = 'data')['Csf1', ])
tmp.srt$Csf1r_score <- scale(GetAssayData(tmp.srt, assay = 'ST', layer = 'data')['Csf1r', ])
tmp.srt$Colocal_Csf1_Csf1r <- GetColocalProb(tmp.srt, meta_features = c('Csf1_score', 'Csf1r_score'))
tmp.srt$Colocal_FBYapHi_ResMP_bin <- ifelse(tmp.srt$Colocal_FBYapHi_ResMP > -log10(0.05), 'Yes', 'No')
tmp.srt$Colocal_Csf1_Csf1r_bin <- ifelse(tmp.srt$Colocal_Csf1_Csf1r > -log10(0.05), 'Yes', 'No')
tmp.srt <- tmp.srt[, tmp.srt$Zone %in% c('RA', 'LA')]
p7.1 <- CountCellBarPlot(tmp.srt[, tmp.srt$Sample == 'LatsCKO_Veh'],
                         stack.by = 'Colocal_FBYapHi_ResMP_bin', group.by = 'Zone', cols = c('white', 'red3')) +
    labs(title = 'LatsCKO+Veh FBYapHi+ResMP')
p7.2 <- CountCellBarPlot(tmp.srt[, tmp.srt$Sample == 'LatsCKO_Inh'],
                         stack.by = 'Colocal_FBYapHi_ResMP_bin', group.by = 'Zone', cols = c('white', 'red3')) +
    labs(title = 'LatsCKO+Inh FBYapHi+ResMP')
p7.3 <- CountCellBarPlot(tmp.srt[, tmp.srt$Sample == 'LatsCKO_Veh'],
                         stack.by = 'Colocal_Csf1_Csf1r_bin', group.by = 'Zone', cols = c('white', 'red3'))+
    labs(title = 'LatsCKO+Veh Csf1+Csf1r')
p7.4 <- CountCellBarPlot(tmp.srt[, tmp.srt$Sample == 'LatsCKO_Inh'],
                         stack.by = 'Colocal_Csf1_Csf1r_bin', group.by = 'Zone', cols = c('white', 'red3'))+
    labs(title = 'LatsCKO+Inh Csf1+Csf1r')
p7 <- wrap_plots(p7.1, p7.2, p7.3, p7.4, ncol = 2) &
    theme_Publication(aspect.ratio = 1)
p7
PlotPDF('3.07.bar.inh_vs_veh_coloclalization_count', 3, 3)
p7
dev.off()


##~~ Plot 8 Chord Cell Chat LatsCKO FB-Immune   ####
ra_cont_lats.cc_list <- readRDS('analysis/PART26.full_ra_cell_type.cellchat_list.rds')
ra_cont_lats.cont.cc <- ra_cont_lats.cc_list[[1]]
ra_cont_lats.lats.cc <- ra_cont_lats.cc_list[[2]]

tmp.srt <- DropMetaLevels(lr_fb.srt[, lr_fb.srt$Group3 %in% c('Control_RA', 'LatsCKO_RA')])
Idents(tmp.srt) <- 'Group3'
fb_deg <- FindAllMarkers(tmp.srt, only.pos = T)
fb_deg <- fb_deg[fb_deg$p_val_adj < 0.05 & fb_deg$avg_log2FC > .25, ]
lats_lig <- intersect(ra_cont_lats.lats.cc@LR$LRsig$ligand, fb_deg$gene[fb_deg$cluster == 'LatsCKO_RA'])
lats_net <- ra_cont_lats.lats.cc@LR$LRsig[ra_cont_lats.lats.cc@LR$LRsig$ligand %in% lats_lig, ]
View(lats_net[! duplicated(lats_net$pathway_name), ])
poi <- c(
    'TGFb',
    'BMP',
    'GDF',
    'CSF',
    'SEMA3',
    'PERIOSTIN',
    'ANGPTL'
)
PlotPDF('3.08.chord.ra_latscko_fb_to_imm_specific', 6, 6)
netVisual_chord_gene(ra_cont_lats.lats.cc,
                     sources.use = 'FB', targets.use = 'Immune',
                     slot.name = "net",
                     signaling = poi,
                     title.name = 'LatsCKO RA - FB to Immune',
                     scale = F,
                     link.border = T,
                     big.gap = 35, small.gap = 5,
                     lab.cex = 0.5, legend.pos.x = 10, legend.pos.y = 30, annotationTrackHeight = 0.05)

dev.off()


##~~ Plot 10 Track Csf1 activation   ####
tmp.srt <- DropMetaLevels(inh_fb.srt[, ! inh_fb.srt$Group1 %in% c('Control_RA_Inh_2' ,'LatsCKO_RA_Inh_2')])
DefaultAssay(tmp.srt) <- 'ATAC'
p10 <- CoveragePlot(
    object = tmp.srt,
    group.by = 'Group1',
    assay = 'ATAC',
    region = 'Csf1',
    #region = 'chr3-107727909-107773608', ## Csf1
    extend.upstream = 0,
    extend.downstream = 60e3, window = 5e2
)
p10
PlotPDF('3.10.track.csf1', 12, 5)
p10
dev.off()

##~~ Plot 11 Track Igf1r activations   ####
tmp.srt <- DropMetaLevels(inh_fb.srt[, ! inh_fb.srt$Group1 %in% c('Control_RA_Inh_2' ,'LatsCKO_RA_Inh_2')])
DefaultAssay(tmp.srt) <- 'ATAC'
p11 <- CoveragePlot(
    object = tmp.srt,
    group.by = 'Group1',
    assay = 'ATAC',
    #region = 'Igf1r',
    region = 'chr7-67850000-68030000',
    #extend.upstream = 5e4,
    #extend.downstream = 1e4,
    window = 5e2, ymax = 100,
)
p11
PlotPDF('3.11.track.igf1r', 12, 5)
p11
dev.off()


##~~ Plot 12 Box Macrophage glycolysis   ####
tmp.srt <- AddModuleScore2(inh.srt,
                           features = list(gl1.mu,
                                           gl2.mu,
                                           gl3.mu,
                                           gl4.mu,
                                           gl5.mu,
                                           gl6.mu,
                                           U(c(gl2.mu, gl3.mu, gl4.mu)),
                                           gl9.mu
                           ),
                           names = c('Yap_CF_Target',  ## gl1.hs
                                     'Hallmark_Glycolysis',  ## gl2.hs
                                     'KEGG_Glycolysis_Glucogenesis',  ## gl3.hs
                                     'Reactome_Glycolysis',  ## gl4.hs
                                     'GO_Glycolysis',  ## gl5.hs
                                     'GO_Glycolysis_Up_in_LatsCKO',  ## gl6.hs
                                     'Union_Glycolysis',  ## U(c(gl2.hs, gl3.hs, gl4.hs, gl5.hs)
                                     'OCP_Marker'  ## gl9.hs
                           ),
                           return_z = T)
mp_inh.srt <- tmp.srt[, tmp.srt$Cell_state %in% c('MP1_Recruit', 'MP2_Resident', 'MP3_Prol')]
a <- split(mp_inh.srt$Reactome_Glycolysis, mp_inh.srt$Group2)
pval_1v2 <- t.test(a$Control_RA_Veh, a$LatsCKO_RA_Veh)$p.value # 0.000000000004521874
pval_2v4 <- t.test(a$LatsCKO_RA_Veh, a$LatsCKO_RA_Inh)$p.value # 0.7456699
pval_3v4 <- t.test(a$Control_RA_Inh, a$LatsCKO_RA_Inh)$p.value # 0.0000000000000000000134866

p12 <- BoxPlot(mp_inh.srt, feature =  'Reactome_Glycolysis', group.by = 'Group2', cols = Color_sample_inh,
               split.by = 'Cell_state') +
    theme_Publication(aspect.ratio = 1) +
        RotatedAxis()
p12
PlotPDF('3.12.box.mp_glycolysis', 9, 3)
p12
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####  Figure 5  ####
##~~ Plot 1 Box snRNA LatsCKO veh vs inh FB, Yap OCP and glycolysis score  ####
x <- intersect(gl9.mu, VariableFeatures(inh.srt))
tmp.srt <- AddModuleScore2(inh.srt,
                           features = list(gl1.mu,
                                           gl2.mu,
                                           gl3.mu,
                                           gl4.mu,
                                           gl5.mu,
                                           gl6.mu,
                                           U(c(gl2.mu, gl3.mu, gl4.mu)),
                                           x
                           ),
                           names = c('Yap_CF_Target',  ## gl1.hs
                                     'Hallmark_Glycolysis',  ## gl2.hs
                                     'KEGG_Glycolysis_Glucogenesis',  ## gl3.hs
                                     'Reactome_Glycolysis',  ## gl4.hs
                                     'GO_Glycolysis',  ## gl5.hs
                                     'GO_Glycolysis_Up_in_LatsCKO',  ## gl6.hs
                                     'Union_Glycolysis',  ## U(c(gl2.hs, gl3.hs, gl4.hs, gl5.hs)
                                     'OCP_Marker'  ## gl9.hs
                           ),
                           return_z = T)
tmp2.srt <- DropMetaLevels(tmp.srt[, tmp.srt$Cell_state %in% c('FB4_YapHi')])

a <- split(tmp2.srt$Yap_CF_Target, tmp2.srt$Group2)
a_pval_2v4 <- t.test(a$LatsCKO_RA_Veh, a$LatsCKO_RA_Inh)$p.value # 0.0004295836
b <- split(tmp2.srt$OCP_Marker, tmp2.srt$Group2)
b_pval_2v4 <- t.test(b$LatsCKO_RA_Veh, b$LatsCKO_RA_Inh)$p.value # 0.02192241
c <- split(tmp2.srt$Reactome_Glycolysis, tmp2.srt$Group2)
c_pval_2v4 <- t.test(c$LatsCKO_RA_Veh, c$LatsCKO_RA_Inh)$p.value # 0.01088195

p1.1 <- BoxPlot(tmp2.srt, feature =  'Yap_CF_Target', group.by = 'Group2', cols = Color_sample_inh[c(2, 4)])
p1.2 <- BoxPlot(tmp2.srt, feature =  'OCP_Marker', group.by = 'Group2', cols = Color_sample_inh[c(2, 4)])
p1.3 <- BoxPlot(tmp2.srt, feature =  'Reactome_Glycolysis', group.by = 'Group2', cols = Color_sample_inh[c(2, 4)])
p1 <- wrap_plots(p1.1, p1.2, p1.3, ncol = 1) &
    theme_Publication(aspect.ratio = 3)
p1
PlotPDF('4.01.box.fb4_yap_ocp_glycolysis', 3, 9)
p1
dev.off()

##~~ Plot 2 Box snRNA LatsCKO veh vs inh FB, Yap OCP and glycolysis score in steady state CF  ####
x <- intersect(gl9.mu, VariableFeatures(inh.srt))
tmp.srt <- AddModuleScore2(inh.srt,
                           features = list(gl1.mu,
                                           gl2.mu,
                                           gl3.mu,
                                           gl4.mu,
                                           gl5.mu,
                                           gl6.mu,
                                           U(c(gl2.mu, gl3.mu, gl4.mu)),
                                           x
                           ),
                           names = c('Yap_CF_Target',  ## gl1.hs
                                     'Hallmark_Glycolysis',  ## gl2.hs
                                     'KEGG_Glycolysis_Glucogenesis',  ## gl3.hs
                                     'Reactome_Glycolysis',  ## gl4.hs
                                     'GO_Glycolysis',  ## gl5.hs
                                     'GO_Glycolysis_Up_in_LatsCKO',  ## gl6.hs
                                     'Union_Glycolysis',  ## U(c(gl2.hs, gl3.hs, gl4.hs, gl5.hs)
                                     'OCP_Marker'  ## gl9.hs
                           ),
                           return_z = T)
tmp2.srt <- DropMetaLevels(tmp.srt[, ! tmp.srt$Cell_state %in% c('FB4_YapHi')])

a <- split(tmp2.srt$Yap_CF_Target, tmp2.srt$Group2)
a_pval_1v2 <- t.test(a$Control_RA_Veh, a$LatsCKO_RA_Veh)$p.value # 0.6287861
a_pval_2v4 <- t.test(a$LatsCKO_RA_Veh, a$LatsCKO_RA_Inh)$p.value # 0.1802767
a_pval_3v4 <- t.test(a$Control_RA_Inh, a$LatsCKO_RA_Inh)$p.value # 0.00001855734
b <- split(tmp2.srt$OCP_Marker, tmp2.srt$Group2)
b_pval_1v2 <- t.test(b$Control_RA_Veh, b$LatsCKO_RA_Veh)$p.value # 0.0001590439
b_pval_2v4 <- t.test(b$LatsCKO_RA_Veh, b$LatsCKO_RA_Inh)$p.value # 0.0000004195964
b_pval_3v4 <- t.test(b$Control_RA_Inh, b$LatsCKO_RA_Inh)$p.value # 0.000004166929
c <- split(tmp2.srt$Reactome_Glycolysis, tmp2.srt$Group2)
c_pval_1v2 <- t.test(c$Control_RA_Veh, c$LatsCKO_RA_Veh)$p.value # 0.000000000000003189291
c_pval_2v4 <- t.test(c$LatsCKO_RA_Veh, c$LatsCKO_RA_Inh)$p.value # 0.00000000000000004757691
c_pval_3v4 <- t.test(c$Control_RA_Inh, c$LatsCKO_RA_Inh)$p.value # 0.00000000000000000000000000000000000000000000236012

p2.1 <- BoxPlot(tmp2.srt, feature =  'Yap_CF_Target', group.by = 'Group2', cols = Color_sample_inh)
p2.2 <- BoxPlot(tmp2.srt, feature =  'OCP_Marker', group.by = 'Group2', cols = Color_sample_inh)
p2.3 <- BoxPlot(tmp2.srt, feature =  'Reactome_Glycolysis', group.by = 'Group2', cols = Color_sample_inh)
p2 <- wrap_plots(p2.1, p2.2, p2.3, ncol = 1) &
    theme_Publication(aspect.ratio = 3)
p2
PlotPDF('4.02.box.fb123_yap_ocp_glycolysis', 3, 12)
p2
dev.off()


tmp3.srt <- DropMetaLevels(tmp.srt[, ! tmp.srt$Cell_state %in% c('FB4_YapHi')])


##~~ Plot 3 ST LatsCKO veh vs inh  Sox9 and Col2a1  ####
p3 <- FeaturePlotST_Dark(st.srt, features = c('Sox9', 'Col2a1'), pt.sizes = pt.size.wh,
                         minvals = c(1, 1),
                         maxvals = c(3, 3),
                         ncol = 6, asp = 1) &
    NoLegend()
p3
PlotPDF('4.03.st.sox9_col2a1', 28, 8)
p3
dev.off()


##~~ Plot 4 Box snRNA LatsCKO veh vs inh FB4 Sox9 and Col2a1  ####
p4.1 <- BoxPlot(st.srt[, st.srt$Zone == 'RA' & st.srt$Sample %in% c('LatsCKO_Veh', 'LatsCKO_Inh')],
               feature = c('Sox9'),
               group.by = 'Sample', cols = Color_sample_inh[c(2, 4)])
p4.2 <- BoxPlot(st.srt[, st.srt$Zone == 'RA' & st.srt$Sample %in% c('LatsCKO_Veh', 'LatsCKO_Inh')],
                feature = c('Col2a1'),
                group.by = 'Sample', cols = Color_sample_inh[c(2, 4)])
p4 <- wrap_plots(p4.1, p4.2, ncol = 2) &
    theme_Publication(aspect.ratio = 2) &
        RotatedAxis()
p4
PlotPDF('4.04.box.st_ra_sox9_col2a1', 4, 4)
p4
dev.off()


##~~ Plot 4 Violin snRNA LatsCKO veh vs inh FB4 Sox9 and Col2a1  ####
p4 <- VlnPlot2(st.srt[, st.srt$Zone == 'RA' & st.srt$Sample %in% c('LatsCKO_Veh', 'LatsCKO_Inh')],
               features = c('Sox9', 'Col2a1'),
               group.by = 'Sample', pt.size = -1, assay = 'SCT', adjust = 1.2, cols = Color_sample_inh[c(2, 4)]) &
    theme_Publication(aspect.ratio = 2)
p4
PlotPDF('4.04.vln.st_ra_sox9_col2a1', 4, 4)
p4
dev.off()


##~~ Plot 5 ST Control vs LatsCKO Vegll3, Col11a1, Col12a1  ####
p5 <- FeaturePlotST_Dark(st.srt, features = c('Vgll3', 'Col11a1', 'Col12a1'), pt.sizes = pt.size.wh,
                         minvals = rep(0, 3),
                         maxvals = rep(3, 3),
                         ncol = 6, asp = 1) &
    NoLegend()
p5
PlotPDF('4.05.st.vgll3_col11a1_col12a1', 28, 12)
p5
dev.off()


##~~ Plot 6 Feature snRNA Control vs LatsCKO Vegll3, Col11a1, Col12a1  ####
p6 <- FeaturePlot3(lr_fb.srt, features = c('Vgll3', 'Col11a1', 'Col12a1'), reduction = 'sub_clean_RNA_umap') &
    scale_color_distiller(palette = 'RdYlBu') &
    theme_Publication(aspect.ratio = 1)
p6
PlotPDF('4.06.feat_density.vgll3_col11a1_col12a1', 4, 4)
p6
dev.off()

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----





####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Not Used  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~  Plot 1  ST Feature - ST Glycolysis Score  ####
tmp.srt <- AddModuleScore2(st.srt, features = list(U(c(gl2.mu, gl3.mu, gl4.mu))),
                           names = 'Glycolysis', return_z = T)

data <- tmp.srt@meta.data
data$X <- tmp.srt@reductions$Spatial@cell.embeddings[,1]
data$Y <- tmp.srt@reductions$Spatial@cell.embeddings[,2]
data$Glycolysis[data$Glycolysis > 3] <- 3
data$Glycolysis[data$Glycolysis < -1] <- -1

p1 <- ggplot(data) +
    geom_point(aes(x = X, y = Y, color = Glycolysis), size = 0.7) +
    scale_color_distiller(limits = c(-1, 3), palette = 'RdYlBu') +
    facet_wrap(~Sample, scales = 'free') +
    labs(color = 'Glycolysis Genes\nExp. Z-score') +
    theme_classic() +
    NoAxes() +
    theme(aspect.ratio = 1,  strip.background = element_blank())
p1
PlotPDF('4.01.st_feat.glycolysis', 7, 3)
p1
dev.off()


##~~  Plot 2  ST Feature - ST YAP Score  ####
tmp.srt <- AddModuleScore2(st.srt, features = list(gl1.mu), names = 'YAP_activity', return_z = T)

data <- tmp.srt@meta.data
data$X <- tmp.srt@reductions$Spatial@cell.embeddings[,1]
data$Y <- tmp.srt@reductions$Spatial@cell.embeddings[,2]
data$YAP_activity[data$YAP_activity > 3] <- 3
data$YAP_activity[data$YAP_activity < -1] <- -1

p2 <- ggplot(data) +
    geom_point(aes(x = X, y = Y, color = YAP_activity), size = 0.7) +
    scale_color_distiller(limits = c(-1, 3), palette = 'RdYlBu') +
    facet_wrap(~Sample, scales = 'free') +
    labs(color = 'Yap Targets\nExp. Z-score') +
    theme_classic() +
    NoAxes() +
    theme(aspect.ratio = 1,  strip.background = element_blank())
p2
PlotPDF('4.02.st_feat.yap_activity', 7, 3)
p2
dev.off()

##~~  Plot 3  Scatter - ST Glycolysis + Yap Correlation  ####
tmp.srt <- AddModuleScore2(st.srt, features = list(gl1.mu, U(c(gl2.mu, gl3.mu, gl4.mu))),
                           names = c('YAP_activity', 'Glycolysis'), return_z = T)
tmp.srt$Coloc_Glyco_YapAct <- GetColocalProb(tmp.srt, meta_features = c('Glycolysis', 'YAP_activity'))

data <- tmp.srt@meta.data
data$X <- tmp.srt@reductions$Spatial@cell.embeddings[,1]
data$Y <- tmp.srt@reductions$Spatial@cell.embeddings[,2]
data$Coloc_Glyco_YapAct[data$Coloc_Glyco_YapAct > 4] <- 4

p3 <- ggplot(data) +
    geom_point(aes(x = X, y = Y, color = Coloc_Glyco_YapAct), size = 0.7) +
    scale_color_distiller(limits = c(0, 4), palette = 'RdYlBu') +
    facet_wrap(~Sample, scales = 'free') +
    labs(color = '-Log10 p', caption = 'Yap Activity and Glycolysis Colocalization P Value') +
    theme_classic() +
    NoAxes() +
    theme(aspect.ratio = 1,  strip.background = element_blank())
p3
PlotPDF('4.03.st_feat.colocal_yap_activity_glycolysis', 7, 3)
p3
dev.off()

tmp.srt <- AddModuleScore2(full.srt, features = list(gl1.mu, U(c(gl2.mu, gl3.mu, gl4.mu))), assay = 'RNA',
                           names = c('YAP_activity', 'Glycolysis'), return_z = T)
val <- split(tmp.srt[, tmp.srt$Cell_type == 'FB']$Glycolysis, tmp.srt[, tmp.srt$Cell_type == 'FB']$Group1)
t.test(val[[1]], val[[3]])$p.value

BoxPlot(tmp.srt[, tmp.srt$Cell_type == 'FB'], feature = 'YAP_activity', group.by = 'Group1') +
    BoxPlot(tmp.srt[, tmp.srt$Cell_type == 'FB'], feature = 'Glycolysis', group.by = 'Group1')


ggplot(tmp.srt[, tmp.srt$Cell_type == 'FB' & tmp.srt$Group1 == 'LatsCKO + Vehicle']@meta.data) +
    geom_point(aes(x = YAP_activity, y = Glycolysis, color = Cell_state)) +
    geom_smooth(aes(x = YAP_activity, y = Glycolysis), formula = y ~ x, level = 0.99) +
    theme_classic() +
    theme(aspect.ratio = 1)


tmp.srt <- AddModuleScore2(sn.srt, features = list(gl1.mu, U(c(gl2.mu, gl3.mu, gl4.mu))), assay = 'RNA',
                           names = c('YAP_activity', 'Glycolysis'), return_z = T)
ggplot(tmp.srt[, tmp.srt$Cell_type == 'CF' & tmp.srt$Sample == 'Mut RA2']@meta.data) +
    geom_point(aes(x = YAP_activity, y = Glycolysis, color = clustering_level3)) +
    geom_smooth(aes(x = YAP_activity, y = Glycolysis), formula = y ~ x, level = 0.99) +
    theme_classic() +
    theme(aspect.ratio = 1)



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Test  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inh_fb.srt <- AddModuleScore2(inh_fb.srt, features = list(Yap_target$Yap_target_cf_cutrun), names = 'Yap_score',
                              return_z = T)


tmp.srt <- inh_fb.srt[, unlist(DownsampleByMeta(inh_fb.srt, meta_var = 'Group2', down_to_min_group = T, random = T))]
data <- data.frame(tmp.srt@meta.data)
data$X <- tmp.srt@reductions$sub_clean_RNA_umap@cell.embeddings[, 1]
data$Y <- tmp.srt@reductions$sub_clean_RNA_umap@cell.embeddings[, 2]

ggplot(data) +
    geom_point(aes(x = X, y = Y, color = Yap_score)) +
    scale_color_distiller(palette = 'RdYlBu') +
    facet_wrap(~Group2) +
    theme_classic() +
    theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
VlnPlot2(tmp.srt, feature = 'Yap_score', group.by = 'Group2')
BoxPlot(tmp.srt[, tmp.srt$Cell_state == 'FB4_YapHi'], feature = 'Yap_score', group.by = 'Group2')

FeaturePlot2(inh_fb.srt, features = 'Yap_score', split.by = 'Group2', reduction = 'sub_clean_RNA_umap')


