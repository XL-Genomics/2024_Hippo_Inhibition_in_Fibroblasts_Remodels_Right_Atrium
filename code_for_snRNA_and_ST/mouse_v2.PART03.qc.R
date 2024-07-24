####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Part <- 'PART03_QC'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2019_atrial_latscko_crtsai/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scATACseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('local', Ver, Part, 'mouse', Project_dir = '2019_atrial_latscko_crtsai', Data_drive = 'bree')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load sample information  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample_meta.df <- read.csv(paste0(Docu_dir, 'mouse_sample_meta.csv'))
studies <- U(sample_meta.df$study)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
srt <- readRDS('integrated/PART02.merged.cbn_atac.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----



####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  QC RNA (4 Multi + 10 snRNA Samples)  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n_sample <- 14
vars <- c('nCount_CBN', 'nFeature_CBN', 'pct_mito_CBN')

cutoffs <- data.frame(name = 1:n_sample)

p <- VlnPlot(srt, features = vars, ncol = 3, log = F, pt.size = 0, cols = mycol_20) + NoLegend()
p
PlotPDF('01.1.vln.meta_pre_filter', 16, 4)
p
dev.off()

nCount_CBN_cutoffs <- rep(NA, 2*n_sample)
nFeature_CBN_cutoffs <- nCount_CBN_cutoffs
pct_mito_CBN_cutoffs <- nCount_CBN_cutoffs[1:n_sample]
max_mito <- 5
min_nCount_CBN <- 500
min_nFeature_CBN <- 200
for(i in 1:n_sample){
        meta <- srt@meta.data[srt$Group1 == levels(srt$Group1)[i], ]
        nCount_CBN_cutoffs[c(i, i+n_sample)] <- GetOutlier(meta$nCount_CBN, iqr_multiplier = 3)
        nFeature_CBN_cutoffs[c(i, i+n_sample)] <- GetOutlier(meta$nFeature_CBN, iqr_multiplier = 3)
        pct_mito_CBN_cutoffs[i] <- GetOutlier(meta$pct_mito_CBN, iqr_multiplier = 3)[2]
        if(nCount_CBN_cutoffs[i] < min_nCount_CBN){nCount_CBN_cutoffs[i] <- min_nCount_CBN}
        if(nFeature_CBN_cutoffs[i] < min_nFeature_CBN){nFeature_CBN_cutoffs[i] <- min_nFeature_CBN}
        if(pct_mito_CBN_cutoffs[i] > max_mito){pct_mito_CBN_cutoffs[i] <- max_mito}
}

p <- list()
for(i in 1:n_sample){
        meta <- srt@meta.data[srt$Group1 == levels(srt$Group1)[i], ]
        p[[i]] <- ggplot(meta, aes(x = Group1, y = nCount_CBN)) +
                geom_violin(fill = mycol_30[3]) +
                geom_hline(yintercept = c(nCount_CBN_cutoffs[i],
                                          nCount_CBN_cutoffs[i+n_sample]),
                           color = mycol_14[1], size = 1) +
                labs(title = levels(meta$Group1)[i],
                     subtitle = paste0(round(nCount_CBN_cutoffs[i], 0) ,
                                       ' < RNA Count < ',
                                       round(nCount_CBN_cutoffs[i+n_sample], 0))) +
                scale_y_continuous(limits = c(0, 4e4))
}
p2 <- wrap_plots(p, ncol = 4) &
        NoLegend() &
        theme_classic() &
        theme(aspect.ratio = 2, axis.title = element_blank())
p2
PlotPDF('01.2.vln_merged_cutoff_nCount_CBN', 12, 18)
p2
dev.off()


p <- list()
for(i in 1:n_sample){
        meta <- srt@meta.data[srt$Group1 == levels(srt$Group1)[i], ]
        p[[i]] <- ggplot(meta, aes(x = Group1, y = nFeature_CBN)) +
                geom_violin(fill = mycol_30[3]) +
                geom_hline(yintercept = c(nFeature_CBN_cutoffs[i],
                                          nFeature_CBN_cutoffs[i+n_sample]),
                           color = mycol_14[1], size = 1) +
                labs(title = levels(meta$Group1)[i],
                     subtitle = paste0(round(nFeature_CBN_cutoffs[i], 0) ,
                                       ' < Gene Count < ',
                                       round(nFeature_CBN_cutoffs[i+n_sample], 0))) +
                scale_y_continuous(limits = c(0, 12e3))
}
p2 <- wrap_plots(p, ncol = 4) &
        NoLegend() &
        theme_classic() &
        theme(aspect.ratio = 2, axis.title = element_blank())
p2
PlotPDF('01.3.vln_merged_cutoff_nFeature_CBN', 12, 18)
p2
dev.off()


p <- list()
for(i in 1:n_sample){
        meta <- srt@meta.data[srt$Group1 == levels(srt$Group1)[i], ]
        p[[i]] <- ggplot(meta, aes(x = Group1, y = pct_mito_CBN)) +
                geom_violin(fill = mycol_30[3]) +
                geom_hline(yintercept = pct_mito_CBN_cutoffs[i], color = mycol_14[1], size = 1) +
                labs(title = levels(meta$Group1)[i],
                     subtitle = paste0('Mito% < ', round(pct_mito_CBN_cutoffs[i], 2))) +
                scale_y_continuous(limits = c(0, 5))
}
p2 <- wrap_plots(p, ncol = 4) &
        NoLegend() &
        theme_classic() &
        theme(aspect.ratio = 2, axis.title = element_blank())
p2
PlotPDF('01.4.vln_merged_cutoff_mito_CBN', 12, 18)
p2
dev.off()


## Save cutoff values
cutoffs$nCount_max <- nCount_CBN_cutoffs[(1:n_sample)+n_sample]
cutoffs$nCount_min <- nCount_CBN_cutoffs[(1:n_sample)]
cutoffs$nFeature_max <- nFeature_CBN_cutoffs[(1:n_sample)+n_sample]
cutoffs$nFeature_min <- nFeature_CBN_cutoffs[(1:n_sample)]
cutoffs$mito_max <- pct_mito_CBN_cutoffs


srt$LowQ_CBN <- F
for(i in 1:n_sample){
        srt$LowQ_CBN[srt$Group1 == levels(srt$Group1)[i] &
                             (srt$pct_mito_CBN>pct_mito_CBN_cutoffs[i] |
                                      srt$nCount_CBN<nCount_CBN_cutoffs[i] |
                                      srt$nCount_CBN>nCount_CBN_cutoffs[i+n_sample] |
                                      srt$nFeature_CBN<nFeature_CBN_cutoffs[i] |
                                      srt$nFeature_CBN>nFeature_CBN_cutoffs[i+n_sample])] <- T
}
toss_by_rna <- Cells(srt)[srt$LowQ_CBN]
L(toss_by_rna) ## 18369
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Generate QC ATAC Metrics (4 Multiome samples)  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
n_sample <- 4

DefaultAssay(srt) <- 'ATAC'
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(srt) <- annotations


atac.srt <- DropMetaLevels(srt[, srt$Group1 %in% c('Control_RA_Veh',
                                                   'LatsCKO_RA_Veh',
                                                   'Control_RA_Inh_1',
                                                   'LatsCKO_RA_Inh_1')])
levels(atac.srt$Group1)

# compute nucleosome signal score per cell
atac.srt <- NucleosomeSignal(atac.srt, assay = 'ATAC')
VlnPlot2(atac.srt, features = 'nucleosome_signal', group.by = 'Group1', y.max = 1.5)

# compute TSS enrichment score per cell
# srt <- TSSEnrichment(srt, fast = F, assay = 'ATAC')
atac.srt$TSS.enrichment <- 0 ## Not run due to Seurat v5 incompatibility

total_frag_list <- list()
samples <- levels(atac.srt$Group1)
# add fraction of reads in peaks
path <- c('/Volumes/shire/data/scmulti/2023_Csf1r_CTsai/matrix/2023_csf1r_CTsai_flox_V/outs/atac_fragments.tsv.gz',
          '/Volumes/shire/data/scmulti/2023_Csf1r_CTsai//matrix/2023_csf1r_CTsai_CKO_V/outs/atac_fragments.tsv.gz',
          '/Volumes/shire/data/scmulti/2023_Csf1r_CTsai/matrix/2023_csf1r_CTsai_flox_G/outs/atac_fragments.tsv.gz',
          '/Volumes/shire/data/scmulti/2023_Csf1r_CTsai/matrix/2023_csf1r_CTsai_CKO_G/outs/atac_fragments.tsv.gz'
) ## !! Note the order of the fragment file must match the order of Group1 levels !!
for(i in 1:n_sample) {
        total_frag_list[[i]] <- CountFragments(path[i], cells = atac.srt$Cell_barcode[atac.srt$Group1 == samples[i]])
        rownames(total_frag_list[[i]]) <- colnames(atac.srt)[atac.srt$Group1 == samples[i]]
}
x <- total_frag_list[[1]]
for(i in 2:n_sample){x <- rbind(x, total_frag_list[[i]])}
O(x$CB, as.vector(atac.srt$Cell_barcode))
x <- x[Cells(atac.srt), ]
atac.srt$Fragments <- x$frequency_count
VlnPlot2(atac.srt, features = 'Fragments', group.by = 'Group1', log = T)

atac.srt <- FRiP(atac.srt, assay = 'ATAC', total.fragments = 'Fragments')
VlnPlot2(atac.srt, features = 'FRiP', group.by = 'Group1')

atac.srt$blacklist_fraction <- FractionCountsInRegion(atac.srt, assay = 'ATAC', regions = blacklist_mm10)
VlnPlot2(atac.srt, features = 'blacklist_fraction', group.by = 'Group1')

## Add meta data back to full seurat
srt$Nucleosome_signal_ATAC <- NA
srt$Nucleosome_pct_ATAC <- NA
srt$TSS_enrich_ATAC <- NA
srt$Fragments_ATAC <- NA
srt$Reads_in_peak_pct_ATAC <- NA
srt$Blacklist_pct_ATAC <- NA

srt@meta.data[colnames(atac.srt),  c("Nucleosome_signal_ATAC", 
                                     "Nucleosome_pct_ATAC",
                                     "TSS_enrich_ATAC",
                                     "Fragments_ATAC", 
                                     "Reads_in_peak_pct_ATAC",
                                     "Blacklist_pct_ATAC")] <-
        atac.srt@meta.data[, c("nucleosome_signal", 
                               "nucleosome_percentile",
                               "TSS.enrichment",
                               "Fragments", 
                               "FRiP",
                               "blacklist_fraction")]
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Set QC Threshold for ATAC  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
atac.srt <- DropMetaLevels(srt[, srt$Group1 %in% c('Control_RA_Veh',
                                                   'LatsCKO_RA_Veh',
                                                   'Control_RA_Inh_1',
                                                   'LatsCKO_RA_Inh_1')])

vars <- c('Nucleosome_signal_ATAC',
          'TSS_enrich_ATAC', 
          'Blacklist_pct_ATAC', 
          'Reads_in_peak_pct_ATAC', 
          'Fragments_ATAC')
p <- VlnPlot2(atac.srt, features = vars, ncol = 2, log = F, pt.size = 0, cols = mycol_10) + NoLegend()
p
PlotPDF('02.1.vln.meta_pre_filter', 8, 12)
print(p)
dev.off()

FRiP_cutoffs <- rep(0, n_sample*2)
fragments_cutoffs <- rep(0, n_sample*2)
max_fragments <- 100e3
min_fragments <- 100
min_FRiP <- 0.15
max_Blacklist_pct <- 0.1
max_Nucleosome_signal <- 2
min_TSS_enrich <- -1 ## not used 

for(i in 1:n_sample){
        meta <- atac.srt@meta.data[atac.srt$Group1 == levels(atac.srt$Group1)[i], ]
        FRiP_cutoffs[c(i, i+n_sample)] <- GetOutlier(meta$Reads_in_peak_pct_ATAC, iqr_multiplier = 3)
        if(FRiP_cutoffs[i] < min_FRiP){FRiP_cutoffs[i] <- min_FRiP}
        fragments_cutoffs[c(i, i+n_sample)] <- GetOutlier(meta$Fragments_ATAC, iqr_multiplier = 3)
        if(fragments_cutoffs[i] < min_fragments){fragments_cutoffs[i] <- min_fragments}
        if(fragments_cutoffs[i+n_sample] > max_fragments){fragments_cutoffs[i+n_sample] <- max_fragments}
}

p <- list()
for(i in 1:n_sample){
        meta <- atac.srt@meta.data[atac.srt$Group1 == levels(atac.srt$Group1)[i], ]
        p[[i]] <- ggplot(meta, aes(x = Group1, y = Reads_in_peak_pct_ATAC)) +
                geom_violin(fill = mycol_30[3]) +
                geom_hline(yintercept = FRiP_cutoffs[i], color = mycol_14[1], size = 1) +
                labs(title = levels(meta$Group1)[i],
                     subtitle = paste0(round(FRiP_cutoffs[i], 2) , ' < % Reads in Peak')) +
                scale_y_continuous(limits = c(0, 1))
}
p2 <- wrap_plots(p, ncol = 4) &
        NoLegend() &
        theme_classic() &
        theme(aspect.ratio = 2, axis.title = element_blank())
PlotPDF('02.2.vln_merged_cutoff_FRiP_ATAC', 12, 8)
p2
dev.off()

p <- list()
for(i in 1:n_sample){
        meta <- atac.srt@meta.data[atac.srt$Group1 == levels(atac.srt$Group1)[i], ]
        p[[i]] <- ggplot(meta, aes(x = Group1, y = log1p(Fragments_ATAC))) +
                geom_violin(fill = mycol_30[3]) +
                geom_hline(yintercept = c(log1p(fragments_cutoffs[i]),
                                          log1p(fragments_cutoffs[i+n_sample])),
                           color = mycol_14[1], size = 1) +
                labs(title = levels(meta$Group1)[i],
                     subtitle = paste0(round(fragments_cutoffs[i], 0) ,
                                       ' < Fragments < ',
                                       round(fragments_cutoffs[i+n_sample], 0))) +
                scale_y_continuous(limits = c(1, 15))
}
p2 <- wrap_plots(p, ncol = 4) &
        NoLegend() &
        theme_classic() &
        theme(aspect.ratio = 2, axis.title = element_blank())
PlotPDF('02.3.vln_merged_cutoff_Fragments_ATAC', 12, 8)
p2
dev.off()

p <- list()
for(i in 1:n_sample){
        meta <- atac.srt@meta.data[atac.srt$Group1 == levels(atac.srt$Group1)[i], ]
        p[[i]] <- ggplot(meta, aes(x = Group1, y = Nucleosome_signal_ATAC)) +
                geom_violin(fill = mycol_30[3]) +
                geom_hline(yintercept = max_Nucleosome_signal, color = mycol_14[1], size = 1) +
                labs(title = levels(meta$Group1)[i],
                     subtitle = 'Nucleosome signal < 2') +
                scale_y_continuous(limits = c(0, 2.5))
}
p2 <- wrap_plots(p, ncol = 4) &
        NoLegend() &
        theme_classic() &
        theme(aspect.ratio = 2, axis.title = element_blank())
PlotPDF('02.4.vln_merged_cutoff_Nucleosome_signal_ATAC', 12, 8)
p2
dev.off()

# p <- list()
# for(i in 1:n_sample){
#         meta <- atac.srt@meta.data[atac.srt$Group1 == levels(atac.srt$Group1)[i], ]
#         p[[i]] <- ggplot(meta, aes(x = Group1, y = TSS_enrich_ATAC)) +
#                 geom_violin(fill = mycol_30[3]) +
#                 geom_hline(yintercept = min_TSS_enrich, color = mycol_14[1], size = 1) +
#                 labs(title = levels(meta$Group1)[i],
#                      subtitle = 'TSS enrichment > 2') +
#                 scale_y_continuous(limits = c(0, 30))
# }
# p2 <- wrap_plots(p, ncol = 4) &
#         NoLegend() &
#         theme_classic() &
#         theme(aspect.ratio = 2, axis.title = element_blank())
# PlotPDF('02.5.vln_merged_cutoff_TSS_enrich_ATAC', 12, 8)
# p2
# dev.off()

p <- list()
for(i in 1:n_sample){
        meta <- atac.srt@meta.data[atac.srt$Group1 == levels(atac.srt)[i], ]
        p[[i]] <- ggplot(meta, aes(x = Group1, y = Blacklist_pct_ATAC)) +
                geom_violin(fill = mycol_30[3]) +
                geom_hline(yintercept = max_Blacklist_pct, color = mycol_14[1], size = 1) +
                labs(title = levels(meta$Group1)[i],
                     subtitle = 'Blacklist % < 5%') +
                scale_y_continuous(limits = c(0, 0.15))
}
p2 <- wrap_plots(p, ncol = 4) &
        NoLegend() &
        theme_classic() &
        theme(aspect.ratio = 2, axis.title = element_blank())
PlotPDF('02.6.vln_merged_cutoff_Blacklist_pct_ATAC', 12, 8)
p2
dev.off()

atac.srt$LowQ_ATAC <- F
for(i in 1:n_sample){
        atac.srt$LowQ_ATAC[atac.srt$Group1 == levels(atac.srt$Group1)[i] &
                                   (atac.srt$Nucleosome_signal_ATAC > max_Nucleosome_signal |
                                            atac.srt$TSS_enrich_ATAC < min_TSS_enrich |
                                            atac.srt$Blacklist_pct_ATAC > max_Blacklist_pct |
                                            atac.srt$Fragments_ATAC < fragments_cutoffs[i] |
                                            atac.srt$Fragments_ATAC > fragments_cutoffs[i+n_sample] |
                                            atac.srt$Reads_in_peak_pct_ATAC < FRiP_cutoffs[i])] <- T
}
srt$LowQ_ATAC <- F
srt$LowQ_ATAC[colnames(atac.srt)] <- atac.srt$LowQ_ATAC

## Save cutoff values
cutoffs[, c('nucleosome_signal_max', 
            'tss_enrich_min', 
            'blacklist_pct_max', 
            'fragments_min', 
            'fragments_max', 
            'reads_in_peak_pct_min')] <- NA
cutoffs$nucleosome_signal_max[1:n_sample] <- max_Nucleosome_signal
cutoffs$tss_enrich_min[1:n_sample] <- min_TSS_enrich
cutoffs$blacklist_pct_max[1:n_sample] <- max_Blacklist_pct
cutoffs$fragments_min[1:n_sample] <- fragments_cutoffs[(1:n_sample)]
cutoffs$fragments_max[1:n_sample] <- fragments_cutoffs[(1:n_sample)+n_sample]
cutoffs$reads_in_peak_pct_min[1:n_sample] <- FRiP_cutoffs[(1:n_sample)]
cutoffs$name <- NULL
cutoffs$Sample_id <- sample_meta.df$Sample_id

toss_by_atac <- colnames(srt)[srt$LowQ_ATAC] ## 1879
L(toss_by_atac)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(srt@meta.data, 'integrated/PART03.merged.pre_filter.srt_meta.rds')
saveRDS(cutoffs, 'analysis/PART03.cell_filtered.low_quality_cutoffs.df.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Remove cells  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
O(toss_by_atac, toss_by_rna) ## 1879   18369    458
L(union(toss_by_rna, toss_by_atac)) ## 19790
L(union(toss_by_rna, toss_by_atac))/ncol(srt) ## 17.2%

p <- DimPlot2(srt, group.by = c('LowQ_CBN', 'LowQ_ATAC'), cols = c('grey80', 'red'))
p
PlotPDF('03.umap.pc_umap_low_qual_cells', 12, 6)
p
dev.off()

clean.srt <- srt[, ! colnames(srt) %in% union(toss_by_rna, toss_by_atac)]
DefaultAssay(clean.srt) <- 'CBN'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(clean.srt, 'integrated/PART03.merged.flt.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
