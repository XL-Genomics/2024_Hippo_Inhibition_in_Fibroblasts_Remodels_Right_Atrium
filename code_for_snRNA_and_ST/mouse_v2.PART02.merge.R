####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Part <- 'PART02_Merge'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2019_atrial_latscko_crtsai/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scATACseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('local', Ver, Part, 'mouse', Project_dir = '2019_atrial_latscko_crtsai', Data_drive = 'bree')

plan("multisession", workers = 4)
options(future.globals.maxSize = 100*1024^3) # for 100 Gb RAM
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Sample Information  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sample_meta.df <- read.csv(paste0(Docu_dir, 'mouse_sample_meta.csv'))
studies <- U(sample_meta.df$Study)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####  ~~~~~~~~~~~~~~~~~~~~~~~  ####
####  Process snMultiome data  ####
####  ~~~~~~~~~~~~~~~~~~~~~~~  ####


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####   Load Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
raw1.srt <- readRDS('individual/1.2023_Csf1r_CTsai.raw.srt.rds')
cbn1.srt <- readRDS('individual/1.2023_Csf1r_CTsai.cbn.srt.rds')

srt_raw_list <- list(raw1.srt[, raw1.srt$library == '2023_csf1r_CTsai_flox_V'],
                     raw1.srt[, raw1.srt$library == '2023_csf1r_CTsai_flox_G'],
                     raw1.srt[, raw1.srt$library == '2023_csf1r_CTsai_CKO_V'],
                     raw1.srt[, raw1.srt$library == '2023_csf1r_CTsai_CKO_G'])
srt_cbn_list <- list(cbn1.srt[, cbn1.srt$library == '2023_csf1r_CTsai_flox_V'],
                     cbn1.srt[, cbn1.srt$library == '2023_csf1r_CTsai_flox_G'],
                     cbn1.srt[, cbn1.srt$library == '2023_csf1r_CTsai_CKO_V'],
                     cbn1.srt[, cbn1.srt$library == '2023_csf1r_CTsai_CKO_G'])

n_sample <- 4
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Reduce ATAC Peaks  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Create a unified set of peaks to quantify in each dataset
combined.peaks <- Signac::reduce(x = c(srt_raw_list[[1]]@assays$ATAC@ranges,
                                       srt_raw_list[[2]]@assays$ATAC@ranges,
                                       srt_raw_list[[3]]@assays$ATAC@ranges,
                                       srt_raw_list[[4]]@assays$ATAC@ranges))

## Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths < 10e3 & peakwidths > 20]
combined.peaks ## 109741

## Create Fragment objects
path <- c(
        '/Volumes/shire/data/scmulti/2023_Csf1r_CTsai//matrix/2023_csf1r_CTsai_flox_V/outs/atac_fragments.tsv.gz',
        '/Volumes/shire/data/scmulti/2023_Csf1r_CTsai/matrix/2023_csf1r_CTsai_flox_G/outs/atac_fragments.tsv.gz',
        '/Volumes/shire/data/scmulti/2023_Csf1r_CTsai//matrix/2023_csf1r_CTsai_CKO_V/outs/atac_fragments.tsv.gz',
        '/Volumes/shire/data/scmulti/2023_Csf1r_CTsai/matrix/2023_csf1r_CTsai_CKO_G/outs/atac_fragments.tsv.gz'
)

frag_list <- list()
for(i in 1:n_sample) {
        frag_list[[i]] <- CreateFragmentObject(path = path[i],
                                               cells = as.vector(srt_raw_list[[i]]$orig.name))
}

## Quantify peaks in each dataset
counts_list <- list()
for(i in 1:n_sample) {
        counts_list[[i]] <- FeatureMatrix(fragments = frag_list[[i]],
                                          features = combined.peaks,
                                          cells = as.vector(srt_raw_list[[i]]$orig.name))
}

## Create new ATAC object
new_atac_list <- list()
for(i in 1:n_sample) {
        atac_assay <- CreateChromatinAssay(counts_list[[i]], fragments = frag_list[[i]], genome = 'mm10')
        meta <- srt_raw_list[[i]]@meta.data
        rownames(meta) <- meta$orig.name
        meta <- meta[colnames(atac_assay), ]
        new_atac_list[[i]] <- CreateSeuratObject(atac_assay, assay = 'ATAC', meta.data = meta)
        new_atac_list[[i]] <- RenameCells(new_atac_list[[i]],
                                          new.names = paste0(new_atac_list[[i]]$study,
                                                             '_',
                                                             new_atac_list[[i]]$library,
                                                             ':',
                                                             new_atac_list[[i]]$orig.name),
                                          for.merge = F)
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Merge ATAC Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
merge_atac.srt <- merge(x = new_atac_list[[1]], y = new_atac_list[2:n_sample])
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
Annotation(merge_atac.srt) <- annotations
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Merge Cell Ranger GEX Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
merge_raw_gex.srt <- merge(x = srt_raw_list[[1]], y = srt_raw_list[2:n_sample])
merge_raw_gex.srt <- RenameAssays(merge_raw_gex.srt, 'RNA' = 'RAW')
DefaultAssay(merge_raw_gex.srt) <- 'RAW'
merge_raw_gex.srt <- DietSeurat(merge_raw_gex.srt, assays = 'RAW')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Combine CellRanger GEX and ATAC Assays  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
identical(as.vector(Cells(merge_raw_gex.srt)), as.vector(Cells(merge_atac.srt))) ## TRUE

srt <- merge_atac.srt
srt[['RAW']] <- merge_raw_gex.srt@assays$RAW
DefaultAssay(srt) <- 'RAW'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Merge CellBender GEX Data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
merge_cbn_gex.srt <- merge(x = srt_cbn_list[[1]], y = srt_cbn_list[2:n_sample])
merge_cbn_gex.srt <- RenameAssays(merge_cbn_gex.srt, 'RNA' = 'CBN')
merge_cbn_gex.srt <- JoinLayers(merge_cbn_gex.srt)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Filter for Consensus Cells of CellRanger and CellBender  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
O(Cells(merge_cbn_gex.srt), Cells(srt)) ## 36037  32975  32852
consensus <- intersect(Cells(merge_cbn_gex.srt), Cells(srt))

srt <- srt[, consensus]
merge_cbn_gex.srt <- merge_cbn_gex.srt[, consensus]

cbn_count <- GetAssayData(merge_cbn_gex.srt, assay = 'CBN', layer = 'counts')
cbn_count <- cbn_count[, Cells(srt)]
identical(Cells(srt, assay = 'RAW'), colnames(cbn_count)) ## TRUE

srt[['CBN']] <- CreateAssayObject(counts = cbn_count, min.cells = 0, min.features = 0)
DefaultAssay(srt) <- 'CBN'
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Add Meta Data from Sample Sheet  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meta <- srt@meta.data
meta <- inner_join(x = meta, y = sample_meta.df, join_by(library == Name_on_disk))
identical(Cells(srt), paste0(meta$study, '_', meta$library, ':', meta$orig.name)) ## TRUE
rownames(meta) <- Cells(srt)
meta$orig.ident <- NULL
colnames(meta)[colnames(meta) == 'study'] <- 'Study'
colnames(meta)[colnames(meta) == 'library'] <- 'Library'
colnames(meta)[colnames(meta) == 'orig.name'] <- 'Orig.name'
srt@meta.data <- meta

Idents(srt) <- srt$Group1

multi.srt <- srt
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####  ~~~~~~~~~~~~~~~~~~~~~~~  ####
####  Process snRNA data       ####
####  ~~~~~~~~~~~~~~~~~~~~~~~  ####


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Merge raw   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rds_name <- paste0(1:3, '.', studies, '.raw.srt.rds')
data.list <- list()
ncells <- 0
for(i in 2:L(rds_name)) {
        data.list[[i]] <- readRDS(paste0('individual/', rds_name[i]))
        message(data.list[[i]]$study[1], ' ', ncol(data.list[[i]]))
        ncells <- ncells + ncol(data.list[[i]])
}
merged.raw.srt <- merge(data.list[[2]], data.list[3:L(data.list)])
c(ncells, ncol(merged.raw.srt)) ## 82300 82300
rm(data.list)
gc()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Merge cbn datasets   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rds_name <- paste0(1:3, '.', studies, '.cbn.srt.rds')
data.list <- list()
ncells <- 0
for(i in 2:L(rds_name)) {
        data.list[[i]] <- readRDS(paste0('individual/', rds_name[i]))
        message(data.list[[i]]$study[1], ' ', ncol(data.list[[i]]))
        ncells <- ncells + ncol(data.list[[i]])
}
merged.cbn.srt <- merge(data.list[[2]], data.list[3:L(data.list)])
c(ncells, ncol(merged.cbn.srt)) ## 85423 85423
rm(data.list)
gc()

## Intersect datasets
O(Cells(merged.raw.srt), Cells(merged.cbn.srt)) ## common 81955
consensus_cell <- intersect(Cells(merged.raw.srt), Cells(merged.cbn.srt))
merged.cbn.srt <- merged.cbn.srt[, consensus_cell]
consensus_gene <- intersect(rownames(merged.raw.srt), rownames(merged.cbn.srt))
merged.cbn.srt <- DietSeurat(merged.cbn.srt, features = consensus_gene)
merged.cbn.srt <- JoinLayers(merged.cbn.srt)

merged.cbn.srt <- RenameAssays(merged.cbn.srt, assay.name = 'RNA', new.assay.name = 'CBN')
merged.cbn.srt <- NormalizeData(merged.cbn.srt)

rm(merged.raw.srt)
gc()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Add Meta Data from Sample Sheet  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meta <- merged.cbn.srt@meta.data
meta <- inner_join(x = meta, y = sample_meta.df, join_by(library == Name_on_disk))
identical(Cells(merged.cbn.srt), paste0(meta$study, '_', meta$library, ':', meta$orig.name)) ## TRUE
rownames(meta) <- Cells(merged.cbn.srt)
meta$orig.ident <- NULL
colnames(meta)[colnames(meta) == 'study'] <- 'Study'
colnames(meta)[colnames(meta) == 'library'] <- 'Library'
colnames(meta)[colnames(meta) == 'orig.name'] <- 'Orig.name'
merged.cbn.srt@meta.data <- meta

sn.srt <- merged.cbn.srt
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----





####  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  ####
####  Merge snMultiome and snRNA data       ####
####  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  ####

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
srt <- merge(sn.srt, multi.srt)
srt <- JoinLayers(srt)
srt <- DietSeurat(srt, layers = c('counts', 'data'), assays = c('CBN', 'ATAC'))


srt$Group1 <- revalue(srt$Sample_id, 
                      c('P01_S001' = 'Control_RA_Veh',
                        'P01_S002' = 'Control_RA_Inh_1',
                        'P01_S003' = 'LatsCKO_RA_Veh',
                        'P01_S004' = 'LatsCKO_RA_Inh_1',
                        'P02_S001' = 'Control_RA_1',
                        'P02_S002' = 'Control_RA_2',
                        'P02_S003' = 'LatsCKO_RA_1',
                        'P02_S004' = 'LatsCKO_RA_2',
                        'P02_S005' = 'Control_LA_1',
                        'P02_S006' = 'Control_LA_2',
                        'P02_S007' = 'LatsCKO_LA_1',
                        'P02_S008' = 'LatsCKO_LA_2',
                        'P03_S001' = 'Control_RA_Inh_2',
                        'P03_S002' = 'LatsCKO_RA_Inh_2'
                      ))
srt$Group1 <- factor(srt$Group1, levels = c(
        'Control_LA_1',
        'Control_LA_2',
        'LatsCKO_LA_1',
        'LatsCKO_LA_2',
        'Control_RA_1',
        'Control_RA_2',
        'LatsCKO_RA_1',
        'LatsCKO_RA_2',
        'Control_RA_Veh',
        'LatsCKO_RA_Veh',
        'Control_RA_Inh_1',
        'LatsCKO_RA_Inh_1',
        'Control_RA_Inh_2',
        'LatsCKO_RA_Inh_2'
))

srt$Group2 <- revalue(srt$Group1, 
                      c('Control_LA_1' = 'Control_LA',
                        'Control_LA_2' = 'Control_LA',
                        'LatsCKO_LA_1' = 'LatsCKO_LA',
                        'LatsCKO_LA_2' = 'LatsCKO_LA',
                        'Control_RA_1' = 'Control_RA',
                        'Control_RA_2' = 'Control_RA',
                        'LatsCKO_RA_1' = 'LatsCKO_RA',
                        'LatsCKO_RA_2' = 'LatsCKO_RA',
                        'Control_RA_Veh' = 'Control_RA_Veh',
                        'LatsCKO_RA_Veh' = 'LatsCKO_RA_Veh',
                        'Control_RA_Inh_1' = 'Control_RA_Inh',
                        'LatsCKO_RA_Inh_1' = 'LatsCKO_RA_Inh',
                        'Control_RA_Inh_2' = 'Control_RA_Inh',
                        'LatsCKO_RA_Inh_2' = 'LatsCKO_RA_Inh'))
srt$Group2 <- factor(srt$Group2, levels = c('Control_LA', 'LatsCKO_LA',
                                            'Control_RA', 'LatsCKO_RA',
                                            'Control_RA_Veh', 'LatsCKO_RA_Veh',
                                            'Control_RA_Inh', 'LatsCKO_RA_Inh'))

srt$Group3 <- factor(revalue(srt$Group2, c('Control_RA_Veh' = 'Control_RA', 'LatsCKO_RA_Veh' = 'LatsCKO_RA')), 
                     levels = c('Control_LA', 'LatsCKO_LA', 
                                'Control_RA', 'LatsCKO_RA',
                                'Control_RA_Inh', 'LatsCKO_RA_Inh'))

Table(srt$Group1, srt$Group2)
Table(srt$Group1, srt$Group3)


####  Check sex labeling  ####
tmp.srt <- AddModuleScore2(srt, features = list(chrX_genes$chrX_only_genes, chrY_genes$chrY_only_genes),
                           return_z = T, names = c('ChrX', 'ChrY'))
p1 <- BoxPlot(tmp.srt, feature = 'ChrX', group.by = 'Group1', cols = mycol_20) +
        scale_y_continuous(limits = c(-3, 3))
p2 <- BoxPlot(tmp.srt, feature = 'ChrY', group.by = 'Group1', cols = mycol_20) +
        scale_y_continuous(limits = c(-3, 3))
p3 <- ggplot(srt@meta.data) + geom_point(aes(x = Group1, y = Sex, color = Group1), size = 5) +
        scale_color_manual(values = mycol_20) +
        theme_classic() +
        RotatedAxis()
p <- wrap_plots(p1, p2, p3, ncol = 1)
PlotPDF('1.1.box.sex_genes_per_sample', 6, 12)
p
dev.off()

####  Recalculate cell cycle scores  ####
srt <- AddModuleScore2(srt, features = list(s.genes, g2m.genes), names = c('S.Score', 'G2M.Score'), return_z = T)

####  Recalculate pct mito  ####
srt <- PercentageFeatureSet(srt, pattern = '^mt-', col.name = 'pct_mito_CBN', assay = 'CBN')

####  Remove duplicated  ####
srt$Study.1 <- NULL
srt$nCount_RNA <- NULL
srt$nFeature_RNA <- NULL
colnames(srt@meta.data)[colnames(srt@meta.data) == 'Orig.name'] <- 'Cell_barcode' 

Idents(srt) <- 'Group1'

####  Check genotype and condition labeling  ####
p <- DotPlot2(srt, features = c('eGFP', 'creERT')) |
        ggplot(srt@meta.data) + geom_point(aes(y = Group1, x = Genotype_l, color = Group1), size = 5) +
        scale_color_manual(values = mycol_20) +
        scale_y_discrete(limits = rev) +
        theme_classic() +
        RotatedAxis()
p
PlotPDF('1.2.dot.genotype_condition_per_sample', 12, 8)
p
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Get PC UMAP Embedding   ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
srt <- ProcessSrt_std(srt, assay = 'CBN', do.umap = T, npcs = 30)
p1 <- DimPlot2(srt, cols = mycol_20)
p1

dsp.srt <- srt[, unlist(DownsampleByMeta(srt, meta_var = 'Group1', down_to_min_group = T, random = T))]
p2 <- DimPlot2(dsp.srt, group.by = 'Group1', cols = mycol_20)
p3 <- DimPlot2(dsp.srt, group.by = 'Group2', cols = mycol_10)
p4 <- DimPlot2(dsp.srt, group.by = 'Group3', cols = mycol_10)

PlotPDF('2.1.umap.pc_umap_full_data', 15, 15)
p1 + p2 + p3 + p4
dev.off()

srt <- DietSeurat(srt, assays = c('CBN', 'ATAC'), layers = c('counts', 'data'), dimreducs = c('pca', 'umap'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(srt, 'integrated/PART02.merged.cbn_atac.srt.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
