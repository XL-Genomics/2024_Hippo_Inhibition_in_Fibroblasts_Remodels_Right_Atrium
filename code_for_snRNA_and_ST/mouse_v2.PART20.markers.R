####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Part <- 'PART20_Markers'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2019_atrial_latscko_crtsai/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scATACseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('local', Ver, Part, 'mouse', Project_dir = '2019_atrial_latscko_crtsai', Data_drive = 'bree')

plan("multisession", workers = 24)
options(future.globals.maxSize = 100*1024^3) # for 100 Gb RAM
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load data  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
srt <- readRDS('integrated/PART19.consolidated.srt.rds')
srt <- srt[, srt$Non_ambiguous]
srt <- DropMetaLevels(srt)

srt@misc$marker <- list()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Cell type markers  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compute Cell_type RNA markers
Idents(srt) <- 'Cell_type'
Cell_type_marker <- FindAllMarkers(srt,
                                   assay = 'RNA',
                                   test.use = 'wilcox',
                                   only.pos = T,
                                   logfc.threshold = 1,
                                   random.seed = 505)
Cell_type_marker <- Cell_type_marker[Cell_type_marker$p_val_adj < 0.01, ]
Table(Cell_type_marker$cluster)

## Compute Cell_type ATAC markers
Cell_type_peaks <- FindAllMarkers(srt,
                                  assay = 'ATAC',
                                  only.pos = T,
                                  logfc.threshold = 1,
                                  test.use = 'LR',
                                  return.thresh = 0.01,
                                  latent.vars = 'nCount_ATAC')
Cell_type_peaks <- Cell_type_peaks[Cell_type_peaks$p_val_adj < 0.01, ]
Table(Cell_type_peaks$cluster)

srt@misc$marker$Cell_type_rna_marker <- Cell_type_marker
srt@misc$marker$Cell_type_atac_marker <- Cell_type_peaks
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(srt@misc$marker, 'analysis/PART20.markers.srt_misc.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Cell state markers  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Compute Cell state RNA markers
Idents(srt) <- 'Cell_state'

mk_list <- list()
ct <- c('CM', 'FB', 'EC', 'Immune', 'EpiC', 'Mural')
for(i in 1:L(ct)) {
    message(ct[i])
    sub.srt <- srt[, srt$Cell_type == ct[i]]
    mk_list[[ct[i]]] <- FindAllMarkers(sub.srt,
                                       assay = 'RNA',
                                       test.use = 'wilcox',
                                       only.pos = T,
                                       random.seed = 505)
    mk_list[[ct[i]]] <- mk_list[[ct[i]]][mk_list[[ct[i]]]$p_val_adj < 0.01, ]
}

## Compute Cell state ATAC markers
peak_list <- list()
ct <- c('CM', 'FB', 'EC', 'Immune', 'EpiC', 'Mural')
for(i in 1:L(ct)) {
    message(ct[i])
    sub.srt <- srt[, srt$Cell_type == ct[i]]
    peak_list[[ct[i]]] <- FindAllMarkers(sub.srt,
                                         assay = 'ATAC',
                                         only.pos = T,
                                         test.use = 'LR',
                                         random.seed = 505,
                                         latent.vars = 'nCount_ATAC')
    peak_list[[ct[i]]] <- peak_list[[ct[i]]][peak_list[[ct[i]]]$p_val_adj < 0.01, ]
}

srt@misc$marker$Cell_state_rna_marker <- mk_list
srt@misc$marker$Cell_state_atac_marker <- peak_list
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Link DAP with genes  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
DefaultAssay(srt) <- 'ATAC'

## For Cell_type
x <- Cell_type_peaks$gene
y <- ClosestFeature(srt, regions = x)
Cell_type_peaks <- cbind(Cell_type_peaks, y)

## For Cell_state
ct <- c('CM', 'FB', 'EC', 'Immune', 'Mural')
for(i in 1:L(ct)) {
    x <- peak_list[[ct[i]]]$gene
    y <- ClosestFeature(srt, regions = x)
    peak_list[[ct[i]]] <- cbind(peak_list[[ct[i]]], y)
}

srt@misc$marker$Cell_type_atac_marker <- Cell_type_peaks
srt@misc$marker$Cell_state_atac_marker <- peak_list
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
saveRDS(srt@misc$marker, 'analysis/PART20.markers.srt_misc.rds') ## Update
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
