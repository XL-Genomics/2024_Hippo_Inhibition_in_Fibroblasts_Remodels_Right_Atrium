####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Part <- 'PART04_Doublet'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2019_atrial_latscko_crtsai/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scATACseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('local', Ver, Part, 'mouse', Project_dir = '2019_atrial_latscko_crtsai', Data_drive = 'bree')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library('reticulate')
use_condaenv(condaenv = "scrublet", conda = "/Users/Felix/Conda/anaconda3/bin/python")
scr <- import('scrublet')
plt <- import("matplotlib.pyplot")
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Seurat  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
merged.dlt.srt <- readRDS('integrated/PART03.merged.flt.srt.rds')
merged.dlt.srt <- DietSeurat(merged.dlt.srt, assays = 'CBN', dimreducs = names(merged.dlt.srt@reductions))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###  Identify doublets for each dataset  (linear processing) ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
doublet_rate <- 0.15 ## Assuming 15% doublet formation rate
all_Doublet_SC <- c()
all_Scrublet <- c()

all_samples <- levels(merged.dlt.srt$Group1)
for(i in 1:L(all_samples)){
        gc()
        message(paste0('Processing ', all_samples[i], ' ...'))
        tmp.srt <- merged.dlt.srt[, merged.dlt.srt$Group1 == all_samples[i]]
        results <- GetDoublet(srt_obj = tmp.srt, doublet_rate = doublet_rate, dimN.var.toal = 0.85)
        all_Doublet_SC <- c(all_Doublet_SC, results[[1]])
        ## plot umap
        PlotPDF(paste0('1.', str_pad(i, pad = 0, width = 2), '.', all_samples[i], '.doublets_found'), 10, 5)
        print(results[[2]])
        dev.off()
        all_Scrublet <- c(all_Scrublet, results[[3]])
        ## plot scrublet histogram
        PlotPDF(paste0('1.', str_pad(i, pad = 0, width = 2), '.', all_samples[i], '.scrublet_hist'), 8, 4)
        print(plt$show(results[[4]]$plot_histogram()[[1]]))
        dev.off()
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        saveRDS(all_Doublet_SC, 'analysis/PART04.cell_filtered.scrublet_doublets.rds')
        saveRDS(all_Scrublet,   'analysis/PART04.cell_filtered.scrublet_score.rds')
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
}

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Evaluate doublets  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## save doublets to the main Seurat
merged.dlt.srt$Doublet_SC <- F
merged.dlt.srt$Doublet_SC[all_Doublet_SC] <- T
merged.dlt.srt$Doublet_SC_score <- NA
merged.dlt.srt$Doublet_SC_score[names(all_Scrublet)] <- all_Scrublet


p <- DimPlot2(merged.dlt.srt, cols = 'grey75', reduction = 'umap', raster = F, pt.size = 0.001,
             cells.highlight = all_Doublet_SC, cols.highlight = 'red', sizes.highlight = 0.001) +
        labs(title = paste0('Total cells: ', ncol(merged.dlt.srt), '  Scrublet Doublets found: ', L(all_Doublet_SC))) +
        NoLegend() +
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
p
PlotPDF('2.01.merged_filtered.scrublet_doublets', 10, 10)
p
dev.off()

p <- FeaturePlot2(merged.dlt.srt, reduction = 'umap', raster = F, features = 'Doublet_SC_score') +
        NoLegend() +
        theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())
p
PlotPDF('2.02.merged_filtered.scrublet_score', 10, 10)
p
dev.off()

p <- VlnPlot2(merged.dlt.srt, features = 'Doublet_SC_score', group.by = 'Group1', cols = mycol_20) +
        NoLegend()
p
PlotPDF('2.03.merged_filtered.scrublet_score_per_sample', 5, 5)
p
dev.off()

p1 <- CountCellBarPlot(merged.dlt.srt, group.by = 'Group1', stack.by = 'Doublet_SC',
                       cols = mycol_10, percentage = F)
p2 <- CountCellBarPlot(merged.dlt.srt, group.by = 'Group1', stack.by = 'Doublet_SC',
                       cols = mycol_10, percentage = T)
wrap_plots(p1, p2, ncol = 1)
PlotPDF('2.04.merged_filtered.pct_sc_dlt_study', 10, 15, onefile = T)
wrap_plots(p1, p2, ncol = 1)
dev.off()
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(merged.dlt.srt@meta.data, 'integrated/PART04.merged.dlt.srt_meta.rds') ## Create new
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
