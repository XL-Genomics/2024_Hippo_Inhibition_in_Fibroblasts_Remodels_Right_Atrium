####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Part <- 'PART26_CellChat'

Code_dir <- paste0('~/Documents/Bioinformatics/project/2019_atrial_latscko_crtsai/code/mouse_v', Ver, '/')

source(Sys.readlink(paste0(Code_dir, 'src/bioinformatics.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/scATACseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/stRNAseq.R')))
source(Sys.readlink(paste0(Code_dir, 'src/geneset_mouse.R')))
source(paste0(Code_dir, 'mouse_v', Ver, '.helper_functions.R'))

InitiateProject('local', Ver, Part, 'mouse', Project_dir = '2019_atrial_latscko_crtsai', Data_drive = 'bree')

suppressMessages(library('CellChat'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Seurat  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
full.srt <- readRDS('integrated/PART19.consolidated.srt.rds')
full.srt <- DropMetaLevels(full.srt[, full.srt$Non_ambiguous])
full.srt <- DietSeurat(full.srt, assays = 'RNA', dimreducs = names(full.srt@reductions))

##~~  Clean RNA data of LA+RA Cont+LatsCKO ####
lr.srt <- DropMetaLevels(full.srt[, full.srt$Study == '2023_Latscko_CTsai'])
lr_cont_la.srt <- DropMetaLevels(lr.srt[, lr.srt$Group2 == 'Control_LA'])
lr_cont_ra.srt <- DropMetaLevels(lr.srt[, lr.srt$Group2 == 'Control_RA'])
lr_lats_la.srt <- DropMetaLevels(lr.srt[, lr.srt$Group2 == 'LatsCKO_LA'])
lr_lats_ra.srt <- DropMetaLevels(lr.srt[, lr.srt$Group2 == 'LatsCKO_RA'])

##~~  Clean RNA data of RA Csf1r Inhibitor Cont+LatsCKO ####
inh.srt <- DropMetaLevels(full.srt[, full.srt$Study == '2023_Csf1r_CTsai'])
inh_cont_veh.srt <- DropMetaLevels(inh.srt[, inh.srt$Group2 == 'Control_RA_Veh'])
inh_lats_veh.srt <- DropMetaLevels(inh.srt[, inh.srt$Group2 == 'LatsCKO_RA_Veh'])
inh_cont_inh.srt <- DropMetaLevels(inh.srt[, inh.srt$Group2 == 'Control_RA_Inh'])
inh_lats_inh.srt <- DropMetaLevels(inh.srt[, inh.srt$Group2 == 'LatsCKO_RA_Inh'])
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Prepare CellChat database  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~---
## Adding FGF13 interactions per NicheNet Database
DoCellChat <- function(srt, group.by, assay = 'RNA', secreted_only = F, trim = 0.1){
    CellChatDB <- CellChatDB.mouse
    if(secreted_only){CellChatDB <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")}
    ####  Build CellChat Object  ####
    cch <- createCellChat(object = srt, group.by = group.by, assay = assay)
    cch <- addMeta(cch, meta = srt@meta.data)
    groupSize <- as.numeric(table(cch@idents))
    cch@DB <- CellChatDB
    ####  Preprocessing the expression data for cell-cell communication analysis  ####
    cch <- subsetData(cch) # subset the expression data of signaling genes for saving computation cost
    cch <- identifyOverExpressedGenes(cch)
    cch <- identifyOverExpressedInteractions(cch)
    cch <- projectData(cch, PPI.human)
    ####  Inference of cell-cell communication network  ####
    cch <- computeCommunProb(cch, type = "truncatedMean", trim = trim)
    # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
    cch <- filterCommunication(cch)
    # CellChat computes the communication probability on signaling pathway level
    cch <- computeCommunProbPathway(cch)
    cch <- aggregateNet(cch)
    return(cch)
}
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  CellChat workflow for Control vs Lats in LA and RA  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----
## 1. RA Control vs LatsCKO
ra_ctrl_cell_type.cch <- DoCellChat(lr_cont_ra.srt, group.by = 'Cell_type', secreted_only = T, trim = 0.1)
ra_lats_cell_type.cch <- DoCellChat(lr_lats_ra.srt, group.by = 'Cell_type', secreted_only = T, trim = 0.1)
full_ra_cell_type.list <- list(ra_ctrl_cell_type.cch, ra_lats_cell_type.cch)
full_ra_cell_type.cch <- mergeCellChat(list(Control = ra_ctrl_cell_type.cch,
                                            LatsCKO = ra_lats_cell_type.cch),
                                       add.names = c('Control', 'LatsCKO'))

ra_ctrl_cell_state.cch <- DoCellChat(lr_cont_ra.srt, group.by = 'Cell_state', secreted_only = T, trim = 0.1)
ra_lats_cell_state.cch <- DoCellChat(lr_lats_ra.srt, group.by = 'Cell_state', secreted_only = T, trim = 0.1)
full_ra_cell_state.list <- list(ra_ctrl_cell_state.cch, ra_lats_cell_state.cch)
full_ra_cell_state.cch <- mergeCellChat(list(Control = ra_ctrl_cell_state.cch,
                                             LatsCKO = ra_lats_cell_state.cch),
                                        add.names = c('Control', 'LatsCKO'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(full_ra_cell_type.cch, 'analysis/PART26.full_ra_cell_type.cellchat.rds')
saveRDS(full_ra_cell_type.list, 'analysis/PART26.full_ra_cell_type.cellchat_list.rds')
saveRDS(full_ra_cell_state.cch, 'analysis/PART26.full_ra_cell_state.cellchat.rds')
saveRDS(full_ra_cell_state.list, 'analysis/PART26.full_ra_cell_state.cellchat_list.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## 2. Control LA vs RA
ra_ctrl_cell_type.cch <- DoCellChat(lr_cont_ra.srt, group.by = 'Cell_type', secreted_only = T, trim = 0.1)
la_ctrl_cell_type.cch <- DoCellChat(lr_cont_la.srt, group.by = 'Cell_type', secreted_only = T, trim = 0.1)
full_cont_cell_type.list <- list(ra_ctrl_cell_type.cch, la_ctrl_cell_type.cch)
full_cont_cell_type.cch <- mergeCellChat(list(LA = la_ctrl_cell_type.cch,
                                              RA = ra_ctrl_cell_type.cch),
                                         add.names = c('LA', 'RA'))

ra_ctrl_cell_state.cch <- DoCellChat(lr_cont_ra.srt, group.by = 'Cell_state', secreted_only = T, trim = 0.1)
la_ctrl_cell_state.cch <- DoCellChat(lr_cont_la.srt, group.by = 'Cell_state', secreted_only = T, trim = 0.1)
full_cont_cell_state.list <- list(ra_ctrl_cell_state.cch, la_ctrl_cell_state.cch)
full_cont_cell_state.cch <- mergeCellChat(list(LA = la_ctrl_cell_state.cch,
                                               RA = ra_ctrl_cell_state.cch),
                                          add.names = c('LA', 'RA'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(full_cont_cell_type.cch, 'analysis/PART26.full_cont_cell_type.cellchat.rds')
saveRDS(full_cont_cell_type.list, 'analysis/PART26.full_cont_cell_type.cellchat_list.rds')
saveRDS(full_cont_cell_state.cch, 'analysis/PART26.full_cont_cell_state.cellchat.rds')
saveRDS(full_cont_cell_state.list, 'analysis/PART26.full_cont_cell_state.cellchat_list.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## 3. RA WT Veh Control vs Csf1r Inhibitor

ctrl_veh_cell_type.cch <- DoCellChat(inh_cont_veh.srt, group.by = 'Cell_type', secreted_only = T, trim = 0.1)
ctrl_inh_cell_type.cch <- DoCellChat(inh_cont_inh.srt, group.by = 'Cell_type', secreted_only = T, trim = 0.1)
full_ctrl_cell_type.list <- list(ctrl_veh_cell_type.cch, ctrl_inh_cell_type.cch)
full_ctrl_cell_type.cch <- mergeCellChat(list(Veh = ctrl_veh_cell_type.cch,
                                              Inh = ctrl_inh_cell_type.cch),
                                         add.names = c('Veh', 'Inh'))

ctrl_veh_cell_state.cch <- DoCellChat(inh_cont_veh.srt, group.by = 'Cell_state', secreted_only = T, trim = 0.1)
ctrl_inh_cell_state.cch <- DoCellChat(inh_cont_inh.srt, group.by = 'Cell_state', secreted_only = T, trim = 0.1)
full_ctrl_cell_state.list <- list(ctrl_veh_cell_state.cch, ctrl_inh_cell_state.cch)
full_ctrl_cell_state.cch <- mergeCellChat(list(Veh = ctrl_veh_cell_state.cch,
                                               Inh = ctrl_inh_cell_state.cch),
                                          add.names = c('Veh', 'Inh'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(full_ctrl_cell_type.cch, 'analysis/PART26.full_ctrl_cell_type.cellchat.rds')
saveRDS(full_ctrl_cell_type.list, 'analysis/PART26.full_ctrl_cell_type.cellchat_list.rds')
saveRDS(full_ctrl_cell_state.cch, 'analysis/PART26.full_ctrl_cell_state.cellchat.rds')
saveRDS(full_ctrl_cell_state.list, 'analysis/PART26.full_ctrl_cell_state.cellchat_list.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## 4. RA LatsCKO Veh vs Csf1r Inhibitor
lats_veh_cell_type.cch <- DoCellChat(inh_lats_veh.srt, group.by = 'Cell_type', secreted_only = T, trim = 0.1)
lats_inh_cell_type.cch <- DoCellChat(inh_lats_inh.srt, group.by = 'Cell_type', secreted_only = T, trim = 0.1)
full_lats_cell_type.list <- list(lats_veh_cell_type.cch, lats_inh_cell_type.cch)
full_lats_cell_type.cch <- mergeCellChat(list(Veh = lats_veh_cell_type.cch,
                                              Inh = lats_inh_cell_type.cch),
                                         add.names = c('Veh', 'Inh'))

lats_veh_cell_state.cch <- DoCellChat(inh_lats_veh.srt, group.by = 'Cell_state', secreted_only = T, trim = 0.1)
lats_inh_cell_state.cch <- DoCellChat(inh_lats_inh.srt, group.by = 'Cell_state', secreted_only = T, trim = 0.1)
full_lats_cell_state.list <- list(lats_veh_cell_state.cch, lats_inh_cell_state.cch)
full_lats_cell_state.cch <- mergeCellChat(list(Veh = lats_veh_cell_state.cch,
                                               Inh = lats_inh_cell_state.cch),
                                          add.names = c('Veh', 'Inh'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(full_lats_cell_type.cch, 'analysis/PART26.full_lats_cell_type.cellchat.rds')
saveRDS(full_lats_cell_type.list, 'analysis/PART26.full_lats_cell_type.cellchat_list.rds')
saveRDS(full_lats_cell_state.cch, 'analysis/PART26.full_lats_cell_state.cellchat.rds')
saveRDS(full_lats_cell_state.list, 'analysis/PART26.full_lats_cell_state.cellchat_list.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## 5. RA Inhibitor Cont vs LatsCKO
inh_cont_cell_type.cch <- DoCellChat(inh_cont_inh.srt, group.by = 'Cell_type', secreted_only = T, trim = 0.1)
inh_lats_cell_type.cch <- DoCellChat(inh_lats_inh.srt, group.by = 'Cell_type', secreted_only = T, trim = 0.1)
full_inh_cell_type.list <- list(inh_cont_cell_type.cch, inh_lats_cell_type.cch)
full_inh_cell_type.cch <- mergeCellChat(list(Control = inh_cont_cell_type.cch,
                                             LatsCKO = inh_lats_cell_type.cch),
                                        add.names = c('Control', 'LatsCKO'))

inh_cont_cell_state.cch <- DoCellChat(inh_cont_inh.srt, group.by = 'Cell_state', secreted_only = T, trim = 0.1)
inh_lats_cell_state.cch <- DoCellChat(inh_lats_inh.srt, group.by = 'Cell_state', secreted_only = T, trim = 0.1)
full_inh_cell_state.list <- list(inh_cont_cell_state.cch, inh_lats_cell_state.cch)
full_inh_cell_state.cch <- mergeCellChat(list(Control = inh_cont_cell_state.cch,
                                              LatsCKO = inh_lats_cell_state.cch),
                                         add.names = c('Control', 'LatsCKO'))

####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
saveRDS(full_inh_cell_type.cch, 'analysis/PART26.full_inh_cell_type.cellchat.rds')
saveRDS(full_inh_cell_type.list, 'analysis/PART26.full_inh_cell_type.cellchat_list.rds')
saveRDS(full_inh_cell_state.cch, 'analysis/PART26.full_inh_cell_state.cellchat.rds')
saveRDS(full_inh_cell_state.cch, 'analysis/PART26.full_inh_cell_state.cellchat_list.rds')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~