####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ver <- '2'
Part <- 'PART01_Data_Collection'

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
studies <- U(sample_meta.df$Study)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Data 1 - 2023_Csf1r_CTsai Multiome  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
study <- '2023_Csf1r_CTsai'
library <- c(
        '2023_csf1r_CTsai_flox_V',
        '2023_csf1r_CTsai_flox_G',
        '2023_csf1r_CTsai_CKO_V',
        '2023_csf1r_CTsai_CKO_G'
)
gh5 <- paste0('/Volumes/shire/data/scmulti/', study, '/matrix/', library, '/outs/filtered_feature_bc_matrix.h5')
afrag <- paste0('/Volumes/shire/data/scmulti/', study, '/matrix/', library, '/outs/atac_fragments.tsv.gz')
merge.srt <- MakeDataset(
        study = study,
        library = library,
        mode = 'multiome',
        matrix_dir = paste0('/Volumes/shire/data/scmulti/', study, '/matrix/', library),
        gex_filtered_h5_file = gh5,
        atac_frag_file = afrag
)
saveRDS(merge.srt, paste0('individual/1.', study, '.raw.srt.rds'))

merge.srt <- MakeDataset(
        study = study,
        library = library,
        mode = 'cellbender',
        matrix_dir = paste0('/Volumes/shire/data/scmulti/', study, '/matrix/cellbender_v1/', library),
        gex_filtered_h5_file = NULL,
        atac_frag_file = NULL
)
saveRDS(merge.srt, paste0('individual/1.', study, '.cbn.srt.rds'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Data 2 - 2023_Latscko_CTsai snRNA  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
study <- '2023_Latscko_CTsai'
library <- c(
        '2023_Latscko_CTsai_G1',
        '2023_Latscko_CTsai_G2',
        '2023_Latscko_CTsai_G3',
        '2023_Latscko_CTsai_G4',
        '2023_Latscko_CTsai_G5',
        '2023_Latscko_CTsai_G6',
        '2023_Latscko_CTsai_G7',
        '2023_Latscko_CTsai_G8'
)

merge.srt <- MakeDataset(
        study = study,
        library = library,
        mode = '10x',
        matrix_dir = paste0('/Volumes/shire/data/scrnaseq/', study, '/matrix/', library)
)
saveRDS(merge.srt, paste0('individual/2.', study, '.raw.srt.rds'))

merge.srt <- MakeDataset(
        study = study,
        library = library,
        mode = 'cellbender',
        matrix_dir = paste0('/Volumes/shire/data/scrnaseq/', study, '/matrix/cellbender_v1/', library)
)
saveRDS(merge.srt, paste0('individual/2.', study, '.cbn.srt.rds'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Data 3 - 2023_Csf1r_CTsai snRNA  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
study <- '2023_Csf1r_CTsai'
library <- c(
        '2023_Csf1r_CTsai_F_GW',
        '2023_Csf1r_CTsai_KO_GW'
)

merge.srt <- MakeDataset(
        study = study,
        library = library,
        mode = '10x',
        matrix_dir = paste0('/Volumes/shire/data/scrnaseq/', study, '/matrix/', library)
)
saveRDS(merge.srt, paste0('individual/3.', study, '.raw.srt.rds'))

merge.srt <- MakeDataset(
        study = study,
        library = library,
        mode = 'cellbender',
        matrix_dir = paste0('/Volumes/shire/data/scrnaseq/', study, '/matrix/cellbender_v1/', library)
)
saveRDS(merge.srt, paste0('individual/3.', study, '.cbn.srt.rds'))
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~----