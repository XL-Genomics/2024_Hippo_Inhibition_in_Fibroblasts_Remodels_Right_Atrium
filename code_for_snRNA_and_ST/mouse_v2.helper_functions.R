####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  scMultiomics of Lats1/2 CKO Mouse With Csf1r Inhibitor
####  2024-01-01 by Xiao LI (Texas Heart Institute, US)
####
####  Helper functions for data processing
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Load Libraries  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
suppressMessages(library('hdf5r'))
suppressMessages(library('data.table'))
library('miloR')
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate global variables  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
markers_lvl1 <- c('Tnnt2',  'Nppa',           # CM
                  'Col1a1', 'Tcf21',          # CF
                  'Pecam1',  'Cdh5',          # EC
                  'Prox1',   'Flt4',          # LEC
                  'Acta2',   'Rgs5',          # SMC
                  'Pdgfrb', 'Vtn',            # Pericyte
                  'Upk3b',  'Wt1',            # Epicardium
                  'Npr3',   'H19',            # Endocardium
                  'Plp1',   'Nrn1',           # Schwann
                  'Nrsn1',  'Npy',            # Neuronal
                  'Tenm4',  'Car3',           # Adipocyte
                  'Ptprc',  'H2-D1',          # Immune
                  'Hba-a1', 'Hbb-bh1',        # RBC
                  'Mki67',  'Top2a')          # Mitotic
markers_lvl1_min <- markers_lvl1[seq(1, L(markers_lvl1), 2)]
markers_lvl2_early <- c('Epcam',  'Krt8',     # Mesoderm progen
                        'Nkx2-5', 'Isl1',     # Mesoderm
                        'Foxa2', 'Pax9',      # Endoderm
                        'Sox2',  'Wnt6')      # Ectoderm
markers_lvl2_immune <- c('Adgre1', 'Fcgr1',   # Mf, Mono
                         'Xcr1',   'Cd209a',  # Dc
                         'Cd79a',  'Ms4a1',   # B
                         'Cd3e',   'Nkg7',    # T+NK
                         'S100a9', 'S100a8')  # Granulocyte
markers_lvl2_bm <- c('Ccl5',   # ILC
                     'Vpreb1', # B
                     'Pf4',  # megakaryocytes
                     'Car1',   # erythrocytes
                     'Prss34',    # basophils
                     'Prg2')  # Geosinophils
markers_lvl2_tc <- c('Cd4',
                     'Trdc', #
                     'Cd8a', #
                     'Foxp3',  # Treg
                     'Ifngr1', # Treg, Th1
                     'Slamf6',   # Tfh
                     'Pdcd1',    # Cd8_exhausted
                     'Ccr7', # Cd8_mem (ccr7+Pdcd1+)
                     'Ctla4',
                     'Xcl1',
                     'Gzmb',
                     'Tox')  # Tfh
markers_lvl2_cm <- c('Nr2f1', 'Cav1', # Atrium
                     'Isl1', 'Tcn', # Outflow tract
                     'Myl2', 'Mpped2', # Ventrical
                     'Shox2', 'Pitx2'
)

options(scipen = 999)
options(future.globals.maxSize= 2000*1024^2)
set.seed(505)
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####  Initiate global functions  ####
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
InitiateProject <- function(Machine = 'local', Ver, Part, Catagory = 'mouse', Project_dir, Data_drive){
        if (Machine == 'local'){
                Data_dir <<- paste0('/Volumes/', Data_drive, '/project/', Project_dir, '/')
                File_dir <<- paste0('~/Documents/Bioinformatics/project/', Project_dir, '/')
        } else if (Machine == 'lorien') {
                Data_dir <<- paste0('/moria/', Project_dir, '/')
                File_dir <<- Data_dir
        } else if (Machine == 'bobbyd') {
                Data_dir <<- paste0('~/work/', Project_dir, '/')
                File_dir <<- Data_dir
        } else {stop('Machine not defined')}
        ####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        Rdata_dir <<- paste0(Data_dir, 'rdata/', Catagory, '_v', Ver)
        Docu_dir <<- paste0(File_dir, 'doc/', Catagory, '_v', Ver, '/')
        Plot_dir <<- paste0(File_dir, 'plot/', Catagory, '_v', Ver, '/', Part, '/')
        Meta_dir <<- paste0(File_dir, 'meta/', Catagory, '_v', Ver, '/', Part, '/')
        dir.create(file.path(Rdata_dir), showWarnings = F, recursive = T)
        setwd(Rdata_dir)
        dir.create(file.path('./individual'), showWarnings = F, recursive = T)
        dir.create(file.path('./integrated'), showWarnings = F, recursive = T)
        dir.create(file.path('./analysis'), showWarnings = F, recursive = T)
        dir.create(file.path('./external'), showWarnings = F, recursive = T)
        dir.create(file.path('./tmp'), showWarnings = F, recursive = T)
        dir.create(file.path(Plot_dir), showWarnings = F, recursive = T)
        dir.create(file.path(Meta_dir), showWarnings = F, recursive = T)
        dir.create(file.path(Docu_dir), showWarnings = F, recursive = T)
}

MakeSrt <- function(mode,
                    matrix_dir,
                    study,
                    library,
                    gex_filtered_h5_file = NULL, ## For Multiome only
                    atac_frag_file = NULL ## For Multiome only
                    ) {
        if (mode == '10x') {
                matrix <- paste0(matrix_dir, '/outs/filtered_feature_bc_matrix/')
                srt <- CreateSeuratObject(counts = Read10X(data.dir = matrix),
                                          min.cells = 1, min.features = 1, project = study)
        }
        else if (mode == 'cellbender') {
                matrix <- paste0(matrix_dir, '/cellbender_filtered.h5')
                srt <- CreateSeuratObject(counts = ReadCB_h5(matrix), min.cells = 1, min.features = 1, project = study)
        }
        else if (mode == 'matrix') {
                matrix <- read.table(gzfile(paste0(matrix_dir[1], '.matrix.csv.gz')), header = T, sep = ',')
                srt <- CreateSeuratObject(counts = matrix, min.cells = 1, min.features = 1, project = study)
        }
        else if (mode == 'multiome') {
                # load both modalities
                inputdata.10x <- Read10X_h5(gex_filtered_h5_file)
                # extract ATAC data
                atac_counts <- inputdata.10x$Peaks
                # only use peaks in standard chromosomes
                grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
                grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
                atac_counts <- atac_counts[as.vector(grange.use), ]
                annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
                seqlevelsStyle(annotations) <- 'UCSC'
                genome(annotations) <- "mm10"
                chrom_assay <- CreateChromatinAssay(
                        counts = atac_counts,
                        sep = c(":", "-"),
                        genome = 'mm10',
                        fragments = atac_frag_file,
                        min.cells = 1,
                        min.features = 1,
                        annotation = annotations
                )
                srt <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC")
                # extract RNA data
                rna_counts <- inputdata.10x$`Gene Expression`
                # Create Seurat object
                srt[["RNA"]] <- CreateAssayObject(counts = rna_counts)
                DefaultAssay(srt) <- 'RNA'
        }
        srt$study <- study
        srt$library <- library
        srt$orig.name <- Cells(srt)
        srt$orig.ident <- NULL
        srt <- RenameCells(srt, new.names = paste0(srt$study,
                                                   '_',
                                                   srt$library,
                                                   ':',
                                                   srt$orig.name
                                                   ), for.merge = F)
        return(srt)
}

MakeDataset <- function(study,
                        library,
                        mode,
                        matrix_dir,
                        gex_filtered_h5_file,
                        atac_frag_file
                        ){
        srt.list <- list()
        for(i in 1:L(library)) {
                message('Processing sample: ', i)
                sample_meta_sub.df <- sample_meta.df[sample_meta.df$Study == study &
                                                             sample_meta.df$Name_on_disk == library[i], ]
                message('Sample metadata found')
                srt.list[[i]] <- MakeSrt(mode = mode,
                                         matrix_dir = matrix_dir[i],
                                         study = study,
                                         library = sample_meta_sub.df$Name_on_disk,
                                         gex_filtered_h5_file = gex_filtered_h5_file[i],
                                         atac_frag_file = atac_frag_file[i]
                )
                print('Seurat generated...')
        }
        if(L(srt.list) > 1) {
                merge.srt <- merge(srt.list[[1]], srt.list[2:L(srt.list)])
        } else {
                merge.srt <- srt.list[[1]]
        }
        return(merge.srt)
}



ProcessSrt_std <- function(srt_obj, var.toal = 0.75, assay = 'RNA', do.umap = T, npcs = 50, ...) {
        srt.out <- srt_obj %>%
                NormalizeData() %>%
                CellCycleScoring(s.features = s.genes,
                                 g2m.features = g2m.genes,
                                 set.ident = F) %>%
                PercentageFeatureSet(pattern = '^mt-', col.name = paste0('pct_mito_', assay), assay = assay)
        srt.out@meta.data[, paste0('pct_mito_', assay)][is.nan(srt.out@meta.data[, paste0('pct_mito_', assay)])] <- 0
        srt.out <- srt.out %>%
                FindVariableFeatures() %>%
                ScaleData(vars.to.regress = c(paste0( c('nFeature_', 'pct_mito_'), assay), 'S.Score', 'G2M.Score'),
                          ...) %>%
                RunPCA(verbose = F, seed.use = 505, npcs = npcs)
        dimN <- FindDimNumber(srt_obj = srt.out, var.toal = var.toal, reduction = 'pca')
        if(do.umap){srt.out <- RunUMAP(srt.out, dims = 1:dimN, seed.use = 505)}
        return(srt.out)
}

ProcessSrt_sct <- function(srt_obj, var.toal = 0.75, do.umap = T, assay = 'RNA') {
        srt.out <- srt_obj %>%
                NormalizeData() %>%
                CellCycleScoring(s.features = s.genes,
                                 g2m.features = g2m.genes,
                                 set.ident = F) %>%
                PercentageFeatureSet(pattern = '^mt-', col.name = paste0('pct_mito_', assay), assay = assay)
        srt.out@meta.data[, paste0('pct_mito_', assay)][is.nan(srt.out@meta.data[, paste0('pct_mito_', assay)])] <- 0
        srt.out <- SCTransform(srt.out,
                               assay = assay,
                               method = "glmGamPoi",
                               seed.use = 505,
                               return.only.var.genes = F,
                               vars.to.regress = c(paste0( c('nFeature_', 'pct_mito_'), assay), 'S.Score', 'G2M.Score')
                               ) %>%
                RunPCA(verbose = F, seed.use = 505)
        dimN <- FindDimNumber(srt_obj = srt.out, var.toal = var.toal, reduction = 'pca')
        if(do.umap){srt.out <- RunUMAP(srt.out, dims = 1:dimN, seed.use = 505)}
        return(srt.out)
}

ProcessSrt_hmn <- function(srt_obj, haromize.by = 'sample', harmony_max_iter = 30, var.toal = 0.75, assay = 'RNA',
                           do.umap = T, harmony_early_stop = T, ...) {
        if(harmony_early_stop){
                srt.out <- srt_obj %>%
                        RunHarmony(group.by.vars = haromize.by,
                                   max.iter.harmony = harmony_max_iter,
                                   plot_convergence = T,
                                   assay.use = assay, ...)
        } else {
                srt.out <- srt_obj %>%
                        RunHarmony(group.by.vars = haromize.by,
                                   epsilon.cluster = -Inf, epsilon.harmony = -Inf,
                                   plot_convergence = T,
                                   max.iter.harmony = harmony_max_iter,
                                   assay.use = assay, ...)
        }
        if(do.umap){
                dimN <- FindDimNumber(srt_obj = srt.out, var.toal = var.toal, reduction = 'harmony')
                srt.out <- RunUMAP(srt.out, dims = 1:dimN, seed.use = 505, reduction = 'harmony',
                                   reduction.name = 'hmn_umap', reduction.key = 'hmnumap_', ...)
        }
        return(srt.out)
}

ProcessSrt_clust <- function(srt_obj, resolution = 0.2, reduction = 'harmony', var.toal = 0.75, ...){
        dimN <- FindDimNumber(srt_obj = srt_obj, var.toal = var.toal, reduction = reduction)
        srt.out <- srt_obj %>%
                FindNeighbors(dims = 1:dimN, reduction = reduction, force.recalc = T) %>%
                FindClusters(resolution = resolution, ...)
}

IntersectSrts <- function(ref.srt, que.srt, ref_data.name, que_data.name) {
        ## Identify shared cells
        shared_cells <- intersect(Cells(ref.srt), Cells(que.srt))
        title <- paste0(ref_data.name, ' has: ', ncol(ref.srt), '... ',
                        que_data.name, ' has: ', ncol(que.srt), '... ',
                        'Shared cells: ',        L(shared_cells))
        ref.srt <- ref.srt %>%
                ProcessSrt_std() %>%
                ProcessSrt_hmn()
        p1 <- DimPlot(ref.srt,
                      cells.highlight = shared_cells,
                      cols.highlight = 'grey60', cols = 'red', pt.size = 0.1, sizes.highlight = 0.1, raster = F) +
                labs(title = paste0(ref_data.name, ' has: ', ncol(ref.srt), '... Shared cells: ', L(shared_cells))) +
                theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank(),
                      legend.position = 'bottom')
        que.srt <- que.srt %>%
                ProcessSrt_std() %>%
                ProcessSrt_hmn()
        p2 <- DimPlot(que.srt,
                      cells.highlight = shared_cells,
                      cols.highlight = 'grey60', cols = 'red', pt.size = 0.1, sizes.highlight = 0.1, raster = F) +
                labs(title = paste0(que_data.name, ' has: ', ncol(que.srt), '... Shared cells: ', L(shared_cells))) +
                theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank(),
                      legend.position = 'bottom')
        print(title)
        ## Create new intersect srt from query data
        que.srt <- que.srt[, shared_cells]
        print('New seurat object generated:')
        print(que.srt)

        ## Match meta data
        ref.meta <- ref.srt@meta.data
        que.meta <- ref.meta[Cells(que.srt),
                             c('pub', 'sample', 'orig.name', 'method', 'platform', 'protocol', 'processed',
                               'tissue', 'enrichment', 'preparation', 'sex', 'age', 'genotype_s', 'genotype_l',
                               'condition', 'strain', 'replicate', 'group', 'batch')]
        print(paste0(nrow(que.meta), ' cells found in reference metadata'))
        que.srt@meta.data <- cbind(que.srt@meta.data, que.meta)
        que.srt$processed <- 'CellBender'
        print(head(que.srt@meta.data, n = 3))
        return(list(p1, p2, que.srt))
}

## A Cellbender bug workaround:
ReadCB_h5 <- function(filename, use.names = TRUE, unique.features = TRUE) {
        if (!requireNamespace('hdf5r', quietly = TRUE)) {
                stop("Please install hdf5r to read HDF5 files")
        }
        if (!file.exists(filename)) {
                stop("File not found")
        }
        infile <- hdf5r::H5File$new(filename = filename, mode = 'r')
        genomes <- names(x = infile)
        output <- list()
        if (hdf5r::existsGroup(infile, 'matrix')) {
                # cellranger version 3
                message('CellRanger version 3+ format H5')
                if (use.names) {
                        feature_slot <- 'features/name'
                } else {
                        feature_slot <- 'features/id'
                }
        } else {
                message('CellRanger version 2 format H5')
                if (use.names) {
                        feature_slot <- 'gene_names'
                } else {
                        feature_slot <- 'genes'
                }
        }
        for (genome in genomes) {
                counts <- infile[[paste0(genome, '/data')]]
                indices <- infile[[paste0(genome, '/indices')]]
                indptr <- infile[[paste0(genome, '/indptr')]]
                shp <- infile[[paste0(genome, '/shape')]]
                features <- infile[[paste0(genome, '/', feature_slot)]][]
                barcodes <- infile[[paste0(genome, '/barcodes')]]
                sparse.mat <- sparseMatrix(
                        i = indices[] + 1,
                        p = indptr[],
                        x = as.numeric(x = counts[]),
                        dims = shp[],
                        repr = 'T'
                )
                if (unique.features) {
                        features <- make.unique(names = features)
                }
                rownames(x = sparse.mat) <- features
                colnames(x = sparse.mat) <- barcodes[]
                sparse.mat <- as(object = sparse.mat, Class = "CsparseMatrix")
                # Split v3 multimodal
                if (infile$exists(name = paste0(genome, '/features'))) {
                        types <- infile[[paste0(genome, '/features/feature_type')]][]
                        types.unique <- unique(x = types)
                        if (length(x = types.unique) > 1) {
                                message("Genome ",
                                        genome, "
                                        has multiple modalities, returning a list of matrices for this genome")
                                sparse.mat <- sapply(
                                        X = types.unique,
                                        FUN = function(x) {
                                                return(sparse.mat[which(x = types == x), ])
                                        },
                                        simplify = FALSE,
                                        USE.NAMES = TRUE
                                )
                        }
                }
                output[[genome]] <- sparse.mat
        }
        infile$close_all()
        if (length(x = output) == 1) {
                return(output[[genome]])
        } else{
                return(output)
        }
}

QuickCheck <- function(srt_obj, markers, ...) {
        p1 <- VlnPlot(srt_obj,
                       features = c('nFeature_RNA', 'nCount_RNA', 'pct_mito'),
                       group.by = 'batch', ncol = 3, pt.size = -1) &
                theme(aspect.ratio = 0.5)
        p1 <- wrap_plots(list((
                p1[[1]] + theme(axis.text.x = element_blank())),
                (p1[[2]] + theme(axis.text.x = element_blank())),
                p1[[3]]), ncol = 1)
        p2.1 <- DimPlot(srt_obj, group.by = "batch", reduction = 'hmn_umap', raster = F) +
                labs(x = "UMAP 1", y = "UMAP 2", title = paste0(ncol(srt_obj), " Cells x ", nrow(srt_obj), " Genes")) +
                guides(col = guide_legend(ncol = 1)) +
                theme(aspect.ratio = 1,
                      axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom")
        p2.2 <- DimPlot(srt_obj, group.by = "batch", reduction = 'umap', raster = F) +
                labs(x = "UMAP 1", y = "UMAP 2", title = paste0(ncol(srt_obj), " Cells x ", nrow(srt_obj), " Genes")) +
                guides(col = guide_legend(ncol = 1)) +
                theme(aspect.ratio = 1,
                      axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "bottom")
        p3 <- FeaturePlot2(srt_obj,
                           reduction = 'hmn_umap',
                           features = intersect(markers, rownames(srt_obj)),
                           ncol = ceiling(L(intersect(markers, rownames(srt_obj))) / 4))
        return(list(p2.1, p2.2, p1, p3))
}
SaveH5ad <- function(srt, path, name, assay, raw_count_only = F, h5seurat_keep = F, verbose = T){
        for(i in 1:ncol(srt@meta.data)){srt@meta.data[, i] <- as.vector(srt@meta.data[, i])} ## convert factor to vectors
        srt <- DietSeurat(srt, scale.data = F,
                          assays = assay,
                          dimreducs = names(srt@reductions),
                          graphs = names(srt@graphs))
        if(raw_count_only){srt <- SetAssayData(srt, slot = 'data', new.data = GetAssayData(srt, slot = 'counts'))}
        if(verbose){message('Raw matrix:')
                print(GetAssayData(srt, slot = 'counts')[1:20, 1:10])
                message('Data matrix:')
                print(GetAssayData(srt, slot = 'data')[1:20, 1:10])
                message('Scaled Data matrix:')}
        if(sum(dim(GetAssayData(srt, assay = assay, slot = 'scale.data')))==0){message('No scaled data slot')} else{
                print(GetAssayData(srt, assay = assay, slot = 'scale.data')[1:20, 1:10])}
        SaveH5Seurat(object = srt, filename = paste0(path, '/', name, '.h5Seurat'), overwrite = T, verbose = verbose)
        Convert(paste0(path, '/', name, '.h5Seurat'), dest = "h5ad", assay = assay, overwrite = T, verbose = verbose)
        if(!h5seurat_keep){system(paste0('rm ', path, '/', name, '.h5Seurat'))}
}

Subcluster <- function(plot_num, srt, celltype, dimN, resolution, features = NULL, fast_marker = F,
                       regress = c('nFeature_CBN', 'pct_mito_CBN')){
        ## Seurat processing
        message('\n', 'Processing Seurat...')
        srt <- srt |>
                FindVariableFeatures(verbose = F) |>
                ScaleData(verbose = F, vars.to.regress = regress) |>
                RunPCA(seed.use = 505, verbose = F) |>
                RunHarmony(group.by.vars = 'sample', max.iter.harmony = 50, assay.use = 'CBN', verbose = F) |>
                RunUMAP(reduction = 'harmony', dims = 1:dimN, verbose = F)
        srt <- DietSeurat(srt, dimreducs = c('umap', 'harmony'))
        gc()
        p <- wrap_plots(list(DimPlot2(srt, group.by = 'sample', raster = F, cols = mycol_10),
                             DimPlot2(srt, group.by = 'genotype', raster = F, cols = mycol_10)
                             ), ncol = 2)
        PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.meta'), 30, 20)
        print(p)
        dev.off()
        if(is.null(features)){ features <- markers_lvl1}
        p <- FeaturePlot2(srt, features = U(features), ncol = ceiling(LU(features)^0.5), reduction = 'umap')
        PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.feature'),
                ceiling(LU(features)^0.5)*2,
                floor(LU(features)^0.5)*2)
        print(p)
        dev.off()
        ## Subcluster evaluation
        message('\n', 'Clustering with ', dimN, ' components, at ', resolution, ' resolution...')
        srt <- srt |>
                FindNeighbors(reduction = 'harmony', dims = 1:dimN) |>
                FindClusters(resolution = resolution)
        for(k in resolution){
                Idents(srt) <- paste0('CBN_snn_res.', k)
                levels(srt) <- str_sort(levels(srt), numeric = T)
                if(fast_marker){
                        marker <- FindAllMarkers(srt, assay = 'CBN', only.pos = T, return.thresh = 0.0001, max.cells.per.ident = 5e3)
                } else {
                        marker <- FindAllMarkers(srt, assay = 'CBN', only.pos = T, return.thresh = 0.0001)
                }
                srt@misc$marker[[paste0('CBN_snn_res.', k)]] <- marker
                cols <- mycol_20
                if(LU(srt@active.ident)>20){cols <- mycol_40}
                p <- DimPlot2(srt, label = T, repel = T, raster = F, cols = cols) + labs(title = paste0('Louvain_', k))
                PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.', k, 'res.umap'), 10, 10)
                print(p)
                dev.off()
                p <- MarkerHeatmap(srt, marker.df = marker, top = 20)
                PlotPDF(paste0(plot_num, '.', celltype, '.', dimN, 'dims.', k, 'res.heatmap'),
                        20,
                        L(levels(srt))*20/8)
                print(p)
                dev.off()
                gc()
        }
        return(srt)
}

GetDoublet <- function(srt_obj, doublet_rate, dimN.var.toal){
        ## Scrublet (run via reticulate)
        mtx <- GetAssayData(srt_obj, layer = "counts", assay = 'CBN')
        mtx <- t(mtx)
        scrub_model <- scr$Scrublet(mtx, expected_doublet_rate = doublet_rate)
        rst <- scrub_model$scrub_doublets(min_gene_variability_pctl = dimN.var.toal*100,
                                          n_prin_comps = 30L,
                                          min_counts = 2, min_cells = 3)
        rst[[2]] <- scrub_model$call_doublets(threshold = 0.3) ## adjusted based on histogram
        sc_doublets <- Cells(srt_obj)[rst[[2]]]
        sc_singlets <- Cells(srt_obj)[!rst[[2]]]
        srt_obj$Scrublet_doublet <- 'Singlet'
        srt_obj$Scrublet_doublet[rst[[2]]] <- 'Doublet'
        Scrublet <- rst[[1]]
        names(Scrublet) <- Cells(srt_obj)
        
        p2 <- DimPlotSplit(srt_obj, split_by = 'Scrublet_doublet', split_order = c('Singlet', 'Doublet'),
                           cols.highlight = mycol_14[c(2, 1)], ncol = 2)
        p2[[1]] <- p2[[1]] + labs(title = paste0('Srub Singlet: ', L(sc_singlets), ' Cells'))
        p2[[2]] <- p2[[2]] + labs(title = paste0('Srub Doublet: ', L(sc_doublets), ' Cells'))
        p <- wrap_plots(
                p2[[1]],
                p2[[2]],
                ncol = 2)
        return(list(
                sc_doublets,
                p,
                Scrublet,
                scrub_model
        ))
}

RunMilo <- function(seurat, k = 20, d = 50, alpha = 0.1, prop = 0.2, umap_dim, pc_dim = 'scVI'){
        sce <- as.SingleCellExperiment(seurat)
        milo.meta <- seurat@meta.data
        milo.obj <- Milo(sce)
        reducedDim(milo.obj, "UMAP") <- reducedDim(sce, str_to_upper(umap_dim))
        milo.obj <- buildGraph(milo.obj, k = k, d = d, reduced.dim = str_to_upper(pc_dim))
        milo.obj <- makeNhoods(milo.obj, k = k, d = d, refined = TRUE, prop = prop, reduced_dims = str_to_upper(pc_dim))
        milo.obj <- calcNhoodDistance(milo.obj, d = d, reduced.dim = str_to_upper(pc_dim))
        milo.obj <- countCells(milo.obj, samples = "Group1", meta.data = milo.meta)
        
        milo.design <- as.data.frame(xtabs(~ Group2 + Group1, data = milo.meta))
        milo.design <- milo.design[milo.design$Freq > 0, ]
        milo.design <- distinct(milo.design)
        rownames(milo.design) <- milo.design$Group1
        
        milo.res <- testNhoods(milo.obj, design = ~ Group2, design.df = milo.design, reduced.dim = str_to_upper(pc_dim))
        milo.obj <- buildNhoodGraph(milo.obj)
        return(list(milo.obj, milo.res))
}

## Modified version of Milo plotting functions
library(ggraph)
plotNhoodGraph.my <- function (x, layout = "UMAP", colour_by = NA, subset.nhoods = NULL,
                               edge_alpha = 0.1, edge_color = 'grey90',
                               size_range = c(0.5, 3), node_stroke = 0.1, ...)
{
        nh_graph <- nhoodGraph(x)
        if (!is.null(subset.nhoods)) {
                nh_graph <- igraph::induced_subgraph(nh_graph,
                                                     vids = which(as.numeric(V(nh_graph)$name) %in%
                                                                          unlist(nhoodIndex(x)[subset.nhoods])))
        }
        nh_graph <- permute(nh_graph, order(vertex_attr(nh_graph)$size,
                                            decreasing = TRUE))
        if (is.character(layout)) {
                redDim <- layout
                layout <- reducedDim(x, redDim)[as.numeric(vertex_attr(nh_graph)$name),
                ]
                if (!any(class(layout) %in% c("matrix"))) {
                        warning("Coercing layout to matrix format")
                        layout <- as(layout, "matrix")
                }
        }
        if (!is.na(colour_by)) {
                if (colour_by %in% colnames(colData(x))) {
                        col_vals <- colData(x)[as.numeric(vertex_attr(nh_graph)$name),
                                               colour_by]
                        if (!is.numeric(col_vals)) {
                                col_vals <- as.character(col_vals)
                        }
                        V(nh_graph)$colour_by <- col_vals
                }
                else {
                        stop(colour_by, "is not a column in colData(x)")
                }
        }
        else {
                V(nh_graph)$colour_by <- V(nh_graph)$size
                colour_by <- "Nhood size"
        }
        if (colour_by %in% c("logFC")) {
                plot.g <- simplify(nh_graph)
                pl <- ggraph(simplify(nh_graph), layout = layout) +
                        geom_edge_link0(aes(width = weight), edge_colour = edge_color, edge_alpha = edge_alpha) +
                        geom_node_point(aes(fill = colour_by, size = size), shape = 21, stroke = node_stroke) +
                        scale_size(range = size_range, name = "Nhood size") +
                        scale_edge_width(range = c(0.2, 3), name = "overlap size") +
                        theme_classic(base_size = 14) +
                        theme(axis.line = element_blank(),
                              axis.text = element_blank(),
                              axis.ticks = element_blank(),
                              axis.title = element_blank())
        }
        if (is.numeric(V(nh_graph)$colour_by)) {
                pl <- pl + scale_fill_gradient2(name = colour_by)
        }
        else {
                mycolors <- colorRampPalette(brewer.pal(11, "RdBu"))(length(unique(V(nh_graph)$colour_by)))
                pl <- pl + scale_fill_manual(values = mycolors, name = colour_by, na.value = "white")
        }
        pl
}
plotNhoodGraphDA.my <- function (x, milo_res, alpha = 0.05, res_column = "logFC", ...)
{
        if (is.character(layout)) {
                if (!layout %in% names(reducedDims(x))) {
                        stop(layout, "is not in readucedDim(x) - choose a different layout")
                }
        }
        signif_res <- milo_res
        signif_res[signif_res$SpatialFDR > alpha, res_column] <- 0
        colData(x)[res_column] <- NA
        colData(x)[unlist(nhoodIndex(x)[signif_res$Nhood]), res_column] <- signif_res[, res_column]
        plotNhoodGraph.my(x, colour_by = res_column, ...)
}

DoCellChat <- function(srt, group.by, assay = 'RNA', secreted_only = F, trim = 0.1){
        CellChatDB <- CellChatDB.human
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
