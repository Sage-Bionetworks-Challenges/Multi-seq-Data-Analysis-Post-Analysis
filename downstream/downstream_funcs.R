aggregate_output_by_mean<- function(data_path_list, new_filename) {

  data_list <- lapply(data_path_list, function(f) {
    fread(f, verbose = FALSE) %>% tibble::column_to_rownames("V1")
  })
  shared_cells <- Reduce(intersect, lapply(data_list, colnames))
  shared_genes <- Reduce(intersect, lapply(data_list, rownames))

  data_list <- lapply(data_list, function(f) {
    f <- f[shared_genes, shared_cells]
    f %>% tibble::rownames_to_column("V1") %>%
      as.data.table() %>%
      melt.data.table(id.vars = 1) %>% 
        setnames(c("genes", "cells", "values"))
  })
      
  new_data <- rbindlist(data_list, use.names = TRUE) %>%
    .[, .(values = mean(values)), by = .(cells, genes)] %>%
    dcast(genes ~ cells, value.var = "values") %>%
    tibble::column_to_rownames("genes")

  fwrite(new_data, new_filename, row.names = TRUE)
  return(new_data)
}

# Prepare seurat object from 10x style matrices
matrices_to_seurat <- function(data_path, nickname, dataset, cells, data_type, condition, npcs = 10){
  # create output directory if missing
  
  if(!file.exists('seurat_objects')){
    dir.create('seurat_objects')
  }
  
  if(!file.exists('plots')){
    dir.create('plots')
  }
  
  # Build output file name
  output_filename <- paste0('seurat_objects/', nickname, '.rds')
  
  # Check if file already exists and/or was recorded in metadata
  # if(file.exists(output_filename)){
  #   print('Error: RDS file already exists with this name')
  # }else{
    # print('Processing Seurat object - ETA ~5min')
    
    # read in matrix and convert to seurat object
    # raw <- Read10X(data.dir = data_folder)
    raw <- data.table::fread(data_path, data.table = FALSE) %>% tibble::column_to_rownames("V1")
    so <- CreateSeuratObject(counts = raw,
                             project = nickname)
    
    # Subset to just cells in reference/10X 'filtered' list
    so <- subset(so,
                 cells = cells)
    
    # create cell_barcode metadata
    
    so@meta.data$cell_barcode <- row.names(so@meta.data)

    # Process seurat object with standard workflow
    so <- NormalizeData(so, verbose = FALSE) %>%
      FindVariableFeatures(nfeatures = 3000, verbose = FALSE) %>%
      ScaleData(verbose = FALSE) %>%
      RunPCA(verbose = FALSE)
    suppressWarnings(
      so <- so %>% RunUMAP(reduction = 'pca', dims = 1:npcs, verbose = FALSE)
    )
    so <- so %>%
      RunUMAP(reduction = 'pca', dims = 1:npcs, verbose = FALSE) %>%
      FindNeighbors(reduction = 'pca', dims = 1:npcs, verbose = FALSE) %>%
      FindClusters(resolution = seq(from = 0.1, to = 1, by = 0.1), verbose = FALSE)
    
    # Save umap colored by nCount_UMI
    # curr_umap <- FeaturePlot(so, 
    #                          features = 'nCount_RNA')+
    #   ggtitle(nickname)+
    #   coord_equal()
    
    # ggsave(filename = paste0('plots/', nickname, '_umap_by_counts.pdf'),
    #        plot = curr_umap)
    
    # Save metadata for seurat object
    
    meta <- tibble(file = data_path,
                   nickname = nickname,
                   dataset = dataset,
                   data_type = data_type,
                   condition = condition,
                   date_created = Sys.time(),
                   output_file = output_filename,
                   npcs = npcs)
    
    
    if(!file.exists('so_metadata.csv')){
      # Save seurat object & metadata
      saveRDS(so, file = output_filename)
      write_csv(meta, 'so_metadata.csv')
    }else{
      prior_meta <- read_csv('so_metadata.csv')
      
      if(meta$nickname %in% prior_meta$nickname){
        prior_meta <- prior_meta[-which(prior_meta$nickname == meta$nickname), ]
      }

      new_meta <- rbind(prior_meta,
                        meta)
      
      # Save seurat object & metadata
      saveRDS(so, file = output_filename)
      write_csv(new_meta, 'so_metadata.csv')
    
    return(so)
  }
  
}



compare_so_to_ref <- function(query,
                              reference){
  # match genes/cells
  shared_genes <-  intersect(rownames(query), rownames(reference))
  shared_cells <-  intersect(colnames(query), colnames(reference))
  reference@assays$RNA@data <- reference@assays$RNA@data[shared_genes, shared_cells]
  reference@assays$RNA@counts <- reference@assays$RNA@counts[shared_genes, shared_cells]
  query@assays$RNA@data <- query@assays$RNA@data[shared_genes, shared_cells]
  query@assays$RNA@counts <- query@assays$RNA@counts[shared_genes, shared_cells]


  if(sum(colnames(query) == colnames(reference))/ncol(reference) != 1){
    print('Error - cell barcodes do not line up between query and reference data set')
  }
  
  # Clustering similarity (NMI and ARI)
  paired_cluster_resolutions <- intersect(colnames(query@meta.data),
                                          colnames(reference@meta.data)) %>%
    grep(., pattern = '^RNA_snn_res.', value = TRUE)
  
  
  for(i in paired_cluster_resolutions){
    curr_nmi <- aricode::NMI(c1 = query@meta.data[query@meta.data$cell_barcode %in% shared_cells, colnames(query@meta.data)==i],
                             c2 = reference@meta.data[reference@meta.data$cell_barcode %in% shared_cells, colnames(reference@meta.data) == i])
    
    curr_ari <- aricode::ARI(c1 = query@meta.data[query@meta.data$cell_barcode %in% shared_cells, colnames(query@meta.data)==i],
                             c2 = reference@meta.data[reference@meta.data$cell_barcode %in% shared_cells, colnames(reference@meta.data) == i])
    
    clust_tibble <- tibble(query_name = query@project.name,
                           reference_name = reference@project.name,
                           resolution = as.numeric(str_remove(i, pattern = '^RNA_snn_res.')),
                           ari = curr_ari,
                           nmi = curr_nmi,
                           n_clusters.query = length(unique(query@meta.data[query@meta.data$cell_barcode %in% shared_cells, colnames(query@meta.data)==i])),
                           n_clusters.ref = length(unique(reference@meta.data[query@meta.data$cell_barcode %in% shared_cells, colnames(reference@meta.data) == i])))
    
    if(i == paired_cluster_resolutions[1]){
      cluster_metrics <- clust_tibble
    }else{
      cluster_metrics <- rbind(cluster_metrics,
                               clust_tibble)
    }
  }
  
  # KNN overlap
  unique_barcodes <- colnames(reference)
  for(i in unique_barcodes){
    shared_nn <- sum(reference@graphs$RNA_nn[i,] & query@graphs$RNA_nn[i,])
    union_nn <- sum(reference@graphs$RNA_nn[i,] | query@graphs$RNA_nn[i,])
    
    curr_knn_tibble <- tibble(query_name = query@project.name,
                              reference_name = reference@project.name,
                              barcode = i,
                              shared_nn = shared_nn,
                              union_nn = union_nn,
                              jaccard_shared = shared_nn/union_nn)
    
    if(i == unique_barcodes[1]){
      knn_tibble <- curr_knn_tibble
    }else{
      knn_tibble <- rbind(knn_tibble,
                          curr_knn_tibble)
    }
  }
  
  knn_summary <- tibble(query_name = query@project.name,
                        reference_name = reference@project.name,
                        knn_jaccard_mean = mean(knn_tibble$jaccard_shared))
  
  # Drop out comparison
  
  tp <- sum(reference@assays$RNA@counts == 0 & query@assays$RNA@counts == 0)
  fn <- sum(reference@assays$RNA@counts == 0 & query@assays$RNA@counts != 0)
  fp <- sum(reference@assays$RNA@counts != 0 & query@assays$RNA@counts == 0)
  tn <- sum(reference@assays$RNA@counts != 0 & query@assays$RNA@counts != 0)
  
  f1 <- 2*tp/(2*tp+fp+fn)
  
  drop_outs <- tibble(query_name = query@project.name,
                      reference_name = reference@project.name,
                      tp = tp,
                      fp = fp,
                      tn = tn,
                      fn = fn,
                      f1 = f1)
  
  metrics <- list(cluster_metrics = cluster_metrics,
                  knn_tibble = knn_tibble,
                  knn_summary = knn_summary,
                  drop_outs = drop_outs)
  
  return(metrics)
}



lfc_comparison <- function(comparsion){
  
  perturb_nickname = comparsion[1]
  control_nickname = comparsion[2]
  
  print(str_glue("Running DE analysis on {toString(perturb_nickname)} vs {toString(control_nickname)}"))
  so_prtb <- readRDS(file.path("seurat_objects",  str_glue("{perturb_nickname}.rds")))
  so_ctrl <- readRDS(file.path("seurat_objects",  str_glue("{control_nickname}.rds")))
  
  # Add condition metadata
  so_prtb@meta.data$condition <- 'perturbation'
  so_ctrl@meta.data$condition <- 'control'
  
  # Merge to a single seurat object
  so_merge <- merge(x = so_prtb,
                    y = so_ctrl)
  
  Idents(so_merge) <- 'condition'
  
  # Normalize merged seurat object
  so_merge <- NormalizeData(so_merge, verbose = FALSE)

  # Find markers with default wilcox test
  lfc <- FindMarkers(so_merge,
                      ident.1 = 'perturbation',
                      ident.2 = 'control',
                      verbose = FALSE) %>%
    mutate(gene = row.names(.)) %>%
    mutate(ident.1 = perturb_nickname) %>%
    mutate(ident.2 = control_nickname)
  
  return(lfc)
}
