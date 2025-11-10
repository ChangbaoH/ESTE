##### scRNA_Analysis

# Note on Example Code
# The example code provided here has been reorganized for clarity and testing purposes.
# Please review it carefully and modify as necessary for your specific use case.

stop()

temp <- list.files("./data/Plot/figure6/numbat/")
sample_Array <- as.data.frame(array(NA, dim = c(100, 2)))
index_temp <- 1
for (i in 1:length(temp)) {
  temp1 <- list.files(paste("./data/Plot/figure6/numbat/", temp[i], sep = ""))
  for (j in temp1) {
    sample_Array[index_temp,1] <- temp[i]
    sample_Array[index_temp,2] <- j
    index_temp <- index_temp + 1
  }
}
sample_Array <- na.omit(sample_Array)

## clone_seg
{
  save_file_temp <- "./data/Plot/figure6/clone_seg/"
  files <- sample_Array

  temp1 <- mclapply(1:dim(files)[1], function(sample_index){
    cancer <- files[sample_index,1]
    sample <- files[sample_index,2]

    save_file <- paste("./data/Plot/figure6/numbat/",cancer,"/",sample,"/",sep = "")
    nb <- try(Numbat$new(out_dir = paste("./data/Plot/figure6/numbat/",cancer,"/",sample,"/",sep = "") ), silent = T)
    if(any(class(nb) == "try-error")){
      return(NULL)
    }

    temp <- nb$clone_post[,2:3]
    temp <- as.data.frame(unique(temp))
    temp <- temp[order(temp[,1]),]
    cluster_label <- c()
    segments <- c()
    for (i in 2:dim(temp)[1]) {
      temp1 <- setdiff(strsplit(temp[i,2],"[,]")[[1]], segments)
      cluster_label <- c(cluster_label, rep(i-1,length(temp1)))
      segments <- c(segments, temp1)
    }
    y_breaks <- seq(0, 1, length.out = length(segments) + 1)

    colors <- c("#B09C85FF","#74512D","#91D1C2FF","#BCBD22FF","#FEB941","#F39C12FF","#F39B7FFF",
                "#E64B35FF","#BB0021FF","#8E44ADFF","#AD88C6","#2A6EBBFF","#0E46A3","#3C5488FF")
    if(length(unique(cluster_label)) <= 7){
      color_index <- (1:length(unique(cluster_label)))*2-1
    }else{
      color_index <- (1:7)*2-1
      temp <- length(unique(cluster_label)) - 7
      color_index <- c(color_index, (1:temp)*2)
      color_index <- sort(color_index)
    }
    colors <- colors[color_index]
    seg_df <- data.frame(
      segment = segments,
      cluster = cluster_label,
      ymin = y_breaks[-length(y_breaks)],
      ymax = y_breaks[-1]
    )
    seg_df$cluster <- factor(seg_df$cluster, unique(seg_df$cluster))
    seg_df_temp <- seg_df[,2:4]
    seg_df_temp$"y_ticks" <- 0
    for (i in 1:dim(seg_df_temp)[1]) {
      seg_df_temp[i,4] <- mean(c(min(seg_df_temp[which(seg_df_temp$cluster == seg_df_temp[i,1]),2]),
                                 max(seg_df_temp[which(seg_df_temp$cluster == seg_df_temp[i,1]),3])))
    }
    seg_df_temp <- seg_df_temp[,c(1,4)]
    seg_df_temp <- unique(seg_df_temp)
    seg_df_temp <- rbind(seg_df_temp[1,],seg_df_temp)
    seg_df_temp[,1] <- as.numeric(seg_df_temp$cluster)
    seg_df_temp[1,] <- c(0,0)
    seg_df_temp[,1] <- seg_df_temp$cluster + 1
    seg_df_temp[,1] <- paste("clone",seg_df_temp[,1],sep = "")

    p_left <- ggplot(seg_df, aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax, fill = cluster)) +
      geom_rect(color = NA) +
      scale_fill_manual(values = colors) +
      scale_x_continuous(name   = NULL,breaks = c(0.5),labels = "NULL",limits = c(0, 1),expand = c(0,0)) +
      # scale_y_continuous(name   = NULL,breaks = seg_df_temp$y_ticks,labels = seg_df_temp$cluster,limits = c(0, 1),expand = c(0,0)) +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      geom_text(data = seg_df,
                aes(x = 0.5,y = (ymin + ymax) / 2,label = segment),color = "white",size = 6,family = "Arial"
      )+
      theme_bw(base_size = 30) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=30, color = "black", family = "Arial"),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(t = 0.2, r = 0, b = 0.2, l = 0.2, unit = "cm"),
        panel.border = element_blank(),
        legend.position = "none"
      )


    joint_post <- nb$joint_post %>% as.data.frame()
    joint_post <- joint_post[,c(1,2,3,4,64,68,69)]
    joint_post <- joint_post[which(joint_post[,3] %in% segments),]
    sc_geno <- as.data.frame(array(0, dim = c(length(unique(joint_post[,1])),length(unique(joint_post[,3])))))
    rownames(sc_geno) <- unique(joint_post[,1])
    colnames(sc_geno) <- unique(joint_post[,3])
    joint_post <- joint_post[which(joint_post[,7] != "neu"),]
    joint_post <- joint_post[which(joint_post[,4] == joint_post[,7]),]
    for (i in 1:dim(joint_post)[1]) {
      sc_geno[joint_post[i,1],joint_post[i,3]] <- 1
    }

    sc_geno <- sc_geno[,seg_df$segment,drop = F]
    for (i in 1:dim(sc_geno)[2]) {
      sc_geno <- sc_geno[order(sc_geno[,i],decreasing = T),,drop = F]
    }
    sc_geno$"num" <- rowSums(sc_geno)
    sc_geno <- sc_geno[order(sc_geno$num, decreasing = F),]
    temp <- as.data.frame(nb$clone_post[,1:2])
    rownames(temp) <- temp[,1]
    temp <- temp[rownames(sc_geno),]
    sc_geno$"clone" <- temp[,2]
    sc_geno <- sc_geno[order(sc_geno$clone, decreasing = F),]

    x_ticks <-  sc_geno
    temp <- c()
    pos_lines_cell <- c()
    for (i in unique(x_ticks$clone)) {
      if(i == 1){
        temp <- floor(max(which(x_ticks$clone == i))/2)
      }else{
        temp <- c(temp, floor((max(which(x_ticks$clone == i)) - max(which(x_ticks$clone == i-1)))/2) +
                    max(which(x_ticks$clone == i-1)))
        pos_lines_cell <- c(pos_lines_cell, min(which(x_ticks$clone == i)))
      }
    }
    x_ticks <- x_ticks[temp,]
    x_ticks$clone <- paste("clone",x_ticks$clone,sep = "")

    sc_geno <- sc_geno[,seg_df$segment,drop = F]

    data_df <- as.data.frame(sc_geno, stringsAsFactors = FALSE)
    data_df$cell <- rownames(data_df)
    cells <- data_df$cell

    long_list <- lapply(segments, function(seg) {data.frame(cell = cells,segment = seg, value   = sc_geno[, seg])})
    data_df <- do.call(rbind, long_list)
    rownames(data_df) <- NULL

    data_df <- merge(data_df,by = "segment",
                     data.frame(segment = segments,ymin = y_breaks[-length(y_breaks)],ymax = y_breaks[-1],stringsAsFactors = FALSE),)
    data_df$y_mid <- (data_df$ymin + data_df$ymax) / 2
    data_df$tile_height <- data_df$ymax - data_df$ymin

    data_df$cell_f <- factor(data_df$cell, levels = cells)
    data_df$value_f <- factor(data_df$value, levels = c(0,1))

    p_right <- ggplot() +
      geom_tile(data = data_df,
                aes(x      = cell_f,y      = y_mid,fill   = value_f,width  = 1.0,height = tile_height),color = NA) +
      scale_fill_manual(values = c("0" = "white", "1" = "#FFF176"),
      ) +
      scale_x_discrete(name   = NULL,breaks = rownames(x_ticks),labels = x_ticks$clone,expand = c(0,0)) +
      scale_y_continuous(name   = NULL,breaks = (y_breaks[-1] + y_breaks[-length(y_breaks)]) / 2,
                         labels = segments,limits = c(0, 1),expand = c(0,0)) +
      geom_vline(xintercept = pos_lines_cell,color = "red",linetype = "solid",size = 0.5) +
      theme_bw(base_size = 30) +
      theme(axis.title.x = element_blank(),
            axis.title.y = element_text(size=40, color = "black", family = "Arial"),
            panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
            axis.text.x = element_text(size=30, color = "black", family = "Arial"),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            legend.title = element_text(size=40, color = "black", family = "Arial"),
            legend.text = element_text(size=40, color = "black", family = "Arial"),
            strip.text.y = element_text(size=30, color = "black", family = "Arial"),
            # strip.text.y = element_blank(),
            plot.title = element_text(hjust = 0.5, color = "black", family = "Arial"),
            plot.margin = margin(t = 0.2, r = 0, b = 0.2, l = 0.2, unit = "cm"),
            axis.line.x = element_line(color = "black"),
            axis.line.y = element_line(color = "black"),
            panel.border = element_blank(),
            legend.position = "none"
      )

    save_file <- paste("./data/Plot/figure6/clone_seg/",cancer,"/",sample,"/",sep = "")
    if(!dir.exists(save_file)){
      dir.create(save_file,recursive = T)
    }
    saveRDS(sc_geno,paste(save_file, "sc_geno.rds" , sep = ""))

    pdf(paste(save_file, "clone_seg.pdf" , sep = ""), width = 1440/96, height = 960/96)
    p <- plot_grid(p_left, p_right, rel_widths = c(1,18))
    print(p)
    dev.off()
  }, mc.cores = length(sample_Array), mc.cleanup = TRUE)

}

## seurat + numbat
{
  temp1 <- mclapply(1:length(sample_Array), function(sample_index){
    cancer <- files[sample_index,1]
    sample <- files[sample_index,2]

    data_file <- paste("./Data/Cancer/scRNA/",cancer,"/",sample,"/",sep = "")
    temp <- Read10X(paste(data_file,sample,"/outs/filtered_feature_bc_matrix/",sep=""))
    # rownames(temp)[1:5]
    scRNA <- CreateSeuratObject(
      temp,
      project = "scRNA tutorial",
      min.cells = 3,
      min.features = 200)

    scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
    scRNA[["percent.rb"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")
    HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    HB.genes <- CaseMatch(HB.genes, rownames(scRNA))
    scRNA[["percent.hb"]]<-PercentageFeatureSet(scRNA, features=HB.genes)

    qc <- data.frame(
      nFeature = scRNA$nFeature_RNA,
      nCount   = scRNA$nCount_RNA,
      pctMT    = scRNA$percent.mt,
      pctRB    = scRNA$percent.rb,
      pctHB    = scRNA$percent.hb
    )

    # low_q  <- 0.01
    high_q <- 0.99
    minGene <- 200
    # minGene <- quantile(qc$nFeature, low_q)
    maxGene <- quantile(qc$nFeature, high_q)
    maxUMI  <- quantile(qc$nCount,   high_q)
    # pctMT  <- median(qc$pctMT) + 3 * mad(qc$pctMT)
    pctRB  <- median(qc$pctRB) + 3 * mad(qc$pctRB)
    pctHB  <- median(qc$pctHB) + 3 * mad(qc$pctHB)
    # if(pctMT == 0){pctMT <- 10}
    if(pctRB == 0){pctRB <- 10}
    if(pctHB == 0){pctHB <- 10}
    scRNA <- RunMiQC(scRNA, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.75,
                     model.slot = "flexmix_model")
    scRNA <- subset(scRNA, miQC.keep == "keep")
    scRNA <- subset(
      scRNA,
      subset = nFeature_RNA > minGene &
        nFeature_RNA < maxGene &
        nCount_RNA   < maxUMI  &
        # percent.mt   < pctMT   &
        percent.rb   < pctRB   &
        percent.hb   < pctHB
    )
    save_file <- paste("./data/Plot/figure6/seurat/",cancer,"/",sample,"/",sep = "")
    if(!dir.exists(save_file)){
      dir.create(save_file,recursive = T)
    }
    saveRDS(scRNA, paste(save_file,"scRNA.rds",sep = ""))
    counts_mat <- GetAssayData(scRNA, layer = "counts", assay = "RNA")

    saveRDS(counts_mat, paste(save_file,"counts_mat.rds",sep = ""))

    df_allele <- read_tsv(paste(data_file,"mydata/",sample,"/",sample,"_allele_counts.tsv.gz",sep = ""))
    df_allele <- as.data.frame(df_allele)

    save_file <- paste("./data/Plot/figure6/numbat/",cancer,"/",sample,"/",sep = "")
    if(!dir.exists(save_file)){
      dir.create(save_file,recursive = T)
    }

    out = run_numbat(
      counts_mat, # gene x cell integer UMI count matrix
      ref_hca, # reference expression profile, a gene x cell type normalized expression level matrix
      df_allele, # allele dataframe generated by pileup_and_phase script
      genome = "hg38",
      t = 1e-5,
      ncores = 5,
      plot = TRUE,
      out_dir = save_file
    )

  }, mc.cores = length(sample_Array), mc.cleanup = TRUE)

}

## este
{
  temp1 <- mclapply(1:dim(files)[1], function(sample_index){
    cancer <- files[sample_index,1]
    sample <- files[sample_index,2]

    save_file <- paste("./data/Plot/figure6/clone_seg/",cancer,"/",sample,"/",sep = "")
    sc_geno <- readRDS(paste(save_file,"sc_geno.rds",sep = ""))

    save_file <- paste("./data/Plot/figure6/este/",cancer,"/",sample,"/",sep = "")
    if(!dir.exists(save_file)){
      dir.create(save_file,recursive = T)
    }

    isF <- matrix(as.integer(c(1)),nrow = 1, ncol = 1)
    setD <- matrix(as.integer(1),nrow = 1, ncol = 2)
    setD[1,1] <- as.integer(0)
    setD[1,2] <- as.integer(dim(sc_geno)[1]-1)
    eventD <- matrix(as.integer(c(1,dim(sc_geno)[2])), ncol = 2, byrow = T)
    sc_geno <- cbind(rep(as.integer(1),dim(sc_geno)[1]),sc_geno)
    sc_geno <- as.data.frame(sc_geno)
    for (i in 1:dim(sc_geno)[2]) {
      sc_geno[,i] <- as.integer(sc_geno[,i])
    }
    sc_geno <- as.matrix(sc_geno)

    multi_thrds <- 1
    vote_thrds <- 1L
    ncores1 <- 6
    ncores2 <- 8
    poset_np <- min(dim(sc_geno)[2]-1, 8)

    if(!file.exists(paste(save_file,"eps.rds",sep = ""))){
      eps <- NULL
      eps <- estimate_Epsilon_ForMulti(pat = sc_geno, isF = isF, setD = setD, eventD = eventD, multi_thrds = multi_thrds,
                                       threshold1 = as.integer(8), n_p = as.integer(8),
                                       T = 10.0, N_iter = 200L, thrds=as.integer(ncores1))

      saveRDS(eps, paste(save_file,"eps.rds",sep = ""))
    }
    eps <- readRDS(paste(save_file,"eps.rds",sep = ""))

    if(!file.exists(paste(save_file,"este_fit.rds",sep = ""))){
      poset <- find_Poset_ForVote(pat = sc_geno, isF = isF, eps = eps, setD = setD, eventD = eventD,
                                  Fine_Tune_Num = as.integer(2L), vote_size = 5L, vote_threshold = 0.5,
                                  vote_thrds = vote_thrds,
                                  threshold = 0.0001, threshold2 = 0.01, n_p = as.integer(poset_np), is_update_eps = FALSE,
                                  T = 10, N_iter = 200L, thrds = as.integer(ncores2))

      fit <- estimate_Lambda(pat = sc_geno, poset = poset, isF = isF, eps = eps, setD = setD, eventD = eventD,
                             lambdaS = 1.0, L = 100L, sampling = 'add-remove',
                             maxIter = 500L, updateStepSize=20L, tol=0.001, maxLambda=1e6,
                             neighborhoodDist=1L, is_update_eps = FALSE,
                             thrds=as.integer(40))
      fit$lambdaS <- 1.0
      fit$poset <- poset

      saveRDS(fit, paste(save_file,"este_fit.rds",sep = ""))
    }

  }, mc.cores = dim(files)[1], mc.cleanup = TRUE)

  ## time_T
  time_Help_For_SC <- function(genotype, este_fit, cores = 10){

    rateT <- rateTimeHelp(este_fit$lambdaS, este_fit$lambda, este_fit$poset)

    gammaClusterRes <- try(GammaCluster(rateT$LB, c(0.8,0.2)),silent = T)
    if(class(gammaClusterRes) == "try-error"){
      gammaClusterRes <- list()
      gammaClusterRes$responsibilities <- as.data.frame(matrix(1, length(rateT$LB), 1))
    }

    topo_path <- topological_Sort(este_fit$poset)
    parent_Set <- parent_Set_Help(este_fit$poset)

    temp <- mclapply(1:dim(genotype)[1], function(index_G){
      Geno <- unlist(genotype[index_G,])
      posetCluster <- list()
      for (j in length(topo_path):1) {
        e <- topo_path[j]
        if(Geno[e] == 1){
          if(length(posetCluster)==0){
            posetCluster[[1]] <- c(e)
          }else {
            isOk <- FALSE
            for (k in 1:length(posetCluster)) {
              for (l in posetCluster[[k]]) {
                for (d in 1:length(parent_Set)) {
                  if(length(parent_Set[[d]])>0 & e%in%parent_Set[[d]] & l%in%parent_Set[[d]]){
                    posetCluster[[k]] <- c(posetCluster[[k]], topo_path[j])
                    isOk <- TRUE
                    break
                  }
                }
                if(isOk)
                  break
              }
              if(isOk)
                break
            }
            if(!isOk){
              posetCluster[[length(posetCluster)+1]] <- c(topo_path[j])
            }
          }
        }
      }

      age_Rate <- NULL

      if(length(posetCluster) == 0){
        age_Rate <- log(2)*min(rateT$LB)

      }else{
        for (j in 1:length(posetCluster)) {
          posetCluster[[j]] <- which(rateT$LB == max(rateT$LB[posetCluster[[j]]]))
        }

        posetClusterWeight <- sapply(1:length(posetCluster), function(j){
          gammaClusterRes$responsibilities[posetCluster[[j]],1]
        })
        posetClusterWeight <- posetClusterWeight/sum(posetClusterWeight)

        age_Temp <- c()
        for (j in 1:length(posetCluster)) {
          time_L <- rateT$LB[posetCluster[[j]][1]]
          time_U <- rateT$UB[posetCluster[[j]][1]]
          if (time_U != Inf){
            age_Temp <- c(age_Temp, time_L + (time_U - time_L)*log(2))
          }else{
            age_Temp <- c(age_Temp,  time_L + log(2))
          }
        }

        age_Rate <- sum(age_Temp*posetClusterWeight)

      }

      return(age_Rate)
    }, mc.cores = cores, mc.cleanup = TRUE)

    sc_Time <- as.data.frame(array(NA, dim = c(dim(genotype)[1],2)))
    colnames(sc_Time) <- c("cell", "rate_T")
    rownames(sc_Time) <- sc_Time[,1] <- rownames(genotype)
    sc_Time[,2] <- unlist(temp)

    return(sc_Time)

  }

  for (sample_index in 1:dim(files)[1]) {
    cancer <- files[sample_index,1]
    sample <- files[sample_index,2]

    save_file <- paste("./data/Plot/figure6/clone_seg/",cancer,"/",sample,"/",sep = "")
    sc_geno <- readRDS(paste(save_file,"sc_geno.rds",sep = ""))
    Genotype <- sc_geno

    save_file <- paste("./data/Plot/figure6/este/",cancer,"/",sample,"/",sep = "")
    res <- readRDS(paste(save_file,"este_fit.rds",sep = ""))
    eps_temp <- rep(res$epsilon1[1,1], dim(sc_geno)[2])
    Donor_Pair_Genotype <- find_most_Compatible_Genotype_by_Flipping(sc_geno, eps_temp, res$poset, max_iter = 100000, cores = 30)
    rownames(Donor_Pair_Genotype) <- rownames(sc_geno)
    colnames(Donor_Pair_Genotype) <- colnames(sc_geno)

    saveRDS(Donor_Pair_Genotype, paste(save_file,"pair_sc_geno.rds",sep = ""))
    sc_geno_Sampling <- readRDS(paste(save_file,"pair_sc_geno.rds",sep = ""))

    sc_Time <- time_Help_For_SC(sc_geno_Sampling, res, 50)

    sc_Time <- na.omit(sc_Time)
    saveRDS(sc_Time, paste(save_file,"sc_Time.rds",sep = ""))

  }

}

## gene expression time
{
  ## seurat
  temp1 <- mclapply(1:dim(files)[1], function(sample_index){
    cancer <- files[sample_index,1]
    sample <- files[sample_index,2]

    save_file <- paste("./data/Plot/figure6/seurat/",cancer,"/",sample,"/",sep = "")
    scRNA <- readRDS(paste(save_file,"scRNA.rds",sep = ""))

    scRNA <- NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
    scRNA <- FindVariableFeatures(scRNA,selection.method = "vst",nfeatures = 2000)
    all.genes <- rownames(scRNA)
    scRNA <- ScaleData(scRNA, features = all.genes)
    scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
    scRNA <- JackStraw(scRNA, num.replicate = 100, dims = 50)
    scRNA <- ScoreJackStraw(scRNA, dims = 1:50)

    pvals <- scRNA@reductions$pca@jackstraw$overall.p.values
    pc.p <- pvals[,2]
    cutoff_idx <- which(pc.p >= 0.01)[1]
    if (!is.na(cutoff_idx)) {
      n_pcs <- cutoff_idx - 1
    } else {
      n_pcs <- length(pc.p)
    }

    dims_Use <- 1:n_pcs
    scRNA <- FindNeighbors(scRNA, dims = dims_Use)
    scRNA <- FindClusters(scRNA, resolution = 0.5)
    scRNA <- RunUMAP(scRNA, dims = dims_Use)
    scRNA <- RunTSNE(scRNA, dims = dims_Use)
    saveRDS(scRNA, paste(save_file,"scRNA_dr.rds",sep = ""))

  }, mc.cores = dim(files)[1], mc.cleanup = TRUE)

  ## gene exrichment
  temp1 <- mclapply(1:dim(files)[1], function(sample_index){
    cancer <- files[sample_index,1]
    sample <- files[sample_index,2]

    save_file <- paste("./data/Plot/figure6/este/",cancer,"/",sample,"/",sep = "")
    sc_Time <- readRDS(paste(save_file,"sc_Time.rds",sep = ""))

    save_file <- paste("./data/Plot/figure6/seurat/",cancer,"/",sample,"/",sep = "")
    counts_mat <- readRDS(paste(save_file,"counts_mat.rds",sep = ""))
    counts_mat <- as.matrix(counts_mat)

    sc_Time <- sc_Time[intersect(colnames(counts_mat),rownames(sc_Time)),]
    counts_mat <- counts_mat[,rownames(sc_Time)]

    cor_res <- as.data.frame(array(NA, dim = c(dim(counts_mat)[1],2)))
    colnames(cor_res) <- c("R_Value","P_Value")
    rownames(cor_res) <- rownames(counts_mat)
    for (i in 1:dim(counts_mat)[1]) {
      x <- counts_mat[i,]
      index <- which(x>0)
      if(length(index) >= 30){
        x <- x[index]
        y <- sc_Time[names(index),2]
        temp <- cor.test(y,x,method = "spearman")
        cor_res[i,1] <- temp$estimate
        cor_res[i,2] <- temp$p.value
      }
    }

    cor_res <- cbind(rownames(cor_res), cor_res)
    colnames(cor_res)[1] <- "cell"
    cor_res <- na.omit(cor_res)

    save_file <- paste("./data/Plot/figure6/geneExpression/",cancer,"/",sample,"/",sep = "")
    if(!dir.exists(save_file)){
      dir.create(save_file, recursive = T)
    }
    saveRDS(cor_res, paste(save_file,"sc_cor_Res.rds",sep = ""))

    temp1 <- cor_res[cor_res[,2]>=0,]
    temp1 <- temp1[order(temp1[,2],decreasing = T),]
    temp1 <- temp1[!duplicated(temp1[,1]), ]

    temp2 <- cor_res[cor_res[,2]<0,]
    temp2 <- temp2[order(temp2[,2],decreasing = F),]
    temp2 <- temp2[!duplicated(temp2[,1]), ]
    temp2 <- temp2[order(temp2[,2],decreasing = T),]
    temp1 <- rbind(temp1,temp2)

    gene_rank <- as.numeric(temp1[,2])
    names(gene_rank) <- temp1[,1]
    gene_rank <- sort(gene_rank,decreasing = T)

    gmt_dir <- paste("./data/Genom/msigdb_v2024.1.Hs_files_to_download_locally/")
    gmt_files <- list.files(gmt_dir)
    gmt_files <- gmt_files[unlist(sapply(gmt_files, function(x){
      grepl("symbols",x)
    }))]

    for (i in 1:length(gmt_files)) {
      custom_gmt <- read.gmt(paste(gmt_dir,gmt_files[i],sep = ""))
      gsea_result <- NULL
      gsea_result <- GSEA(
        geneList = gene_rank,
        TERM2GENE = custom_gmt,
        verbose = FALSE
      )
      res_cancer_file6 <- paste(save_file,gsub(".gmt","",gmt_files[i]),"/",sep = "")
      if(!dir.exists(res_cancer_file6)){
        dir.create(res_cancer_file6, recursive = T)
      }
      saveRDS(gsea_result, paste(res_cancer_file6,"gsea_result.rds",sep = ""))
    }

  }, mc.cores = dim(files)[1], mc.cleanup = TRUE)
}









