##### pancancer Immune fitness Inference

# Note on Example Code
# The example code provided here has been reorganized for clarity and testing purposes.
# Please review it carefully and modify as necessary for your specific use case.

stop()

## data
{
  library(readr)
  library(readxl)

  ICGC_drivers <- readRDS("./data/ICGC/ICGC_drivers1.rds")
  res_file <- "./data/Results/"

  ################## TCIA TSNADB
  TSNADB_Read <- function(file_T){
    data <- read.table(file_T, fill = T, sep = "\t")
    data[1,20:21] <- c("frequency_in_tissue", "frequency_in_all")
    colnames(data) <- data[1,]
    data <- data[-1,]
    return(data)
  }

  #### TCIA
  TCIA_NEO <- as.data.frame(read_tsv("./data/ImmuneInference/TCIA/TCIA-NeoantigensData.tsv"))
  unique(TCIA_NEO[,2])

  ####
  TSNADB2_SNV_All <- read.table("./data/ImmuneInference/TSNADB/SNV-derived.txt", fill = T, header = T, sep = "\t")
  TSNADB2_SNV_All[,1] <- unlist(sapply(TSNADB2_SNV_All[,3],function(x){
    strsplit(x,"_")[[1]][1]
  }))

  TSNADB2_INDEL_All <- read.table("./data/ImmuneInference/TSNADB/INDEL-derived.txt", fill = T, header = T, sep = "\t")
  TSNADB2_INDEL_All[,1] <- unlist(sapply(TSNADB2_INDEL_All[,3],function(x){
    strsplit(x,"_")[[1]][1]
  }))

  TSNADB2_FUSION_All <- read.table("./data/ImmuneInference/TSNADB/Fusion-derived.txt", fill = T, header = T, sep = "\t")
  TSNADB2_FUSION_All[,1] <- unlist(sapply(TSNADB2_FUSION_All[,3],function(x){
    gsub("--","-",x)
  }))


  ## Cancer Info
  {
    cancer_conver <- list("Breast-AdenoCA" = "Breast",
                          "CNS-GBM" = "Brain",
                          "ColoRect-AdenoCA" = "Colorectal",
                          "Liver-HCC" = "Liver",
                          "Lung-AdenoCA" = "Lung",
                          "Lung-SCC" =  "Lung",
                          "Prost-AdenoCA" = "Prostate" ,
                          "Skin-Melanoma" = "Skin" ,
                          "Uterus-AdenoCA" = "Uterus")

    TCIA_Cancer_array <- list("Breast-AdenoCA" = "BRCA",
                              "CNS-GBM" = "GBM",
                              "ColoRect-AdenoCA" = NULL,
                              "Liver-HCC" = "LIHC",
                              "Lung-AdenoCA" = "LUAD",
                              "Lung-SCC" =  "LUSC",
                              "Prost-AdenoCA" = "PRAD",
                              "Skin-Melanoma" = "SKCM" ,
                              "Uterus-AdenoCA" = "UCEC")

  }

  ## Neoantigens Info
  for (cancer_index in 1:length(ICGC_drivers)) {
    cancer_temp <- names(ICGC_drivers)[cancer_index]
    if(is.null(cancer_conver[[cancer_temp]]) | is.null(TCIA_Cancer_array[[cancer_temp]])){
      next
    }
    cat(cancer_index,"  ",cancer_temp,"\n")

    res_cancer_file <- paste(res_file,cancer_temp,"/Immune/",sep = "")
    if(!dir.exists(res_cancer_file)){
      dir.create(res_cancer_file)
    }

    TCIA_NEO_Temp <- TCIA_NEO[which(TCIA_NEO[,2] == TCIA_Cancer_array[[cancer_temp]]),]
    File1 <- cancer_conver[[cancer_temp]]

    TSNADB2_SNV_All_Temp <- TSNADB2_SNV_All[which(TSNADB2_SNV_All[,2] == File1),]
    TSNADB2_INDEL_All_Temp <- TSNADB2_INDEL_All[which(TSNADB2_INDEL_All[,2] == File1),]
    TSNADB2_FUSION_All_Temp <- TSNADB2_FUSION_All[which(TSNADB2_FUSION_All[,2] == File1),]

    TCGA_neoantigens <- as.data.frame(read_excel("./data/ImmuneInference/1-s2.0-S2211124718303954-mmc6.xlsx",
                                            sheet = "Fusion neoantigens", skip = 1))
    TCGA_neoantigens <- TCGA_neoantigens[TCGA_neoantigens[,2] == TCIA_Cancer_array[[cancer_temp]],]
    TCGA_neoantigens[,1] <- unlist(sapply(TCGA_neoantigens[,1], function(x){
      substr(x,1,12)
    }))

    patient <- unique(c(TCIA_NEO_Temp[,1], TCGA_neoantigens[,1]))

    temp <- mclapply(1:length(patient), function(i){
      patient_T <- patient[i]
      TCIA_NEO_Temp_P <- TCIA_NEO_Temp[which(TCIA_NEO_Temp[,1] == patient_T),,drop = F]
      TCGA_neoantigens_P <- TCGA_neoantigens[which(TCGA_neoantigens[,1] == patient_T),,drop = F]

      event <- c()
      Genes_T <- c()
      Fusion_T <- c()
      if(dim(TCIA_NEO_Temp_P)[1] > 0){
        Genes_T <- unique(TCIA_NEO_Temp_P[,3])
        event <- c(event, Genes_T)
      }
      if(dim(TCGA_neoantigens_P)[1] > 0){
        Fusion_T <- unique(TCGA_neoantigens_P[,4])
        event <- c(event, Fusion_T)
      }

      res <- array(NA, dim = c(length(event), 3))
      rownames(res) <- event
      colnames(res) <- c("peptide_N", "TSNADB2_peptide", "kind")

      for (j in 1:dim(res)[1]) {
        event <- rownames(res)[j]
        if(length(intersect(event, Genes_T)) > 0){
          peptide <- TCIA_NEO_Temp_P[which(TCIA_NEO_Temp_P[,3] == event),4]
          res[j,1] <- length(peptide)
          res[j,2] <- length(unique(TSNADB2_SNV_All_Temp[which(TSNADB2_SNV_All_Temp[,1] == event & TSNADB2_SNV_All_Temp[,5] %in% peptide), 5])) +
            length(unique(TSNADB2_INDEL_All_Temp[which(TSNADB2_INDEL_All_Temp[,1] == event & TSNADB2_INDEL_All_Temp[,5] %in% peptide), 5]))
          res[j,3] <- "Gene"
        }else if(length(intersect(event, Fusion_T)) > 0){
          peptide <- TCGA_neoantigens_P[which(TCGA_neoantigens_P[,4] == event),6]
          res[j,1] <- length(peptide)
          res[j,2] <- length(unique(TCGA_neoantigens_P[which(TCGA_neoantigens_P[,4] == event & TCGA_neoantigens_P[,6] %in% peptide), 6]))
          res[j,3] <- "Fusion"
        }
      }

      res <- as.data.frame(res)
      for (j in 1:2) {
        res[,j] <- as.integer(res[,j])
      }
      return(list(patient = patient_T, res = res))
    }, mc.cores = min(20, length(patient)), mc.cleanup = TRUE)

    saveRDS(temp, paste(res_cancer_file,"NeoRes.rds",sep = ""))
  }

  ## Cancer Rate Info
  for (cancer_index in 1:length(ICGC_drivers)) {
    cancer_temp <- names(ICGC_drivers)[cancer_index]
    if(is.null(cancer_conver[[cancer_temp]]) | is.null(TCIA_Cancer_array[[cancer_temp]])){
      next
    }
    cat(cancer_index,"  ",cancer_temp,"\n")

    res_cancer_file <- paste(res_file,cancer_temp,"/Immune/",sep = "")
    if(!dir.exists(res_cancer_file)){
      dir.create(res_cancer_file)
    }

    res_cancer_file1 <- paste(res_file,cancer_temp,"/time/",sep = "")
    Donor_Pair_Geno_Last <- readRDS(paste(res_cancer_file1, "Donor_Pair_Geno_Last.rds", sep = ""))
    sample_Patient <- array(NA, dim = c(length(Donor_Pair_Geno_Last), 3))
    colnames(sample_Patient) <- c("sample", "Patient", "censored_Rate")
    for (i in 1:length(Donor_Pair_Geno_Last)) {
      if(!is.null(Donor_Pair_Geno_Last[[i]])){
        sample_Patient[i,1] <- names(Donor_Pair_Geno_Last)[i]
        sample_Patient[i,3] <- Donor_Pair_Geno_Last[[i]]$age_Rate_Total
        if(substr(sample_Patient[i,1],1,4) == "TCGA"){
          sample_Patient[i, 2] <- substr(sample_Patient[i,1],1,12)
        }
      }
    }
    sample_Patient <- na.omit(sample_Patient)
    Neo <- readRDS(paste(res_cancer_file,"NeoRes.rds",sep = ""))
    patient_T <- c()
    for (i in 1:length(Neo)) {
      patient_T <- c(patient_T, Neo[[i]]$patient)
    }
    sample_Patient <- sample_Patient[which(sample_Patient[,2] %in% patient_T),]
    sample_Patient <- sample_Patient[,c(2,3)]
    saveRDS(sample_Patient, paste(res_cancer_file,"sample_Patient_Rate.rds",sep = ""))

    rateT <- readRDS(paste(res_cancer_file1, "rateT.rds", sep = ""))
    driver_genes <- ICGC_drivers[[cancer_index]]$driver_gene

    fitness_Inference_Input_Gene <- list()
    min_Rate <- 0.1
    for (kind in c("TSNADB2_peptide")) {
      total_T <- 0
      for (i in 1:length(Neo)) {
        if(Neo[[i]]$patient %in% sample_Patient[,1]){
          temp <- array(NA, dim = c(dim(Neo[[i]]$res)[1], 6))
          colnames(temp) <- c("Patient", "Gene", "censored_Rate", "gene_LB", "gene_real_Rate", "TSNADB2_peptide")
          temp[,1] <- Neo[[i]]$patient
          temp[,2] <- rownames(Neo[[i]]$res)
          temp[,3] <- sample_Patient[which(sample_Patient[,1] == Neo[[i]]$patient), 2]
          temp[,6] <- Neo[[i]]$res[, "TSNADB2_peptide"]

          index <- c()
          for (j in 1:dim(temp)[1]) {
            if(temp[j,2] %in% driver_genes){
              temp[j,"gene_LB"] <- rateT$LB[which(driver_genes == temp[j,2])]
              if(as.numeric(temp[j,"gene_LB"]) < as.numeric(temp[j,"censored_Rate"]) ){
                temp[j,"gene_real_Rate"] <- as.numeric(temp[j,"censored_Rate"]) - as.numeric(temp[j,"gene_LB"])
                index <- c(index, j)
              }
            }else{
              index <- c(index, j)
            }
          }
          temp <- temp[index,,drop = FALSE]
          if( dim(temp)[1] > 0 & length(which(as.numeric(temp[, kind]) > 0))/dim(temp)[1] >= min_Rate){
            if(total_T == 0){
              fitness_Inference_Input_Gene1 <- temp
            }else{
              fitness_Inference_Input_Gene1 <- rbind(fitness_Inference_Input_Gene1, temp)
            }
            total_T <- total_T + 1
          }
        }
      }

      fitness_Inference_Input_Gene[[kind]] <- fitness_Inference_Input_Gene1[,c("Patient", "Gene", "censored_Rate", "gene_LB", "gene_real_Rate", kind)]
    }

    saveRDS(fitness_Inference_Input_Gene, paste(res_cancer_file,"fitness_Inference_Input_Gene.rds",sep = ""))

  }

}


library(parallel)
library(mvtnorm)
library(MCMCpack)
library(MASS)

library(devtools)
load_all()

ICGC_drivers <- readRDS("./data/ICGC/ICGC_drivers1.rds")
res_file <- "./data/Results/"

# stage1 MCMC
{
  cancer_index <- 1
  cancer_temp <- names(ICGC_drivers)[cancer_index]
  cat(cancer_index,"  ",cancer_temp,"\n")
  res_cancer_file <- paste(res_file,cancer_temp,"/Immune/",sep = "")
  if(!dir.exists(res_cancer_file)){
    dir.create(res_cancer_file)
  }
  fitness_Inference_Input_Gene_List <- readRDS(paste(res_cancer_file,"/fitness_Inference_Input_Gene.rds",sep = ""))
  sample_Patient <- readRDS(paste(res_cancer_file,"/sample_Patient_Rate.rds",sep = ""))
  driver_genes <- ICGC_drivers[[cancer_index]]$driver_gene


  for (hh in 1:length(fitness_Inference_Input_Gene_List)) {
    if(names(fitness_Inference_Input_Gene_List)[hh] != "TSNADB2_peptide"){
      next
    }
    res_cancer_file1 <- paste(res_cancer_file,names(fitness_Inference_Input_Gene_List)[hh],"/",sep = "")
    if(!dir.exists(res_cancer_file1)){
      dir.create(res_cancer_file1)
    }
    files_temp <- list.files(res_cancer_file1)
    if(length(files_temp) < 10){
      seed_temp <- sample(1:99999, 10 - length(files_temp))

      fitness_Inference_Input_Gene <- fitness_Inference_Input_Gene_List[[hh]]
      fitness_Inference_Input_Gene <- fitness_Inference_Input_Gene[,c("Patient", "Gene", "gene_LB", names(fitness_Inference_Input_Gene_List)[hh])]
      if(names(fitness_Inference_Input_Gene_List)[hh] != "TSNADB2_peptide"){
        fitness_Inference_Input_Gene <- fitness_Inference_Input_Gene[
          which(unlist(sapply(fitness_Inference_Input_Gene[,2], function(x){
            if(grepl("-",x)){return(TRUE)}else{return(FALSE)}
          })) == FALSE) ,]
      }

      if(cancer_temp %in% c("CNS-GBM","Liver-HCC","Prost-AdenoCA")){
        min_Patient <- 2
      }else if(cancer_temp %in% c("Breast-AdenoCA","Lung-SCC")){
        min_Patient <- 3
      }else if(cancer_temp %in% c("Lung-AdenoCA","Skin-Melanoma","Uterus-AdenoCA")){
        min_Patient <- 6
      }


      ## stage1 model function
      {
        fitness_Inference_Patient_Gene_Prior <- function(fitness_Inference_Input_Gene,
                                                         sample_Patient,
                                                         driver_genes,
                                                         min_Patient = 4,
                                                         sigma_p_square = 0.0001,
                                                         neighbor_size = 10,
                                                         history_num = 100,
                                                         history_skip = 10,
                                                         history_size = 25000,
                                                         knn_used_in_history_size = 20000,
                                                         gene_history_num = 100,
                                                         gene_history_size = 25000,
                                                         gene_used_in_history_size = 20000,
                                                         max_iter = 40000,
                                                         cellular_clock = "censored_Rate",
                                                         target_recept_rate = 0.234,
                                                         f_b_scaling_update_strength = 10,
                                                         r_update_sigma = 0.01,
                                                         show_size = 1000,
                                                         seed = NULL){

          if(!is.null(seed)){
            seed <- sample(1:99999, 1)
          }
          set.seed(seed)

          fitness_Inference_Input_Gene <- as.data.frame(fitness_Inference_Input_Gene)
          driver_gene_Temp <- intersect(unique(fitness_Inference_Input_Gene[,2]), driver_genes)
          gene_Temp <- setdiff(unique(fitness_Inference_Input_Gene[,2]), driver_gene_Temp)
          num_Temp <- unlist(sapply(1:length(gene_Temp), function(i){
            return(length(which(fitness_Inference_Input_Gene[,2] == gene_Temp[i])))
          }))

          gene_Temp <- gene_Temp[which( unlist(sapply(1:length(gene_Temp), function(i){
            if((grepl("-",gene_Temp[i]) & num_Temp[i] >= 2) | (!grepl("-",gene_Temp[i]) & num_Temp[i] >= min_Patient)){
              return(TRUE)
            }else{return(FALSE)}})) == TRUE) ]
          gene_Temp <- gene_Temp[order(gene_Temp, decreasing = F)]

          driver_gene_Temp_t_min <- c()
          for (i in 1:length(driver_gene_Temp)) {
            temp <- sample_Patient[which(sample_Patient[,1] %in% unique(fitness_Inference_Input_Gene[which(fitness_Inference_Input_Gene[,2] == driver_gene_Temp[i]), 1])),
                                   cellular_clock]
            temp <- min(as.numeric(temp))
            driver_gene_Temp_t_min <- c(driver_gene_Temp_t_min, temp)
            fitness_Inference_Input_Gene[which(fitness_Inference_Input_Gene[,2] == driver_gene_Temp[i]), 3] <-
              as.numeric(fitness_Inference_Input_Gene[which(fitness_Inference_Input_Gene[,2] == driver_gene_Temp[i]), 3])/temp
          }
          names(driver_gene_Temp_t_min) <- driver_gene_Temp

          temp <- which(!is.na(fitness_Inference_Input_Gene[,3]) & as.numeric(fitness_Inference_Input_Gene[,3]) >= 1)
          if(length(temp) > 0){
            fitness_Inference_Input_Gene <- fitness_Inference_Input_Gene[-temp,]
          }
          driver_gene_Temp <- intersect(unique(fitness_Inference_Input_Gene[,2]), driver_genes)
          driver_gene_Temp_t_min <- driver_gene_Temp_t_min[driver_gene_Temp]

          gene_Total <- c(gene_Temp, driver_gene_Temp)

          fitness_Inference_Input_Gene <- fitness_Inference_Input_Gene[which(fitness_Inference_Input_Gene[,2] %in% gene_Total),]
          sample_Patient <- sample_Patient[which(sample_Patient[,1] %in% unique(fitness_Inference_Input_Gene[,1])) ,]
          sample_Patient <- as.data.frame(sample_Patient)
          sample_Patient[,2] <- as.numeric(sample_Patient[,2])

          sample_Patient <- sample_Patient[order(sample_Patient[,cellular_clock], decreasing = F),]
          rownames(sample_Patient) <- 1:dim(sample_Patient)[1]
          fitness_Inference_Input_Gene <- cbind(rep(NA, dim(fitness_Inference_Input_Gene)[1]), fitness_Inference_Input_Gene)
          fitness_Inference_Input_Gene <- as.data.frame(fitness_Inference_Input_Gene)
          for (i in 1:dim(sample_Patient)[1]) {
            fitness_Inference_Input_Gene[which(fitness_Inference_Input_Gene[,2] == sample_Patient[i,1]),1] <- i
          }

          fitness_Inference_Input_Gene[,1] <- as.numeric(fitness_Inference_Input_Gene[,1])
          fitness_Inference_Input_Gene[,4] <- as.numeric(fitness_Inference_Input_Gene[,4])
          fitness_Inference_Input_Gene[,5] <- as.numeric(fitness_Inference_Input_Gene[,5])
          fitness_Inference_Input_Gene <- fitness_Inference_Input_Gene[order(fitness_Inference_Input_Gene[,3], decreasing = F),]
          fitness_Inference_Input_Gene <- fitness_Inference_Input_Gene[order(fitness_Inference_Input_Gene[,1], decreasing = F),]
          fitness_Inference_Input_Gene <- fitness_Inference_Input_Gene[,-1]
          rownames(fitness_Inference_Input_Gene) <- 1:dim(fitness_Inference_Input_Gene)[1]
          rownames(sample_Patient) <- sample_Patient[,1]
          sample_Patient <- sample_Patient[,-1,drop = F]
          scaleNum <- max(sample_Patient[, cellular_clock])
          sample_Patient[,"cellular_clock"] <- sample_Patient[, cellular_clock]/scaleNum
          driver_gene_Temp_t_min <- driver_gene_Temp_t_min/scaleNum
          # for (i in 1:dim(fitness_Inference_Input_Gene)[1]) {
          #   if(!is.na(fitness_Inference_Input_Gene[i,3])){
          #     fitness_Inference_Input_Gene[i,3] <- (fitness_Inference_Input_Gene[i,3]/scaleNum)/sample_Patient[fitness_Inference_Input_Gene[i,1],"cellular_clock"]
          #   }
          # }

          ########  initializaton  #############
          P_Total <- dim(sample_Patient)[1]
          patients <- rownames(sample_Patient)

          ### patient_List fp_bp_Sampling zpk_Sampling
          data_List <- vector("list", dim(sample_Patient)[1])
          names(data_List) <- patients

          fp_bp_Sampling <- array(NA, dim = c(max_iter, P_Total, 2))
          fp_bp_Sampling[1,,] <- 0
          accepted_fp_bp_Sampling <- array(NA, dim = c(max_iter, P_Total))
          accepted_fp_bp_Sampling[1,] <- TRUE
          fp_bp_scaling_Sampling <- array(NA, dim = c(max_iter, P_Total))

          for (i in 1:dim(sample_Patient)[1]) {
            data_List[[i]] <- list()
            ### data & zpk need update
            if(length(which(fitness_Inference_Input_Gene[,1] == patients[i] &
                            fitness_Inference_Input_Gene[,2] %in% gene_Temp)) > 0){
              data_List[[i]]$data_nu <- fitness_Inference_Input_Gene[which(fitness_Inference_Input_Gene[,1] == patients[i] &
                                                                             fitness_Inference_Input_Gene[,2] %in% gene_Temp), c(2,4),drop = F]
              ### zpk
              data_List[[i]]$zpk_nu_Sampling <- array(NA, dim = c(max_iter, dim(data_List[[i]]$data_nu)[1]))
              data_List[[i]]$zpk_nu_Sampling[1,] <- 1
            }else{
              data_List[[i]]$data_nu <- NULL
              data_List[[i]]$zpk_nu_Sampling <- NULL
            }
            ### data & zpk not need update
            if (length(which(fitness_Inference_Input_Gene[,1] == patients[i] &
                             fitness_Inference_Input_Gene[,2] %in% driver_gene_Temp)) > 0){
              data_List[[i]]$data_nnu <- fitness_Inference_Input_Gene[which(fitness_Inference_Input_Gene[,1] == patients[i] &
                                                                              fitness_Inference_Input_Gene[,2] %in% driver_gene_Temp),c(2,3,4),drop = F]
              data_List[[i]]$zpk_nnu_Sampling <- array(NA, dim = c(max_iter, dim(data_List[[i]]$data_nnu)[1]))
              data_List[[i]]$zpk_nnu_Sampling[1,] <- 1
            }else{
              data_List[[i]]$data_nnu <- NULL
              data_List[[i]]$zpk_nnu_Sampling <- NULL
            }
            ### patient neighbors
            half_neighbor_size <- floor(neighbor_size/2)
            start <- max(1, i - half_neighbor_size)
            end <- min(P_Total, start + neighbor_size)
            data_List[[i]]$neighbor <- start:end
            ### fp_bp prior
            data_List[[i]]$prior <- list(mu0 = c(0, 0),
                                         Lambda = matrix(c(1e-4, 0, 0, 1e-4), nrow = 2),
                                         kappa0 = 1,
                                         nu = 1)
            ### mu_fp_bp sigma_fp_bp
            data_List[[i]]$mu_fp_bp_Sampling <- array(NA, dim = c(max_iter, 2))
            data_List[[i]]$mu_fp_bp_Sampling[1,] <- data_List[[i]]$prior$mu0
            data_List[[i]]$sigma_fp_bp_Sampling <- array(NA, dim = c(max_iter, 2, 2))
            data_List[[i]]$sigma_fp_bp_Sampling[1,,] <- data_List[[i]]$prior$Lambda
            ### fp_bp_scaling
            fp_bp_scaling_Sampling[1,i] <- sigma_p_square
            ### fp_bp_accept_rate
            data_List[[i]]$fp_bp_accept_rate <- rep(NA, max_iter)
            data_List[[i]]$fp_bp_accept_rate[1] <- 1
          }

          ### gene_List rk_Sampling scaling mhStepSize
          gene_List <- vector("list", length(gene_Total))
          names(gene_List) <- gene_Total

          for (i in 1:length(gene_Temp)) {
            gene_List[[i]] <- list()
            ### patient
            gene_List[[i]]$patient <- c()
            ### prior
            gene_List[[i]]$prior <- c(0.5, 0.1)
            ### mu nu
            gene_List[[i]]$mu_Sampling <- rep(NA, max_iter)
            gene_List[[i]]$mu_Sampling[1] <- 0.5
            gene_List[[i]]$nu_Sampling <- rep(NA, max_iter)
            gene_List[[i]]$nu_Sampling[1] <- 0.1
          }

          for (i in 1:length(data_List)) {
            data_nu <- data_List[[i]]$data_nu
            if(!is.null(data_nu)){
              for (j in 1:dim(data_nu)[1]) {
                gene <- data_nu[j,1]
                gene_List[[gene]]$patient <- c(gene_List[[gene]]$patient, patients[i])
              }
            }
          }

          for (i in 1:length(gene_Temp)) {
            ### data
            gene_List[[i]]$data <- array(NA, dim = c(max_iter, length(gene_List[[i]]$patient)))
            colnames(gene_List[[i]]$data) <- gene_List[[i]]$patient
            gene_List[[i]]$data[1,] <- 0.5

            gene_List[[i]]$accepted_rpk_Sampling <- array(NA, dim = c(max_iter, length(gene_List[[i]]$patient)))
            colnames(gene_List[[i]]$accepted_rpk_Sampling) <- gene_List[[i]]$patient
            gene_List[[i]]$accepted_rpk_Sampling[1,] <- TRUE

            gene_List[[i]]$t_min <- min(sample_Patient[gene_List[[i]]$patient,"cellular_clock"])
          }

          ######## pi
          f <- function(x, a){x/(1-exp(-x)) - a}
          temp <- array(NA, dim = c(length(gene_Total), 6))
          for (i in 1:length(gene_Total)) {
            temp_g <- fitness_Inference_Input_Gene[which(fitness_Inference_Input_Gene[,"Gene"] == gene_Total[i]),4]
            temp[i,1] <- length(which(temp_g == 0))
            temp[i,2] <- length(temp_g) - temp[i,1]
            temp[i,3] <- temp[i,1]/length(temp_g)

            if(temp[i,2] > 0){
              temp[i,4] <- mean(temp_g[which(temp_g>0)])
              lambda <- uniroot(function(x){f(x, a = temp[i,4])}, c(-1, 1000))$root
              p_0 <- exp(-lambda)
              temp[i,5] <- p_0
              temp[i,6] <- max(0, temp[i,3] - p_0)
              # temp[i,6] <- max(0, (temp[i,3] - p_0)/(1 - p_0))
            }
          }
          pi_temp <- temp[,6]
          pi_temp <- pi_temp[!is.na(pi_temp) & pi_temp > 0]
          mu_Temp <- mean(pi_temp)
          sig_Temp <- var(pi_temp)

          pi_Sampling <- rep(NA, max_iter)
          pi_Sampling[1] <- mu_Temp

          pi_beta <- (mu_Temp*(1 - mu_Temp)/sig_Temp - 1)*(1 - mu_Temp)
          pi_alpha <- pi_beta*mu_Temp/(1 - mu_Temp)

          sigmap <- matrix(c(1,0,0,1), nrow = 2)

          if (history_size <= history_num*history_skip){
            cat("history_size <=", history_num*history_skip, "\n")
            return(NULL)
          }
          if (history_size <= knn_used_in_history_size){
            cat("history_size <=", knn_used_in_history_size , "\n")
            return(NULL)
          }
          history_size_T <- history_num*history_skip
          index_Temp <- unlist(sapply(1:history_num, function(i){
            i*history_skip
          }))
          knn_skip <- floor(knn_used_in_history_size/history_size_T)
          knn_size <- history_size_T*knn_skip


          if (gene_history_size <= gene_used_in_history_size){
            cat("gene_history_size <=", gene_used_in_history_size , "\n")
            return(NULL)
          }
          gene_history_size_T <- gene_history_num*history_skip
          gene_index_Temp <- unlist(sapply(1:gene_history_num, function(i){
            i*history_skip
          }))
          gene_skip <- floor(gene_used_in_history_size/gene_history_size_T)
          gene_size <- gene_history_size_T*gene_skip


          tp_T <- sample_Patient[,"cellular_clock"]

          tp_T <- (tp_T - mean(tp_T))/sd(tp_T)
          ########  iterations  ############
          for (iter in 2:max_iter) {
            if(iter%%show_size == 0){
              cat(iter,"\n")
            }

            temp <- vector("list", P_Total)
            for (p in 1:P_Total) {
              tp <- sample_Patient[p,"cellular_clock"]

              fp_bp <- fp_bp_Sampling[iter - 1,p,]
              accepted_fp_bp <- FALSE
              fp_bp_accept_rate <- NA

              mu_fp_bp <- data_List[[p]]$mu_fp_bp_Sampling[iter - 1,]
              sigma_fp_bp <- data_List[[p]]$sigma_fp_bp_Sampling[iter - 1,,]
              lambda_fp_bp <- solve(sigma_fp_bp)
              scaling <- fp_bp_scaling_Sampling[iter - 1,p]

              gene <- np <- zp <- rp_scaling <- rp <- t_min_arr <- accepted_rp <- c()

              ### update rk
              if(!is.null(data_List[[p]]$data_nu)){
                data_nu <- data_List[[p]]$data_nu
                gene <- data_nu[, 1]
                np <- data_nu[, 2]
                zp <- data_List[[p]]$zpk_nu_Sampling[iter - 1,]
                rp <- accepted_rp <- rep(NA, length(gene))
                t_min_arr <- c()

                for (i in 1:dim(data_nu)[1]) {
                  t_min <- gene_List[[gene[i]]]$t_min
                  t_min_arr <- c(t_min_arr, t_min)
                  rpk <- gene_List[[gene[i]]]$data[iter - 1, patients[p]]
                  rpk_prim <- -1
                  while (rpk_prim<=0 | rpk_prim>=1) {
                    rpk_prim <- rnorm(1, rpk, r_update_sigma)
                  }

                  lambda_prim <- exp(fp_bp[1]*(tp - t_min*rpk_prim) + fp_bp[2])
                  lambda <- exp(fp_bp[1]*(tp - t_min*rpk) + fp_bp[2])

                  tmp_prim <- (zp[i] == 1)*dpois(np[i], lambda_prim, log = T)
                  tmp <- (zp[i] == 1)*dpois(np[i], lambda, log = T)

                  mu_rk <- gene_List[[gene[i]]]$mu_Sampling[iter - 1]
                  nu_rk <- gene_List[[gene[i]]]$nu_Sampling[iter - 1]
                  alpha_rk <- mu_rk*nu_rk
                  beta_rk <- nu_rk - alpha_rk

                  tmp_1 <- (alpha_rk - 1)*log(rpk_prim) + (beta_rk - 1)*log(rpk_prim)
                  tmp_2 <- (alpha_rk - 1)*log(rpk) + (beta_rk - 1)*log(rpk)

                  accept_Rate <- min(1, exp(tmp_prim + tmp_1 - tmp - tmp_2))

                  if (is.null(accept_Rate) | is.na(accept_Rate)){
                    accept_Rate <- 0
                  }

                  if(runif(1,0,1) <= accept_Rate){
                    rp[i] <- rpk_prim
                    accepted_rp[i] <- TRUE
                  }else{
                    rp[i] <- rpk
                    accepted_rp[i] <- FALSE
                  }

                }

              }

              if(!is.null(data_List[[p]]$data_nnu)){
                data_nnu <- data_List[[p]]$data_nnu
                rp <- c(rp, data_nnu[,2])
                t_min_arr <- c(t_min_arr, driver_gene_Temp_t_min[data_nnu[,1]])
                np <- c(np, data_nnu[,3])
                zp <- c(zp, data_List[[p]]$zpk_nnu_Sampling[iter - 1,])
              }

              ### update fp & bp
              fp_bp_prim <- rmvnorm(1, fp_bp, sigmap*scaling)

              lambda_prim <- exp(fp_bp_prim[1]*(tp - t_min_arr*rp) + fp_bp_prim[2])
              lambda <- exp(fp_bp[1]*(tp - t_min_arr*rp) + fp_bp[2])

              tmp_prim <- sum((zp == 1)*dpois(np, lambda_prim, log = T))
              tmp <- sum((zp == 1)*dpois(np, lambda, log = T))

              tmp_1 <- -0.5*(fp_bp_prim - mu_fp_bp)%*%lambda_fp_bp%*%t(fp_bp_prim - mu_fp_bp)
              tmp_2 <- -0.5*(fp_bp - mu_fp_bp)%*%lambda_fp_bp%*%(fp_bp - mu_fp_bp)

              accept_Rate <- min(1, exp(tmp_prim + tmp_1 - tmp - tmp_2))
              if (is.null(accept_Rate) | is.na(accept_Rate)){
                accept_Rate <- 0
              }
              fp_bp_accept_rate <- accept_Rate

              if(runif(1,0,1) <= accept_Rate){
                fp_bp <- fp_bp_prim
                accepted_fp_bp <- TRUE
              }

              scaling <- scaling*exp(f_b_scaling_update_strength*(accept_Rate - target_recept_rate)/(log(iter)))

              ### updata zpk
              pi <- pi_Sampling[iter - 1]
              lambda <- exp(fp_bp[1]*(tp - t_min_arr*rp) + fp_bp[2])
              for (i in 1:length(np)) {
                if(np[i] == 0){
                  r_tmp <- pi/(pi+(1-pi)*dpois(0,lambda[i],log=F))
                  zp[i] <- 1*(runif(1,0,1)>r_tmp)
                }
              }

              ### update pi parameters
              alpha_p <- sum((zp == 0)*(np==0))
              beta_p <- sum(zp == 1)

              temp[[p]] <- list(gene = gene, zp = zp, rp = rp, accepted_rp = accepted_rp, fp_bp_accept_rate = fp_bp_accept_rate,
                                fp_bp = fp_bp, accepted_fp_bp = accepted_fp_bp, scaling = scaling, alpha_p = alpha_p,
                                beta_p = beta_p)

            }

            ### update fp_bp accepted_fp_bp fp_bp_scaling zp rp accepted_rp
            alpha_Total <- 0
            beta_Total <- 0
            for (i in 1:length(temp)) {
              alpha_Total <- alpha_Total + temp[[i]]$alpha_p
              beta_Total <- beta_Total + temp[[i]]$beta_p

              fp_bp_Sampling[iter, i,] <- temp[[i]]$fp_bp
              accepted_fp_bp_Sampling[iter, i] <- temp[[i]]$accepted_fp_bp
              fp_bp_scaling_Sampling[iter, i] <- temp[[i]]$scaling
              data_List[[i]]$fp_bp_accept_rate[iter] <- temp[[i]]$fp_bp_accept_rate

              length_gene <- length(temp[[i]]$gene)

              if(length_gene > 0){
                data_List[[i]]$zpk_nu_Sampling[iter, ] <- temp[[i]]$zp[1:length_gene]
                for (j in 1:length_gene) {
                  gene_List[[temp[[i]]$gene[j]]]$data[iter, patients[i]] <- temp[[i]]$rp[j]
                  gene_List[[temp[[i]]$gene[j]]]$accepted_rpk_Sampling[iter, patients[i]] <- temp[[i]]$accepted_rp[j]
                }
              }

              if(length_gene < length(temp[[i]]$zp)){
                data_List[[i]]$zpk_nnu_Sampling[iter,] <- temp[[i]]$zp[(length_gene + 1):length(temp[[i]]$zp)]
              }
            }

            ### update mu_fp_bp sigma_fp_bp prior
            if (iter >= (history_size + history_size_T - 1)){
              temp <- vector("list", P_Total)
              for (p in 1:P_Total) {
                data_Temp <- fp_bp_Sampling[(iter - history_size_T + 1):iter, p,]
                data_Temp <- data_Temp[index_Temp,]
                par_Temp <- PosteriorDraw_Mvnormal_for_1(data_List[[p]]$prior, data_Temp)
                prior_Temp <- list(mu0 = colMeans(data_Temp),
                                   Lambda = (dim(data_Temp)[1] - 1)*var(data_Temp),
                                   kappa0 = dim(data_Temp)[1],
                                   nu = dim(data_Temp)[1])
                temp[[p]] <- list(par_Temp = par_Temp, prior_Temp = prior_Temp)
              }

            }else if(iter > history_size){
              temp <- vector("list", P_Total)
              for (p in 1:P_Total) {
                data_Temp <- fp_bp_Sampling[history_size:iter, p,]

                data_Temp_knn <- rbind(knn_cache[[p]], data_Temp)
                data_Temp_knn <- data_Temp_knn[(1 + dim(data_Temp)[1]):(dim(data_Temp_knn)[1]),]

                data_Temp <- data_Temp_knn[index_Temp,]
                par_Temp <- PosteriorDraw_Mvnormal_for_1(data_List[[p]]$prior, data_Temp)
                prior_Temp <- list(mu0 = colMeans(data_Temp),
                                   Lambda = (dim(data_Temp)[1] - 1)*var(data_Temp),
                                   kappa0 = dim(data_Temp)[1],
                                   nu = dim(data_Temp)[1])
                temp[[p]] <- list(par_Temp = par_Temp, prior_Temp = prior_Temp)
              }

            }else if(iter == history_size){
              index_Temp1 <- unlist(sapply(1:history_size_T, function(i){
                (i - 1)*knn_skip + 1
              }))

              knn_cache <- list()
              for (p in 1:P_Total) {
                data_Temp <- fp_bp_Sampling[(iter - knn_size + 1):iter, p,]
                knn_cache[[p]] <- data_Temp[index_Temp1,,drop = F]
              }

              temp <- vector("list", P_Total)
              for (p in 1:P_Total) {
                data_Temp <- fp_bp_Sampling[iter, p,,drop = F]

                data_Temp_knn <- rbind(knn_cache[[p]], data_Temp)
                data_Temp_knn <- data_Temp_knn[(1 + dim(data_Temp)[1]):(dim(data_Temp_knn)[1]),]

                data_Temp <- data_Temp_knn[index_Temp,]
                par_Temp <- PosteriorDraw_Mvnormal_for_1(data_List[[p]]$prior, data_Temp)
                prior_Temp <- list(mu0 = colMeans(data_Temp),
                                   Lambda = (dim(data_Temp)[1] - 1)*var(data_Temp),
                                   kappa0 = dim(data_Temp)[1],
                                   nu = dim(data_Temp)[1])

                temp[[p]] <- list(par_Temp = par_Temp, prior_Temp = prior_Temp)
              }

            }else{
              fp_T <- fp_bp_Sampling[iter,,1]
              bp_T <- fp_bp_Sampling[iter,,2]

              fp_T <- (fp_T - mean(fp_T))/sd(fp_T)
              bp_T <- (bp_T - mean(bp_T))/sd(bp_T)

              temp <- vector("list", P_Total)
              for (p in 1:P_Total) {
                dis_T <- ((tp_T - tp_T[p])^2 + (fp_T - fp_T[p])^2 + (bp_T - bp_T[p])^2)
                data_Temp <- fp_bp_Sampling[iter, order(dis_T, decreasing = F)[1:(neighbor_size + 1)],]
                par_Temp <- PosteriorDraw_Mvnormal_for_1(data_List[[p]]$prior, data_Temp)
                prior_Temp <- list(mu0 = colMeans(data_Temp),
                                   Lambda = (dim(data_Temp)[1] - 1)*var(data_Temp),
                                   kappa0 = dim(data_Temp)[1],
                                   nu = dim(data_Temp)[1])
                temp[[p]] <- list(par_Temp = par_Temp, prior_Temp = prior_Temp)
              }

            }

            for (i in 1:P_Total) {
              data_List[[i]]$mu_fp_bp_Sampling[iter,] <- temp[[i]]$par_Temp$mu
              data_List[[i]]$sigma_fp_bp_Sampling[iter,,] <- temp[[i]]$par_Temp$sig
              data_List[[i]]$prior <- temp[[i]]$prior_Temp
            }

            ### update mu_rp sigma_rp prior
            if (iter >= (gene_history_size + gene_history_size_T - 1)){
              temp <- vector("list", length(gene_Temp))
              for (i in 1:length(gene_Temp)) {
                data_Temp <- gene_List[[i]]$data[(iter - gene_history_size_T + 1):iter,,drop = F]
                data_Temp <- data_Temp[gene_index_Temp,,drop = F]
                data_Temp <- as.numeric(data_Temp)

                par_Temp <- MetropolisHastings_Beta(gene_List[[i]]$prior, data_Temp, r_update_sigma,
                                                    c(gene_List[[i]]$mu_Sampling[iter - 1],
                                                      gene_List[[i]]$nu_Sampling[iter - 1]))

                mu_Temp <- mean(data_Temp)
                sig_Temp <- var(data_Temp)

                beta <- (mu_Temp*(1 - mu_Temp)/sig_Temp - 1)*(1 - mu_Temp)
                alpha <- beta*mu_Temp/(1 - mu_Temp)

                prior_Temp <- c(mu_Temp, alpha + beta)

                temp[[i]] <- list(par_Temp = par_Temp, prior_Temp = prior_Temp)
              }

            }else if(iter > gene_history_size){
              temp <- vector("list", length(gene_Temp))
              for (i in 1:length(gene_Temp)) {
                data_Temp <- gene_List[[i]]$data[gene_history_size:iter,,drop = F]

                data_Temp_cache <- rbind(gene_cache[[i]], data_Temp)
                data_Temp_cache <- data_Temp_cache[(1 + dim(data_Temp)[1]):(dim(data_Temp_cache)[1]),]
                data_Temp <- as.numeric(data_Temp_cache[gene_index_Temp,])

                par_Temp <- MetropolisHastings_Beta(gene_List[[i]]$prior, data_Temp, r_update_sigma,
                                                    c(gene_List[[i]]$mu_Sampling[iter - 1],
                                                      gene_List[[i]]$nu_Sampling[iter - 1]))
                mu_Temp <- mean(data_Temp)
                sig_Temp <- var(data_Temp)

                beta <- (mu_Temp*(1 - mu_Temp)/sig_Temp - 1)*(1 - mu_Temp)
                alpha <- beta*mu_Temp/(1 - mu_Temp)

                prior_Temp <- c(mu_Temp, alpha + beta)

                temp[[i]] <- list(par_Temp = par_Temp, prior_Temp = prior_Temp)
              }

            }else if(iter == gene_history_size){
              index_Temp1 <- unlist(sapply(1:gene_history_size_T, function(i){
                (i - 1)*gene_skip + 1
              }))

              gene_cache <- list()
              for (i in 1:length(gene_Temp)) {
                data_Temp <- gene_List[[i]]$data[(iter - gene_size + 1):iter,,drop = F]
                gene_cache[[i]] <- data_Temp[index_Temp1,,drop = F]
              }

              temp <- vector("list", length(gene_Temp))
              for (i in 1:length(gene_Temp)) {
                data_Temp <- gene_List[[i]]$data[iter,,drop = F]

                data_Temp_cache <- rbind(gene_cache[[i]], data_Temp)
                data_Temp_cache <- data_Temp_cache[(1 + dim(data_Temp)[1]):(dim(data_Temp_cache)[1]),]
                data_Temp <- as.numeric(data_Temp_cache[gene_index_Temp,])

                par_Temp <- MetropolisHastings_Beta(gene_List[[i]]$prior, data_Temp, r_update_sigma,
                                                    c(gene_List[[i]]$mu_Sampling[iter - 1],
                                                      gene_List[[i]]$nu_Sampling[iter - 1]))
                mu_Temp <- mean(data_Temp)
                sig_Temp <- var(data_Temp)

                beta <- (mu_Temp*(1 - mu_Temp)/sig_Temp - 1)*(1 - mu_Temp)
                alpha <- beta*mu_Temp/(1 - mu_Temp)

                prior_Temp <- c(mu_Temp, alpha + beta)

                temp[[i]] <- list(par_Temp = par_Temp, prior_Temp = prior_Temp)
              }

            }else{
              temp <- vector("list", length(gene_Temp))
              for (i in 1:length(gene_Temp)) {
                data_Temp <- gene_List[[i]]$data[iter,]

                par_Temp <- MetropolisHastings_Beta(gene_List[[i]]$prior, data_Temp, r_update_sigma,
                                                    c(gene_List[[i]]$mu_Sampling[iter - 1],
                                                      gene_List[[i]]$nu_Sampling[iter - 1]))
                mu_Temp <- mean(data_Temp)
                sig_Temp <- var(data_Temp)

                beta <- (mu_Temp*(1 - mu_Temp)/sig_Temp - 1)*(1 - mu_Temp)
                alpha <- beta*mu_Temp/(1 - mu_Temp)

                prior_Temp <- c(mu_Temp, alpha + beta)

                temp[[i]] <- list(par_Temp = par_Temp, prior_Temp = prior_Temp)
              }
            }

            for (i in 1:length(gene_Temp)) {
              gene_List[[i]]$mu_Sampling[iter] <- temp[[i]]$par_Temp[1]
              gene_List[[i]]$nu_Sampling[iter] <- temp[[i]]$par_Temp[2]
              gene_List[[i]]$prior <- temp[[i]]$prior_Temp
            }

            ### update pi
            pi_Sampling[iter] <- rbeta(1, pi_alpha + alpha_Total, pi_beta + beta_Total)

          }


          return(list(fitness_Inference_Input_Gene = fitness_Inference_Input_Gene,
                      sample_Patient = sample_Patient,
                      gene_Temp = gene_Temp,
                      driver_gene_Temp = driver_gene_Temp,
                      fp_bp_Sampling = fp_bp_Sampling,
                      accepted_fp_bp_Sampling = accepted_fp_bp_Sampling,
                      fp_bp_scaling_Sampling = fp_bp_scaling_Sampling,
                      data_List = data_List,
                      gene_List = gene_List,
                      pi_alpha = pi_alpha,
                      pi_beta = pi_beta,
                      pi_Sampling = pi_Sampling,
                      scaleNum = scaleNum,
                      driver_gene_Temp_t_min = driver_gene_Temp_t_min
          ))

        }

        # You can also use this function by este::fitness_Inference_Patient_Gene_Prior.

      }

      temp <- mclapply(1:length(seed_temp), function(i){
        seed_this <- seed_temp[i]
        set.seed(seed_this)
        res <-  fitness_Inference_Patient_Gene_Prior(fitness_Inference_Input_Gene = fitness_Inference_Input_Gene,
                                                     sample_Patient = sample_Patient,
                                                     driver_genes = driver_genes,
                                                     min_Patient = min_Patient,
                                                     sigma_p_square = 0.0001,
                                                     neighbor_size = 10,
                                                     history_num = 100,
                                                     history_skip = 10,
                                                     history_size = 25000,
                                                     knn_used_in_history_size = 20000,
                                                     gene_history_num = 100,
                                                     gene_history_size = 25000,
                                                     gene_used_in_history_size = 20000,
                                                     max_iter = 40000,
                                                     cellular_clock = "censored_Rate",
                                                     target_recept_rate = 0.234,
                                                     f_b_scaling_update_strength = 10,
                                                     r_update_sigma = 0.01,
                                                     show_size = 1000,
                                                     seed = seed_this)

        saveRDS(res, paste(res_cancer_file1,"res_",seed_this,".rds",sep = ""))
        return(NULL)
      }, mc.cores = length(seed_temp), mc.cleanup = TRUE)

    }
  }
}


# stage2 MCMC
{
  library(parallel)
  library(mvtnorm)
  library(MCMCpack)
  library(MASS)

  ICGC_drivers <- readRDS("./ICGC_drivers1.rds")
  res_file <- "./data/Results/"

  mclapply(1:length(ICGC_drivers), function(cancer_index){
    cancer_temp <- names(ICGC_drivers)[cancer_index]
    if(is.null(cancer_conver[[cancer_temp]]) | is.null(TCIA_Cancer_array[[cancer_temp]])){
      next
    }
    cat(cancer_index,"  ",cancer_temp,"\n")
    res_cancer_file <- paste(res_file,cancer_temp,"/Immune/",sep = "")
    type <- c("TSNADB2_peptide")
    index_temp <- c(1:500)*10

    if(dir.exists(res_cancer_file)){
      for (i in 1:length(type)) {
        res_cancer_file1 <- paste(res_cancer_file,type[i],"/",sep = "")
        files_temp <- list.files(res_cancer_file1)
        prior_data <- list(data_info = list(), data_list = list())
        if(length(files_temp) > 0){
          for (j in 1:length(files_temp)) {
            temp <- readRDS(paste(res_cancer_file1,files_temp[j],sep = ""))

            if(j == 1){
              fitness_Inference_Input_Gene <- temp$fitness_Inference_Input_Gene
              sample_Patient <- temp$sample_Patient
              driver_gene_Temp <- temp$driver_gene_Temp
              gene_Temp <- temp$gene_Temp
              gene_Total <- c(gene_Temp, driver_gene_Temp)
              driver_gene_Temp_t_min <- temp$driver_gene_Temp_t_min
              scaleNum <- temp$scaleNum
              pi_alpha <- temp$pi_alpha
              pi_beta <- temp$pi_beta

              prior_data$data_info <- list(fitness_Inference_Input_Gene = fitness_Inference_Input_Gene, sample_Patient = sample_Patient,
                                           driver_gene_Temp = driver_gene_Temp, gene_Temp = gene_Temp, gene_Total = gene_Total,
                                           driver_gene_Temp_t_min = driver_gene_Temp_t_min, scaleNum = scaleNum,
                                           pi_alpha = pi_alpha, pi_beta = pi_beta)

            }

            f_b_prior_data <- vector("list", length(temp$data_List))
            names(f_b_prior_data) <- names(temp$data_List)
            for (k in 1:length(temp$data_List)) {
              data_temp <- temp$fp_bp_Sampling[35000:40000,k,]
              f_b_prior_data[[k]] <- data_temp[index_temp,]
            }

            gene_prior_data <- vector("list", length(temp$gene_List))
            names(gene_prior_data) <- names(temp$gene_List)
            for (k in 1:length(temp$gene_List)) {
              data_temp <- temp$gene_List[[names(gene_prior_data)[k]]]
              data_temp <- data_temp$data[35000:40000,,drop = F]
              gene_prior_data[[k]] <- data_temp[index_temp,]
            }

            data_temp <- temp$pi_Sampling[35000:40000]
            pi_prior <- data_temp[index_temp]

            prior_data$data_list[[j]] <- list(f_b_prior_data = f_b_prior_data, gene_prior_data = gene_prior_data, pi_prior = pi_prior)

          }

          saveRDS(prior_data, paste(res_cancer_file,type[i],"_prior_data.rds",sep = ""))

        }

      }
    }


    return(NULL)
  }, mc.cores = 10, mc.cleanup = TRUE)

  temp <- as.data.frame(array(NA, dim = c(2000, 3)))
  index <- 1
  for (cancer_index in 1:length(ICGC_drivers)) {
    cancer_temp <- names(ICGC_drivers)[cancer_index]
    if(is.null(cancer_conver[[cancer_temp]]) | is.null(TCIA_Cancer_array[[cancer_temp]])){
      next
    }
    cat(cancer_index,"  ",cancer_temp,"\n")

    res_cancer_file <- paste(res_file,cancer_temp,"/Immune/",sep = "")
    if(!dir.exists(res_cancer_file)){
      dir.create(res_cancer_file)
    }

    res_cancer_file1 <- paste(res_cancer_file,"mcmcRes/",sep = "")
    if(!dir.exists(res_cancer_file1)){
      dir.create(res_cancer_file1)
    }

    type <- c("TSNADB2_peptide")

    for (hh in 1:length(type)) {
      if(file.exists(paste(res_cancer_file,type[hh],"_prior_data.rds",sep = "") )){
        res_cancer_file2 <- paste(res_cancer_file1,type[hh],"/",sep = "")
        if(!dir.exists(res_cancer_file2)){
          dir.create(res_cancer_file2)
        }
        files_temp <- list.files(res_cancer_file2)
        # prior_data <- readRDS(paste(res_cancer_file,type[hh],"_prior_data.rds",sep = ""))

        if(length(files_temp) < 10){
          seed_temp <- sample(1:99999, 10 - length(files_temp))
          for (i in 1:length(seed_temp)) {
            temp[index,1] <- cancer_temp
            temp[index,2] <- type[hh]
            temp[index,3] <- seed_temp[i]
            index <- index + 1
          }
        }
      }
    }
  }
  temp <- na.omit(temp)

  ## stage2 model function
  {
    fitness_Inference_Patient_Gene <- function(prior_data,
                                               sigma_p_square = 0.0001,
                                               max_iter = 10000,
                                               cellular_clock = "censored_Rate",
                                               target_recept_rate = 0.234,
                                               f_b_scaling_update_strength = 10,
                                               r_update_sigma = 0.01,
                                               show_size = 1000,
                                               seed = NULL){

      if(!is.null(seed)){
        seed <- sample(1:99999, 1)
      }
      set.seed(seed)

      fitness_Inference_Input_Gene <- prior_data$data_info$fitness_Inference_Input_Gene
      sample_Patient <- prior_data$data_info$sample_Patient
      driver_gene_Temp <- prior_data$data_info$driver_gene_Temp
      gene_Temp <- prior_data$data_info$gene_Temp
      gene_Total <- prior_data$data_info$gene_Total
      driver_gene_Temp_t_min <- prior_data$data_info$driver_gene_Temp_t_min
      scaleNum <- prior_data$data_info$scaleNum

      ########  initializaton  #############
      P_Total <- dim(sample_Patient)[1]
      patients <- rownames(sample_Patient)

      ### patient_List fp_bp_Sampling zpk_Sampling
      data_List <- vector("list", dim(sample_Patient)[1])
      names(data_List) <- patients

      fp_bp_Sampling <- array(NA, dim = c(max_iter, P_Total, 2))
      accepted_fp_bp_Sampling <- array(NA, dim = c(max_iter, P_Total))
      accepted_fp_bp_Sampling[1,] <- TRUE
      fp_bp_scaling_Sampling <- array(NA, dim = c(max_iter, P_Total))

      for (i in 1:P_Total) {
        temp <- prior_data$data_list[[1]]$f_b_prior_data[[patients[i]]]
        for (j in 2:length(prior_data$data_list)) {
          temp <- rbind(temp, prior_data$data_list[[j]]$f_b_prior_data[[patients[i]]])
        }
        fp_bp_Sampling[1,i,] <- colMeans(temp)
      }

      for (i in 1:dim(sample_Patient)[1]) {
        data_List[[i]] <- list()
        ### data & zpk need update
        if(length(which(fitness_Inference_Input_Gene[,1] == patients[i] &
                        fitness_Inference_Input_Gene[,2] %in% gene_Temp)) > 0){
          data_List[[i]]$data_nu <- fitness_Inference_Input_Gene[which(fitness_Inference_Input_Gene[,1] == patients[i] &
                                                                         fitness_Inference_Input_Gene[,2] %in% gene_Temp), c(2,4),drop = F]
          ### zpk
          data_List[[i]]$zpk_nu_Sampling <- array(NA, dim = c(max_iter, dim(data_List[[i]]$data_nu)[1]))
          data_List[[i]]$zpk_nu_Sampling[1,] <- 1
        }else{
          data_List[[i]]$data_nu <- NULL
          data_List[[i]]$zpk_nu_Sampling <- NULL
        }
        ### data & zpk not need update
        if (length(which(fitness_Inference_Input_Gene[,1] == patients[i] &
                         fitness_Inference_Input_Gene[,2] %in% driver_gene_Temp)) > 0){
          data_List[[i]]$data_nnu <- fitness_Inference_Input_Gene[which(fitness_Inference_Input_Gene[,1] == patients[i] &
                                                                          fitness_Inference_Input_Gene[,2] %in% driver_gene_Temp),c(2,3,4),drop = F]
          data_List[[i]]$zpk_nnu_Sampling <- array(NA, dim = c(max_iter, dim(data_List[[i]]$data_nnu)[1]))
          data_List[[i]]$zpk_nnu_Sampling[1,] <- 1
        }else{
          data_List[[i]]$data_nnu <- NULL
          data_List[[i]]$zpk_nnu_Sampling <- NULL
        }

        ### mu_fp_bp sigma_fp_bp
        Lambda_temp <- var(prior_data$data_list[[1]]$f_b_prior_data[[patients[i]]])
        for (j in 2:length(prior_data$data_list)) {
          Lambda_temp <- Lambda_temp + var(prior_data$data_list[[j]]$f_b_prior_data[[patients[i]]])
        }
        Lambda_temp <- Lambda_temp/length(prior_data$data_list)

        data_List[[i]]$mu_fp_bp_Sampling <- fp_bp_Sampling[1,i,]
        data_List[[i]]$sigma_fp_bp_Sampling <- Lambda_temp

        ### fp_bp_scaling
        fp_bp_scaling_Sampling[1,i] <- sigma_p_square
        ### fp_bp_accept_rate
        data_List[[i]]$fp_bp_accept_rate <- rep(NA, max_iter)
        data_List[[i]]$fp_bp_accept_rate[1] <- 1
      }

      ### gene_List rk_Sampling scaling mhStepSize
      gene_List <- vector("list", length(gene_Total))
      names(gene_List) <- gene_Total

      for (i in 1:length(gene_Temp)) {
        gene_List[[i]] <- list()
        ### patient
        gene_List[[i]]$patient <- c()
        ### alpha beta
        mu_Temp <- mean(as.numeric(prior_data$data_list[[1]]$gene_prior_data[[names(gene_List)[i]]]))
        sig_Temp <- var(as.numeric(prior_data$data_list[[1]]$gene_prior_data[[names(gene_List)[i]]]))
        for (j in 2:length(prior_data$data_list)) {
          mu_Temp <- mu_Temp + mean(as.numeric(prior_data$data_list[[j]]$gene_prior_data[[names(gene_List)[i]]]))
          sig_Temp <- sig_Temp + var(as.numeric(prior_data$data_list[[j]]$gene_prior_data[[names(gene_List)[i]]]))
        }
        mu_Temp <- mu_Temp/length(prior_data$data_list)
        sig_Temp <- sig_Temp/length(prior_data$data_list)

        gene_List[[i]]$beta <- (mu_Temp*(1 - mu_Temp)/sig_Temp - 1)*(1 - mu_Temp)
        gene_List[[i]]$alpha <- gene_List[[i]]$beta*mu_Temp/(1 - mu_Temp)

      }

      for (i in 1:length(data_List)) {
        data_nu <- data_List[[i]]$data_nu
        if(!is.null(data_nu)){
          for (j in 1:dim(data_nu)[1]) {
            gene <- data_nu[j,1]
            gene_List[[gene]]$patient <- c(gene_List[[gene]]$patient, patients[i])
          }
        }
      }

      for (i in 1:length(gene_Temp)) {
        ### data
        gene_List[[i]]$data <- array(NA, dim = c(max_iter, length(gene_List[[i]]$patient)))
        colnames(gene_List[[i]]$data) <- gene_List[[i]]$patient
        mu_Temp <- mean(as.numeric(prior_data$data_list[[1]]$gene_prior_data[[names(gene_List)[i]]]))
        for (j in 2:length(prior_data$data_list)) {
          mu_Temp <- mu_Temp + mean(as.numeric(prior_data$data_list[[j]]$gene_prior_data[[names(gene_List)[i]]]))
        }
        gene_List[[i]]$data[1,] <- mu_Temp/length(prior_data$data_list)

        gene_List[[i]]$accepted_rpk_Sampling <- array(NA, dim = c(max_iter, length(gene_List[[i]]$patient)))
        colnames(gene_List[[i]]$accepted_rpk_Sampling) <- gene_List[[i]]$patient
        gene_List[[i]]$accepted_rpk_Sampling[1,] <- TRUE

        gene_List[[i]]$t_min <- min(sample_Patient[gene_List[[i]]$patient,"cellular_clock"])
      }

      ######## pi
      mu_Temp <- mean(prior_data$data_list[[1]]$pi_prior)
      sig_Temp <- var(prior_data$data_list[[1]]$pi_prior)
      for (j in 2:length(prior_data$data_list)) {
        mu_Temp <- mu_Temp + mean(prior_data$data_list[[j]]$pi_prior)
        sig_Temp <- sig_Temp + var(prior_data$data_list[[j]]$pi_prior)
      }
      mu_Temp <- mu_Temp/length(prior_data$data_list)
      sig_Temp <- sig_Temp/length(prior_data$data_list)

      pi_Sampling <- rep(NA, max_iter)
      pi_Sampling[1] <- mu_Temp

      pi_beta <- (mu_Temp*(1 - mu_Temp)/sig_Temp - 1)*(1 - mu_Temp)
      pi_alpha <- pi_beta*mu_Temp/(1 - mu_Temp)

      sigmap <- matrix(c(1,0,0,1), nrow = 2)

      ########  iterations  ############
      for (iter in 2:max_iter) {
        if(iter%%show_size == 0){
          cat(iter,"\n")
        }

        temp <- vector("list", P_Total)
        for (p in 1:P_Total) {
          tp <- sample_Patient[p,"cellular_clock"]

          fp_bp <- fp_bp_Sampling[iter - 1,p,]
          accepted_fp_bp <- FALSE
          fp_bp_accept_rate <- NA

          mu_fp_bp <- data_List[[p]]$mu_fp_bp_Sampling
          sigma_fp_bp <- data_List[[p]]$sigma_fp_bp_Sampling
          lambda_fp_bp <- solve(sigma_fp_bp)
          scaling <- fp_bp_scaling_Sampling[iter - 1,p]

          gene <- np <- zp <- rp_scaling <- rp <- t_min_arr <- accepted_rp <- c()

          ### update rk
          if(!is.null(data_List[[p]]$data_nu)){
            data_nu <- data_List[[p]]$data_nu
            gene <- data_nu[, 1]
            np <- data_nu[, 2]
            zp <- data_List[[p]]$zpk_nu_Sampling[iter - 1,]
            rp <- accepted_rp <- rep(NA, length(gene))
            t_min_arr <- c()

            for (i in 1:dim(data_nu)[1]) {
              t_min <- gene_List[[gene[i]]]$t_min
              t_min_arr <- c(t_min_arr, t_min)
              rpk <- gene_List[[gene[i]]]$data[iter - 1, patients[p]]
              rpk_prim <- -1
              while (rpk_prim<=0 | rpk_prim>=1) {
                rpk_prim <- rnorm(1, rpk, r_update_sigma)
              }

              lambda_prim <- exp(fp_bp[1]*(tp - t_min*rpk_prim) + fp_bp[2])
              lambda <- exp(fp_bp[1]*(tp - t_min*rpk) + fp_bp[2])

              tmp_prim <- (zp[i] == 1)*dpois(np[i], lambda_prim, log = T)
              tmp <- (zp[i] == 1)*dpois(np[i], lambda, log = T)

              alpha_rk <- gene_List[[gene[i]]]$alpha
              beta_rk <- gene_List[[gene[i]]]$beta

              tmp_1 <- (alpha_rk - 1)*log(rpk_prim) + (beta_rk - 1)*log(rpk_prim)
              tmp_2 <- (alpha_rk - 1)*log(rpk) + (beta_rk - 1)*log(rpk)

              accept_Rate <- min(1, exp(tmp_prim + tmp_1 - tmp - tmp_2))

              if (is.null(accept_Rate) | is.na(accept_Rate)){
                accept_Rate <- 0
              }

              if(runif(1,0,1) <= accept_Rate){
                rp[i] <- rpk_prim
                accepted_rp[i] <- TRUE
              }else{
                rp[i] <- rpk
                accepted_rp[i] <- FALSE
              }

            }

          }

          if(!is.null(data_List[[p]]$data_nnu)){
            data_nnu <- data_List[[p]]$data_nnu
            rp <- c(rp, data_nnu[,2])
            t_min_arr <- c(t_min_arr, driver_gene_Temp_t_min[data_nnu[,1]])
            np <- c(np, data_nnu[,3])
            zp <- c(zp, data_List[[p]]$zpk_nnu_Sampling[iter - 1,])
          }

          ### update fp & bp
          fp_bp_prim <- rmvnorm(1, fp_bp, sigmap*scaling)

          lambda_prim <- exp(fp_bp_prim[1]*(tp - t_min_arr*rp) + fp_bp_prim[2])
          lambda <- exp(fp_bp[1]*(tp - t_min_arr*rp) + fp_bp[2])

          tmp_prim <- sum((zp == 1)*dpois(np, lambda_prim, log = T))
          tmp <- sum((zp == 1)*dpois(np, lambda, log = T))

          tmp_1 <- -0.5*(fp_bp_prim - mu_fp_bp)%*%lambda_fp_bp%*%t(fp_bp_prim - mu_fp_bp)
          tmp_2 <- -0.5*(fp_bp - mu_fp_bp)%*%lambda_fp_bp%*%(fp_bp - mu_fp_bp)

          accept_Rate <- min(1, exp(tmp_prim + tmp_1 - tmp - tmp_2))
          if (is.null(accept_Rate) | is.na(accept_Rate)){
            accept_Rate <- 0
          }
          fp_bp_accept_rate <- accept_Rate

          if(runif(1,0,1) <= accept_Rate){
            fp_bp <- fp_bp_prim
            accepted_fp_bp <- TRUE
          }

          scaling <- scaling*exp(f_b_scaling_update_strength*(accept_Rate - target_recept_rate)/(log(iter)))

          ### updata zpk
          pi <- pi_Sampling[iter - 1]
          lambda <- exp(fp_bp[1]*(tp - t_min_arr*rp) + fp_bp[2])
          for (i in 1:length(np)) {
            if(np[i] == 0){
              r_tmp <- pi/(pi+(1-pi)*dpois(0,lambda[i],log=F))
              zp[i] <- 1*(runif(1,0,1)>r_tmp)
            }
          }

          ### update pi parameters
          alpha_p <- sum((zp == 0)*(np==0))
          beta_p <- sum(zp == 1)

          temp[[p]] <- list(gene = gene, zp = zp, rp = rp, accepted_rp = accepted_rp, fp_bp_accept_rate = fp_bp_accept_rate,
                            fp_bp = fp_bp, accepted_fp_bp = accepted_fp_bp, scaling = scaling, alpha_p = alpha_p,
                            beta_p = beta_p)

        }

        ### update fp_bp accepted_fp_bp fp_bp_scaling zp rp accepted_rp
        alpha_Total <- 0
        beta_Total <- 0
        for (i in 1:length(temp)) {
          alpha_Total <- alpha_Total + temp[[i]]$alpha_p
          beta_Total <- beta_Total + temp[[i]]$beta_p

          fp_bp_Sampling[iter, i,] <- temp[[i]]$fp_bp
          accepted_fp_bp_Sampling[iter, i] <- temp[[i]]$accepted_fp_bp
          fp_bp_scaling_Sampling[iter, i] <- temp[[i]]$scaling
          data_List[[i]]$fp_bp_accept_rate[iter] <- temp[[i]]$fp_bp_accept_rate

          length_gene <- length(temp[[i]]$gene)

          if(length_gene > 0){
            data_List[[i]]$zpk_nu_Sampling[iter, ] <- temp[[i]]$zp[1:length_gene]
            for (j in 1:length_gene) {
              gene_List[[temp[[i]]$gene[j]]]$data[iter, patients[i]] <- temp[[i]]$rp[j]
              gene_List[[temp[[i]]$gene[j]]]$accepted_rpk_Sampling[iter, patients[i]] <- temp[[i]]$accepted_rp[j]
            }
          }

          if(length_gene < length(temp[[i]]$zp)){
            data_List[[i]]$zpk_nnu_Sampling[iter,] <- temp[[i]]$zp[(length_gene + 1):length(temp[[i]]$zp)]
          }
        }

        ### update pi
        pi_Sampling[iter] <- rbeta(1, pi_alpha + alpha_Total, pi_beta + beta_Total)

      }


      return(list(fitness_Inference_Input_Gene = fitness_Inference_Input_Gene,
                  sample_Patient = sample_Patient,
                  gene_Temp = gene_Temp,
                  driver_gene_Temp = driver_gene_Temp,
                  fp_bp_Sampling = fp_bp_Sampling,
                  accepted_fp_bp_Sampling = accepted_fp_bp_Sampling,
                  fp_bp_scaling_Sampling = fp_bp_scaling_Sampling,
                  data_List = data_List,
                  gene_List = gene_List,
                  pi_alpha = pi_alpha,
                  pi_beta = pi_beta,
                  pi_Sampling = pi_Sampling,
                  scaleNum = scaleNum,
                  driver_gene_Temp_t_min = driver_gene_Temp_t_min
      ))

    }

    # You can also use this function by este::fitness_Inference_Patient_Gene

  }

  mclapply(1:dim(temp)[1], function(index){
    cancer_temp <- temp[index,1]
    cat(cancer_index,"  ",cancer_temp,"\n")

    res_cancer_file <- paste(res_file,cancer_temp,"/Immune/",sep = "")
    if(!dir.exists(res_cancer_file)){
      dir.create(res_cancer_file)
    }

    res_cancer_file1 <- paste(res_cancer_file,"mcmcRes/",sep = "")
    if(!dir.exists(res_cancer_file1)){
      dir.create(res_cancer_file1)
    }

    type <- temp[index,2]
    seed_this <- as.integer(temp[index,3])

    if(file.exists(paste(res_cancer_file,type,"_prior_data.rds",sep = "") )){
      res_cancer_file2 <- paste(res_cancer_file1,type,"/",sep = "")
      if(!dir.exists(res_cancer_file2)){
        dir.create(res_cancer_file2)
      }

      prior_data <- readRDS(paste(res_cancer_file,type,"_prior_data.rds",sep = ""))
      set.seed(seed_this)
      res <- fitness_Inference_Patient_Gene(prior_data = prior_data,
                                            sigma_p_square = 0.0001,
                                            max_iter = 10000,
                                            cellular_clock = "censored_Rate",
                                            target_recept_rate = 0.234,
                                            f_b_scaling_update_strength = 10,
                                            r_update_sigma = 0.01,
                                            show_size = 1000,
                                            seed = seed_this)

      saveRDS(res, paste(res_cancer_file2,"res_",seed_this,".rds",sep = ""))
    }

    return(NULL)
  }, mc.cores = 10, mc.cleanup = TRUE)


}





