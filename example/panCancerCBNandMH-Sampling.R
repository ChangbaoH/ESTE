## pan Cancer CBN and MH-Sampling

# Note on Example Code
# The example code provided here has been reorganized for clarity and testing purposes.
# Please review it carefully and modify as necessary for your specific use case.
stop()

################# Genotypes CBN Infer este
{
  library("este")
  library("parallel")
  ICGC_drivers <- readRDS("./data/ICGC/ICGC_drivers1.rds")
  res_file <- "./data/Results/"

  for (cancer_index in 1:length(ICGC_drivers)) {
    cancer_temp <- names(ICGC_drivers)[cancer_index]
    res_cancer_file <- paste(res_file,cancer_temp,"/",sep = "")
    res_cancer_file <- paste(res_cancer_file, "este/",sep = "")

    if(!dir.exists(res_cancer_file)){
      dir.create(res_cancer_file)
    }

    Genotype <- read.csv(paste(res_cancer_file, "Genotype.csv", sep = ""), check.names = F)
    Genotype <- Genotype[order(Genotype[,"dataSet"]),]
    write.csv(Genotype, paste(res_cancer_file, "Genotype_sort.csv", sep = ""), row.names = F)
    Genotype <- read.csv(paste(res_cancer_file, "Genotype_sort.csv", sep = ""), check.names = F)

    cat("-----------------------------------\n")
    cat(cancer_temp, ": ", dim(Genotype)[1],"\n")
    dataSet_info <- readRDS(paste(res_cancer_file, "dataSet_info.rds", sep = ""))

    for (i in 1:length(unique(Genotype[,"dataSet"]))) {
      temp <- unique(Genotype[,"dataSet"])[i]
      for (j in 1:length(dataSet_info)) {
        if(dataSet_info[[j]]$dataset_index == temp){
          cat(dataSet_info[[j]]$dataSet, ": ",length(which(Genotype[,"dataSet"] == temp)), "\n")
        }
      }
    }

    freq_table <- table(Genotype[,"dataSet"])
    eps_res <- as.data.frame(array(-1, dim = c(length(unique(Genotype[,"dataSet"])),3)))
    eps_res[,3] <- unique(Genotype[,"dataSet"])
    Genotype <- Genotype[which(Genotype[,"dataSet"] %in% names(freq_table[freq_table > 90]) ),]

    driver_chr <- ICGC_drivers[[cancer_index]]$driver_chr
    driver_genes <- ICGC_drivers[[cancer_index]]$driver_gene

    dataSet_unique <- unique(Genotype[,"dataSet"])

    mat <- as.matrix(Genotype[, c(driver_genes, driver_chr)])
    isF <- matrix(as.integer(c(1)) ,nrow = length(dataSet_unique), ncol = 2)
    setD <- as.data.frame(array(NA, dim = c(length(dataSet_unique),2)))
    for (l in 1:length(dataSet_unique)) {
      setD[l,1] <- min(which(Genotype[,"dataSet"] == dataSet_unique[l]) - 1)
      setD[l,2] <- max(which(Genotype[,"dataSet"] == dataSet_unique[l]) - 1)
    }
    for (l in 1:dim(setD)[2]) {
      setD[,l] <- as.integer(setD[,l])
    }
    setD <- as.matrix(setD)

    eventD <- as.data.frame(array(NA, dim = c(2,2)))
    eventD[1,1] <- 1
    eventD[1,2] <- length(driver_genes)
    eventD[2,1] <- length(driver_genes) + 1
    eventD[2,2] <- length(driver_genes) + length(driver_chr)
    for (l in 1:dim(eventD)[2]) {
      eventD[,l] <- as.integer(eventD[,l])
    }
    eventD <- as.matrix(eventD)

    Geno <- cbind(rep(as.integer(1),dim(mat)[1]),mat)
    Geno <- as.data.frame(Geno)
    for (l in 1:dim(Geno)[2]) {
      Geno[,l] <- as.integer(Geno[,l])
    }
    Geno <- as.matrix(Geno)


    poset_np <- min(length(c(driver_genes, driver_chr)), 8)

    eps <- estimate_Epsilon_ForMulti(pat = Geno, isF = isF, setD = setD, eventD = eventD, multi_thrds = 1,
                                     threshold1 = as.integer(8), n_p = as.integer(8),
                                     T = 10.0, N_iter = 200L, thrds=as.integer( 10 ))
    saveRDS(eps, paste(res_cancer_file, "este_eps.rds", sep = ""))


    for (i in 1:dim(eps)[1]) {
      eps_res[which(eps_res[,3] == dataSet_unique[i]),1:2] <- eps[i,1:2]
    }

    poset <- find_Poset_ForVote(pat = Geno, isF = isF, eps = eps, setD = setD, eventD = eventD,
                                Fine_Tune_Num = as.integer(2L), vote_size = 5L, vote_threshold = 0.5,
                                vote_thrds = 1L, threshold = 0.0001, threshold2 = 0.01, n_p = as.integer(poset_np), is_update_eps = FALSE,
                                T = 10, N_iter = 200L, thrds = as.integer(20))

    fit <- estimate_Lambda(pat = Geno, poset = poset, isF = isF, eps = eps, setD = setD, eventD = eventD,
                           lambdaS = 1.0, L = 100L, sampling = 'add-remove',
                           maxIter = 500L, updateStepSize=20L, tol=0.001, maxLambda=1e6,
                           neighborhoodDist=1L, is_update_eps = FALSE,
                           thrds=as.integer(20))
    fit$lambdaS <- 1.0
    fit$poset <- poset

    for (i in 1:dim(eps_res)[1]) {
      if(eps_res[i,1] < 0){
        Genotype <- read.csv(paste(res_cancer_file, "Genotype_sort.csv", sep = ""), check.names = F)
        Genotype <- Genotype[which(Genotype[,"dataSet"] %in% eps_res[i,3] ),]

        driver_chr <- ICGC_drivers[[cancer_index]]$driver_chr
        driver_genes <- ICGC_drivers[[cancer_index]]$driver_gene

        dataSet_unique <- unique(Genotype[,"dataSet"])

        mat <- as.matrix(Genotype[, c(driver_genes, driver_chr)])
        isF <- matrix(as.integer(c(1)) ,nrow = length(dataSet_unique), ncol = 2)
        setD <- as.data.frame(array(NA, dim = c(length(dataSet_unique),2)))
        for (l in 1:length(dataSet_unique)) {
          setD[l,1] <- min(which(Genotype[,"dataSet"] == dataSet_unique[l]) - 1)
          setD[l,2] <- max(which(Genotype[,"dataSet"] == dataSet_unique[l]) - 1)
        }
        for (l in 1:dim(setD)[2]) {
          setD[,l] <- as.integer(setD[,l])
        }
        setD <- as.matrix(setD)

        eventD <- as.data.frame(array(NA, dim = c(2,2)))
        eventD[1,1] <- 1
        eventD[1,2] <- length(driver_genes)
        eventD[2,1] <- length(driver_genes) + 1
        eventD[2,2] <- length(driver_genes) + length(driver_chr)
        for (l in 1:dim(eventD)[2]) {
          eventD[,l] <- as.integer(eventD[,l])
        }
        eventD <- as.matrix(eventD)

        Geno <- cbind(rep(as.integer(1),dim(mat)[1]),mat)
        Geno <- as.data.frame(Geno)
        for (l in 1:dim(Geno)[2]) {
          Geno[,l] <- as.integer(Geno[,l])
        }
        Geno <- as.matrix(Geno)

        temp <- estimate_Epsilon_based_on_Poset_and_Lambda(pat = Geno, poset = fit$poset, isF = isF, eps = matrix(colMeans(fit$epsilon1), nrow = 1),
                                                           setD = setD, eventD = eventD, lambda = fit$lambda,
                                                           lambdaS = 1.0, L = 100L, sampling = 'add-remove',
                                                           maxIter = 500L, updateStepSize=20L, tol=0.001, maxLambda=1e6,
                                                           neighborhoodDist=1L, thrds=as.integer(5))

        eps_res[i,1:2] <- unique(as.numeric(temp$epsilon2))
      }
    }

    fit$eps_res <- eps_res

    saveRDS(fit, paste(res_cancer_file, "este_fit.rds", sep = ""))

  }


}



################# Genotypes CBN Infer mccbn
{
  library("este")
  library("parallel")
  library("mccbn")

  ICGC_drivers <- readRDS("./data/ICGC/ICGC_drivers1.rds")
  res_file <- "./data/Results/"

  for (cancer_index in 1:length(ICGC_drivers)) {
    cat(cancer_index,"\n")
    cancer_temp <- names(ICGC_drivers)[cancer_index]
    res_cancer_file <- paste(res_file,cancer_temp,"/",sep = "")
    res_cancer_file <- paste(res_cancer_file, "este/",sep = "")

    Genotype <- read.csv(paste(res_cancer_file, "Genotype_sort.csv", sep = ""), check.names = F)
    driver_chr <- ICGC_drivers[[cancer_index]]$driver_chr
    driver_genes <- ICGC_drivers[[cancer_index]]$driver_gene
    h_cbn_data <- as.matrix(Genotype[, c(driver_genes, driver_chr)])

    posets = candidate_posets(h_cbn_data, rep(1, dim(h_cbn_data)[1]), 0.9)
    poset0 = posets[[length(posets)]]
    fit = adaptive.simulated.annealing(poset0, h_cbn_data, L=100, max.iter.asa=10, seed=10L, thrds=10L)
    saveRDS(fit, paste(res_cancer_file, "mccbn_fit.rds", sep = ""))
  }



}



################# Genotypes MH Sampling
{
  library("MASS")
  library("fitdistrplus")
  library("matrixcalc")
  library("quadprog")
  library("Matrix")
  library("mclust")
  library("data.table")
  library("readr")
  library("devtools")
  load_all()

  ICGC_Cli_patient <- as.data.frame(fread("./data/ICGC/pancan_pcawg_2020/data_clinical_patient.txt", skip = 4, fill = TRUE, sep = "\t"))

  ICGC_drivers <- readRDS("./data/ICGC/ICGC_drivers1.rds")
  res_file <- "./data/Results/"

  sample <- read.csv("./data/ICGC/pcawg_sample_sheet.csv")
  time <- read.csv("./data/ICGC/time.csv")
  cli <- read.csv("./data/ICGC/cli.csv",row.names = 1)

  chr1 <- read.csv("./data/Genom/chromData.csv", row.names = 1)
  chr <- read.csv("./data/Genom/GRCH37_chrBand1.csv")

  ## cBioPortal Cli
  {
    age_string <- unique(c("Diagnosis Age","Age at Which Sequencing was Reported (Years)",
                           "AGE","Age at First Diagnosis", "Age at Sequencing",
                           "Age At Diagnosis", "Age at Sampling", "Age at Diagnosis",
                           "Age", "Age At Surgery","Age at Surgery/Biopsy","Age (yrs)",
                           "Age at Resection","Age (Sampling)","Patient Age At Diagnosis",
                           "Age at Sample Collection","Age At Procurement","Age (tumor resected)",
                           "Age at Initial Diagnosis","Age At Initial Pathological Diagnosis","Age at SCC"))

    for (cancer_index in 1:length(ICGC_drivers)) {

      cancer_temp <- names(ICGC_drivers)[cancer_index]
      res_cancer_file <- paste(res_file,cancer_temp,"/time/",sep = "")
      if(!dir.exists(res_cancer_file)){
        dir.create(res_cancer_file)
      }

      cBioCli <- as.data.frame(read.table(paste("./data/cBioPortal/cli/", cancer_temp, "/combined_study_clinical_data.tsv", sep = ""),
                                          comment.char = "#", quote = "", header = T, fill = TRUE, check.names = FALSE, sep = "\t"))

      age_string_temp <- intersect(age_string, colnames(cBioCli))

      if(cancer_temp == "CNS-GBM"){
        for (i in 1:dim(cBioCli)[1]) {
          temp <- na.omit(as.numeric(cBioCli[i,"Age at Initial Diagnosis"])/365)
          if(length(temp) > 0){
            for (j in 1:length(age_string_temp)) {
              if( !is.na(as.numeric(cBioCli[i, age_string_temp[j]])) & max(abs(as.numeric(cBioCli[i, age_string_temp[j]]) - temp), 1) == 1 ){
                cBioCli[i, age_string_temp[j]] <- temp
              }
            }
          }
        }

      }

      cBioCli$"AGE_cache" <- rep(NA, dim(cBioCli)[1])
      for (i in 1:dim(cBioCli)[1]) {
        temp <- c()
        for (j in 1:length(age_string_temp)) {
          temp1 <- na.omit(as.numeric(cBioCli[i,age_string_temp[j]]))
          if(length(temp1) > 0){
            temp <- c(temp, mean(as.numeric(temp1)))
          }
        }

        if(length(temp) > 0){
          cBioCli[i,"AGE_cache"] <- min(temp)
        }
      }
      cBioCli <- cBioCli[,c("Patient ID","AGE_cache")]
      colnames(cBioCli)[2] <- "AGE"
      saveRDS(cBioCli, paste(res_cancer_file, "cBioCli_filter.rds",sep = ""))
    }

  }

  for (cancer_index in 1:length(ICGC_drivers)) {
    cancer_temp <- names(ICGC_drivers)[cancer_index]
    res_cancer_file <- paste(res_file,cancer_temp,"/este/",sep = "")
    Genotype <- read.csv(paste(res_cancer_file, "Genotype_sort.csv", sep = ""), check.names = F)
    cat("-----------------------------------\n")
    cat(cancer_temp, ": ", dim(Genotype)[1],"\n")
    dataSet_info <- readRDS(paste(res_cancer_file, "dataSet_info.rds", sep = ""))

    for (i in 1:length(dataSet_info)) {
      if(all(dataSet_info[[i]]$dataSet == "ICGC") ){
        data_Set_ICGC <- dataSet_info[[i]]$dataset_index
      }
    }

    for (i in 1:length(unique(Genotype[,"dataSet"]))) {
      temp <- unique(Genotype[,"dataSet"])[i]
      cat(i, ": ", temp,"  ", min(which(Genotype[,"dataSet"] == temp))," ",max(which(Genotype[,"dataSet"] == temp)),"\n")
    }

    driver_chr <- ICGC_drivers[[cancer_index]]$driver_chr
    driver_genes <- ICGC_drivers[[cancer_index]]$driver_gene

    chr3 <- as.data.frame(array(NA, dim = c(length(driver_chr),10)))
    colnames(chr3) <- c(colnames(chr1),"evnet","evnetKind")
    chr3[,9] <- driver_chr
    for (i in 1:dim(chr3)[1]) {
      temp <- substr(driver_chr[i], 2, nchar(driver_chr[i]))
      chr3[i, 1:8] <- chr1[which(chr1[,8] == temp),1:8]
      chr3[i, 10] <- substr(driver_chr[i], 1, 1)
    }

    chr3 <- chr3[order(chr3[,3],decreasing = F),]
    chr3 <- chr3[order(chr3[,2],decreasing = F),]
    chr3 <- chr3[order(as.integer(chr3[,5]),decreasing = F),]

    write.csv(chr3, "./data/Genom/chromData3.csv",row.names = F)
    chr3 <- read.csv("./data/Genom/chromData3.csv")

    ##### time_Donor_Geno
    time_Donor_Geno <- Genotype[which(Genotype$dataSet == data_Set_ICGC),]
    rownames(time_Donor_Geno) <- time_Donor_Geno[,"sample"]
    time_Donor_Geno <- as.matrix(time_Donor_Geno[,c(driver_genes, driver_chr)])

    res <- readRDS(paste(res_cancer_file, "este_fit.rds", sep = ""))
    event <- c(driver_genes,driver_chr)

    rateT <- rateTimeHelp(res$lambdaS, res$lambda, res$poset)
    for (i in 1:length(event)) {
      cat(event[i],"   ",rateT$LB[i],"   ",rateT$UB[i],"\n" ,sep = "" )
    }

    gammaClusterRes <- GammaCluster(rateT$LB, c(0.9,0.1))

    res_cancer_file <- paste(res_file,cancer_temp,"/time/",sep = "")
    if(!dir.exists(res_cancer_file)){
      dir.create(res_cancer_file)
    }

    ##### time_data
    time_temp <- time[which(time[,13] == cancer_temp & !(time[,10] %in% c("WGD","DoubleGain"))),]
    time_temp <- cbind(time_temp,array(NA,dim = c(dim(time_temp)[1],5)))
    colnames(time_temp)[14:18] <- c("chrN", "event", "num", "donor", "age")

    for (i in 1:dim(time_temp)[1]) {
      time_temp[i,17] <- sample[which(sample[,6] == time_temp[i,1]),4]
      time_temp[i,18] <- cli[time_temp[i,1],1]
      if(is.na(time_temp[i,18])){
        time_temp[i,18] <- ICGC_Cli_patient[which(ICGC_Cli_patient[,1] == time_temp[i,17]),"AGE"]
      }
      time_temp[i,14] <- paste("chr",time_temp[i,2],sep = "")
      if(time_temp[i,10]=="SingleGain"){
        np <- "+"
      }else{
        np <- "-"
      }
      temp <- chr3[which(chr3[, 1] == time_temp[i,14] & (chr3[, 2] <= time_temp[i,4] & chr3[, 3] >= time_temp[i,3]) & chr3[, 10] == np),9]

      time_temp[i,15] <- temp[1]
      time_temp[i,16] <- length(temp)
      if(length(temp) > 1){
        for (j in 2:length(temp)) {
          time_temp[i,15] <- paste(time_temp[i,15], "|", temp[j], sep = "")
        }
      }
    }

    time_temp1 <- time_temp[which(time_temp[,16] > 0),]


    # time_GMM <- GaussCluster(time_temp1$time, 3, c(0.2,0.6,0.2), c(0,0.5,1), c(1e-4,0.1,1e-4))
    # snv_GMM <- GaussCluster(time_temp1$no.snvs,3,c(0.2,0.6,0.2),c(2,17,100),c(1,100,10000))
    # saveRDS(time_GMM, paste(res_cancer_file, "time_GMM.rds", sep = ""))
    # saveRDS(snv_GMM, paste(res_cancer_file, "snv_GMM.rds", sep = ""))
    # time_GMM <- readRDS(paste(res_cancer_file, "time_GMM.rds", sep = ""))
    # snv_GMM <- readRDS(paste(res_cancer_file, "snv_GMM.rds", sep = ""))

    time_GMM <- Mclust(time_temp1$time, G=3, verbose=FALSE)
    snv_GMM <- Mclust(time_temp1$no.snvs, G=3, verbose=FALSE)
    saveRDS(time_GMM, paste(res_cancer_file, "time_GMM.rds", sep = ""))
    saveRDS(snv_GMM, paste(res_cancer_file, "snv_GMM.rds", sep = ""))
    time_GMM <- readRDS(paste(res_cancer_file, "time_GMM.rds", sep = ""))
    snv_GMM <- readRDS(paste(res_cancer_file, "snv_GMM.rds", sep = ""))

    noSelectTimeIndex1 <- which(time_GMM$classification != 2)
    noSelectTimeIndex2 <- which(snv_GMM$classification == 1)

    filter <- time_temp1[-intersect(noSelectTimeIndex1, noSelectTimeIndex2),]
    filter_donor <- unique(filter$donor)

    filter_num <- array(NA,dim = c(length(filter_donor),length(driver_chr)))
    rownames(filter_num) <- filter_donor
    colnames(filter_num) <- driver_chr
    filter_time <- array(NA,dim = c(length(filter_donor),length(driver_chr)+1))
    rownames(filter_time) <- filter_donor
    colnames(filter_time) <- c(driver_chr, "age")
    for (i in 1:dim(filter_time)[1]) {
      filter_time[i, "age"] <- filter[which(filter$donor == rownames(filter_time)[i])[1], 18]
    }
    filter_time <- filter_time[which(!is.na(filter_time[,dim(filter_time)[2]])),]
    filter_donor <- rownames(filter_time)
    filter_snvN <- filter_time[,-(length(driver_chr)+1)]
    filter_num <- filter_num[rownames(filter_time),]

    for (i in 1:dim(filter_time)[1]) {
      temp <- filter[which(filter$donor == rownames(filter_time)[i]),]
      temp1 <- array(0, dim = c(length(driver_chr),3))
      rownames(temp1) <- driver_chr
      for (j in 1:dim(temp)[1]) {
        eventTemp <- strsplit(x=temp[j, 15],split='[|]')[[1]]
        temp1[eventTemp, 1] <- temp1[eventTemp, 1] + temp[j, 6]
        temp1[eventTemp, 2] <- temp1[eventTemp, 2] + 1
        temp1[eventTemp, 3] <- temp1[eventTemp, 1] + temp[j, 9]
      }
      temp1 <- temp1[which(temp1[,2] > 0),,drop = FALSE]
      if(dim(temp1)[1] > 0){
        for (j in 1:dim(temp1)[1]) {
          filter_time[i, rownames(temp1)[j]] <- temp1[j,1]/temp1[j,2]
          filter_num[i, rownames(temp1)[j]] <- temp1[j,3]/temp1[j,2]
        }
      }
    }

    cancer_time <- array(NA, dim = c(dim(filter_time)[1], dim(filter_time)[2]+2))
    rownames(cancer_time) <- rownames(filter_time)
    colnames(cancer_time) <- c(colnames(filter_time), "healthTime", "sampleTime")
    cancer_time[,1:(dim(cancer_time)[2]-2)] <- filter_time[,1:(dim(cancer_time)[2]-2)]
    filter_snvN <- t(rbind(filter_snvN, rep(2,dim(filter_snvN)[2])))

    saveRDS(rateT, paste(res_cancer_file, "rateT.rds", sep = ""))
    saveRDS(gammaClusterRes, paste(res_cancer_file, "gammaClusterRes.rds", sep = ""))
    saveRDS(res, paste(res_cancer_file, "este_fit.rds", sep = ""))
    saveRDS(time_temp, paste(res_cancer_file, "time_temp.rds", sep = ""))
    saveRDS(time_temp1, paste(res_cancer_file, "time_temp1.rds", sep = ""))
    saveRDS(filter, paste(res_cancer_file, "filter.rds", sep = ""))
    saveRDS(filter_time, paste(res_cancer_file, "filter_time.rds", sep = ""))
    saveRDS(filter_num, paste(res_cancer_file, "filter_num.rds", sep = ""))
    saveRDS(cancer_time, paste(res_cancer_file, "cancer_time.rds", sep = ""))
    saveRDS(filter_snvN, paste(res_cancer_file, "filter_snvN.rds", sep = ""))




    ##### create most compatible genotyp for true genotyp
    cBioCli <- readRDS(paste(res_cancer_file, "cBioCli_filter.rds",sep = ""))
    cBioCli <- cBioCli[which(!is.na(cBioCli[,"AGE"])),]
    cBioCli[,"AGE"] <- round(cBioCli[,"AGE"], 2)
    cBioCli <- cBioCli[which(cBioCli[,"AGE"] > 0),]

    genotype_other <- Genotype[which(Genotype$dataSet != data_Set_ICGC),]
    dornor_other <- intersect(cBioCli[,"Patient ID"], genotype_other[,"sample"])
    cBioCli <- cBioCli[which(cBioCli[,"Patient ID"] %in% dornor_other),]
    genotype_other <- genotype_other[which(genotype_other[,"sample"] %in% dornor_other),]

    used_Genotype_For_Sampling <- Genotype[which(Genotype$dataSet == data_Set_ICGC),]
    rownames(used_Genotype_For_Sampling) <- used_Genotype_For_Sampling[,"sample"]
    used_Genotype_For_Sampling <- as.matrix(used_Genotype_For_Sampling[, c(driver_genes, driver_chr)])

    temp <- genotype_other
    rownames(temp) <- temp[, "sample"]
    temp <- as.matrix(temp[, c(driver_genes, driver_chr)])
    used_Genotype_For_Sampling <- rbind(used_Genotype_For_Sampling, temp)

    Donor_Pair_Genotype <- as.data.frame(array(0, dim = dim(used_Genotype_For_Sampling)))
    rownames(Donor_Pair_Genotype) <- rownames(used_Genotype_For_Sampling)
    colnames(Donor_Pair_Genotype) <- colnames(used_Genotype_For_Sampling)


    for (i in 1:dim(res$eps_res)[1]) {
      donor_temp <- Genotype[which(Genotype[,"dataSet"] ==  res$eps_res[i,3]),"sample"]
      donor_temp <- intersect(donor_temp, rownames(used_Genotype_For_Sampling))
      if( length(donor_temp) > 0 ){
        used_Genotype_For_Sampling_temp <- used_Genotype_For_Sampling[donor_temp,,drop = FALSE]
        eps_temp <- c(rep(res$eps_res[i,1], length(driver_genes)), rep(res$eps_res[i,2], length(driver_chr)) )

        Donor_Pair_Geno_temp <- find_most_Compatible_Genotype_by_Flipping(used_Genotype_For_Sampling_temp,
                                                                          eps_temp, res$poset, max_iter = 100000, cores = 10)
        for (j in 1:dim(Donor_Pair_Geno_temp)[1]) {
          Donor_Pair_Genotype[donor_temp[j],] <- Donor_Pair_Geno_temp[donor_temp[j],]
        }

      }
    }


    saveRDS(Donor_Pair_Genotype, paste(res_cancer_file, "Donor_Pair_Genotype.rds", sep = ""))

    ###################### filter
    Donor_Pair_Genotype <- readRDS(paste(res_cancer_file, "Donor_Pair_Genotype.rds", sep = ""))

    filter_Donor_Pair_Genotype <- donor_Pair_Genotype_Filter(Donor_Pair_Genotype, cancer_time)

    saveRDS(filter_Donor_Pair_Genotype, paste(res_cancer_file, "filter_Donor_Pair_Genotype.rds", sep = ""))

    filter_Donor_Pair_Genotype <- readRDS(paste(res_cancer_file, "filter_Donor_Pair_Genotype.rds", sep = ""))
    topo_path <- topological_Sort(res$poset)
    parent_Set <- parent_Set_Help(res$poset)

    cBioCli <- readRDS(paste(res_cancer_file, "cBioCli_filter.rds",sep = ""))
    cBioCli <- cBioCli[which(!is.na(cBioCli[,"AGE"])),]
    cBioCli[,"AGE"] <- round(cBioCli[,"AGE"], 2)
    cBioCli <- cBioCli[which(cBioCli[,"AGE"] > 0),]
    cBioCli <- cBioCli[which(cBioCli[,"Patient ID"] %in% intersect(cBioCli[,"Patient ID"], Genotype[,"sample"])), ]

    Donor_Pair_Geno_Last <- MH_Pretreatment_based_on_flipping(filter_Donor_Pair_Genotype, time_temp, cancer_time, filter_num, rateT,
                                                              driver_genes, TRUE, cBioCli, topo_path, parent_Set, gammaClusterRes)


    saveRDS(topo_path, paste(res_cancer_file, "topo_path.rds", sep = ""))
    saveRDS(parent_Set, paste(res_cancer_file, "parent_Set.rds", sep = ""))
    saveRDS(Donor_Pair_Geno_Last, paste(res_cancer_file, "Donor_Pair_Geno_Last.rds", sep = ""))

  }

  ## MCMC
  library("parallel")
  library("fitdistrplus")

  res_file <- "./data/Results/"
  ICGC_drivers <- readRDS("./data/ICGC/ICGC_drivers1.rds")

  for (cancer_index in 1:length(ICGC_drivers)) {
    cancer_temp <- names(ICGC_drivers)[cancer_index]
    res_cancer_file <- paste(res_file,cancer_temp,"/time/",sep = "")
    res_cancer_file1 <- paste(res_file,cancer_temp,"/time1/",sep = "")

    if(!dir.exists(res_cancer_file1)){
      dir.create(res_cancer_file1)
    }

    Donor_Pair_Geno_Last <- readRDS(paste(res_cancer_file, "Donor_Pair_Geno_Last.rds", sep = ""))

    temp <- MH_Sampling_based_on_flipping(Donor_Pair_Geno_Last, saveFile = paste(res_cancer_file1,"MH_Res_Last.rds",sep = ""),
                                  isWeight = TRUE, iterations = 5, g.sf = 10, maxSampleIter = 50,
                                  errorTolerance = 0.01, errorToleranceGrowthRate = 1.5, cores1 = 8,
                                  cores2 = 10, show_size = 1000, seed = 1234)

    saveRDS(temp, paste(res_cancer_file1,"MH_Res_Last.rds",sep = ""))

  }


}




