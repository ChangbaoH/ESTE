## pan cancer data

# Note on Example Code
# The example code provided here has been reorganized for clarity and testing purposes.
# Please review it carefully and modify as necessary for your specific use case.
stop()

library("readr")
library("data.table")
library("rtracklayer")
library("dplyr")
library("parallel")
library("liftOver")
library("org.Hs.eg.db")

library("devtools")
library("parallel")
load_all()

#################  all ICGC cancer
{
  ICGC_time <- read.csv(.Date("./data/ICGC/time.csv"))

  ICGC_Cancer <- as.data.frame(array(NA, dim = c(length(unique(ICGC_time$histology_abbreviation)), 2)))
  colnames(ICGC_Cancer) <- c("cancer","sampleNum")
  ICGC_Cancer$cancer <-  unique(ICGC_time$histology_abbreviation)
  for (i in 1:dim(ICGC_Cancer)[1]) {
    ICGC_Cancer[i,2] <- length(unique(ICGC_time[which(ICGC_time$histology_abbreviation == ICGC_Cancer[i,1]),1]))
  }
  ICGC_Cancer <- ICGC_Cancer[order(ICGC_Cancer[,1]),]
  ICGC_Cancer <- ICGC_Cancer[which(!ICGC_Cancer[,1] %in% c("Bone-Benign","Breast-DCIS","Breast-LobularCA","Cervix-AdenoCA",
                                                           "Myeloid-AML","Myeloid-MDS","Myeloid-MPN","Thy-AdenoCA")),]
  ICGC_Cancer <- ICGC_Cancer[order(ICGC_Cancer[,2],decreasing = T),]

}

#################  cBioPortal
{
  ICGC_Cancer <- ICGC_Cancer[which(ICGC_Cancer[,1] %in% c("Breast-AdenoCA","CNS-GBM","ColoRect-AdenoCA","Liver-HCC","Lung-AdenoCA",
                                                          "Lung-SCC","Prost-AdenoCA","Skin-Melanoma","Uterus-AdenoCA")),]
  ICGC_Cancer <- ICGC_Cancer[order(ICGC_Cancer[,1]),]
  rownames(ICGC_Cancer) <- 1:dim(ICGC_Cancer)[1]
}

#################  drivers ICGC
{
  drivers <- vector("list",dim(ICGC_Cancer)[1])
  names(drivers) <- ICGC_Cancer[,1]
  ICGC_drivers <- read.csv("./data/ICGC/ICGC_driver.csv")
  unique(ICGC_drivers[,1])
  unused_cancer <- c()
  for (i in 1:length(drivers)) {
    temp <- names(drivers)[i]
    if(temp %in% ICGC_drivers[,1]){
      drivers_temp <- ICGC_drivers[which(ICGC_drivers[,1] == temp),2]
      list_temp <- list(driver_chr = c(), driver_gene = c())
      for (j in 1:length(drivers_temp)) {
        if(substr(drivers_temp[j],1,5) %in% c("gain_","loss_","homde")){
          drivers_temp1 <- gsub("gain_","+",drivers_temp[j])
          drivers_temp1 <- gsub("loss_","-",drivers_temp1)
          drivers_temp1 <- gsub("homdel_","-",drivers_temp1)
          list_temp$driver_chr <- c(list_temp$driver_chr, drivers_temp1)
        }else if(drivers_temp[j] != "WGD"){
          list_temp$driver_gene <- c(list_temp$driver_gene, drivers_temp[j])
        }
      }
      drivers[[i]] <- list_temp
    }else{
      unused_cancer <- c(unused_cancer, temp)
    }
  }
  ICGC_Cancer <- ICGC_Cancer[which(!ICGC_Cancer[,1] %in% unused_cancer),]
  drivers_temp <- list()
  for (i in 1:length(drivers)) {
    if(! names(drivers)[i] %in% unused_cancer){
      drivers_temp[[length(drivers_temp) + 1]] <- drivers[[i]]
      names(drivers_temp)[length(drivers_temp)] <- names(drivers)[i]
    }
  }
  write.csv(ICGC_Cancer, "./data/ICGC/ICGC_Cancer.csv",row.names = F)
  saveRDS(drivers_temp, "./data/ICGC/ICGC_drivers.rds")
  ICGC_Cancer <- read.csv("./data/ICGC/ICGC_Cancer.csv")
  ICGC_drivers <- readRDS("./data/ICGC/ICGC_drivers.rds")

}

#################  driver genes
{
  driver_gene_dataset <- as.data.frame(array(NA, dim = c(length(ICGC_drivers), 2)))
  colnames(driver_gene_dataset) <- c("ICGC","compendium")
  driver_gene_dataset[,1] <- ICGC_Cancer[,1]
  driver_gene_dataset[,2] <- c("BRCA","GB|GBM","COAD|COADREAD|READ","HCC","LUAD","LUSC","PRAD","MEL|SKCM","UCEC")

  # A compendium of mutational cancer driver genes | IntOGen-Drivers-20240920
  compendium_drivers <- as.data.frame(read_tsv("./data/driverGenes/Compendium_Cancer_Genes.tsv"))
  unique(compendium_drivers$CANCER_TYPE)
  cohort <- read_tsv("./data/driverGenes/cohorts.tsv")
  cohort <- unique(cohort[,2:3])
  write.csv(cohort, "./data/driverGenes/cohort_cancer.csv",row.names = F)
  cohort <- read.csv("./data/driverGenes/cohort_cancer.csv")
  for (i in 1:length(ICGC_drivers)) {
    temp <- names(ICGC_drivers)[i]
    temp1 <- driver_gene_dataset[which(driver_gene_dataset[,1] == temp), 2]
    temp1 <- strsplit(temp1, "[|]")[[1]]
    temp1 <- compendium_drivers[which(compendium_drivers[,4] %in% temp1 & compendium_drivers[,12] == "TRUE"),]
    temp1 <- unique(temp1[,1])
    ICGC_drivers[[i]]$driver_gene <- unique(c(ICGC_drivers[[i]]$driver_gene, temp1))
  }
  saveRDS(ICGC_drivers, "./data/ICGC/ICGC_drivers.rds")
  ICGC_drivers <- readRDS("./data/ICGC/ICGC_drivers.rds")
}


################# ICGC infomation dele
# donor <- read.csv("./data/ICGC/donor.all_projects.csv")
sample <- read.csv("./data/ICGC/pcawg_sample_sheet.csv")
time <- read.csv("./data/ICGC/time.csv")
cli <- read.csv("./data/ICGC/cli.csv",row.names = 1)
# unique(time$histology_abbreviation)
rm(donor)
################# ICGC time
chr <- read.csv("./data/Genom/chromData.csv",row.names = 1)
mclapply(names(ICGC_drivers), function(cancer_temp){
  file1 <- paste("./data/ICGC/",cancer_temp,"/",sep = "")
  if(!dir.exists(file1)){
    dir.create(file1)
  }

  cancer_chr_data <- time[which(time[,13] == cancer_temp & !(time[,10] %in% c("WGD","DoubleGain"))),]
  cancer_chr_data <- cbind(cancer_chr_data,array(NA,dim = c(dim(cancer_chr_data)[1],10)))
  colnames(cancer_chr_data)[14:23] <- c("chrN1","chrN2","event1","event2","chrS1","chrE1","chrS2","chrE2","donor","age")

  cancer_chr_data$chrS1 <- as.integer(cancer_chr_data$chrS1)
  cancer_chr_data$chrE1 <- as.integer(cancer_chr_data$chrE1)
  cancer_chr_data$chrS2 <- as.integer(cancer_chr_data$chrS2)
  cancer_chr_data$chrE2 <- as.integer(cancer_chr_data$chrE2)
  for (i in 1:dim(cancer_chr_data)[1]) {
    cancer_chr_data[i,22] <- sample[which(sample[,6] == cancer_chr_data[i,1]),4]
    cancer_chr_data[i,23] <- cli[cancer_chr_data[i,1],1]
    cancer_chr_data[i,14] <- paste("chr",cancer_chr_data[i,2],sep = "")
    temp <- chr[which(chr[,1] == cancer_chr_data[i,14]),]
    l <- temp[which(temp[,2] <= cancer_chr_data[i,3] & temp[,3] >= cancer_chr_data[i,4]),,drop = F]
    if (dim(l)[1] > 0){
      cancer_chr_data[i,15] <- paste(temp[1,6],l[1,4],sep = "")
      cancer_chr_data[i,20:21] <- l[1,2:3]
    }
    mp <- max(temp[which(temp[,7] == "p"),3])
    mq <- max(temp[which(temp[,7] == "q"),3])
    if (cancer_chr_data[i,4] <= mp){
      cancer_chr_data[i,14] <- paste(temp[1,6],"p",sep = "")
      cancer_chr_data[i,18:19] <- c(as.integer(1),as.integer(mp))
    }else{
      cancer_chr_data[i,14] <- paste(temp[1,6],"q",sep = "")
      cancer_chr_data[i,18:19] <- c(as.integer(mp+1),as.integer(mq))
    }
    if(cancer_chr_data[i,10] %in% c("SingleGain")){
      cancer_chr_data[i,16] <- paste("+",cancer_chr_data[i,14],sep = "")
      if(!is.na(cancer_chr_data[i,15])){
        cancer_chr_data[i,17] <- paste("+",cancer_chr_data[i,15],sep = "")
      }
    }else{
      cancer_chr_data[i,16] <- paste("-",cancer_chr_data[i,14],sep = "")
      if(!is.na(cancer_chr_data[i,15])){
        cancer_chr_data[i,17] <- paste("-",cancer_chr_data[i,15],sep = "")
      }
    }
  }
  write.csv(cancer_chr_data, paste(file1, "cancer_chr_data.csv", sep = ""), row.names = F)
}, mc.cores = length(ICGC_drivers), mc.cleanup = TRUE)
rm(cancer_temp, cancer_chr_data)

ICGC_Cli_patient <- as.data.frame(fread("./data/ICGC/pancan_pcawg_2020/data_clinical_patient.txt", skip = 4, fill = TRUE, sep = "\t"))
ICGC_Cli_sample <- as.data.frame(fread("./data/ICGC/pancan_pcawg_2020/data_clinical_sample.txt", skip = 4, fill = TRUE, sep = "\t"))

################# GRCH37 Gene ID
{
  ID_37 <- as.data.frame(rtracklayer::import('./data/Genom/Homo_sapiens.GRCh37.75.gtf.gz'))
  ID_37 <- ID_37[which(ID_37[,7] == "gene"),]
  ID_37 <- ID_37[,c(1,2,3,4,5,10,11)]
  ID_37 <- ID_37[which(ID_37[,1] %in% c(1:22,"X","Y")),]
  ID_37[,1] <- sapply(ID_37[,1], function(x){
    x <- gsub("X","23",x)
    x <- gsub("Y","24",x)
  })
  ID_37 <- distinct(ID_37,gene_id,.keep_all = T)
  write.csv(ID_37,"./data/Genom/ID_37.csv",row.names = F)
  rm(ID_37)
}

ID <- read.csv("./data/Genom/ID_37.csv")
for (i in 1:3) {
  ID[,i] <- as.integer(ID[,i])
}

################# driver chr
min_variance_partition <- function(nums, target) {
  if(sum(nums) <= target){
    return(list(nums))
  }
  find_partitions <- function(nums, target) {
    n <- length(nums)
    all_partitions <- list()
    if (n == 0) return(list(list()))
    if (sum(nums) <= target){
      return(list(list(nums)))
    }
    for (i in 1:n) {
      group_sum <- sum(nums[1:i])
      if (group_sum >= target) {
        if (i == n) {
          all_partitions <- append(all_partitions, list(list(nums[1:i])))
        } else {
          remaining <- nums[(i+1):n]
          sub_partitions <- find_partitions(remaining, target)
          for (partition in sub_partitions) {
            all_partitions <- append(all_partitions, list(c(list(nums[1:i]), partition)))
          }
        }
      }
    }
    return(all_partitions)
  }

  partitions <- find_partitions(nums, target)

  max_group_count <- max(sapply(partitions, length))
  max_partitions <- Filter(function(p) length(p) == max_group_count, partitions)

  variances <- sapply(max_partitions, function(partition) {
    group_sums <- sapply(partition, sum)
    if(length(group_sums) == 1){
      return(0)
    }else{
      var(group_sums)
    }
  })

  best_partition <- max_partitions[[which.min(variances)]]
  return(best_partition)
}
chr_deal <- function(r,r1){
  chr <- read.csv("./data/Genom/GRCH37_chrBand.csv")
  chr <- cbind(chr, array(NA, dim = c(dim(chr)[1], 2)))
  colnames(chr) <- c("chr","chromStart","chromEnd","name","nChr","arm","bandLength")
  chr$nChr <- sapply(chr$chr, function(x){
    x <- gsub("chr","",x)
    x <- gsub("X","23",x)
    x <- gsub("Y","24",x)
  })
  for (i in 1:dim(chr)[1]) {
    if(grepl("p",chr[i,4])){
      chr[i,6] <- "p"
    }else{
      chr[i,6] <- "q"
    }
    chr[i,7] <- as.integer(chr[i,3])  - as.integer(chr[i,2])
  }
  chr$"name2" <- sapply(chr$name, function(x){
    strsplit(x,"[.]")[[1]][1]
  })
  chr$"name2" <- paste(chr[,5], chr$"name2", sep = "")
  write.csv(chr, "./data/Genom/GRCH37_chrBand1.csv", row.names = F)


  chr <- read.csv("./data/Genom/GRCH37_chrBand.csv")
  chr <- cbind(chr, array(NA, dim = c(dim(chr)[1], 3)))
  colnames(chr) <- c("chr","chromStart","chromEnd","name","nChr","arm","bandLength","chrName")
  chr$nChr <- sapply(chr$chr, function(x){
    x <- gsub("chr","",x)
    x <- gsub("X","23",x)
    x <- gsub("Y","24",x)
  })
  chr$name <- sapply(chr$name, function(x){
    strsplit(x,"[.]")[[1]][1]
  })
  chr$name <- paste(chr$nChr, chr$name, sep = "")
  chr_Band <- as.data.frame(array(NA,dim = c(length(unique(chr[,"name"])),8)))
  colnames(chr_Band) <- c(colnames(chr))
  chr_Band[,4] <- unique(chr[,"name"])
  for (i in 1:dim(chr_Band)[1]) {
    chr_Band[i,1] <- unique(chr[which(chr[,4] == chr_Band[i,4]),1])
    chr_Band[i,2] <- min(chr[which(chr[,4] == chr_Band[i,4]),2])
    chr_Band[i,3] <- max(chr[which(chr[,4] == chr_Band[i,4]),3])
    chr_Band[i,5] <- max(chr[which(chr[,4] == chr_Band[i,4]),5])
    if(grepl("p",chr_Band[i,4])){
      chr_Band[i,6] <- "p"
    }else{
      chr_Band[i,6] <- "q"
    }
    chr_Band[i,7] <- as.integer(chr_Band[i,3])  - as.integer(chr_Band[i,2])
  }
  for (i in c(2,3,5,7)) {
    chr_Band[,i] <- as.integer(chr_Band[,i])
  }

  chr <- chr_Band
  mean_L <- mean(as.numeric(chr_Band[,7]))

  for (i in 1:24) {
    temp <- chr_Band[which(chr_Band[,5] == i & chr_Band[,6] == "p"),,drop = FALSE]
    temp <- temp[order(as.numeric(temp[,2]), decreasing = F),]
    temp1 <- min_variance_partition(temp[,7], r*mean_L)
    index <- 1
    for (j in 1:length(temp1)) {
      temp[min(which(is.na(temp[,8]))):(min(which(is.na(temp[,8]))) + length(temp1[[j]]) - 1),8] <-
        rep(paste(temp[1,5], temp[1,6], index, sep = ""), length(temp1[[j]][1]))
      index <- index + 1
    }
    for (j in 1:dim(temp)[1]) {
      chr[which(chr[,4] == temp[j,4]),8] <- temp[j,8]
    }

    temp <- chr_Band[which(chr_Band[,5] == i & chr_Band[,6] == "q"),,drop = FALSE]
    temp <- temp[order(as.numeric(temp[,2]), decreasing = F),]
    temp1 <- min_variance_partition(temp[,7], r*mean_L)
    index <- 1
    for (j in 1:length(temp1)) {
      temp[min(which(is.na(temp[,8]))):(min(which(is.na(temp[,8]))) + length(temp1[[j]]) - 1),8] <-
        rep(paste(temp[1,5], temp[1,6], index, sep = ""), length(temp1[[j]][1]))
      index <- index + 1
    }
    for (j in 1:dim(temp)[1]) {
      chr[which(chr[,4] == temp[j,4]),8] <- temp[j,8]
    }
  }

  chr_Band <- as.data.frame(array(NA,dim = c(length(unique(chr[,"chrName"])),8)))
  colnames(chr_Band) <- colnames(chr)
  chr_Band[,8] <- unique(chr[,"chrName"])
  for (i in 1:dim(chr_Band)[1]) {
    chr_Band[i,c(1,4,5,6)] <- chr[which(chr[,"chrName"] == chr_Band[i,8])[1],c(1,4,5,6)]
    chr_Band[i,2] <- min(chr[which(chr[,8] == chr_Band[i,8]),2])
    chr_Band[i,3] <- max(chr[which(chr[,8] == chr_Band[i,8]),3])
    chr_Band[i,7] <- as.integer(chr_Band[i,3]) - as.integer(chr_Band[i,2])
  }

  for (i in 1:24) {
    index_C <- max(which(chr_Band[,5] == i & chr_Band[,6] == "p"))
    if(index_C > 1 & chr_Band[index_C-1,5] == chr_Band[index_C,5] & chr_Band[index_C-1,6] == chr_Band[index_C,6] &
       chr_Band[index_C,7] <= r1*r*mean_L ){
      chr_Band[index_C,8] <- chr_Band[index_C-1,8]
    }
    index_C <- max(which(chr_Band[,5] == i & chr_Band[,6] == "q"))
    if(index_C > 1 & chr_Band[index_C-1,5] == chr_Band[index_C,5] & chr_Band[index_C-1,6] == chr_Band[index_C,6] &
       chr_Band[index_C,7] <= r1*r*mean_L){
      chr_Band[index_C,8] <- chr_Band[index_C-1,8]
    }
  }

  chr <- chr_Band
  chr_Band <- as.data.frame(array(NA,dim = c(length(unique(chr[,"chrName"])),8)))
  colnames(chr_Band) <- colnames(chr)
  chr_Band[,8] <- unique(chr[,"chrName"])
  for (i in 1:dim(chr_Band)[1]) {
    chr_Band[i,c(1,4,5,6)] <- chr[which(chr[,"chrName"] == chr_Band[i,8])[1],c(1,4,5,6)]
    chr_Band[i,2] <- min(chr[which(chr[,8] == chr_Band[i,8]),2])
    chr_Band[i,3] <- max(chr[which(chr[,8] == chr_Band[i,8]),3])
    chr_Band[i,7] <- as.integer(chr_Band[i,3]) - as.integer(chr_Band[i,2])
  }

  return(chr_Band)
}

chr <- chr_deal(2,0.75)
write.csv(chr, "./data/Genom/chromData.csv")
chr1 <- read.csv("./data/Genom/chromData.csv", row.names = 1)
chr <- read.csv("./data/Genom/GRCH37_chrBand1.csv")

chr1 <- chr1[which(chr1[,5] < 23),]
for (cancer_temp in names(ICGC_drivers)) {
  file1 <- paste("./data/ICGC/",cancer_temp,"/",sep = "")
  if(!dir.exists(file1)){
    dir.create(file1)
  }
  driver_chr <- ICGC_drivers[[cancer_temp]]$driver_chr

  for (i in 1:length(ICGC_drivers[[cancer_temp]]$driver_chr)) {
    kind <- substr(driver_chr[i],1,1)
    remain <- substr(driver_chr[i],2, nchar(driver_chr[i]))

    if(remain %in% ID[,7]){
      chr_index <- ID[which(ID[,7] == remain), 1]
      chr_start <- ID[which(ID[,7] == remain), 2]
      chr_end <- ID[which(ID[,7] == remain), 3]
    }else{
      if_Dot <- FALSE
      for (j in 1:nchar(remain)) {
        if(substr(remain, j, j) == "."){
          if_Dot = TRUE
          break
        }
      }
      if(if_Dot){
        index_C <- 1
        while (!(substr(remain,index_C,index_C) %in% c("p","q"))) {
          index_C <- index_C+1
        }
        chr_index <- substr(remain,1,index_C-1)
        chr_start <- chr[which(chr[,5] == substr(remain,1,index_C-1) & chr[,4] == substr(remain,index_C,nchar(remain))),2]
        chr_end <- chr[which(chr[,5] == substr(remain,1,index_C-1) & chr[,4] == substr(remain,index_C,nchar(remain))),3]
      }else{
        if(grepl("p",remain) | grepl("q",remain)){
          index_C <- 1
          while (!(substr(remain,index_C,index_C) %in% c("p","q"))) {
            index_C <- index_C+1
          }
          chr_index <- substr(remain,1,index_C-1)
          if(substr(remain, nchar(remain), nchar(remain)) %in% c("p","q")){
            chr_start <- chr[which(chr[,5] == substr(remain, 1, nchar(remain)-1) &
                                     chr[,6] == substr(remain, nchar(remain), nchar(remain))),2]
            chr_end <- chr[which(chr[,5] == substr(remain, 1, nchar(remain)-1) &
                                   chr[,6] == substr(remain, nchar(remain), nchar(remain))),3]
          }else{
            chr_start <- chr[which(chr[,8] == remain ),2]
            chr_end <- chr[which(chr[,8] == remain ),3]
          }
        }else{
          chr_index <- remain
          chr_start <- chr[which(chr[,5] == remain ),2]
          chr_end <- chr[which(chr[,5] == remain ),3]
        }
      }
    }
    driver_chr[i] <- NA
    for (j in 1:length(chr_start)) {
      remain <- chr1[which(chr1[,5] == chr_index & chr_start[j] >= chr1[,2] & chr_end[j] <= chr1[,3]),8]
      driver_chr <- c(driver_chr, paste(kind, remain, sep = ""))
    }
  }
  driver_chr <- as.character(na.omit(driver_chr))
  driver_chr <- unique(driver_chr)
  driver_chr <- driver_chr[which(nchar(driver_chr) > 1)]

  ICGC_drivers[[cancer_temp]]$driver_chr_row <- ICGC_drivers[[cancer_temp]]$driver_chr
  ICGC_drivers[[cancer_temp]]$driver_chr <- driver_chr
}
saveRDS(ICGC_drivers, "./data/ICGC/ICGC_drivers1.rds")
ICGC_drivers <- readRDS("./data/ICGC/ICGC_drivers1.rds")


################# ICGC Genotype
ICGC_Mutation <- as.data.frame(fread("./data/ICGC/pancan_pcawg_2020/data_mutations.txt", skip = 2, fill = TRUE, sep = "\t"))

mclapply(names(ICGC_drivers), function(cancer_temp){
  file1 <- paste("./data/ICGC/",cancer_temp,"/",sep = "")

  cancer_chr_data <- read.csv(paste(file1, "cancer_chr_data.csv", sep = ""))

  donor_temp <- intersect(cancer_chr_data[,22], sample[,4])

  ICGC_Cli_sample_temp <- ICGC_Cli_sample[which(ICGC_Cli_sample[,8] == cancer_temp),]

  donor_temp <- intersect(donor_temp, unique(ICGC_Cli_sample_temp[which(!is.na(ICGC_Cli_sample_temp[,8])),2]))
  sample_temp <- ICGC_Cli_sample_temp[which(ICGC_Cli_sample_temp[,2] %in% donor_temp & !is.na(ICGC_Cli_sample_temp[,8])), 1]
  sample_temp <- intersect(sample_temp, sample[,8])
  donor_temp <- ICGC_Cli_sample_temp[which(ICGC_Cli_sample_temp[,1] %in% sample_temp & !is.na(ICGC_Cli_sample_temp[,8])), 2]
  sample_temp <- ICGC_Cli_sample_temp[which(ICGC_Cli_sample_temp[,2] %in% donor_temp & !is.na(ICGC_Cli_sample_temp[,8])), 1]

  driver_chr <- ICGC_drivers[[cancer_temp]]$driver_chr
  driver_genes <- ICGC_drivers[[cancer_temp]]$driver_gene

  Genotype <- array(NA,dim = c(length(donor_temp) ,length(driver_genes)+length(driver_chr)+2))
  Genotype[,dim(Genotype)[2]] <- donor_temp
  colnames(Genotype) <- c(driver_genes,driver_chr,"dataSet","sample")
  Genotype[,dim(Genotype)[2]-1] <- "ICGC"

  for (i in 1:dim(Genotype)[1]) {
    sample_temp1 <- ICGC_Cli_sample_temp[which(ICGC_Cli_sample_temp[,2] == Genotype[i, dim(Genotype)[2]] & !is.na(ICGC_Cli_sample_temp[,8])), 1]
    ICGC_Mutation_temp <- ICGC_Mutation[which(ICGC_Mutation[,17] %in% sample_temp1), 1]
    Genotype[i,which(driver_genes %in% ICGC_Mutation_temp)] <- 1
    Genotype[i,which(!(driver_genes %in% ICGC_Mutation_temp))] <- 0

    for (j in paste("./data/ICGC/cna/",sample[which(sample[,4] == Genotype[i,dim(Genotype)[2]]),6],".consensus.20170119.somatic.cna.annotated.txt",sep = "")) {
      if(file.exists(j)){
        chr_temp <- as.data.frame(fread(j))[,1:4]
        for (k in (length(driver_genes)+1):(length(driver_genes)+length(driver_chr))) {
          chrEvent_temp <- colnames(Genotype)[k]
          chrEvent_temp1 <- substr(chrEvent_temp, 2, nchar(chrEvent_temp))
          chr_start <- as.integer(chr1[which(chr1[,8] == chrEvent_temp1) ,2])
          chr_end <- as.integer(chr1[which(chr1[,8] == chrEvent_temp1) ,3])
          chr_name <- chr1[which(chr1[,8] == chrEvent_temp1) ,5]

          if(substr(chrEvent_temp, 1, 1) == "-"){
            if (length(which(chr_temp[,1] == chr_name & ((chr_temp[,2] >= chr_start & chr_temp[,2] <= chr_end)|(
              chr_temp[,3] >= chr_start & chr_temp[,3] <= chr_end)) & chr_temp[,4] <= 1))>0){
              Genotype[i,k] <- 1
            }else{
              Genotype[i,k] <- 0
            }
          }else{
            if (length(which(chr_temp[,1] == chr_name & ((chr_temp[,2] >= chr_start & chr_temp[,2] <= chr_end)|(
              chr_temp[,3] >= chr_start & chr_temp[,3] <= chr_end)) & chr_temp[,4] >= 3))>0){
              Genotype[i,k] <- 1
            }else{
              Genotype[i,k] <- 0
            }
          }
        }
        break
      }
    }
  }
  write.csv(Genotype, paste(file1, "Genotype.csv", sep = ""), row.names = FALSE)

}, mc.cores = length(ICGC_drivers), mc.cleanup = TRUE)

### create results dir
{
  res_file <- "./data/Results/"
  if(!dir.exists(res_file)){
    dir.create(res_file)
  }
  for (cancer_index in 1:length(ICGC_drivers)) {
    cancer_temp <- names(ICGC_drivers)[cancer_index]
    res_cancer_file <- paste(res_file,cancer_temp,"/",sep = "")
    if(!dir.exists(res_cancer_file)){
      dir.create(res_cancer_file)
    }
  }

  ID <- read.csv("./data/Genom/ID_37.csv")
  chr1 <- read.csv("./data/Genom/chromData.csv", row.names = 1)

  cBioPortal_file <- "./data/cBioPortal/"
  ch <- import.chain('./data/Genom/hg38ToHg19.over.chain')
}

### parallel create Genotypes
{
  library(future)
  library(future.apply)
  plan(multisession, workers = length(ICGC_drivers))

  undeal_cancer_res <- future_lapply(1:length(ICGC_drivers), function(cancer_index) {

    undeal_cancer <- list()

    cancer_temp <- names(ICGC_drivers)[cancer_index]
    cancer_file <- paste(cBioPortal_file,cancer_temp,"/",sep = "")
    cancer_files <- list.files(cancer_file)
    res_cancer_file <- paste(res_file,cancer_temp,"/",sep = "")

    driver_chr <- ICGC_drivers[[cancer_index]]$driver_chr
    driver_genes <- ICGC_drivers[[cancer_index]]$driver_gene
    chr_Gene_index <- list()
    for (i in 1:length(driver_chr)) {
      chrEvent_temp1 <- substr(driver_chr[i], 2, nchar(driver_chr[i]))
      chr_start <- chr1[which(chr1[,8] == chrEvent_temp1) ,2]
      chr_end <- chr1[which(chr1[,8] == chrEvent_temp1) ,3]
      chr_name <- chr1[which(chr1[,8] == chrEvent_temp1) ,5]
      ID_temp <- ID[which(ID[,1] == chr_name),]
      chr_Gene_index[[i]] <- ID_temp[which(ID_temp[,2] >= chr_start & ID_temp[,3] <= chr_end),7]
      names(chr_Gene_index)[i] <- driver_chr[i]
    }
    saveRDS(chr_Gene_index, paste(res_cancer_file, "chr_Gene_index.rds", sep = ""))
    chr_Gene_index <- readRDS(paste(res_cancer_file, "chr_Gene_index.rds", sep = ""))

    for (cancer_file_temp_index in 1:length(cancer_files)) {
      cancer_file_temp <- cancer_files[cancer_file_temp_index]
      cancer_dir_temp <- paste(cancer_file, gsub(".tar.gz", "", cancer_file_temp), "/", sep = "")

      if(grepl(".tar.gz", cancer_file_temp) & !dir.exists(cancer_dir_temp)){
        system(paste("tar -zxvf", paste(cancer_file, cancer_file_temp, sep = ""), "-C", cancer_file, sep = " "))
      }

      if(grepl(".tar.gz", cancer_file_temp) & dir.exists(cancer_dir_temp) ){
        cat(cancer_temp, ":", gsub(".tar.gz", "", cancer_file_temp), "\n")
        ## Genotype
        mutation <- read.table(paste(cancer_dir_temp, "data_mutations.txt", sep = ""), comment.char = "#", quote = "",
                               header = T, fill = TRUE, check.names = FALSE, sep = "\t")

        cli_sample <- read.table(paste(cancer_dir_temp, "data_clinical_sample.txt", sep = ""), comment.char = "#", quote = "",
                                 header = T, fill = TRUE, check.names = FALSE, sep = "\t")

        if("CANCER_TYPE_DETAILED" %in% colnames(cli_sample) & "SAMPLE_ID" %in% colnames(cli_sample) &
           "Tumor_Sample_Barcode" %in% colnames(mutation)){
          # sample_temp <- intersect(mutation[,"Tumor_Sample_Barcode"], colnames(cnaseg))
          sample_temp <- mutation[,"Tumor_Sample_Barcode"]
          sample_temp <- unique(intersect(sample_temp, cli_sample[,"SAMPLE_ID"]))
          donor_temp <- unique(cli_sample[which(cli_sample[,"SAMPLE_ID"] %in% sample_temp) ,"PATIENT_ID"])
          sample_temp <- unique(cli_sample[which(cli_sample[,"PATIENT_ID"] %in% donor_temp) ,"SAMPLE_ID"])

          Genotype <- array(NA,dim = c(length(donor_temp) ,length(driver_genes)+length(driver_chr)+3))
          colnames(Genotype) <- c(driver_genes,driver_chr,"dataSet","sample","cancer_type")
          Genotype[, "sample"] <- donor_temp
          Genotype[, "dataSet"] <- gsub(".tar.gz", "", cancer_file_temp)

          for (i in 1:dim(Genotype)[1]) {
            sample_temp1 <- cli_sample[which(cli_sample[,"PATIENT_ID"] == Genotype[i, "sample"]), "SAMPLE_ID"]
            Mutation_temp <- mutation[which(mutation[, "Tumor_Sample_Barcode"] %in% sample_temp1), "Hugo_Symbol"]
            Genotype[i,which(driver_genes %in% Mutation_temp)] <- 1
            Genotype[i,which(!(driver_genes %in% Mutation_temp))] <- 0
            Genotype[i, "cancer_type"] <- paste(cli_sample[which(cli_sample[,"PATIENT_ID"] == Genotype[i, "sample"]),
                                                           "CANCER_TYPE_DETAILED"], collapse = "|")
          }

          if(file.exists(paste(cancer_dir_temp, "data_cna_hg19.seg", sep = "")) | file.exists(paste(cancer_dir_temp, "data_cna_hg38.seg", sep = ""))){
            if(file.exists(paste(cancer_dir_temp, "data_cna_hg38.seg", sep = ""))){
              cnaseg <- read.table(paste(cancer_dir_temp, "data_cna_hg38.seg", sep = ""), comment.char = "#", quote = "",
                                   header = T, fill = TRUE, check.names = FALSE, sep = "\t")
              cnaseg$"chr_temp" <- paste("chr",cnaseg[,2],sep = "")
              hg38.gr <- GRanges(seqnames = cnaseg[,"chr_temp"], ranges=IRanges(start=cnaseg[,3], end = cnaseg[,4]))
              hg19.gr <- liftOver(hg38.gr, ch)
              cnaseg$"start_temp" <- cnaseg[,3]
              cnaseg$"end_temp" <- cnaseg[,4]

              hg19.gr_results <- mclapply(1:length(hg19.gr), function(j){
                hg19.gr_temp <- hg19.gr[[j]]
                hg19.gr_temp_start <- hg19.gr_temp@ranges@start[1]
                hg19.gr_temp_end <- hg19.gr_temp@ranges@start[length(hg19.gr_temp@ranges)] +
                  hg19.gr_temp@ranges@width[length(hg19.gr_temp@ranges)] - 1
                return(c(hg19.gr_temp_start, hg19.gr_temp_end))
              }, mc.cores = 10, mc.cleanup = TRUE)

              cnaseg[,3] <- sapply(hg19.gr_results, function(x){
                x[1]
              })
              cnaseg[,4] <- sapply(hg19.gr_results, function(x){
                x[2]
              })
              cnaseg <- cnaseg[,1:6]

            }else{
              cnaseg <- read.table(paste(cancer_dir_temp, "data_cna_hg19.seg", sep = ""), comment.char = "#", quote = "",
                                   header = T, fill = TRUE, check.names = FALSE, sep = "\t")
            }
            for (j in c(3,4)) {
              cnaseg[,j] <- as.integer(cnaseg[,j])
            }
            cnaseg[,6] <- round(2*(2^(cnaseg[,6])))

            for (i in 1:dim(Genotype)[1]) {
              if(length(which(cnaseg[,1] %in% cli_sample[which(cli_sample[,"PATIENT_ID"] == Genotype[i, "sample"]),
                                                         "SAMPLE_ID"])) > 0){
                cnaseg_temp  <- cnaseg[which(cnaseg[,1] %in% cli_sample[which(cli_sample[,"PATIENT_ID"] == Genotype[i, "sample"]),
                                                                        "SAMPLE_ID"]), c(2,3,4,6)]

                for (k in (length(driver_genes)+1):(length(driver_genes)+length(driver_chr))) {
                  chrEvent_temp <- colnames(Genotype)[k]
                  chrEvent_temp1 <- substr(chrEvent_temp, 2, nchar(chrEvent_temp))
                  chr_start <- as.integer(chr1[which(chr1[,8] == chrEvent_temp1) ,2])
                  chr_end <- as.integer(chr1[which(chr1[,8] == chrEvent_temp1) ,3])
                  chr_name <- chr1[which(chr1[,8] == chrEvent_temp1) ,5]

                  if(substr(chrEvent_temp, 1, 1) == "-"){
                    if (length(which(cnaseg_temp[,1] == chr_name & ((cnaseg_temp[,2] >= chr_start & cnaseg_temp[,2] <= chr_end)|(
                      cnaseg_temp[,3] >= chr_start & cnaseg_temp[,3] <= chr_end)) & cnaseg_temp[,4] <= 1))>0){
                      Genotype[i,k] <- 1
                    }else{
                      Genotype[i,k] <- 0
                    }
                  }else{
                    if (length(which(cnaseg_temp[,1] == chr_name & ((cnaseg_temp[,2] >= chr_start & cnaseg_temp[,2] <= chr_end)|(
                      cnaseg_temp[,3] >= chr_start & cnaseg_temp[,3] <= chr_end)) & cnaseg_temp[,4] >= 3)) > 0){
                      Genotype[i,k] <- 1
                    }else{
                      Genotype[i,k] <- 0
                    }
                  }
                }
              }
            }

          }else if(file.exists(paste(cancer_dir_temp, "data_cna.txt", sep = ""))){
            cnaseg <- read.table(paste(cancer_dir_temp, "data_cna.txt", sep = ""), comment.char = "#", quote = "",
                                 header = T, fill = TRUE, check.names = FALSE, sep = "\t")

            if(!"Hugo_Symbol" %in% colnames(cnaseg)){
              ENTREZID <- select(org.Hs.eg.db, keys = as.character(cnaseg[,"Entrez_Gene_Id"]),
                                 columns = "SYMBOL", keytype = "ENTREZID")
              cnaseg$"Hugo_Symbol" <- ENTREZID[,"SYMBOL"]
            }

            for (i in 1:dim(Genotype)[1]) {
              if(length(intersect(colnames(cnaseg), cli_sample[which(cli_sample[,"PATIENT_ID"] == Genotype[i, "sample"]),
                                                               "SAMPLE_ID"])) > 0){

                cnaseg_temp  <- cnaseg[, c("Hugo_Symbol", intersect(colnames(cnaseg), cli_sample[
                  which(cli_sample[,"PATIENT_ID"] == Genotype[i, "sample"]),"SAMPLE_ID"]))]
                for (k in (length(driver_genes)+1):(length(driver_genes)+length(driver_chr))) {
                  chr_genes <- chr_Gene_index[[colnames(Genotype)[k]]]
                  chr_genes_cna <- as.integer(unlist(cnaseg_temp[which(cnaseg_temp[,1] %in% chr_genes), 2:dim(cnaseg_temp)[2]]))
                  chr_genes_cna <- na.omit(chr_genes_cna)
                  if( length(chr_genes_cna) > 0 ){
                    uniqv <- unique(chr_genes_cna)
                    chr_genes_cna <- uniqv[which.max(tabulate(match(chr_genes_cna, uniqv)))]
                    chrEvent_temp <- colnames(Genotype)[k]

                    if(substr(chrEvent_temp, 1, 1) == "-"){
                      if (chr_genes_cna < 0){
                        Genotype[i,k] <- 1
                      }else{
                        Genotype[i,k] <- 0
                      }
                    }else{
                      if (chr_genes_cna > 0){
                        Genotype[i,k] <- 1
                      }else{
                        Genotype[i,k] <- 0
                      }
                    }
                  }
                }
              }
            }
          }

          Genotype <- na.omit(Genotype)
          write.csv(Genotype, paste(res_cancer_file, gsub(".tar.gz", "", cancer_file_temp),"_Genotype.csv", sep = ""), row.names = FALSE)
        }else{
          undeal_cancer <- c(undeal_cancer, gsub(".tar.gz", "", cancer_file_temp))
        }
      }
    }

    return(undeal_cancer)

  })

  plan(sequential)
  print(undeal_cancer_res)

}

### Genotypes Filter
{
  for (cancer_index in 1:length(ICGC_drivers)) {
    cancer_temp <- names(ICGC_drivers)[cancer_index]
    cancer_file <- paste(cBioPortal_file,cancer_temp,"/",sep = "")
    res_cancer_file <- paste(res_file,cancer_temp,"/",sep = "")
    cancer_files <- list.files(res_cancer_file)
    index <- 1
    for (i in 1:length(cancer_files)) {
      if( grepl("_Genotype.csv", cancer_files[i])){
        temp <- read.csv(paste(res_cancer_file, cancer_files[i], sep = ""),check.names = F)
        if(index ==1){
          Genotype <- temp
          index <- index + 1
        }else{
          Genotype <- rbind(Genotype, temp)
        }
      }
    }
    write.csv(Genotype, paste(res_cancer_file, "Genotype_temp.csv", sep = ""), row.names = F)
  }

  for (cancer_index in 1:length(ICGC_drivers)) {
    cancer_temp <- names(ICGC_drivers)[cancer_index]
    res_cancer_file <- paste(res_file,cancer_temp,"/",sep = "")
    Genotype <- read.csv(paste(res_cancer_file, "Genotype_temp.csv", sep = ""), check.names = F)

    temp <- Genotype[, "cancer_type", drop = F]
    temp <- unique(temp)

    if(cancer_temp == "Breast-AdenoCA"){
      Genotype <- Genotype[which(!Genotype[,"cancer_type"] %in% c("Breast Mixed Ductal and Lobular Carcinoma",
                                                                  "Breast Mixed Ductal and Lobular Carcinoma|Breast Mixed Ductal and Lobular Carcinoma|Breast Mixed Ductal and Lobular Carcinoma",
                                                                  "Breast Mixed Ductal and Lobular Carcinoma|Breast Mixed Ductal and Lobular Carcinoma","Breast",
                                                                  "Metaplastic Breast Cancer","Paget Disease of the Nipple","Adenoid Cystic Breast Cancer",
                                                                  "Solid Papillary Carcinoma of the Breast","Malignant Phyllodes Tumor of the Breast",
                                                                  "Breast Mixed Ductal and Lobular Carcinoma|Breast Mixed Ductal and Lobular Carcinoma",
                                                                  "Adenoid Cystic Breast Cancer","Breast|Breast|Breast","Adenoid Cystic Breast Cancer")),]
    }else if(cancer_temp == "CNS-GBM"){
      Genotype <- Genotype[which(!Genotype[,"cancer_type"] %in% c("Oligodendroglioma|Oligodendroglioma",
                                                                  "Diffuse Astrocytoma|Diffuse Astrocytoma","Anaplastic Astrocytoma|Diffuse Astrocytoma",
                                                                  "Anaplastic Oligodendroglioma|Oligodendroglioma","Anaplastic Astrocytoma|Anaplastic Astrocytoma",
                                                                  "Anaplastic Astrocytoma|Anaplastic Astrocytoma|Anaplastic Astrocytoma|Diffuse Astrocytoma",
                                                                  "Diffuse Astrocytoma|Anaplastic Astrocytoma|Anaplastic Astrocytoma",
                                                                  "Diffuse Astrocytoma|Diffuse Astrocytoma|Diffuse Astrocytoma","Diffuse Glioma",
                                                                  "Diffuse Glioma|Diffuse Glioma","Anaplastic Oligodendroglioma|Anaplastic Oligodendroglioma",
                                                                  "Glioma","Anaplastic Astrocytoma","Anaplastic Astrocytoma|Anaplastic Astrocytoma|Anaplastic Astrocytoma",
                                                                  "Diffuse Astrocytoma","Anaplastic Oligodendroglioma","Oligodendroglioma","Oligoastrocytoma|Oligoastrocytoma",
                                                                  "Diffuse Glioma|Diffuse Glioma|Diffuse Glioma|Diffuse Glioma")),]

    }else if(cancer_temp == "ColoRect-AdenoCA"){
      Genotype <- Genotype[which(!Genotype[,"cancer_type"] %in% c("Small Bowel Cancer",
                                                                  "Duodenal Adenocarcinoma","Small Intestinal Carcinoma")),]
    }else if(cancer_temp == "Liver-HCC"){
      Genotype <- Genotype[which(!Genotype[,"cancer_type"] %in% c("Fibrolamellar Carcinoma")),]
    }else if(cancer_temp == "Lung-AdenoCA"){
      Genotype <- Genotype[which(!Genotype[,"cancer_type"] %in% c("Small Cell Lung Cancer|Small Cell Lung Cancer",
                                                                  "Combined Small Cell Lung Carcinoma|Large Cell Neuroendocrine Carcinoma",
                                                                  "Non-Small Cell Lung Cancer|Non-Small Cell Lung Cancer|Non-Small Cell Lung Cancer|Non-Small Cell Lung Cancer|Non-Small Cell Lung Cancer|Non-Small Cell Lung Cancer|Non-Small Cell Lung Cancer|Non-Small Cell Lung Cancer|Large Cell Neuroendocrine Carcinoma",
                                                                  "Non-Small Cell Lung Cancer|Large Cell Neuroendocrine Carcinoma","Non-Small Cell Lung Cancer|Non-Small Cell Lung Cancer",
                                                                  "Non-Small Cell Lung Cancer|Malignant Tumor","Non-Small Cell Lung Cancer|Adenocarcinoma, NOS",
                                                                  "Non-Small Cell Lung Cancer|Cancer of Unknown Primary, NOS","Non-Small Cell Lung Cancer|Combined Small Cell Lung Carcinoma",
                                                                  "Lung Squamous Cell Carcinoma|Non-Small Cell Lung Cancer",
                                                                  "Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma","Non-Small Cell Lung Cancer|Poorly Differentiated Non-Small Cell Lung Cancer",
                                                                  "Non-Small Cell Lung Cancer|Poorly Differentiated Carcinoma, NOS","Lung Squamous Cell Carcinoma|Squamous Cell Carcinoma, NOS",
                                                                  "Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma",
                                                                  "Lung Squamous Cell Carcinoma|Poorly Differentiated Non-Small Cell Lung Cancer",
                                                                  "Non-Small Cell Lung Cancer|Non-Small Cell Lung Cancer|Non-Small Cell Lung Cancer|Pleomorphic Carcinoma of the Lung",
                                                                  "Non-Small Cell Lung Cancer|Cancer of Unknown Primary","Non-Small Cell Lung Cancer|Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma",
                                                                  "Non-Small Cell Lung Cancer|Combined Small Cell Lung Carcinoma|Large Cell Neuroendocrine Carcinoma",
                                                                  "Non-Small Cell Lung Cancer|Non-Small Cell Lung Cancer|Non-Small Cell Lung Cancer",
                                                                  "Non-Small Cell Lung Cancer|Sarcoma, NOS","Non-Small Cell Lung Cancer|Prostate Adenocarcinoma",
                                                                  "Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma",
                                                                  "Lung Squamous Cell Carcinoma|Thymic Carcinoma|Thymic Carcinoma|Lung Squamous Cell Carcinoma",
                                                                  "Lung Squamous Cell Carcinoma|Breast Invasive Ductal Carcinoma|Breast Invasive Ductal Carcinoma",
                                                                  "Non-Small Cell Lung Cancer|Non-Small Cell Lung Cancer|Large Cell Neuroendocrine Carcinoma",
                                                                  "Non-Small Cell Lung Cancer|Poorly Differentiated Carcinoma, NOS|Cancer of Unknown Primary",
                                                                  "Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma|Non-Small Cell Lung Cancer",
                                                                  "Non-Small Cell Lung Cancer|Cancer of Unknown Primary, NOS|Cancer of Unknown Primary",
                                                                  "Lung Squamous Cell Carcinoma","Non-Small Cell Lung Cancer",
                                                                  "Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma|Oral Cavity Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma",
                                                                  "Non-Small Cell Lung Cancer|Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma",
                                                                  "Non-Small Cell Lung Cancer|Lung Neuroendocrine Tumor","Large Cell Neuroendocrine Carcinoma",
                                                                  "Non-Small Cell Lung Cancer|Lung Squamous Cell Carcinoma")),]


    }else if(cancer_temp == "Lung-SCC"){
      Genotype <- Genotype[which(Genotype[,"cancer_type"] %in% c(
        "Lung Squamous Cell Carcinoma|Non-Small Cell Lung Cancer",
        "Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma","Non-Small Cell Lung Cancer|Poorly Differentiated Non-Small Cell Lung Cancer",
        "Non-Small Cell Lung Cancer|Poorly Differentiated Carcinoma, NOS","Lung Squamous Cell Carcinoma|Squamous Cell Carcinoma, NOS",
        "Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma",
        "Lung Squamous Cell Carcinoma|Poorly Differentiated Non-Small Cell Lung Cancer",
        "Non-Small Cell Lung Cancer|Cancer of Unknown Primary","Non-Small Cell Lung Cancer|Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma",
        "Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma",
        "Lung Squamous Cell Carcinoma|Thymic Carcinoma|Thymic Carcinoma|Lung Squamous Cell Carcinoma",
        "Lung Squamous Cell Carcinoma|Breast Invasive Ductal Carcinoma|Breast Invasive Ductal Carcinoma",
        "Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma|Non-Small Cell Lung Cancer",
        "Lung Squamous Cell Carcinoma","Non-Small Cell Lung Cancer",
        "Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma|Oral Cavity Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma",
        "Non-Small Cell Lung Cancer|Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma|Lung Squamous Cell Carcinoma",
        "Non-Small Cell Lung Cancer|Lung Squamous Cell Carcinoma")),]

    }else if(cancer_temp == "Prost-AdenoCA"){
      Genotype <- Genotype[which(!Genotype[,"cancer_type"] %in% c("Prostate Neuroendocrine Carcinoma","Prostate Small Cell Carcinoma",
                                                                  "Prostate Squamous Cell Carcinoma","Prostate Small Cell Carcinoma|Prostate Small Cell Carcinoma",
                                                                  "Prostate Neuroendocrine Carcinoma|Prostate Neuroendocrine Carcinoma")),]
    }else if(cancer_temp == "Skin-Melanoma"){
      Genotype <- Genotype[which(!Genotype[,"cancer_type"] %in% c("[Not Available]","[Not Available]|[Not Available]",
                                                                  "Pheochromocytoma")),]
    }else if(cancer_temp == "Uterus-AdenoCA"){
      Genotype <- Genotype[which(!Genotype[,"cancer_type"] %in% c("Uterine Clear Cell Carcinoma",
                                                                  "Uterine Undifferentiated Carcinoma","Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor",
                                                                  "Uterine Mixed Endometrial Carcinoma","Uterine Neuroendocrine Carcinoma",
                                                                  "Uterine Dedifferentiated Carcinoma",
                                                                  "Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor|Uterine Carcinosarcoma/Uterine Malignant Mixed Mullerian Tumor")),]
    }

    Genotype <- Genotype[,1:(dim(Genotype)[2]-1)]

    write.csv(Genotype, paste(res_cancer_file, "Genotype_temp1.csv", sep = ""), row.names = F)

  }

  for (cancer_index in 1:length(ICGC_drivers)) {
    cancer_temp <- names(ICGC_drivers)[cancer_index]
    cat(cancer_temp,"\n")
    res_cancer_file <- paste(res_file,cancer_temp,"/",sep = "")
    Genotype <- read.csv(paste(res_cancer_file, "Genotype_temp1.csv", sep = ""), check.names = F)
    ICGC_Genotype <- read.csv(paste("./data/ICGC/", cancer_temp, "/Genotype.csv", sep = ""), check.names = F)

    if(all(colnames(Genotype) == colnames(ICGC_Genotype))){
      Genotype <- rbind(Genotype, ICGC_Genotype)

      Genotype[,"sample"] <- sapply(Genotype[,"sample"], function(x){
        if(grepl("TCGA",x)){
          return(substr(x,1,12))
        }else{
          return(x)
        }
      })

      driver_chr <- ICGC_drivers[[cancer_index]]$driver_chr
      driver_genes <- ICGC_drivers[[cancer_index]]$driver_gene

      Genotype1 <- as.data.frame(array(NA, dim = c(length(unique(Genotype[,"sample"])), dim(Genotype)[2])))
      colnames(Genotype1) <- colnames(Genotype)
      Genotype1[, "sample"] <- unique(Genotype[,"sample"])

      Genotype_res <- mclapply(1:dim(Genotype1)[1], function(i){
        cancer_sample <- Genotype1[i,"sample"]
        Genotype_temp <- Genotype[which(Genotype[,"sample"] == cancer_sample), c(driver_genes, driver_chr),drop = F]
        if(dim(Genotype_temp)[1] > 1){
          Genotype_temp1 <- colSums(Genotype_temp)/dim(Genotype_temp)[1]
          Genotype_temp1[Genotype_temp1 >= 0.5] <- 1
          Genotype_temp1[Genotype_temp1 < 0.5] <- 0
          Genotype_temp <- Genotype_temp1
        }

        return(list(Genotype_temp = Genotype_temp, sample = cancer_sample,
                    datsSet = paste(Genotype[which(Genotype[,"sample"] == cancer_sample), "dataSet"],collapse = "|")))
      }, mc.cores = 20, mc.cleanup = TRUE)

      for (i in 1:dim(Genotype1)[1]) {
        Genotype1[i, c(driver_genes, driver_chr)] <- Genotype_res[[i]]$Genotype_temp
        Genotype1[i, "dataSet"] <- Genotype_res[[i]]$datsSet
      }

      dataSet_map <- array(NA, dim = c(length(unique(Genotype[,"dataSet"])),2))
      dataSet_map[,1] <- unique(Genotype[,"dataSet"])
      rownames(dataSet_map) <- dataSet_map[,1]
      dataSet_list <- list()
      temp <- unique(Genotype1$dataSet)
      for (i in 1:length(temp)) {
        temp1 <- strsplit(temp[i], "[|]")[[1]]
        if(length(dataSet_list) == 0){
          dataSet_list[[length(dataSet_list) + 1]] <- temp1
        }else{
          if_union <- FALSE
          for (j in 1:length(dataSet_list)) {
            if(length(intersect(dataSet_list[[j]], temp1)) > 0){
              if_union <- TRUE
              dataSet_list[[j]] <- unique(c(dataSet_list[[j]], temp1))
              break
            }
          }
          if(!if_union){
            dataSet_list[[length(dataSet_list) + 1]] <- temp1
          }
        }
      }
      tcga_temp <- c()
      dataSet_info <- list()
      for (i in 1:length(dataSet_list)) {
        if(any(grepl("tcga", dataSet_list[[i]]))){
          tcga_temp <- unique(c(tcga_temp, dataSet_list[[i]]))
        }else{
          dataSet_info[[length(dataSet_info) + 1]] <-
            list(dataSet = dataSet_list[[i]], dataset_index = paste("dataSet", length(dataSet_info) + 1, sep = ""))
        }
      }
      if(length(tcga_temp) > 0){
        dataSet_info[[length(dataSet_info) + 1]] <-
          list(dataSet = tcga_temp, dataset_index = paste("dataSet", length(dataSet_info) + 1, sep = ""))
      }
      saveRDS(dataSet_info, paste(res_cancer_file, "dataSet_info.rds", sep = ""))

      for (i in 1:length(dataSet_info)) {
        for (temp in dataSet_info[[i]]$dataSet) {
          Genotype1[which(grepl(temp, Genotype1[,"dataSet"])), "dataSet"] <- dataSet_info[[i]]$dataset_index
        }
      }

      Genotype1 <- na.omit(Genotype1)
      write.csv(Genotype1, paste(res_cancer_file, "Genotype.csv", sep = ""), row.names = F)

    }else{
      cat(cancer_temp, " have errors!\n")
    }
  }


}









