## baselineTestExample createData, ct_cbn

# Note on Example Code
# The example code provided here has been reorganized for clarity and testing purposes.
# Please review it carefully and modify as necessary for your specific use case.
stop()

library("devtools")
library("parallel")
# clean_dll()
# document()
load_all()

## ct-cbn file
ct_cbn_file <- "~/program/ct-cbn-0.1.04b/h-cbn"

## generate simulate data
simulateFile <- "./data/simulateGenotypeData/"
while (!dir.exists(simulateFile)){
  dir.create(simulateFile, recursive = T)
}
# ceiling(sqrt(dim(mat)[1]*0.3))
numEventPerKind <- c(6, 8, 12, 16, 20)
numEventKind <- c(2, 3)
numSet <- c(2, 3)
epslion <- c(0.05, 0.1, 0.2)
rep_num <- 10
for (i in numEventPerKind) {
  numSamplePerSet <- ceiling((i*10/3)^2)
  for (j in numEventKind) {
    k <- j
    file1 <- paste(simulateFile,i,"-",j,"-",k,"/",sep = "")
    while (!dir.exists(file1)){
      dir.create(file1, recursive = T)
    }
    for (eps in epslion) {
      file2 <- paste(file1,eps,"/",sep = "")
      while (!dir.exists(file2)) {
        dir.create(file2, recursive = T)
      }
      for (rep in 1:rep_num) {
        file3 <- paste(file2,"rep",rep,".rds",sep = "")
        while (!file.exists(file3)) {
          temp <- simulation_Data_Generate( rep(i,j), rep(numSamplePerSet, k), baseEpsilon = eps)
          saveRDS(temp, file3)
        }
      }
    }
  }
}


## convert numEventPerKind, numEventKind, epslion to exp_res
exp_res <- array(NA, dim = c(length(numEventPerKind)*length(numEventKind)*length(epslion)*rep_num, 5))
colnames(exp_res) <- c("numP","eventKind","setKind","eps","rep")
index <- 1
for (i in numEventPerKind) {
  for (j in numEventKind) {
    for (eps in epslion) {
      for (rep in 1:rep_num) {
        exp_res[index,1] <- i
        exp_res[index,2] <- j
        exp_res[index,3] <- j
        exp_res[index,4] <- eps
        exp_res[index,5] <- rep
        index <- index + 1
      }
    }
  }
}

## create simulate res file
simulateResFile <- "./data/simulateResData/"
while (!dir.exists(simulateResFile)){
  dir.create(simulateResFile, recursive = T)
}

thrds <- 40

## parallel test h-cbn
temp1 <- mclapply(1:dim(exp_res)[1], function(index){
  i <- exp_res[index, 1]
  j <- exp_res[index, 2]
  k <- exp_res[index, 3]
  eps <- exp_res[index, 4]
  rep <- exp_res[index, 5]

  numSamplePerSet <- ceiling((i*10/3)^2)
  file1 <- paste(simulateFile,i,"-",j,"-",k,"/",sep = "")
  file2 <- paste(file1,eps,"/",sep = "")
  file3 <- paste(file2,"rep",rep,".rds",sep = "")

  res_file1 <- paste(simulateResFile,i,"-",j,"-",k,"/",sep = "")
  res_file2 <- paste(res_file1,eps,"/",sep = "")
  res_file3 <- paste(res_file2, "temp", rep, "/",sep = "")
  res_file4 <- paste(res_file3, "test", "/",sep = "")

  while (!dir.exists(res_file3)){
    dir.create(res_file3, recursive = T)
  }
  while (!dir.exists(res_file4)){
    dir.create(res_file4, recursive = T)
  }

  data <- readRDS(file3)
  h_cbn_data <- data$obs_events

  write(c(dim(h_cbn_data)[1], dim(h_cbn_data)[2] + 1), paste(res_file3, "test.pat", sep=""))
  write.table(cbind(rep(1, dim(h_cbn_data)[1]), h_cbn_data),
              paste(res_file3, "test.pat", sep=""),
              row.names=FALSE, col.names=FALSE, append=TRUE)

  poset_temp <- array(NA, dim = c(dim(h_cbn_data)[2] - 1,2))
  colsum_Data <- array(NA, dim = c(dim(h_cbn_data)[2],2))
  colsum_Data[,1] <- 1:dim(colsum_Data)[1]
  colsum_Data[,2] <- colSums(h_cbn_data)
  colsum_Data <- colsum_Data[order(colsum_Data[,2],decreasing = T),]
  colsum_Data_sort_index <- colsum_Data[,1]
  for (a in 1:(length(colsum_Data_sort_index) - 1)) {
    poset_temp[a,1] <- colsum_Data_sort_index[a]
    poset_temp[a,2] <- colsum_Data_sort_index[a + 1]
  }
  write(dim(h_cbn_data)[2], paste(res_file3, "test.poset", sep=""))
  write.table(poset_temp,paste(res_file3, "test.poset", sep=""), row.names=FALSE, col.names=FALSE, append=TRUE)
  write(0, paste(res_file3, "test.poset", sep=""), append=TRUE)

  system(paste("OMP_NUM_THREADS=1 ", ct_cbn_file, " -f ", gsub("test/","test",res_file4), "-w > ", gsub("test/","res.txt",res_file4)))

}, mc.cores = thrds, mc.cleanup = TRUE)


