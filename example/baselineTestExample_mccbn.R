## baselineTestExample_mccbn

# Note on Example Code
# The example code provided here has been reorganized for clarity and testing purposes.
# Please review it carefully and modify as necessary for your specific use case.
stop()

library("mccbn")
library("parallel")
## generate simulate data
simulateFile <- "./data/simulateGenotypeData/"
# ceiling(sqrt(dim(mat)[1]*0.3))
numEventPerKind <- c(6, 8, 12, 16, 20)
numEventKind <- c(2, 3)
numSet <- c(2, 3)
epslion <- c(0.05, 0.1, 0.2)
rep_num <- 10

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

simulateResFile <- "./data/simulateResData/"

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

  posets = candidate_posets(h_cbn_data, rep(1, dim(h_cbn_data)[1]), 0.9)
  poset0 = posets[[length(posets)]]
  fit = adaptive.simulated.annealing(poset0, h_cbn_data, L=100, max.iter.asa=10, seed=10L, thrds=4L)
  saveRDS(fit, paste(res_file3, "fit.rds", sep = ""))

}, mc.cores = 12, mc.cleanup = TRUE)


#########################################################################
## data scaling test
library("mccbn")
library("parallel")

## generate simulate data
simulateFile <- "./data/simulateGenotypeData1/"
while (!dir.exists(simulateFile)){
  dir.create(simulateFile, recursive = T)
}
# ceiling(sqrt(dim(mat)[1]*0.3))
numEventPerKind <- c(6, 8, 12, 16, 20)
numEventKind <- c(2)
numSet <- c(2)
epslion <- c(0.05)
rep_num <- 10
data_scaling <- c(1, 4, 8)

exp_res <- array(NA, dim = c(length(numEventPerKind)*length(data_scaling)*length(numEventKind)*length(epslion)*rep_num, 6))
colnames(exp_res) <- c("numP","eventKind","setKind","eps","scaling","rep")
index <- 1
for (i in numEventPerKind) {
  for (j in numEventKind) {
    for (s in data_scaling) {
      for (eps in epslion) {
        for (rep in 1:rep_num) {
          exp_res[index,1] <- i
          exp_res[index,2] <- j
          exp_res[index,3] <- j
          exp_res[index,4] <- eps
          exp_res[index,5] <- s
          exp_res[index,6] <- rep
          index <- index + 1
        }
      }
    }
  }
}

simulateResFile <- "./data/simulateResData1_mccbn/"

## parallel test mccbn
{
  for (index in 1:dim(exp_res)[1]) {
    i <- exp_res[index, 1]
    j <- exp_res[index, 2]
    k <- exp_res[index, 3]
    eps <- exp_res[index, 4]
    s <- exp_res[index, 5]
    rep <- exp_res[index, 6]

    numSamplePerSet <- s*ceiling((i*10/3)^2)
    file1 <- paste(simulateFile,i,"-",j,"-",k,"-",s,"/",sep = "")
    file2 <- paste(file1,eps,"/",sep = "")
    file3 <- paste(file2,"rep",rep,".rds",sep = "")

    res_file1 <- paste(simulateResFile,i,"-",j,"-",k,"-",s,"/",sep = "")
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

    posets = candidate_posets(h_cbn_data, rep(1, dim(h_cbn_data)[1]), 0.9)
    poset0 = posets[[length(posets)]]
    fit = adaptive.simulated.annealing(poset0, h_cbn_data, L=100, max.iter.asa=10, seed=10L, thrds=10L)
    saveRDS(fit, paste(res_file3, "fit.rds", sep = ""))
  }

}




