## baselineTestExample_este

# Note on Example Code
# The example code provided here has been reorganized for clarity and testing purposes.
# Please review it carefully and modify as necessary for your specific use case.
stop()

library("este")
library("parallel")

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
{
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

    if(!file.exists(paste(res_file3, "este_fit.rds", sep = ""))){
      data <- readRDS(file3)
      mat <- as.matrix(data$obs_events)
      isF <- matrix(as.integer(c(1)) ,nrow = k, ncol = j)
      setD <- as.data.frame(array(NA, dim = c(k,2)))
      for (l in 1:k) {
        setD[l,2] <- numSamplePerSet*l - 1
        setD[l,1] <- setD[l,2] - numSamplePerSet + 1
      }
      for (l in 1:dim(setD)[2]) {
        setD[,l] <- as.integer(setD[,l])
      }
      setD <- as.matrix(setD)

      numEvent_temp <- rep(i, j)
      eventD <- as.data.frame(array(NA, dim = c(j,2)))
      for (l in 1:j) {
        range_temp <- indexToRange(numEvent_temp,l)
        eventD[l,] <- range_temp
      }
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

      ncores1 <- 4
      ncores2 <- 5

      poset_np <- min(i*j, 8)

      eps <- estimate_Epsilon_ForMulti(pat = Geno, isF = isF, setD = setD, eventD = eventD, multi_thrds = 1,
                                       threshold1 = as.integer(8), n_p = as.integer(8),
                                       T = 10.0, N_iter = 200L, thrds=as.integer(ncores1))


      poset <- find_Poset_ForVote(pat = Geno, isF = isF, eps = eps, setD = setD, eventD = eventD,
                                  Fine_Tune_Num = as.integer(2L), vote_size = 5L, vote_threshold = 0.5,
                                  vote_thrds = 1L, threshold = 0.0001, threshold2 = 0.01, n_p = as.integer(poset_np), is_update_eps = FALSE,
                                  T = 10, N_iter = 200L, thrds = as.integer(ncores2))

      fit <- estimate_Lambda(pat = Geno, poset = poset, isF = isF, eps = eps, setD = setD, eventD = eventD,
                             lambdaS = 1.0, L = 100L, sampling = 'add-remove',
                             maxIter = 500L, updateStepSize=20L, tol=0.001, maxLambda=1e6,
                             neighborhoodDist=1L, is_update_eps = FALSE,
                             thrds=as.integer(ncores2))
      fit$lambdaS <- 1.0
      fit$poset <- poset

      saveRDS(fit, paste(res_file3, "este_fit.rds", sep = ""))

    }

  }, mc.cores = 10, mc.cleanup = TRUE)

}


#########################################################
## data scaling test
library("devtools")
library("parallel")

load_all()

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
for (i in numEventPerKind) {
  for (s in data_scaling) {
    numSamplePerSet <- s*ceiling((i*10/3)^2)
    for (j in numEventKind) {
      k <- j
      file1 <- paste(simulateFile,i,"-",j,"-",k,"-",s,"/",sep = "")
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


}

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

simulateResFile <- "./data/simulateResData1/"

## parallel test h-cbn
temp1 <- mclapply(1:dim(exp_res)[1], function(index){
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
  res_file5 <- paste(res_file3, "cache", "/",sep = "")

  while (!dir.exists(res_file3)){
    dir.create(res_file3, recursive = T)
  }
  while (!dir.exists(res_file4)){
    dir.create(res_file4, recursive = T)
  }


  if(!file.exists(paste(res_file3, "este_fit.rds", sep = ""))){
    while (!dir.exists(res_file5)){
      dir.create(res_file5, recursive = T)
    }

    exp_res[index, 6] <- "starting"
    write.csv(exp_res, paste(simulateFile, "exp_res.csv", sep = ""))

    data <- readRDS(file3)
    mat <- as.matrix(data$obs_events)
    isF <- matrix(as.integer(c(1)) ,nrow = k, ncol = j)
    setD <- as.data.frame(array(NA, dim = c(k,2)))
    for (l in 1:k) {
      setD[l,2] <- numSamplePerSet*l - 1
      setD[l,1] <- setD[l,2] - numSamplePerSet + 1
    }
    for (l in 1:dim(setD)[2]) {
      setD[,l] <- as.integer(setD[,l])
    }
    setD <- as.matrix(setD)

    numEvent_temp <- rep(i, j)
    eventD <- as.data.frame(array(NA, dim = c(j,2)))
    for (l in 1:j) {
      range_temp <- indexToRange(numEvent_temp,l)
      eventD[l,] <- range_temp
    }
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

    ncores1 <- 6
    ncores2 <- 10

    poset_np <- min(i*j, 8)

    if(!file.exists(paste(res_file5, "eps.rds", sep = ""))){
      eps <- estimate_Epsilon_ForMulti(pat = Geno, isF = isF, setD = setD, eventD = eventD, multi_thrds = 1,
                                       threshold1 = as.integer(8), n_p = as.integer(8),
                                       T = 10.0, N_iter = 200L, thrds=as.integer(ncores1))
      saveRDS(eps, paste(res_file5, "eps.rds", sep = ""))
    }
    eps <- readRDS(paste(res_file5, "eps.rds", sep = ""))

    poset <- find_Poset_ForVote(pat = Geno, isF = isF, eps = eps, setD = setD, eventD = eventD,
                                Fine_Tune_Num = as.integer(2L), vote_size = 5L, vote_threshold = 0.5,
                                vote_thrds = 1L, threshold = 0.0001, threshold2 = 0.01, n_p = as.integer(poset_np), is_update_eps = FALSE,
                                T = 10, N_iter = 200L, thrds = as.integer(ncores2))


    fit <- estimate_Lambda(pat = Geno, poset = poset, isF = isF, eps = eps, setD = setD, eventD = eventD,
                           lambdaS = 1.0, L = 100L, sampling = 'add-remove',
                           maxIter = 500L, updateStepSize=20L, tol=0.001, maxLambda=1e6,
                           neighborhoodDist=1L, is_update_eps = FALSE,
                           thrds=as.integer(ncores2))

    fit$lambdaS <- 1.0
    fit$poset <- poset

    saveRDS(fit, paste(res_file3, "este_fit.rds", sep = ""))

  }else{
    exp_res[index, 6] <- "finished"
    write.csv(exp_res, paste(simulateFile, "exp_res.csv", sep = ""))
  }

}, mc.cores = 6, mc.cleanup = TRUE)


