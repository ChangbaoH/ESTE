## baseline test

#' @title Generate simulation data
#' @export
#' @description Generate simulation data
#' @param numEventArray A array describe the event number in each event type
#' @param numSampleArray A array describe the set number in each set type
#' @param lambdaSampling The value of lambdas. Default: 1
#' @param lambdaSamplingScaling The value of scaling which used to sample lambdas. Default: 3
#' @param baseEpsilon The value of baseEpsilon. Default: 0.05
#' @param epsilonSamplingRate The value of rate which used to sample epsilon. Default: 0.5
#' @param graphDensity The value which used to control the sparsity of poset. Default: 0.2
simulation_Data_Generate <- function(numEventArray, numSampleArray, lambdaSampling = 1, lambdaSamplingScaling = 3,
                                     baseEpsilon = 0.05, epsilonSamplingRate = 0.5, graphDensity = 0.2){

  numEventKind <- length(numEventArray)
  numSet <- length(numSampleArray)

  if(epsilonSamplingRate >= 1){
    cat("Error: epsilonSamplingRate >= 1\n")
    return(NULL)
  }

  # eps = matrix(runif(numEventKind*numSet, baseEpsilon - baseEpsilon*epsilonSamplingRate, baseEpsilon + baseEpsilon*epsilonSamplingRate),
  #              ncol = numEventKind, byrow = T)

  eps_min <- baseEpsilon - baseEpsilon*epsilonSamplingRate
  eps_max <- baseEpsilon + baseEpsilon*epsilonSamplingRate
  eps <- c(eps_min, eps_min + (1:(numEventKind*numSet - 1))*(eps_max - eps_min)/(numEventKind*numSet - 1))
  eps = matrix(eps, ncol = numEventKind, byrow = T)

  lambda_s = lambdaSampling
  mean_time_s = 1/lambda_s

  n <- sum(numEventArray)
  N <- sum(numSampleArray)

  poset <- random_Poset(n, graph_density=graphDensity)

  mean_times = runif(n, mean_time_s/lambdaSamplingScaling, lambdaSamplingScaling*mean_time_s)
  lambdas = 1/mean_times
  T_events <- matrix(0, N, n)
  T_sampling <- rexp(N, lambda_s)

  for (i in 1:n) {
    T_events[, i] <- rexp(N, lambdas[i])
  }
  T_sum_events <- matrix(0, N, n)
  topo_path = topological_Sort(poset)

  for (e in topo_path) {
    parents <- which(poset[, e] == 1)
    if (length(parents) == 0) {
      T_sum_events[, e] = T_events[, e]
    }else if (length(parents) == 1) {
      T_sum_events[, e] = T_events[, e] + T_sum_events[,parents]
    }else {
      T_sum_events[, e] = T_events[, e] + apply(T_sum_events[,parents], 1, max)
    }
  }
  obs_events <- hidden_genotypes <- matrix(0, N, n)

  for (i in 1:n) {
    index_i <- targetToIndex(numEventArray, i)

    indexes <- which(T_sum_events[, i] <= T_sampling)
    hidden_genotypes[indexes, i] <- 1
    for (j in 1:numSet) {
      range_temp <- indexToRange(numSampleArray, j)
      index_temp <- indexes[which(range_temp[1] <= indexes & indexes <= range_temp[2])]
      obs_events[index_temp, i] <- rbinom(length(index_temp), 1, 1-eps[j,index_i])
    }

    indexes <- which(T_sum_events[, i] > T_sampling)
    for (j in 1:numSet) {
      range_temp <- indexToRange(numSampleArray, j)
      index_temp <- indexes[which(range_temp[1] <= indexes & indexes <= range_temp[2])]
      obs_events[index_temp, i] <- rbinom(length(index_temp), 1, eps[j,index_i])
    }
  }

  eps_obs <- eps

  for (i in 1:numEventKind) {
    for (j in 1:numSet) {
      range_i <- indexToRange(numEventArray, i)
      range_j <- indexToRange(numSampleArray, j)

      eps_obs[j,i] <- sum(abs(hidden_genotypes[range_j[1]:range_j[2],range_i[1]:range_i[2]] -
                            obs_events[range_j[1]:range_j[2],range_i[1]:range_i[2]]))/(
                              (range_j[2] - range_j[1] + 1)*(range_i[2] - range_i[1] + 1))

    }
  }

  return(list(eps = eps, eps_obs = eps_obs, poset = poset, lambdas = lambdas, T_sampling = T_sampling, T_events = T_events,
              T_sum_events = T_sum_events, hidden_genotypes = hidden_genotypes, obs_events = obs_events))

}



#' @title Generate simulation time data
#' @export
#' @description Generate simulation data
#' @param numEventArray A array describe the event number in each event type
#' @param numSampleArray A array describe the set number in each set type
#' @param initialNum The initial number of independent sample
#' @param lambdaSampling The value of lambdas. Default: 1
#' @param sample_event_rate The sample event rate for independent
#' @param lambdaSamplingScaling The value of scaling which used to sample lambdas. Default: 3
#' @param baseEpsilon The value of baseEpsilon. Default: 0.05
#' @param epsilonSamplingRate The value of rate which used to sample epsilon. Default: 0.5
#' @param graphDensity The value which used to control the sparsity of poset. Default: 0.2
simulation_Time_Data_Generate <- function(numEventArray, numSampleArray, initialNum, sample_event_rate, lambdaSampling = 1, lambdaSamplingScaling = 3,
                                          baseEpsilon = 0.05, epsilonSamplingRate = 0.5, graphDensity = 0.2){

  if(initialNum >= sum(numSampleArray)){
    cat("initialNum >= ", sum(numSampleArray), sep = "\n")
    return(NULL)
  }

  numEventKind <- length(numEventArray)
  numSet <- length(numSampleArray)

  if(epsilonSamplingRate >= 1){
    cat("Error: epsilonSamplingRate >= 1\n")
    return(NULL)
  }

  eps_min <- baseEpsilon - baseEpsilon*epsilonSamplingRate
  eps_max <- baseEpsilon + baseEpsilon*epsilonSamplingRate
  eps <- c(eps_min, eps_min + (1:(numEventKind*numSet - 1))*(eps_max - eps_min)/(numEventKind*numSet - 1))
  eps = matrix(eps, ncol = numEventKind, byrow = T)

  lambda_s = lambdaSampling
  mean_time_s = 1/lambda_s

  n <- sum(numEventArray)
  N <- initialNum

  poset <- random_Poset(n, graph_density=graphDensity)
  mean_times = runif(n, mean_time_s/lambdaSamplingScaling, lambdaSamplingScaling*mean_time_s)
  lambdas = 1/mean_times
  T_events <- matrix(0, N, n)
  T_sampling <- rexp(N, lambda_s)
  for (i in 1:n) {
    T_events[, i] <- rexp(N, lambdas[i])
  }
  topo_path = topological_Sort(poset)
  N <- sum(numSampleArray) - N
  T_events <- rbind(T_events, matrix(0, N, n))
  T_sampling <- c(T_sampling, rexp(N, lambda_s))

  for (i in (initialNum + 1):(initialNum + N)) {
    index <- 0
    while (!(index > 0 & index <= n)) {
      index <- floor(min(sample(1:n,1)*sample_event_rate, n))
    }
    temp <- topo_path[1:index]
    T_events_temp <- T_events[sample(1:(i-1),1), ]
    if(index < n){
      temp <- topo_path[(index + 1):n]
      for (j in 1:length(temp)) {
        T_events_temp[temp[j]] <- rexp(1, lambdas[temp[j]])
      }
    }
    T_events[i,] <- T_events_temp
    T_sampling[i] <- rexp(1, lambda_s)
  }

  N <- dim(T_events)[1]

  T_sum_events <- rbind(matrix(0, N, n))
  for (e in topo_path) {
    parents <- which(poset[, e] == 1)
    if (length(parents) == 0) {
      T_sum_events[, e] = T_events[, e]
    }else if (length(parents) == 1) {
      T_sum_events[, e] = T_events[, e] + T_sum_events[,parents]
    }else {
      T_sum_events[, e] = T_events[, e] + apply(T_sum_events[,parents], 1, max)
    }
  }

  obs_events <- hidden_genotypes <- matrix(0, N, n)

  for (i in 1:n) {
    index_i <- targetToIndex(numEventArray, i)

    indexes <- which(T_sum_events[, i] <= T_sampling)
    hidden_genotypes[indexes, i] <- 1
    for (j in 1:numSet) {
      range_temp <- indexToRange(numSampleArray, j)
      index_temp <- indexes[which(range_temp[1] <= indexes & indexes <= range_temp[2])]
      obs_events[index_temp, i] <- rbinom(length(index_temp), 1, 1-eps[j,index_i])
    }

    indexes <- which(T_sum_events[, i] > T_sampling)
    for (j in 1:numSet) {
      range_temp <- indexToRange(numSampleArray, j)
      index_temp <- indexes[which(range_temp[1] <= indexes & indexes <= range_temp[2])]
      obs_events[index_temp, i] <- rbinom(length(index_temp), 1, eps[j,index_i])
    }
  }

  eps_obs <- eps

  for (i in 1:numEventKind) {
    for (j in 1:numSet) {
      range_i <- indexToRange(numEventArray, i)
      range_j <- indexToRange(numSampleArray, j)

      eps_obs[j,i] <- sum(abs(hidden_genotypes[range_j[1]:range_j[2],range_i[1]:range_i[2]] -
                                obs_events[range_j[1]:range_j[2],range_i[1]:range_i[2]]))/(
                                  (range_j[2] - range_j[1] + 1)*(range_i[2] - range_i[1] + 1))

    }
  }

  return(list(eps = eps, eps_obs = eps_obs, poset = poset, lambdas = lambdas, T_sampling = T_sampling, T_events = T_events,
              T_sum_events = T_sum_events, hidden_genotypes = hidden_genotypes, obs_events = obs_events))

}
