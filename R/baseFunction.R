## base function

#' @title Estimate epsilon
#' @export
#' @description Estimate epsilon for multi set and event
#' @param pat A matrix for observed genotypes. The col represents the event
#' index and the row represents a observation whose entries indicate whether
#' an event has been observed "1" or not "0". Data for different datasets
#' and event sets is represented using a block matrix
#'
#' @param isF A matrix indicates whether the data in dataset[i] and eventset[j]
#' has been filled "True" or not "False"
#'
#' @param isCE A matrix indicates whether calculate epslion of dataset[i] and
#' eventset[j],"True" or not "False"
#'
#' @param eps The initial matrix of epslion
#'
#' @param setD The description of dataset index
#'
#' @param eventD The description of eventset index
#'
#' @param multi_thrds The number of threads for epsilon estimation. Default: 1
#'
#' @param threshold The threshold of epslion estimate. Default: 0.05
#'
#' @param threshold1 The threshold of epslion estimate method selection. Default: 7
#'
#' @param T The temperature of simulated annealing algorithm. Default: 10.0
#'
#' @param N_iter The max iter number of simulated annealing algorithm. Default: 150
#'
#' @param lambdaS Rate of the sampling process. Default: 1.0
#'
#' @param thrds The number of threads for parallel execution. Default: 1
#'
#' @param R The loop number of guess Function
#'
#' @return A matrix of estimated epsilon
estimate_Epsilon_ForMulti <- function(pat, isF, isCE = NULL, eps =NULL, setD, eventD, multi_thrds = 1L,
                                      threshold = 0.05, threshold1 = 7, n_p = 0L, T = 10.0, N_iter = 150L, lambdaS = 1.0, thrds=1L,
                                      R = 1L){

  dataSetN = nrow(isF)
  eventSetN = ncol(isF)

  if (is.null(isCE)){
    isCE = isF
  }

  if (is.null(eps)){
    eps = matrix(as.numeric(0),nrow = dataSetN,ncol = eventSetN)
  }

  total <- dim(setD)[1] * dim(eventD)[1]
  eps_res <- matrix(as.numeric(0),nrow = dataSetN,ncol = eventSetN)

  temp1 <- list()

  if(multi_thrds <= 1){
    for (l in 1:total) {
      l_temp <- l - 1
      i <- floor(l_temp/dim(eventD)[1]) + 1
      j <- l_temp%%dim(eventD)[1] + 1

      pat_temp <- as.matrix(pat[c(setD[i,1]:setD[i,2])+1, c(1,(eventD[j,1]:eventD[j,2])+1) ,drop = FALSE])
      isF_temp <- as.matrix(isF[i,j,drop = FALSE])
      isCE_temp <- as.matrix(isCE[i,j,drop = FALSE])
      eps_temp <- as.matrix(eps[i,j,drop = FALSE])
      setD_temp <- as.matrix(setD[i,,drop = FALSE])
      setD_temp <- setD_temp - setD_temp[1,1]
      eventD_temp <- as.matrix(eventD[j,,drop = FALSE])
      eventD_temp <- eventD_temp - eventD_temp[1,1] + as.integer(1)

      temp <- estimate_Epsilon(pat = pat_temp, isF = isF_temp, isCE = isCE_temp, eps = eps_temp,
                               setD = setD_temp, eventD = eventD_temp, threshold = threshold,
                               threshold1 = threshold1, n_p = n_p,T = T, N_iter = N_iter, lambdaS = lambdaS,
                               thrds = thrds, R = R)

      temp1[[l]] <- list(temp = temp, i = i, j = j)
    }

  }else{
    temp1 <- mclapply(1:total, function(l){
      l_temp <- l - 1
      i <- floor(l_temp/dim(eventD)[1]) + 1
      j <- l_temp%%dim(eventD)[1] + 1

      pat_temp <- as.matrix(pat[c(setD[i,1]:setD[i,2])+1, c(1,(eventD[j,1]:eventD[j,2])+1) ,drop = FALSE])
      isF_temp <- as.matrix(isF[i,j,drop = FALSE])
      isCE_temp <- as.matrix(isCE[i,j,drop = FALSE])
      eps_temp <- as.matrix(eps[i,j,drop = FALSE])
      setD_temp <- as.matrix(setD[i,,drop = FALSE])
      setD_temp <- setD_temp - setD_temp[1,1]
      eventD_temp <- as.matrix(eventD[j,,drop = FALSE])
      eventD_temp <- eventD_temp - eventD_temp[1,1] + as.integer(1)

      temp <- estimate_Epsilon(pat = pat_temp, isF = isF_temp, isCE = isCE_temp, eps = eps_temp,
                               setD = setD_temp, eventD = eventD_temp, threshold = threshold,
                               threshold1 = threshold1, n_p = n_p,T = T, N_iter = N_iter, lambdaS = lambdaS,
                               thrds = thrds, R = R)

      return(list(temp = temp, i = i, j = j))

    }, mc.cores = multi_thrds, mc.cleanup = TRUE)

  }

  for (l in 1:length(temp1)) {
    eps_res[temp1[[l]]$i, temp1[[l]]$j] <- temp1[[l]]$temp
  }

  return(eps_res)

}




#' @title Find poset set fastly for multiple data and event sets
#' @export
#'
#' @description Fast poset set estimation for the hidden conjunctive Bayesian
#' network model (H-CBN), First determine the partial order relationship
#' between any two events, and then perform fine-tune
#'
#' @param pat A matrix for observed genotypes. The col represents the event
#' index and the row represents a observation whose entries indicate whether
#' an event has been observed "1" or not "0". Data for different datasets
#' and event sets is represented using a block matrix
#'
#' @param isF A matrix indicates whether the data in dataset[i] and eventset[j]
#' has been filled "True" or not "False"
#'
#' @param isCE A matrix indicates whether calculate epslion of dataset[i] and
#' eventset[j],"True" or not "False"
#'
#' @param eps The matrix of epslion
#'
#' @param setD The description of dataset index
#'
#' @param eventD The description of eventset index
#'
#' @param Fine_Tune_Num Fine tune number for more accurate poset, experience shows that
#' the best number is 2. Default: 2
#'
#' @param vote_size The number of vote that take a vote to find a consensus poset. Default: 1
#'
#' @param vote_threshold The threshold in the vote. Default: 0.5
#'
#' @param vote_thrds The number of threads in the vote. Default: 1
#'
#' @param thrds The number of threads for parallel execution. Default: 1
#'
#' @param threshold The threshold for determining partial order relationships, it can filter
#' out more false positives 1 with the larger value. experience shows that it is appropriate to
#' threshold <= 0.01, if you set threshold <= 0, the program will search the threshold with gmm,
#' we recommend setting it to a small value(such as 0.0001) when Fine_Tune_Num >= 2. Default: 0.0001
#'
#' @param threshold2 The threshold for determining partial order relationships considering the tolerance,
#' it is appropriate to make threshold2 >= threshold. Default: 0.01
#'
#' @param gmm_min_cluster The min cluster number of gmm. Default: 2
#'
#' @param gmm_convenience_threshold The value used for stop increasing cluster number, it can
#' increase the cluster number with the smaller value, experience shows that it is appropriate to
#' make gmm_convenience_threshold <= 0.005. Default: 0.001
#'
#' @param gmm_filter_threshold The value used for searching threshold after gmm, it can increase
#' the threshold with the larger value, experience shows that it is appropriate to
#' make gmm_convenience_threshold <= gmm_filter_threshold <= 0.0075. Default: 0.002
#'
#' @param n_p The vector contain the number of events during fine-tune. Default: [6, 6]
#'
#' @param is_update_eps whether update eps in fine tune, experience
#'  shows that it is appropriate to make is_update_eps = FALSE. Default: FALSE
#'
#' @param T The temperature of simulated annealing algorithm. Default: 10.0
#'
#' @param N_iter The max iter number of simulated annealing algorithm. Default: 150
#'
#' @param lambdaS Rate of the sampling process. Default: 1.0
#'
#' @param R The loop number of guess Function
#'
#' @return A matrix of poset
find_Poset_ForVote <- function(pat, isF, isCE = NULL, eps, setD, eventD, Fine_Tune_Num=2L, vote_size = 1L,
                              vote_threshold = 0.5, vote_thrds = 1L, thrds=1L, threshold = 0.0001, threshold2 = 0.01,
                              gmm_min_cluster = 2, gmm_convenience_threshold = 0.001, gmm_filter_threshold = 0.002,
                              n_p = 6L, is_update_eps = FALSE, T = 10.0, N_iter = 150L, lambdaS = 1.0, R = 1L) {

  if(vote_size <= 1){
    poset <- find_Poset(pat = pat, isF = isF, isCE = isCE, eps = eps, setD = setD, eventD = eventD,
                        Fine_Tune_Num = Fine_Tune_Num, thrds = thrds, threshold = threshold, threshold2 = threshold2,
                        gmm_min_cluster = gmm_min_cluster, gmm_convenience_threshold = gmm_convenience_threshold,
                        gmm_filter_threshold = gmm_filter_threshold, n_p = n_p, is_update_eps = is_update_eps,
                        T = T, N_iter = N_iter, lambdaS = lambdaS, R = R)
    return(poset)
  }else{
    temp <- list()
    if(vote_thrds <= 1){
      for (l in 1:vote_size) {
        poset <- find_Poset(pat = pat, isF = isF, isCE = isCE, eps = eps, setD = setD, eventD = eventD,
                            Fine_Tune_Num = Fine_Tune_Num, thrds = thrds, threshold = threshold, threshold2 = threshold2,
                            gmm_min_cluster = gmm_min_cluster, gmm_convenience_threshold = gmm_convenience_threshold,
                            gmm_filter_threshold = gmm_filter_threshold, n_p = n_p, is_update_eps = is_update_eps,
                            T = T, N_iter = N_iter, lambdaS = lambdaS, R = R)
        temp[[l]] <- poset
      }
    }else{
      temp <- mclapply(1:vote_size, function(l){
        poset <- find_Poset(pat = pat, isF = isF, isCE = isCE, eps = eps, setD = setD, eventD = eventD,
                            Fine_Tune_Num = Fine_Tune_Num, thrds = thrds, threshold = threshold, threshold2 = threshold2,
                            gmm_min_cluster = gmm_min_cluster, gmm_convenience_threshold = gmm_convenience_threshold,
                            gmm_filter_threshold = gmm_filter_threshold, n_p = n_p, is_update_eps = is_update_eps,
                            T = T, N_iter = N_iter, lambdaS = lambdaS, R = R)

        return(poset)
      }, mc.cores = vote_thrds, mc.cleanup = TRUE)
    }

    for (i in 1:vote_size) {
      if(i == 1){
        poset <- temp[[i]]
      }else{
        poset_temp <- temp[[i]]
        for (m in 1:dim(poset_temp)[1]) {
          for (n in 1:dim(poset_temp)[2]) {
            if(poset_temp[m,n] == 1){
              poset[m,n] <- poset[m,n] + 1
            }
          }
        }
      }
    }

    poset <- poset/vote_size
    poset[poset <= vote_threshold] <- 0
    poset[poset > vote_threshold] <- 1
    poset <- as.data.frame(poset)
    for (i in 1:dim(poset)[2]) {
      poset[,i] <- as.integer(poset[,i])
    }
    poset <- as.matrix(poset)

    return(poset)
  }

}




#' @title Estimate Lambda for multiple data and event sets
#' @export
#'
#' @description Lambda estimation for the hidden conjunctive Bayesian
#' network model (H-CBN), the estimation based on mccbn
#' (https://github.com/cbg-ethz/MC-CBN)
#'
#' @param pat A matrix for observed genotypes. The col represents the event
#' index and the row represents a observation whose entries indicate whether
#' an event has been observed "1" or not "0". Data for different datasets
#' and event sets is represented using a block matrix. pat must fill fully.
#'
#' @param poset An adjacency matrix represents partial order sets
#'
#' @param isF A matrix indicates whether the data in dataset[i] and eventset[j]
#' has been filled "True" or not "False"
#'
#' @param eps The matrix of epslion
#'
#' @param setD The description of dataset index
#'
#' @param eventD The description of eventset index
#'
#' @param Multi_size The number of Lambda estimation. Default: 1
#'
#' @param Multi_thrds The number of threads in the Lambda estimation. Default: 1
#'
#' @param lambdaS Rate of the sampling process. Default: 1.0
#'
#' @param L Number of samples to be drawn from the proposal in the E-step.
#' Default: 100
#'
#' @param sampling sampling scheme to generate hidden genotypes, OPTIONS:
#' "forward" - generate occurrence times according to current rate
#'             parameters, and, from them, generate the hidden genotypes
#' "add-remove" - generate genotypes, from observed genotypes, using a two-steps
#'                proposal. First, pick a move uniformly at random: either to
#'                add or to remove an event. Events are chosen to be removed
#'                with probability proportional to their rates, and to be added
#'                with an inverse probability. Second, make genotypes compatible
#'                with the poset by either adding or removing all events
#'                incompatible with the poset
#' "backward" - enumerate all genotypes with Hamming distance
#' "bernoulli" - generate genotypes from a Bernoulli distribution with success
#'               probability p = epsilon
#' "pool" - generate a pool of compatible genotypes according to current rate
#'          parameters and sample K observations proportional to their Hamming
#'          distance
#'
#' @param maxIter the maximum number of EM iterations. Default: 100
#'
#' @param updateStepSize number of EM steps after which the number of samples,
#'  L is doubled. L is increased, if the difference in the parameter estimates
#'  between such consecutive batches is greater than the tolerance level, tol.
#'  Default: 20
#'
#' @param tol convergence tolerance for the error rate and the rate parameters.
#' The EM runs until the difference between the average estimates in the last
#' two batches is smaller than tol, or until max.iter is reached. Default: 0.001
#'
#' @param maxLambda an optional upper bound on the value of the rate
#' parameters. Default: 1e6
#'
#' @param neighborhood.dist an integer value indicating the Hamming distance
#' between the observation and the samples generated by "backward" sampling.
#' This option is used if sampling is set to "backward". Default: 1
#'
#' @param thrds number of threads for parallel execution. Default: 1
#'
#' @param is_update_eps whether update eps in lambda estimate program, experience
#'  shows that it is appropriate to make is_update_eps = FALSE. Default: FALSE
#'
#' @param seed seed for reproducibility
estimate_Lambda_ForMulti <- function(pat, poset, isF, eps, setD, eventD, Multi_size = 1L, Multi_thrds = 1L,
                                     thrds=1L, lambdaS = 1.0, L = 100L, sampling = 'add-remove',
                                     maxIter=100L, updateStepSize=20L, tol=0.001, maxLambda=1e6,
                                     neighborhoodDist=1L, is_update_eps = FALSE, seed=NULL) {
  if(Multi_size <= 1){
    fit <- estimate_Lambda(pat = pat, poset = poset, isF = isF, eps = eps, setD = setD, eventD = eventD,
                           thrds = thrds, lambdaS = lambdaS, L = L, sampling = sampling,
                           maxIter=maxIter, updateStepSize=updateStepSize, tol=tol, maxLambda=maxLambda,
                           neighborhoodDist=neighborhoodDist, is_update_eps = is_update_eps, seed=seed)
    return(fit)
  }else{
    temp <- list()
    if(Multi_thrds <= 1){
      for (l in 1:Multi_size) {
        fit <- estimate_Lambda(pat = pat, poset = poset, isF = isF, eps = eps, setD = setD, eventD = eventD,
                               thrds = thrds, lambdaS = lambdaS, L = L, sampling = sampling,
                               maxIter=maxIter, updateStepSize=updateStepSize, tol=tol, maxLambda=maxLambda,
                               neighborhoodDist=neighborhoodDist, is_update_eps = is_update_eps, seed=seed)
        temp[[l]] <- fit
      }
    }else{
      temp <- mclapply(1:Multi_size, function(l){
        fit <- estimate_Lambda(pat = pat, poset = poset, isF = isF, eps = eps, setD = setD, eventD = eventD,
                               thrds = thrds, lambdaS = lambdaS, L = L, sampling = sampling,
                               maxIter=maxIter, updateStepSize=updateStepSize, tol=tol, maxLambda=maxLambda,
                               neighborhoodDist=neighborhoodDist, is_update_eps = is_update_eps, seed=seed)

        return(fit)
      }, mc.cores = Multi_thrds, mc.cleanup = TRUE)
    }

    for (i in 1:Multi_size) {
      if(i == 1){
        fit <- temp[[i]]
        epsilon1_mean <- fit$epsilon1
        epsilon2_mean <- fit$epsilon2
        lambda_mean <- fit$lambda
      }else{
        fit <- temp[[i]]
        epsilon1_mean <- epsilon1_mean + fit$epsilon1
        epsilon2_mean <- epsilon2_mean + fit$epsilon2
        lambda_mean <- lambda_mean + fit$lambda

      }
    }

    return(list(epsilon1 = epsilon1_mean/Multi_size, epsilon2 = epsilon2_mean/Multi_size, lambda = lambda_mean/Multi_size))
  }

}





#' @title Target to index
#' @export
#' @description Find the max index which sum(arr|index|) >= target
#' @param arr A array
#' @param target A number
#' @returns A number
targetToIndex <- function(arr, target){

  if(length(arr) == 1){
    return(1)
  }else{
    for (i in 2:length(arr)) {
      arr[i] <- arr[i] + arr[i-1]
    }
    index <- min(which(arr >= target))

    return(index)
  }

}


#' @title Index To Target
#' @export
#' @description Find the range condition on arr and index
#' @param arr A array
#' @param index A number
#' @returns A range array
indexToRange <- function(arr, index){

  if(index == 1){
    return(c(1,arr[1]))
  }else if(index <= length(arr)){
    for (i in 2:index) {
      arr[i] <- arr[i] + arr[i-1]
    }
    return(c(arr[index-1]+1, arr[index]))
  }

  return(NULL)
}


#' @title Topological Sort
#' @export
#' @description Perform topological sort with poset
#' @param poset A matrix for poset
#' @returns A sorted list
topological_Sort <- function(poset) {
  p = ncol(poset)
  sorted.list = rep(0, p)
  for(i in 1:p) {
    sorted.list[i] = min.degree = which.min(apply(poset, 2, sum))
    poset[min.degree, ] = 0
    poset[, min.degree] = Inf
  }

  sorted.list
}


#' @title Random Poset
#' @export
#' @description Generate a random poset
#' @param p The num of event
#' @param graph_density The density of graph
#' @param trans_reduced Add the transitive reduced poset to the list
#' @returns A poset
random_Poset <- function(p, graph_density=0.15, trans_reduced = TRUE){
  random_Pposets(nr_pos=1, nr_muts=p , ldenses=graph_density, trans_reduced = trans_reduced)[[1]]
}

trans_reduction <- function(A) {

  old_names = dimnames(A)

  colnames(A) =  rownames(A) = paste("n", 1:nrow(A), sep='')
  R <- as.relation(A)

  RT = transitive_reduction(R)
  res = relation_incidence(RT)
  res = res[rownames(A),colnames(A)]

  dimnames(res) = old_names
  res
}

random_Pposets <- function(nr_pos, nr_muts , ldenses, trans_reduced = TRUE){

  if(length(ldenses) == 1 ) {
    ldenses = rep(ldenses, nr_pos)
  } else if(length(ldenses) < nr_pos ) {
    stop( paste("Invalid value for ldenses ", ldenses) )
  }

  if(length(nr_muts) == 1 ) {
    nr_muts = rep(nr_muts, nr_pos)
  } else if(length(nr_muts) < nr_pos ) {
    stop( paste("Invalid value for nr_muts ", nr_edges) )
  }

  nr_edges = ceiling(choose(nr_muts, 2)*ldenses)

  # Initialisation
  all_posets = rep(0, nr_pos)
  pos = 1

  while (pos <= nr_pos){

    poset <- matrix(0, nr_muts[pos], nr_muts[pos])
    poset[upper.tri(poset)][sample(choose(nr_muts[pos], 2), nr_edges[pos])] <- 1

    # added by Hesam. Add the transitive reduced poset to the list
    if(trans_reduced) {
      poset = trans_reduction(poset)
    }

    all_posets[pos] = list(poset)
    pos = pos + 1
  }
  all_posets
}

#' @title Rate Time Help
#' @export
#' @description Convert the lambda to the time bound based on lambdaS
#' @param lambda_s The sampling lambda in CBN
#' @param lambdas The lambda of each event in CBN
#' @param poset The poset matrix in CBN
#' @returns A list that contain the time's low and up bound of each event
rateTimeHelp <- function(lambda_s, lambdas, poset){
  T_events <- lambda_s/lambdas
  T_sum_events <- rep(0, length(lambdas))
  T_sum_events1 <- rep(0, length(lambdas))

  topo_path = topological_Sort(poset)
  for (e in topo_path) {
    parents <- which(poset[, e] == 1)
    if (length(parents) == 0) {
      T_sum_events[e] = T_events[e]
    }else {
      T_sum_events[e] = T_events[e] + max(T_sum_events[parents])
    }
  }

  for (e in topo_path) {
    childrens <- which(poset[e, ] == 1)
    if (length(childrens) == 0) {
      T_sum_events1[e] = Inf
    }else {
      T_sum_events1[e] = min(T_sum_events[childrens])
    }
  }

  return(list(LB=T_sum_events, UB=T_sum_events1))
}




Genotype_Generate <- function(sampling_N_Total, max_One_N, poset, lambda_s = 1, lambdas, cores = 1){
  topo_path = topological_Sort(poset)
  n <- length(lambdas)
  loop <- ceiling(sampling_N_Total/max_One_N)
  for (l in 1:loop) {
    if(l == loop){
      sampling_N <- sampling_N_Total - max_One_N*(loop-1)
    }else{
      sampling_N <- max_One_N
    }

    temp <- mclapply(1:cores, function(i){
      if(i == cores){
        N <- sampling_N - (cores-1)*ceiling(sampling_N/cores)
      }else{
        N <- ceiling(sampling_N/cores)
      }

      T_sampling <- rexp(N, lambda_s)
      T_events <- matrix(0, N, n)
      for (i in 1:n) {
        T_events[, i] <- rexp(N, lambdas[i])
      }
      T_sum_events <- matrix(0, N, n)

      for (e in topo_path) {
        parents <- which(poset[, e] == 1)
        if (length(parents) == 0) {
          T_sum_events[, e] = T_events[, e]
        }
        else if (length(parents) == 1) {
          T_sum_events[, e] = T_events[, e] + T_sum_events[,parents]
        }
        else {
          T_sum_events[, e] = T_events[, e] + apply(T_sum_events[,parents], 1, max)
        }
      }
      hidden_genotypes = matrix(0, N, n)
      for (i in 1:n) {
        indexes =  which(T_sum_events[, i] <= T_sampling)
        hidden_genotypes[indexes, i] = 1
      }
      return(unique(hidden_genotypes))
    }, mc.cores = cores, mc.cleanup = TRUE)

    ans <- temp[[1]]
    for (i in 2:length(temp)) {
      ans <- rbind(ans, temp[[i]])
    }
    ans <- unique(ans)

    if (l == 1){
      res <- ans
    }else{
      res <- unique(rbind(res, ans))
    }
  }

  return(res)
}



Donor_Pair_Geno_Gerate <- function(time_Donor_Geno, Genotype_Sampling, eps, genotype_len, N, cores = 1){
  eps1 <- 1-eps

  temp <- mclapply(1:dim(time_Donor_Geno)[1], function(i){
    time_Donor_Geno_Temp <- time_Donor_Geno[i,]

    g_e_matrix <- abs(sweep(Genotype_Sampling, 2, time_Donor_Geno_Temp, "-"))
    g_t_matrix <- 1 - g_e_matrix
    log_error_value_vector <- rowSums(log(g_t_matrix * eps1 + g_e_matrix * eps))

    sort_Temp <- order(log_error_value_vector, decreasing = T)[1:N]

    # for (j in 1:dim(Genotype_Sampling_Temp)[1]) {
    #   g_e <- abs(Genotype_Sampling_Temp[j, 1:genotype_len] - time_Donor_Geno[i,])
    #   g_t <- (1-g_e)
    #   Genotype_Sampling_Temp[j, genotype_len + 1] <- sum(log(g_t*eps1 + g_e*eps))
    # }

    return(cbind(Genotype_Sampling[sort_Temp,], log_error_value = log_error_value_vector[sort_Temp]))

  }, mc.cores = cores, mc.cleanup = TRUE)

  for (i in 1:length(temp)) {
    names(temp)[[i]] <- rownames(time_Donor_Geno)[i]
  }
  return(temp)
}


#' @title Find most compatible genotype by sampling
#' @export
#' @description Find most compatible genotype by sampling, first create compatible genotype by Genotype_Generate,
#' then find top k compatible genotype by Donor_Pair_Geno_Gerate
#' @param sampling_loop The number of outer sampling loop
#' @param sampling_Num_Total The number of sampling compatible genotype in each outer sampling loop,
#' the total number = sampling_loop * sampling_Num_Total
#' @param max_One_Num The max sampling number in each innerloop
#' @param poset The poset matrix in CBN
#' @param lambda_s The sampling lambda in CBN
#' @param lambdas The lambda of each event in CBN
#' @param drivers The drivers
#' @param time_Donor_Geno The true genotype
#' @param eps The epsilon in CBN
#' @param max_Select_Num The max number of compatible genotype of each outer sampling loop
#' @param cache_dir The cache dir used for save temp files
#' @param cores1 The parallel cores in Genotype_Generate
#' @param cores2 The parallel cores in Donor_Pair_Geno_Gerate
#' @returns The list of donor pair compatible genotype
find_most_Compatible_Genotype_by_Sampling <- function(sampling_loop, sampling_Num_Total, max_One_Num, poset, lambda_s, lambdas,
                                                      drivers, time_Donor_Geno, eps,
                                                      max_Select_Num, cache_dir, cores1 = 1, cores2 = 1){
  cache_dir <- paste(cache_dir, "/", sep = "")
  cache_dir <- gsub("//","/",cache_dir)

  genotype_len <- length(eps)

  for (s_loop in 1:sampling_loop) {
    # cat("\nGenotype_Sampling  ", s_loop)
    if(!file.exists(paste(cache_dir, "Genotype_Sampling_", s_loop, ".rds", sep = ""))){
      Genotype_Sampling <- Genotype_Generate(sampling_Num_Total, max_One_Num, poset, lambda_s, lambdas, cores1)
      # Genotype_Sampling <- cbind(Genotype_Sampling,rep(NA, dim(Genotype_Sampling)[1]))
      # colnames(Genotype_Sampling) <- c(driver_genes, driver_chr, "log_error_value")
      colnames(Genotype_Sampling) <- c(driver_genes, driver_chr)
      rownames(Genotype_Sampling) <- paste("sampling_Geno_", s_loop, "_", 1:dim(Genotype_Sampling)[1],sep = "")

      saveRDS(Genotype_Sampling, paste(cache_dir, "Genotype_Sampling_", s_loop, ".rds", sep = ""))
    }else{
      Genotype_Sampling <- readRDS(paste(cache_dir, "Genotype_Sampling_", s_loop, ".rds", sep = ""))
    }

    # cat("   Donor_Pair_Geno  ", s_loop, "\n")
    temp <- Donor_Pair_Geno_Gerate(time_Donor_Geno, Genotype_Sampling, eps, genotype_len, max_Select_Num, cores2)
    if (s_loop == 1){
      if(file.exists(paste(cache_dir,"Donor_Pair_Geno_", s_loop, "_temp.rds", sep = ""))){
        Donor_Pair_Geno <- readRDS(paste(cache_dir,"Donor_Pair_Geno_", s_loop, "_temp.rds", sep = ""))
      }else{
        for (ll in 1:length(temp)) {
          Donor_Pair_Geno[[ll]] <- temp[[ll]]
        }
      }
    }else{
      for (ll in 1:length(Donor_Pair_Geno)) {
        Donor_Pair_Geno[[ll]] <- rbind(Donor_Pair_Geno[[ll]], temp[[ll]])
      }
    }

    saveRDS(Donor_Pair_Geno, paste(cache_dir,"Donor_Pair_Geno_", s_loop, "_temp.rds", sep = "") )
  }

  return(Donor_Pair_Geno)
}

#' @title Find most compatible genotype by flipping
#' @export
#' @description Find most compatible genotype by flipping compatible event, it is faster than Sampling
#' @param genotype The true genotype that from the same dataset
#' @param eps The epsilon of each event
#' @param poset The poset matrix in CBN
#' @param max_iter The max loop iter for each sample
#' @param cores The parallel cores
#' @returns The matrix of donor pair compatible genotype
find_most_Compatible_Genotype_by_Flipping <- function(genotype, eps, poset, max_iter = 10000, cores = 1){

  if(dim(poset)[1]!=length(eps) | dim(poset)[2]!=length(eps)){
    cat("Error: Please check length(eps) and dim(poset)!\n")
  }

  parent_Child_Set <- list()
  topo_path = as.integer(topological_Sort(poset))
  for (e in topo_path) {
    parent_Child_Set[[e]] <- list(par = c(), chi = c())
    parent_Child_Set[[e]]$par <- which(poset[, e] == 1)
    parent <- parent_Child_Set[[e]]$par
    if(length(parent) > 0){
      for (i in parent) {
        if(length(parent_Child_Set[[i]]$par) > 0 ){
          parent_Child_Set[[e]]$par <- unique(c(parent_Child_Set[[e]]$par, parent_Child_Set[[i]]$par))
        }
      }
    }
  }

  for (e in topo_path) {
    child <- setdiff(1:length(topo_path), e)
    if(length(parent_Child_Set[[e]]$par) > 0){
      child <- setdiff(child, parent_Child_Set[[e]]$par)
    }
    if(length(child) > 0){
      for (i in child) {
        if(length(parent_Child_Set[[i]]$par)>0 & e%in%parent_Child_Set[[i]]$par ){
          parent_Child_Set[[e]]$chi <- unique(c(parent_Child_Set[[e]]$chi, i))
        }
      }
    }
  }

  n <- length(eps)

  topo_path1 <- topo_path
  for (i in 1:n) {
    topo_path1[i] <- which(topo_path == i)
  }

  temp <- mclapply(1:dim(genotype)[1], function(i){
    is_Deal <- TRUE
    index <- 1

    while (is_Deal) {
      if(index == 1){
        geno_T <- genotype[i,]
      }else if(index > max_iter){
        geno_T <- NULL
        break
      }else{
        geno_T[diff_T[1,6]] <- 1-geno_T[diff_T[1,6]]
      }

      diff_T <- array(0, dim = c(n,6))
      for (j in 1:n) {
        if(length(parent_Child_Set[[j]]$par)>0){
          diff_T[j,1] <- length(which(geno_T[parent_Child_Set[[j]]$par] != 1))
        }
        if(length(parent_Child_Set[[j]]$chi)>0){
          diff_T[j,2] <- length(which(geno_T[parent_Child_Set[[j]]$chi] != 0))
        }
        if(geno_T[j] == 1){
          diff_T[j,3] <- diff_T[j,1] - diff_T[j,2]
        }else{
          diff_T[j,3] <- diff_T[j,2] - diff_T[j,1]
        }
      }

      diff_T[,4] <- eps
      diff_T[,5] <- topo_path1
      diff_T[,6] <- 1:n
      diff_T <- diff_T[order(diff_T[,5], decreasing = T),]
      diff_T <- diff_T[order(diff_T[,4], decreasing = T),]
      diff_T <- diff_T[order(diff_T[,3], decreasing = T),]

      if(diff_T[1,3] > 0){
        is_Deal <- TRUE
      }else{
        is_Deal <- FALSE
      }
      index <- index + 1

    }

    return(geno_T)
  }, mc.cores = cores, mc.cleanup = TRUE)

  genotype <- cbind(genotype, rep(NA, dim(genotype)[1]))

  for (i in 1:length(temp)) {
    if(is.null(temp[[i]])){
      genotype[i,n+1] <- FALSE
    }else{
      genotype[i,n+1] <- TRUE
      genotype[i,1:n] <- temp[[i]]
    }
  }
  genotype <- genotype[which(genotype[,n+1] == TRUE),1:n]
  return(genotype)
}



#' @title Donor Pair Genotype Filter
#' @export
#' @description Filter donor pair genotype that used for time estimate
#' @param Donor_Pair_Genotype The donor pair genotype
#' @param cancer_time The cancer time
#' @param method The method used to find most compatible genotype. Default: flipping
#' @param rate_Event_Save The rate used to filter genotype. Default: 0.01
#' @returns A filter donor pair genotype list that contain Genotype, is_Have_mustE, is_Used_Time_Estimate
donor_Pair_Genotype_Filter <- function(Donor_Pair_Genotype, cancer_time, method = "flipping", rate_Event_Save = 0.01){

  filter_Donor_Pair_Genotype <- list()

  if(method == "sampling"){
    genotype_len <- dim(Donor_Pair_Genotype[[1]])[2]-1
    for (i in 1:length(Donor_Pair_Genotype)) {
      donor <- names(Donor_Pair_Genotype)[i]

      Geno_Temp <- Donor_Pair_Genotype[[i]]
      Geno_Temp <- Geno_Temp[rownames(unique(Geno_Temp[,1:genotype_len])),]

      colnames(Geno_Temp)[dim(Geno_Temp)[2]] <- "Normalized_P"
      log_values <- Geno_Temp[,"Normalized_P"]
      max_log_value <- max(log_values)
      adjusted_log_values <- log_values - max_log_value
      exp_values <- exp(adjusted_log_values)
      Geno_Temp[,"Normalized_P"] <- exp_values / sum(exp_values)
      Geno_Temp <- Geno_Temp[order(Geno_Temp[,"Normalized_P"], decreasing = T),]
      Geno_Temp <- Geno_Temp[which(Geno_Temp[,"Normalized_P"] >=
                                     max(Geno_Temp[,"Normalized_P"]*rate_Event_Save)),]
      if(class(Geno_Temp)[1] == "numeric"){
        Geno_Temp <- matrix(Geno_Temp, nrow = 1)
        colnames(Geno_Temp) <- c(colnames(Donor_Pair_Genotype[[i]])[-1], "Normalized_P")
      }

      Geno_Temp[,"Normalized_P"] <- Geno_Temp[,"Normalized_P"]/sum(Geno_Temp[,"Normalized_P"])


      if (donor %in% rownames(cancer_time)){
        mustE <- colnames(cancer_time)[which(!is.na(cancer_time[donor,]))]
        mustE <- setdiff(mustE, c("age","healthTime", "sampleTime"))

        is_mustE <- c()

        for (j in 1:dim(Geno_Temp)[1]) {
          if(sum(Geno_Temp[j, mustE]) > 0){
            is_mustE <- c(is_mustE, TRUE)
          }else{
            is_mustE <- c(is_mustE, FALSE)
          }
        }

        if (any(is_mustE) == TRUE){
          filter_Donor_Pair_Genotype[[i]] <- list(Genotype = Geno_Temp, is_Have_mustE = is_mustE,
                                                  is_Used_Time_Estimate = TRUE)
        }else{
          filter_Donor_Pair_Genotype[[i]] <- list(Genotype = Geno_Temp, is_Used_Time_Estimate = FALSE)
        }
      }else{
        filter_Donor_Pair_Genotype[[i]] <- list(Genotype = Geno_Temp, is_Used_Time_Estimate = FALSE)
      }
      names(filter_Donor_Pair_Genotype)[i] <- donor
    }

    return(filter_Donor_Pair_Genotype)

  }else if(method == "flipping"){

    for (i in 1:dim(Donor_Pair_Genotype)[1]) {
      donor <- rownames(Donor_Pair_Genotype)[i]

      Geno_Temp <- Donor_Pair_Genotype[i,]

      if (donor %in% rownames(cancer_time)){
        mustE <- colnames(cancer_time)[which(!is.na(cancer_time[donor,]))]
        mustE <- setdiff(mustE, c("age","healthTime", "sampleTime"))

        is_mustE <- FALSE
        if(sum(Geno_Temp[mustE]) > 0){
          is_mustE <- TRUE
        }

        if (any(is_mustE) == TRUE){
          filter_Donor_Pair_Genotype[[i]] <- list(Genotype = Geno_Temp, is_Have_mustE = is_mustE,
                                                  is_Used_Time_Estimate = TRUE)
        }else{
          filter_Donor_Pair_Genotype[[i]] <- list(Genotype = Geno_Temp, is_Used_Time_Estimate = FALSE)
        }

      }else{
        filter_Donor_Pair_Genotype[[i]] <- list(Genotype = Geno_Temp, is_Used_Time_Estimate = FALSE)
      }
      names(filter_Donor_Pair_Genotype)[i] <- donor
    }

    return(filter_Donor_Pair_Genotype)


  }else{
    return(NULL)
  }

}



#' @title parent Set Help
#' @export
#' @description find parent set for each event
#' @param poset The poset matrix in CBN
#' @returns A list contain the parent set of each event
parent_Set_Help <- function(poset){
  parent_Set <- list()
  topo_path = topological_Sort(poset)
  index_C <- 1
  for (e in topo_path) {
    parent_Set[[index_C]] <- which(poset[, e] == 1)
    parent <- parent_Set[[index_C]]
    for (i in parent) {
      if(!is.null(parent_Set[[paste(i, "_" ,sep = "")]])){
        parent_Set[[index_C]] <- c(parent_Set[[index_C]], parent_Set[[paste(i, "_" ,sep = "")]] )
      }
    }
    parent_Set[[index_C]] <- unique(parent_Set[[index_C]])
    names(parent_Set)[index_C] <-  paste(e, "_" ,sep = "")
    index_C <- index_C + 1
  }
  return(parent_Set)
}



#' @title parent Set Help1
#' @export
#' @description find parent set for each event
#' @param poset The poset matrix in CBN
#' @returns A list contain the parent set of each event
parent_Set_Help1 <- function(poset){
  topo_path <- as.integer(topological_Sort(poset))
  child_Set <- list()
  for (i in length(topo_path):1) {
    e <- topo_path[i]
    child_Set[[e]] <- list(child = c())
    child_Set[[e]]$child <- which(poset[e,] == 1)
    child <- child_Set[[e]]$child
    if(length(child) > 0){
      for (j in child) {
        if( length(child_Set[[j]]$child) > 0){
          child_Set[[e]]$child <- unique(c(child_Set[[e]]$child, child_Set[[j]]$child))
        }
      }
    }
  }
  return(child_Set)
}



#' @title Gamma Cluster
#' @export
#' @description Gamma Cluster of 1 dim data based on GMM
#' @param data The vector of input data
#' @param rate_cluster The rate of cluster
#' @param max_iter The max update iters
#' @param tol The tolerance
#' @returns A list that contain the estimation results of weights, shape params, rate params, responsibilities, cluster
GammaCluster <- function(data,rate_cluster,max_iter = 10000,tol = 1e-8){

  data_temp <- sort(data)
  n <- length(data)
  k <- length(rate_cluster)
  weights <- rate_cluster/sum(rate_cluster)


  if(k > 1){
    shape_params <- c()
    rate_params <- c()
    boundaries <- cumsum(c(0, weights)) * n
    for (i in 1:length(weights)) {
      start <- ceiling(boundaries[i]) + 1
      end <- floor(boundaries[i + 1])
      gamma_fit <- fitdist(data_temp[start:end],"gamma",method = "mle")
      shape_params <- c(shape_params, gamma_fit$estimate[1])
      rate_params <- c(rate_params, gamma_fit$estimate[2])
    }
    e_step <- function(data, weights, shape_params, rate_params) {
      responsibilities <- matrix(0, n, k)
      for (j in 1:k) {
        temp <- dgamma(data, shape = shape_params[j], rate = rate_params[j])
        temp[temp == 0] <- 1e-10
        responsibilities[, j] <- weights[j] * temp
      }
      responsibilities <- responsibilities / rowSums(responsibilities)
      return(responsibilities)
    }
    m_step <- function(data, responsibilities) {
      weights <- colMeans(responsibilities)
      shape_params <- numeric(k)
      rate_params <- numeric(k)

      for (j in 1:k) {
        weighted_data <- responsibilities[, j] * data
        mean_data <- sum(weighted_data) / sum(responsibilities[, j])
        var_data <- sum(responsibilities[, j] * (data - mean_data)^2) / sum(responsibilities[, j])

        if(var_data == 0){
          return(NULL)
        }

        shape_params[j] <- mean_data^2 / var_data
        rate_params[j] <- mean_data / var_data
      }

      return(list(weights = weights, shape_params = shape_params, rate_params = rate_params))
    }
    log_likelihoods <- numeric(max_iter)
    for (iter in 1:max_iter) {
      # E Step
      responsibilities <- e_step(data, weights, shape_params, rate_params)

      # M Step
      params <- m_step(data, responsibilities)
      if(is.null(params)){
        cluster <- apply(responsibilities, 1, which.max)
        return(list(weights = weights, shape_params = shape_params, rate_params = rate_params, responsibilities = responsibilities,
                    cluster = cluster)  )
      }
      sumshape <- sum((params$shape_params - shape_params)^2)
      sumrate <- sum((params$rate_params - rate_params)^2)
      if((sumshape <= tol) & (sumrate <= tol)){
        break
      }

      weights <- params$weights
      shape_params <- params$shape_params
      rate_params <- params$rate_params

      # log_likelihood
      log_likelihood <- sum(log(rowSums(sapply(1:k, function(j) {
        weights[j] * dgamma(data, shape = shape_params[j], rate = rate_params[j])
      }))))

      log_likelihoods[iter] <- log_likelihood
    }

    cluster <- apply(responsibilities, 1, which.max)

    return(list(weights = weights, shape_params = shape_params, rate_params = rate_params, responsibilities = responsibilities,
                cluster = cluster)  )

  }else{
    return(list())
  }

}




#' @title genotype SolveQP
#' @export
#' @description find the best z1 and z2 by SolveQP based on genotype
#' @param posetCluster The list of different subpathway
#' @param age The age of the donor
#' @param donor The id of the donor
#' @param cancer_time The cancer time
#' @param Geno The compatible genotype
#' @param driver_genes_len The length(driver_genes)
#' @param filter_num The number matrix of chr time event
#' @param rateT The list of low and up bound rate of each event
#' @param isWeight If isWeight = true, use filter num to weighted average
#' @returns A array of the best z1 and z2
genotype_SolveQP <- function(posetCluster, age, donor, cancer_time, Geno, driver_genes_len, filter_num, rateT, posetClusterWeight,
                             isWeight = TRUE){

  ti <- cancer_time[donor,1:(dim(cancer_time)[2]-3)]
  indexS <- intersect(which(!is.na(ti)), which(Geno == 1)-driver_genes_len)
  ti <- ti[indexS]
  num <- filter_num[donor,1:(dim(filter_num)[2])][indexS]
  indexS <- indexS + driver_genes_len
  ri <- rateT$LB[indexS]

  Amat_L <- matrix(0,nrow = 2, ncol = 2)
  Amat_L[1,1] <- -1
  Amat_L[2,1] <- 1
  bvec_L <- c(-age,age)
  Amat_L_array <- rep(0, length(posetCluster))
  Weight1 <- posetClusterWeight
  Weight2 <- 0

  for (j in 1:length(posetCluster)) {
    Amat_L_array[j] <- -rateT$LB[posetCluster[[j]][1]]
    if(!is.infinite(rateT$UB[posetCluster[[j]][1]])){
      Amat_L[2,2] <- Amat_L[2,2] + rateT$UB[posetCluster[[j]][1]]*posetClusterWeight[j]
      Weight2 <- Weight2 + posetClusterWeight[j]
    }
  }
  if(Weight2 == 0){
    Amat_L[2,2] <- Inf
    Amat_L[1,2] <- sum(Amat_L_array*Weight1)/sum(Weight1)
  }else{
    Amat_L[2,2] <- Amat_L[2,2]/Weight2
    while ( abs(sum(Amat_L_array*Weight1)/sum(Weight1)) >=  Amat_L[2,2] ) {
      delta <- rep(0, length(posetCluster))
      for (j in 1:length(posetCluster)) {
        Weight2 <- Weight1
        Weight2[j] <- 0
        delta[j] <- abs(sum(Amat_L_array*Weight1)/sum(Weight1)) - abs(sum(Amat_L_array*Weight2)/sum(Weight2))
      }
      Weight1[which(delta == max(delta))] <- 0
    }
    Amat_L[1,2] <- sum(Amat_L_array*Weight1)/sum(Weight1)
  }


  Amat <- rbind(Amat_L, matrix(c(1,0,0,1),nrow = 2, byrow = T))
  indexS <- which(!(Amat[,2] == Inf | Amat[,2] == -Inf))
  Amat <- Amat[indexS,]
  Amat <- t(Amat)
  bvec <- c(bvec_L, 0, 0)
  bvec <- bvec[indexS]

  Dmat <- matrix(NA, nrow = 2, ncol = 2)
  Dmat[1,1] <- 1
  if (isWeight){
    Dmat[1,2] <- Dmat[2,1] <- sum(ri*num/sum(num))
    Dmat[2,2] <- sum(ri*ri*num/sum(num))
    dvec <- c(sum(ti*num/sum(num)),sum(ri*ti*num/sum(num)))*2*age
  }else{
    Dmat[1,2] <- Dmat[2,1] <- mean(ri)
    Dmat[2,2] <- mean(ri*ri)
    dvec <- c(mean(ti),mean(ri*ti))*2*age
  }

  Dmat <- Dmat*2

  if(!is.positive.definite(Dmat)){
    Dmat <- as.matrix(nearPD(Dmat)$mat)
  }

  sol <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec)
  sol_Temp <- sol$solution[1:2]
  sol_Temp[which(sol_Temp <= 0)] <- 0.0000001

  errors_MSE <- (0.5*sol_Temp%*%Dmat - dvec)%*%sol_Temp
  Amat <- t(Amat)
  Dmat <- Dmat*0.5


  return(list(sol_last = sol_Temp, Amat = Amat, bvec = bvec,
              Dmat = Dmat, dvec = dvec, errors_MSE = errors_MSE))

}




#' @title MH Pretreatment based on flipping
#' @export
#' @description MH Pretreatment based on flipping
#' @param filter_Donor_Pair_Genotype The filter donor pair genotype
#' @param time_temp The time temp
#' @param cancer_time The cancer time
#' @param filter_num The number matrix of chr time event
#' @param rateT The list of low and up bound rate of each event
#' @param driver_genes The driver genes
#' @param isWeight If isWeight = true, use filter num to weighted average
#' @param cBioCli The cBioCli
#' @param topo_path The topo path
#' @param parent_Set The parent set
#' @returns A array of the best z1 and z2
MH_Pretreatment_based_on_flipping <- function(filter_Donor_Pair_Genotype, time_temp, cancer_time, filter_num, rateT, driver_genes, isWeight,
                            cBioCli, topo_path, parent_Set, gammaClusterRes){

  Donor_Pair_Geno_Last <- vector("list", length(filter_Donor_Pair_Genotype))

  for (i in 1:length(filter_Donor_Pair_Genotype)) {
    genotype_len <- dim(filter_Donor_Pair_Genotype[[i]]$Genotype)[2]

    if (filter_Donor_Pair_Genotype[[i]]$is_Used_Time_Estimate == TRUE){
      donor <- names(filter_Donor_Pair_Genotype)[i]
      age <- as.numeric(cancer_time[donor,"age"])
      if(is.na(age)){
        next
      }

      Genotype <- filter_Donor_Pair_Genotype[[i]]$Genotype
      Genotype <- cbind(Genotype, array(NA, dim = c(dim(Genotype)[1],3)))
      colnames(Genotype)[(genotype_len + 1):(genotype_len + 3)] <- c("is_Have_mustE","healthTime", "sampleTime")
      Genotype[, genotype_len + 1] <- filter_Donor_Pair_Genotype[[i]]$is_Have_mustE

      filter_Donor_Pair_Genotype[[i]]$is_Have_mustE <- NULL

      Geno <- unlist(Genotype[1,1:genotype_len])
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

      for (j in 1:length(posetCluster)) {
        posetCluster[[j]] <- which(rateT$LB == max(rateT$LB[posetCluster[[j]]]))
      }

      posetClusterWeight <- sapply(1:length(posetCluster), function(j){
        gammaClusterRes$responsibilities[posetCluster[[j]],1]
      })
      posetClusterWeight <- posetClusterWeight/sum(posetClusterWeight)

      sol_res <- genotype_SolveQP(posetCluster, age, donor, cancer_time, Geno, length(driver_genes),
                                  filter_num, rateT, posetClusterWeight, isWeight = TRUE)

      if( any(as.numeric(sol_res$sol_last) <= 0.0001) ){
        next
      }

      Genotype[1, (genotype_len + 2):(genotype_len + 3)] <- sol_res$sol_last

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

      Donor_Pair_Geno_Last[[i]] <- list(Genotype = Genotype, age = age,
                                        is_Used_Time_Estimate = TRUE,
                                        sol_res = sol_res, age_Rate_Total = age_Rate, posetClusterWeight = posetClusterWeight
      )
      names(Donor_Pair_Geno_Last)[i] <- donor
    }
  }

  z1 <- z2 <- c()
  for (i in 1:length(Donor_Pair_Geno_Last)) {
    if(!is.null(Donor_Pair_Geno_Last[[i]])){
      z1 <- c(z1, Donor_Pair_Geno_Last[[i]]$sol_res$sol_last[1])
      z2 <- c(z2, Donor_Pair_Geno_Last[[i]]$sol_res$sol_last[2])
    }
  }
  z1[which(z1 <= 0.000001)] <- 0.000001
  z2[which(z2 <= 0.000001)] <- 0.000001
  z1_z2 <- mean(z1)/mean(z2)

  for (i in 1:length(filter_Donor_Pair_Genotype)) {
    genotype_len <- dim(filter_Donor_Pair_Genotype[[i]]$Genotype)[2]
    if (is.null(Donor_Pair_Geno_Last[[i]])){
      donor <- names(filter_Donor_Pair_Genotype)[i]
      if(donor %in% cBioCli[,"Patient ID"]){
        age <- mean(as.numeric(cBioCli[which(cBioCli[,"Patient ID"] == donor),"AGE"]) )
      }else{
        age <- as.numeric(unique(time_temp[which(time_temp[,"donor"] == donor), "age"]))
      }
      if(is.na(age)){
        next
      }

      Genotype <- filter_Donor_Pair_Genotype[[i]]$Genotype
      Genotype <- cbind(Genotype, array(NA, dim = c(dim(Genotype)[1],3)))
      colnames(Genotype)[(genotype_len + 1):(genotype_len + 3)] <- c("is_Have_mustE","healthTime", "sampleTime")

      Geno <- unlist(Genotype[1,1:genotype_len])
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

      if(length(posetCluster) == 0){
        posetClusterWeight <- NULL

        r_new <- z1_z2 + log(2)*min(rateT$LB)
        sol_Temp <- age/r_new
        Genotype[1, (genotype_len + 2):(genotype_len + 3)] <- c(z1_z2*sol_Temp, sol_Temp)

        Amat <- matrix(NA,nrow = 2, ncol = 2)
        Amat[1,1] <- -1
        Amat[1,2] <- 0
        Amat[2,1] <- 1
        Amat[2,2] <- min(rateT$LB)

        bvec <- c(-age,age)
        Amat <- rbind(Amat, matrix(c(1,0,0,1),nrow = 2, byrow = T))
        bvec <- c(bvec, 0, 0)

        sol_res <- list(sol_last = Genotype[1, (genotype_len + 2):(genotype_len + 3)], Amat = Amat, bvec = bvec)
        age_Rate <- log(2)*min(rateT$LB)

      }else{
        for (j in 1:length(posetCluster)) {
          posetCluster[[j]] <- which(rateT$LB == max(rateT$LB[posetCluster[[j]]]))
        }

        posetClusterWeight <- sapply(1:length(posetCluster), function(j){
          gammaClusterRes$responsibilities[posetCluster[[j]],1]
        })
        posetClusterWeight <- posetClusterWeight/sum(posetClusterWeight)

        r_new <- c()
        for (j in 1:length(posetCluster)) {
          time_L <- rateT$LB[posetCluster[[j]][1]]
          time_U <- rateT$UB[posetCluster[[j]][1]]
          if (time_U != Inf){
            r_new <- c(r_new, time_L + (time_U - time_L)*log(2))
          }else{
            r_new <- c(r_new, time_L + log(2))
          }
        }
        r_new <- z1_z2 + sum(r_new*posetClusterWeight)
        sol_Temp <- age/r_new

        Genotype[1, (genotype_len + 2):(genotype_len + 3)] <- c(z1_z2*sol_Temp, sol_Temp)


        Amat_L <- matrix(0,nrow = 2, ncol = 2)
        Amat_L[1,1] <- -1
        Amat_L[2,1] <- 1
        bvec_L <- c(-age,age)
        Amat_L_array <- rep(0, length(posetCluster))
        Weight1 <- posetClusterWeight
        Weight2 <- 0

        for (j in 1:length(posetCluster)) {
          Amat_L_array[j] <- -rateT$LB[posetCluster[[j]][1]]
          if(!is.infinite(rateT$UB[posetCluster[[j]][1]])){
            Amat_L[2,2] <- Amat_L[2,2] + rateT$UB[posetCluster[[j]][1]]*posetClusterWeight[j]
            Weight2 <- Weight2 + posetClusterWeight[j]
          }
        }
        if(Weight2 == 0){
          Amat_L[2,2] <- Inf
          Amat_L[1,2] <- sum(Amat_L_array*Weight1)/sum(Weight1)
        }else{
          Amat_L[2,2] <- Amat_L[2,2]/Weight2

          while ( abs(sum(Amat_L_array*Weight1)/sum(Weight1)) >=  Amat_L[2,2] ) {
            delta <- rep(0, length(posetCluster))
            for (j in 1:length(posetCluster)) {
              Weight2 <- Weight1
              Weight2[j] <- 0
              delta[j] <- abs(sum(Amat_L_array*Weight1)/sum(Weight1)) - abs(sum(Amat_L_array*Weight2)/sum(Weight2))
            }
            Weight1[which(delta == max(delta))] <- 0
          }
          Amat_L[1,2] <- sum(Amat_L_array*Weight1)/sum(Weight1)
        }



        Amat <- rbind(Amat_L, matrix(c(1,0,0,1),nrow = 2, byrow = T))
        indexS <- which(!(Amat[,2] == Inf | Amat[,2] == -Inf))
        Amat <- Amat[indexS,]
        bvec <- c(bvec_L, 0, 0)
        bvec <- bvec[indexS]

        sol_res <- list(sol_last = Genotype[1, (genotype_len + 2):(genotype_len + 3)], Amat = Amat, bvec = bvec)

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

      Donor_Pair_Geno_Last[[i]] <- list(Genotype = Genotype, age = age,
                                        is_Used_Time_Estimate = FALSE,
                                        sol_res = sol_res, age_Rate_Total = age_Rate, posetClusterWeight = posetClusterWeight
      )
      names(Donor_Pair_Geno_Last)[i] <- donor
    }
  }

  return(Donor_Pair_Geno_Last)


}



#' @title Gamma Parameter GS
#' @export
#' @description Gamma Parameter Gibbs Sampling
#' @param curr.alpha The current alpha
#' @param curr.beta The current beta
#' @param log.hyper.p The log hyperparameter of p
#' @param hyper.q The hyperparameter of q
#' @param hyper.r The hyperparameter of r
#' @param hyper.s The hyperparameter of s
#' @param x The current sampling data
#' @param gamma.sf The scaling perparameter for exploration
#' @returns A array of the new.alpha, new.beta, accepted
Gamma_Parameter_GS <- function(curr.alpha, curr.beta, log.hyper.p, hyper.q, hyper.r, hyper.s, x, gamma.sf) {
  alpha.proposed <- rgamma(n=1, curr.alpha * gamma.sf, gamma.sf)
  beta.proposed <- rgamma(n=1, curr.beta * gamma.sf, gamma.sf)

  proposal.dist.ratio <- dgamma(x = curr.alpha, shape = alpha.proposed * gamma.sf, rate = gamma.sf) *
    dgamma(x = curr.beta, shape = beta.proposed * gamma.sf, rate = gamma.sf) /
    dgamma(x = alpha.proposed, shape = curr.alpha * gamma.sf, rate = gamma.sf) /
    dgamma(x = beta.proposed, shape = curr.beta * gamma.sf, rate = gamma.sf)

  log.post.p <- log.hyper.p + sum(log(x))
  post.q <- hyper.q + sum(x)
  post.r <- hyper.r + length(x)
  post.s <- hyper.s + length(x)

  log.density.ratio <- (alpha.proposed - curr.alpha) * log.post.p -
    post.q * (beta.proposed - curr.beta) -
    post.r * (lgamma(alpha.proposed) - lgamma(curr.alpha)) +
    post.s * (alpha.proposed * log(beta.proposed) - curr.alpha * log(curr.beta))

  importance.ratio <- proposal.dist.ratio * exp(log.density.ratio)
  if (runif(n=1) < importance.ratio) {
    new.alpha <- alpha.proposed
    new.beta <- beta.proposed
    accepted <- TRUE
  } else {
    new.alpha <- curr.alpha
    new.beta <- curr.beta
    accepted <- FALSE
  }

  return(c(new.alpha, new.beta, accepted))
}



#' @title MH Sampling based on flipping
#' @export
#' @description MH Sampling based on flipping
#' @param Donor_Pair_Geno_Last The result from MH_Pretreatment_based_on_flipping
#' @param saveFile The file used to save result
#' @param isWeight If isWeight = true, use filter num to weighted average
#' @param iterations The sampling iterations
#' @param g.sf The scaling perparameter for exploration
#' @param maxSampleIter The max sample iterations of reject sampling
#' @param errorTolerance The error tolerance of reject sampling
#' @param errorToleranceGrowthRate The error tolerance growth rate of reject sampling
#' @param cores1 The parallel cores of sampling for chr time contain data
#' @param cores2 The parallel cores of sampling for no chr time contain data
#' @param show_size The show size
#' @param seed The seed
#' @returns A list of alpha, beta, lambda, accepted, chr_Time_Contain, no_chr_Time_Contain, time_Sampling
MH_Sampling_based_on_flipping <- function(Donor_Pair_Geno_Last, saveFile, isWeight = TRUE, iterations = 10000, g.sf = 10, maxSampleIter = 100,
                        errorTolerance = 0.01, errorToleranceGrowthRate = 1.1, cores1 = 6, cores2 = 10, show_size = 1000, seed = 1234){

  chr_Time_Contain <- list()
  no_chr_Time_Contain <- list()

  for (i in 1:length(Donor_Pair_Geno_Last)) {
    if (!is.null(Donor_Pair_Geno_Last[[i]])){
      if(Donor_Pair_Geno_Last[[i]]$is_Used_Time_Estimate){
        chr_Time_Contain[[length(chr_Time_Contain)+1]] <- list(donor = names(Donor_Pair_Geno_Last)[i],
                                                               sol_last = as.numeric(Donor_Pair_Geno_Last[[i]]$sol_res$sol_last),
                                                               Dmat = Donor_Pair_Geno_Last[[i]]$sol_res$Dmat,
                                                               dvec = Donor_Pair_Geno_Last[[i]]$sol_res$dvec,
                                                               Amat = Donor_Pair_Geno_Last[[i]]$sol_res$Amat,
                                                               bvec = Donor_Pair_Geno_Last[[i]]$sol_res$bvec,
                                                               errors_MSE = as.numeric(Donor_Pair_Geno_Last[[i]]$sol_res$errors_MSE),
                                                               age_Rate_Total = Donor_Pair_Geno_Last[[i]]$age_Rate_Total)
      }else{
        no_chr_Time_Contain[[length(no_chr_Time_Contain)+1]] <- list(donor = names(Donor_Pair_Geno_Last)[i],
                                                                     sol_last = as.numeric(Donor_Pair_Geno_Last[[i]]$sol_res$sol_last),
                                                                     Amat = Donor_Pair_Geno_Last[[i]]$sol_res$Amat,
                                                                     bvec = Donor_Pair_Geno_Last[[i]]$sol_res$bvec,
                                                                     age_Rate_Total = Donor_Pair_Geno_Last[[i]]$age_Rate_Total)
      }
    }
  }

  N_chr <- length(chr_Time_Contain)
  N_age <- length(no_chr_Time_Contain)
  N <- N_chr + N_age

  set.seed(seed)
  cores1 <- as.integer(cores1)
  cores2 <- as.integer(cores2)

  if(N_chr < cores1 | N_age < cores2){
    cat("Warring!\ncores1 should smaller than", N_chr,"\n")
    cat("Warring!\ncores2 should smaller than", N_age,"\n")
    return(NULL)
  }else{
    cat("cores1:  ", cores1, "\n")
    cat("cores2:  ", cores2, "\n")
  }

  time_Sampling <- array(data = NA, dim=c(iterations, N, 2))
  for (i in 1:N_chr) {
    time_Sampling[1, i,] <- chr_Time_Contain[[i]]$sol_last
  }
  for (i in 1:N_age) {
    time_Sampling[1, i + N_chr,] <- no_chr_Time_Contain[[i]]$sol_last
  }

  alpha <- beta <- lambda <- accepted <- rep(NA, iterations)

  gamma_fit <- fitdist(time_Sampling[1,,1],"gamma",method = "mle")
  alpha[1] <- gamma_fit$estimate[1]
  beta[1] <- gamma_fit$estimate[2]
  accepted[1] <- TRUE
  lambda[1] <- 1/mean(time_Sampling[1,,2])

  log.hyperparameter.p <- sum(log(time_Sampling[1,,1]))
  hyperparameter.q <- sum(time_Sampling[1,,1])
  hyperparameter.r <- hyperparameter.s <- N

  for (iter in 2:iterations) {
    if(show_size>0 & iter%%show_size==0){
      cat(iter, "   ", mean(lambda[floor(iter-0.9*show_size):(iter-1)]), "\n")
      saveRDS(list(alpha = alpha, beta = beta, lambda = lambda, accepted = accepted,
                   chr_Time_Contain = chr_Time_Contain, no_chr_Time_Contain = no_chr_Time_Contain,
                   time_Sampling = time_Sampling, iter = iter), saveFile)
    }

    temp <- mclapply(1:N_chr, function(i){
      Dmat <- chr_Time_Contain[[i]]$Dmat
      dvec <- chr_Time_Contain[[i]]$dvec
      Amat <- chr_Time_Contain[[i]]$Amat
      bvec <- chr_Time_Contain[[i]]$bvec
      index <- which(bvec != 0)
      Amat <- Amat[index,]
      bvec <- bvec[index]

      total_errors_MSE <- chr_Time_Contain[[i]]$errors_MSE
      isAccepte <- FALSE
      delta_total_errors_MSE <- abs(total_errors_MSE) * errorTolerance
      delta_bvec <- abs(bvec) * errorTolerance

      while (!isAccepte) {
        for (j in 1:maxSampleIter) {
          time_T <- c(rgamma(1, alpha[iter-1], beta[iter-1]), rexp(1, lambda[iter-1]))
          if(((time_T%*%Dmat - dvec)%*%time_T <= total_errors_MSE) & all(Amat%*%time_T >= bvec)){
            isAccepte <- TRUE
            break
          }
        }
        if(!isAccepte){
          total_errors_MSE <- total_errors_MSE + delta_total_errors_MSE
          bvec <- bvec - delta_bvec
          delta_total_errors_MSE <- delta_total_errors_MSE * errorToleranceGrowthRate
          delta_bvec <- delta_bvec * errorToleranceGrowthRate
        }
      }
      return(time_T)
    }, mc.cores = cores1, mc.cleanup = TRUE)

    for (i in 1:N_chr) {
      time_Sampling[iter, i,] <- temp[[i]]
    }

    time_Chr <- c(mean(time_Sampling[iter, 1:N_chr,1]), mean(time_Sampling[iter, 1:N_chr,2]))

    temp <- mclapply(1:N_age, function(i){
      Amat <- no_chr_Time_Contain[[i]]$Amat
      bvec <- no_chr_Time_Contain[[i]]$bvec
      index <- which(bvec != 0)
      Amat <- Amat[index,]
      bvec <- c(bvec[index], -Amat%*%time_Chr)
      Amat <- rbind(Amat, -Amat)
      isAccepte <- FALSE
      delta_bvec <- abs(bvec) * errorTolerance

      while (!isAccepte) {
        for (j in 1:maxSampleIter) {
          time_T <- c(rgamma(1, alpha[iter-1], beta[iter-1]), rexp(1, lambda[iter-1]))
          if(all(Amat%*%time_T >= bvec)){
            isAccepte <- TRUE
            break
          }
        }
        if(!isAccepte){
          bvec <- bvec - delta_bvec
          delta_bvec <- delta_bvec * errorToleranceGrowthRate
        }
      }
      return(time_T)
    }, mc.cores = cores2, mc.cleanup = TRUE)


    for (i in 1:N_age) {
      time_Sampling[iter, i + N_chr,] <- temp[[i]]
    }

    update <- Gamma_Parameter_GS(alpha[iter-1], beta[iter-1], log.hyperparameter.p, hyperparameter.q, hyperparameter.r, hyperparameter.s,
                                 time_Sampling[iter,,1], g.sf)
    alpha[iter] <- update[1]
    beta[iter] <- update[2]
    accepted[iter] <- update[3]

    lambda[iter] <- rgamma(1, shape = 0.01 + N, rate = 0.01 + sum(time_Sampling[iter,,2]))

  }

  saveRDS(list(alpha = alpha, beta = beta, lambda = lambda, accepted = accepted,
               chr_Time_Contain = chr_Time_Contain, no_chr_Time_Contain = no_chr_Time_Contain,
               time_Sampling = time_Sampling), saveFile)

  return(list(alpha = alpha, beta = beta, lambda = lambda, accepted = accepted,
              chr_Time_Contain = chr_Time_Contain, no_chr_Time_Contain = no_chr_Time_Contain,
              time_Sampling = time_Sampling))


}


#' @title Get Results From Sampling
#' @export
#' @description get res from sampling
#' @param MH_Res The result of MH-Sampling
#' @param start The start of useful data
#' @param end The end of useful data
#' @param skip_size The skip size
#' @returns A list of alpha, beta, lambda, donor, censored_Time, time_Sampling
get_Res_From_Sampling <- function(MH_Res, start, end, skip_size){
  alpha <- MH_Res$alpha[start]
  beta <- MH_Res$beta[start]
  lambda <- MH_Res$lambda[start]
  donor <- rep(NA, length(MH_Res$chr_Time_Contain) + length(MH_Res$no_chr_Time_Contain))
  censored_Time <- list()
  waiting_time <- list()
  evolve_time <- list()

  censored_Rate <- rep(NA, length(MH_Res$chr_Time_Contain) +length(MH_Res$no_chr_Time_Contain))

  N_chr <- length(MH_Res$chr_Time_Contain)
  index_S <- start

  for (j in 1:length(MH_Res$chr_Time_Contain)) {
    donor[j] <- MH_Res$chr_Time_Contain[[j]]$donor
    censored_Rate[j] <- MH_Res$chr_Time_Contain[[j]]$age_Rate_Total
    waiting_time[[j]] <- MH_Res$time_Sampling[start,j,1]
    evolve_time[[j]] <- MH_Res$time_Sampling[start,j,2]
    censored_Time[[j]] <- MH_Res$chr_Time_Contain[[j]]$age_Rate_Total * MH_Res$time_Sampling[start,j,2]
  }

  for (j in 1:length(MH_Res$no_chr_Time_Contain)) {
    donor[j + N_chr] <- MH_Res$no_chr_Time_Contain[[j]]$donor
    censored_Rate[j + N_chr] <- MH_Res$no_chr_Time_Contain[[j]]$age_Rate_Total
    waiting_time[[j + N_chr]] <- MH_Res$time_Sampling[start,j + N_chr,1]
    evolve_time[[j + N_chr]] <- MH_Res$time_Sampling[start,j + N_chr,2]
    censored_Time[[j + N_chr]] <- MH_Res$no_chr_Time_Contain[[j]]$age_Rate_Total * MH_Res$time_Sampling[start,j + N_chr,2]
  }

  index <- floor((end - start)/skip_size)
  for (i in 1:index) {
    alpha <- c(alpha, MH_Res$alpha[start + i*skip_size])
    beta <- c(beta, MH_Res$beta[start + i*skip_size])
    lambda <- c(lambda, MH_Res$lambda[start + i*skip_size])
    index_S <- c(index_S, start + i*skip_size)

    for (j in 1:length(MH_Res$chr_Time_Contain)) {
      waiting_time[[j]] <- c(waiting_time[[j]], MH_Res$time_Sampling[start + i*skip_size,j,1])
      evolve_time[[j]] <- c(evolve_time[[j]], MH_Res$time_Sampling[start + i*skip_size,j,2])
      censored_Time[[j]] <- c(censored_Time[[j]], MH_Res$chr_Time_Contain[[j]]$age_Rate_Total *
                                MH_Res$time_Sampling[start + i*skip_size,j,2])
    }

    for (j in 1:length(MH_Res$no_chr_Time_Contain)) {
      waiting_time[[j + N_chr]] <- c(waiting_time[[j + N_chr]], MH_Res$time_Sampling[start + i*skip_size,j + N_chr,1])
      evolve_time[[j + N_chr]] <- c(evolve_time[[j + N_chr]], MH_Res$time_Sampling[start + i*skip_size,j + N_chr,2])
      censored_Time[[j + N_chr]] <- c(censored_Time[[j + N_chr]], MH_Res$no_chr_Time_Contain[[j]]$age_Rate_Total *
                                        MH_Res$time_Sampling[start + i*skip_size,j + N_chr,2])
    }

  }

  time_Sampling <- MH_Res$time_Sampling[index_S,,]

  waiting_time_res <- sapply(1:length(waiting_time), function(i){
    mean(waiting_time[[i]])
  })
  evolve_time_res <- sapply(1:length(evolve_time), function(i){
    mean(evolve_time[[i]])
  })
  censored_Time_res <- sapply(1:length(censored_Time), function(i){
    mean(censored_Time[[i]])
  })

  return(list(alpha = mean(alpha), beta = mean(beta), lambda = mean(lambda),
              donor = donor, waiting_time = waiting_time_res, evolve_time = evolve_time_res,
              censored_Time = censored_Time_res, time_Sampling = time_Sampling))
}



#' @title Posterior Parameters Mvnormal
#' @export
#' @description sample posterior parameters from mvnormal distribution
#' @param priorParameters The prior parameters
#' @param x The data
#' @returns A list of mu_n, t_n, kappa_n, nu_n
PosteriorParameters_Mvnormal <- function(priorParameters, x) {
  nrow_x <- nrow(x)
  colmean_x <- colMeans(x)

  kappa0 <- priorParameters$kappa0
  mu0 <- priorParameters$mu0

  kappa_n <- kappa0 + nrow_x
  nu_n <- priorParameters$nu + nrow_x

  mu_n <- (kappa0 * mu0 + nrow_x * colmean_x)/(nrow_x + kappa0)

  sum_squares <- (nrow_x - 1) * var(x)
  sum_squares[is.na(sum_squares)] <- 0
  t_n <- priorParameters$Lambda +
    sum_squares +
    ((kappa0 * nrow_x)/(kappa0 + nrow_x)) * ((mu0 - colmean_x) %*% t(mu0 - colmean_x))

  return(list(mu_n = mu_n, t_n = t_n, kappa_n = kappa_n, nu_n = nu_n))
}



#' @title Posterior Draw Mvnormal for one
#' @export
#' @description sample one posterior from mvnormal distribution with simple solution for numerical error
#' @param priorParameters The prior parameters
#' @param x The data
#' @returns A list of mu, sig
PosteriorDraw_Mvnormal_for_1 <- function(priorParameters, x) {
  if (!is.matrix(x)) {
    x <- matrix(x, ncol = length(x))
  }

  post_parameters <- PosteriorParameters_Mvnormal(priorParameters, x)
  A <- post_parameters$t_n
  eigen_values <- eigen(A)$values
  if (any(eigen_values < 1e-10)){
    eigen_vectors <- eigen(A)$vectors
    eigen_values[which(eigen_values < 1e-10)] <- 1e-10
    A <- eigen_vectors %*% diag(eigen_values) %*% t(eigen_vectors)
  }

  sig <- rWishart(1, post_parameters$nu_n, A)[,,1]

  mu <- mvtnorm::rmvnorm(1,
                         post_parameters$mu_n,
                         ginv(post_parameters$kappa_n * sig))

  return(list(mu = mu, sig = sig/post_parameters$kappa_n^2))
}



#' @title Likelihood Beta
#' @export
#' @description calculate the likelihood of Beta distribution
#' @param x The data
#' @param param The parameters
#' @returns A number of likelihood
Likelihood_Beta <- function(x, param) {
  mu <- param[1]
  nu <- param[2]
  return(dbeta(x, mu * nu, (1 - mu) * nu))
}



#' @title Prior Density Beta
#' @export
#' @description calculate the prior density of Beta distribution
#' @param param The param
#' @param priorParameters The priorParameters
#' @returns A number of prior density
PriorDensity_Beta <- function(priorParameters, param) {

  muDensity <- dunif(param[1], 0, 1)
  nuDensity <- dgamma(1/param[2], priorParameters[1], priorParameters[2])

  return(as.numeric(muDensity * nuDensity))
}



#' @title Mh Parameter Proposal Beta
#' @export
#' @description MH parameter proposal of Beta distribution
#' @param update_sigma The update sigma
#' @param old_params The old params
#' @returns The new params
MhParameterProposal_Beta <- function(update_sigma, old_params) {

  new_params <- old_params

  new_params[1] <- old_params[1] + rnorm(1, 0, update_sigma)

  if (new_params[1] > 1 | new_params[1] < 0) {
    new_params[1] <- old_params[1]
  }

  new_params[2] <- abs(old_params[2] + rnorm(1, 0, update_sigma))

  return(new_params)
}



#' @title MetropolisHastings Beta
#' @export
#' @description Metropolis-Hastings sampling of Beta distribution
#' @param priorParameters The prior parameters
#' @param x The data
#' @param update_sigma The update sigma
#' @param start_pos The start position
#' @returns The param
MetropolisHastings_Beta <- function(priorParameters, x, update_sigma, start_pos) {

  old_param <- start_pos

  old_prior <- log(PriorDensity_Beta(priorParameters, old_param))
  old_Likelihood <- sum(log(Likelihood_Beta(x, old_param)))

  prop_param <- MhParameterProposal_Beta(update_sigma, old_param)
  new_prior <- log(PriorDensity_Beta(priorParameters, prop_param))
  new_Likelihood <- sum(log(Likelihood_Beta(x, prop_param)))

  accept_prob <- min(1, exp(new_prior + new_Likelihood - old_prior - old_Likelihood))

  if (is.na(accept_prob) | !length(accept_prob) ) {
    accept_prob <- 0
  }

  if (runif(1) < accept_prob) {
    old_param <- prop_param
  }

  return(old_param)
}



#' @title Fitness Prior Inference Patient Gene Stage1
#' @export
#' @description Inference fitness for patients with gene neoantigen data
#' @param fitness_Inference_Input_Gene The fitness inference gene input data
#' @param sample_Patient The patient info
#' @param driver_genes The driver genes in CBN that don't need to update the parameters
#' @param min_Patient The genes that will be used at least min_Patient, default 4
#' @param sigma_p_square The initial variance of the proposed distribution in the adaptive MCMC, default 0.0001
#' @param neighbor_size The neighbor size in knn prior stage, default 10
#' @param history_num The size of history data in history prior stage, default 100
#' @param history_skip The thining num in history prior stage, default 10
#' @param history_size The total size before history prior stage, default 25000
#' @param knn_used_in_history_size The size of knn prior stage used for history prior stage, default 20000
#' @param gene_history_num The gene size of history data in history prior stage, default 100
#' @param gene_history_size The total gene size before history prior stage, default 25000
#' @param gene_used_in_history_size The gene size of knn prior stage used for history prior stage, default 20000
#' @param max_iter The max sampling number, default 40000
#' @param cellular_clock The colname that used as cellular_clock in sample_Patient, default "censored_Rate"
#' @param target_recept_rate The target recept rate in the adaptive MCMC, default 0.234
#' @param f_b_scaling_update_strength The scaling update strength in the adaptive MCMC, default 10
#' @param r_update_sigma The sigma of the proposed distribution in the gene update stage, default 0.01
#' @param show_size The show size, default 1000
#' @param seed The seed
#' @returns The list of fitness_Inference_Input_Gene, sample_Patient, gene_Temp, driver_gene_Temp,
#'          fp_bp_Sampling, accepted_fp_bp_Sampling, fp_bp_scaling_Sampling, data_List, gene_List,
#'          pi_alpha, pi_beta, pi_Sampling, scaleNum, driver_gene_Temp_t_min
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




#' @title Fitness Inference Patient Gene Stage2
#' @export
#' @description Inference fitness for patients with gene neoantigen data and fixed prior
#' @param prior_data The prior data from stage1
#' @param sigma_p_square The initial variance of the proposed distribution in the adaptive MCMC, default 0.0001
#' @param max_iter The max sampling number, default 10000
#' @param cellular_clock The colname that used as cellular_clock in sample_Patient, default "censored_Rate"
#' @param target_recept_rate The target recept rate in the adaptive MCMC, default 0.234
#' @param f_b_scaling_update_strength The scaling update strength in the adaptive MCMC, default 10
#' @param r_update_sigma The sigma of the proposed distribution in the gene update stage, default 0.01
#' @param show_size The show size, default 1000
#' @param seed The seed
#' @returns The list of fitness_Inference_Input_Gene, sample_Patient, gene_Temp, driver_gene_Temp,
#'          fp_bp_Sampling, accepted_fp_bp_Sampling, fp_bp_scaling_Sampling, data_List, gene_List,
#'          pi_alpha, pi_beta, pi_Sampling, scaleNum, driver_gene_Temp_t_min
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

