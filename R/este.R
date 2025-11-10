#' @title Estimate epsilon fastly for multiple data and event sets
#' @export
#'
#' @description Fast epsilon estimation for the hidden conjunctive Bayesian
#' network model (H-CBN), use a subset of the events to estimate epsilon when
#' the number of events is greater than 7
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
#' @param eps The initial matrix of epslion
#'
#' @param setD The description of dataset index
#'
#' @param eventD The description of eventset index
#'
#' @param threshold The threshold of epslion estimate. Default: 0.05
#'
#' @param threshold1 The threshold of epslion estimate method selection. Default: 7
#'
#' @param n_p The max num for estimate epsilon be. Default: the value of threshold1
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
#' @return A value of epsilon
estimate_Epsilon <- function(pat, isF, isCE = NULL, eps =NULL, setD, eventD,
        threshold = 0.05, threshold1 = 7L, n_p = 0L, T = 10.0, N_iter = 150L, lambdaS = 1.0, thrds=1L,
        R = 1L) {

  #---------------------------------------------------------
  if (!is.matrix(pat) | !is.integer(pat)){
    cat("pat must be an integer matrix!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.matrix(isF) | !is.integer(isF)){
    cat("isF must be an integer matrix!\n")
    break
  }
  dataSetN = nrow(isF)
  eventSetN = ncol(isF)

  #---------------------------------------------------------
  if (!is.null(isCE)){
    if (!is.matrix(isCE) | !is.integer(isCE)){
      cat("isCE must be an integer matrix!\n")
      break
    }else{
      if ((nrow(isCE) != dataSetN) | (ncol(isCE) != eventSetN)){
        cat("dim of isCE is error!\n")
        break
      }
    }
  }else{
    isCE = isF
  }

  #---------------------------------------------------------
  if (!is.null(eps)){
    if (!is.matrix(eps) | !is.numeric(eps)){
      cat("eps must be a numeric matrix!\n")
      break
    }else{
      if ((nrow(eps) != dataSetN) | (ncol(eps) != eventSetN)){
        cat("dim of eps is error!\n")
        break
      }
    }
  }else{
    eps = matrix(as.numeric(0),nrow = dataSetN,ncol = eventSetN)
  }

  #---------------------------------------------------------
  if (!is.matrix(setD) | !is.integer(setD)){
    cat("setD must be an integer matrix!\n")
    break
  }else{
    if ((nrow(setD) != dataSetN) | (ncol(setD) != 2)){
      cat("dim of setD is error!\n")
      break
    }
  }

  #---------------------------------------------------------
  if (!is.matrix(eventD) | !is.integer(eventD)){
    cat("eventD must be an integer matrix!\n")
    break
  }else{
    if ((nrow(eventD) != eventSetN) | (ncol(eventD) != 2)){
      cat("dim of eventD is error!\n")
      break
    }
  }

  #---------------------------------------------------------
  if (!is.numeric(threshold)){
    cat("threshold must be numeric!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.integer(threshold1)){
    cat("threshold1 must be integer!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.integer(n_p)){
    cat("n_p must be integer!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.numeric(T)){
    cat("T must be numeric!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.integer(N_iter)){
    cat("N_iter must be an integer!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.numeric(lambdaS)){
    cat("lambdaS must be numeric!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.integer(thrds)){
    cat("thrds must be an integer!\n")
    break
  }else{
    thrds = min(thrds,max(1,parallel::detectCores()-2));
    cat("Parallel Threads: ",thrds,"\n")
  }

  #---------------------------------------------------------
  if (!is.integer(R)){
    cat("R must be an integer!\n")
    break
  }else{
    R = min(R,1);
  }

  #---------------------------------------------------------
  .Call('_este_estimate_epsilon_wrapper', PACKAGE = 'este', pat, isF, isCE, eps,
        setD, eventD, threshold, threshold1, n_p, T, N_iter, lambdaS, thrds, R)

}



#' @title Estimate epsilon based on poset and lambda
#' @export
#'
#' @description Estimate epsilon based on poset and lambda for the hidden
#' conjunctive Bayesian network model (H-CBN), the estimation based on mccbn
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
#' @param lambda The lambda of each event
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
#' @param seed seed for reproducibility
#'
#' @return A list of epsilon1, epsilon2, lambda
estimate_Epsilon_based_on_Poset_and_Lambda <- function(pat, poset, isF, eps, setD, eventD, lambda, lambdaS = 1.0,
                            L = 100L, sampling = c('forward','add-remove',
                                                   'backward','bernoulli','pool'), maxIter=100L,
                            updateStepSize=20L, tol=0.001, maxLambda=1e6,
                            neighborhoodDist=1L, thrds=1L, seed=NULL) {

  #---------------------------------------------------------
  if (!is.matrix(pat) | !is.integer(pat)){
    cat("pat must be an integer matrix!\n")
    break
  }
  E <- ncol(pat)-1

  #---------------------------------------------------------
  if (!is.matrix(poset) | !is.integer(poset)){
    cat("poset must be an integer matrix!\n")
    break
  }else{
    if ((nrow(poset) != ncol(poset)) | (ncol(poset) != E)){
      cat("dim of poset is error!\n")
      break
    }
  }

  #---------------------------------------------------------
  if (!is.matrix(isF) | !is.integer(isF)){
    cat("isF must be an integer matrix!\n")
    break
  }
  dataSetN = nrow(isF)
  eventSetN = ncol(isF)

  #---------------------------------------------------------
  if (!is.matrix(eps) | !is.numeric(eps)){
    cat("eps must be a numeric matrix!\n")
    break
  }else{
    if ((nrow(eps) != dataSetN) | (ncol(eps) != eventSetN)){
      cat("dim of eps is error!\n")
      break
    }
  }

  #---------------------------------------------------------
  if (!is.matrix(setD) | !is.integer(setD)){
    cat("setD must be an integer matrix!\n")
    break
  }else{
    if ((nrow(setD) != dataSetN) | (ncol(setD) != 2)){
      cat("dim of setD is error!\n")
      break
    }
  }

  #---------------------------------------------------------
  if (!is.matrix(eventD) | !is.integer(eventD)){
    cat("eventD must be an integer matrix!\n")
    break
  }else{
    if ((nrow(eventD) != eventSetN) | (ncol(eventD) != 2)){
      cat("dim of eventD is error!\n")
      break
    }
  }

  #---------------------------------------------------------
  if (!is.numeric(lambda)){
    cat("lambda must be numeric vector!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.numeric(lambdaS)){
    cat("lambdaS must be numeric!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.integer(L)){
    cat("L must be an integer!\n")
    break
  }

  #---------------------------------------------------------
  sampling <- match.arg(sampling)
  if (!is.integer(maxIter)){
    cat("maxIter must be an integer!\n")
    break
  }
  if (!is.integer(updateStepSize)){
    cat("updateStepSize must be an integer!\n")
    break
  }
  if (updateStepSize > maxIter)
    updateStepSize <- as.integer(maxIter / 5)

  #---------------------------------------------------------
  if (!is.numeric(tol)){
    cat("tol must be numeric!\n")
    break
  }
  if (!is.numeric(maxLambda)){
    cat("maxLambda must be numeric!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.integer(neighborhoodDist)){
    cat("neighborhoodDist must be a integer!\n")
    break
  }
  if (!is.integer(thrds)){
    cat("thrds must be an integer!\n")
    break
  }else{
    thrds = min(thrds,max(1,parallel::detectCores()-2));
    cat("Parallel Threads: ",thrds,"\n")
  }

  #---------------------------------------------------------
  if (!is.null(seed)){
    if (!is.integer(seed)){
      cat("seed must be a integer!\n")
      break
    }
  }else{
    seed <- sample.int(3e4, 1)
  }

  #---------------------------------------------------------
  .Call('_este_estimate_epsilon_based_on_poset_and_lambda_wrapper', PACKAGE = 'este', pat, poset, isF,
        eps, setD, eventD, lambda, lambdaS, L, sampling, maxIter, updateStepSize,
        tol, maxLambda, neighborhoodDist, thrds, seed)

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
find_Poset <- function(pat, isF, isCE = NULL, eps, setD, eventD, Fine_Tune_Num=2L, thrds=1L, threshold = 0.0001, threshold2 = 0.01,
                       gmm_min_cluster = 2, gmm_convenience_threshold = 0.001, gmm_filter_threshold = 0.002,
             n_p = 6L, is_update_eps = FALSE, T = 10.0, N_iter = 150L, lambdaS = 1.0, R = 1L) {

  #---------------------------------------------------------
  if (!is.matrix(pat) | !is.integer(pat)){
    cat("pat must be an integer matrix!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.matrix(isF) | !is.integer(isF)){
    cat("isF must be an integer matrix!\n")
    break
  }
  dataSetN = nrow(isF)
  eventSetN = ncol(isF)

  #---------------------------------------------------------
  if (!is.null(isCE)){
    if (!is.matrix(isCE) | !is.integer(isCE)){
      cat("isCE must be an integer matrix!\n")
      break
    }else{
      if ((nrow(isCE) != dataSetN) | (ncol(isCE) != eventSetN)){
        cat("dim of isCE is error!\n")
        break
      }
    }
  }else{
    isCE = isF
  }

  #---------------------------------------------------------
  if (!is.matrix(eps) | !is.numeric(eps)){
    cat("eps must be a numeric matrix!\n")
    break
  }else{
    if ((nrow(eps) != dataSetN) | (ncol(eps) != eventSetN)){
      cat("dim of eps is error!\n")
      break
    }
  }

  #---------------------------------------------------------
  if (!is.matrix(setD) | !is.integer(setD)){
    cat("setD must be an integer matrix!\n")
    break
  }else{
    if ((nrow(setD) != dataSetN) | (ncol(setD) != 2)){
      cat("dim of setD is error!\n")
      break
    }
  }

  #---------------------------------------------------------
  if (!is.matrix(eventD) | !is.integer(eventD)){
    cat("eventD must be an integer matrix!\n")
    break
  }else{
    if ((nrow(eventD) != eventSetN) | (ncol(eventD) != 2)){
      cat("dim of eventD is error!\n")
      break
    }
  }

  #---------------------------------------------------------
  if (!is.integer(Fine_Tune_Num)){
    cat("Fine_Tune_Num must be integer!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.integer(thrds)){
    cat("thrds must be an integer!\n")
    break
  }else{
    thrds = min(thrds,max(1,parallel::detectCores()-2));
    cat("Parallel Threads: ",thrds,"\n")
  }

  #---------------------------------------------------------
  if (!is.numeric(threshold)){
    cat("threshold must be numeric!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.numeric(threshold2)){
    cat("threshold2 must be numeric!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.integer(n_p)){
    cat("n_p must be integer!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.logical(is_update_eps)){
    cat("is_update_eps must be logical!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.numeric(T)){
    cat("T must be numeric!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.integer(N_iter)){
    cat("N_iter must be an integer!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.numeric(lambdaS)){
    cat("lambdaS must be numeric!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.integer(R)){
    cat("R must be an integer!\n")
    break
  }else{
    R = min(R,1);
  }

  #---------------------------------------------------------
  .Call('_este_find_poset_wrapper', PACKAGE = 'este', pat, thrds, isF, isCE, eps,
        setD, eventD, Fine_Tune_Num, threshold, threshold2, gmm_min_cluster,
        gmm_convenience_threshold, gmm_filter_threshold, n_p, is_update_eps,
        T, N_iter, lambdaS, R)
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
#'
#' @return A list of epsilon1, epsilon2, lambda
estimate_Lambda <- function(pat, poset, isF, eps, setD, eventD, lambdaS = 1.0,
                  L = 100L, sampling = c('forward','add-remove',
                  'backward','bernoulli','pool'), maxIter=100L,
                  updateStepSize=20L, tol=0.001, maxLambda=1e6,
                  neighborhoodDist=1L, thrds=1L, is_update_eps = FALSE, seed=NULL) {

  #---------------------------------------------------------
  if (!is.matrix(pat) | !is.integer(pat)){
    cat("pat must be an integer matrix!\n")
    break
  }
  E <- ncol(pat)-1

  #---------------------------------------------------------
  if (!is.matrix(poset) | !is.integer(poset)){
    cat("poset must be an integer matrix!\n")
    break
  }else{
    if ((nrow(poset) != ncol(poset)) | (ncol(poset) != E)){
      cat("dim of poset is error!\n")
      break
    }
  }

  #---------------------------------------------------------
  if (!is.matrix(isF) | !is.integer(isF)){
    cat("isF must be an integer matrix!\n")
    break
  }
  dataSetN = nrow(isF)
  eventSetN = ncol(isF)

  #---------------------------------------------------------
  if (!is.matrix(eps) | !is.numeric(eps)){
    cat("eps must be a numeric matrix!\n")
    break
  }else{
    if ((nrow(eps) != dataSetN) | (ncol(eps) != eventSetN)){
      cat("dim of eps is error!\n")
      break
    }
  }

  #---------------------------------------------------------
  if (!is.matrix(setD) | !is.integer(setD)){
    cat("setD must be an integer matrix!\n")
    break
  }else{
    if ((nrow(setD) != dataSetN) | (ncol(setD) != 2)){
      cat("dim of setD is error!\n")
      break
    }
  }

  #---------------------------------------------------------
  if (!is.matrix(eventD) | !is.integer(eventD)){
    cat("eventD must be an integer matrix!\n")
    break
  }else{
    if ((nrow(eventD) != eventSetN) | (ncol(eventD) != 2)){
      cat("dim of eventD is error!\n")
      break
    }
  }

  #---------------------------------------------------------
  if (!is.numeric(lambdaS)){
    cat("lambdaS must be numeric!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.integer(L)){
    cat("L must be an integer!\n")
    break
  }

  #---------------------------------------------------------
  sampling <- match.arg(sampling)
  if (!is.integer(maxIter)){
    cat("maxIter must be an integer!\n")
    break
  }
  if (!is.integer(updateStepSize)){
    cat("updateStepSize must be an integer!\n")
    break
  }
  if (updateStepSize > maxIter)
    updateStepSize <- as.integer(maxIter / 5)

  #---------------------------------------------------------
  if (!is.numeric(tol)){
    cat("tol must be numeric!\n")
    break
  }
  if (!is.numeric(maxLambda)){
    cat("maxLambda must be numeric!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.integer(neighborhoodDist)){
    cat("neighborhoodDist must be a integer!\n")
    break
  }
  if (!is.integer(thrds)){
    cat("thrds must be an integer!\n")
    break
  }else{
    thrds = min(thrds,max(1,parallel::detectCores()-2));
    cat("Parallel Threads: ",thrds,"\n")
  }

  #---------------------------------------------------------
  if (!is.logical(is_update_eps)){
    cat("is_update_eps must be logical!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.null(seed)){
    if (!is.integer(seed)){
      cat("seed must be a integer!\n")
      break
    }
  }else{
    seed <- sample.int(3e4, 1)
  }

  #---------------------------------------------------------
  .Call('_este_estimate_lambda_wrapper', PACKAGE = 'este', pat, poset, isF,
        eps, setD, eventD, lambdaS, L, sampling, maxIter, updateStepSize,
        tol, maxLambda, neighborhoodDist, thrds, is_update_eps, seed)


}





#' @title Large scale continuous time conjunctive Bayesian network (lsctCBN)
#' for multiple data and event sets
#' @export
#'
#' @description Determine whether the genotype is compatible with the poset
#'
#' @param genotype A matrix
#'
#' @param poset A matrix
is_Compatible <- function(genotype, poset) {

  #---------------------------------------------------------
  if (!is.array(genotype)){
    cat("genotype must be an array!\n")
    break
  }

  #---------------------------------------------------------
  if (!is.matrix(poset)){
    cat("poset must be an matrix!\n")
    break
  }

  .Call('_este_isCompatible_wrapper', PACKAGE = 'este', genotype, poset)

}






#' @title Large scale continuous time conjunctive Bayesian network (lsctCBN)
#' for multiple data and event sets
#' @export
#'
#' @description Sample time for censor age
#'
#' @param data A data matrix
#'
#' @param alpha Alpha parameter of exp distribution
#'
#' @param beta Beta parameter of exp distribution
#'
#' @param lambda Parameter of exp distribution
#'
#' @param maxSampleIter Max ssample iter
#'
#' @param rateE1 Tolerated amplification rate
#'
#' @param thrds Parallel thrds
#'
sample_Age_T <- function(data, alpha, beta, lambda, maxSampleIter, rateE1, thrds, seed) {

  .Call('_este_sample_Age_T_wrapper', PACKAGE = 'este', data, alpha, beta, lambda, maxSampleIter, rateE1, thrds, seed)

}




#' @title Large scale continuous time conjunctive Bayesian network (lsctCBN)
#' for multiple data and event sets
#' @export
#'
#' @description Sample time for censor age (lightweight)
#'
#' @param data A data matrix
#'
#' @param alpha Alpha parameter of exp distribution
#'
#' @param beta Beta parameter of exp distribution
#'
#' @param lambda Parameter of exp distribution
#'
#' @param maxSampleIter Max ssample iter
#'
#' @param rateE1 Tolerated amplification rate
#'
#' @param thrds Parallel thrds
#'
sample_Age_TLW <- function(data, alpha, beta, lambda, thrds) {

  .Call('_este_sample_Age_TLW_wrapper', PACKAGE = 'este', data, alpha, beta, lambda, thrds)

}











