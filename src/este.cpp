/** title: Wrapper of este
 *
 * @author Hu Changbao
 * @email 1437894182@qq.com
 */
#include <Rcpp.h>
#include "interfaceFun.h"

using namespace Rcpp;
using namespace Eigen;
using namespace std;


//-----------------------------------------------------------------------------------------
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
Rcpp::NumericMatrix estimate_epsilon_wrapper(Rcpp::IntegerMatrix& ipat,
                  Rcpp::IntegerMatrix& iisF, Rcpp::IntegerMatrix& iisCE,
                  Rcpp::NumericMatrix& ieps, Rcpp::IntegerMatrix& isetD,
                  Rcpp::IntegerMatrix& ieventD, Rcpp::NumericVector ithreshold,
                  Rcpp::IntegerVector ithreshold1, Rcpp::IntegerVector in_p,
                  Rcpp::NumericVector iT, Rcpp::IntegerVector iN_iter,
                  Rcpp::NumericVector ilambdaS, Rcpp::IntegerVector ithrds,
                  Rcpp::IntegerVector iR) {

  int N = ipat.nrow();
  int n = ipat.ncol()-1;

  dataD dD;
  /*      genotype poset setNum eventNum     */
  dD.pat = as<MyIntMatrix>(ipat);
  dD.poset.resize(n+1,n+1);
  dD.poset.setZero();
  dD.setNum = iisF.nrow();
  dD.eventNum = iisF.ncol();
  /*      isF isCE epsilon poset_temp      */
  dD.isF = as<MyBoolMatrix>(iisF);
  dD.isCE = as<MyBoolMatrix>(iisCE);
  dD.epsilon = as<MyDoubleMatrix>(ieps);
  dD.poset_temp.resize(n+1,n+1);
  dD.poset_temp.setZero();
  /*      setD eventD colEvent     */
  for (int i = 0; i < isetD.nrow(); ++i) {
    dD.setD.push_back(doubleD(isetD(i,0),isetD(i,1)));
  }
  for (int i = 0; i < ieventD.nrow(); ++i) {
    dD.eventD.push_back(doubleD(ieventD(i,0),ieventD(i,1)));
  }
  dD.colEvent.push_back(-1);
  for (int i = 0; i < dD.eventNum; ++i) {
    for (int j = dD.eventD[i].F; j <= dD.eventD[i].S; ++j) {
      dD.colEvent.push_back(i);
    }
  }
  /*      GENOTYPE     */
  std::vector<int**> GENOTYPE(30);
  GENOTYPE[3] = precompute_binary(3);
  /*      parallel     */
  int thrds = as<int>(ithrds);
  #ifdef _OPENMP
    omp_set_num_threads(thrds);
    omp_set_dynamic(false);
  #endif
  /*      threshold T N_iter lambdaS R      */
  double threshold = max(0.5,as<double>(ithreshold));
  int threshold1 = as<int>(ithreshold1);
  int n_p = as<int>(in_p);
  double T = as<double>(iT);
  int N_iter = as<int>(iN_iter);
  double lambdaS = as<double>(ilambdaS);
  int R = as<int>(iR);
  /*      estimate epsilon     */
  estimate_epsilon_mix(dD, thrds, GENOTYPE, threshold, threshold1, n_p, T, N_iter, lambdaS, R);
  /*      GENOTYPE Free    */
  for (int i = 0; i < GENOTYPE.size(); ++i) {
    if (GENOTYPE[i] != NULL){
      free_int_Matrix(GENOTYPE[i], pow2(i));
    }
  }

  return wrap(dD.epsilon);
}



//-----------------------------------------------------------------------------------------
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
Rcpp::List estimate_epsilon_based_on_poset_and_lambda_wrapper(Rcpp::IntegerMatrix& ipat,
                                   Rcpp::IntegerMatrix& iposet, Rcpp::IntegerMatrix& iisF,
                                   Rcpp::NumericMatrix& ieps, Rcpp::IntegerMatrix& isetD,
                                   Rcpp::IntegerMatrix& ieventD, Rcpp::NumericVector& ilambda,
                                   Rcpp::NumericVector ilambdaS,
                                   Rcpp::IntegerVector iL, Rcpp::CharacterVector isampling,
                                   Rcpp::IntegerVector imaxIter, Rcpp::IntegerVector iupdateStepSize,
                                   Rcpp::NumericVector itol, Rcpp::NumericVector imaxLambda,
                                   Rcpp::IntegerVector ineighborhoodDist, Rcpp::IntegerVector ithrds,
                                   Rcpp::IntegerVector iseed) {


  /*      obs poset eps setD    */
  int rowN = 0;
  for (int i = 0; i < iisF.nrow(); ++i) {
    bool is_Fill = true;
    for (int j = 0; j < iisF.ncol(); ++j) {
      if (iisF(i,j) == 0){
        is_Fill = false;
        break;
      }
    }
    if (is_Fill){
      rowN ++;
    } else{
      break;
    }
  }
  int N_total = 0;
  int n = iposet.nrow();
  MyDoubleMatrix eps(rowN,n);
  std::vector<doubleD> setD(rowN);

  std::vector<int> colEvent;
  colEvent.push_back(-1);
  for (int i = 0; i < ieventD.nrow(); ++i) {
    for (int j = ieventD(i,0); j <= ieventD(i,1); ++j) {
      colEvent.push_back(i);
    }
  }

  for (int i = 0; i < rowN; ++i) {
    N_total += isetD(i,1) - isetD(i,0) +1;
    for (int j = 0; j < n; ++j) {
      eps(i,j) = ieps(i,colEvent[j+1]);
      setD[i] = doubleD(isetD(i,0),isetD(i,1));
    }
  }

  MyBoolMatrix obs;
  obs.resize(N_total,n);
  for (int i = 0; i < N_total; ++i) {
    for (int j = 0; j < obs.cols(); ++j) {
      obs(i,j) = ipat(i,j+1);
    }
  }

  MyIntMatrix poset = as<MyIntMatrix>(iposet);

  /*      eventD     */
  std::vector<doubleD> eventD(ieventD.nrow());
  for (int i = 0; i < ieventD.nrow(); ++i) {
    eventD[i] = doubleD(ieventD(i,0),ieventD(i,1));
  }

  /*      L sampling seed thrds lambdaS      */
  unsigned int L = as<unsigned int>(iL);
  std::string sampling = as<std::string>(isampling);
  int seed = as<int>(iseed);
  int thrds = as<int>(ithrds);
  double lambdaS = as<double>(ilambdaS);

  /*      max_iter update_step_size tol max_lambda neighborhood_dist      */
  int maxIter = as<int>(imaxIter);
  int updateStepSize = as<int>(iupdateStepSize);
  double tol = as<double>(itol);
  double maxLambda = as<double>(imaxLambda);
  unsigned int neighborhoodDist = as<unsigned int>(ineighborhoodDist);

  #ifdef _OPENMP
    omp_set_num_threads(thrds);
    omp_set_dynamic(false);
  #endif

  VectorXd lambdaRes(poset.rows());
  for (int i = 0; i < poset.rows(); ++i) {
    lambdaRes(i) = ilambda[i];
  }


  /*      estimate lambda     */
  estimate_epsilon_based_on_poset_and_lambda(lambdaRes, poset, obs, eps, setD, eventD, L, sampling, seed,
                  thrds, lambdaS, maxIter, updateStepSize, tol, maxLambda, neighborhoodDist);

  return Rcpp::List::create(Named("epsilon1")=ieps, Named("epsilon2")=eps, Named("lambda")=lambdaRes);

}


//-----------------------------------------------------------------------------------------
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
Rcpp::IntegerMatrix find_poset_wrapper(Rcpp::IntegerMatrix& ipat, Rcpp::IntegerVector ithrds,
                  Rcpp::IntegerMatrix& iisF, Rcpp::IntegerMatrix& iisCE,
                  Rcpp::NumericMatrix& ieps, Rcpp::IntegerMatrix& isetD,
                  Rcpp::IntegerMatrix& ieventD, Rcpp::IntegerVector iFine_Tune_Num,
                  Rcpp::NumericVector ithreshold, Rcpp::NumericVector ithreshold2,
                  Rcpp::IntegerVector igmm_min_cluster,
                  Rcpp::NumericVector igmm_convenience_threshold,
                  Rcpp::NumericVector igmm_filter_threshold,
                  Rcpp::IntegerVector in_p, Rcpp::LogicalVector iis_update_eps,
                  Rcpp::NumericVector iT,
                  Rcpp::IntegerVector iN_iter, Rcpp::NumericVector ilambdaS,
                  Rcpp::IntegerVector iR) {

  int N = ipat.nrow();
  int n = ipat.ncol()-1;

  dataD dD;
  /*      genotype poset setNum eventNum     */
  dD.pat = as<MyIntMatrix>(ipat);
  dD.poset.resize(n+1,n+1);
  dD.poset.setZero();
  dD.setNum = iisF.nrow();
  dD.eventNum = iisF.ncol();
  /*      isF isCE epsilon poset_temp      */
  dD.isF = as<MyBoolMatrix>(iisF);
  dD.isCE = as<MyBoolMatrix>(iisCE);
  dD.epsilon = as<MyDoubleMatrix>(ieps);
  dD.poset_temp.resize(n+1,n+1);
  dD.poset_temp.setZero();
  /*      setD eventD colEvent     */
  for (int i = 0; i < isetD.nrow(); ++i) {
    dD.setD.push_back(doubleD(isetD(i,0),isetD(i,1)));
  }
  for (int i = 0; i < ieventD.nrow(); ++i) {
    dD.eventD.push_back(doubleD(ieventD(i,0),ieventD(i,1)));
  }
  dD.colEvent.push_back(-1);
  for (int i = 0; i < dD.eventNum; ++i) {
    for (int j = dD.eventD[i].F; j <= dD.eventD[i].S; ++j) {
      dD.colEvent.push_back(i);
    }
  }
  /*      GENOTYPE     */
  std::vector<int**> GENOTYPE(30);
  GENOTYPE[3] = precompute_binary(3);
  /*      parallel     */
  int thrds = as<int>(ithrds);
  #ifdef _OPENMP
    omp_set_num_threads(thrds);
    omp_set_dynamic(false);
  #endif
  /*      Fine_Tune_Num threshold n_p T N_iter lambdaS R      */
  int Fine_Tune_Num = as<int>(iFine_Tune_Num);
  double threshold = as<double>(ithreshold);
  double threshold2 = as<double>(ithreshold2);
  int gmm_min_cluster = as<int>(igmm_min_cluster);
  double gmm_convenience_threshold = as<double>(igmm_convenience_threshold);
  double gmm_filter_threshold = as<double>(igmm_filter_threshold);
  // std::vector<int> n_p = as<std::vector<int>>(in_p);
  int n_p = as<int>(in_p);
  bool is_update_eps = as<bool>(iis_update_eps);
  double T = as<double>(iT);
  int N_iter = as<int>(iN_iter);
  double lambdaS = as<double>(ilambdaS);
  int R = as<int>(iR);
  /*      find poset     */
  find_poset_mix(dD, thrds, GENOTYPE, Fine_Tune_Num, threshold, threshold2, gmm_min_cluster, gmm_convenience_threshold,
                 gmm_filter_threshold, n_p, is_update_eps, T, N_iter, lambdaS, R);
  /*      GENOTYPE Free    */
  for (int i = 0; i < GENOTYPE.size(); ++i) {
    if (GENOTYPE[i] != NULL){
      free_int_Matrix(GENOTYPE[i], pow2(i));
    }
  }

  return wrap(dD.poset.bottomRightCorner(n, n));
}


//-----------------------------------------------------------------------------------------
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
Rcpp::List estimate_lambda_wrapper(Rcpp::IntegerMatrix& ipat,
          Rcpp::IntegerMatrix& iposet, Rcpp::IntegerMatrix& iisF,
          Rcpp::NumericMatrix& ieps, Rcpp::IntegerMatrix& isetD,
          Rcpp::IntegerMatrix& ieventD, Rcpp::NumericVector ilambdaS,
          Rcpp::IntegerVector iL, Rcpp::CharacterVector isampling,
          Rcpp::IntegerVector imaxIter, Rcpp::IntegerVector iupdateStepSize,
          Rcpp::NumericVector itol, Rcpp::NumericVector imaxLambda,
          Rcpp::IntegerVector ineighborhoodDist, Rcpp::IntegerVector ithrds,
          Rcpp::LogicalVector iis_update_eps, Rcpp::IntegerVector iseed) {


  /*      obs poset eps setD    */
  int rowN = 0;
  for (int i = 0; i < iisF.nrow(); ++i) {
    bool is_Fill = true;
    for (int j = 0; j < iisF.ncol(); ++j) {
      if (iisF(i,j) == 0){
        is_Fill = false;
        break;
      }
    }
    if (is_Fill){
      rowN ++;
    } else{
      break;
    }
  }
  int N_total = 0;
  int n = iposet.nrow();
  MyDoubleMatrix eps(rowN,n);
  std::vector<doubleD> setD(rowN);

  std::vector<int> colEvent;
  colEvent.push_back(-1);
  for (int i = 0; i < ieventD.nrow(); ++i) {
    for (int j = ieventD(i,0); j <= ieventD(i,1); ++j) {
      colEvent.push_back(i);
    }
  }

  for (int i = 0; i < rowN; ++i) {
    N_total += isetD(i,1) - isetD(i,0) +1;
    for (int j = 0; j < n; ++j) {
      eps(i,j) = ieps(i,colEvent[j+1]);
      setD[i] = doubleD(isetD(i,0),isetD(i,1));
    }
  }

  MyBoolMatrix obs;
  obs.resize(N_total,n);
  for (int i = 0; i < N_total; ++i) {
    for (int j = 0; j < obs.cols(); ++j) {
      obs(i,j) = ipat(i,j+1);
    }
  }

  MyIntMatrix poset = as<MyIntMatrix>(iposet);

  /*      eventD     */
  std::vector<doubleD> eventD(ieventD.nrow());
  for (int i = 0; i < ieventD.nrow(); ++i) {
    eventD[i] = doubleD(ieventD(i,0),ieventD(i,1));
  }

  /*      L sampling seed thrds lambdaS      */
  unsigned int L = as<unsigned int>(iL);
  std::string sampling = as<std::string>(isampling);
  int seed = as<int>(iseed);
  int thrds = as<int>(ithrds);
  double lambdaS = as<double>(ilambdaS);
  bool is_update_eps = as<bool>(iis_update_eps);

  /*      max_iter update_step_size tol max_lambda neighborhood_dist      */
  int maxIter = as<int>(imaxIter);
  int updateStepSize = as<int>(iupdateStepSize);
  double tol = as<double>(itol);
  double maxLambda = as<double>(imaxLambda);
  unsigned int neighborhoodDist = as<unsigned int>(ineighborhoodDist);

  #ifdef _OPENMP
    omp_set_num_threads(thrds);
    omp_set_dynamic(false);
  #endif

  // cout << "rowN: " << rowN << "N_total: " << N_total << "n: " << n << endl;
  // cout << "obs:" << obs.rows() << "   " << obs.cols() << endl;
  // cout << "poset: " << endl << poset << endl;
  // cout << "eps:" << endl << eps << endl;
  // for (int i = 0; i < setD.size(); ++i) {
  //   cout << setD[i].F << "  " << setD[i].S << endl;
  // }
  // for (int i = 0; i < eventD.size(); ++i) {
  //   cout << eventD[i].F << "  " << eventD[i].S << endl;
  // }
  // cout << "L: " << L << "   sampling: " << sampling << endl;
  // cout << "seed: " << seed << "   thrds: " << thrds << endl;
  // cout << "lambdaS: " << lambdaS << "   max_iter: " << maxIter << endl;
  // cout << "update_step_size: " << updateStepSize << "   tol: " << tol << endl;
  // cout << "max_lambda: " << maxLambda << "   neighborhood_dist: " << neighborhoodDist << endl;

  VectorXd lambdaRes;
  /*      estimate lambda     */
  lambda_estimate(lambdaRes, poset, obs, eps, setD, eventD, L, sampling, is_update_eps, seed,
                  thrds, lambdaS, maxIter, updateStepSize, tol, maxLambda, neighborhoodDist);

  return Rcpp::List::create(Named("epsilon1")=ieps, Named("epsilon2")=eps, Named("lambda")=lambdaRes);

}




// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
Rcpp::IntegerVector isCompatible_wrapper(Rcpp::IntegerVector& igenotype, Rcpp::IntegerMatrix& iposet) {

  std::vector<int> genotype = as<std::vector<int>>(igenotype);
  MyIntMatrix poset = as<MyIntMatrix>(iposet);

  int p = poset.rows();
  for (int i = 0; i < p; ++i) {
    if (genotype[i] == 1) {
      for (int j = 0; j < p; ++j){
        if(poset(j,i) == 1 && genotype[j] == 0){
          return wrap(0);
        }
      }
    }
  }

  return wrap(1);
}







// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
Rcpp::NumericMatrix sample_Age_T_wrapper(Rcpp::NumericMatrix& idata,
                    Rcpp::NumericVector ialpha, Rcpp::NumericVector ibeta,
                    Rcpp::NumericVector ilambda, Rcpp::IntegerVector imaxSampleIter,
                    Rcpp::NumericVector irateE1, Rcpp::IntegerVector ithrds, Rcpp::IntegerVector iseed) {

  MyDoubleMatrix data = as<MyDoubleMatrix>(idata);
  double alpha = as<double>(ialpha);
  double beta = as<double>(ibeta);
  double lambda = as<double>(ilambda);

  int maxSampleIter = as<int>(imaxSampleIter);
  double rateE1 = as<double>(irateE1);

  int N = data.rows();
  MyDoubleMatrix res;
  res.resize(N,2);

  /*      parallel     */
  int thrds = as<int>(ithrds);
  #ifdef _OPENMP
    omp_set_num_threads(thrds);
    omp_set_dynamic(false);
  #endif
  int seed = as<int>(iseed);
  std::mt19937 gen(seed);
  std::gamma_distribution<double> gamma_dist(alpha, 1/beta);
  std::exponential_distribution<double> exp_dist(lambda);
  double z1_proposed = 0;
  double z2_proposed = 0;
  double t1 = 0;
  double t2 = 0;

  #pragma omp parallel for
  for (int i = 0;i < N; i++){
    bool isAccepte = false;
    double a1 = data(i, 0);
    double a2 = data(i, 1);
    double age = data(i, 2);
    double low_B = data(i, 3);
    double up_B = data(i, 4);
    double step_V = (up_B + low_B)*rateE1/2;

    for (int j = 0; j < maxSampleIter; ++j){
      z1_proposed = gamma_dist(gen);
      z2_proposed = exp_dist(gen);
      t1 = z1_proposed + z2_proposed*a1;
      t2 = z1_proposed + z2_proposed*a2;
      if (t1 <= age && age <= t2 && low_B <= t1 && t2 <= up_B){
        isAccepte = true;
        break;
      }
    }

    while (!isAccepte){
      up_B = up_B + step_V;
      low_B = low_B - step_V;
      for (int j = 0; j < maxSampleIter; ++j){
        z1_proposed = gamma_dist(gen);
        z2_proposed = exp_dist(gen);
        t1 = z1_proposed + z2_proposed*a1;
        t2 = z1_proposed + z2_proposed*a2;
        if (t1 <= age && age <= t2 && low_B <= t1 && t2 <= up_B){
          isAccepte = true;
          break;
        }
      }
    }
    res(i,0) = z1_proposed;
    res(i,1) = z2_proposed;
  }

  return wrap(res);
}





// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::export]]
Rcpp::NumericMatrix sample_Age_TLW_wrapper(Rcpp::NumericMatrix& idata,
                                         Rcpp::NumericVector ialpha, Rcpp::NumericVector ibeta,
                                         Rcpp::NumericVector ilambda, Rcpp::IntegerVector ithrds) {

  MyDoubleMatrix data = as<MyDoubleMatrix>(idata);
  double alpha = as<double>(ialpha);
  double beta = as<double>(ibeta);
  double lambda = as<double>(ilambda);

  int N = data.rows();
  MyDoubleMatrix res;
  res.resize(N,2);

  /*      parallel     */
  int thrds = as<int>(ithrds);
  #ifdef _OPENMP
    omp_set_num_threads(thrds);
    omp_set_dynamic(false);
  #endif
  std::random_device rd;
  std::mt19937 gen(rd());
  std::gamma_distribution<double> gamma_dist(alpha, 1/beta);
  std::exponential_distribution<double> exp_dist(lambda);
  double z1_proposed = 0;
  double z2_proposed = 0;
  double t1 = 0;
  double t2 = 0;

  #pragma omp parallel for
  for (int i = 0;i < N; i++){
    bool isAccepte = false;
    double a1 = data(i, 0);
    double a2 = data(i, 1);
    double age = data(i, 2);

    while (!isAccepte){
      z1_proposed = gamma_dist(gen);
      z2_proposed = exp_dist(gen);
      t1 = z1_proposed + z2_proposed*a1;
      t2 = z1_proposed + z2_proposed*a2;
      if (t1 <= age && age <= t2){
        isAccepte = true;
        break;
      }
    }
    res(i,0) = z1_proposed;
    res(i,1) = z2_proposed;
  }

  return wrap(res);
}



