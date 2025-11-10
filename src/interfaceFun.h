#ifndef INTERFACEFUN_H
#define INTERFACEFUN_H

#include "lambdaCal.h"
#include "gmm.h"

using namespace std;
using namespace Eigen;


//------------------------------------------------------------------------------------------------------------------------
/*                  copy_Matrix                       */
void copy_Matrix(int** a,dataD& dD,int N,int row,int l,int* s){
  int shift = dD.setD[row].F;
  for (int i = 1; i < l+1; ++i) {
    int t = s[i-1];
    for (int j = 0; j < N; ++j) {
      a[j][i] = dD.pat(j+shift,t);
    }
  }
}

/*                  model_data_init                       */
void model_data_init(model *M_SA,int** Data_SA,data* D_SA,int N){
  for (int i=0;i<3;i++)
  {
    M_SA[i].n = 2;
    M_SA[i].P = get_int_matrix(3, 3);
  }
  M_SA[0].P[1][2]=1;
  M_SA[1].P[2][1]=1;
  for (int i=0;i<3;i++)
  {
    M_SA[i].lin_ext = get_int_array(2);
    M_SA[i].J_P = bfs_order_ideals(M_SA[i].P, M_SA[i].n+1, &(M_SA[i].m), M_SA[i].lin_ext);
    parents(&M_SA[i]);
    children(&M_SA[i]);
  }
  for (int i = 0; i < N; ++i) {
    Data_SA[i][0] = 1;
  }
  for (int i = 0; i < 4; ++i) {
    D_SA[i].Q = get_int_matrix(3, 3);
    D_SA[i].J_Q = get_int_array(8);
    D_SA[i].g = get_int_array(3);
    D_SA[i].t = get_double_array(3);
  }
}

/*                  estimate_epsilon when n <= 7                       */
double estimate_epsilon_se(dataD& dD,int N,int n,int row,int col, int** GENOTYPE,std::vector<doubleD>& ECount,
                           int T=10,int N_iter=150,double S=1.0,int R = 1){
  // if (n > 7){
  //   cout << "The event numbers should small than 7！" << endl;
  //   return -1;
  // }
  bool isFree = false;
  if (ECount.size() == 0){
    for (int i = dD.eventD[col].F; i <= dD.eventD[col].S; ++i) {
      doubleD temp;
      temp.F = i;
      temp.S = 0;
      for (int j = dD.setD[row].F; j <= dD.setD[row].S; ++j) {
        if (dD.pat(j,i)==1){
          temp.S ++;
        }
      }
      ECount.push_back(temp);
    }
    std::sort(ECount.begin(),ECount.end(),[](doubleD a,doubleD b){
      return a.S > b.S;
    });
    isFree = true;
  }
  double eps = MIN(dD.epsilon(row,col), 1.0);
  double epsilon = 0;
  double total_loglik_i = 0;
  int* copyIndex = get_int_array(n);

  /*                      Initial                      */
  int N_u;
  int** data_e = get_int_matrix(N,n+1);
  for (int i = 0; i < N; ++i) {
    data_e[i][0] = 1;
  }
  for (int i = 0; i < n; ++i) {
    copyIndex[i] = ECount[i].F;
  }
  copy_Matrix(data_e,dD,N,row,n,copyIndex);

  model M_i;
  M_i.n = n;
  M_i.lin_ext = get_int_array(M_i.n);
  M_i.P = get_int_matrix(M_i.n+1,M_i.n+1);
  for (int i = 1; i < n; ++i) {
    M_i.P[i][i+1] = 1;
  }

  double* lambda_i = get_double_array(M_i.n+1);
  lambda_i[0] = S;

  int* pat_idx = get_int_array(N);
  data* D = new_make_data_set(data_e, N, n, &N_u, pat_idx);

  /*  single run: */
  select_poset(eps, &M_i, lambda_i, D, N_u, R, GENOTYPE);
  /* Compute variables */
  double* Prob = get_double_array(M_i.m);
  double** condprob = get_double_matrix(M_i.m, N_u);
  double* lambda_exit = get_double_array(M_i.m);
  compute_lambda_exit(lambda_i, &M_i, lambda_exit);
  int* lattice_index = get_int_array(pow2(M_i.n+1));
  for (int i=0; i < M_i.m; i++)
    lattice_index[M_i.J_P[i]]=i;
  compute_all_prob(lambda_i, &M_i, lambda_exit, Prob, lattice_index,GENOTYPE);

  /* Estimate epsilon | lambda */
  epsilon=0.5; // initial value
  total_loglik_i = EM_epsilon(&M_i, D, N_u, Prob, condprob, &epsilon,GENOTYPE);
  free(Prob);
  free(lambda_exit);
  free(lattice_index);
  free_double_Matrix(condprob,M_i.m);

  epsilon = 0.00001;
  total_loglik_i = EM_EM(&M_i, D, N_u, lambda_i, &epsilon, GENOTYPE);

  /* Local search */
  local_search(&M_i, D, N_u, lambda_i, &epsilon, total_loglik_i, T, N_iter, GENOTYPE);

  std::unordered_set<int> set_temp;
  for (int i = 1; i < n+1; ++i) {
    set_temp.insert(copyIndex[i-1]);
    for (int j = 1; j < n+1; ++j) {
      if (M_i.P[i][j] == 1){
        dD.poset_temp(copyIndex[i-1],copyIndex[j-1]) = 1;
      }
    }
  }
  dD.event_temp.push_back(set_temp);

  epsilon = 0.00001;
  total_loglik_i = EM_EM(&M_i, D, N_u, lambda_i, &epsilon, GENOTYPE);

  /*               need free or clear                 */
  free_data(D,N_u,n);
  free(pat_idx);
  free(lambda_i);
  free_model(&M_i);
  free_int_Matrix(data_e, N);
  free(copyIndex);
  if (isFree){
    std::vector<doubleD>().swap(ECount);
  }

  return epsilon;
}

/*                  estimate_epsilon when n > 7                       */
double estimate_epsilon_be(dataD& dD, int N, int n,int row,int col, int** GENOTYPE, int** GENOTYPE3, int n_p,
                           double threshold = 0.05,int T=10,int N_iter=150, double S=1.0,int R = 1){

  int N_u;
  double eps = MIN(dD.epsilon(row,col), 1.0);
  std::vector<doubleD> ECount;
  double epsilon = 0;
  double total_loglik_i = 0;
  // int n_p = 6;
  int* copyIndex = get_int_array(n_p);
  int* copyIndex1 = get_int_array(2);

  for (int i = dD.eventD[col].F; i <= dD.eventD[col].S; ++i) {
    doubleD temp;
    temp.F = i;
    temp.S = 0;
    for (int j = dD.setD[row].F; j <= dD.setD[row].S; ++j) {
      if (dD.pat(j,i)==1){
        temp.S ++;
      }
    }
    ECount.push_back(temp);
  }
  std::sort(ECount.begin(),ECount.end(),[](doubleD a,doubleD b){
    return a.S > b.S;
  });

  /*                      Calculate the initial epsilon                    */
  epsilon = estimate_epsilon_se(dD,N,n_p,row,col,GENOTYPE,ECount);

  /*                           estimate_epsilon                           */
  int** data_e = get_int_matrix(N,n_p+1);
  for (int i = 0; i < N; ++i) {
    data_e[i][0] = 1;
  }
  double* lambda_i = get_double_array(n_p+1);
  lambda_i[0] = S;
  int* pat_idx = get_int_array(N);

  int count_SA = 20;
  double* epsilon_SA = get_double_array(count_SA);
  int epsilon_count_size = -1;

  model M_SA[3];
  int** Data_SA = get_int_matrix(N,3);
  data* D_SA = (data*)calloc(4, sizeof(data));
  model_data_init(M_SA,Data_SA,D_SA,N);

  std::unordered_set<int> selectE;
  double* lambda = get_double_array(3);
  std::vector<doubleD> PC;
  std::unordered_map<int,int> PI;

  for (int i = 0; i < count_SA; ++i) {
    for (int j = 0; j < n-1; ++j) {
      if (selectE.size() == n_p){
        break;
      }
      /*  Look for events that used to estimate e  */
      for (int k = j+1; k < min((int) (1.3334*n_p),(int) (0.75*n)); ++k) {
        if (selectE.size() == n_p){
          break;
        }
        if (selectE.find(ECount[k].F)==selectE.end()){
          copyIndex1[0] = ECount[j].F;
          copyIndex1[1] = ECount[k].F;
          copy_Matrix(Data_SA,dD,N,row,2,copyIndex1);
          /*  Compute D_SA.count D_SA.g  */
          memset(pat_idx,0,N*sizeof(int));
          N_u = make_data_set(D_SA,Data_SA, N, 2, pat_idx);
          /*  Compute is_compatible  */
          memset(lambda,0.0,3*sizeof(double));
          lambda[0] = S;
          select_poset1(&M_SA[2], lambda, D_SA, N_u, R, GENOTYPE3);
          double e_temp2 = EM_EM1(&M_SA[2], D_SA, N_u, lambda, &epsilon, GENOTYPE3);
          memset(lambda,0.0,3*sizeof(double));
          lambda[0] = S;
          select_poset1(&M_SA[0], lambda, D_SA, N_u, R, GENOTYPE3);
          double e_temp0 = EM_EM1(&M_SA[0], D_SA, N_u, lambda, &epsilon, GENOTYPE3);

          if (e_temp0 > e_temp2){
            if (selectE.size() <= n_p - 2){
              selectE.insert(ECount[j].F);
              selectE.insert(ECount[k].F);
              doubleD temp = {ECount[j].F,ECount[k].F};
              PC.push_back(temp);
              PI.insert({ECount[j].F,-1});
              PI.insert({ECount[k].F,-1});
            }

          } else{
            memset(lambda,0.0,3*sizeof(double));
            lambda[0] = S;
            select_poset1(&M_SA[1], lambda, D_SA, N_u, R, GENOTYPE3);
            e_temp0 = EM_EM1(&M_SA[1], D_SA, N_u, lambda, &epsilon, GENOTYPE3);
            if (e_temp0 > e_temp2){
              if (selectE.size() <= n_p - 2){
                selectE.insert(ECount[j].F);
                selectE.insert(ECount[k].F);
                doubleD temp = {ECount[k].F,ECount[j].F};
                PC.push_back(temp);
                PI.insert({ECount[j].F,-1});
                PI.insert({ECount[k].F,-1});
              }

            }
          }
        }

      }
    }

    if (selectE.size() < n_p){
      for (int j = 0; j < n-1; ++j) {
        if (selectE.size() == n_p){
          break;
        }
        if (selectE.find(ECount[j].F)==selectE.end()){
          selectE.insert(ECount[j].F);
          PI.insert({ECount[j].F,-1});
        }
      }
    }


    N_u = 0;
    for (unordered_set<int>::iterator it = selectE.begin();it != selectE.end();it ++) {
      copyIndex[N_u] = *it;
      PI.find(*it)->second = N_u+1;
      N_u ++;
    }
    copy_Matrix(data_e,dD,N,row,n_p,copyIndex);

    model M_i;
    M_i.n = n_p;
    M_i.lin_ext = get_int_array(M_i.n);
    M_i.P = get_int_matrix(M_i.n+1,M_i.n+1);
    for (int j = 0; j < PC.size(); ++j) {
      M_i.P[PI.find(PC[j].F)->second][PI.find(PC[j].S)->second] = 1;
    }

    memset(pat_idx,0,N*sizeof(int));
    data* D = new_make_data_set(data_e, N, n_p, &N_u, pat_idx);

    /*  single run: */
    select_poset(eps, &M_i, lambda_i, D, N_u, R, GENOTYPE);
    /* Compute variables */
    epsilon = 0.00001;
    total_loglik_i = EM_EM(&M_i, D, N_u, lambda_i, &epsilon, GENOTYPE);
    /* Local search */
    local_search(&M_i, D, N_u, lambda_i, &epsilon, total_loglik_i, T, N_iter, GENOTYPE);


    std::unordered_set<int> set_temp;
    for (int j = 1; j < n_p + 1; ++j) {
      set_temp.insert(copyIndex[j-1]);
      for (int k = 1; k < n_p + 1; ++k) {
        if (M_i.P[j][k] == 1){
          dD.poset_temp(copyIndex[j-1],copyIndex[k-1]) = 1;
        }
      }
    }
    dD.event_temp.push_back(set_temp);



    epsilon = 0.00001;
    total_loglik_i = EM_EM(&M_i, D, N_u, lambda_i, &epsilon, GENOTYPE);
    /*               need free or clear                 */
    free_data(D,N_u,n_p);
    free_model(&M_i);
    memset(lambda_i,0.0,(n_p+1)*sizeof (double));
    lambda_i[0] = S;

    selectE.clear();
    PC.clear();
    PI.clear();

    bool isBreak = FALSE;

    if (isnormal(epsilon)==1){
      if (epsilon_count_size > -1){
        for (int i = 0; i <= epsilon_count_size; ++i) {
          if (abs(epsilon_SA[i]-epsilon)/epsilon_SA[i] <= threshold){
            epsilon = (epsilon+epsilon_SA[i])/2;
            isBreak = TRUE;
            break;
          }
        }
      }
      epsilon_count_size ++;
      epsilon_SA[epsilon_count_size] = epsilon;
    }
    if (isBreak == TRUE){
      break;
    }
  }

  /*  free  */
  free(copyIndex);
  free(copyIndex1);
  std::unordered_map<int,int>().swap(PI);
  std::vector<doubleD>().swap(PC);
  free(lambda);

  std::unordered_set<int>().swap(selectE);
  free_data(D_SA,4,2);
  free_int_Matrix(Data_SA,N);
  for (int i = 0; i < 3; ++i) {
    free_model(&M_SA[i]);
  }
  free(epsilon_SA);
  free(pat_idx);
  free(lambda_i);

  free_int_Matrix(data_e, N);
  std::vector<doubleD>().swap(ECount);

  return epsilon;
}

/*                  estimate_epsilon interface                       */
void estimate_epsilon_mix(dataD& dD, int thrds, std::vector<int**>& GENOTYPE, double threshold = 0.05,
                          int threshold1 = 7, int n_p = 0, int T = 10,int N_iter=150, double S=1.0,int R=1){

  if(n_p <= 0){
    n_p = threshold1;
  }

  for (int i = 0; i < dD.eventNum; ++i) {
    int n = dD.eventD[i].S-dD.eventD[i].F+1;
    if (n>threshold1){
      if (GENOTYPE[n_p + 1] == NULL){
        GENOTYPE[n_p + 1] = precompute_binary(n_p + 1);
      }
    }else{
      if (GENOTYPE[n+1] == NULL){
        GENOTYPE[n+1] = precompute_binary(n+1);
      }
    }
  }
  int thrdsN = dD.eventNum*dD.setNum;

  // #pragma omp parallel for
  for (int t = 0; t < thrdsN; ++t) {
    int i = t / dD.eventNum;
    int j = t % dD.eventNum;
    int n = dD.eventD[j].S-dD.eventD[j].F+1;
    if (dD.isCE(i,j)){
      int N = dD.setD[i].S-dD.setD[i].F+1;
      if (n>threshold1){
        dD.epsilon(i,j) = estimate_epsilon_be(dD,N,n,i,j, GENOTYPE[n_p + 1],GENOTYPE[3], n_p,
                   threshold,T,N_iter,S,R);
      }else{
        std::vector<doubleD> ECount;
        dD.epsilon(i,j) = estimate_epsilon_se(dD,N,n,i,j,GENOTYPE[n+1],ECount, T,N_iter,S,R);
      }
    }
  }


}


void MCEM_hcbn_with_fixed_lambda(Model& model, MatrixXb& obs, unsigned int L, std::string& sampling, ControlEM& control_EM,
                                 const unsigned int thrds, Context& ctx) {

  // Initialization and instantiation of variables
  const vertices_size_type p = model._size; // Number of mutations / events
  int rowN = model._epsilon.rows();
  int eN = model._eventD.size();

  const unsigned int N = obs.rows();         // Number of observations / genotypes
  double N_eff = 0.0;
  unsigned int K = 0;
  unsigned int update_step_size = control_EM.update_step_size;

  MyDoubleMatrix avg_eps(rowN,p);
  MyDoubleMatrix avg_eps_current(rowN,p);

  bool tol_comparison = true;
  MyDoubleMatrix expected_dist = MyDoubleMatrix::Zero(rowN,eN);

  MatrixXd expected_Tdiff = MatrixXd::Zero(N, p);
  VectorXd Tdiff_colsum(p);
  MatrixXd Tdiff_pool;
  VectorXd scale_cumulative;


  if (sampling == "add-remove") {
    scale_cumulative.resize(p);
    if (model.get_update_node_idx())
      model.update_node_idx();
  } else if (sampling == "pool") {
    K = p * L;
    Tdiff_pool.resize(K, p);

  }

  for (unsigned int iter = 0; iter < control_EM.max_iter; ++iter) {
    MatrixXd T_pool;
    /*    update eps & lambda    */
    if (iter == update_step_size) {
      avg_eps_current /= control_EM.update_step_size;
      if (tol_comparison) {
        if (matrixDivision1(avg_eps, avg_eps_current, control_EM.tol))
          break;
        // L *= 2;
      }
      avg_eps = avg_eps_current;

      update_step_size += control_EM.update_step_size;
      tol_comparison = true;

      /* Restart averaging */
      avg_eps_current.fill(0.0);
    }

    /* E step
     * Conditional expectation for the sufficient statistics per observation and event
     * scale_cumulative: cumulative E(lambda)
     */
    if (sampling == "add-remove") {
      scale_cumulative = scale_path_to_mutation(model);
    } else if (sampling == "pool") {
      /* All threads share the same pool of mutation times */
      T_pool.resize(K, p);
      /* Generate observation times from a given poset and given rates */
      T_pool = sample_times(K, model, Tdiff_pool, ctx.rng);
    }

    expected_dist.fill(0.0);

    /*     parallel    */
    //        #ifdef _OPENMP
    //            omp_set_num_threads(thrds);
    //        #endif
    /*     generate thrds rng     */
    auto rngs = ctx.get_auxiliary_rngs(thrds);
    #pragma omp declare reduction(matrix_add: MyDoubleMatrix:omp_out=add_matrix(omp_out,omp_in)) initializer(omp_priv = omp_orig)
    {
      #pragma omp parallel for reduction(matrix_add:expected_dist) schedule(static)
      for (unsigned int i = 0; i < N; ++i) {
        int rowI = 0;
        for (int j = 0; j < model._setD.size(); ++j) {
          if (i>=model._setD[j].F && i<=model._setD[j].S){
            rowI = j;
            break;
          }
        }

        MyIntMatrix d_pool;
        if (sampling == "pool") {
          VectorXd T_sampling(K);
          /*    generate genotype    */
          MatrixXb genotype_pool =
          generate_genotypes(T_pool, model, T_sampling, (*rngs)[omp_get_thread_num()]);
          d_pool = hamming_dist_mat1(genotype_pool, obs.row(i), model);
        }

        DataImportanceSampling importance_sampling = importance_weight(
          obs.row(i), L, model, rowI, sampling, scale_cumulative, d_pool, Tdiff_pool,
          control_EM.neighborhood_dist, (*rngs)[omp_get_thread_num()]);
        double aux = importance_sampling.w.sum();
        if (aux > 0) {
          /* Only consider observations with at least one feasible sample */
          //                    int L_eff = L;
          /*                    if (sampling == "backward")
           L_eff = (importance_sampling.w.array() > 0).count();*/
          for (int j = 0; j < eN; ++j) {
            expected_dist(rowI,j) += importance_sampling.w.dot(importance_sampling.dist.col(j).cast<double>()) / aux;
          }
          expected_Tdiff.row(i) = (importance_sampling.Tdiff.transpose() * importance_sampling.w) / aux;
        }
      }
    }

    /* M-step */
    for (int i = 0; i < rowN; ++i) {
      for (int j = 0; j < eN; ++j) {
        expected_dist(i,j) = expected_dist(i,j)/((model._setD[i].S-model._setD[i].F+1)*(model._eventD[j].S-model._eventD[j].F+1));
      }
    }

    MyDoubleMatrix epsTemp(rowN,p);

    for (int i = 0; i < eN; ++i) {
      for (int j = model._eventD[i].F; j <= model._eventD[i].S; ++j) {
        for (int k = 0; k < rowN; ++k) {
          epsTemp(k,j-1) = expected_dist(k,i);
        }
      }
    }
    model._epsilon = epsTemp;

    RowVectorXd weights(N);
    weights.fill(1.0);
    Tdiff_colsum = weights*expected_Tdiff;
    avg_eps_current += model._epsilon;
    if (iter + 1 == control_EM.max_iter) {
      unsigned int num_iter = control_EM.max_iter - update_step_size +
        control_EM.update_step_size;
      avg_eps_current /= num_iter;
    }

  }

  model._epsilon = avg_eps_current;

}


/*                  estimate epsilon based on poset and lambda                       */
void estimate_epsilon_based_on_poset_and_lambda(VectorXd& lambdaRes, MyIntMatrix& poset, MyBoolMatrix& obs, MyDoubleMatrix& eps,
                                                std::vector<doubleD>& setD,
                     std::vector<doubleD>& eventD, unsigned int L, std::string& sampling, int seed,
                     int thrds, double lambda_s=1.0, int max_iter=100L, int update_step_size=20L,
                     double tol=0.001, double max_lambda=1e6, unsigned int neighborhood_dist=1){

  const auto p = poset.rows(); // Number of mutations / events
  VectorXd ilambda(p+1);
  ilambda(0) = lambda_s;
  for (int i = 1; i < p + 1; ++i) {
    ilambda(i) = lambdaRes(i-1);
  }

  VectorXd times;
  times.resize(obs.rows());
  bool sampling_times_available = false;

  edge_container edge_list = adjacency_mat2list(poset);
  Model M(edge_list, p, lambda_s);
  M.set_lambda(ilambda, max_lambda);
  M._epsilon = eps;
  M._setD = setD;
  M._eventD = eventD;
  M.has_cycles();
  if (M.cycle){
    cout << "Poset has cycle, please remove some edges for the poset and re-enter it!";
    return;
  }
  M.topological_sort();

  ControlEM control_EM(max_iter, update_step_size, tol, max_lambda, neighborhood_dist);

  /* Call the underlying C++ function */
  Context ctx(seed);

  MCEM_hcbn_with_fixed_lambda(M, obs, L, sampling, control_EM, thrds, ctx);

  lambdaRes = M._lambda;
  eps = M._epsilon;
}


//------------------------------------------------------------------------------------------------------------------------
// /*                  model_data_init                       */
void model_data_init1(model *M_SA){
  for (int i=0;i<3;i++)
  {
    M_SA[i].n = 2;
    M_SA[i].P = get_int_matrix(3, 3);
  }
  M_SA[0].P[1][2]=1;
  M_SA[1].P[2][1]=1;
  for (int i=0;i<3;i++)
  {
    M_SA[i].lin_ext = get_int_array(2);
    M_SA[i].J_P = bfs_order_ideals(M_SA[i].P, M_SA[i].n+1, &(M_SA[i].m), M_SA[i].lin_ext);
    parents(&M_SA[i]);
    children(&M_SA[i]);
  }
}

// /*                  partial ordering relation judgment                       */
std::vector<double> find_poset(dataD& dD, int col1,int col2, int** GENOTYPE3, double threshold = 0.0, double S=1.0,int R = 1){

  int rowN = 0;
  int N_total = 0;
  for (int i = 0; i < dD.setNum; ++i) {
    if (dD.isF(i,dD.colEvent[col1])== true && dD.isF(i,dD.colEvent[col2])== true){
      rowN ++;
    }else{
      break;
    }
  }

  std::vector<int> N(rowN);
  std::vector<int> N_u(rowN);
  std::vector<double> epsilon1(rowN);
  std::vector<double> epsilon2(rowN);
  std::vector<double*> lambda(rowN);
  std::vector<int**> Data_SA(rowN);
  std::vector<data1*> D_SA(rowN);

  double* union_lambda = get_double_array(3);
  union_lambda[0] = S;
  int* copyIndex = get_int_array(2);
  copyIndex[0] = col1;
  copyIndex[1] = col2;

  for (int i = 0; i < rowN; ++i) {
    int N_temp = dD.setD[i].S - dD.setD[i].F +1;
    N_total += N_temp;
    N[i]=N_temp;
    epsilon1[i]=dD.epsilon(i,dD.colEvent[col1]);
    epsilon2[i]=dD.epsilon(i,dD.colEvent[col2]);
    double* lambda_temp = get_double_array(3);
    lambda_temp[0] = S;
    lambda[i]=lambda_temp;

    int** Data_SA_temp = get_int_matrix(N_temp,3);
    for (int i = 0; i < N_temp; ++i) {
      Data_SA_temp[i][0] = 1;
    }
    copy_Matrix(Data_SA_temp,dD,N_temp,i,2,copyIndex);
    Data_SA[i]=Data_SA_temp;

    data1* D_SA_temp = (data1*)calloc(4, sizeof(data1));
    for (int j = 0; j < 4; ++j) {
      D_SA_temp[j].Q = get_int_matrix(3, 3);
      D_SA_temp[j].J_Q = get_int_array(8);
      D_SA_temp[j].g = get_int_array(3);
      D_SA_temp[j].t = get_double_array(3);
      D_SA_temp[j].compatible_state = get_int_array(8);
    }
    D_SA[i]=D_SA_temp;
  }

  model M_SA[3];
  model_data_init1(M_SA);

  /*  Compute D_SA.count D_SA.g  */
  for (int i = 0; i < rowN; ++i) {
    N_u[i] = make_data_set1(D_SA[i],Data_SA[i], N[i], 2);
  }

  /*  Compute is_compatible & */
  for (int i = 0; i < rowN; ++i) {
    select_poset2(&M_SA[2], lambda[i], D_SA[i], N_u[i], R, GENOTYPE3);
  }
  for (int i = 1; i < 3; ++i) {
    for (int j = 0; j < rowN; ++j) {
      union_lambda[i] += lambda[j][i]*N[j];
    }
    union_lambda[i] = union_lambda[i]/N_total;
  }
  double e_temp2 = EM_EM2(&M_SA[2], D_SA, N, N_total, N_u,  lambda, union_lambda,
                          epsilon1, epsilon2,rowN, GENOTYPE3);


  for (int i = 0; i < rowN; ++i) {
    memset(lambda[i],0.0,3*sizeof(double));
    lambda[i][0] = S;
    select_poset2(&M_SA[0], lambda[i], D_SA[i], N_u[i], R, GENOTYPE3);
  }
  memset(union_lambda,0.0,3*sizeof(double));
  union_lambda[0] = S;
  for (int i = 1; i < 3; ++i) {
    for (int j = 0; j < rowN; ++j) {
      union_lambda[i] += lambda[j][i]*N[j];
    }
    union_lambda[i] = union_lambda[i]/N_total;
  }
  double e_temp0 = EM_EM2(&M_SA[0], D_SA, N, N_total, N_u,  lambda, union_lambda,
                          epsilon1, epsilon2, rowN, GENOTYPE3);



  for (int i = 0; i < rowN; ++i) {
    memset(lambda[i],0.0,3*sizeof(double));
    lambda[i][0] = S;
    select_poset2(&M_SA[1], lambda[i], D_SA[i], N_u[i], R, GENOTYPE3);
  }
  memset(union_lambda,0.0,3*sizeof(double));
  union_lambda[0] = S;
  for (int i = 1; i < 3; ++i) {
    for (int j = 0; j < rowN; ++j) {
      union_lambda[i] += lambda[j][i]*N[j];
    }
    union_lambda[i] = union_lambda[i]/N_total;
  }
  double e_temp1 = EM_EM2(&M_SA[1], D_SA, N, N_total, N_u,  lambda, union_lambda,
                          epsilon1, epsilon2, rowN, GENOTYPE3);


  std::vector<double> res;
  res.push_back(e_temp0);
  res.push_back(e_temp1);
  res.push_back(e_temp2);


  /*         free         */
  free(copyIndex);
  free(union_lambda);
  for (int i = 0; i < rowN; ++i) {
    free_data1(D_SA[i],4,2);
    free_int_Matrix(Data_SA[i],N[i]);
    free(lambda[i]);
  }
  std::vector<data1*>().swap(D_SA);
  std::vector<int**>().swap(Data_SA);
  std::vector<double*>().swap(lambda);
  std::vector<double>().swap(epsilon2);
  std::vector<double>().swap(epsilon1);
  std::vector<int>().swap(N_u);
  std::vector<int>().swap(N);

  return res;
}

// /*                  Determine whether the partial order relationship between events i and j has been fine-tuned                      */
bool is_in_SA_mix(std::vector<std::unordered_set<int>>& event_temp, int col1, int col2){
  for (int i = 0; i < event_temp.size(); ++i) {
    if (event_temp[i].find(col1)!=event_temp[i].end() && event_temp[i].find(col2)!=event_temp[i].end()){
      return true;
    }
  }
  return false;
}

// /*                  Determine whether the partial order relationship between events i and j has been fine-tuned                       */
bool is_in_SA(dataD& dD,int col1,int col2){
  return is_in_SA_mix(dD.event_temp, col1,col2);
}

// /*                  poset fine-tune                       */
// void poset_SA(dataD& dD,std::unordered_set<int>& event_Use, MyIntMatrix& poset_temp, std::vector<std::unordered_set<int>>& event_temp,
//               std::vector<int**>& GENOTYPE, double T=10,int N_iter= 150,double S=1.0,int R = 1){
//
//   int n = event_Use.size();
//   if (GENOTYPE[n+1] == NULL){
//     GENOTYPE[n+1] = precompute_binary(n+1);
//   }
//   std::vector<int> event;
//   for (auto i = event_Use.begin(); i != event_Use.end(); ++i) {
//     event.push_back(*i);
//   }
//   event.shrink_to_fit();
//   std::sort(event.begin(),event.end());
//
//   int rowN = 0;
//   int N_total = 0;
//   for (int i = 0; i < dD.setNum; ++i) {
//     bool is_Fill = true;
//     for (int j = 0; j < n; ++j) {
//       if (dD.isF(i,dD.colEvent[event[j]])== false){
//         is_Fill = false;
//         break;
//       }
//     }
//     if (is_Fill){
//       rowN ++;
//     } else{
//       break;
//     }
//   }
//
//   //------------------------------------------------------------------------------------------------------------------
//   std::vector<int> N(rowN);
//   std::vector<int> N_u(rowN);
//   MyDoubleMatrix epsilon;
//   epsilon.resize(rowN,n);
//   int temp0 = dD.colEvent[event[0]];
//   int temp1 = 0;
//   std::vector<std::vector<int>> eventCol;
//   eventCol.push_back(std::vector<int>());
//   for (int i = 0; i < n; ++i) {
//     if (dD.colEvent[event[i]] == temp0){
//       eventCol[temp1].push_back(i);
//     } else{
//       eventCol[temp1].shrink_to_fit();
//       std::vector<int> vector_temp;
//       vector_temp.push_back(i);
//       eventCol.push_back(vector_temp);
//       temp0 = dD.colEvent[event[i]];
//       temp1 ++;
//     }
//   }
//   eventCol.shrink_to_fit();
//   MyDoubleMatrix union_epsilon;
//   union_epsilon.resize(rowN,eventCol.size());
//   for (int i = 0; i < rowN; ++i) {
//     for (int j = 0; j < eventCol.size(); ++j) {
//       union_epsilon(i,j) = dD.epsilon(i, dD.colEvent[event[eventCol[j][0]]]);
//     }
//   }
//
//   std::vector<double*> lambda(rowN);
//   std::vector<int**> Data_SA(rowN);
//   std::vector<data2*> D_SA(rowN);
//   double* union_lambda = get_double_array(n+1);
//   union_lambda[0] = S;
//   int* copyIndex = get_int_array(n);
//   for (int i = 0; i < n; ++i) {
//     copyIndex[i] = event[i];
//   }
//   for (int i = 0; i < rowN; ++i) {
//     int N_temp = dD.setD[i].S - dD.setD[i].F +1;
//     N_total += N_temp;
//     N[i]=N_temp;
//     for (int j = 0; j < n; ++j) {
//       epsilon(i,j) = dD.epsilon(i,dD.colEvent[event[j]]);
//     }
//
//     double* lambda_temp = get_double_array(n+1);
//     lambda_temp[0] = S;
//     lambda[i]=lambda_temp;
//
//     int** Data_SA_temp = get_int_matrix(N_temp,n+1);
//     for (int i = 0; i < N_temp; ++i) {
//       Data_SA_temp[i][0] = 1;
//     }
//     copy_Matrix(Data_SA_temp,dD,N_temp,i,n,copyIndex);
//     Data_SA[i] = Data_SA_temp;
//
//     int* pat_idx = get_int_array(N_temp);
//     data2* D_SA_temp = new_make_data_set1(Data_SA_temp, N_temp, n, &N_u[i], pat_idx);
//     D_SA[i] = D_SA_temp;
//
//     free(pat_idx);
//   }
//
//   double total_loglik_i = 0;
//
//   //------------------------------------------------------------------------------------------------------------------
//   /*    Initial Model    */
//   model M_i;
//   M_i.n = n;
//   M_i.lin_ext = get_int_array(M_i.n);
//   M_i.P = get_int_matrix(M_i.n+1,M_i.n+1);
//   for (int i = 0; i < n; ++i) {
//     for (int j = 0; j < n; ++j) {
//       if (dD.poset(event[i],event[j]) == 1){
//         M_i.P[i+1][j+1] = 1;
//       }
//     }
//   }
//
//   /*  single run: */
//   select_poset3(&M_i, lambda, D_SA, N_u, R, rowN, GENOTYPE[n+1]);
//
//   /*   update union_lambda   */
//   for (int i = 1; i <= n; ++i) {
//     union_lambda[i] = 0;
//     for (int j = 0; j < rowN; ++j) {
//       union_lambda[i] += lambda[j][i]*N[j];
//     }
//     union_lambda[i] = union_lambda[i]/N_total;
//   }
//
//   /* Compute variables */
//   double* Prob = get_double_array(M_i.m);
//   std::vector<double**> condprob(rowN);
//   for (int i = 0; i < rowN; ++i) {
//     condprob[i] = get_double_matrix(M_i.m, N_u[i]);
//   }
//   double* lambda_exit = get_double_array(M_i.m);
//   compute_lambda_exit(union_lambda, &M_i, lambda_exit);
//   int* lattice_index = get_int_array(pow2(M_i.n+1));
//   for (int i=0; i < M_i.m; i++)
//     lattice_index[M_i.J_P[i]]=i;
//   compute_all_prob(union_lambda, &M_i, lambda_exit, Prob, lattice_index,GENOTYPE[n+1]);
//
//   /*   Estimate epsilon   */
//   total_loglik_i = EM_epsilon4(&M_i, D_SA, N_u, Prob, condprob, epsilon, union_epsilon, eventCol, rowN, GENOTYPE[n+1]);
//
//   /*    free variables    */
//   free(lattice_index);
//   free(lambda_exit);
//   for (int i = 0; i < condprob.size(); ++i) {
//     free_double_Matrix(condprob[i],M_i.m);
//   }
//   std::vector<double**>().swap(condprob);
//   free(Prob);
//
//   //------------------------------------------------------------------------------------------------------------------
//   total_loglik_i = EM_EM4(&M_i, D_SA, N, N_total, N_u, lambda, union_lambda, epsilon, union_epsilon, eventCol, rowN, GENOTYPE[n+1]);
//
//
//   /*    cout << "----------------------" << endl;
//    for (int i = 0; i < rowN; ++i) {
//    for (int j = 0; j < n+1; ++j) {
//    cout << lambda[i][j] << "  ";
//    lambda[i][0] = S;
//    }
//    cout << endl;
//    }
//    for (int i = 0; i < n+1; ++i) {
//    cout << union_lambda[i] << "  ";
//    }
//    union_lambda[0] = S;
//    cout << endl;*/
//   /* Local search */
//   int m_temp = local_search3(&M_i, D_SA, N, N_total, N_u, lambda, union_lambda, epsilon, union_epsilon, eventCol, total_loglik_i, rowN, T, N_iter,
//                              GENOTYPE[n+1]);
//
//   std::unordered_set<int> set_temp;
//   for (int i = 1; i < n+1; ++i) {
//     set_temp.insert(copyIndex[i-1]);
//     for (int j = 1; j < n+1; ++j) {
//       if (M_i.P[i][j] == 1){
//         poset_temp(copyIndex[i-1],copyIndex[j-1]) = 1;
//       }
//     }
//   }
//   event_temp.push_back(set_temp);
//
//   M_i.lin_ext = get_int_array(M_i.n);  // a linear extension of the poset
//   M_i.J_P = bfs_order_ideals(M_i.P, M_i.n+1, &(M_i.m), M_i.lin_ext);
//
//   /*         free         */
//   free(copyIndex);
//   free(union_lambda);
//   for (int i = 0; i < rowN; ++i) {
//     free_data2(D_SA[i],N_u[i],n,m_temp);
//     free_int_Matrix(Data_SA[i],N[i]);
//     free(lambda[i]);
//   }
//   std::vector<data2*>().swap(D_SA);
//   std::vector<int**>().swap(Data_SA);
//   std::vector<double*>().swap(lambda);
//   epsilon.resize(0,0);
//   std::vector<int>().swap(N_u);
//   std::vector<int>().swap(N);
//   std::vector<int>().swap(event);
// }

// /*                  poset fine-tune                       */
void poset_SA(dataD& dD,std::unordered_set<int>& event_Use, MyIntMatrix& poset_temp, std::unordered_set<int>& event_temp,
              std::vector<int**>& GENOTYPE, double T=10,int N_iter= 150,double S=1.0,int R = 1){

  int n = event_Use.size();

  // Must be a synchronization code block.
  global_mutex_GENOTYPE.lock();
  if (GENOTYPE[n+1] == NULL){
    GENOTYPE[n+1] = precompute_binary(n+1);
  }
  global_mutex_GENOTYPE.unlock();

  std::vector<int> event;
  for (auto i = event_Use.begin(); i != event_Use.end(); ++i) {
    event.push_back(*i);
  }
  event.shrink_to_fit();
  std::sort(event.begin(),event.end());

  int rowN = 0;
  int N_total = 0;
  for (int i = 0; i < dD.setNum; ++i) {
    bool is_Fill = true;
    for (int j = 0; j < n; ++j) {
      if (dD.isF(i,dD.colEvent[event[j]])== false){
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

  //------------------------------------------------------------------------------------------------------------------
  std::vector<int> N(rowN);
  std::vector<int> N_u(rowN);
  MyDoubleMatrix epsilon;
  epsilon.resize(rowN,n);
  int temp0 = dD.colEvent[event[0]];
  int temp1 = 0;
  std::vector<std::vector<int>> eventCol;
  eventCol.push_back(std::vector<int>());
  for (int i = 0; i < n; ++i) {
    if (dD.colEvent[event[i]] == temp0){
      eventCol[temp1].push_back(i);
    } else{
      eventCol[temp1].shrink_to_fit();
      std::vector<int> vector_temp;
      vector_temp.push_back(i);
      eventCol.push_back(vector_temp);
      temp0 = dD.colEvent[event[i]];
      temp1 ++;
    }
  }
  eventCol.shrink_to_fit();
  MyDoubleMatrix union_epsilon;
  union_epsilon.resize(rowN,eventCol.size());
  for (int i = 0; i < rowN; ++i) {
    for (int j = 0; j < eventCol.size(); ++j) {
      union_epsilon(i,j) = dD.epsilon(i, dD.colEvent[event[eventCol[j][0]]]);
    }
  }

  std::vector<double*> lambda(rowN);
  std::vector<int**> Data_SA(rowN);
  std::vector<data2*> D_SA(rowN);
  double* union_lambda = get_double_array(n+1);
  union_lambda[0] = S;
  int* copyIndex = get_int_array(n);
  for (int i = 0; i < n; ++i) {
    copyIndex[i] = event[i];
  }
  for (int i = 0; i < rowN; ++i) {
    int N_temp = dD.setD[i].S - dD.setD[i].F +1;
    N_total += N_temp;
    N[i]=N_temp;
    for (int j = 0; j < n; ++j) {
      epsilon(i,j) = dD.epsilon(i,dD.colEvent[event[j]]);
    }

    double* lambda_temp = get_double_array(n+1);
    lambda_temp[0] = S;
    lambda[i]=lambda_temp;

    int** Data_SA_temp = get_int_matrix(N_temp,n+1);
    for (int i = 0; i < N_temp; ++i) {
      Data_SA_temp[i][0] = 1;
    }
    copy_Matrix(Data_SA_temp,dD,N_temp,i,n,copyIndex);
    Data_SA[i] = Data_SA_temp;

    int* pat_idx = get_int_array(N_temp);
    data2* D_SA_temp = new_make_data_set1(Data_SA_temp, N_temp, n, &N_u[i], pat_idx);
    D_SA[i] = D_SA_temp;

    free(pat_idx);
  }

  double total_loglik_i = 0;

  //------------------------------------------------------------------------------------------------------------------
  /*    Initial Model    */
  model M_i;
  M_i.n = n;
  M_i.lin_ext = get_int_array(M_i.n);
  M_i.P = get_int_matrix(M_i.n+1,M_i.n+1);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (dD.poset(event[i],event[j]) == 1){
        M_i.P[i+1][j+1] = 1;
      }
    }
  }

  /*  single run: */
  select_poset3(&M_i, lambda, D_SA, N_u, R, rowN, GENOTYPE[n+1]);

  /*   update union_lambda   */
  for (int i = 1; i <= n; ++i) {
    union_lambda[i] = 0;
    for (int j = 0; j < rowN; ++j) {
      union_lambda[i] += lambda[j][i]*N[j];
    }
    union_lambda[i] = union_lambda[i]/N_total;
  }

  /* Compute variables */
  double* Prob = get_double_array(M_i.m);
  std::vector<double**> condprob(rowN);
  for (int i = 0; i < rowN; ++i) {
    condprob[i] = get_double_matrix(M_i.m, N_u[i]);
  }
  double* lambda_exit = get_double_array(M_i.m);
  compute_lambda_exit(union_lambda, &M_i, lambda_exit);
  int* lattice_index = get_int_array(pow2(M_i.n+1));
  for (int i=0; i < M_i.m; i++)
    lattice_index[M_i.J_P[i]]=i;
  compute_all_prob(union_lambda, &M_i, lambda_exit, Prob, lattice_index,GENOTYPE[n+1]);

  /*   Estimate epsilon   */
  total_loglik_i = EM_epsilon4(&M_i, D_SA, N_u, Prob, condprob, epsilon, union_epsilon, eventCol, rowN, GENOTYPE[n+1]);

  /*    free variables    */
  free(lattice_index);
  free(lambda_exit);
  for (int i = 0; i < condprob.size(); ++i) {
    free_double_Matrix(condprob[i],M_i.m);
  }
  std::vector<double**>().swap(condprob);
  free(Prob);

  //------------------------------------------------------------------------------------------------------------------
  total_loglik_i = EM_EM4(&M_i, D_SA, N, N_total, N_u, lambda, union_lambda, epsilon, union_epsilon, eventCol, rowN, GENOTYPE[n+1]);


  /*    cout << "----------------------" << endl;
   for (int i = 0; i < rowN; ++i) {
   for (int j = 0; j < n+1; ++j) {
   cout << lambda[i][j] << "  ";
   lambda[i][0] = S;
   }
   cout << endl;
   }
   for (int i = 0; i < n+1; ++i) {
   cout << union_lambda[i] << "  ";
   }
   union_lambda[0] = S;
   cout << endl;*/
  /* Local search */
  int m_temp = local_search3(&M_i, D_SA, N, N_total, N_u, lambda, union_lambda, epsilon, union_epsilon, eventCol, total_loglik_i, rowN, T, N_iter,
                             GENOTYPE[n+1]);

  std::unordered_set<int> set_temp;
  for (int i = 1; i < n+1; ++i) {
    set_temp.insert(copyIndex[i-1]);
    for (int j = 1; j < n+1; ++j) {
      if (M_i.P[i][j] == 1){
        poset_temp(copyIndex[i-1],copyIndex[j-1]) = 1;
      }
    }
  }

  event_temp = set_temp;

  M_i.lin_ext = get_int_array(M_i.n);  // a linear extension of the poset
  M_i.J_P = bfs_order_ideals(M_i.P, M_i.n+1, &(M_i.m), M_i.lin_ext);

  /*         free         */
  free(copyIndex);
  free(union_lambda);
  for (int i = 0; i < rowN; ++i) {
    free_data2(D_SA[i],N_u[i],n,m_temp);
    free_int_Matrix(Data_SA[i],N[i]);
    free(lambda[i]);
  }
  std::vector<data2*>().swap(D_SA);
  std::vector<int**>().swap(Data_SA);
  std::vector<double*>().swap(lambda);
  epsilon.resize(0,0);
  std::vector<int>().swap(N_u);
  std::vector<int>().swap(N);
  std::vector<int>().swap(event);

}



// /*                  poset fine-tune with fixed eps                      */
void poset_SA_with_fixed_eps(dataD& dD,std::unordered_set<int>& event_Use, MyIntMatrix& poset_temp, std::unordered_set<int>& event_temp,
              std::vector<int**>& GENOTYPE, double T=10,int N_iter= 150,double S=1.0,int R = 1){

  int n = event_Use.size();

  // Must be a synchronization code block.
  global_mutex_GENOTYPE.lock();
  if (GENOTYPE[n+1] == NULL){
    GENOTYPE[n+1] = precompute_binary(n+1);
  }
  global_mutex_GENOTYPE.unlock();

  std::vector<int> event;
  for (auto i = event_Use.begin(); i != event_Use.end(); ++i) {
    event.push_back(*i);
  }
  event.shrink_to_fit();
  std::sort(event.begin(),event.end());

  int rowN = 0;
  int N_total = 0;
  for (int i = 0; i < dD.setNum; ++i) {
    bool is_Fill = true;
    for (int j = 0; j < n; ++j) {
      if (dD.isF(i,dD.colEvent[event[j]])== false){
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

  //------------------------------------------------------------------------------------------------------------------
  std::vector<int> N(rowN);
  std::vector<int> N_u(rowN);
  MyDoubleMatrix epsilon;
  epsilon.resize(rowN,n);
  int temp0 = dD.colEvent[event[0]];
  int temp1 = 0;
  std::vector<std::vector<int>> eventCol;
  eventCol.push_back(std::vector<int>());
  for (int i = 0; i < n; ++i) {
    if (dD.colEvent[event[i]] == temp0){
      eventCol[temp1].push_back(i);
    } else{
      eventCol[temp1].shrink_to_fit();
      std::vector<int> vector_temp;
      vector_temp.push_back(i);
      eventCol.push_back(vector_temp);
      temp0 = dD.colEvent[event[i]];
      temp1 ++;
    }
  }
  eventCol.shrink_to_fit();
  MyDoubleMatrix union_epsilon;
  union_epsilon.resize(rowN,eventCol.size());
  for (int i = 0; i < rowN; ++i) {
    for (int j = 0; j < eventCol.size(); ++j) {
      union_epsilon(i,j) = dD.epsilon(i, dD.colEvent[event[eventCol[j][0]]]);
    }
  }

  std::vector<double*> lambda(rowN);
  std::vector<int**> Data_SA(rowN);
  std::vector<data2*> D_SA(rowN);
  double* union_lambda = get_double_array(n+1);
  union_lambda[0] = S;
  int* copyIndex = get_int_array(n);
  for (int i = 0; i < n; ++i) {
    copyIndex[i] = event[i];
  }
  for (int i = 0; i < rowN; ++i) {
    int N_temp = dD.setD[i].S - dD.setD[i].F +1;
    N_total += N_temp;
    N[i]=N_temp;
    for (int j = 0; j < n; ++j) {
      epsilon(i,j) = dD.epsilon(i,dD.colEvent[event[j]]);
    }

    double* lambda_temp = get_double_array(n+1);
    lambda_temp[0] = S;
    lambda[i]=lambda_temp;

    int** Data_SA_temp = get_int_matrix(N_temp,n+1);
    for (int i = 0; i < N_temp; ++i) {
      Data_SA_temp[i][0] = 1;
    }
    copy_Matrix(Data_SA_temp,dD,N_temp,i,n,copyIndex);
    Data_SA[i] = Data_SA_temp;

    int* pat_idx = get_int_array(N_temp);
    data2* D_SA_temp = new_make_data_set1(Data_SA_temp, N_temp, n, &N_u[i], pat_idx);
    D_SA[i] = D_SA_temp;

    free(pat_idx);
  }

  double total_loglik_i = 0;

  //------------------------------------------------------------------------------------------------------------------
  /*    Initial Model    */
  model M_i;
  M_i.n = n;
  M_i.lin_ext = get_int_array(M_i.n);
  M_i.P = get_int_matrix(M_i.n+1,M_i.n+1);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (dD.poset(event[i],event[j]) == 1){
        M_i.P[i+1][j+1] = 1;
      }
    }
  }

  /*  single run: */
  select_poset3(&M_i, lambda, D_SA, N_u, R, rowN, GENOTYPE[n+1]);

  /*   update union_lambda   */
  for (int i = 1; i <= n; ++i) {
    union_lambda[i] = 0;
    for (int j = 0; j < rowN; ++j) {
      union_lambda[i] += lambda[j][i]*N[j];
    }
    union_lambda[i] = union_lambda[i]/N_total;
  }

  /* Compute variables */
  double* Prob = get_double_array(M_i.m);
  std::vector<double**> condprob(rowN);
  for (int i = 0; i < rowN; ++i) {
    condprob[i] = get_double_matrix(M_i.m, N_u[i]);
  }
  double* lambda_exit = get_double_array(M_i.m);
  compute_lambda_exit(union_lambda, &M_i, lambda_exit);
  int* lattice_index = get_int_array(pow2(M_i.n+1));
  for (int i=0; i < M_i.m; i++)
    lattice_index[M_i.J_P[i]]=i;
  compute_all_prob(union_lambda, &M_i, lambda_exit, Prob, lattice_index,GENOTYPE[n+1]);

  /*   Estimate epsilon   */
  // total_loglik_i = EM_epsilon4(&M_i, D_SA, N_u, Prob, condprob, epsilon, union_epsilon, eventCol, rowN, GENOTYPE[n+1]);

  /*    free variables    */
  free(lattice_index);
  free(lambda_exit);
  for (int i = 0; i < condprob.size(); ++i) {
    free_double_Matrix(condprob[i],M_i.m);
  }
  std::vector<double**>().swap(condprob);
  free(Prob);

  //------------------------------------------------------------------------------------------------------------------
  total_loglik_i = EM_EM5(&M_i, D_SA, N, N_total, N_u, lambda, union_lambda, epsilon, union_epsilon, eventCol, rowN, GENOTYPE[n+1]);


  /*    cout << "----------------------" << endl;
   for (int i = 0; i < rowN; ++i) {
   for (int j = 0; j < n+1; ++j) {
   cout << lambda[i][j] << "  ";
   lambda[i][0] = S;
   }
   cout << endl;
   }
   for (int i = 0; i < n+1; ++i) {
   cout << union_lambda[i] << "  ";
   }
   union_lambda[0] = S;
   cout << endl;*/
  /* Local search */
  int m_temp = local_search4(&M_i, D_SA, N, N_total, N_u, lambda, union_lambda, epsilon, union_epsilon, eventCol, total_loglik_i, rowN, T, N_iter,
                             GENOTYPE[n+1]);

  // std::cout << "union_epsilon:" << std::endl << union_epsilon << std::endl;

  std::unordered_set<int> set_temp;
  for (int i = 1; i < n+1; ++i) {
    set_temp.insert(copyIndex[i-1]);
    for (int j = 1; j < n+1; ++j) {
      if (M_i.P[i][j] == 1){
        poset_temp(copyIndex[i-1],copyIndex[j-1]) = 1;
      }
    }
  }

  event_temp = set_temp;

  M_i.lin_ext = get_int_array(M_i.n);  // a linear extension of the poset
  M_i.J_P = bfs_order_ideals(M_i.P, M_i.n+1, &(M_i.m), M_i.lin_ext);

  /*         free         */
  free(copyIndex);
  free(union_lambda);
  for (int i = 0; i < rowN; ++i) {
    free_data2(D_SA[i],N_u[i],n,m_temp);
    free_int_Matrix(Data_SA[i],N[i]);
    free(lambda[i]);
  }
  std::vector<data2*>().swap(D_SA);
  std::vector<int**>().swap(Data_SA);
  std::vector<double*>().swap(lambda);
  epsilon.resize(0,0);
  std::vector<int>().swap(N_u);
  std::vector<int>().swap(N);
  std::vector<int>().swap(event);

}




// // // /*                  post_delete_poset                       */
// void post_delete_poset(dataD& dD, std::vector<int**>& GENOTYPE, int n_p=6, double T=10,int N_iter= 150,double S=1.0,int R = 1){
//   int n = dD.poset.rows();
//   std::vector<std::unordered_set<int>> parent(n);
//   std::vector<std::unordered_set<int>> children(n);
//
//   for (int i = 1; i < n; ++i) {
//     for (int j = 1; j < n; ++j) {
//       if (dD.poset(i,j) == 1){
//         parent[j].insert(i);
//         children[i].insert(j);
//       }
//     }
//   }
//
//   bool is_first = true;
//   std::vector<int> leaf_node;
//   std::mutex m;
//
//   while (is_first | leaf_node.size() > 0){
//     if(is_first){
//       for (int i = 1; i < n; ++i) {
//         if(children[i].size() == 0 && parent[i].size() != 0){
//           leaf_node.push_back(i);
//         }
//       }
//       is_first = false;
//     }
//
//     #pragma omp parallel for
//     for (int l = 0; l < leaf_node.size(); ++l) {
//       MyIntMatrix poset_temp;
//       poset_temp.resize(n,n);
//       poset_temp.setZero();
//       std::vector<std::unordered_set<int>> event_temp;
//
//       std::unordered_set<int> event_Use;
//       std::deque<int> event_D;
//       event_D.push_back(leaf_node[l]);
//
//       while (event_Use.size() < n_p){
//         if (event_D.size() > 0){
//           int temp = event_D.front();
//           event_D.pop_front();
//           event_Use.insert(temp);
//           //                    Must be a synchronization code block.
//           m.lock();
//           if (parent[temp].size() > 0){
//             auto i = parent[temp].begin();
//             while (i != parent[temp].end() && (event_Use.size()+event_D.size())<n_p){
//               event_D.push_back(*i);
//               int temp_1 = *i;
//               i++;
//               parent[temp].erase(temp_1);
//               children[temp_1].erase(temp);
//             }
//           }
//           m.unlock();
//         }else{
//           break;
//         }
//       }
//
//       poset_SA(dD,event_Use, poset_temp, event_temp, GENOTYPE,T, N_iter, S, R);
//
//       for (int i = 1; i < n; ++i) {
//         for (int j = 1; j < n; ++j) {
//           if (dD.poset(i,j)==1 && poset_temp(i,j)==0 && is_in_SA_mix(event_temp,i,j)){
//             dD.poset(i,j) = 0;
//           }
//         }
//       }
//     }
//
//     leaf_node.clear();
//     for (int i = 1; i < n; ++i) {
//       if(children[i].size() == 0 && parent[i].size() != 0){
//         leaf_node.push_back(i);
//       }
//     }
//   }
// }


// // /*                  post_delete_poset                       */
void post_delete_poset(dataD& dD, std::vector<int**>& GENOTYPE,
                       int n_p=6, bool is_update_eps = false, double T=10,int N_iter= 150,double S=1.0,int R = 1){
  int n = dD.poset.rows();
  //    std::vector<std::unordered_set<int>> parent(n);
  //    std::vector<std::unordered_set<int>> children(n);

  MyIntMatrix parent;
  parent.resize(n,n);
  parent.setZero();

  MyIntMatrix children;
  children.resize(n,n);
  children.setZero();

  for (int i = 1; i < n; ++i) {
    for (int j = 1; j < n; ++j) {
      if (dD.poset(i,j) == 1){
        parent(j,i) = 1;
        children(i,j) = 1;
      }
    }
  }

  std::deque<int> leaf_node;

  for (int i = 1; i < n; ++i) {
    if((children.row(i).array() == 0).all() &&  !((parent.row(i).array() == 0).all()) ){
      leaf_node.push_back(i);
    }
  }
  //    for (int i = 1; i <= n; ++i) {
  //        if(!((dD.poset.row(i).array() == 0).all()) && (children.row(i).array() == 0).all()){
  //            leaf_node.push_back(i);
  //        }
  //    }

  int parallel_index = 0;
  std::vector<std::unordered_set<int>> event_list_parallel;
  std::vector<MyIntMatrix> poset_temp_parallel;
  std::vector<std::unordered_set<int>> event_temp_parallel;

  while (leaf_node.size() > 0){
    std::unordered_set<int> event_Use;
    std::deque<int> event_D;
    int temp = leaf_node.front();
    leaf_node.pop_front();
    event_D.push_back(temp);
    std::unordered_set<int> event_Use_other;

    while (event_Use.size() < n_p){
      if(event_D.size() == 0){
        break;
      }
      int temp1 = event_D.front();
      event_D.pop_front();
      event_Use.insert(temp1);
      if ( !((parent.row(temp1).array() == 0).all()) ){
        for (int i = 1; i < n; ++i) {
          if(parent(temp1, i) == 1){
            event_D.push_back(i);
          }
        }
      }
      if(event_D.size() == 0 && event_Use.size() < n_p){
        for (const auto& element : event_Use) {
          for (int i = 1; i < n; ++i) {
            if(event_Use.count(i) == 0 && dD.poset(element,i) == 1){
              event_D.push_back(i);
              event_Use_other.insert(i);
            }
          }
        }
      }
    }

    std::unordered_set<int> event_Use_diff;
    for (const auto& element : event_Use) {
      if (event_Use_other.count(element) == 0) {
        event_Use_diff.insert(element);
      }
    }

    for (const auto& element : event_Use_diff) {
      if(element == temp){
        for (const auto& element1 : event_Use_diff) {
          if(element1 != element && parent(element, element1) == 1 ){
            parent(element, element1) = 0;
            children(element1, element) = 0;
          }
        }
      }else{
        if( (children.row(element).array() == 0).all() ){
          for (const auto& element1 : event_Use_diff) {
            if(element1 != element && parent(element, element1) == 1 ){
              parent(element, element1) = 0;
              children(element1, element) = 0;
            }
          }
        }
      }
    }

    for (int i = 1; i < n; ++i) {
      if((children.row(i).array() == 0).all() && !((parent.row(i).array() == 0).all()) ){
        if (std::find(leaf_node.begin(),leaf_node.end(),i) == leaf_node.end()){
          leaf_node.push_back(i);
        }
      }
    }

    event_list_parallel.push_back(event_Use);

    //        MyIntMatrix poset_temp;
    //        poset_temp.resize(n,n);
    //        poset_temp.setZero();
    //        poset_temp_parallel.push_back(poset_temp);
    //
    //        std::unordered_set<int> event_temp;
    //        event_temp_parallel.push_back(event_temp);
  }

  std::vector<std::unordered_set<int>> event_list_parallel_unique;

  for (int l = 0; l < event_list_parallel.size(); ++l) {
    if(l == 0){
      event_list_parallel_unique.push_back(event_list_parallel[l]);
      MyIntMatrix poset_temp;
      poset_temp.resize(n,n);
      poset_temp.setZero();
      poset_temp_parallel.push_back(poset_temp);
      std::unordered_set<int> event_temp;
      event_temp_parallel.push_back(event_temp);
    }else{
      bool is_same = false;
      for (int i = 0; i < event_list_parallel_unique.size(); ++i) {
        if(event_list_parallel[l] == event_list_parallel_unique[i]){
          is_same = true;
          break;
        }
      }
      if(!is_same){
        event_list_parallel_unique.push_back(event_list_parallel[l]);
        MyIntMatrix poset_temp;
        poset_temp.resize(n,n);
        poset_temp.setZero();
        poset_temp_parallel.push_back(poset_temp);
        std::unordered_set<int> event_temp;
        event_temp_parallel.push_back(event_temp);
      }

    }
  }

  #pragma omp parallel for
  for (int l = 0; l < event_list_parallel_unique.size(); ++l) {
    if(is_update_eps){
      poset_SA(dD, event_list_parallel_unique[l], poset_temp_parallel[l], event_temp_parallel[l], GENOTYPE,T, N_iter, S, R);
    }else{
      poset_SA_with_fixed_eps(dD, event_list_parallel_unique[l], poset_temp_parallel[l], event_temp_parallel[l], GENOTYPE,T, N_iter, S, R);
    }
  }

  //
  for (int l = 0; l < event_list_parallel_unique.size(); ++l) {
    for (int i = 1; i < n; ++i) {
      for (int j = 1; j < n; ++j) {
        if (dD.poset(i,j)==1 && poset_temp_parallel[l](i,j)==0 &&
            event_temp_parallel[l].find(i) != event_temp_parallel[l].end() &&
            event_temp_parallel[l].find(j) != event_temp_parallel[l].end()){
          dD.poset(i,j) = 0;
        }
      }
    }
  }



}



// /*                  delete_error_poset                       */
void poset_streamlining(dataD& dD){
  std::vector<int> eventNum(dD.poset.cols());

  for (int i = 0; i < dD.pat.rows(); ++i) {
    for (int j = 1; j < eventNum.size(); ++j) {
      if (dD.pat(i,j) == 1){
        eventNum[j] ++;
      }
    }
  }
  for (int i = 1; i < eventNum.size(); ++i) {
    for (int j = 1; j < eventNum.size(); ++j) {
      if (dD.poset(i,j) == 1 && (eventNum[i] <= eventNum[j]) ){
        dD.poset(i,j) = 0;
      }
    }
  }
}



// /*                  isCyclic Kahn                       */
bool isCyclicKahn(MyIntMatrix& adjMatrix) {
  int V = adjMatrix.rows();
  vector<int> inDegree(V, 0);

  for (int i = 0; i < V; ++i) {
    for (int j = 0; j < V; ++j) {
      if (adjMatrix(i,j) == 1) {
        inDegree[j]++;
      }
    }
  }

  std::queue<int> q;
  for (int i = 0; i < V; ++i) {
    if (inDegree[i] == 0) {
      q.push(i);
    }
  }

  int count = 0;

  while (!q.empty()) {
    int u = q.front();
    q.pop();
    count++;

    for (int v = 0; v < V; ++v) {
      if (adjMatrix(u,v) == 1) {
        inDegree[v]--;
        if (inDegree[v] == 0) {
          q.push(v);
        }
      }
    }
  }

  return (count != V);
}



// /*                  poset to DAG                      */
void poset_to_DAG(dataD& dD){
  std::vector<int> Count_temp(dD.poset.rows());
  for (int i = 1; i < dD.poset.rows(); ++i) {
    int temp = 0;
    for (int j = 0; j < dD.pat.rows(); ++j) {
      if (dD.pat(j,i)==1){
        temp ++;
      }
    }
    Count_temp[i] = temp;
  }

  std::vector<threeD> ECount;
  for (int i = 1; i < dD.poset.rows(); ++i) {
    for (int j = 1; j < dD.poset.rows(); ++j) {
      if(dD.poset(i,j) == 1){
        threeD temp;
        temp.F = i;
        temp.S = j;
        temp.T = Count_temp[j] - Count_temp[i];
        ECount.push_back(temp);
      }
    }
  }

  std::sort(ECount.begin(),ECount.end(),[](threeD a,threeD b){
    return a.T > b.T;
  });

  for (int i = 0; i < ECount.size(); ++i) {
    if(isCyclicKahn(dD.poset)){
      dD.poset(ECount[i].F, ECount[i].S) = 0;
    }else{
      break;
    }
  }

}





// /*                  search threshold with gmm                     */
double threshold_gmm(MyDoubleMatrix& loglik_matrix, int gmm_min_cluster = 2,double gmm_convenience_threshold = 0.0025,
                     double gmm_filter_threshold = 0.005){
  std::vector<double> data;

  for (int i = 0; i < loglik_matrix.rows(); ++i) {
    double e_temp0 = loglik_matrix(i,0);
    double e_temp1 = loglik_matrix(i,1);
    double e_temp2 = loglik_matrix(i,2);

    double val_1 = max(abs(e_temp0 - e_temp2), 0.000001);
    double val_2 = max(abs(e_temp1 - e_temp2), 0.000001);

    double sub_val_1_2 = min(val_1, val_2)/max(val_1, val_2);

    if ((e_temp0 - e_temp1) > 0 && (e_temp0 - e_temp2) > 0){
      data.push_back(sub_val_1_2);
    } else if ((e_temp1 - e_temp0) > 0 && (e_temp1 - e_temp2) > 0){
      data.push_back(sub_val_1_2);
    }
  }

  // for (int i = 0; i < data.size(); ++i) {
  //   cout << data[i] << endl;
  // }

  if (data.size() <= 4){
    return 0;
  }
  double threshold = -1.0;
  std::sort(data.begin(), data.end());
  std::vector<std::vector<double>> data_cluster;
  if(gmm_min_cluster < 2){
    gmm_min_cluster = 2;
  }
  int cluster_num = gmm_min_cluster;
  try {
    while (threshold < 0 && cluster_num <= int (data.size()/2)){
      std::vector<std::vector<double>> data_temp(cluster_num);
      int data_size = int (data.size()/cluster_num);
      for (int i = 0; i < data.size(); ++i) {
        int index_temp = min(int (i/data_size), cluster_num - 1);
        data_temp[index_temp].push_back(data[i]);
      }

      std::vector<double> mean_temp(cluster_num);
      std::vector<double> var_temp(cluster_num);
      for (int i = 0; i < cluster_num; ++i) {
        mean_temp[i] = std::accumulate(data_temp[i].begin(),data_temp[i].end(),0.0)/data_temp[i].size();
        double var_temp_2 = 0.0;
        for (double value : data_temp[i]) {
          var_temp_2 += std::pow(value - mean_temp[i], 2);
        }
        var_temp[i] = max(0.000001, var_temp_2/data_temp[i].size());
      }

      GaussianMixtureModel gmm(mean_temp, var_temp);
      cout << "---------------------------" << endl;
      std::cout << "cluster_num: " << cluster_num << endl;
      // std::cout << "Initial parameters: " << cluster_num << std::endl;
      // gmm.printParameters();
      std::vector<std::vector<double>> cluster_res = gmm.runEM(data);
      std::cout << "GMM parameters: " << cluster_num << std::endl;
      gmm.printParameters();
      std::cout << "cluster_res: " << cluster_num << std::endl;
      for (int i = 0; i < cluster_res.size(); ++i) {
        cout << i << ": " << cluster_res[i][0] << "  " << cluster_res[i][1] << endl;
      }
      if (gmm.mu[0] <= gmm_convenience_threshold){
        for (int i = 0; i < cluster_num; ++i) {
          if(gmm.mu[i] <= gmm_filter_threshold){
            threshold = cluster_res[i][1];
          }
        }
        break;
      }

      // cout << "---------------------------" << endl;
      cluster_num += 1;
    }
  }  catch (const std::exception& e) {
    std::cout << "Catch exception in GMM : " << e.what() << std::endl;
  }

  cout << "---------------------------" << endl;
  if (threshold <= 0){
    threshold = 0.0075;
  }

  return threshold;
}


/*                  find poset interface                       */
// void find_poset_mix(dataD& dD, int thrds, std::vector<int**>& GENOTYPE, int Fine_Tune_Num = 2, double threshold = 0,
//                     int gmm_min_cluster = 2, double gmm_convenience_threshold = 0.005,
//                     double gmm_filter_threshold = 0.0075, int n_p = 6, bool is_update_eps = false,
//                     double T=10,int N_iter= 150, double S=1.0,int R = 1){
//
//   // for (int i = 1; i < dD.poset.cols(); ++i) {
//   //   for (int j = i+1; j < dD.poset.cols(); ++j) {
//   //     cout << i << "  " << j << endl;
//   //     find_poset(dD,i,j, GENOTYPE[3], threshold, S, R);
//   //   }
//   // }
//
//   std::vector<doubleD> eventS;
//   for (int i = 1; i < dD.poset.cols(); ++i) {
//     for (int j = i+1; j < dD.poset.cols(); ++j) {
//       eventS.push_back(doubleD(i,j));
//     }
//   }
//
//   MyDoubleMatrix loglik_matrix;
//   loglik_matrix.resize(eventS.size(), 3);
//
//   int thrdsN = eventS.size();
//   int loop = (thrdsN/thrds) + (thrdsN%thrds>0?1:0);
//   for (int l = 0; l < loop; ++l) {
//     #pragma omp parallel for
//     for (int t = l*thrds; t < min((l+1)*thrds,thrdsN); ++t) {
//       int i = eventS[t].F;
//       int j = eventS[t].S;
//       std::vector<double> res_temp = find_poset(dD,i,j, GENOTYPE[3], S, R);
//       for (int res_in = 0; res_in < res_temp.size(); ++res_in) {
//         loglik_matrix(t,res_in) = res_temp[res_in];
//       }
//     }
//   }
//
//   if(threshold <= 0){
//     threshold = threshold_gmm(loglik_matrix, gmm_min_cluster, gmm_convenience_threshold, gmm_filter_threshold);
//   }
//
//   cout << "---------------------------" << endl;
//   std::cout << "threshold: " << threshold << endl;
//
//   for (int i = 0; i < eventS.size(); ++i) {
//     int col1 = eventS[i].F;
//     int col2 = eventS[i].S;
//     double e_temp0 = loglik_matrix(i,0);
//     double e_temp1 = loglik_matrix(i,1);
//     double e_temp2 = loglik_matrix(i,2);
//
//     double val_1 = max(abs(e_temp0 - e_temp2), 0.000001);
//     double val_2 = max(abs(e_temp1 - e_temp2), 0.000001);
//
//     double sub_val_1_2 = min(val_1, val_2)/max(val_1, val_2);
//
//     if ((e_temp0 - e_temp1) > 0 && (e_temp0 - e_temp2) > 0 && sub_val_1_2 > threshold ){
//       dD.poset(col1,col2) = 1;
//     } else if ((e_temp1 - e_temp0) > 0 && (e_temp1 - e_temp2) > 0 && sub_val_1_2 > threshold ){
//       dD.poset(col2,col1) = 1;
//     }
//
//   }
//
//   while(Fine_Tune_Num > 0){
//     std::cout << "n_p: " << n_p << endl;
//     post_delete_poset(dD, GENOTYPE, n_p, is_update_eps, T,N_iter,S,R);
//     Fine_Tune_Num -= 1;
//   }
//
//   // if(is_Fine_Tune){
//   //   std::cout << "n_p: " << n_p << endl;
//   //   post_delete_poset(dD, GENOTYPE, n_p,T,N_iter,S,R);
//   //   std::cout << "n_p: " << n_p << endl;
//   //   post_delete_poset(dD, GENOTYPE, n_p,T,N_iter,S,R);
//   // }
//
//   poset_streamlining(dD);
//
// }
//

void find_poset_mix(dataD& dD, int thrds, std::vector<int**>& GENOTYPE, int Fine_Tune_Num = 2, double threshold = 0, double threshold2 = 0.01,
                    int gmm_min_cluster = 2, double gmm_convenience_threshold = 0.005,
                    double gmm_filter_threshold = 0.0075, int n_p = 6, bool is_update_eps = false,
                    double T=10,int N_iter= 150, double S=1.0,int R = 1){


  std::vector<doubleD> eventS;
  for (int i = 1; i < dD.poset.cols(); ++i) {
    for (int j = i+1; j < dD.poset.cols(); ++j) {
      eventS.push_back(doubleD(i,j));
    }
  }

  MyDoubleMatrix loglik_matrix;
  loglik_matrix.resize(eventS.size(), 3);

  int thrdsN = eventS.size();
  int loop = (thrdsN/thrds) + (thrdsN%thrds>0?1:0);
  for (int l = 0; l < loop; ++l) {
  #pragma omp parallel for
    for (int t = l*thrds; t < min((l+1)*thrds,thrdsN); ++t) {
      int i = eventS[t].F;
      int j = eventS[t].S;
      std::vector<double> res_temp = find_poset(dD,i,j, GENOTYPE[3], S, R);
      for (int res_in = 0; res_in < res_temp.size(); ++res_in) {
        loglik_matrix(t,res_in) = res_temp[res_in];
      }
    }
  }

  if(threshold <= 0){
    threshold = threshold_gmm(loglik_matrix, gmm_min_cluster, gmm_convenience_threshold, gmm_filter_threshold);
  }

  cout << "---------------------------" << endl;
  std::cout << "threshold: " << threshold << endl;
  std::cout << "threshold2: " << threshold2 << endl;


  for (int i = 0; i < eventS.size(); ++i) {
    int col1 = eventS[i].F;
    int col2 = eventS[i].S;
    double e_temp0 = loglik_matrix(i,0);
    double e_temp1 = loglik_matrix(i,1);
    double e_temp2 = loglik_matrix(i,2);

    double val_1 = max(abs(e_temp0 - e_temp2), 0.000001);
    double val_2 = max(abs(e_temp1 - e_temp2), 0.000001);

    double sub_val_1_2 = min(val_1, val_2)/max(val_1, val_2);

    if ((e_temp0 - e_temp1) > 0 && (e_temp0 - e_temp2) > 0 && sub_val_1_2 > threshold ){
      dD.poset(col1,col2) = 1;
    }else if((e_temp1 - e_temp0) > 0 && (e_temp1 - e_temp2) > 0 && sub_val_1_2 > threshold ){
      dD.poset(col2,col1) = 1;
    }else if((e_temp0 - e_temp1) > 0 && (e_temp0 - e_temp2) <= 0 && sub_val_1_2 < threshold2){
      dD.poset(col1,col2) = 1;
    }else if((e_temp1 - e_temp0) > 0 && (e_temp1 - e_temp2) <= 0 && sub_val_1_2 < threshold2){
      dD.poset(col2,col1) = 1;
    }

  }

  poset_to_DAG(dD);

  for (int i = 0; i < Fine_Tune_Num; ++i) {
    post_delete_poset(dD, GENOTYPE, n_p, is_update_eps, T,N_iter,S,R);
  }


  // if(is_Fine_Tune){
  //   std::cout << "n_p: " << n_p << endl;
  //   post_delete_poset(dD, GENOTYPE, n_p,T,N_iter,S,R);
  //   std::cout << "n_p: " << n_p << endl;
  //   post_delete_poset(dD, GENOTYPE, n_p,T,N_iter,S,R);
  // }

  poset_streamlining(dD);

}



//------------------------------------------------------------------------------------------------------------------------
/*                  MCEM hcbn                       */
void MCEM_hcbn(Model& model, MatrixXb& obs, unsigned int L, std::string& sampling, ControlEM& control_EM,
               const unsigned int thrds, Context& ctx) {

    // Initialization and instantiation of variables
    const vertices_size_type p = model._size; // Number of mutations / events
    int rowN = model._epsilon.rows();
    int eN = model._eventD.size();

    const unsigned int N = obs.rows();         // Number of observations / genotypes
    double N_eff = 0.0;
    unsigned int K = 0;
    unsigned int update_step_size = control_EM.update_step_size;
    VectorXd avg_lambda = VectorXd::Zero(p);
    VectorXd avg_lambda_current = VectorXd::Zero(p);

    MyDoubleMatrix avg_eps(rowN,p);
    MyDoubleMatrix avg_eps_current(rowN,p);

    bool tol_comparison = true;
    MyDoubleMatrix expected_dist = MyDoubleMatrix::Zero(rowN,eN);

    MatrixXd expected_Tdiff = MatrixXd::Zero(N, p);
    VectorXd Tdiff_colsum(p);
    MatrixXd Tdiff_pool;
    VectorXd scale_cumulative;


    if (sampling == "add-remove") {
        scale_cumulative.resize(p);
        if (model.get_update_node_idx())
            model.update_node_idx();
    } else if (sampling == "pool") {
        K = p * L;
        Tdiff_pool.resize(K, p);

    }

    for (unsigned int iter = 0; iter < control_EM.max_iter; ++iter) {
        MatrixXd T_pool;
        /*    update eps & lambda    */
        if (iter == update_step_size) {
            avg_lambda_current /= control_EM.update_step_size;
            avg_eps_current /= control_EM.update_step_size;
            if (tol_comparison) {
                if (matrixDivision1(avg_eps, avg_eps_current, control_EM.tol) &&
                    ((avg_lambda - avg_lambda_current).array().abs() <= control_EM.tol).all())
                    break;
                // L *= 2;
            }
            avg_lambda = avg_lambda_current;
            avg_eps = avg_eps_current;

            update_step_size += control_EM.update_step_size;
            tol_comparison = true;

            /* Restart averaging */
            avg_lambda_current = VectorXd::Zero(p);
            avg_eps_current.fill(0.0);
        }

        /* E step
         * Conditional expectation for the sufficient statistics per observation and event
         * scale_cumulative: cumulative E(lambda)
         */
        if (sampling == "add-remove") {
            scale_cumulative = scale_path_to_mutation(model);
        } else if (sampling == "pool") {
            /* All threads share the same pool of mutation times */
            T_pool.resize(K, p);
            /* Generate observation times from a given poset and given rates */
            T_pool = sample_times(K, model, Tdiff_pool, ctx.rng);
        }

        expected_dist.fill(0.0);

        /*     generate thrds rng     */
        auto rngs = ctx.get_auxiliary_rngs(thrds);
        #pragma omp declare reduction(matrix_add: MyDoubleMatrix:omp_out=add_matrix(omp_out,omp_in)) initializer(omp_priv = omp_orig)
        {
            #pragma omp parallel for reduction(matrix_add:expected_dist) schedule(static)
            for (unsigned int i = 0; i < N; ++i) {
                int rowI = 0;
                for (int j = 0; j < model._setD.size(); ++j) {
                    if (i>=model._setD[j].F && i<=model._setD[j].S){
                        rowI = j;
                        break;
                    }
                }

                MyIntMatrix d_pool;
                if (sampling == "pool") {
                  VectorXd T_sampling(K);
                  /*    generate genotype    */
                  MatrixXb genotype_pool =
                              generate_genotypes(T_pool, model, T_sampling, (*rngs)[omp_get_thread_num()]);
                  d_pool = hamming_dist_mat1(genotype_pool, obs.row(i), model);
                }

                DataImportanceSampling importance_sampling = importance_weight(
                        obs.row(i), L, model, rowI, sampling, scale_cumulative, d_pool, Tdiff_pool,
                        control_EM.neighborhood_dist, (*rngs)[omp_get_thread_num()]);

                double aux = importance_sampling.w.sum();
                if (aux > 0) {
                    for (int j = 0; j < eN; ++j) {
                        expected_dist(rowI,j) += importance_sampling.w.dot(importance_sampling.dist.col(j).cast<double>()) / aux;
                    }
                    expected_Tdiff.row(i) = (importance_sampling.Tdiff.transpose() * importance_sampling.w) / aux;
                }
            }
        }

        /* M-step */
        for (int i = 0; i < rowN; ++i) {
            for (int j = 0; j < eN; ++j) {
                expected_dist(i,j) = expected_dist(i,j)/((model._setD[i].S-model._setD[i].F+1)*(model._eventD[j].S-model._eventD[j].F+1));
            }
        }
        MyDoubleMatrix epsTemp(rowN,p);
        // cout << expected_dist << endl;
        for (int i = 0; i < eN; ++i) {
            for (int j = model._eventD[i].F; j <= model._eventD[i].S; ++j) {
                for (int k = 0; k < rowN; ++k) {
                    epsTemp(k,j-1) = expected_dist(k,i);
                }
            }
        }
        model._epsilon = epsTemp;

        RowVectorXd weights(N);
        weights.fill(1.0);
        Tdiff_colsum = weights*expected_Tdiff;
        model.set_lambda((Tdiff_colsum/N).array().inverse(), control_EM.max_lambda);
        avg_lambda_current +=  model._lambda;
        avg_eps_current += model._epsilon;

        if (iter + 1 == control_EM.max_iter) {
            unsigned int num_iter = control_EM.max_iter - update_step_size +
                                    control_EM.update_step_size;
            avg_lambda_current /= num_iter;
            avg_eps_current /= num_iter;
        }

    }

    model._lambda = avg_lambda_current;
    model._epsilon = avg_eps_current;

}


void MCEM_hcbn_with_fixed_eps(Model& model, MatrixXb& obs, unsigned int L, std::string& sampling, ControlEM& control_EM,
               const unsigned int thrds, Context& ctx) {

  // Initialization and instantiation of variables
  const vertices_size_type p = model._size; // Number of mutations / events
  int rowN = model._epsilon.rows();
  int eN = model._eventD.size();

  const unsigned int N = obs.rows();         // Number of observations / genotypes
  double N_eff = 0.0;
  unsigned int K = 0;
  unsigned int update_step_size = control_EM.update_step_size;
  VectorXd avg_lambda = VectorXd::Zero(p);
  VectorXd avg_lambda_current = VectorXd::Zero(p);

  bool tol_comparison = true;
  MyDoubleMatrix expected_dist = MyDoubleMatrix::Zero(rowN,eN);

  MatrixXd expected_Tdiff = MatrixXd::Zero(N, p);
  VectorXd Tdiff_colsum(p);
  MatrixXd Tdiff_pool;
  VectorXd scale_cumulative;


  if (sampling == "add-remove") {
    scale_cumulative.resize(p);
    if (model.get_update_node_idx())
      model.update_node_idx();
  } else if (sampling == "pool") {
    K = p * L;
    Tdiff_pool.resize(K, p);

  }

  for (unsigned int iter = 0; iter < control_EM.max_iter; ++iter) {
    MatrixXd T_pool;
    /*    update eps & lambda    */
    if (iter == update_step_size) {
      avg_lambda_current /= control_EM.update_step_size;

      if (tol_comparison) {
        if (((avg_lambda - avg_lambda_current).array().abs() <= control_EM.tol).all())
          break;
        // L *= 2;
      }
      avg_lambda = avg_lambda_current;

      update_step_size += control_EM.update_step_size;
      tol_comparison = true;

      /* Restart averaging */
      avg_lambda_current = VectorXd::Zero(p);
    }

    /* E step
     * Conditional expectation for the sufficient statistics per observation and event
     * scale_cumulative: cumulative E(lambda)
     */
    if (sampling == "add-remove") {
      scale_cumulative = scale_path_to_mutation(model);
    } else if (sampling == "pool") {
      /* All threads share the same pool of mutation times */
      T_pool.resize(K, p);
      /* Generate observation times from a given poset and given rates */
      T_pool = sample_times(K, model, Tdiff_pool, ctx.rng);
    }

    expected_dist.fill(0.0);

    /*     generate thrds rng     */
    auto rngs = ctx.get_auxiliary_rngs(thrds);
    #pragma omp declare reduction(matrix_add: MyDoubleMatrix:omp_out=add_matrix(omp_out,omp_in)) initializer(omp_priv = omp_orig)
    {
      #pragma omp parallel for reduction(matrix_add:expected_dist) schedule(static)
      for (unsigned int i = 0; i < N; ++i) {
        int rowI = 0;
        for (int j = 0; j < model._setD.size(); ++j) {
          if (i>=model._setD[j].F && i<=model._setD[j].S){
            rowI = j;
            break;
          }
        }

        MyIntMatrix d_pool;
        if (sampling == "pool") {
          VectorXd T_sampling(K);
          /*    generate genotype    */
          MatrixXb genotype_pool =
          generate_genotypes(T_pool, model, T_sampling, (*rngs)[omp_get_thread_num()]);
          d_pool = hamming_dist_mat1(genotype_pool, obs.row(i), model);
        }

        DataImportanceSampling importance_sampling = importance_weight(
          obs.row(i), L, model, rowI, sampling, scale_cumulative, d_pool, Tdiff_pool,
          control_EM.neighborhood_dist, (*rngs)[omp_get_thread_num()]);

        double aux = importance_sampling.w.sum();
        if (aux > 0) {
          for (int j = 0; j < eN; ++j) {
            expected_dist(rowI,j) += importance_sampling.w.dot(importance_sampling.dist.col(j).cast<double>()) / aux;
          }
          expected_Tdiff.row(i) = (importance_sampling.Tdiff.transpose() * importance_sampling.w) / aux;
        }
      }
    }

    /* M-step */
    for (int i = 0; i < rowN; ++i) {
      for (int j = 0; j < eN; ++j) {
        expected_dist(i,j) = expected_dist(i,j)/((model._setD[i].S-model._setD[i].F+1)*(model._eventD[j].S-model._eventD[j].F+1));
      }
    }


    RowVectorXd weights(N);
    weights.fill(1.0);
    Tdiff_colsum = weights*expected_Tdiff;
    model.set_lambda((Tdiff_colsum/N).array().inverse(), control_EM.max_lambda);
    avg_lambda_current +=  model._lambda;


    if (iter + 1 == control_EM.max_iter) {
      unsigned int num_iter = control_EM.max_iter - update_step_size +
        control_EM.update_step_size;
      avg_lambda_current /= num_iter;
    }

  }

  model._lambda = avg_lambda_current;

}



/*                  estimate lambda                       */
void lambda_estimate(VectorXd& lambdaRes, MyIntMatrix& poset, MyBoolMatrix& obs, MyDoubleMatrix& eps, std::vector<doubleD>& setD,
                     std::vector<doubleD>& eventD, unsigned int L, std::string& sampling, bool is_Update_eps, int seed,
                     int thrds, double lambda_s=1.0, int max_iter=100L, int update_step_size=20L,
                     double tol=0.001, double max_lambda=1e6, unsigned int neighborhood_dist=1){

  const auto p = poset.rows(); // Number of mutations / events
  VectorXd ilambda(p+1);
  uniform_real_distribution<double> u(lambda_s/3, lambda_s*3);
  default_random_engine e(time(NULL));
  ilambda(0) = lambda_s;
  for (int i = 1; i < p + 1; ++i) {
    ilambda(i) = u(e);
  }

  VectorXd times;
  times.resize(obs.rows());
  bool sampling_times_available = false;

  edge_container edge_list = adjacency_mat2list(poset);
  Model M(edge_list, p, lambda_s);
  M.set_lambda(ilambda, max_lambda);
  M._epsilon = eps;
  M._setD = setD;
  M._eventD = eventD;
  M.has_cycles();
  if (M.cycle){
    cout << "Poset has cycle, please remove some edges for the poset and re-enter it!";
    return;
  }
  M.topological_sort();

  ControlEM control_EM(max_iter, update_step_size, tol, max_lambda, neighborhood_dist);

  /* Call the underlying C++ function */
  Context ctx(seed);

  if(is_Update_eps){
    MCEM_hcbn(M, obs, L, sampling, control_EM, thrds, ctx);
  }else{
    MCEM_hcbn_with_fixed_eps(M, obs, L, sampling, control_EM, thrds, ctx);
  }

  lambdaRes = M._lambda;
  eps = M._epsilon;
}
















#endif //INTERFACEFUN_H
