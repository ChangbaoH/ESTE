#ifndef ESTE_LAMBDACAL_H
#define ESTE_LAMBDACAL_H
#include <Rcpp.h>

#include "ct-cbn.h"
#include <random>
#include <cmath>
#include <random>
#include <vector>
#include <memory>
#include "rng_utils.hpp"
#include <exception>
#include <string>
#include <cassert>
#include <thread>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/topological_sort.hpp>
#include <boost/math/special_functions/binomial.hpp>

using Eigen::Map;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::RowVectorXd;
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;
typedef Eigen::Matrix<bool, 1, Eigen::Dynamic> RowVectorXb;
typedef Map<VectorXd> MapVecd;
typedef Map<VectorXi> MapVeci;
typedef Map<MatrixXi> MapMati;
typedef Map<MatrixXd> MapMatd;
typedef Map<RowVectorXd> MapRowVecd;

struct Event {
    unsigned int event_id;
};

/* Define Graph
 * edge: hash_setS   node:vecS   kind:bidirectionalS()   VertexProperties:Event
 */
typedef boost::adjacency_list<boost::hash_setS, boost::vecS, boost::bidirectionalS, Event> Poset;
typedef boost::graph_traits<Poset>::vertices_size_type vertices_size_type;
typedef boost::graph_traits<Poset>::vertex_descriptor Node;
typedef std::pair<vertices_size_type, vertices_size_type> Edge;
typedef std::vector<Node> node_container;
typedef std::vector<Edge> edge_container;


class Context {
public:
    typedef ::rng_type rng_type;
    typedef boost::ptr_vector<rng_type> rng_vector_type;

    rng_type rng;

    Context(int seed): rng(seed) {}

    /*    num_rngs:thrds   reserve some space  */
    std::unique_ptr<rng_vector_type> get_auxiliary_rngs(int num_rngs) {
        std::unique_ptr<rng_vector_type> rngs(new rng_vector_type());
        rngs->reserve(num_rngs);
        for (int t = 0; t < num_rngs; ++t)
            rngs->push_back(new rng_type(rng()));
        return std::move(rngs);
    }
};



class Model {
public:
    Poset poset;                // Adjacency list encoding cover releations
    node_container topo_path;   // A topological ordering of the nodes
    bool cycle;
    VectorXd _lambda;
    float _lambda_s;
    MyDoubleMatrix _epsilon;
    std::vector<doubleD> _setD;
    std::vector<doubleD> _eventD;
    std::vector<int> _node_idx;
    bool _update_node_idx;
    vertices_size_type _size;


    /* Constructor using the edge iterator constructor */
    Model(const edge_container& edge_list, unsigned int p, float lambda_s=1.0) :
            poset(edge_list.begin(), edge_list.end(), p), cycle(false),
            _lambda(p), _lambda_s(lambda_s), _node_idx(p), _size(p) {
        topo_path.reserve(p);
        for (unsigned int i = 0; i < p; ++i) {
            poset[i].event_id = i;
            _node_idx[i] = i;
            _update_node_idx = false;
        }
    }


    void set_lambda(const Eigen::Ref<const VectorXd>& lambda,
                    const float max_lambda);

    std::vector<node_container> get_direct_predecessors() const;

    void has_cycles();

    void topological_sort();

    void update_node_idx();

    inline const int get_node_idx(const unsigned int i) const;

    inline bool get_update_node_idx() const;
};



class DataImportanceSampling {
public:
    VectorXd w;
    MyIntMatrix dist;
    MatrixXd Tdiff;

    /*Parametrized Constructor */
    DataImportanceSampling(int L, int p ,int eN) : w(L), dist(L, eN),
                                                             Tdiff(L, p) {
        w.setConstant(1.0);
    }

};



/* Class containing customisable options for the EM algorithm */
class ControlEM {
public:
    unsigned int max_iter;
    unsigned int update_step_size; // increase L every 'update_step_size' to reach a desirable 'tol'
    double tol;                    // convergence tolerance
    float max_lambda;
    unsigned int neighborhood_dist;

    ControlEM(unsigned int max_iter=100, unsigned int update_step_size=20,
              double tol=0.001, float max_lambda=1e6,
              unsigned int neighborhood_dist=1) :
            max_iter(max_iter), update_step_size(update_step_size), tol(tol),
            max_lambda(max_lambda), neighborhood_dist(neighborhood_dist) {}
};



/**
 * Exception thrown when a cycle is found.
 */
class not_acyclic_exception : public std::exception {
public:
    not_acyclic_exception(
            const char *error = "Poset does not correspond to an acyclic graph") {
        errorMessage = error;
    }

    const char *what() const noexcept {
        return errorMessage.c_str();
    }

private:
    std::string errorMessage;
};


/*                                        backward_sampling.h                                        */
int n_choose_k(unsigned int n, unsigned int k);

void neighbors(const unsigned int p, unsigned int k,
               Eigen::Ref<MatrixXb> output_mat);

void coin_tossing(MatrixXb& output_mat, const Model& model, int rowI, Context::rng_type& rng);


/*                                        compatible_observations.h                                        */
bool is_compatible(const RowVectorXb& genotype, const Model& model);


/*                                        lambdaCal.h                                        */
edge_container adjacency_mat2list(const MatrixXi& poset);



/*                                        add_remove.h                                        */
VectorXd scale_path_to_mutation(const Model& model);

RowVectorXb draw_sample(const RowVectorXb& genotype, const Model& model,
                        const unsigned int move, const VectorXd& remove_weight,
                        const VectorXd& add_weight, double& q_choice,
                        const int idx_remove, const int idx_add,
                        bool compatible);

MatrixXd generate_mutation_times(const MatrixXb& obs, const Model& model, VectorXd& dens, VectorXd& sampling_time,
                                 Context::rng_type& rng);

VectorXd cbn_density_log(const MatrixXd& time, const VectorXd& lambda);


// ----------------------------------------------------------------------------------------------------------------------
// Model implementation

/*                                         lambdaCal.h                                        */
/*                                        model function                                        */
void Model::set_lambda(const Eigen::Ref<const VectorXd>& lambda, const float max_lambda) {
  _lambda = lambda;
  for (unsigned int j = 0; j < lambda.size(); ++j)
    if (lambda[j] > max_lambda || !std::isfinite(lambda[j]))
      _lambda[j] = max_lambda;
}

//' @description Obtain (direct) predecessors/parents per node
 std::vector<node_container> Model::get_direct_predecessors() const {

   std::vector<node_container> parents(_size);
   /* Loop through nodes in topological order */
   for (node_container::const_reverse_iterator v = topo_path.rbegin();
        v != topo_path.rend(); ++v) {
     /* Loop through (direct) predecessors/parents of node v */
     boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
     for (boost::tie(in_begin, in_end) = boost::in_edges(*v, poset);
          in_begin != in_end; ++in_begin) {
       Node u = boost::get(boost::get(&Event::event_id, poset), source(*in_begin, poset));
       parents[poset[*v].event_id].push_back(u);
     }
   }
   return parents;
 }

void Model::has_cycles() {
  /* Check for cycles using strongly connected components as proxy */
  std::vector<vertices_size_type> component(_size);
  vertices_size_type num_scc = boost::strong_components(
    poset,
    boost::make_iterator_property_map(
      component.begin(), boost::get(boost::vertex_index, poset)));
  if (num_scc != _size)
    cycle = true;
}

void Model::topological_sort() {
  topo_path.clear();
  boost::topological_sort(poset, std::back_inserter(topo_path));
}

void Model::update_node_idx() {
  for (unsigned int i = 0; i < _size; ++i)
    _node_idx[poset[i].event_id] = i;
  _update_node_idx = false;
}

const int Model::get_node_idx(const unsigned int i) const {
  return _node_idx[i];
}

bool Model::get_update_node_idx() const {
  return _update_node_idx;
}


/*                                        backward_sampling.h                                        */
int n_choose_k(unsigned int n, unsigned int k) {
  if (k > n)
    return 0;
  return boost::math::binomial_coefficient<double>(n, k);
}

void neighbors(const unsigned int p, unsigned int k,
               Eigen::Ref<MatrixXb> output_mat) {

  if (p >= 1) {
    unsigned int nrows = n_choose_k(p - 1, k);
    unsigned int nrows_alt = n_choose_k(p - 1, k - 1);

    /* Mutation as observed */
    if (nrows > 0)
      neighbors(p - 1, k, output_mat.block(0, 1, nrows, p - 1));

    /* Change the observed mutation */
    if (nrows_alt > 0) {
      output_mat.block(nrows, 0, nrows_alt, 1).setConstant(!output_mat(nrows, 0));
      neighbors(p - 1, k - 1, output_mat.block(nrows, 1, nrows_alt, p - 1));
    }
  }
}

void runif_std(const unsigned int N, std::vector<double>& output,
               Context::rng_type& rng) {

  std::uniform_real_distribution<double> distribution(0, 1);
  for (unsigned int i = 0; i < N; ++i)
    output[i] = distribution(rng);
}

void coin_tossing(MatrixXb& output_mat, const Model& model, int rowI, Context::rng_type& rng) {

  const unsigned int p = output_mat.cols();
  const unsigned int N = output_mat.rows();
  std::vector<double> result(p * N);
  runif_std(result.size(), result, rng);

  MatrixXd mask(p, N);
  mask = Map<MatrixXd>(result.data(), N, p);

  for (unsigned int i = 0; i < p; ++i){
    int temp = 0;
    for (int j = 0; j < model._eventD.size(); ++j) {
      if (i>=(model._eventD[j].F-1) && i<=(model._eventD[j].S-1)){
        temp = j;
        break;
      }
    }
    output_mat.col(i) = (mask.col(i).array() < model._epsilon(rowI,temp)).
    select(!output_mat(0, i), output_mat.col(i));
  }

}


/*                                        compatible_observations.h                                        */
bool is_compatible(const RowVectorXb& genotype, const Model& model) {

  auto id = boost::get(&Event::event_id, model.poset);
  /* Loop through nodes in reverse topological order */
  for (node_container::const_iterator v = model.topo_path.begin(); v != model.topo_path.end(); ++v) {
    if (genotype[model.poset[*v].event_id]) {
      /* Loop through (direct) predecessors/parents of node v */
      boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
      for (boost::tie(in_begin, in_end) = boost::in_edges(*v, model.poset); in_begin != in_end; ++in_begin)
        if (!genotype[boost::get(id, source(*in_begin, model.poset))])
          return false;
    }
  }
  return true;
}



/*                                        lambdaCal.h                                        */
//' @description Extract cover relations from an adjacency matrix
 edge_container adjacency_mat2list(const MatrixXi& poset) {
   const auto p = poset.rows();
   edge_container edge_list;
   for (unsigned int i = 0; i < p; ++ i)
     for (unsigned int j = 0; j < p; ++j)
       if (poset(i, j) == 1)
         edge_list.push_back(Edge(i, j));

       return edge_list;
 }


/*                                        add_remove.h                                        */
/*    return cumulative E(lambda_i)    */
VectorXd scale_path_to_mutation(const Model& model) {

  VectorXd scale_cumulative = VectorXd::Zero(model._size);
  /*   get the node's property map of graph  */
  auto id = boost::get(&Event::event_id, model.poset);
  /* Loop through nodes in topological order */
  for (node_container::const_reverse_iterator v = model.topo_path.rbegin(); v != model.topo_path.rend(); ++v) {
    /*     scale_cumulative:E(lambda)     */
    scale_cumulative[model.poset[*v].event_id] = 1.0 / model._lambda[model.poset[*v].event_id];
    double max_scale = -1.0;
    /* Loop through (direct) predecessors/parents of node v */
    boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
    for (boost::tie(in_begin, in_end) = boost::in_edges(*v, model.poset); in_begin != in_end; ++in_begin)
      max_scale = std::max(max_scale, scale_cumulative[boost::get(id, source(*in_begin, model.poset))]);

    if (max_scale > -1.0)
      scale_cumulative[model.poset[*v].event_id] += max_scale;
  }
  return scale_cumulative;
}

void add_all(RowVectorXb& genotype, const Model& model) {

  auto id = boost::get(&Event::event_id, model.poset);
  /* Loop through nodes in reverse topological order */
  for (node_container::const_iterator v = model.topo_path.begin();
       v != model.topo_path.end(); ++v) {
    if (genotype[model.poset[*v].event_id]) {
      /* Loop through (direct) predecessors/parents of node v */
      boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
      for (boost::tie(in_begin, in_end) = boost::in_edges(*v, model.poset);
           in_begin != in_end; ++in_begin) {
        Node u = boost::get(id, source(*in_begin, model.poset));
        if (!genotype[u])
          genotype[u] = true;
      }
    }
  }
}

void add_parents(RowVectorXb& genotype, const Node v, const Model& model) {

  assert(genotype[model.poset[v].event_id]);
  auto id = boost::get(&Event::event_id, model.poset);
  /* Loop through predecessors/parents of node v */
  boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
  for (boost::tie(in_begin, in_end) = boost::in_edges(v, model.poset); in_begin != in_end; ++in_begin) {
    Node u = boost::get(id, source(*in_begin, model.poset));
    if (!genotype[u]) {
      genotype[u] = true;
      add_parents(genotype, source(*in_begin, model.poset), model);
    }
  }
}

void remove_all(RowVectorXb& genotype, const Model& model) {

  auto id = boost::get(&Event::event_id, model.poset);
  /* Loop through nodes in topological order */
  for (node_container::const_reverse_iterator v = model.topo_path.rbegin();
       v != model.topo_path.rend(); ++v) {
    if (!genotype[model.poset[*v].event_id]) {
      /* Loop through (direct) successor/children of node v */
      boost::graph_traits<Poset>::out_edge_iterator out_begin, out_end;
      for (boost::tie(out_begin, out_end) = out_edges(*v, model.poset);
           out_begin != out_end; ++out_begin) {
        Node w = boost::get(id, target(*out_begin, model.poset));
        if (genotype[w])
          genotype[w] = false;
      }
    }
  }
}

void remove_children(RowVectorXb& genotype, const Node v, const Model& model) {

  assert(!genotype[model.poset[v].event_id]);
  auto id = boost::get(&Event::event_id, model.poset);
  /* Loop through (direct) successor/children of node v */
  boost::graph_traits<Poset>::out_edge_iterator out_begin, out_end;
  for (boost::tie(out_begin, out_end) = out_edges(v, model.poset);
       out_begin != out_end; ++out_begin) {
    Node w = boost::get(id, target(*out_begin, model.poset));
    if (genotype[w]) {
      genotype[w] = false;
      remove_children(genotype, target(*out_begin, model.poset), model);
    }
  }
}

RowVectorXb draw_sample(const RowVectorXb& genotype, const Model& model,
                        const unsigned int move, const VectorXd& remove_weight,
                        const VectorXd& add_weight, double& q_choice,
                        const int idx_remove, const int idx_add,
                        bool compatible) {

  RowVectorXb sample = genotype;
  switch (move) {
  case 0 :
    /* Pick an event to be added */
    q_choice = add_weight[idx_add] / add_weight.sum();
    sample[idx_add] = true;
    if (compatible)
      add_parents(sample, model.get_node_idx(idx_add), model);
    else
      add_all(sample, model);
    break;
  case 1 :
    /* Pick an event to be removed */
    q_choice = remove_weight[idx_remove] / remove_weight.sum();
    sample[idx_remove] = false;
    if (compatible)
      remove_children(sample, model.get_node_idx(idx_remove), model);
    else
      remove_all(sample, model);
    break;
  case 2 :
    /* Stand-still */
    q_choice = 1;
    break;
  }
  return sample;
}

/* Probability density function of an exponential distribution in log scale */
VectorXd dexp_log(const VectorXd& time, double rate) {
  VectorXd ret = std::log(rate) - (rate * time).array();
  return ret;
}

/* Cumulative distribution function of an exponential distribution in log scale */
VectorXd pexp_log(VectorXd& time, double rate) {
  // FIXME: -(-rate * time).array().expm1().array().log()
  VectorXd ret = (1 - (-rate * time).array().exp()).array().log();
  return ret;
}

MatrixXd generate_mutation_times(const MatrixXb& obs, const Model& model, VectorXd& dens, VectorXd& sampling_time,
                                 Context::rng_type& rng) {

  unsigned int N = obs.rows();
  unsigned int p = obs.cols();
  MatrixXd time_events = MatrixXd::Zero(N, p);
  MatrixXd time_events_sum = MatrixXd::Zero(N, p);
  MatrixXd cutoff = MatrixXd::Zero(N, p);

  /* Generate sampling times sampling_time ~ Exp(lambda_{s}) */
  sampling_time = rexp(N, model._lambda_s, rng);

  std::vector<node_container> parents = model.get_direct_predecessors();
  /* Loop through nodes in topological order */
  for (node_container::const_reverse_iterator v = model.topo_path.rbegin();
       v != model.topo_path.rend(); ++v) {
    /* Alternative vec1.cwiseMax(vec2)
     * see https://eigen.tuxfamily.org/dox/group__QuickRefPage.html)
     */
    VectorXd time_parents_max = VectorXd::Zero(N);
    if (!parents[model.poset[*v].event_id].empty()) {
      MatrixXd aux1 = MatrixXd::Zero(N, parents[model.poset[*v].event_id].size());

      for (unsigned int u = 0; u < parents[model.poset[*v].event_id].size(); ++u)
        aux1.col(u) = time_events_sum.col(parents[model.poset[*v].event_id][u]);
      time_parents_max = aux1.rowwise().maxCoeff();
    }

    /* if x = 1, Z ~ TExp(lambda, 0, sampling_time - time{max parents})
     * if x = 0, Z ~ TExp(lambda, 0, inf)
     */
    VectorXd cutoff = obs.col(model.poset[*v].event_id).select(sampling_time - time_parents_max,
                              std::numeric_limits<double>::infinity());
    VectorXd time = rtexp(N, model._lambda[model.poset[*v].event_id], cutoff, rng);

    VectorXd aux2 = obs.col(model.poset[*v].event_id).select(time_parents_max,
                            time_parents_max.cwiseMax(sampling_time));
    time_events_sum.col(model.poset[*v].event_id) = aux2 + time;
    time_events.col(model.poset[*v].event_id) = time_events_sum.col(model.poset[*v].event_id) - time_parents_max;

    dens += dexp_log(time, model._lambda[model.poset[*v].event_id]) -
      pexp_log(cutoff, model._lambda[model.poset[*v].event_id]);
  }
  return time_events;
}


VectorXd cbn_density_log(const MatrixXd& time, const VectorXd& lambda) {
  unsigned int nrows = time.rows();
  unsigned int ncols = time.cols();
  VectorXd ret(nrows);
  ret.setConstant(lambda.array().log().sum());
  //    ret -= time * lambda;
  for (unsigned int i = 0; i < nrows; ++i){
    for (unsigned j = 0; j < ncols; ++j) {
      ret[i] -= time(i,j)*lambda(j);
    }
  }
  return ret;
}



// ----------------------------------------------------------------------------------------------------------------------
// Sampling implementation

bool matrixDivision1(MyDoubleMatrix& avg_eps, MyDoubleMatrix& avg_eps_current, double t){
  for (int i = 0; i < avg_eps.rows(); ++i) {
    for (int j = 0; j < avg_eps.cols(); ++j) {
      if (fabs(avg_eps(i,j) - avg_eps_current(i,j)) > t){
        return false;
      }
    }
  }
  return true;
}



MyDoubleMatrix add_matrix(MyDoubleMatrix& X, MyDoubleMatrix& Y) {
  for (int i = 0; i < Y.rows(); i++) {
    for (int j = 0; j < Y.cols(); j++) {
      X(i,j) += Y(i,j);
    }
  }
  return X;
}



//' @noRd
     //' @param N number of samples to be drawn
     std::vector<int> rdiscrete_std(const unsigned int N, const VectorXd& weights, Context::rng_type& rng) {
       std::discrete_distribution<int> distribution(weights.data(),weights.data() + weights.size());
       std::vector<int> ret(N);
       for (unsigned int i = 0; i < N; ++i)
         ret[i] = distribution(rng);
       return ret;
     }


VectorXd log_bernoulli_process(const VectorXd& dist, const double eps,
                               const unsigned int p) {

  const unsigned int L = dist.size();
  VectorXd log_prob = VectorXd::Zero(L);

  if (eps == 0) {
    /* NOTE: If all observations are compatible with poset (which can happen
     * because noisy observations can be compatible), then eps can be 0, as
     * well as dist.
     */
    for (unsigned int i = 0; i < L; ++i) {
      if (dist[i] != 0) {
        double log_eps = std::log(eps + DBL_EPSILON);
        double log1m_eps = std::log1p(- eps - DBL_EPSILON);
        log_prob[i] = log_eps * dist[i] + log1m_eps * (p - dist[i]);
      } else {
        log_prob[i] = 0;
      }
    }
  } else {
    double log_eps = std::log(eps);
    double log1m_eps = std::log1p(-eps);//log(1+x)
    log_prob = (log_eps * dist).array() + log1m_eps * (p - dist.array());
  }
  return log_prob;
}


double log_bernoulli_process(const unsigned int dist, const double eps,
                             const unsigned int p) {
  double log_prob = 0.0;

  if (eps == 0) {
    /* NOTE: If all observations are compatible with poset (which can happen
     * because noisy observations can be compatible), then eps can be 0, as
     * well as dist.
     */
    if (dist != 0)
      log_prob = std::log(eps + DBL_EPSILON) * dist +
        std::log1p(- eps - DBL_EPSILON) * (p - dist);
  } else {
    log_prob = std::log(eps) * dist + std::log1p(-eps) * (p - dist);
  }
  return log_prob;
}


//' Generate observations from a given poset and given rates
     //'
     //' @noRd
     //' @param N number of samples
     //' @return returns matrix containing observations
     MatrixXb sample_genotypes(const unsigned int N, const Model& model, MatrixXd& T_events,
                               VectorXd& T_sampling, Context::rng_type& rng) {
       /* Initialization and instantiation of variables */
       const vertices_size_type p = model._size;  // Number of mutations / events
       MatrixXd T_events_sum;
       MatrixXb obs;
       T_events_sum.setZero(N, p);
       obs.setConstant(N, p, false);
       auto id = boost::get(&Event::event_id, model.poset);

       /* Generate occurence times T_events_{j} ~ Exp(lambda_{j}) */
       for (unsigned int j = 0; j < p; ++j)
         T_events.col(j) = rexp(N, model._lambda[j], rng);

       /* Use sampling times when available */
       T_sampling = rexp(N, model._lambda_s, rng);

       VectorXd T_max = VectorXd::Zero(N);
       /* Loop through nodes in topological order */
       for (node_container::const_reverse_iterator v = model.topo_path.rbegin(); v != model.topo_path.rend(); ++v) {
         T_max = VectorXd::Zero(N);
         /* Loop through parents of node v */
         boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
         for (boost::tie(in_begin, in_end) = boost::in_edges(*v, model.poset); in_begin != in_end; ++in_begin) {
           Node u = boost::get(id, source(*in_begin, model.poset));
           T_max = (T_events_sum.col(u).array() > T_max.array()).select(T_events_sum.col(u), T_max);
         }
         T_events_sum.col(model.poset[*v].event_id) = T_events.col(model.poset[*v].event_id) + T_max;
         obs.col(model.poset[*v].event_id) =
           (T_events_sum.col(model.poset[*v].event_id).array() <= T_sampling.array()).select(true, obs.col(model.poset[*v].event_id));
       }
       return obs;
     }


//' Generate observation times from a given poset and given rates
          //'
          //' @noRd
          //' @param N number of samples
          //' @return returns matrix containing occurrence times
          MatrixXd sample_times(const unsigned int N, const Model& model, MatrixXd& T_events, Context::rng_type& rng) {
            /* Initialization and instantiation of variables */
            const vertices_size_type p = model._size;  // Number of mutations / events
            MatrixXd T_events_sum;
            T_events_sum.setZero(N, p);
            /*      generate node id     */
            auto id = boost::get(&Event::event_id, model.poset);

            /* Generate occurence times T_events_{j} ~ Exp(lambda_{j}) */
            for (unsigned int j = 0; j < p; ++j)
              T_events.col(j) = rexp(N, model._lambda[j], rng);

            VectorXd T_max = VectorXd::Zero(N);
            /* Loop through nodes in topological order */
            for (node_container::const_reverse_iterator v = model.topo_path.rbegin(); v != model.topo_path.rend(); ++v) {
              T_max = VectorXd::Zero(N);
              /* Loop through parents of node v */
              boost::graph_traits<Poset>::in_edge_iterator in_begin, in_end;
              for (boost::tie(in_begin, in_end) = boost::in_edges(*v, model.poset); in_begin != in_end; ++in_begin) {
                Node u = boost::get(id, source(*in_begin, model.poset));
                T_max = (T_events_sum.col(u).array() > T_max.array()).select(T_events_sum.col(u), T_max);
              }
              T_events_sum.col(model.poset[*v].event_id) = T_events.col(model.poset[*v].event_id) + T_max;
            }
            return T_events_sum;
          }


MatrixXb generate_genotypes(const MatrixXd& T_events_sum, const Model& model, VectorXd& T_sampling, Context::rng_type& rng) {
  /* Initialization and instantiation of variables */
  const unsigned int N = T_events_sum.rows(); // Number of genotypes to be drawn
  const vertices_size_type p = model._size;  // Number of mutations / events
  MatrixXb obs;
  obs.setConstant(N, p, false);
  T_sampling = rexp(N, model._lambda_s, rng);

  /* Loop through nodes in topological order */
  for (node_container::const_reverse_iterator v = model.topo_path.rbegin(); v != model.topo_path.rend(); ++v)
    obs.col(model.poset[*v].event_id) =
      (T_events_sum.col(model.poset[*v].event_id).array() <= T_sampling.array()).select(true, obs.col(model.poset[*v].event_id));

  return obs;
}


MatrixXb generate_genotypes111(const MatrixXd& T_events_sum, const Model& model, VectorXd& T_sampling) {
  /* Initialization and instantiation of variables */
  const unsigned int N = T_events_sum.rows(); // Number of genotypes to be drawn
  const vertices_size_type p = model._size;  // Number of mutations / events
  MatrixXb obs;

  return obs;
}


//' Compute Hamming distance between two vectors
   //'
   //' @noRd
   //' @return returns Hamming distance
   int hamming_dist(const VectorXi& x, const VectorXi& y) {
     return (x - y).array().abs().sum();
   }


//' Compute Hamming distance between a matrix and a vector row-wise
   //'
   //' @noRd
   //' @return returns a vector containing the Hamming distance
   VectorXi hamming_dist_mat(const MatrixXb& x, const RowVectorXb& y) {
     const int N = x.rows();
     return (x.array() != y.replicate(N, 1).array()).rowwise().count().cast<int>();
     // return (x.rowwise() - y).array().abs().rowwise().sum();
   }


//' Compute Hamming distance between a matrix and a vector row-wise
   //'
   //' @noRd
   //' @return returns a vector containing the Hamming distance
   MyIntMatrix hamming_dist_mat1(const MatrixXb& x, const RowVectorXb& y, const Model& m) {
     const int N = x.rows();
     MyIntMatrix d_pool = MyIntMatrix::Zero(N,m._eventD.size());
     for (int i = 0; i < N; ++i) {
       for (int j = 0; j < m._eventD.size(); ++j) {
         int temp = 0;
         for (int k = m._eventD[j].F-1; k <= m._eventD[j].S-1; ++k) {
           if (x(i,k) != y[k]){
             temp ++;
           }
         }
         d_pool(i,j) = temp;
       }
     }
     return d_pool;
     // return (x.rowwise() - y).array().abs().rowwise().sum();
   }




//' Compute importance weights and (expected) sufficient statistics by
   //' importance sampling
   //'
   //' @noRd
   //' @param genotype a p-dimesional vector corresponding to an observed
   //' genotype - e.g., mutated (1) and non-mutated (0) genes
   //' @param L number of samples
   //' @param model
   //' @param time sampling time
   //' @param sampling variable indicating which proposal to use
   //' @return returns importance weights and (expected) sufficient statistics
   DataImportanceSampling importance_weight(
       const RowVectorXb& genotype, unsigned int L, const Model& model, const int rowI, const std::string& sampling,
       const VectorXd& scale_cumulative, const MyIntMatrix& dist_pool, const MatrixXd& Tdiff_pool,
       const unsigned int neighborhood_dist, Context::rng_type& rng) {

     /* Initialization and instantiation of variables */
     const vertices_size_type p = model._size; // Number of mutations / events
     int eN = model._eventD.size();
     unsigned int reps = 0;
     unsigned int L_aux = 1;
     std::vector<int> nrows(neighborhood_dist+1);
     unsigned int nrows_sum = 0;

     if (sampling == "backward") {
       /* Using the leading and the first k order terms in X */
       for (unsigned int i = 0; i <= neighborhood_dist; ++i) {
         nrows[i] = n_choose_k(p, i);
         nrows_sum += nrows[i];
       }
       reps = L / nrows_sum;
       L = reps * nrows_sum;
     } else if (sampling == "bernoulli") {
       reps = 1;
       L_aux = std::max(L / reps, L_aux);
       L = reps * L_aux;
     }
     DataImportanceSampling importance_sampling(L, p, eN);

     if (sampling == "forward") {
       /* Generate L samples from poset with parameters lambda and lambda_s.
        * In particular, epsilon is zero (default value) - because the idea is to
        * generate samples of X (true genotype)
        */
       MatrixXb samples(L, p);
       VectorXd T_sampling(L);
       samples = sample_genotypes(L, model, importance_sampling.Tdiff, T_sampling,rng);
       importance_sampling.dist = hamming_dist_mat1(samples, genotype, model);
       // importance_sampling.w = pow(model.get_epsilon(), importance_sampling.dist.array()) *
       //   pow(1 - model.get_epsilon(), p - importance_sampling.dist.array());
       for (int i = 0; i < eN; ++i) {
         int length = model._eventD[i].S - model._eventD[i].F + 1;
         for (int j = 0; j < L; ++j) {
           importance_sampling.w[j] *= pow(model._epsilon(rowI,i),1.0*importance_sampling.dist(j,i))
           * pow(1-model._epsilon(rowI,i),1.0*(length-importance_sampling.dist(j,i)));
         }
       }
     }else if (sampling == "add-remove") {
       MatrixXb samples(L, p);
       VectorXd log_prob_Y_X(L);
       VectorXd log_prob_X(L);
       VectorXd log_proposal(L);
       /* Pick a move: add or remove with equal probability
        * Remove: choose event x with prob. ~ sum_{j \in max_path to x}(1/lambda_j)
        * (not applicable if genotype corresponds to the wild type)
        * Add: choose event x with an inverse probability (not applicable if the
        * genotype corresponds to the resistant genotype)
        * 'Stand-still': remain unperturbed (applicable if genotype is compatible
        * with the poset)
        */
       bool compatible = is_compatible(genotype, model);
       unsigned int mutations = genotype.count();
       bool wild_type = (mutations == 0);
       bool resistant_type = (mutations == p);
       std::vector<int> move(L);
       if (compatible) {
         if (wild_type) {
           /* possible moves: add (move 0) or stand-still (move 2) */
           Eigen::Vector3d weights(0.5, 0.0, 0.5);
           move = rdiscrete_std(L, weights, rng);
           log_proposal.setConstant(std::log(0.5));
         } else if (resistant_type) {
           /* possible moves: remove (move 1) or stand-still (move 2) */
           Eigen::Vector3d weights(0.0, 0.5, 0.5);
           move = rdiscrete_std(L, weights, rng);
           log_proposal.setConstant(std::log(0.5));
         } else {
           Eigen::Vector3d weights = Eigen::Vector3d::Ones();
           move = rdiscrete_std(L, weights, rng);
           log_proposal.setConstant(-std::log(3));
         }
       } else {
         /* possible moves: add (move 0) or remove (move 1) */
         Eigen::Vector2d weights = Eigen::Vector2d::Constant(0.5);
         move = rdiscrete_std(L, weights, rng);
         log_proposal.setConstant(std::log(0.5));
       }
       double q_choice;
       VectorXd remove_weight = genotype.select(scale_cumulative.transpose(), 0);
       VectorXd add_weight = scale_cumulative.array().inverse();
       add_weight = genotype.select(0, add_weight.transpose());
       std::vector<int> idx_remove = rdiscrete_std(L, remove_weight, rng);
       std::vector<int> idx_add = rdiscrete_std(L, add_weight, rng);
       for (unsigned int l = 0; l < L; ++l) {
         samples.row(l) = draw_sample(genotype, model, move[l], remove_weight,
                     add_weight, q_choice, idx_remove[l],
                                                     idx_add[l], compatible);
         log_proposal[l] += std::log(q_choice);
       }
       importance_sampling.dist = hamming_dist_mat1(samples, genotype, model);
       VectorXd T_sampling(L);
       /* Generate mutation times based on samples */
       importance_sampling.Tdiff = generate_mutation_times(samples, model, log_proposal, T_sampling, rng);
       log_prob_Y_X.fill(0.0);
       for (int i = 0; i < eN; ++i) {
         int temp = model._eventD[i].S-model._eventD[i].F+1;
         log_prob_Y_X += log_bernoulli_process(importance_sampling.dist.col(i).cast<double>(),
                                               model._epsilon(rowI,i), temp)*temp/(1.0*p);
       }
       log_prob_X = cbn_density_log(importance_sampling.Tdiff, model._lambda);
       importance_sampling.w = (log_prob_Y_X + log_prob_X - log_proposal).array().exp();
     } else if (sampling == "pool") {
       // VectorXd q_prob = pow(model.get_epsilon(), dist_pool.array()) *
       //   pow(1 - model.get_epsilon(), p - dist_pool.array());
       VectorXd q_prob(dist_pool.rows());
       q_prob.fill(1.0);
       for (int i = 0; i < eN; ++i) {
         int length = model._eventD[i].S - model._eventD[i].F + 1;
         for (int j = 0; j < q_prob.rows(); ++j) {
           q_prob(j) *= pow(model._epsilon(rowI,i),1.0*dist_pool(j,i))
           * pow(1.0-model._epsilon(rowI,i),1.0*(length-dist_pool(j,i)));
         }
       }
       /* In the unlikely event that q_prob is 0 for all samples, default to
        * random sampling
        */
       bool random = false;
       if (q_prob.sum() == 0) {
         q_prob.setConstant(1);
         random = true;
       }
       double q_prob_sum = q_prob.sum();
       q_prob /= q_prob_sum;

       /* Draw L samples with replacement and with weights q_prob */
       std::vector<int> idxs_sample = rdiscrete_std(L, q_prob, rng);
       int idx;
       for (unsigned int l = 0; l < L; ++l) {
         idx = idxs_sample[l];
         importance_sampling.dist.row(l) = dist_pool.row(idx);
         importance_sampling.Tdiff.row(l) = Tdiff_pool.row(idx);
       }
       if (random){
         importance_sampling.w.fill(0.0);
         for (int i = 0; i < eN; ++i) {
           int temp = model._eventD[i].S-model._eventD[i].F+1;
           importance_sampling.w += log_bernoulli_process(importance_sampling.dist.col(i).cast<double>(),
                                                          model._epsilon(rowI,i),temp)*temp/(1.0*p);
         }
       }
       else{
         importance_sampling.w.setConstant(q_prob_sum / dist_pool.rows());
       }
     }else if (sampling == "backward") {
       VectorXd log_prob_Y_X(L);
       VectorXd log_prob_X(L);
       VectorXd log_proposal = VectorXd::Zero(L);

       MatrixXb samples = genotype.replicate(nrows_sum, 1);
       MyIntMatrix dist = MyIntMatrix::Zero(nrows_sum,eN);
       VectorXd log_prob_Y_X_aux = VectorXd::Zero(nrows_sum);
       for (int i = 0; i < eN; ++i) {
         int temp = model._eventD[i].S-model._eventD[i].F+1;
         log_prob_Y_X_aux[0] += log_bernoulli_process(0.0,
                                                      model._epsilon(rowI,i),temp)*temp/(1.0*p);
       }
       nrows_sum = 1;
       for (unsigned int i = 1; i <= neighborhood_dist; ++i) {
         neighbors(p, i, samples.block(nrows_sum, 0, nrows[i], p));
         double constTemp = 0;
         for (int j = 0; j < eN; ++j) {
           int temp = model._eventD[i].S-model._eventD[i].F+1;
           dist.col(j).segment(nrows_sum, nrows[i]).setConstant(i);
           constTemp += log_bernoulli_process(i,
                                              model._epsilon(rowI,j),temp)*temp/(1.0*p);
         }
         log_prob_Y_X_aux.segment(nrows_sum, nrows[i]).setConstant(constTemp);
         nrows_sum += nrows[i];
       }

       MatrixXb samples_rep = samples.replicate(reps, 1);
       importance_sampling.dist = dist.replicate(reps, 1);
       log_prob_Y_X = log_prob_Y_X_aux.replicate(reps, 1);
       VectorXd T_sampling(L);

       /* Generate mutation times based on samples */
       importance_sampling.Tdiff =
       generate_mutation_times(samples_rep, model, log_proposal, T_sampling, rng);
       /* Account for incompatible samples */
       Eigen::Matrix<bool, Eigen::Dynamic, 1> incompatible_samples = (importance_sampling.Tdiff.array() < 0.0).rowwise().any();
       log_proposal = log_proposal.array() - std::log(nrows_sum - incompatible_samples.count() / reps);
       log_prob_X = cbn_density_log(importance_sampling.Tdiff, model._lambda);
       importance_sampling.w = (log_prob_Y_X + log_prob_X - log_proposal).array().exp();

       /* Downweight samples that are not feasible / incompatible with current poset */
       importance_sampling.w = incompatible_samples.select(0, importance_sampling.w);
     }else if (sampling == "bernoulli") {
       VectorXd log_prob_Y_X(L);
       VectorXd log_prob_X(L);
       VectorXd log_proposal = VectorXd::Zero(L);

       MatrixXb samples = genotype.replicate(L_aux, 1);
       coin_tossing(samples, model, rowI, rng);

       MyIntMatrix dist = MyIntMatrix::Zero(L_aux,eN);
       VectorXd log_prob_Y_X_aux = VectorXd::Zero(L_aux);
       dist = hamming_dist_mat1(samples, genotype, model);

       MatrixXb samples_rep = samples.replicate(reps, 1);
       importance_sampling.dist = dist.replicate(reps, 1);

       VectorXd T_sampling(L);
       /* Generate mutation times based on samples */
       importance_sampling.Tdiff =
       generate_mutation_times(samples_rep, model, log_proposal, T_sampling, rng);

       log_prob_X = cbn_density_log(importance_sampling.Tdiff, model._lambda);

       importance_sampling.w = (log_prob_X - log_proposal).array().exp();
       /* Account for incompatible samples and downweight them */
       Eigen::Matrix<bool, Eigen::Dynamic, 1> incompatible_samples = (importance_sampling.Tdiff.array() < 0.0).rowwise().any();
       importance_sampling.w = incompatible_samples.select(0, importance_sampling.w);
     }

     return importance_sampling;
   }




#endif //ESTE_LAMBDACAL_H
