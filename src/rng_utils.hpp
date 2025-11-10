#ifndef RNG_UTILS_HPP
#define RNG_UTILS_HPP
#include <Rcpp.h>

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

// [[Rcpp::depends(BH)]]
#include <boost/core/noncopyable.hpp>
#include <limits>

#include <random>
#include <algorithm>
#include <functional>
#include <array>

using namespace Rcpp;
using Eigen::Map;
using Eigen::VectorXd;
using namespace std;

template <typename RNG_TYPE>
VectorXd rtexp(const unsigned int N, const double rate, const VectorXd& cutoff,
               RNG_TYPE& rng);

template <typename RNG_TYPE>
VectorXd rexp(const unsigned int N, const double rate, RNG_TYPE& rng) {
  double inf = std::numeric_limits<double>::infinity();
  VectorXd aux(N);
  aux.setConstant(inf);
  return rtexp(N, rate, aux, rng);
}

class StdRng: private boost::noncopyable {
public:
  StdRng(unsigned int seed) {
    std::ranlux24_base seeder_rng = std::ranlux24_base(seed);

    // adapted from https://stackoverflow.com/a/15509942
    std::array<int, std::mt19937::state_size> seed_data;
    std::generate_n(seed_data.data(), seed_data.size(), std::ref(seeder_rng));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));

    _rng.seed(seq);
  }

  typedef std::mt19937::result_type result_type;

  result_type operator()() {
    return _rng();
  }

  static constexpr result_type min() {
    return std::mt19937::min();
  }

  static constexpr result_type max() {
    return std::mt19937::max();
  }

  std::mt19937& get_rng() {
    return _rng;
  }

private:
  std::mt19937 _rng;
};

typedef StdRng rng_type;

template <>
VectorXd rtexp<StdRng>(const unsigned int N, const double rate,
                       const VectorXd& cutoff, StdRng& rng);

template <>
VectorXd rtexp<StdRng>(const unsigned int N, const double rate,
                       const VectorXd& cutoff, StdRng& rng) {
  std::uniform_real_distribution<double> distribution(0.0, 1.0);
  VectorXd result(N);

  for (unsigned int i = 0; i < N; ++i) {
    double u = distribution(rng.get_rng());
    result[i] = -std::log1p(u * std::expm1(-cutoff[i] * rate)) / rate;
  }

  return result;
}


#endif
