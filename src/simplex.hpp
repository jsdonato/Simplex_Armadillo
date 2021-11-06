//#include "combination.hpp"
#include <optional>
#include <armadillo>
#include <algorithm>

class Simplex {
public:
  Simplex(arma::mat A_, arma::mat b_, arma::mat c_);

  void Run();

  arma::mat getPrimalSolution() { return x_bar.value(); }

  arma::mat getDualSolution() { return y_bar.value(); }

  double getOptimalValue() { return arma::as_scalar(c * x_bar.value()); }

private:
  const arma::mat A;
  const arma::mat b;
  const arma::mat c; 
  std::optional<arma::mat> x_bar;
  std::optional<arma::mat> y_bar;
  std::vector<size_t> beta;
  std::vector<size_t> eta; 

  void initializePartitions();

  //For a positive integer n and positive integer k, this function generates all subsets of {0,1,...,n-1}
  //of size k and their respective complement with respect to {0,1,...,n-1}.  This is used to find a partition
  //beta (and eta) such that A_beta is invertible.
  std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> combPair(size_t n, size_t k);
  void makeCombPairsHelper(std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& result,
		         std::vector<size_t>& temp, std::vector<size_t>& super_set, size_t n, size_t k, size_t start);



};
