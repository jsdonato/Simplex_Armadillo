#include "simplex.hpp"

Simplex::Simplex(arma::mat A_, arma::mat b_, arma::mat c_) :
	         A(A_), b(b_), c(c_) {
  if (A.n_cols >= A.n_rows) {
    initializePartitions();
  }
  else {
    std::cout << "ERROR: A.n_rows > A.n_cols\n";
    exit(0);
  }
}

void Simplex::initializePartitions() {
 auto result = combPair(A.n_cols, A.n_rows);
 for (const auto &p : result) {
   if (arma::det(A.cols(arma::conv_to<arma::uvec>::from(p.first))) != 0) {
     beta = p.first; 
     eta = p.second; 
     return;
   }
  }
}

void Simplex::Run() {
  x_bar = arma::mat(A.n_cols, 1); 
 
  while (true) {
    //Find primal solution x_bar and dual solution y_bar associated with 
    //partitions beta and eta.
    arma::mat A_beta = A.cols(arma::conv_to<arma::uvec>::from(beta)); 
    arma::mat A_beta_inv = A_beta.i();
    y_bar = c.cols(arma::conv_to<arma::uvec>::from(beta)) * A_beta_inv;
    arma::mat x_bar_beta = A_beta_inv * b;
    x_bar.zeros();
    for (size_t i = 0; i < beta.size(); i++) { x_bar(beta[i], 0) = x_bar_beta(i, 0); }
    
    //If c_bar_eta = (c-y_barA)_eta>=0 then stop since this implies
    //that x_bar and y_bar are optimal.
    arma::mat c_bar = c - y_bar * A;
    arma::mat c_bar_eta = c_bar.cols(arma::conv_to<arma::uvec>::from(eta));
    bool status1 = arma::all(arma::vectorise(c_bar_eta) >= 0.0);
    if (status1) return;

    //Choose a non basic index eta_j s.t. c_bar_eta_j < 0.
    size_t eta_j;
    size_t j;
    for (size_t i = 0; i < c_bar_eta.n_rows; i++) {
      if (c_bar_eta(0, i) < 0.0) {
        eta_j = eta[i];
	j = i;
	break;
      }
    }

    //If A_bar_eta_j = (A_beta)^-1 * A_eta_j <= 0 then stop since this implies that
    //the primal problem (P) is unbounded and the dual problem (D) is infeasible.
    arma::mat A_bar_eta_j = A_beta_inv * A.col(eta_j);
    bool status2 = arma::all(arma::vectorise(A_bar_eta_j) <= 0.0);
    if (status2) {
      x_bar.reset();
      y_bar.reset();
      return;
    }

    //Let i^* = argmin_{i : A_bar_eta_j(i) > 0}{x_bar_beta(i)/A_bar_eta_j(i)}
    //and replace eta_j with beta_i^*.
    size_t i_star;
    double min = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < A_bar_eta_j.n_rows; i++) {
      if (A_bar_eta_j(i, 0) > 0 && x_bar_beta(i, 0) / A_bar_eta_j(i, 0) < min) {
        min = x_bar_beta(i, 0) / A_bar_eta_j(i, 0);
	i_star = i;

      }
    }
    std::swap(beta[i_star], eta[j]);

  }

}

void Simplex::makeCombPairsHelper(std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>>& result,
		         std::vector<size_t>& temp, std::vector<size_t>& super_set, size_t n, size_t k, size_t start){
  if (k == 0) {
    std::vector<size_t> diff(n - k);
    std::vector<size_t>::iterator it;
    it = std::set_difference(super_set.begin(), super_set.end(), temp.begin(), temp.end(), diff.begin());
    diff.resize(it-diff.begin());
    result.push_back(std::make_pair(temp, diff));
    return;
  }
 
  for (int i = start; i < n; ++i) {
    temp.push_back(i);
    makeCombPairsHelper(result, temp, super_set, n, k - 1, i + 1);
    temp.pop_back();
  }
}

std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> Simplex::combPair(size_t n, size_t k) {
  std::vector<std::pair<std::vector<size_t>, std::vector<size_t>>> result;
  std::vector<size_t> temp;
  std::vector<size_t> super_set;
  for (size_t i = 0; i < n; i++) { super_set.push_back(i); }

  makeCombPairsHelper(result, temp, super_set, n, k, 0);

  return result;


}
