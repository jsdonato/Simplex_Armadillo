#include "src/simplex.hpp"
#include <armadillo>
#include <iostream>

int main() {
  arma::mat A = {{1, 1, 1, 0}, 
	         {2, 1, 0, 1}};
  arma::mat b = arma::colvec({12, 16});
  arma::mat c = {-40, -30, 0, 0};
  
  Simplex s(A, b, c);
  s.Run();
  std::cout << "OPTIMAL VALUE:\n" << s.getOptimalValue() << std::endl;
  std::cout << "PRIMAL SOLUTION:\n" << s.getPrimalSolution() << std::endl;   
  std::cout << "DUAL SOLUTION:\n" << s.getDualSolution() << std::endl; 
}
