# Simplex Armadillo
Library with implementation of simplex algorithm which uses the Armadillo C++ library.  The implementation of this algorithm is based on the description and derivation provided in the following text (https://github.com/jon77lee/JLee_LinearOptimizationBook/blob/master/JLee.4.0.pdf) which professor Jon Lee utilized and wrote to teach MATH 561 and the University of Michigan.

## Description
The algorithms solves the following two problems as the same time.
```
min cx                     max yb
    Ax=b  (P)              yA<=c  (D)
    x>=0
```
Where `A \in R^mxn`, `b \in R^mx1`, and `c \in R^1xn`.

## Features
| Feature | Description |
|------------|------------|
| `Simplex(arma::mat A, arma::mat b, arma::mat c)` | Creates simplex object |
| `void Run()` | Runs the simplex algorithm on the problem described by matricies `A`, `b`, and `c` inputed to `Simplex` constructor. |
| `arma::mat getPrimalSolution()` | Gets solution associated with `(P)` if solution is not found then it returns an empty `arma::mat` object. |
| `arma::mat getDualSolution()` | Gets solution associated with `(D)` if solution is not found then it returns an empty `arma::mat` object. |
| `double getOptimalValue()` | Returns the optimal value of `cx` and `yb` |

