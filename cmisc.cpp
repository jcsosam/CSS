#include <RcppArmadillo.h>
#include <Rmath.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include "mycfunctions.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace R;
using namespace Rcpp;
using namespace std;


// [[Rcpp::export]]
mat a_iter (double I, double K, const mat &Xmcmc, const mat &beta, const mat &U, const mat &V)
{
        int i, ii;
        double N = Xmcmc.n_rows, P = Xmcmc.n_cols;  // N = I*I
        mat a(N, I, fill::zeros);
        for (int j = 0; j < I; j++) {
                for (int n = 0; n < N; n++) {
                        i = n2i(n, I), ii = n2ii(n, I);
                        if (i != ii) {
                                a.at(n, j) = 0.0;
                                for (int p = 0; p < P; p++) a.at(n, j) += Xmcmc.at(n, p) * beta.at(j, p);
                                for (int k = 0; k < K; k++) a.at(n, j) += U.at(i, j*K + k) * V.at(ii, j*K + k);
                        }
                }
        }
        return a;
}