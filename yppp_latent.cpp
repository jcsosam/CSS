#include <RcppArmadillo.h>
#include <Rmath.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "fcdlatent.h"
#include "mycfunctions.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace R;
using namespace Rcpp;
using namespace std;

umat na_ind(const mat &B)
{
        uword n = B.n_rows, p = B.n_cols, i, j;
        umat out(n, p, fill::zeros);
        for (i = 0; i < n; i++) {
                for (j = 0; j < p; j++) {
                        if (NumericVector::is_na(B.at(i, j))) {
                                out.at(i, j) = 1;
                        }
                }
        }
        return out;
}

vec sample_y (const double &I, const double &K, const umat &na_idx, const vec &bet, const mat &U, const mat &V)
{
        int i, ii, idx = 0;
        double N = I*I;
        vec y(accu(na_idx), fill::zeros);
        for (int j = 0; j < I; j++) {
                for (int n = 0; n < N; n++) {
                        if (na_idx.at(n,j) == 1) {
                                i = n2i(n, I), ii = n2ii(n, I);  
                                y[idx] = R::rbinom (1, R::pnorm(linear_predictor_latent(i, ii, j, K, bet, U, V), 0.0, 1.0, 1, 0));
                                ++idx;
                        }
                }
        }
        return y;
}

// [[Rcpp::export]]
vec yppp_latent (mat Yna, const double &K, const int &nburn, const int &nsams, const int &nskip)
{
        double S = nburn + nskip*nsams;  // n MCMC iterations
        double N = Yna.n_rows;           // n dyadic measurements, N = I*I
        double I = Yna.n_cols;           // n actors

        //prior elicitation
        double a_psi = 1.0;
        double b_psi = 1.0;
        double a_vsi = 2.0;
        double b_vsi = (a_vsi - 1.0)/4.0;
        double omesq = 0.25;
        double a_sig = 2.0;
        double b_sig = (a_sig - 1.0)/(2.0*sqrt(2.0*K));
        double a_tau = 2.0;
        double b_tau = (a_tau - 1.0)/(sqrt(2.0*K));
        double kapsq = 1.0/(2.0*sqrt(2.0*K));
        
        //parameter initialization from the prior
        vec zero_K(K, fill::zeros); // zero vector size K
        mat I_K(K, K, fill::eye);   // identity matrix size K
        
        double psi      = R::beta(a_psi, b_psi)  ;
        double vsisq    = rinvgamma(a_vsi, b_vsi);
        double sigsq_u  = rinvgamma(a_sig, b_sig);
        double sigsq_v  = rinvgamma(a_sig, b_sig);
        double tausq_u  = rinvgamma(a_tau, b_tau);
        double tausq_v  = rinvgamma(a_tau, b_tau);  
        double nu = rnorm(0.0, sqrt(omesq));
        uvec gama(I); 
        uvec xi(I); 
        vec bet(I);
        for (int j = 0; j < I; j++) {
                gama[j] = R::rbinom(1, psi);
                xi[j]   = R::rbinom(1, psi);
                bet[j]  = R::rnorm(nu, sqrt(vsisq));
        }         
        mat eta  = rmvnorm(I, zero_K, kapsq*I_K);   
        mat zeta = rmvnorm(I, zero_K, kapsq*I_K);  
        mat U(I, K*I, fill::zeros);
        mat V(I, K*I, fill::zeros);
        for (int i = 0; i < I; i++) {
                for (int j = 0; j < I; j++) {
                        if (i != j) {
                                U(i, span(j*K, (j + 1)*K - 1)) = rmvnorm (1,  eta.row(i).t(), sigsq_u*I_K);
                                V(i, span(j*K, (j + 1)*K - 1)) = rmvnorm (1, zeta.row(i).t(), sigsq_v*I_K);
                        } else {  // i == j
                                // U
                                if (gama(i) == 1) {
                                        U(i, span(i*K, (i + 1)*K - 1)) = rmvnorm (1, eta.row(i).t(), sigsq_u*I_K);
                                } else {
                                        U(i, span(i*K, (i + 1)*K - 1)) = rmvnorm (1, zero_K, tausq_u*I_K);
                                }
                                // V
                                if (xi(i) == 1) {
                                        V(i, span(i*K, (i + 1)*K - 1)) = rmvnorm (1, zeta.row(i).t(), sigsq_v*I_K);
                                } else {
                                        V(i, span(i*K, (i + 1)*K - 1)) = rmvnorm (1, zero_K, tausq_v*I_K);
                                } 
                        }
                }
        }
        mat Z(N, I, fill::zeros);  // no need to initialize Z
        
        // posterior predictive probabilities
        umat na_idx = na_ind(Yna);
        double nmiss = accu(na_idx);
        vec y_ppp(nmiss, fill::zeros);
        vec y(nmiss);
        
        umat Ymcmc(N, I);
        for (int j = 0; j < I; ++j) {
                for (int n = 0; n < N; n++) {
                        if (na_idx.at(n, j) == 0) {
                                Ymcmc.at(n, j) = Yna.at(n, j);
                        } 
                }
        }
        
        //chain
        for (int s = 0; s < S; s++) {
                //update missing values
                y = sample_y(I, K, na_idx, bet, U, V);
                
                int h = 0;
                for (int j = 0; j < I; j++) {
                        for (int n = 0; n < N; n++) {
                                if (NumericVector::is_na(Ymcmc.at(n, j))) { 
                                        Ymcmc.at(n, j) = y[n];
                                        ++h;
                                }
                        }
                }
                
                //update
                sample_Z   (I, K, bet, U, V, Z, Ymcmc);
                sample_U   (I, K, sigsq_u, tausq_u, gama, eta , bet, U, V, Z);
                sample_V   (I, K, sigsq_v, tausq_v, xi,   zeta, bet, U, V, Z);
                sample_eta (I, K, kapsq, sigsq_u, gama, eta,  U);
                sample_zeta(I, K, kapsq, sigsq_v, xi  , zeta, V);
                sample_gama(I, K, psi, sigsq_u, tausq_u, gama, eta,  U);
                sample_xi  (I, K, psi, sigsq_v, tausq_v, xi,   zeta, V);
                sample_psi (I, a_psi, b_psi, psi, gama, xi);
                sample_sigsq_u(I, K, a_sig, b_sig, sigsq_u, gama, eta,  U);
                sample_sigsq_v(I, K, a_sig, b_sig, sigsq_v, xi,   zeta, V);
                sample_tausq_u(I, K, a_tau, b_tau, tausq_u, gama, U);
                sample_tausq_v(I, K, a_tau, b_tau, tausq_v, xi,   V);
                sample_bet(I, K, vsisq, nu, bet, U, V, Z);
                sample_nu(I, omesq, nu, vsisq, bet);
                sample_vsisq(I, a_vsi, b_vsi, vsisq, nu, bet);
                
                // ppp
                if ((s + 1 > nburn) && ((s + 1)%nskip == 0)) y_ppp += y / (double)nsams;
        }
        
        return y_ppp ;
}