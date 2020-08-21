#ifndef FCDLATENT_H
#define FCDLATENT_H

#include <math.h>
#include <RcppArmadillo.h>
#include <mycfunctions.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace R;
using namespace Rcpp;
using namespace std;

/*
 * -----------------------------------------------------------------------------
 * arguments
 * 
 * I          double             n of actors
 * K          double             latent dimension
 * psi        double             agreement probability
 * ksq_eta    double             sender space consensus positions variance
 * sigsq_u    double             sender space latent positions variance (agreement)
 * tausq_u    double             sender space latent positions variance (disagreement)
 * ksq_zeta   double             receiver space consensus positions variance
 * sigsq_v    double             receiver space latent positions variance (agreement)
 * tausq_v    double             receiver space latent positions variance (disagreement)
 * vsisq      double             intercepts variance             
 * nu         double             intercepts mean
 * bet        vec      I x 1     linear predictor intercepts
 * gama       uve      I x 1     sender space agreement indicators
 * eta        mat      I x K     sender space concensus positions
 * U          mat      I x I*K   sender space latent positions
 * xi         uvec     I x 1     reveiver space agreement indicators
 * zeta       mat      I x K     receiver space consensus positions
 * V          mat      I x I*K   receiver space latent positions
 * Z          mat      N x J     auxiliary variables
 * Ymcmc      umat     N x I     data (matrices stored by cols)
 * ----------------------------------------------------------------------------- 
 * samplers
 * 
 * sample_Z
 * sample_U
 * sample_V
 * sample_eta
 * sample_zeta
 * sample_bet
 * sample_nu
 * sample_gama
 * sample_xi
 * sample_psi
 * sample_vsisq
 * sample_sigsq_u
 * sample_sigsq_v
 * sample_tausq_u
 * sample_tausq_v
 * -----------------------------------------------------------------------------
*/

void sample_Z (const double &I, const double &K, const vec &bet, const mat &U, const mat &V, mat &Z, const umat &Ymcmc)
{
        // return: N x J matrix
        int i, ii;
        double N = I*I;
        for (int j = 0; j < I; j++) {
                for (int n = 0; n < N; n++) {
                        i  = n2i (n, I);
                        ii = n2ii(n, I);
                        if (i != ii) {
                                if (Ymcmc.at(n, j) == 1) {
                                        Z.at(n, j) = rtruncnorm(linear_predictor_latent(i, ii, j, K, bet, U, V), 'P');
                                } else {
                                        Z.at(n, j) = rtruncnorm(linear_predictor_latent(i, ii, j, K, bet, U, V), 'N');
                                }
                        }
                }
        }
}

void sample_U (const double &I, const double &K, const double &sigsq_u, const double &tausq_u, const uvec &gama, const mat &eta, const vec &bet, mat &U, const mat &V, const mat &Z)
{
        // return: I x K*I matrix
        double var;
        rowvec tmp(K), v(K), mu(K);
        mat C(K, K), sigma(K, K);
        for (int j = 0; j < I; j++) {
                for (int i = 0; i < I; i++) {
                        tmp = zeros<rowvec>(K);
                        C   = zeros<mat>(K, K); 
                        for (int ii = 0; ii < I; ii++) {
                                if (ii != i) {
                                        v = V(ii, span(j*K, (j + 1)*K - 1));
                                        C += v.t() * v;
                                        tmp += (Z.at(i2n(i, ii, I), j) - bet[j]) * v;
                                }
                        }
                        if (i == j) {
                                if (gama[i] == 0) {
                                        var = tausq_u;
                                        mu  = zeros<rowvec>(K);
                                } else{  // gama[i] == 1
                                        var = sigsq_u;
                                        mu  = eta.row(i);
                                }
                        } else{  // i != j
                                var = sigsq_u;
                                mu  = eta.row(i);
                        }
                        //for (int k = 0; k < K; k++) C(k, k) += 1.0/var;
                        C.diag() += 1.0/var;
                        sigma = inv_sympd(C);
                        U(i, span(j*K, (j + 1)*K - 1)) = rmvnorm(1, sigma * ((mu/var + tmp).t()), sigma);
                }
        }
}

void sample_V (const double &I, const double &K, const double &sigsq_v, const double &tausq_v, const uvec &xi, const mat &zeta, const vec &bet, const mat &U, mat &V, const mat &Z)
{
        // return: I x K*I matrix
        double var;
        rowvec tmp(K), u(K), mu(K);
        mat C(K, K), sigma(K, K);
        for (int j = 0; j < I; j++) {
                for (int ii = 0; ii < I; ii++) {
                        tmp = zeros<rowvec>(K);
                        C   = zeros<mat>(K, K); 
                        for (int i = 0; i < I; i++) {
                                if (ii != i) {
                                        u = U(i, span(j*K, (j+1)*K - 1));
                                        C += u.t() * u;
                                        tmp += (Z.at(i2n(i, ii, I), j) - bet[j]) * u;
                                }
                        }
                        if (ii == j) {
                                if (xi(ii) == 0) {
                                        var = tausq_v;
                                        mu  = zeros<rowvec>(K);
                                } else{  // xi(ii) == 1
                                        var = sigsq_v;
                                        mu  = zeta.row(ii);
                                }
                        } else{  // ii != j
                                var = sigsq_v;
                                mu  = zeta.row(ii);
                        }
                        //for (int k = 0; k < K; k++) C(k, k) += 1.0/var;
                        C.diag() += 1.0/var;
                        sigma = inv_sympd(C);
                        V(ii, span(j*K, (j+1)*K - 1)) = rmvnorm(1, sigma * ((mu/var + tmp).t()), sigma);
                }
        }
}

void sample_eta (const double &I, const double &K, const double &ksq_eta, const double &sigsq_u, const uvec &gama, mat &eta, const mat &U)
{
        // return: I x K matrix
        mat sigma0(K, K, fill::zeros), sigma1(K, K, fill::zeros);
        sigma0.diag() += 1.0/(1.0/ksq_eta + (I - 1.0)/sigsq_u);
        sigma1.diag() += 1.0/(1.0/ksq_eta + I/sigsq_u);
        //for (int k = 0; k < K; k++) {
        //     sigma0(k, k) = 1.0/(1.0/ksq_eta + (I - 1.0)/sigsq_u);
        //     sigma1(k, k) = 1.0/(1.0/ksq_eta + I/sigsq_u);
        //}
        rowvec tmp(K);
        for (int i = 0; i < I; i++) {
                mat sigma = sigma0;
                tmp = zeros<rowvec>(K);
                for (int j = 0; j < I; j++) {
                        if (j != i) {
                                tmp += U(i, span(j*K, (j+1)*K - 1));
                        }
                }
                if (gama[i] == 1) {
                        sigma = sigma1;
                        tmp  += U(i, span(i*K, (i + 1)*K - 1));
                }
                eta.row(i) = rmvnorm(1, sigma * (tmp.t()/sigsq_u), sigma);  
        }
}

void sample_zeta (const double &I, const double &K, const double &ksq_zeta, const double &sigsq_v, const uvec &xi, mat &zeta, const mat &V)
{
        // return: I x K matrix
        mat sigma0(K, K, fill::zeros), sigma1(K, K, fill::zeros);
        sigma0.diag() += 1.0/(1.0/ksq_zeta + (I - 1.0)/sigsq_v);
        sigma1.diag() += 1.0/(1.0/ksq_zeta + I/sigsq_v);
        //for (int k = 0; k < K; k++) {
        //        sigma0(k, k) = 1.0/(1.0/ksq_zeta + (I - 1.0)/sigsq_v);
        //        sigma1(k, k) = 1.0/(1.0/ksq_zeta + I/sigsq_v);
        //}
        rowvec tmp(K);
        for (int i = 0; i < I; i++) {
                mat sigma = sigma0;
                tmp = zeros<rowvec>(K);
                for (int j = 0; j < I; j++) {
                        if (j != i) {
                                tmp += V(i, span(j*K, (j+1)*K - 1));
                        }
                }
                if (xi[i] == 1) {
                        sigma = sigma1;
                        tmp  += V(i, span(i*K, (i+1)*K - 1));
                }
                zeta.row(i) = rmvnorm(1, sigma * (tmp.t()/sigsq_v), sigma);  
        }
}

void sample_bet (const double &I, const double &K, const double &vsisq, const double &nu, vec &bet, const mat &U, const mat &V, const mat &Z)
{
        // return: J x 1 matrix
        int i, ii;
        double N = I*I, tmp, latprod, sigsq = 1.0/(1.0/vsisq + I*(I - 1.0));
        for (int j = 0; j < I; j++) {
                tmp = 0.0;
                for (int n = 0; n < N; n++) {
                        i  = n2i (n, I);
                        ii = n2ii(n, I);
                        if (ii != i) {
                                latprod = 0.0;
                                for (int k = 0; k < K; k++) latprod += U.at(i, j*K + k) * V.at(ii, j*K + k);
                                tmp += (Z.at(n, j) - latprod);
                        }
                }
                bet[j] = R::rnorm(sigsq * (nu/vsisq + tmp), sqrt(sigsq));
        }
}

void sample_nu (const double &J, const double &omesq, double &nu, const double &vsisq, const vec &bet)
{
        double sigsq = 1.0/(1.0/omesq + J/vsisq);
        nu = rnorm(sigsq * (accu(bet)/vsisq), sqrt(sigsq));
}

void sample_gama (const double &I, const double &K, const double &psi, const double &sigsq_u, const double &tausq_u, uvec &gama, const mat &eta, const mat &U)
{
        // return: I x 1 vector
        double logsigsq = log(sigsq_u), logtausq = log(tausq_u), logpsi = log(psi), logd1, d1, d0, logp;
        rowvec u(K);
        for (int i = 0; i < I; i++) {
                u = U(i, span(i*K, (i+1)*K - 1));
                logd1 = -(K/2.0)*(std::log(2.0 * M_PI) + logsigsq) - accu(square(u - eta.row(i)))/(2.0*sigsq_u); 
                d1 = exp(logd1);
                d0 = exp(-(K/2.0)*(std::log(2.0 * M_PI) + logtausq) - accu(square(u))/(2.0*tausq_u));
                // \lop p = \log\psi + \log d1 - \log c
                // c = \psi*d1 + (1-\psi)*d0 (normalizing constant)
                logp = logpsi + logd1 - log(psi*d1 + (1.0 - psi)*d0);
                if (log(R::runif(0.0, 1.0)) < logp) {
                        gama[i] = 1; 
                } else {
                        gama[i] = 0;
                }
        }
}

void sample_xi (const double &I, const double &K, const double &psi, const double &sigsq_v, const double &tausq_v, uvec &xi, const mat &zeta, const mat &V)
{
        // return: I x 1 vector
        double logsigsq = log(sigsq_v), logtausq = log(tausq_v), logpsi = log(psi), logd1, d1, d0, logp;
        rowvec v(K);
        for (int i = 0; i < I; i++) {
                v = V(i, span(i*K, (i+1)*K - 1));
                logd1 = -(K/2.0) * (std::log(2.0 * M_PI) + logsigsq) - accu(square(v - zeta.row(i)))/(2.0*sigsq_v); 
                d1 = exp(logd1);
                d0 = exp(-(K/2.0) * (std::log(2.0 * M_PI) + logtausq) - accu(square(v))/(2.0*tausq_v));
                logp = logpsi + logd1 - log(psi*d1 + (1.0 - psi)*d0);
                if (log(R::runif(0.0, 1.0)) < logp) {
                        xi[i] = 1; 
                } else {
                        xi[i] = 0;
                }
        }
}

void sample_psi (const double &I, const double &a_psi, const double &b_psi, double &psi, const uvec &gama, const uvec &xi)
{
        double tmp = accu(gama + xi);
        psi = R::rbeta(a_psi + tmp, b_psi + 2.0*I - tmp);
}

void sample_vsisq (const double &I, const double &a_vsi, const double &b_vsi, double &vsisq, const double &nu, const vec &bet)
{
        double tmp = 0.0;
        for (int j = 0; j < I; j++) tmp += pow(bet[j] - nu, 2.0);
        vsisq = rinvgamma(a_vsi + I/2.0, b_vsi + tmp/2.0);
}

void sample_sigsq_u (const double &I, const double &K, const double &a_sig, const double &b_sig, double &sigsq_u, const uvec &gama, const mat &eta, const mat &U) 
{
        double tmp = 0.0;
        for (int i = 0; i < I; i++) {
                for (int j = 0; j < I; j++) {
                        if (j != i) {
                                for (int k = 0; k < K; k++) tmp += pow(U.at(i, j*K + k) - eta.at(i, k), 2.0);
                        }
                }
                if (gama[i] == 1) {
                        for (int k = 0; k < K; k++) tmp += pow(U(i, i*K + k) - eta(i, k), 2.0);
                }
        }
        sigsq_u = rinvgamma(a_sig + (I*(I - 1.0) + accu(gama))*K/2.0, b_sig + tmp/2.0);
} 

void sample_sigsq_v (const double &I, const double &K, const double &a_sig, const double &b_sig, double &sigsq_v, const uvec &xi, const mat &zeta, const mat &V) 
{
        double tmp = 0.0;
        for (int i = 0; i < I; i++) {
                for (int j = 0; j < I; j++) {
                        if(j != i) {
                                for (int k = 0; k < K; k++)
                                        tmp += pow(V.at(i, j*K + k) - zeta.at(i, k), 2.0);
                        }
                }
                if (xi[i] == 1) {
                        for (int k = 0; k < K; k++) tmp += pow(V.at(i, i*K + k) - zeta.at(i, k), 2.0);
                }
        }
        sigsq_v = rinvgamma (a_sig + (I*(I - 1.0) + accu(xi))*K/2.0, b_sig + tmp/2.0);
}

void sample_tausq_u (const double &I, const double &K, const double &a_tau, const double &b_tau, double &tausq_u, const uvec &gama, const mat &U) 
{
        double tmp = 0.0;
        for (int i = 0; i < I; i++) {
                if (gama[i] == 0) {
                        for (int k = 0; k < K; k++) tmp += pow(U.at(i, i*K + k), 2.0);
                }
        }
        tausq_u = rinvgamma(a_tau + (I - accu(gama))*K/2.0, b_tau + tmp/2.0);
}

void sample_tausq_v (const double &I, const double &K, const double &a_tau, const double &b_tau, double &tausq_v, const uvec &xi, const mat &V) 
{
        double tmp = 0.0;
        for (int i = 0; i < I; i++) {
                if (xi[i] == 0) {
                        for (int k = 0; k < K; k++) tmp += pow(V.at(i, i*K + k), 2.0);
                }
        }
        tausq_v = rinvgamma(a_tau + (I - accu(xi))*K/2.0, b_tau + tmp/2.0);
}

#endif