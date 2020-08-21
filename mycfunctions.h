#ifndef MYCFUNCTIONS_H
#define MYCFUNCTIONS_H

#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace R;
using namespace Rcpp;
using namespace std;

const double log2pi = std::log(2.0 * M_PI);

double round_to_digits (double value, int digits)
{
        // otherwise it will return 'nan' due to the log10() of zero
        if (value == 0.0) 
                return 0.0;
        double factor = pow(10.0, digits - ceil(log10(fabs(value))));
        return round(value*factor)/factor;   
}

char* mypaste0 (string path, string name)
{
        // pastes path and name together
        stringstream strname;
        strname << path << name << ".txt";
        string fullname = strname.str();
        string::iterator p = fullname.begin();
        char* chr = &(*p);
        return( chr );
}

char* mypaste (double K, string path, string name)
{
        stringstream strname;
        strname << path << name << K << ".txt";
        string fullname = strname.str();
        string::iterator p = fullname.begin();
        char* chr = &(*p);
        return chr;
}

int i2n (const int &i, const int &ii, const double &I)
{
        return ii*(int)I + i;
}

int n2i (const int &n, const double &I)
{
        return n%(int)I;
}

int n2ii (const int &n, const double &I)
{
        return floor((double)n/I);
}

int nj2h (const int &n, const int &j, const double &N)
{
        return  n + j*(int)N;
}

mat vec2mat (const double &I, const rowvec &v)
{
        // returns the matrix version with I rows of a given vector
        int H = v.n_elem/I;
        mat M(I, H);
        for (int h = 0; h < H; h++) {
                for (int i = 0; i < I; i++) {
                        M(i, h) = v(i + h*I);  // fills by cols
                }
        }
        return M;
}

mat rmvnorm (const int &n, const vec &mu, const mat &sigma) 
{ 
        // random draw from a multivariate normal distribution
        // return: n x k matrix
        // n: number of samples, k: dimension 
        mat Y = randn(n, sigma.n_cols);
        return repmat(mu, 1, n).t() + Y * chol(sigma);
}

double rinvgamma (const double &alpha, const double &beta)
{
        return 1.0/R::rgamma(alpha, 1.0/beta);
}

double rtruncnorm (const double mu, char caso)
{
        // sample from a TN with mean mu and sd 1 (just two cases)
        double aa = R::pnorm(-mu, 0.0, 1.0, 1, 0), uu = R::runif(0.0, 1.0), out = 0.0;
        switch (caso) {
        case 'N' :  // draw from TN (mu, sd = 1, -Inf, 0)
                out = R::qnorm (aa*uu, 0.0, 1.0, 1, 0) + mu;
                if (out == -datum::inf) { 
                        out = -8.21; 
                }
                break;
        case 'P' :  // draw from TN (mu, sd = 1, 0, Inf)
                out = R::qnorm (aa*(1 - uu) + uu, 0.0, 1.0, 1, 0) + mu;
                if (out == datum::inf) { 
                        out = 8.21; 
                }
                break;
        }
        return out;
}

mat rwishart (const double &df, const mat &S) 
{
        uword m = S.n_rows;
        mat Z(m,m);
        for (uword i = 0; i < m; i++) Z(i, i) = sqrt(R::rchisq(df - i));
        for (uword j = 0; j < m; j++) for (uword i = j + 1; i < m; i++) Z(i, j) = R::rnorm(0.0, 1.0);
        mat C = trimatl(Z).t() * chol(S);
        return( C.t()*C );
}

mat riwishart (const double &df, const mat &S)
{
        return( rwishart(df, S.i()).i() );
}

double dmvnorm (const rowvec &x, const rowvec &mean, const mat &sigma, const bool &logd) 
{ 
        // x: k x 1 row vector
        // k: dimension 
        // density from a multivariate normal distribution
        // return: double
        int xdim = x.n_cols;
        arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
        double rootisum = arma::sum(log(rooti.diag()));
        double constants = -(static_cast<double>(xdim)/2.0)*log2pi;
        arma::vec z = rooti*arma::trans(x - mean);    
        double out = constants - 0.5 * arma::sum(z%z) + rootisum;     
        if(logd == false) out = exp(out);
        return out;
}

double dinvgamma (const double &x, const double &a, const double &b, const bool &logd)
{
        double out = a*log(b) - (a + 1.0)*log(x) - b/x - log(R::gammafn(a));
        if (logd == false) out = exp(out);
        return out;
}

double linear_predictor_latent (const int &i, const int &ii, const int &j, const double &K, const vec &beta, const mat &U, const mat &V)
{       
     // linear predictor latent model
     double out = beta[j];
     for (int k = 0; k < K; k++) out += U.at(i, j*K + k) * V.at(ii, j*K + k);
     return out;
}

double linear_predictor_eigen (const int &i, const int &ii, const int &j, const double &K, const vec &beta, const vec &lambda, const mat &U, const mat &V)
{       
     // linear predictor eigen model
     double out = beta[j];
     for (int k = 0; k < K; k++) out += lambda[k] * U.at(i, j*K + k) * V.at(ii, j*K + k);
     return out;
}

double ll_iter_latent (const double &I, const double &K, const vec &beta, const mat &U, const mat &V, const umat &Ymcmc) 
{
        // computes the log-likelihood in a given iteration of the MCMC
        int i, ii;
        double N = I*I, out = 0.0;
        for (int j = 0; j < I; ++j) {
                for (int n = 0; n < N; n++) {
                        i  = n2i (n, I);
                        ii = n2ii(n, I);
                        if (i != ii) {
                                if (Ymcmc.at(n, j) == 1) {
                                        out += R::pnorm(linear_predictor_latent(i, ii, j, K, beta, U, V), 0.0, 1.0, 1, 1);
                                } else {
                                        out += R::pnorm(linear_predictor_latent(i, ii, j, K, beta, U, V), 0.0, 1.0, 0, 1);
                                }
                        }
                }
        }
        return out;
}

double ll_iter_eigen (const double &I, const double &K, const vec &beta, const vec &lambda, const mat &U, const mat &V, const umat &Ymcmc) 
{
     // computes the log-likelihood in a given iteration of the MCMC
     int i, ii;
     double N = I*I, out = 0.0;
     for (int j = 0; j < I; ++j) {
          for (int n = 0; n < N; n++) {
               i  = n2i (n, I);
               ii = n2ii(n, I);
               if (i != ii) {
                    if (Ymcmc.at(n, j) == 1) {
                         out += R::pnorm(linear_predictor_eigen(i, ii, j, K, beta, lambda, U, V), 0.0, 1.0, 1, 1);
                    } else {
                         out += R::pnorm(linear_predictor_eigen(i, ii, j, K, beta, lambda, U, V), 0.0, 1.0, 0, 1);
                    }
               }
          }
     }
     return out;
}

double ll_iter_swartz (const double &I, const double &mu, const vec &al, const vec &be, const vec &ga, const mat &albe, const mat &alga, const umat &Ymcmc) 
{
        // return: double
        int i, ii;
        double N = I*I, out = 0.0;
        for (int j = 0; j < I; j++) {
                for (int n = 0; n < N; n++) {
                        i = n2i(n, I), ii = n2ii (n, I);
                        if (i != ii) {
                                if (Ymcmc.at(n, j) == 1) {
                                        out += R::pnorm(mu + al[i] + be[ii] + ga[j] + albe.at(i, ii) + alga.at(i, j) + alga.at(ii, j), 0.0, 1.0, 1, 1);
                                } else {
                                        out += R::pnorm(mu + al[i] + be[ii] + ga[j] + albe.at(i, ii) + alga.at(i, j) + alga.at(ii, j), 0.0, 1.0, 0, 1);
                                }
                        }
                }
        }
        return out;
}

mat PHI3_mat (const double &phi_1, const double &phi_2, const double &phi_3)
{
        // return: 3 x 3 matrix
        mat out(3, 3);
        out << 1.0   << phi_1 << phi_2 <<  endr
            << phi_1 << 1.0   << phi_3 <<  endr
            << phi_2 << phi_3 << 1.0   <<  endr;
        return out; 
}

mat PHI_mat (const double &phi_1, const double &phi_2, const double &phi_3, const double &phi_4)
{
        // return: 4 x 4 matrix
        mat PHI(4, 4);
        PHI << 1.0   << phi_1 << phi_2 << phi_3 << endr
            << phi_1 << 1.0   << phi_3 << phi_2 << endr
            << phi_2 << phi_3 << 1.0   << phi_4 << endr
            << phi_3 << phi_2 << phi_4 << 1.0   << endr;
        return PHI; 
}

mat SIG_2_mat (const double &sigsq_albe, const double &sigsq_alga, const vec &phi)
{
        // return: 4 x 4 matrix
        mat SIG_2(4, 4);
        // Check whether determinants are positive or not
        double det1 = det(PHI3_mat(phi[0], phi[1], phi[2]));
        double det2 = det(PHI_mat(phi[0], phi[1], phi[2], phi[3]));
        double cons1, cons2;
        if (det1 < 0.0) { cons1 = 0.0; } else { cons1 = 1.0; }
        if (det2 < 0.0) { cons2 = 0.0; } else { cons2 = 1.0; }
        // Create the off-diagonal entries of Sig2
        // If the generated set of phi1, phi2, phi3 and phi4 values don't have both 
        // det1 > 0 and det2 > 0, the off-diagonal elements of Sig2 are all set equal to zero 
        // main diagonal  
        SIG_2.at(0, 0) = sigsq_albe;
        SIG_2.at(1, 1) = sigsq_albe;
        SIG_2.at(2, 2) = sigsq_alga;
        SIG_2.at(3, 3) = sigsq_alga;
        // off-diagonal
        SIG_2.at(0, 1) = cons1 * cons2 * sigsq_albe * phi[0];
        SIG_2.at(1, 0) = SIG_2.at(0, 1);
        SIG_2.at(0, 2) = cons1 * cons2 * sqrt(sigsq_albe * sigsq_alga) * phi[1];
        SIG_2.at(2, 0) = SIG_2.at(0, 2);
        SIG_2.at(0, 3) = cons1 * cons2 * sqrt(sigsq_albe * sigsq_alga) * phi[2];
        SIG_2.at(3, 0) = SIG_2.at(0, 3);
        SIG_2.at(1, 2) = cons1 * cons2 * sqrt(sigsq_albe * sigsq_alga) * phi[2];
        SIG_2.at(2, 1) = SIG_2.at(1, 2);
        SIG_2.at(1, 3) = cons1 * cons2 * sqrt(sigsq_albe * sigsq_alga) * phi[1];
        SIG_2.at(3, 1) = SIG_2.at(1, 3);
        SIG_2.at(2, 3) = cons1 * cons2 * sigsq_alga * phi[3];
        SIG_2.at(3, 2) = SIG_2.at(2, 3);
        return SIG_2;
}

vec phi_prior (const double &a_phi, const double &b_phi)
{
        // return: 4 x 1 vector
        double phi_1, phi_2, phi_3, phi_4;
        mat PHI3(3, 3), PHI(4, 4);
        do {
                do {
                        phi_1 = R::runif(a_phi, b_phi);
                        phi_2 = R::runif(a_phi, b_phi);
                        phi_3 = R::runif(a_phi, b_phi);
                        PHI3  = PHI3_mat(phi_1, phi_2, phi_3); 
                } while (det(PHI3) < 0.0);
                phi_4 = R::runif(a_phi, b_phi);
                PHI  = PHI_mat(phi_1, phi_2, phi_3, phi_4); 
        } while (det(PHI) < 0.0);
        vec phi(4);
        phi << phi_1 << endr
            << phi_2 << endr
            << phi_3 << endr
            << phi_4 << endr;
        return phi;
}

vec mytunning (const int &s, double ntun, const double &ntun0, double delta, const double &del0, const double &eps0, const double &mix, const double &npar)
{
        // tunning paramter calibration
        // return: 2 x 1 vector: [delta, ntun]^T
        double mix_rate = mix/((double)(s + 1)*npar), mix_diff = mix_rate - del0, tmp, cont;
        vec out(2);
        if ((s + 1) % (int)ntun == 0) {
                if (abs(mix_diff) > eps0) {
                        ntun = ntun0;
                        tmp  = delta;
                        cont = 0.0;
                        do {
                                cont += 1.0;
                                tmp   = delta + (0.01/cont) * mix_diff;
                        } while (!(tmp > 0.0 || cont == 1e+6));
                        delta = tmp;
                } else {
                        ntun = ntun + ntun0;
                }
        }
        out[0] = delta;
        out[1] = ntun;
        return out;
}

#endif