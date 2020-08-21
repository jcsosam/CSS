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

/* ----------------------------------------------------------------------------- 
 * arguments
 * 
 * Ymcmc      mat      N x I     data (matrices stored by cols)
 * K          double             latent dimension
 * a_psi      double             psi   hyper-parameter 
 * b_psi      double             psi   hyper-parameter 
 * a_vsi      double             vsisq hyper-parameter 
 * b_vsi      double             vsisq hyper-parameter 
 * omesq      double             omesq hyper-parameter 
 * a_sig      double             sigsq_u/sigsq_v hyper-parameter 
 * b_sig      double             sigsq_u/sigsq_v hyper-parameter 
 * a_tau      double             tausq_u/tausq_v hyper-parameter 
 * b_tau      double             tausq_u/tausq_v hyper-parameter 
 * kapsq      double             eta/zeta hyper-parameter 
 * psi        double             agreement probability
 * sigsq_u    double             sender space latent positions variance (agreement)
 * tausq_u    double             sender space latent positions variance (disagreement)
 * sigsq_v    double             receiver space latent positions variance (agreement)
 * tausq_v    double             receiver space latent positions variance (disagreement)
 * vsisq      double             intercepts variance             
 * nu         double             intercepts mean
 * bet        vec      I x 1     linear predictor intercepts
 * gama       uvec     I x 1     sender space agreement indicators
 * eta        mat      I x K     sender space concensus positions
 * U          mat      I x I*K   sender space latent positions
 * xi         uvec     I x 1     reveiver space agreement indicators
 * zeta       mat      I x K     receiver space consensus positions
 * V          mat      I x I*K   receiver space latent positions
 * Z          mat      N x J     auxiliary variables 
 * nsam       int                n samples
 * nburn      int                burn-in
 * nskip      int                thinning
 * ndisp      int                display info every ndisp iterations
 * path_outs  string             output path
 * ----------------------------------------------------------------------------- 
 */

// [[Rcpp::export]]
void mcmc_ll_latent_cpp (const umat   &Ymcmc,
                           const double &K, 
                           const double &a_psi, 
                           const double &b_psi, 
                           const double &a_vsi, 
                           const double &b_vsi, 
                           const double &omesq, 
                           const double &a_sig, 
                           const double &b_sig, 
                           const double &a_tau, 
                           const double &b_tau, 
                           const double &kapsq,
                           double psi, 
                           double vsisq, 
                           double sigsq_u, 
                           double sigsq_v, 
                           double tausq_u, 
                           double tausq_v, 
                           double nu, 
                           vec    bet, 
                           uvec   gama, 
                           uvec   xi, 
                           mat    eta, 
                           mat    zeta, 
                           mat    U, 
                           mat    V,
                           const int &nburn, 
                           const int &nsams, 
                           const int &nskip, 
                           const int &ndisp, 
                           string path_outs)
{
        int S    = nburn + nskip*nsams;  // n MCMC iterations
        double N = Ymcmc.n_rows;         // n dyadic measurements, N = I*I
        double I = Ymcmc.n_cols;         // n actors
        mat Z(N, I, fill::zeros);        // no need to initialize Z
        
        //write samples: opening files
        char* full;
        string nam;
        nam  = "ll_chain"; full = mypaste0(path_outs, nam); ofstream ll_chain; ll_chain.open(full);
        nam  = "mcmc_log"; full = mypaste0(path_outs, nam); ofstream mcmc_log; mcmc_log.open(full);
        
        // info
        mcmc_log << "///////////////////////////////" << endl <<
                    "           MCMC info           " << endl << 
                    "latent model  for CSS"           << endl <<
                    "***Directed*** perceptions"      << endl << 
                    "latent dimsension, K = " << K    << endl <<
                    "Number of actors,  I = " << I    << endl <<
                    "nburn = " << nburn               << endl <<
                    "nskip = " << nskip               << endl <<
                    "nsams = " << nsams               << endl <<
                    "niter = " << S                   << endl <<
                    "///////////////////////////////" << endl;
        
        ///////////////////////////////// chain ////////////////////////////////
        for (int s = 1; s <= S; s++) {
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

                // save samples
                if ((s > nburn) && (s % nskip == 0)) ll_chain  << ll_iter_latent(I, K, bet, U, V, Ymcmc)  << "\n";
                // display
                if (s % ndisp == 0) mcmc_log << round_to_digits(100.0*s/S, 2.0) << "% completed" << "\n";
        }
        ////////////////////////////// chain end ///////////////////////////////
        
        mcmc_log << "///////////////////////////////" << endl << 
                    "      End of the algorithm     " << endl <<
                    "///////////////////////////////" << endl;
        
        // close files
        ll_chain.close();
        mcmc_log.close();
}  // end MCMC function