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
void mcmc_latent_cpp (const umat   &Ymcmc,
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
        nam  =   "vsisq_chain"; full = mypaste0(path_outs, nam); ofstream   vsisq_chain;   vsisq_chain.open(full);
        nam  =      "nu_chain"; full = mypaste0(path_outs, nam); ofstream      nu_chain;      nu_chain.open(full);
        nam  =     "bet_chain"; full = mypaste0(path_outs, nam); ofstream     bet_chain;     bet_chain.open(full);
        nam  = "sigsq_u_chain"; full = mypaste0(path_outs, nam); ofstream sigsq_u_chain; sigsq_u_chain.open(full);
        nam  = "sigsq_v_chain"; full = mypaste0(path_outs, nam); ofstream sigsq_v_chain; sigsq_v_chain.open(full);
        nam  = "tausq_u_chain"; full = mypaste0(path_outs, nam); ofstream tausq_u_chain; tausq_u_chain.open(full); 
        nam  = "tausq_v_chain"; full = mypaste0(path_outs, nam); ofstream tausq_v_chain; tausq_v_chain.open(full);
        nam  =     "psi_chain"; full = mypaste0(path_outs, nam); ofstream     psi_chain;     psi_chain.open(full);
        nam  =    "gama_chain"; full = mypaste0(path_outs, nam); ofstream    gama_chain;    gama_chain.open(full);
        nam  =      "xi_chain"; full = mypaste0(path_outs, nam); ofstream      xi_chain;      xi_chain.open(full);
        nam  =     "eta_chain"; full = mypaste0(path_outs, nam); ofstream     eta_chain;     eta_chain.open(full);
        nam  =    "zeta_chain"; full = mypaste0(path_outs, nam); ofstream    zeta_chain;    zeta_chain.open(full);
        nam  =       "U_chain"; full = mypaste0(path_outs, nam); ofstream       U_chain;       U_chain.open(full);
        nam  =       "V_chain"; full = mypaste0(path_outs, nam); ofstream       V_chain;       V_chain.open(full);
        nam  =      "ll_chain"; full = mypaste0(path_outs, nam); ofstream      ll_chain;      ll_chain.open(full);
        nam  =      "mcmc_log"; full = mypaste0(path_outs, nam); ofstream      mcmc_log;      mcmc_log.open(full);
        
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
                if ((s > nburn) && (s % nskip == 0)) {
                        // save scalars
                        ll_chain      << ll_iter_latent(I, K, bet, U, V, Ymcmc) << "\n";
                        psi_chain     << psi     << "\n";
                        sigsq_u_chain << sigsq_u << "\n";
                        sigsq_v_chain << sigsq_v << "\n";
                        tausq_u_chain << tausq_u << "\n";
                        tausq_v_chain << tausq_v << "\n";
                        nu_chain      << nu      << "\n";
                        vsisq_chain   << vsisq   << "\n";
                        // save vectors and matrices
                        for (int i = 0; i < I; i++) {
                                gama_chain << gama[i] << " ";
                                xi_chain   << xi[i]   << " ";
                                bet_chain  << bet[i]  << " ";
                        }
                        for (int k = 0; k < K; k++) {
                                for (int i = 0; i < I; i++) {
                                        eta_chain  << eta.at (i, k) << " ";
                                        zeta_chain << zeta.at(i, k) << " ";
                                }
                        }
                        for (int h = 0; h < K*I; ++h) {
                                for (int i = 0; i < I; i++) {
                                        U_chain << U.at(i, h) << " ";
                                        V_chain << V.at(i, h) << " ";
                                }
                        }
                        gama_chain << "\n";
                        xi_chain   << "\n";
                        bet_chain  << "\n";
                        eta_chain  << "\n";
                        zeta_chain << "\n";
                        U_chain    << "\n";
                        V_chain    << "\n";
                }
                // display
                if (s % ndisp == 0) mcmc_log << round_to_digits(100.0*s/S, 2.0) << "% completed" << "\n";
        }
        ////////////////////////////// chain end ///////////////////////////////
        
        mcmc_log << "///////////////////////////////" << endl << 
                    "      End of the algorithm     " << endl <<
                    "///////////////////////////////" << endl;
        
        // close files
        U_chain.close();
        V_chain.close();
        eta_chain.close();
        zeta_chain.close();
        gama_chain.close();
        xi_chain.close();
        psi_chain.close();
        sigsq_u_chain.close();
        sigsq_v_chain.close();
        tausq_u_chain.close();
        tausq_v_chain.close();
        bet_chain.close();
        nu_chain.close();
        vsisq_chain.close();
        ll_chain.close();
        mcmc_log.close();
}  // end MCMC function