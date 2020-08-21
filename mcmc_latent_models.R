##------------------------------------------------------------------------------
##------------------------------------------------------------------------------
##                              LATENT MODEL
##------------------------------------------------------------------------------
##------------------------------------------------------------------------------

prior_init_latent <- function(I, K, a_psi, b_psi, a_vsi, b_vsi, omesq, a_sig, b_sig, a_tau, b_tau, kapsq)
{
        psi     <- rbeta(1, a_psi, b_psi)
        vsisq   <- 1/rgamma(1, a_vsi, b_vsi)
        sigsq_u <- 1/rgamma(1, a_sig, b_sig)
        sigsq_v <- 1/rgamma(1, a_sig, b_sig)
        tausq_u <- 1/rgamma(1, a_tau, b_tau)
        tausq_v <- 1/rgamma(1, a_tau, b_tau)
        nu      <- rnorm(1, 0, sqrt(omesq))
        bet     <- rnorm(I, nu, sqrt(vsisq))
        gama    <- rbinom(I, 1, psi)
        xi      <- rbinom(I, 1, psi)
        eta     <- matrix(rnorm(I*K, 0, sqrt(kapsq)), I, K)
        zeta    <- matrix(rnorm(I*K, 0, sqrt(kapsq)), I, K)
        U       <- matrix(0, I, I*K) 
        V       <- matrix(0, I, I*K)
        for (i in 1:I) {
                for (j in 1:I) {
                        if (i != j) {
                                U[i, ((j-1)*K+1):(j*K)] = rnorm(K, eta [i,], sqrt(sigsq_u))
                                V[i, ((j-1)*K+1):(j*K)] = rnorm(K, zeta[i,], sqrt(sigsq_v))
                        } else {
                                if (gama[i] == 1) {
                                        U[i, ((j-1)*K+1):(j*K)] = rnorm(K, eta[i,],  sqrt(sigsq_u))
                                } else {
                                        U[i, ((j-1)*K+1):(j*K)] = rnorm(K, rep(0,K), sqrt(tausq_u))
                                }
                                if (xi[i] == 1) {
                                        V[i, ((j-1)*K+1):(j*K)] = rnorm(K, zeta[i,], sqrt(sigsq_v))
                                } else {
                                        V[i, ((j-1)*K+1):(j*K)] = rnorm(K, rep(0,K), sqrt(tausq_v))
                                } 
                        }
                }
        }
        list(psi     = psi    ,
             vsisq   = vsisq  ,
             sigsq_u = sigsq_u,
             sigsq_v = sigsq_v,
             tausq_u = tausq_u,
             tausq_v = tausq_v,
             nu      = nu     ,
             bet     = bet    ,
             gama    = gama   ,
             xi      = xi     ,
             eta     = eta    ,
             zeta    = zeta   ,
             U       = U      ,
             V       = V      
        )
}

mcmc_ic_latent <- function(loc, case, model, dataset, K, nburn, nsams, nskip, ndisp, seed)
{
        ## paths ---------------------------------------------------------------
        path_code <- paste0(loc, "code/")
        path_data <- paste0(loc, "data/")
        path_outs <- paste0(loc, "outs/", dataset, "_", case, "_", model, "/")
        if (!dir.exists(path_outs)) dir.create(path_outs)
        path_outs <- paste0(loc, "outs/", dataset, "_", case, "_", model, "/", "K_", K, "/")
        if (!dir.exists(path_outs)) dir.create(path_outs)
        if (length(dir(path_outs)) > 0) unlink(paste0(path_outs, "*"))
        
        ## data ----------------------------------------------------------------
        load(paste0(path_data, "data_", dataset, ".RData"))
        
        ## libs and source -----------------------------------------------------
        suppressPackageStartupMessages(require(Rcpp))
        sourceCpp(paste0(path_code, "mcmc_ic_", model, ".cpp"))
        
        ## hyper-parameter elicitation -----------------------------------------
        a_psi <- 1
        b_psi <- 1
        a_vsi <- 2
        b_vsi <- (a_vsi-1)/4
        omesq <- 0.25
        a_sig <- 2
        b_sig <- (a_sig-1)/(2*sqrt(2*K))
        a_tau <- 2
        b_tau <- (a_tau-1)/(sqrt(2*K))
        kapsq <- 1/(2*sqrt(2*K))
        hyps  <- data.frame(a_psi, b_psi, a_vsi, b_vsi, omesq, a_sig, b_sig, a_tau, b_tau, kapsq)
        nhyps <- length(hyps)
        
        ## prior initialization ------------------------------------------------
        set.seed(seed)
        pars    <- prior_init_latent(I, K, a_psi, b_psi, a_vsi, b_vsi, omesq, a_sig, b_sig, a_tau, b_tau, kapsq)
        npars   <- length(unlist(pars))
        psi     <- pars$psi 
        vsisq   <- pars$vsisq 
        sigsq_u <- pars$sigsq_u 
        sigsq_v <- pars$sigsq_v
        tausq_u <- pars$tausq_u
        tausq_v <- pars$tausq_v
        nu      <- pars$nu
        bet     <- pars$bet
        gama    <- pars$gama
        xi      <- pars$xi
        eta     <- pars$eta
        zeta    <- pars$zeta
        U       <- pars$U
        V       <- pars$V
        
        ## mcmc main------------------------------------------------------------
        ptm <- proc.time()
        set.seed(seed)
        info <- mcmc_ic_latent_cpp(Ymcmc, K, a_psi, b_psi, a_vsi, b_vsi, omesq, a_sig, b_sig, a_tau, b_tau, kapsq, 
                                   psi, vsisq, sigsq_u, sigsq_v, tausq_u, tausq_v, nu, bet, gama, xi, eta, zeta, U, V, 
                                   nburn, nsams, nskip, ndisp, path_outs)
        ptm <- proc.time() - ptm
        eth <- round(as.numeric(ptm[3])/60/60, 3)
        
        ## save ----------------------------------------------------------------
        mcmc_info <- data.frame(I, K, niter = nburn + nskip*nsams, nburn, nsams, nskip, ndisp, ET_hrs = eth, nhyps, npars)
        write.csv(x = mcmc_info, file = paste0(path_outs, "mcmc_info.csv"), row.names = FALSE)
        write.csv(x = hyps, file = paste0(path_outs, "hyperpars.csv"), row.names = FALSE)
        info
}

mcmc_latent <- function(loc, case, model, dataset, K, nburn, nsams, nskip, ndisp, seed)
{
        ## paths ---------------------------------------------------------------
        path_code <- paste0(loc, "code/")
        path_data <- paste0(loc, "data/")
        path_outs <- paste0(loc, "outs/", dataset, "_", case, "_", model, "/")
        if (!dir.exists(path_outs)) dir.create(path_outs)
        path_outs <- paste0(loc, "outs/", dataset, "_", case, "_", model, "/", "K_", K, "/")
        if (!dir.exists(path_outs)) dir.create(path_outs)
        if (length(dir(path_outs)) > 0) unlink(paste0(path_outs, "*"))
        
        ## data ----------------------------------------------------------------
        load(paste0(path_data, "data_", dataset, ".RData"))
        
        ## libs and source -----------------------------------------------------
        suppressPackageStartupMessages(require(Rcpp))
        sourceCpp(paste0(path_code, "mcmc_", model, ".cpp"))
        
        ## hyper-parameter elicitation -----------------------------------------
        a_psi <- 1
        b_psi <- 1
        a_vsi <- 2
        b_vsi <- (a_vsi-1)/4
        omesq <- 0.25
        a_sig <- 2
        b_sig <- (a_sig-1)/(2*sqrt(2*K))
        a_tau <- 2
        b_tau <- (a_tau-1)/(sqrt(2*K))
        kapsq <- 1/(2*sqrt(2*K))
        hyps  <- data.frame(a_psi, b_psi, a_vsi, b_vsi, omesq, a_sig, b_sig, a_tau, b_tau, kapsq)
        nhyps <- length(hyps)
        
        ## prior initialization ------------------------------------------------
        set.seed(seed)
        pars    <- prior_init_latent(I, K, a_psi, b_psi, a_vsi, b_vsi, omesq, a_sig, b_sig, a_tau, b_tau, kapsq)
        npars   <- length(unlist(pars))
        psi     <- pars$psi 
        vsisq   <- pars$vsisq 
        sigsq_u <- pars$sigsq_u 
        sigsq_v <- pars$sigsq_v
        tausq_u <- pars$tausq_u
        tausq_v <- pars$tausq_v
        nu      <- pars$nu
        bet     <- pars$bet
        gama    <- pars$gama
        xi      <- pars$xi
        eta     <- pars$eta
        zeta    <- pars$zeta
        U       <- pars$U
        V       <- pars$V
        
        ## mcmc main------------------------------------------------------------
        ptm <- proc.time()
        set.seed(seed)
        mcmc_latent_cpp(Ymcmc, K, a_psi, b_psi, a_vsi, b_vsi, omesq, a_sig, b_sig, a_tau, b_tau, kapsq, 
                        psi, vsisq, sigsq_u, sigsq_v, tausq_u, tausq_v, nu, bet, gama, xi, eta, zeta, U, V, 
                        nburn, nsams, nskip, ndisp, path_outs)
        ptm <- proc.time() - ptm
        eth <- round(as.numeric(ptm[3])/60/60, 3)
        
        ## save ----------------------------------------------------------------
        mcmc_info <- data.frame(I, K, niter = nburn +nskip*nsams, nburn, nsams, nskip, ndisp, ET_hrs = eth, nhyps, npars)
        write.csv(x = mcmc_info, file = paste0(path_outs, "mcmc_info.csv"), row.names = FALSE)
        write.csv(x = hyps, file = paste0(path_outs, "hyperpars.csv"), row.names = FALSE)
}

mcmc_ll_latent <- function(loc, case, model, dataset, r, K, nburn, nsams, nskip, ndisp, seed)
{
        ## paths ---------------------------------------------------------------
        path_code <- paste0(loc, "code/")
        path_data <- paste0(loc, "data/")
        path_outs <- paste0(loc, "outs/", dataset, "_", case, "_", model, "/")
        if (!dir.exists(path_outs)) dir.create(path_outs)
        path_outs <- paste0(loc, "outs/", case, "_", dataset, "_", model, "/", "rep_", r, "/")
        if (!dir.exists(path_outs)) dir.create(path_outs)
        if (length(dir(path_outs)) > 0) unlink(paste0(path_outs, "*"))
        
        ## data ----------------------------------------------------------------
        load(paste0(path_data, "data_", dataset, ".RData"))
        
        ## libs and source -----------------------------------------------------
        suppressPackageStartupMessages(require(Rcpp))
        sourceCpp(paste0(path_code, "mcmc_ll_", model, ".cpp"))
        
        ## hyper-parameter elicitation -----------------------------------------
        a_psi <- 1
        b_psi <- 1
        a_vsi <- 2
        b_vsi <- (a_vsi-1)/4
        omesq <- 0.25
        a_sig <- 2
        b_sig <- (a_sig-1)/(2*sqrt(2*K))
        a_tau <- 2
        b_tau <- (a_tau-1)/(sqrt(2*K))
        kapsq <- 1/(2*sqrt(2*K))
        hyps  <- data.frame(a_psi, b_psi, a_vsi, b_vsi, omesq, a_sig, b_sig, a_tau, b_tau, kapsq)
        nhyps <- length(hyps)
        
        ## prior initialization ------------------------------------------------
        set.seed(seed)
        pars    <- prior_init_latent(I, K, a_psi, b_psi, a_vsi, b_vsi, omesq, a_sig, b_sig, a_tau, b_tau, kapsq)
        npars   <- length(unlist(pars))
        psi     <- pars$psi 
        vsisq   <- pars$vsisq 
        sigsq_u <- pars$sigsq_u 
        sigsq_v <- pars$sigsq_v
        tausq_u <- pars$tausq_u
        tausq_v <- pars$tausq_v
        nu      <- pars$nu
        bet     <- pars$bet
        gama    <- pars$gama
        xi      <- pars$xi
        eta     <- pars$eta
        zeta    <- pars$zeta
        U       <- pars$U
        V       <- pars$V
        
        ## mcmc main------------------------------------------------------------
        ptm <- proc.time()
        set.seed(seed)
        mcmc_ll_latent_cpp(Ymcmc, K, a_psi, b_psi, a_vsi, b_vsi, omesq, a_sig, b_sig, a_tau, b_tau, kapsq, 
                           psi, vsisq, sigsq_u, sigsq_v, tausq_u, tausq_v, nu, bet, gama, xi, eta, zeta, U, V, 
                           nburn, nsams, nskip, ndisp, path_outs)
        ptm <- proc.time() - ptm
        eth <- round(as.numeric(ptm[3])/60/60, 3)
        
        ## save ----------------------------------------------------------------
        mcmc_info <- data.frame(I, K, niter = nburn + nskip*nsams, nburn, nsams, nskip, ndisp, ET_hrs = eth, nhyps, npars)
        write.csv(x = mcmc_info, file = paste0(path_outs, "mcmc_info.csv"), row.names = FALSE)
        write.csv(x = hyps, file = paste0(path_outs, "hyperpars.csv"), row.names = FALSE)
}

mcmc_ll_latent_no_init <- function(loc, case, model, dataset, K, psi, vsisq, sigsq_u, sigsq_v, tausq_u, tausq_v, nu, bet, gama, xi, eta, zeta, U, V, nburn, nsams, nskip, ndisp, seed)
{
        ## paths ---------------------------------------------------------------
        path_code <- paste0(loc, "code/")
        path_data <- paste0(loc, "data/")
        path_outs <- paste0(loc, "outs/", dataset, "_", case, "_", model, "/")
        if (!dir.exists(path_outs)) dir.create(path_outs)
        path_outs <- paste0(loc, "outs/", dataset, "_", case, "_", model, "/", "K_", K, "/")
        if (!dir.exists(path_outs)) dir.create(path_outs)
        if (length(dir(path_outs)) > 0) unlink(paste0(path_outs, "*"))
        
        ## data ----------------------------------------------------------------
        load(paste0(path_data, "data_", dataset, ".RData"))
        
        ## libs and source -----------------------------------------------------
        suppressPackageStartupMessages(require(Rcpp))
        sourceCpp(paste0(path_code, "mcmc_", model, ".cpp"))
        
        ## hyper-parameter elicitation -----------------------------------------
        a_psi <- 1
        b_psi <- 1
        a_vsi <- 2
        b_vsi <- (a_vsi-1)/4
        omesq <- 0.25
        a_sig <- 2
        b_sig <- (a_sig-1)/(2*sqrt(2*K))
        a_tau <- 2
        b_tau <- (a_tau-1)/(sqrt(2*K))
        kapsq <- 1/(2*sqrt(2*K))
        hyps  <- data.frame(a_psi, b_psi, a_vsi, b_vsi, omesq, a_sig, b_sig, a_tau, b_tau, kapsq)
        nhyps <- length(hyps)
        
        ## prior initialization ------------------------------------------------
        set.seed(seed)
        pars    <- prior_init_latent(I, K, a_psi, b_psi, a_vsi, b_vsi, omesq, a_sig, b_sig, a_tau, b_tau, kapsq)
        npars   <- length(unlist(pars))
        
        ## mcmc main------------------------------------------------------------
        ptm <- proc.time()
        set.seed(seed)
        mcmc_latent_cpp(Ymcmc, K, a_psi, b_psi, a_vsi, b_vsi, omesq, a_sig, b_sig, a_tau, b_tau, kapsq, 
                        psi, vsisq, sigsq_u, sigsq_v, tausq_u, tausq_v, nu, bet, gama, xi, eta, zeta, U, V, 
                        nburn, nsams, nskip, ndisp, path_outs)
        ptm <- proc.time() - ptm
        eth <- round(as.numeric(ptm[3])/60/60, 3)
        
        ## save ----------------------------------------------------------------
        mcmc_info <- data.frame(I, K, niter = nburn +nskip*nsams, nburn, nsams, nskip, ndisp, ET_hrs = eth, nhyps, npars)
        write.csv(x = mcmc_info, file = paste0(path_outs, "mcmc_info.csv"), row.names = FALSE)
        write.csv(x = hyps, file = paste0(path_outs, "hyperpars.csv"), row.names = FALSE)
}