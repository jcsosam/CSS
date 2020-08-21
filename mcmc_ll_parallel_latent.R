##------------------------------------------------------------------------------
## LATENT SPACE MODEL
## Chains parallel computation for a given K
##------------------------------------------------------------------------------

## path ------------------------------------------------------------------------
path_outs <- paste0(loc, "outs/", dataset, "_", case, "_", model, "/")

## libs ------------------------------------------------------------------------
suppressPackageStartupMessages(require(doParallel))

## info ------------------------------------------------------------------------
cat("///////////////////////////////", "\n",
    "           MCMC info           ", "\n",
    "Model: ", model,                  "\n",
    "Dataset: ", dataset,              "\n",
    "K = ", K,                         "\n",
    "Repetitions = ", nrep,            "\n",
    "nburn = ", nburn,                 "\n",
    "nskip = ", nskip,                 "\n",
    "nsams = ", nsams,                 "\n",
    "niter = ", nburn + nskip*nsams,   "\n", 
    "///////////////////////////////", "\n", sep = "")

## set cores and cluster -------------------------------------------------------
cl <- makeCluster(min(detectCores(), nrep))
registerDoParallel(cl)  
info_table <- foreach (r = 1:nrep, .combine = rbind, .inorder = FALSE, .errorhandling = "remove") %dopar% {
        source(paste0(loc, "code/", "mcmc_latent_models.R"))
        mcmc_ll_latent(loc, case, model, dataset, r, K, nburn, nsams, nskip, ndisp, r*seed)
}
stopCluster(cl)

## info ------------------------------------------------------------------------
cat("Parallel computation COMPLETED",  "\n", 
    "///////////////////////////////", "\n", sep = "")

##------------------------------------------------------------------------------