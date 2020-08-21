##------------------------------------------------------------------------------
## LATENT SPACE MODEL
## Chains parallel computation from Kmin to Kmax
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
    "K = ", Kmin, ",...,", Kmax,       "\n",
    "nburn = ", nburn,                 "\n",
    "nskip = ", nskip,                 "\n",
    "nsams = ", nsams,                 "\n",
    "niter = ", nburn + nskip*nsams,   "\n", 
    "///////////////////////////////", "\n", sep = "")

## set cores and cluster -------------------------------------------------------
cl <- makeCluster(min(detectCores(), Kmax-Kmin+1))
registerDoParallel(cl)  
info_table <- foreach (K = Kmin:Kmax, .combine = rbind, .inorder = TRUE, .errorhandling = "remove") %dopar% {
        source(paste0(loc, "code/", "mcmc_latent_models.R"))
        mcmc_ic_latent(loc, case, model, dataset, K, nburn, nsams, nskip, ndisp, seed)
}
stopCluster(cl)

## info ------------------------------------------------------------------------
cat("Parallel computation COMPLETED",  "\n", 
    "///////////////////////////////", "\n", sep = "")

## table -----------------------------------------------------------------------
colnames(info_table) <- c("AIC", "BIC", "DIC1", "DIC2", "WAIC1", "WAIC2")
info_table <- data.frame(K = Kmin:Kmax, info_table)
write.csv(x = info_table, file = paste0(path_outs, "ic.csv"), row.names = FALSE)

##------------------------------------------------------------------------------