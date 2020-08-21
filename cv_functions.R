cv_parallel_butts <- function (group.idx, L, Ymcmc, nburn, nsams, nskip, path_code)
{
        ## libs
        suppressPackageStartupMessages(require(doParallel))
        suppressPackageStartupMessages(require(foreach))
        
        ## CV computation (parallel)
        cl <- makeCluster(min(detectCores(), L))
        registerDoParallel(cl)
        auc_table <- foreach (l = 1:L, .packages = c("Rcpp","AUC"), .combine = "rbind", .inorder = F, .errorhandling = "remove") %dopar% {
                sourceCpp(paste0(path_code, "yppp_butts.cpp"))
                Yna <- Ymcmc
                Yna[group.idx == l] <- NA
                yppp <- c(yppp_butts(Yna, nburn, nsams, nskip))
                y.true <- as.factor(Ymcmc[group.idx == l])
                auc(roc(yppp, y.true))
        }
        stopCluster(cl)
        
        colnames(auc_table) <- c("auc")
        rownames(auc_table) <- c()
        as.data.frame(auc_table)
}

cv_parallel_swartz <- function (group.idx, L, Ymcmc, nburn, nsams, nskip, path_code)
{
        ## libs
        suppressPackageStartupMessages(require(doParallel))
        suppressPackageStartupMessages(require(foreach))
        
        ## CV computation (parallel)
        cl <- makeCluster(min(detectCores(), L))
        registerDoParallel(cl)
        auc_table <- foreach (l = 1:L, .packages = c("Rcpp","AUC"), .combine = "rbind", .inorder = F, .errorhandling = "remove") %dopar% {
                sourceCpp(paste0(path_code, "yppp_swartz.cpp"))
                Yna <- Ymcmc
                Yna[group.idx == l] <- NA
                yppp <- c(yppp_swartz(Yna, nburn, nsams, nskip))
                y.true <- as.factor(Ymcmc[group.idx == l])
                auc(roc(yppp, y.true))
        }
        stopCluster(cl)

        colnames(auc_table) <- c("auc")
        rownames(auc_table) <- c()
        as.data.frame(auc_table)
}

cv_parallel_latent <- function(group.idx, L, Ymcmc, Kmin, Kmax, nburn, nsams, nskip, path_code)
{
        ## libs
        suppressPackageStartupMessages(require(doParallel))
        suppressPackageStartupMessages(require(foreach))
        
        ## indices
        Kl <- as.matrix(expand.grid(Kmin:Kmax, 1:L))
        H <- nrow(Kl)
        
        ## CV computation (parallel)
        cl <- makeCluster(min(detectCores(), H))
        registerDoParallel(cl)
        A <- foreach (h = 1:H, .packages = c("Rcpp","AUC"), .combine = "rbind", .inorder = F) %dopar% {
                sourceCpp(paste0(path_code, "yppp_latent.cpp"))
                Yna <- Ymcmc
                Yna[group.idx == Kl[h,2]] <- NA
                yppp <- c(yppp_latent(Yna, K = Kl[h,1], nburn, nsams, nskip))
                y.true <- as.factor(Ymcmc[group.idx == Kl[h,2]])
                c(Kl[h,1], auc(roc(yppp, y.true)))
        }
        stopCluster(cl)
       
        auc_table <- NULL
        for(K in Kmin:Kmax) auc_table <- cbind(auc_table, A[A[,1] == K, 2])
        auc_table <- as.data.frame(auc_table)
        colnames(auc_table) <- paste0("auc K = ", Kmin:Kmax)
        rownames(auc_table) <- c()
        auc_table
}