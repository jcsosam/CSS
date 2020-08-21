vec2mat <- function (n, v) 
{
        p <- length(v)/n
        matrix(v, n, p, byrow = F)
}

mat.new.order <- function(mat, labs.old, labs.new)
{
        # reorders rows and cols of a given matrix according to labs.new  
        
        n <- nrow(mat)
        colnames(mat) <- rownames(mat) <- labs.old
        out <- matrix(NA, n, n)
        for(i in 1:n) {
                for(j in 1:n) {
                        out[i, j] <- mat[labs.new[i], labs.new[j]]
                }
        }
        colnames(out) <- rownames(out) <- labs.new
        out
}

procus.mat <- function (Z, Z0)
{
        # calculates the trasformation matrix 
        # Procrustes transform 
        # gives rotation, reflection of Z closest to Z0
        # Z0 taget configuration
        # Borg (2005), p. 430
        dec <- svd(t(Z0) %*% Z)
        (dec$v) %*% t(dec$u)
}

inter.matrix.latent <- function(loc, dataset, K, Xmcmc, bet.data, U.data, V.data)
{
        ## paths 
        path_code <- paste0(loc, "code/")
        path_data <- paste0(loc, "data/")
        
        ## depends on cmisc.cpp
        suppressPackageStartupMessages(require(Rcpp))
        sourceCpp(paste0(path_code, "cmisc.cpp"))
        
        ## data 
        load(paste0(path_data, "data_", dataset, ".RData"))
        
        ## interaction probility matrices for every net
        S <- dim(bet.data)[1]  
        I <- dim(bet.data)[2]
        P <- dim(Xmcmc)[2]
        
        TEMP <- array(NA, c(S, I, I, I))
        for (s in 1:S) {
                bet <- matrix(c(bet.data[s, ]), I, P)
                U   <- matrix(c(U.data[s, ]), I, K*I)
                V   <- matrix(c(V.data[s, ]), I, K*I)
                a   <- a_iter(I, K, Xmcmc, bet, U, V)
                p.lower.tail <- pnorm(a)
                for (j in 1:I) {
                        TEMP[s, , , j] <- matrix(c(p.lower.tail[ , j]), I, I)
                        diag(TEMP[s, , , j]) <- 0
                }
                if (s%%floor(0.1*S) == 0) cat(round(100*s/S, 1), "% completed", "\n", sep = "")
        }
        #THETA_{i,i',j} = mean(Phi( xi,i'^T beta_j + u_ij^T v_i'j ))
        THETA.hat <- array(0, c(I, I, I))
        for (j in 1:I) THETA.hat[,,j] <- apply(TEMP[ , , , j], c(2, 3), mean)
        THETA.hat
}

inter.matrix.swartz <- function(I, mu.data, al.data, be.data, ga.data, albe.data, alga.data)
{
        ## calculates the interaction probility matrices for every net
        S <- dim(mu.data)[1]
        
        mu.array   <- mu.data
        al.array   <- al.data
        be.array   <- be.data
        ga.array   <- ga.data
        albe.array <- array(c(albe.data), c(S, I, I))
        alga.array <- array(c(alga.data), c(S, I, I))
        
        TEMP <- array(NA, c(S, I, I, I))
        for (s in 1:S) {
                mu   <- mu.array[s]
                al   <- al.array[s, ]
                be   <- be.array[s, ]
                ga   <- ga.array[s, ]
                albe <- albe.array[s, , ]
                alga <- alga.array[s, , ]
                for (j in 1:I) {
                        for (i in 1:I) {
                                for (ii in (1:I)[-i]) {
                                        TEMP[s, i, ii, j] <- pnorm(mu + al[i] + be[ii] + ga[j] + albe[i, ii] + alga[i, j] + alga[ii, j])
                                }
                        }
                        diag(TEMP[s, , , j]) <- 0
                }
                
                if (s%%floor(0.1*S) == 0) cat(round(100*s/S, 1), "% completed", "\n", sep = "")
        }
        #THETA_{i,i',j} = mean(Phi( a_{i,i',j} ))
        THETA.hat <- array(0, c(I, I, I))
        for (j in 1:I) THETA.hat[ , , j] <- apply(TEMP[ , , , j], c(2, 3), mean)
        THETA.hat
}

ll.chain.rep.plot.latent <- function(loc, case, model, dataset, K, nrep, nsams, skip = 0) 
{
        source(paste0(loc, "code/mycolors.R"))
        
        ### load ###
        path_outs <- paste0(loc, "outs/", dataset, "_", case, "_", model, "/")
        LL <- matrix(NA, nsams - skip, nrep)
        for (r in 1:nrep) LL[ , r] <- as.matrix(read.table(paste0(path_outs, "rep_", r, "/","ll_chain.txt"), skip = skip))
        
        ### plot ###
        pdf(paste0(path_outs,"ll_chains_K_", K, ".pdf"), width = 10, height = 5, pointsize = 15)
        par(mfrow = c(1, 1), mar = c(4, 4, 3, 2) - 1, mgp = c(2, 1, 0), oma = 0.5*rep(1, 4))
        matplot(LL, type = "l", lty = 1, cex.axis = 0.6, col = mycolors[1:nrep], main = paste0(dataset, " data, ", model, " model, K = ", K),
                xlab = "Iteration", ylab = "Log-likelihood", ylim = range(LL) + c(0, 200))
        dev.off()
}

ll.chain.plot.latent <- function(loc, case, model, dataset, Kmin, Kmax, nsams, skip = 0) 
{
        source(paste0(loc, "code/mycolors.R"))
        
        ### load ###
        path_outs <- paste0(loc, "outs/", dataset, "_", case, "_", model, "/")
        LL <- matrix(NA, nsams - skip, Kmax - Kmin + 1)
        for (K in Kmin:Kmax) LL[ , K-Kmin+1] <- as.matrix(read.table(paste0(path_outs, "K_", K, "/","ll_chain.txt"), skip = skip))
        
        ### plot ###
        pdf(paste0(path_outs,"ll_chains.pdf"), width = 10, height = 5, pointsize = 15)
        par(mfrow = c(1, 1), mar = c(4, 4, 3, 2) - 1, mgp = c(2, 1, 0), oma = 0.5*rep(1, 4))
        matplot(LL, type = "l", lty = 1, cex.axis = 0.6, col = mycolors[Kmin:Kmax], main = paste0(dataset, " data, ", model, " model"),
                xlab = "Iteration", ylab = "Log-likelihood", ylim = range(LL) + c(0, 200))
        legend("top", legend = paste0("K = ", Kmin:Kmax), col = mycolors[Kmin:Kmax], fill = mycolors[Kmin:Kmax], 
               text.col = mycolors[Kmin:Kmax], border = mycolors[Kmin:Kmax], horiz = T,  cex = 0.7, bty = "n")
        dev.off()
}

criteria.plot <- function(loc, case, model, dataset, Kmin, Kmax, nsams)
{
        ### load ###
        path_outs  <- paste0(loc, "outs/", dataset, "_", case, "_", model, "/")
        info_table <- read.csv(paste0(path_outs, "ic.csv"))
        
        ### plot ###
        criteria.names <- c("AIC", "BIC", "DIC1", "DIC2", "WAIC1", "WAIC2")
        pdf(paste0(path_outs, "ic.pdf"), width = 15, height = 10, pointsize = 20)
        par(mfrow = c(2, 3), mar = c(4, 4, 3, 2) - 1, mgp = c(2, 1, 0), oma = 0.5*rep(1, 4))
        for (i in 1:length(criteria.names)) {
                plot(Kmin:Kmax, info_table[ , criteria.names[i]], type = "b", lwd = 2, cex = 0.8, col = "royalblue", 
                     xlab = "K", ylab = "Criterion(K)", main = criteria.names[i])
                abline(v = (Kmin:Kmax)[which.min(info_table[ , criteria.names[i]])], col = 2, lwd = 2, lty = 2)
                grid()
        }
        dev.off()
}

select_K <- function(loc, case, model, dataset, Kmin, Kmax, nsams)
{
        ### load ###
        path_outs  <- paste0(loc, "outs/", dataset, "_", case, "_", model, "/")
        info_table <- read.csv(paste0(path_outs, "ic.csv"))
        v <- apply(X = info_table[,-1], MARGIN = 2, FUN = function(x) info_table[which.min(x), 1])
        uv <- unique(v)
        min(uv[which.max(tabulate(match(v, uv)))])
}

mysample <- function (x, ...) {
        if (length(x) == 1) {
                out <- x[1]
        } else {
                out <- sample(x, ...)
        }
        out
}

get_folds <- function(loc, dataset, L)
{
        ## paths
        path_code <- paste0(loc, "code/")
        path_data <- paste0(loc, "data/")
        path_outs <- paste0(loc, "outs/", dataset, "_goodness/")
        
        ## data
        load(paste0(path_data, "data_", dataset, ".RData"))
        
        ## groups formation
        ones.idx <- Ymcmc == 1
        nones <- colSums(ones.idx)
        
        group.idx <- matrix(NA, N, I)
        for (j in 1:I) {
                for (l in 1:L) {
                        tmp <- is.na(group.idx[ , j])
                        tmp[(1:I)*(I+1) - I] <- F  # main diagonal indices
                        tmp[ones.idx[ , j] == F] <- F  # 0s indices
                        if (sum(tmp) > 0) {
                                idx <- mysample(x = which(tmp == T), size = floor(nones[j]/L), replace = F)
                                group.idx[idx, j] <- l
                                rm(tmp, idx)
                        }
                }
                # sample remaining ones for each fold
                r <- nones[j] - floor(nones[j]/L) * L
                r.idx <- rep(c(0, 1), c(L-r, r))[order(runif(L))]
                for (l in 1:L) {
                        if (r.idx[l] == 1) {
                                tmp <- is.na(group.idx[ , j])
                                tmp[(1:I)*(I+1) - I] <- F  # main diagonal indices
                                tmp[ones.idx[ , j] == F] <- F  # 0s indices
                                if (sum(tmp) > 0) {
                                        idx <- mysample(x = which(tmp == T), size = 1, replace = F)
                                        group.idx[idx, j] <- l
                                        rm(tmp, idx)
                                }
                        }
                } 
        }
        
        nones.group <- matrix(NA, L, I)
        for (j in 1:I) {
                for (l in 1:L) {
                        nones.group[l, j] <- sum( group.idx[,j] == l, na.rm = T )
                }
        }
        
        for (j in 1:I) {
                for (l in 1:L) {
                        tmp <- is.na(group.idx[ , j])
                        tmp[(1:I)*(I+1) - I] <- F  # main diagonal indices
                        tmp[ones.idx[ , j] == T] <- F  # 1s indices
                        idx <- sample(x = which(tmp == T), size = I*(I-1)/L - nones.group[l,j], replace = F)
                        group.idx[idx, j] <- l
                        rm(tmp, idx)
                }  
        }
        
        group.idx[is.na(group.idx)] <- 0
        
        group.idx
}

MCC <- function(pred, act) {
        TP <- sum(act == 1 & pred == 1)
        TN <- sum(act == 0 & pred == 0)
        FP <- sum(act == 0 & pred == 1)
        FN <- sum(act == 1 & pred == 0)
        if (any((TP + FP) == 0, (TP + FN) == 0, (TN + FP) == 0, (TN + FN) == 0)) { 
                denom <- 1
        } else {
                denom <- as.double(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
        } 
        ((TP * TN) - (FP * FN))/sqrt(denom)
}

post.summary1 <- function(X, what, r = 2)
{
        X <- as.matrix(X)
        p <- ncol(X)
        out <- NULL
        for (i in 1:p) {
                out <- rbind(out, round(c(mean(X[,i]), 
                                          sd(X[,i]), 
                                          median(X[,i]), 
                                          quantile(X[,i], c(.025,.975))), r))
        }
        rownames(out) <- colnames(X)
        colnames(out) <- c('mean', 'sd', '50%', '2.5%', '97.5%')
        as.data.frame(out)[what]
}