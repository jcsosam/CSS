list.of.packages <- c("Rcpp","RcppArmadillo","doParallel","foreach","igraph","AUC","data.table")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
    install.packages(new.packages)
}
lapply(list.of.packages,require, character.only=TRUE)
