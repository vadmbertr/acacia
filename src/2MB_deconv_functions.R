library(tictoc)
library(MeDeCom)

tpm_norm <- function(dat) {
  TPM = function(counts,lengths) {
    #rownames(counts) = tolower(rownames(counts))
    #names(lengths) = tolower(names(lengths))
    A = intersect(rownames(counts),names(lengths))
    counts = counts[A,]
    lengths = lengths[A]
    rate = counts / lengths
    apply(rate,2,function(x) 1e6*x/sum(x))
  }
  load("../../src/TPM/human_lengths.rda")
  matrix = TPM(counts = as.matrix(dat), lengths = human_lengths)
  if (is.null(rownames(matrix))) {rownames(matrix) = paste("feature",seq(nrow(matrix)))}
  rownames(matrix) = toupper(rownames(matrix))
  return(matrix)
}

run_MtoR_deconvolutionA <- function(met, D, Amet, A=NULL) {
  k = nrow(Amet)
  tic()
  if (met=="ICA") {
    res_baseline <- fastICA::fastICA(X = D,
                            n.comp = k,
                            maxit = 1000,
                            tol = 1e-09, method = "C")
    res <- fastICA::fastICA(X = D,
                            n.comp = k,
                            maxit = 1000,
                            tol = 1e-09, method = "C",
                            w.init = pracma::pinv(Amet %*% res_baseline$K))
    res$names = rownames(D)
    rownames(res$X) = rownames(D)
    weighted.list <- deconica::generate_markers(df = res,
                                                n = 30,
                                                return = "gene.ranked")
    res$A_weighted = t(deconica::get_scores(res$X, weighted.list, summary = "weighted.mean", na.rm = TRUE))
    colnames(res$A_weighted) = colnames(D)
    res$A_weighted = abs(res$A_weighted) %*% diag(1/colSums(abs(res$A_weighted)))
    elt1=t(A)
    elt1=elt1[,sort(colnames(elt1))]
    elt2=t(res$A_weighted)
    row_order <- c(clue::solve_LSAP((cor(elt1,elt2)+1)^2, maximum = T))
    res <- res$A_weighted[row_order,]
  }
  else if (met=="NMF") {
    x = D[rowSums(D)>0,]
    T_fromM = x%*%pracma::pinv(Amet)
    res <- NMF::nmf(x = x, rank = k, method = "snmf/r", seed = NMF::nmfModel(H=Amet, W=T_fromM))
    res <- apply(X = res@fit@H, 2, function(z) {
      z/sum(z)})
    if (any(rowSums(res)==0)) {
      res[rowSums(res)==0,1] <- 1e-5
    }
    elt1=t(A)
    elt1=elt1[,sort(colnames(elt1))]
    elt2=t(res)
    row_order <- c(clue::solve_LSAP((cor(elt1,elt2)+1)^2, maximum = T))
    res <- res[row_order,]
  }
  time_elapsed=toc()
  return(list(res=res,time_elapsed=time_elapsed))
}

run_RtoM_deconvolutionA <- function(met, D, Arna, A=NULL) {
  k = nrow(Arna)
  tic()
  if (met=="ICA") {
    res_baseline <- fastICA::fastICA(X = D,
                            n.comp = k,
                            maxit = 1000,
                            tol = 1e-09, method = "C")
    res <- fastICA::fastICA(X = D,
                            n.comp = k,
                            maxit = 1000,
                            tol = 1e-09, method = "C",
                            w.init = pracma::pinv(Arna %*% res_baseline$K))
    res$names = rownames(D)
    rownames(res$X) = rownames(D)
    weighted.list <- deconica::generate_markers(df = res,
                                                n = 30,
                                                return = "gene.ranked")
    ICA_scores_weighted = t(deconica::get_scores(res$X, weighted.list, summary = "weighted.mean", na.rm = TRUE))
    colnames(ICA_scores_weighted) = colnames(D)
    A_matrix=abs(ICA_scores_weighted) %*% diag(1/colSums(abs(ICA_scores_weighted)))
    elt1=t(A)
    elt1=elt1[,sort(colnames(elt1))]
    elt2=t(A_matrix)
    row_order <- c(clue::solve_LSAP((cor(elt1,elt2)+1)^2, maximum = T))
    res <- A_matrix[row_order,]
  }
  else if (met=="RefFreeEWAS") {
    res = RefFreeEWAS::RefFreeCellMix(Y = D, K = k, mu0 = D%*%pracma::pinv(Arna), verbose = F)
    res$Omega = t(apply(res$Omega, 2, function(x)
      ifelse(x < 0, 0, x))) #explicit NN constraints
    res$Omega = apply(res$Omega, 2, function(x)
      x / sum(x)) #explicit STO constraint
    elt1=t(A)
    elt1=elt1[,sort(colnames(elt1))]
    elt2=t(res$Omega)
    row_order <- c(clue::solve_LSAP((cor(elt1,elt2)+1)^2, maximum = T))
    res <- res$Omega[row_order,]
  }
  else if (met=="MeDeCom") {
    library(parallel)
    hvg = TOAST::findRefinx(D, nmarker = 1e3)
    result <- MeDeCom::runMeDeCom(D[hvg,], Ks=k, lambdas=c(0,10^(-4:-1)), startA=t(Arna), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=6)
    lambda <- result@parameters$lambdas[which.min(result@outputs$`1`$cve)]
    Apred <- MeDeCom::getProportions(result, K=k, lambda=lambda)
    elt1=t(A)
    elt1=elt1[,sort(colnames(elt1))]
    elt2=t(Apred)
    row_order <- c(clue::solve_LSAP((cor(elt1,elt2)+1)^2, maximum = T))
    res <- Apred[row_order,]
  }
  time_elapsed=toc()
  return(list(res=res,time_elapsed=time_elapsed))
}

run_concat_deconvolution <- function(met, D, Tref=NULL, A=NULL) {
  if (!is.null(A)) {
    k = nrow(A)
  }
  tic()
  if (met=="rlr") {
    res <- t(EpiDISH::epidish(beta.m = D, ref.m = as.matrix(Tref), method = "RPC")$estF)
  }
  else if (met=="cibersort") {
    source("../../src/CIBERSORT/CIBERSORT.R")
    res = CIBERSORT(Tref, D) 
    res = t(res[,1:(ncol(res)-3)])
  }
  else if (met=="ica") {
    res <- fastICA::fastICA(X = D,
                            n.comp = k,
                            maxit = 1000,
                            tol = 1e-09, method = "C")
    res$names = rownames(D)
    rownames(res$X) = rownames(D)
    weighted.list <- deconica::generate_markers(df = res,
                                                n = 30,
                                                return = "gene.ranked")
    res$A_weighted = t(deconica::get_scores(res$X, weighted.list, summary = "weighted.mean", na.rm = TRUE))
    colnames(res$A_weighted) = colnames(D)
    res$A_weighted = abs(res$A_weighted) %*% diag(1/colSums(abs(res$A_weighted)))
    elt1=t(A)
    elt1=elt1[,sort(colnames(elt1))]
    elt2=t(res$A_weighted)
    row_order <- c(clue::solve_LSAP((cor(elt1,elt2)+1)^2, maximum = T))
    res <- res$A_weighted[row_order,]
  }
  time_elapsed=toc()
  return(list(res=res,time_elapsed=time_elapsed))
}

run_jdr_deconvolution_sup <- function(met, D_jdr, Tref=NULL) {
  tic()
  if (met=="rlr") {
    res <- t(EpiDISH::epidish(beta.m = D_jdr, ref.m = as.matrix(Tref), method = "RPC")$estF)
    }
  else if (met=="cibersort") {
    source("../../src/CIBERSORT/CIBERSORT.R")
    res = CIBERSORT(Tref, D_jdr) 
    res = t(res[,1:(ncol(res)-3)])
  }
  time_elapsed=toc()
  time_elapsed=time_elapsed$toc-time_elapsed$tic
  return(list(res=res,time_elapsed=time_elapsed))
}

run_jdr_deconvolution_unsup <- function(met, D_jdr, Atrue=NULL) {
  k=nrow(Atrue)
  x=D_jdr-min(D_jdr)
  tic()
  if (met=="ica") {
    res <- fastICA::fastICA(X = D_jdr,
                            n.comp = k,
                            maxit = 1000,
                            tol = 1e-09, method = "C")
    res$names = rownames(D_jdr)
    rownames(res$X) = rownames(D_jdr)
    weighted.list <- deconica::generate_markers(df = res,
                                                n = 10,
                                                return = "gene.ranked")
    res$A_weighted = t(deconica::get_scores(res$X, weighted.list, summary = "weighted.mean", na.rm = TRUE))
    colnames(res$A_weighted) = colnames(D_jdr)
    res = list(A_matrix=abs(res$A_weighted) %*% diag(1/colSums(abs(res$A_weighted))),
               T_matrix=res$S)
  }
  else if (met == "nmf") {
    res <- NMF::nmf(x = x, rank = k, method = "snmf/r", seed = 1)
    res <- list(A_matrix=apply(X = res@fit@H, 2, function(z) {
      z/sum(z)}), #explicit STO constraint
      T_matrix=res@fit@W)
    if (any(rowSums(res$A_matrix)==0)) {
      res$A_matrix[rowSums(res$A_matrix)==0,1] <- 1e-5
    }
  }
  else if (met=="RefFreeEWAS") {
    res = RefFreeEWAS::RefFreeCellMix(Y = D_jdr, K = k, verbose = F)
    res$Omega = t(apply(res$Omega, 2, function(x)
      ifelse(x < 0, 0, x))) #explicit NN constraints
    res$Omega = apply(res$Omega, 2, function(x)
      x / sum(x)) #explicit STO constraint
    res = list(T_matrix=res$Mu, A_matrix=res$Omega)
  }
  else if (met=="EDec") {
    res = t(EDec::run_edec_stage_1(meth_bulk_samples = as.matrix(D_jdr), 
                                   informative_loci = rownames(D_jdr), 
                                   num_cell_types = k)$proportions)
    res = apply(res, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraints
    res=list(A_matrix=res,T_matrix=NULL)
  }
  else if (met=="MeDeCom") {
    res <- MeDeCom::runMeDeCom(D_jdr, Ks=k, lambdas=c(0,10^(-4:-1)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=6)
    lambda <- res@parameters$lambdas[which.min(res@outputs$`1`$cve)]
    A <- getProportions(res, K=k, lambda=lambda)
    T_ <- getLMCs(res, K=k, lambda=lambda)
    res <- list(A_matrix=A,
                T_matrix=T_)
  }
  time_elapsed=toc()
  if (length(time_elapsed)==1) {time_elapsed=time_elapsed} else {time_elapsed=time_elapsed$toc-time_elapsed$tic}
  
  elt1=t(Atrue)
  elt2=t(res$A_matrix)
  elt1=elt1[,sort(colnames(elt1))]
  corr=cor(elt1,elt2)
  corr[is.na(corr)]=0
  row_order <- c(clue::solve_LSAP((corr+1)^2, maximum = T))
  res$A_matrix <- res$A_matrix[row_order,]
  rownames(res$A_matrix) <- colnames(elt1)
  return(list(res=res$A_matrix,time_elapsed=time_elapsed))
}

