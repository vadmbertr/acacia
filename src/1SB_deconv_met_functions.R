library(tictoc)
library(EpiDISH)
library(MeDeCom)
library(EMeth)
run_sup_deconvolution = function(met, dat, ref_profiles, hvg=NULL) {
  tic(met)
  if (met=="rlr") {
    if (!is.null(hvg)) {
      beta.m = dat[hvg,]
      ref.m = as.matrix(ref_profiles[hvg,])
    }
    else {
      beta.m = dat
      ref.m = as.matrix(ref_profiles)
    }
    if (det(t(beta.m)%*%beta.m)==0) {return(list(res=NA,time_elapsed=NA))}
    else {res <- t(epidish(beta.m, ref.m, method = "RPC")$estF)}
  }
  else if (met=="CIBERSORT") {
    if (!is.null(hvg)) {
      beta.m = dat[hvg,]
      ref.m = as.matrix(ref_profiles[hvg,])
    }
    else {
      beta.m = dat
      ref.m = as.matrix(ref_profiles)
    }
    res <- t(epidish(beta.m, ref.m, method = "CBS")$estF)
  }
  else if (met=="houseman") {
    res <- t(epidish(beta.m = dat, ref.m = as.matrix(ref_profiles), method = "CP")$estF)
    res[res<0]<-0
  }
  else if (met=="EMeth") {
    if (!is.null(hvg)) {
      Y = dat[hvg,]
      mu = as.matrix(ref_profiles[hvg,])
    }
    else {
      Y = dat
      mu = as.matrix(ref_profiles)
    }
    res <- t(cv.emeth(Y = Y, eta = rep(0,ncol(dat)), mu = mu, nu=nrow(dat)*(10^seq(-2,1,1)), aber=F)$result)
  }
  time_elapsed=toc()
  return(list(res=res,time_elapsed=time_elapsed))
}

run_unsup_deconvolution = function(met, dat, ref_profiles, option=c("Tmat","Amat"), prop_simu=NULL, hvg=NULL, dist=F) {
  k = ncol(ref_profiles)
  x = dat[rowSums(dat)>0,]
  if (is.null(hvg)) {hvg=seq(nrow(dat))}
  tic(met)
  if (met=="MeDeCom") {
    res <- runMeDeCom(dat[hvg,], Ks=k, lambdas=c(0,10^(-4:-1)), NINIT=10, NFOLDS=10, ITERMAX=300, NCORES=6)
    lambda <- res@parameters$lambdas[which.min(res@outputs$`1`$cve)]
    A <- getProportions(res, K=k, lambda=lambda)
    T_ <- getLMCs(res, K=k, lambda=lambda)
    res <- list(A_matrix=A,
                T_matrix=T_)
  }
  else if (met=="ICA") {
    if (det(t(dat)%*%dat)==0) {return(list(res=NA,time_elapsed=NA))}
    else {
      res = fastICA::fastICA(X = dat, n.comp = k, maxit = 1000, tol = 1e-09)
      res$names = row.names(dat)
      weighted.list <- deconica::generate_markers(df = res,
                                                n = min(30,nrow(dat)),
                                                return = "gene.ranked")
      ICA_scores_weighted = deconica::get_scores(res$X, weighted.list, summary = "weighted.mean", na.rm = TRUE)
      tmp_dat = t(ICA_scores_weighted)
      colnames(tmp_dat) = colnames(dat)
      res = list(A_matrix=abs(tmp_dat) %*% diag(1/colSums(abs(tmp_dat))),
               T_matrix=res$S)
    }
  }
  else if (met=="RefFreeEWAS") {
    res = RefFreeEWAS::RefFreeCellMix(Y = dat, K = k, verbose = F)
    res$Omega = t(apply(res$Omega, 2, function(x)
      ifelse(x < 0, 0, x))) #explicit NN constraints
    res$Omega = apply(res$Omega, 2, function(x)
      x / sum(x)) #explicit STO constraint
    res = list(T_matrix=res$Mu, A_matrix=res$Omega)
  }
  else if (met=="EDec") {
    res = t(EDec::run_edec_stage_1(meth_bulk_samples = as.matrix(dat), 
                                   informative_loci = rownames(dat)[hvg], 
                                   num_cell_types = ncol(ref_profiles))$proportions)
    res = apply(res, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraints
    res=list(A_matrix=res,T_matrix=NULL)
  }
  time_elapsed=toc()
  #time_elapsed=time_elapsed$toc-time_elapsed$tic
  
  if (option=="Tmat") {
    hvg2=TOAST::findRefinx(ref_profiles, nmarker = 1e3)
    elt1=ref_profiles[hvg2,]
    elt2=res$T_matrix[hvg2,]
  }
  else if (option=="Amat") {
    elt1=t(prop_simu)
    elt2=t(res$A_matrix)
  }
  if (dist) {
    row_order <- c(clue::solve_LSAP(dynutils::calculate_distance(t(elt1),t(elt2)), maximum = F))
  }
  else {
    row_order <- c(clue::solve_LSAP((cor(elt1,elt2)+1)^2, maximum = T))
  }
  res$A_matrix <- res$A_matrix[row_order,]
  rownames(res$A_matrix) <- colnames(elt1)
  return(list(res=res$A_matrix,time_elapsed=time_elapsed))
}
