library(tictoc)
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
  rownames(matrix) = toupper(rownames(matrix))
  return(matrix)
}

marker_selection <- function(ref_profiles) {
  marker_list <- lapply(colnames(ref_profiles), function(x) {
    expr_level_1 <- ref_profiles[,x]
    expr_level_2 <- apply(ref_profiles[,colnames(ref_profiles)!=x], 1, function(y) max(y))
    rel_expr_level_1 <- ref_profiles[,x]/rowSums(ref_profiles)
    rel_expr_level_1[is.na(rel_expr_level_1)] <- 0
    marker <- names(rel_expr_level_1)[rel_expr_level_1>0.99]
    return(marker)})
  names(marker_list) = colnames(ref_profiles)
}

run_sup_deconvolution = function(met, dat, ref_profiles) {
  tic(met)
  if (met=="OLSmag") {
    RESULTS = apply(dat, 2, function(x)
      lm(x ~ as.matrix(ref_profiles))$coefficients[-1])
    RESULTS = apply(RESULTS, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    res = apply(RESULTS, 2, function(x)
      x / sum(x)) #explicit STO constraint
    rownames(res) <-
      unlist(lapply(strsplit(rownames(res), ")"), function(x)
        x[2]))
  }
  else if (met=="OLSmagTPM") {
    RESULTS = apply(tpm_norm(dat), 2, function(x)
      lm(x ~ as.matrix(tpm_norm(ref_profiles)))$coefficients[-1])
    RESULTS = apply(RESULTS, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    res = apply(RESULTS, 2, function(x)
      x / sum(x)) #explicit STO constraint
    rownames(res) <-
      unlist(lapply(strsplit(rownames(res), "))"), function(x)
        x[2]))
  }
  else if (met=="ols") {
    res <- t(granulator::deconvolute(m = tpm_norm(dat), sigMatrix = tpm_norm(ref_profiles), methods = met)$proportions$ols_sig1) #TPM norm
    rownames(res) <-  gsub("\\.","_",rownames(res))
    res <- res[colnames(ref_profiles),colnames(dat)]
    res = apply(res, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    res = apply(res, 2, function(x)
      x / sum(x)) #explicit STO constraint
  }
  else if (met=="DCQ") {
    RESULTS = t(ComICS::dcq(reference_data = ref_profiles,
                            mix_data = dat,
                            marker_set = as.data.frame(row.names(as.matrix(ref_profiles))),
                            number_of_repeats = 10)$average)
    RESULTS = apply(RESULTS, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    RESULTS = apply(RESULTS, 2, function(x)
      x / sum(x)) #explicit STO constraint
    res = RESULTS
  }
  else if (met=="elastic_net") {
    RESULTS = apply(dat, 2, function(z)
      coef(glmnet::glmnet(x = ref_profiles,
                          y = z,
                          alpha = 0.2,
                          standardize = TRUE,
                          lambda = glmnet::cv.glmnet(as.matrix(ref_profiles), z)$lambda.1se))[1:ncol(ref_profiles) + 1, ])
    RESULTS = apply(RESULTS, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    RESULTS = apply(RESULTS, 2, function(x)
      x / sum(x)) #explicit STO constraint
    res = RESULTS
  }
  else if (met=="RLR"){ #RLR = robust linear regression
    RESULTS = do.call(cbind.data.frame,lapply(apply(dat,2,function(x) MASS::rlm(x ~ as.matrix(ref_profiles), maxit=100)), function(y) y$coefficients[-1]))
    RESULTS = apply(RESULTS,2,function(x) ifelse(x < 0, 0, x)) #explicit NN constraint
    RESULTS = apply(RESULTS,2,function(x) x/sum(x)) #explicit STO constraint
    rownames(RESULTS) <- unlist(lapply(strsplit(rownames(RESULTS),")"),function(x) x[2]))
    res = RESULTS
  }
  else if (met=="DeconRNASeq") { #NN quadratic programming; lsei function (default: type=1, meaning lsei from quadprog)
    require(pcaMethods)
    res = t(DeconRNASeq::DeconRNASeq(datasets = as.data.frame(dat), signatures = as.data.frame(ref_profiles), proportions = NULL, checksig = FALSE, known.prop = FALSE, use.scale = FALSE, fig = FALSE)$out.all)
    colnames(res) = colnames(dat)
  }
  else if (met=="nnls") {
    res <- t(granulator::deconvolute(m = tpm_norm(dat), sigMatrix = tpm_norm(ref_profiles), methods = met)$proportions$nnls_sig1) #TPM norm
    res <-  res[gsub("_",".",colnames(ref_profiles)),]
    rownames(res) <- colnames(ref_profiles)
    res <- apply(res, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    res <- apply(res, 2, function(x)
      x / sum(x)) #explicit STO constraint
    res <- res[,colnames(dat)]
  }
  else if (met=="svr4") {
    dat_tpm = tpm_norm(dat)
    hvg = TOAST::findRefinx(dat_tpm, nmarker = 1e4)
    res <- t(granulator::deconvolute(m = dat_tpm[hvg,], sigMatrix = tpm_norm(ref_profiles)[hvg,], methods = 'svr', use_cores = 4)$proportions$svr_sig1) #TPM norm
    res <-  res[gsub("_",".",colnames(ref_profiles)),]
    rownames(res) <- colnames(ref_profiles)
    res <- res[,colnames(dat)]
    res = apply(res, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    res = apply(res, 2, function(x)
      x / sum(x)) #explicit STO constraint
  }
  else if (met=="svr3") {
    dat_tpm = tpm_norm(dat)
    hvg = TOAST::findRefinx(dat_tpm, nmarker = 1e3)
    res <- t(granulator::deconvolute(m = dat_tpm[hvg,], sigMatrix = tpm_norm(ref_profiles)[hvg,], methods = 'svr', use_cores = 4)$proportions$svr_sig1) #TPM norm
    res <-  res[gsub("_",".",colnames(ref_profiles)),]
    rownames(res) <- colnames(ref_profiles)
    res <- res[,colnames(dat)]
    res = apply(res, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    res = apply(res, 2, function(x)
      x / sum(x)) #explicit STO constraint
  }
  else if (met=="svr") {
    dat_tpm = tpm_norm(dat)
    res <- t(granulator::deconvolute(m = dat_tpm, sigMatrix = tpm_norm(ref_profiles), methods = 'svr', use_cores = 4)$proportions$svr_sig1) #TPM norm
    res <-  res[gsub("_",".",colnames(ref_profiles)),]
    rownames(res) <- colnames(ref_profiles)
    res <- res[,colnames(dat)]
    res = apply(res, 2, function(x)
      ifelse(x < 0, 0, x)) #explicit NN constraint
    res = apply(res, 2, function(x)
      x / sum(x)) #explicit STO constraint
  }
  else if (met=="nmf/s1") {
    
  }
  else if (met=="nmf/s2") {
    
  }
  else if (met=="CIBERSORT4") {
    source("../../src/CIBERSORT/CIBERSORT.R")
    dat_tpm = tpm_norm(dat)
    hvg = TOAST::findRefinx(dat_tpm, nmarker = 1e4)
    res = CIBERSORT(tpm_norm(ref_profiles)[hvg,], dat_tpm[hvg,]) 
    res = t(res[,1:(ncol(res)-3)])
  }
  else if (met=="CIBERSORT3") {
    source("../../src/CIBERSORT/CIBERSORT.R")
    dat_tpm = tpm_norm(dat)
    hvg = TOAST::findRefinx(dat_tpm, nmarker = 1e3)
    res = CIBERSORT(tpm_norm(ref_profiles)[hvg,], dat_tpm[hvg,]) 
    res = t(res[,1:(ncol(res)-3)])
  }
  else if (met=="CIBERSORT") {
    source("../../src/CIBERSORT/CIBERSORT.R")
    dat_tpm = tpm_norm(dat)
    res = CIBERSORT(tpm_norm(ref_profiles), dat_tpm) 
    res = t(res[,1:(ncol(res)-3)])
  }
  else if (met=="WISP") {
    getWeight = function(data, centro, scaling = c("none", "scale", "center")[1], cutoff_gobalFtest = 0.05, Rsquared_cutoff = 0.5, cutoff_ttest_weights = 0.05, sum_LessThanOne = TRUE) {
      g = intersect(rownames(data), rownames(centro))
      data = data[g, ]
      centro = as.matrix(centro[g, ])
      if (scaling == "scale"){
        data = scale(data, scale = TRUE)
        centro = scale(centro, scale = TRUE)
      }
      if (scaling == "center"){
        data = scale(data, scale = FALSE)
        centro = scale(centro, scale = FALSE)
      }
      nJ = ncol(centro)
      A = centro
      nSubj = ncol(data)
      mixCoef = matrix(0, nSubj, nJ)
      rownames(mixCoef) = colnames(data)
      colnames(mixCoef) = colnames(centro)
      Amat = cbind(rep(-1, nJ), diag(nJ))
      b0vec = c(-1, rep(0, nJ))
      if(sum_LessThanOne==TRUE){
        meq = 0
      } else {
        meq = 1
      }
      output = data.frame(t(apply(data,2, function(y) {
        obs = which(!is.na(y))
        Dmat = t(centro[obs, ]) %*% centro[obs, ]
        diag(Dmat) <- diag(Dmat) + 1e-08
        sc <- norm(Dmat,"2")
        mixCoef = quadprog::solve.QP(Dmat/sc, (t(centro[obs, ]) %*%
                                                 y[obs])/sc, Amat, b0vec, meq = meq)$sol
        B = as.matrix(y)
        coeff = round(mixCoef, 4)
        names(coeff) = paste("weight", colnames(centro), sep = ".")
        vBeta = matrix(coeff)
        dSigmaSq = sum((B - A %*% vBeta)^2)/(nrow(A) - ncol(A))
        dTotalSq = sum((B)^2)/(nrow(A))
        dModelSq = sum((A %*% vBeta)^2)/(ncol(A))
        mVarCovar = try(dSigmaSq * chol2inv(chol(t(A) %*% A)))
        Adjusted.R.squared = round((dTotalSq - dSigmaSq)/dTotalSq,3)
        Ftest = dModelSq/dSigmaSq
        p.Ftest = stats::pf(q = Ftest, df1 = ncol(A), df2 = nrow(A) - ncol(A), lower.tail = FALSE)
        ng = nrow(centro)
        if (!is.character(mVarCovar)){
          vStdErr = sqrt(diag(mVarCovar))
          CI.inf = sapply(1:nJ,function(j){round(coeff[j] - (stats::qt(1-cutoff_ttest_weights/2,ng-nJ)*vStdErr[j]),2)})
          CI.sup = sapply(1:nJ,function(j){round(coeff[j] + (stats::qt(1-cutoff_ttest_weights/2,ng-nJ)*vStdErr[j]),2)})
          tvalue = sapply(1:nJ,function(j){ coeff[j]/vStdErr[j]})
          p.Ttest = sapply(1:nJ,function(j){stats::pt(abs(coeff[j]/vStdErr[j]), df=ng-nJ, lower=FALSE)*2})
        } else { mVarCovar = NA
        vStdErr = CI.inf = CI.sup = p.Ttest = rep(NA,nJ)
        }
        names(CI.inf) = paste("CI.inf.", colnames(centro),sep="")
        names(CI.sup) = paste("CI.sup.", colnames(centro),sep="")
        names(p.Ttest) = paste("Pvalue.", colnames(centro),sep="")
        coeff.filtered = coeff
        coeff.filtered[p.Ttest>cutoff_ttest_weights]=0
        names(coeff.filtered) = paste(names(coeff),".filtered",sep="")
        dist.Obs.Model = round(sqrt(sum(((centro %*% coeff) -y)^2)),2)
        c(coeff, dist.Obs.Model = dist.Obs.Model,Ftest.pvalue = p.Ftest, Adjusted.R.squared = Adjusted.R.squared,CI.inf,CI.sup,p.Ttest,coeff.filtered)})))
      output$topWeightedClass = colnames(centro)[apply(output[,1:nJ], 1, which.max)]
      output$deltaTopWeights = apply(output[, 1:nJ], 1, function(z) abs(diff(sort(z,decreasing = TRUE)[1:2])))
      CI = sapply(1:nJ,function(j) {
        CIinf = sapply(output[,gsub(" ",".",paste("CI.inf.", colnames(centro)[j], sep=""))], function(x) max(0,x))
        CIsup = sapply(output[,gsub(" ",".",paste("CI.sup.", colnames(centro)[j], sep=""))], function(x) min(1,x))
        paste("[",CIinf , ", ", CIsup, "]", sep="")})
      colnames(CI) = colnames(centro)
      output$CI = CI
      output$WARNING = as.factor(c("OK", "LIMIT")[(output$Ftest.pvalue>cutoff_gobalFtest | output$Adjusted.R.squared<Rsquared_cutoff)+1])
      output = output[,-c(grep("CI.inf", colnames(output)), grep("CI.sup", colnames(output)))]
      return(output)
    }
    resW <- getWeight(dat, ref_profiles, scaling = "none")
    res <- t(resW[,grep("filtered",colnames(resW))])
    res <- apply(res, 2, function(x) x/sum(x))
    rownames(res) = colnames(ref_profiles)
  }
  time_elapsed=toc()
  return(list(res=res,time_elapsed=time_elapsed))
}

run_unsup_deconvolution = function(met, dat, ref_profiles, option=c("Tmat","Amat"), dist=F, prop_simu=NULL) {
  k = ncol(ref_profiles)
  x = dat[rowSums(dat)>0,]
  tic(met)
  if (met == "ICA") {
    res <- gedepir::run_deconv(dat,k=k,method="ICA")
  }
  else if (met == "DECONica") {
    res <- fastICA::fastICA(X = dat, n.comp = k, maxit = 1000, tol = 1e-09)
    res$names = row.names(dat)
    weighted.list <- deconica::generate_markers(df = res,
                                                n = 30,
                                                return = "gene.ranked")
    res$A_weighted = t(deconica::get_scores(res$X, weighted.list, summary = "weighted.mean", na.rm = TRUE))
    colnames(res$A_weighted) = colnames(dat)
    #deconica::stacked_proportions_plot(res$A_weighted)
    res$A_weighted = abs(res$A_weighted) %*% diag(1/colSums(abs(res$A_weighted)))
    res=list(A_matrix=res$A_weighted,
             T_matrix=res$S)
  }
  else if (met == "NMF") {
    res <- gedepir::run_deconv(x,k=k,method="NMF")
    res$A_matrix <- apply(res$A_matrix, 2, function(z) {
      z/sum(z)}) #explicit STO constraint
    if (any(rowSums(res$A_matrix)==0)) {
      res$A_matrix[rowSums(res$A_matrix)==0,1] <- 1e-5
    }
  }
  else if (met == "DECONnmf") {
    res <- NMF::nmf(x = x, rank = k, method = "snmf/r", seed = 1)
    res <- list(A_matrix=apply(X = res@fit@H, 2, function(z) {
      z/sum(z)}), #explicit STO constraint
                T_matrix=res@fit@W)
    if (any(rowSums(res$A_matrix)==0)) {
      res$A_matrix[rowSums(res$A_matrix)==0,1] <- 1e-5
    }
  }
  time_elapsed=toc()
  
  if (option=="Tmat") {
    hvg=names(sort(apply(ref_profiles,1,var), decreasing=T))[seq(1e3)]
    elt1=ref_profiles[hvg,]
    elt2=res$T_matrix[hvg,]
  }
  else if (option=="Amat") {
    elt1=t(prop_simu)
    elt2=t(res$A_matrix)
  }
  elt1=elt1[,sort(colnames(elt1))]
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

