library(TOAST)

## Params
block='met'
if (block=="met") {
deconv_methods_sup <- c("rlr","CIBERSORT","houseman")#,"EMeth")
deconv_methods_unsup <- c("RefFreeEWAS","ICA","EDec","MeDeCom")
source(paste0("../../src/1SB_deconv_met_functions.R"))
}
if (block=="rna") {
deconv_methods_sup <- c("DeconRNASeq","nnls","ols","svr",
                        "CIBERSORT","elastic_net","RLR",
                        "WISP")
deconv_methods_unsup <- c("DECONica","DECONnmf")
source(paste0("../../src/1SB_deconv_rna_functions.R"))
}

## Data
if (block=="met") {
data("RA_100samples")
dat <- RA_100samples$Y_raw
T_ref <- RA_100samples$Blood_ref
Atrue <- t(RA_100samples$Prop)
rm(RA_100samples)
}

if (block=="rna") {
data("CBS_PBMC_array")
dat <- CBS_PBMC_array$mixed_all
T_ref <- CBS_PBMC_array$LM_5
Atrue <- t(CBS_PBMC_array$trueProp)
rm(CBS_PBMC_array)
}

## TOAST
K=nrow(Atrue)
deconv_function <- function(dat, k) {
  res <- fastICA::fastICA(X = dat, n.comp = k, maxit = 1000, tol = 1e-09)
  res$names = row.names(dat)
  weighted.list <- deconica::generate_markers(df = res, n = 30, return = "gene.ranked")
  res$A_weighted = t(deconica::get_scores(res$X, weighted.list, summary = "weighted.mean", na.rm = TRUE))
  colnames(res$A_weighted) = colnames(dat)
  A = abs(res$A_weighted) %*% diag(1/colSums(abs(res$A_weighted)))
  return(t(A))
}
if (nrow(dat)<1000) {nMarker=300} else {nMarker=1000}
toast <- TOAST::csDeconv(dat, K, TotalIter = 10, FUN = deconv_function, nMarker=nMarker)$updatedInx
hvg=TOAST::findRefinx(dat, nmarker=length(toast))

## Deconvolution
tmp <- lapply(deconv_methods_sup, function(meth) {
    if (!file.exists(paste0("0_test_TOAST/deconv/deconv_",block,"_sup_",
                            meth,".rds"))) {
      print(paste0("Running ",meth))
      deconv_res=run_sup_deconvolution(meth, dat[toast,], T_ref[toast,])
      saveRDS(deconv_res,paste0("0_test_TOAST/deconv/deconv_",block,"_sup_",
                                meth,".rds"))
    }
  if (!file.exists(paste0("0_test_TOAST/deconv/deconvbl2_",block,"_sup_",
                          meth,".rds"))) {
    print(paste0("Running ",meth))
    deconv_res=run_sup_deconvolution(meth, dat[hvg,], T_ref[hvg,])
    saveRDS(deconv_res,paste0("0_test_TOAST/deconv/deconvbl2_",block,"_sup_",
                              meth,".rds"))
  }
  })
  
tmp <- lapply(deconv_methods_unsup, function(meth) {
  if (!file.exists(paste0("0_test_TOAST/deconv/deconv_",block,"_unsup_",
                          meth,".rds"))) {
    print(paste0("Running ",meth))
    deconv_res=run_unsup_deconvolution(meth, dat[toast,], T_ref[toast,], option="Amat", prop_simu=Atrue)
    saveRDS(deconv_res,paste0("0_test_TOAST/deconv/deconv_",block,"_unsup_",
                              meth,".rds"))
  }
  if (!file.exists(paste0("0_test_TOAST/deconv/deconvbl2_",block,"_unsup_",
                          meth,".rds"))) {
    print(paste0("Running ",meth))
    deconv_res=run_unsup_deconvolution(meth, dat[hvg,], T_ref[hvg,], option="Amat", prop_simu=Atrue)
    saveRDS(deconv_res,paste0("0_test_TOAST/deconv/deconvbl2_",block,"_unsup_",
                              meth,".rds"))
  }
})

