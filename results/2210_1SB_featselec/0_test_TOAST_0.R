library(TOAST)

## Params
block='rna'
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

## Deconvolution
tmp <- lapply(deconv_methods_sup, function(meth) {
    if (!file.exists(paste0("0_test_TOAST/deconv/deconvbl_",block,"_sup_",
                            meth,".rds"))) {
      print(paste0("Running ",meth))
      deconv_res=run_sup_deconvolution(meth, dat, T_ref)
      saveRDS(deconv_res,paste0("0_test_TOAST/deconv/deconvbl_",block,"_sup_",
                                meth,".rds"))
    }
  })
  
tmp <- lapply(deconv_methods_unsup, function(meth) {
  if (!file.exists(paste0("0_test_TOAST/deconv/deconvbl_",block,"_unsup_",
                          meth,".rds"))) {
    print(paste0("Running ",meth))
    deconv_res=run_unsup_deconvolution(meth, dat, T_ref, option="Amat", prop_simu=Atrue)
    saveRDS(deconv_res,paste0("0_test_TOAST/deconv/deconvbl_",block,"_unsup_",
                              meth,".rds"))
  }
})
