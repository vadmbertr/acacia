library(parallel)

## ----
## Set parameters
## ----
block="met"
deconv_methods_sup <- c("rlr","CIBERSORT","houseman")#,"EMeth")
deconv_methods_unsup <- c("RefFreeEWAS","ICA","EDec","MeDeCom")
pattern = paste0("_T_",block,"_ref")
source(paste0("../../src/1SB_deconv_",block,"_functions.R"))


## ----
## load simulated D data, along with A_true and T_ref
## ----
input_path = paste0("../2210_0simu/simulations/",block,"/")
input_path_matrix = sort(list.files(input_path, pattern = "sim"))
input_path_T = sort(list.files(input_path, pattern = paste0("T_",block,"_ref")))
matrix = lapply(input_path_matrix, function(x) readRDS(paste0(input_path,x))$D_met_sim)
Aref = lapply(input_path_matrix, function(x) readRDS(paste0(input_path,x))$A_ref)
T_ref = lapply(input_path_T, function(x) readRDS(paste0(input_path,x)))
T_ref = lapply(T_ref, as.data.frame)
names(matrix) <- input_path_matrix
names(T_ref) <- input_path_T

## ----
## Deconvolution MET ref-based
## ----
tmp=lapply(seq(length(T_ref)), function (lot) {
  print(paste0("Running dataset n°",lot,"/",length(T_ref)))
  lapply(seq(length(matrix)/length(T_ref)), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    print(paste0("Running simulation n°",sim,"/",length(matrix)/length(T_ref)))
    ref_profiles=T_ref[[lot]]
    ref_profiles=ref_profiles[,sort(colnames(ref_profiles))]
    dat=matrix[[(lot-1)*10+sim]]
    hvg=TOAST::findRefinx(dat, nmarker = 1e3)
    lapply(deconv_methods_sup, function(meth) {
      if (!file.exists(paste0("deconv/",block,"/sup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))) {
        print(paste0("Running ",meth))
        deconv_res=run_sup_deconvolution(meth, dat, ref_profiles, hvg)
        A_pred=deconv_res$res
        timing=unname(deconv_res$time_elapsed)
        saveRDS(A_pred,paste0("deconv/",block,"/sup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        saveRDS(timing,paste0("timing/",block,"/sup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_timing_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
      }
    })
  })
})

## ----
## Deconvolution MET ref-free : restricted to 3e4 probes
## ----
tmp=lapply(seq(length(T_ref)), function (lot) {
  print(paste0("Running dataset n°",lot,"/",length(T_ref)))
  lapply(seq(length(matrix)/length(T_ref)), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    print(paste0("Running simulation n°",sim,"/",length(matrix)/length(T_ref)))
    prop_simu=Aref[[(lot-1)*10+sim]]
    dat=matrix[[(lot-1)*10+sim]]
    ref_profiles=T_ref[[lot]]
    if (nrow(dat)>3e4) {
      hvg=TOAST::findRefinx(dat, nmarker = 3e4)
      dat=dat[hvg,]
      ref_profiles=ref_profiles[hvg,]
    }
    hvg=TOAST::findRefinx(dat, nmarker = 1e3)
    lapply(deconv_methods_unsup, function(meth) {
      if (!file.exists(paste0("deconv/",block,"/unsup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))) {
        print(paste0("Running ",meth))
        deconv_res=run_unsup_deconvolution(meth, dat, ref_profiles, "Amat", prop_simu=prop_simu, hvg=hvg)
        A_pred=deconv_res$res
        timing=unname(deconv_res$time_elapsed)
        saveRDS(A_pred,paste0("deconv/",block,"/unsup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        saveRDS(timing,paste0("timing/",block,"/unsup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_timing_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
      }
    })
  })
})

## For datasets with no ground truth, we can compare similarities between methods, as opposed to known similarities
