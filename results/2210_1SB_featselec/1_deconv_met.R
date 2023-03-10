library(pbapply)
set.seed(2022)

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
n_lot = length(T_ref)
n_sims = length(matrix)/length(T_ref)
names(matrix) <- input_path_matrix
names(T_ref) <- input_path_T

if (!file.exists("feat_selec/hvg.rds")) {
  hvg <- pblapply(matrix, function(x) TOAST::findRefinx(x, nmarker = 3e4))
  saveRDS(hvg,"feat_selec/hvg.rds")
}
if (file.exists("feat_selec/hvg.rds")) {hvg <- readRDS("feat_selec/hvg.rds")}
matrix <- pblapply(seq_along(hvg), function(x) matrix[[x]][hvg[[x]],])
T_ref <- lapply(seq_along(hvg), function(x) T_ref[[ceiling(x/10)]][hvg[[x]],])
rm(hvg)

train_targets <- lapply(T_ref, function(x) {
  matrix=diag(ncol(x))
  colnames(matrix)=colnames(x)
  rownames(matrix)=colnames(x)
  as.data.frame(matrix)})

## ----
## Feat selec
## ----
if (!file.exists("feat_selec/TOAST/met.rds")) {
  K=list(4,5,9)
  deconv_function <- function(dat, k) {
    res <- fastICA::fastICA(X = dat, n.comp = k, maxit = 1000, tol = 1e-09)
    res$names = row.names(dat)
    weighted.list <- deconica::generate_markers(df = res, n = 30, return = "gene.ranked")
    res$A_weighted = t(deconica::get_scores(res$X, weighted.list, summary = "weighted.mean", na.rm = TRUE))
    colnames(res$A_weighted) = colnames(dat)
    A = abs(res$A_weighted) %*% diag(1/colSums(abs(res$A_weighted)))
    return(t(A))
  }
  toast <- lapply(seq(n_lot), function(lot)
    pblapply(seq(n_sims), function(sim)
      TOAST::csDeconv(matrix[[(lot-1)*10+sim]], K[[lot]], TotalIter = 10, FUN = deconv_function)$updatedInx))
  saveRDS(toast,"feat_selec/TOAST/met.rds")
}
if (file.exists("feat_selec/TOAST/met.rds")) {toast = readRDS("feat_selec/TOAST/met.rds")}

if (!file.exists("feat_selec/TOAST/met_restricted.rds")) {
  toast_restricted <- lapply(seq(n_lot), function(x)
    lapply(seq(n_sims), function(y)
      TOAST::findRefinx(matrix[[(x-1)*10+y]], nmarker=length(toast[[x]][[y]]))))
  saveRDS(toast_restricted,"feat_selec/TOAST/met_restricted.rds")
}
if (file.exists("feat_selec/TOAST/met_restricted.rds")) {toast_restricted = readRDS("feat_selec/TOAST/met_restricted.rds")}

if (!file.exists("feat_selec/CimpleG/met.rds")) {
  cimpleg <- lapply(seq(11,length(T_ref)), function(x)
    CimpleG::CimpleG(t(T_ref[[x]]), train_targets = train_targets[[x]], target_columns = colnames(train_targets[[x]]), train_only = T)$signatures)
  saveRDS(cimpleg,"feat_selec/CimpleG/met.rds")
}
if (file.exists("feat_selec/CimpleG/met.rds")) {cimpleg = readRDS("feat_selec/CimpleG/met.rds")}

if (!file.exists("feat_selec/CimpleG/met_restricted.rds")) {
  cimpleg_restricted <- lapply(seq(11,length(matrix)), function(x)
    TOAST::findRefinx(matrix[[x]], nmarker=length(cimpleg[[x-10]])))
  saveRDS(cimpleg_restricted,"feat_selec/CimpleG/met_restricted.rds")
}
if (file.exists("feat_selec/CimpleG/met_restricted.rds")) {cimpleg_restricted = readRDS("feat_selec/CimpleG/met_restricted.rds")}

## ----
## Deconvolution MET ref-based TOAST and CimpleG
## ----
tmp=lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    print(paste0("Running simulation n°",sim,"/",n_sims))
    ref_profiles=T_ref[[(lot-1)*10+sim]]
    ref_profiles=ref_profiles[,sort(colnames(ref_profiles))]
    dat=matrix[[(lot-1)*10+sim]]
    lapply(deconv_methods_sup, function(meth) {
      if (lot>1 & !file.exists(paste0("deconv/CimpleG/",block,"/sup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))) {
        print(paste0("Running ",meth))
        deconv_res=run_sup_deconvolution(meth, dat[cimpleg[[(lot-2)*10+sim]],], ref_profiles[cimpleg[[(lot-2)*10+sim]],])
        deconv_resbl2=run_sup_deconvolution(meth, dat[cimpleg_restricted[[(lot-2)*10+sim]],], ref_profiles[cimpleg_restricted[[(lot-2)*10+sim]],])
        saveRDS(deconv_res$res,paste0("deconv/CimpleG/",block,"/sup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        saveRDS(unname(deconv_res$time_elapsed),paste0("timing/CimpleG/",block,"/sup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_timing_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        saveRDS(deconv_resbl2$res,paste0("deconv/CimpleG/",block,"/sup/",
                                      strsplit(input_path_T[lot],pattern)[[1]][1],
                                      "_Apredbl2_",
                                      meth,
                                      "_",sim_txt,
                                      sim,".rds"))
        saveRDS(unname(deconv_resbl2$time_elapsed),paste0("timing/CimpleG/",block,"/sup/",
                                                       strsplit(input_path_T[lot],pattern)[[1]][1],
                                                       "_timingbl2_",
                                                       meth,
                                                       "_",sim_txt,
                                                       sim,".rds"))
      }
      if (!file.exists(paste0("deconv/TOAST/",block,"/sup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))) {
        print(paste0("Running ",meth))
        deconv_res=run_sup_deconvolution(meth, dat[toast[[lot]][[sim]],], ref_profiles[toast[[lot]][[sim]],])
        deconv_resbl2=run_sup_deconvolution(meth, dat[toast_restricted[[lot]][[sim]],], ref_profiles[toast_restricted[[lot]][[sim]],])
        saveRDS(deconv_res$res,paste0("deconv/TOAST/",block,"/sup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        saveRDS(unname(deconv_res$time_elapsed),paste0("timing/TOAST/",block,"/sup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_timing_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        saveRDS(deconv_resbl2$res,paste0("deconv/TOAST/",block,"/sup/",
                                      strsplit(input_path_T[lot],pattern)[[1]][1],
                                      "_Apredbl2_",
                                      meth,
                                      "_",sim_txt,
                                      sim,".rds"))
        saveRDS(unname(deconv_resbl2$time_elapsed),paste0("timing/TOAST/",block,"/sup/",
                                                       strsplit(input_path_T[lot],pattern)[[1]][1],
                                                       "_timingbl2_",
                                                       meth,
                                                       "_",sim_txt,
                                                       sim,".rds"))
      }
    })
  })
})

## ----
## Deconvolution MET ref-free TOAST and CimpleG
## ----
tmp=lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    print(paste0("Running simulation n°",sim,"/",n_sims))
    prop_simu=Aref[[(lot-1)*10+sim]]
    dat=matrix[[(lot-1)*10+sim]]
    ref_profiles=T_ref[[(lot-1)*10+sim]]
    lapply(deconv_methods_unsup, function(meth) {
      if (lot>1 & !file.exists(paste0("deconv/CimpleG/",block,"/unsup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))) {
        print(paste0("Running ",meth))
        deconv_res=run_unsup_deconvolution(meth, dat[cimpleg[[(lot-2)*10+sim]],], ref_profiles[cimpleg[[(lot-2)*10+sim]],], "Amat", prop_simu=prop_simu)
        deconv_resbl2=run_unsup_deconvolution(meth, dat[cimpleg_restricted[[(lot-2)*10+sim]],], ref_profiles[cimpleg_restricted[[(lot-2)*10+sim]],], "Amat", prop_simu=prop_simu)
        saveRDS(deconv_res$res,paste0("deconv/CimpleG/",block,"/unsup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        saveRDS(unname(deconv_res$time_elapsed),paste0("timing/CimpleG/",block,"/unsup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_timing_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        saveRDS(deconv_resbl2$res,paste0("deconv/CimpleG/",block,"/unsup/",
                                      strsplit(input_path_T[lot],pattern)[[1]][1],
                                      "_Apredbl2_",
                                      meth,
                                      "_",sim_txt,
                                      sim,".rds"))
        saveRDS(unname(deconv_resbl2$time_elapsed),paste0("timing/CimpleG/",block,"/unsup/",
                                                       strsplit(input_path_T[lot],pattern)[[1]][1],
                                                       "_timingbl2_",
                                                       meth,
                                                       "_",sim_txt,
                                                       sim,".rds"))
      }
      
      if (!file.exists(paste0("deconv/TOAST/",block,"/unsup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))) {
        print(paste0("Running ",meth))
        deconv_res=run_unsup_deconvolution(meth, dat[toast[[lot]][[sim]],], ref_profiles[toast[[lot]][[sim]],], "Amat", prop_simu=prop_simu)
        deconv_resbl2=run_unsup_deconvolution(meth, dat[toast_restricted[[lot]][[sim]],], ref_profiles[toast_restricted[[lot]][[sim]],], "Amat", prop_simu=prop_simu)
        saveRDS(deconv_res$res,paste0("deconv/TOAST/",block,"/unsup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        saveRDS(unname(deconv_res$time_elapsed),paste0("timing/TOAST/",block,"/unsup/",
                              strsplit(input_path_T[lot],pattern)[[1]][1],
                              "_timing_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        saveRDS(deconv_resbl2$res,paste0("deconv/TOAST/",block,"/unsup/",
                                      strsplit(input_path_T[lot],pattern)[[1]][1],
                                      "_Apredbl2_",
                                      meth,
                                      "_",sim_txt,
                                      sim,".rds"))
        saveRDS(unname(deconv_resbl2$time_elapsed),paste0("timing/TOAST/",block,"/unsup/",
                                                       strsplit(input_path_T[lot],pattern)[[1]][1],
                                                       "_timingbl2_",
                                                       meth,
                                                       "_",sim_txt,
                                                       sim,".rds"))
      }
    })
  })
})

## For datasets with no ground truth, we can compare similarities between methods, as opposed to known similarities
