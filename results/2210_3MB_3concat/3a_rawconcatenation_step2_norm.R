library(pbapply)
source("../../src/2MB_deconv_functions.R")
source("../../src/score_functions.R")

## ----
## Set parameters
## ----
n_hv=5e3
deconv_meth <- c("cibersort","rlr","ica")
score_methods <- c("rmse","mae","pearson")

## ----
## load D, T and Atrue
## ----
path = "../2210_0simu/simulations/rna/"
n_lot = length(list.files(path, pattern = "T_rna_ref"))
n_sims = length(list.files(path, pattern = "sim"))/n_lot
date = strsplit(list.files(path)[1],"_")[[1]][1]
name_data = unname(sapply(list.files(path, pattern = "T_rna_ref"),function(x) rev(strsplit(strsplit(x,"_T_rna")[[1]][1],"_")[[1]])[1]))

Atrue = lapply(sort(list.files(path, pattern = "sim")), function(x) readRDS(paste0("../2210_0simu/simulations/rna/",x))$A_ref)
Atrue = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) Atrue[[10*(lot-1)+sim]])
})
n_celltypes = unique(sapply(Atrue, function(x) nrow(x[[1]])))
n_sample = unique(sapply(Atrue, function(x) ncol(x[[1]])))
rm(path)

D_hv_scale = readRDS(paste0("3a_rawconcatenation_step1/",
                          date,
                          "_DT_hv_scale.rds"))$D
T_hv_scale = readRDS(paste0("3a_rawconcatenation_step1/",
                            date,
                            "_DT_hv_scale.rds"))$Tref

## ----
## Extract hv
## ----
hv_scale_D <- lapply(D_hv_scale, function(x)
  pblapply(x, function(y)
    rownames(y)[TOAST::findRefinx(y, nmarker = n_hv)]))
hv_scale_T <- lapply(T_hv_scale, function(x)
  pblapply(x, function(y)
    rownames(y)[TOAST::findRefinx(y, nmarker = n_hv)]))
sapply(hv_scale_D, function(x)
  sapply(x, function(y) {
    prop=sum(startsWith(y,"ch"),startsWith(y,"cg"),startsWith(y,"rs"))/n_hv
    print(paste0(prop*100,"% are probes"))}))
sapply(hv_scale_T, function(x)
  sapply(x, function(y) {
    prop=sum(startsWith(y,"ch"),startsWith(y,"cg"),startsWith(y,"rs"))/n_hv
    print(paste0(prop*100,"% are probes"))}))

## ----
## Deconvolution
## ----
tmp=lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    print(paste0("Running simulation n°",sim,"/",n_sims))
    lapply(deconv_meth, function(meth) {
      if (!file.exists(paste0("3a_rawconcatenation_step2_norm/deconv/",
                              date,
                              "_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))) {
        print(paste0("Running method ",meth))
        ref_profiles=T_hv_scale[[lot]][[sim]]
        ref_profiles=ref_profiles[,sort(colnames(ref_profiles))]
        deconv_res=run_concat_deconvolution(meth, D_hv_scale[[lot]][[sim]], ref_profiles, Atrue[[lot]][[sim]])
        saveRDS(deconv_res,paste0("3a_rawconcatenation_step2_norm/deconv/",
                                  date,
                                  "_",name_data[lot],
                                  "_Apred_",
                                  meth,
                                  "_",sim_txt,
                                  sim,".rds"))
      }
    })
  })
})
rm(tmp)

## ----
## Score
## ----
score_perf = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth, function(meth) {
      if (file.exists(paste0("3a_rawconcatenation_step2_norm/deconv/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("3a_rawconcatenation_step2_norm/deconv/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$res
        if (!all(is.na(A_pred))) {
        score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue[[lot]][[sim]], A_pred, x))
        score_perf_res_celltype <- sapply(score_methods, function(x)
          sapply(seq(nrow(A_pred)), function(i)
            score_perf(Atrue[[lot]][[sim]][i,], A_pred[i,], x)))
        score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
        score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
        score_perf_res_sample <- sapply(score_methods, function(x)
          sapply(seq(ncol(A_pred)), function(i)
            score_perf(Atrue[[lot]][[sim]][,i], A_pred[,i], x)))
        score_sd_res_sample=apply(score_perf_res_sample, 2, sd)
        score_med_res_sample=apply(score_perf_res_sample, 2, median)
        return(list("perf_g"=score_perf_res,
                    "perf_c"=score_perf_res_celltype,
                    "sd_c"=score_sd_res_celltype,
                    "med_c"=score_med_res_celltype,
                    "perf_s"=score_perf_res_sample,
                    "sd_s"=score_sd_res_sample,
                    "med_s"=score_med_res_sample))
        }
        else {return(list("perf_g"=rep(NA,length(score_methods)),
                          "perf_c"=matrix(NA,ncol=length(score_methods),nrow=n_celltypes[lot]),
                          "sd_c"=rep(NA,length(score_methods)),
                          "med_c"=rep(NA,length(score_methods)),
                          "perf_s"=matrix(NA,ncol=length(score_methods),nrow=n_sample),
                          "sd_s"=rep(NA,length(score_methods)),
                          "med_s"=rep(NA,length(score_methods))))}
      }
      else {print(paste0(meth," not run yet"))}
    })
    names(res)=deconv_meth
    return(res)
  })
})
score_perf_df = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    res1 = lapply(seq(length(deconv_meth)), function(deconv) {
      df = do.call(rbind,score_perf[[lot]][[sim]][[deconv]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(score_methods, each=nrow(df)),
                      "setting"=rep(rownames(df), ncol(df)))
      df$deconv = deconv_meth[deconv]
      return(df)
    })
    res1 = do.call(rbind,res1)
    res1$sim = sim
    return(res1)
  })
  res2 = do.call(rbind,res2)
  res2$dataset = name_data[lot]
  return(res2)
})
score_perf_df = do.call(rbind,score_perf_df)
saveRDS(score_perf_df, paste0("3a_rawconcatenation_step2_norm/perf_scores/",
                              date,
                              "_scores.rds"))
