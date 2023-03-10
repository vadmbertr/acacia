## ----
## Set parameters
## ----
source("../../src/score_functions.R")
block="met" # to vary, either met or rna

## ----
## load A_true and A_pred
## ----
input_path = paste0("../2210_0simu/simulations/",block,"/")
input_path_A_true = sort(list.files(input_path, pattern = "sim"))
input_path_T = sort(list.files(input_path, pattern = paste0("T_",block,"_ref")))
n_lot = length(input_path_T)
n_sims = length(input_path_A_true)/n_lot
date = strsplit(input_path_T[1],"_")[[1]][1]
name_data = unname(sapply(input_path_T,function(x) strsplit(strsplit(x,"_T_")[[1]][1],paste0(date,"_"))[[1]][2]))

A_true = lapply(input_path_A_true, function(x) readRDS(paste0(input_path,x))$A_ref)
A_true = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) A_true[[10*(lot-1)+sim]])
})
n_sample = ncol(A_true[[1]][[1]])
n_celltypes = sapply(A_true, function(x) nrow(x[[1]]))
rm(input_path,input_path_A_true,input_path_T)

input_path_methodssup = paste0("deconv/",block,"/sup/")
input_path_methodssup = sort(list.files(input_path_methodssup))
deconv_methods_sup = unique(sapply(input_path_methodssup, function(x)
  rev(strsplit(x,"_")[[1]])[2]))
deconv_methods_sup[deconv_methods_sup=="net"] = "elastic_net"
rm(input_path_methodssup)

input_path_methodsunsup = paste0("deconv/",block,"/unsup/")
input_path_methodsunsup = sort(list.files(input_path_methodsunsup))
deconv_methods_unsup = unique(sapply(input_path_methodsunsup, function(x)
  rev(strsplit(x,"_")[[1]])[2]))
rm(input_path_methodsunsup)

score_methods <- c("rmse","mae","pearson")

## ----
## Score ref-based
## ----
score_perf_sup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("deconv/",block,"/sup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                             "_",sim_txt,
                              sim,".rds"))) {
        A_pred=readRDS(paste0("deconv/",block,"/sup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        score_perf_res <- sapply(score_methods, function(x) score_perf(A_true[[lot]][[sim]], A_pred, x))
        score_perf_res_celltype <- sapply(score_methods, function(x)
          sapply(seq(nrow(A_pred)), function(i)
            score_perf(A_true[[lot]][[sim]][i,], A_pred[i,], x)))
        score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
        score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
        score_perf_res_sample <- sapply(score_methods, function(x)
          sapply(seq(ncol(A_pred)), function(i)
            score_perf(A_true[[lot]][[sim]][,i], A_pred[,i], x)))
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
      else {print(paste0(meth," not run yet"))}
      })
    names(res)=deconv_methods_sup
    return(res)
    })
  })
score_perf_sup_df = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    res1 = lapply(seq(length(deconv_methods_sup)), function(deconv) {
      df = do.call(rbind,score_perf_sup[[lot]][[sim]][[deconv]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(colnames(df), each=nrow(df)),
                      "setting"=rep(rownames(df), ncol(df)))
      df$deconv = deconv_methods_sup[deconv]
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
score_perf_sup_df = do.call(rbind,score_perf_sup_df)
saveRDS(score_perf_sup_df, paste0("perf_scores/",
                                  date,
                                  "_scores_",
                                  block,
                                  "_sup.rds"))

## ----
## Score ref-free
## ----
score_perf_unsup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("deconv/",block,"/unsup/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("deconv/",block,"/unsup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        score_perf_res <- sapply(score_methods, function(x) score_perf(A_true[[lot]][[sim]], A_pred, x))
        score_perf_res_celltype <- sapply(score_methods, function(x)
          sapply(seq(nrow(A_pred)), function(i)
            score_perf(A_true[[lot]][[sim]][i,], A_pred[i,], x)))
        score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
        score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
        score_perf_res_sample <- sapply(score_methods, function(x)
          sapply(seq(ncol(A_pred)), function(i)
            score_perf(A_true[[lot]][[sim]][,i], A_pred[,i], x)))
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
      else {print(paste0(meth," not run yet"))}
    })
  })
})
score_perf_unsup_df = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    res1 = lapply(seq(length(deconv_methods_unsup)), function(deconv) {
      df = do.call(rbind,score_perf_unsup[[lot]][[sim]][[deconv]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(colnames(df), each=nrow(df)),
                      "setting"=rep(rownames(df), ncol(df)))
      df$deconv = deconv_methods_unsup[deconv]
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
score_perf_unsup_df = do.call(rbind,score_perf_unsup_df)
saveRDS(score_perf_unsup_df, paste0("perf_scores/",
                                  date,
                                  "_scores_",
                                  block,
                                  "_unsup.rds"))
## For datasets with no ground truth, we can compare similarities between methods, as opposed to known similarities
