## ----
## Set parameters
## ----
library(CimpleG)
source("../../src/score_functions.R")
block="met"

## ----
## load
## ----
data(train_data)
data("train_targets")
dat <- t(train_data)
annot <- train_targets[,c(1,2)]
T_ref <- t(WGCNA::collapseRows(t(dat), annot$CELL_TYPE,  rownames(t(dat)),method="Average")$datET)
Atrue <- t(train_targets[,c(3:16)])
rownames(Atrue) <- gsub("CELL_TYPE_","",rownames(Atrue))
colnames(Atrue) <- annot$GSM
Atrue <- Atrue[rowSums(Atrue)>0,]
T_ref <- T_ref[sort(rownames(T_ref)),]
T_ref <- T_ref[,sort(colnames(T_ref))]
Atrue <- Atrue[sort(rownames(Atrue)),]
Atrue <- Atrue[,sort(colnames(Atrue))]
dat <- dat[sort(rownames(dat)),]
dat <- dat[,sort(colnames(dat))]
rm(annot)

input_path_methods = "0_test_CimpleG/deconv/"
input_path_methodssup = sort(list.files(input_path_methods, pattern = paste0("deconv_",block,"_sup")))
deconv_methods_sup = unique(sapply(input_path_methodssup, function(x)
  strsplit(rev(strsplit(x,"_")[[1]])[1],".rds")[[1]][1]))
deconv_methods_sup[deconv_methods_sup=="net"] = "elastic_net"
rm(input_path_methodssup)

input_path_methodsunsup = sort(list.files(input_path_methods, pattern = paste0("deconv_",block,"_unsup")))
deconv_methods_unsup = unique(sapply(input_path_methodsunsup, function(x)
  strsplit(rev(strsplit(x,"_")[[1]])[1],".rds")[[1]][1]))
rm(input_path_methodsunsup,input_path_methods)

score_methods <- c("rmse","mae","pearson")

## ----
## Score ref-based
## ----
score_perf_sup = lapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("0_test_CimpleG/deconv/deconv_",block,"_sup_",
                             meth,".rds"))) {
        A_pred=readRDS(paste0("0_test_CimpleG/deconv/deconv_",block,"_sup_",
                              meth,".rds"))$res
        if (all(is.na(A_pred))) {return(list("perf_g"=rep(NA,length(score_methods)),
                                          "perf_c"=matrix(NA, ncol=length(score_methods), nrow=nrow(Atrue)),
                                          "sd_c"=rep(NA,length(score_methods)),
                                          "med_c"=rep(NA,length(score_methods)),
                                          "perf_s"=matrix(NA, ncol=length(score_methods), nrow=ncol(Atrue)),
                                          "sd_s"=rep(NA,length(score_methods)),
                                          "med_s"=rep(NA,length(score_methods))))}
        else {
          score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue, A_pred, x))
          score_perf_res_celltype <- sapply(score_methods, function(x)
            sapply(seq(nrow(A_pred)), function(i)
              score_perf(Atrue[i,], A_pred[i,], x)))
          score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
          score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
          score_perf_res_sample <- sapply(score_methods, function(x)
            sapply(seq(ncol(A_pred)), function(i)
              score_perf(Atrue[,i], A_pred[,i], x)))
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
        }
      else {print(paste0(meth," not run yet"))}
  })
score_perf_sup_df = lapply(seq(length(deconv_methods_sup)), function(deconv) {
      df = do.call(rbind,score_perf_sup[[deconv]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(nrow(Atrue))),"sd_c","med_c",paste0("perf_s",seq(ncol(Atrue))),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(score_methods, each=nrow(df)),
                      "setting"=rep(rownames(df), ncol(df)),
                      "feat_selec"='CimpleG')
      df$deconv = deconv_methods_sup[deconv]
      return(df)
    })
score_perf_sup_df = do.call(rbind,score_perf_sup_df)
saveRDS(score_perf_sup_df, paste0("0_test_CimpleG/scores/",
                                        "scores_",block,"_sup.rds"))

score_perf_sup_bl = lapply(deconv_methods_sup, function(meth) {
  if (file.exists(paste0("0_test_CimpleG/deconv/deconvbl_",block,"_sup_",
                         meth,".rds"))) {
    A_pred=readRDS(paste0("0_test_CimpleG/deconv/deconvbl_",block,"_sup_",
                          meth,".rds"))$res
    if (all(is.na(A_pred))) {return(list("perf_g"=rep(NA,length(score_methods)),
                                         "perf_c"=matrix(NA, ncol=length(score_methods), nrow=nrow(Atrue)),
                                         "sd_c"=rep(NA,length(score_methods)),
                                         "med_c"=rep(NA,length(score_methods)),
                                         "perf_s"=matrix(NA, ncol=length(score_methods), nrow=ncol(Atrue)),
                                         "sd_s"=rep(NA,length(score_methods)),
                                         "med_s"=rep(NA,length(score_methods))))}
    else {
      score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue, A_pred, x))
      score_perf_res_celltype <- sapply(score_methods, function(x)
        sapply(seq(nrow(A_pred)), function(i)
          score_perf(Atrue[i,], A_pred[i,], x)))
      score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
      score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
      score_perf_res_sample <- sapply(score_methods, function(x)
        sapply(seq(ncol(A_pred)), function(i)
          score_perf(Atrue[,i], A_pred[,i], x)))
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
  }
  else {print(paste0(meth," not run yet"))}
})
score_perf_sup_df_bl = lapply(seq(length(deconv_methods_sup)), function(deconv) {
  df = do.call(rbind,score_perf_sup_bl[[deconv]])
  rownames(df) = c("perf_g",paste0("perf_c",seq(nrow(Atrue))),"sd_c","med_c",paste0("perf_s",seq(ncol(Atrue))),"sd_s","med_s")
  df = data.frame("values"=c(df),
                  "score"=rep(score_methods, each=nrow(df)),
                  "setting"=rep(rownames(df), ncol(df)),
                  "feat_selec"='None')
  df$deconv = deconv_methods_sup[deconv]
  return(df)
})
score_perf_sup_df_bl = do.call(rbind,score_perf_sup_df_bl)
saveRDS(score_perf_sup_df_bl, paste0("0_test_CimpleG/scores/",
                                        "scoresbl_",block,"_sup.rds"))

score_perf_sup_bl2 = lapply(deconv_methods_sup, function(meth) {
  if (file.exists(paste0("0_test_CimpleG/deconv/deconvbl2_",block,"_sup_",
                         meth,".rds"))) {
    A_pred=readRDS(paste0("0_test_CimpleG/deconv/deconvbl2_",block,"_sup_",
                          meth,".rds"))$res
    if (all(is.na(A_pred))) {return(list("perf_g"=rep(NA,length(score_methods)),
                                         "perf_c"=matrix(NA, ncol=length(score_methods), nrow=nrow(Atrue)),
                                         "sd_c"=rep(NA,length(score_methods)),
                                         "med_c"=rep(NA,length(score_methods)),
                                         "perf_s"=matrix(NA, ncol=length(score_methods), nrow=ncol(Atrue)),
                                         "sd_s"=rep(NA,length(score_methods)),
                                         "med_s"=rep(NA,length(score_methods))))}
    else {
      score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue, A_pred, x))
      score_perf_res_celltype <- sapply(score_methods, function(x)
        sapply(seq(nrow(A_pred)), function(i)
          score_perf(Atrue[i,], A_pred[i,], x)))
      score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
      score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
      score_perf_res_sample <- sapply(score_methods, function(x)
        sapply(seq(ncol(A_pred)), function(i)
          score_perf(Atrue[,i], A_pred[,i], x)))
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
  }
  else {print(paste0(meth," not run yet"))}
})
score_perf_sup_df_bl2 = lapply(seq(length(deconv_methods_sup)), function(deconv) {
  df = do.call(rbind,score_perf_sup_bl2[[deconv]])
  rownames(df) = c("perf_g",paste0("perf_c",seq(nrow(Atrue))),"sd_c","med_c",paste0("perf_s",seq(ncol(Atrue))),"sd_s","med_s")
  df = data.frame("values"=c(df),
                  "score"=rep(score_methods, each=nrow(df)),
                  "setting"=rep(rownames(df), ncol(df)),
                  "feat_selec"='None_restricted')
  df$deconv = deconv_methods_sup[deconv]
  return(df)
})
score_perf_sup_df_bl2 = do.call(rbind,score_perf_sup_df_bl2)
saveRDS(score_perf_sup_df_bl2, paste0("0_test_CimpleG/scores/",
                                           "scoresbl2_",block,"_sup.rds"))

## ----
## Score ref-free
## ----
score_perf_unsup = lapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("0_test_CimpleG/deconv/deconv_",block,"_unsup_",
                             meth,".rds"))) {
        A_pred=readRDS(paste0("0_test_CimpleG/deconv/deconv_",block,"_unsup_",
                              meth,".rds"))$res
        if (all(is.na(A_pred))) {return(list("perf_g"=rep(NA,length(score_methods)),
                                             "perf_c"=matrix(NA, ncol=length(score_methods), nrow=nrow(Atrue)),
                                             "sd_c"=rep(NA,length(score_methods)),
                                             "med_c"=rep(NA,length(score_methods)),
                                             "perf_s"=matrix(NA, ncol=length(score_methods), nrow=ncol(Atrue)),
                                             "sd_s"=rep(NA,length(score_methods)),
                                             "med_s"=rep(NA,length(score_methods))))}
        else {
          score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue, A_pred, x))
          score_perf_res_celltype <- sapply(score_methods, function(x)
            sapply(seq(nrow(A_pred)), function(i)
              score_perf(Atrue[i,], A_pred[i,], x)))
          score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
          score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
          score_perf_res_sample <- sapply(score_methods, function(x)
            sapply(seq(ncol(A_pred)), function(i)
              score_perf(Atrue[,i], A_pred[,i], x)))
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
      }
      else {print(paste0(meth," not run yet"))}
  })
score_perf_unsup_df = lapply(seq(length(deconv_methods_unsup)), function(deconv) {
      df = do.call(rbind,score_perf_unsup[[deconv]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(nrow(Atrue))),"sd_c","med_c",paste0("perf_s",seq(ncol(Atrue))),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(score_methods, each=nrow(df)),
                      "setting"=rep(rownames(df), ncol(df)),
                      "feat_selec"='CimpleG')
      df$deconv = deconv_methods_unsup[deconv]
      return(df)
      })
score_perf_unsup_df = do.call(rbind,score_perf_unsup_df)
saveRDS(score_perf_unsup_df, paste0("0_test_CimpleG/scores/",
                                          "scores_",block,"_unsup.rds"))

score_perf_unsup_bl = lapply(deconv_methods_unsup, function(meth) {
  if (file.exists(paste0("0_test_CimpleG/deconv/deconvbl_",block,"_unsup_",
                         meth,".rds"))) {
    A_pred=readRDS(paste0("0_test_CimpleG/deconv/deconvbl_",block,"_unsup_",
                          meth,".rds"))$res
    if (all(is.na(A_pred))) {return(list("perf_g"=rep(NA,length(score_methods)),
                                         "perf_c"=matrix(NA, ncol=length(score_methods), nrow=nrow(Atrue)),
                                         "sd_c"=rep(NA,length(score_methods)),
                                         "med_c"=rep(NA,length(score_methods)),
                                         "perf_s"=matrix(NA, ncol=length(score_methods), nrow=ncol(Atrue)),
                                         "sd_s"=rep(NA,length(score_methods)),
                                         "med_s"=rep(NA,length(score_methods))))}
    else {
      score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue, A_pred, x))
      score_perf_res_celltype <- sapply(score_methods, function(x)
        sapply(seq(nrow(A_pred)), function(i)
          score_perf(Atrue[i,], A_pred[i,], x)))
      score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
      score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
      score_perf_res_sample <- sapply(score_methods, function(x)
        sapply(seq(ncol(A_pred)), function(i)
          score_perf(Atrue[,i], A_pred[,i], x)))
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
  }
  else {print(paste0(meth," not run yet"))}
})
score_perf_unsup_df_bl = lapply(seq(length(deconv_methods_unsup)), function(deconv) {
  df = do.call(rbind,score_perf_unsup_bl[[deconv]])
  rownames(df) = c("perf_g",paste0("perf_c",seq(nrow(Atrue))),"sd_c","med_c",paste0("perf_s",seq(ncol(Atrue))),"sd_s","med_s")
  df = data.frame("values"=c(df),
                  "score"=rep(score_methods, each=nrow(df)),
                  "setting"=rep(rownames(df), ncol(df)),
                  "feat_selec"='None')
  df$deconv = deconv_methods_unsup[deconv]
  return(df)
})
score_perf_unsup_df_bl = do.call(rbind,score_perf_unsup_df_bl)
saveRDS(score_perf_unsup_df_bl, paste0("0_test_CimpleG/scores/",
                                          "scoresbl_",block,"_unsup.rds"))

score_perf_unsup_bl2 = lapply(deconv_methods_unsup, function(meth) {
  if (file.exists(paste0("0_test_CimpleG/deconv/deconvbl2_",block,"_unsup_",
                         meth,".rds"))) {
    A_pred=readRDS(paste0("0_test_CimpleG/deconv/deconvbl2_",block,"_unsup_",
                          meth,".rds"))$res
    if (all(is.na(A_pred))) {return(list("perf_g"=rep(NA,length(score_methods)),
                                         "perf_c"=matrix(NA, ncol=length(score_methods), nrow=nrow(Atrue)),
                                         "sd_c"=rep(NA,length(score_methods)),
                                         "med_c"=rep(NA,length(score_methods)),
                                         "perf_s"=matrix(NA, ncol=length(score_methods), nrow=ncol(Atrue)),
                                         "sd_s"=rep(NA,length(score_methods)),
                                         "med_s"=rep(NA,length(score_methods))))}
    else {
      score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue, A_pred, x))
      score_perf_res_celltype <- sapply(score_methods, function(x)
        sapply(seq(nrow(A_pred)), function(i)
          score_perf(Atrue[i,], A_pred[i,], x)))
      score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
      score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
      score_perf_res_sample <- sapply(score_methods, function(x)
        sapply(seq(ncol(A_pred)), function(i)
          score_perf(Atrue[,i], A_pred[,i], x)))
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
  }
  else {print(paste0(meth," not run yet"))}
})
score_perf_unsup_df_bl2 = lapply(seq(length(deconv_methods_unsup)), function(deconv) {
  df = do.call(rbind,score_perf_unsup_bl2[[deconv]])
  rownames(df) = c("perf_g",paste0("perf_c",seq(nrow(Atrue))),"sd_c","med_c",paste0("perf_s",seq(ncol(Atrue))),"sd_s","med_s")
  df = data.frame("values"=c(df),
                  "score"=rep(score_methods, each=nrow(df)),
                  "setting"=rep(rownames(df), ncol(df)),
                  "feat_selec"='None_restricted')
  df$deconv = deconv_methods_unsup[deconv]
  return(df)
})
score_perf_unsup_df_bl2 = do.call(rbind,score_perf_unsup_df_bl2)
saveRDS(score_perf_unsup_df_bl2, paste0("0_test_CimpleG/scores/",
                                             "scoresbl2_",block,"_unsup.rds"))

