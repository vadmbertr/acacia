## ----
## Set parameters
## ----
source("../../src/2MB_deconv_functions.R")
deconv_meth <- c("cibersort","rlr","ica")
source("../../src/score_functions.R")
score_methods <- c("rmse","mae","pearson")

## ----
## Load data
## ----
input_path = "../2210_0simu/simulations/rna/"
input_path_matrix = sort(list.files(input_path, pattern = "sim"))
date = strsplit(input_path_matrix[1],"_")[[1]][1]
Atrue = lapply(input_path_matrix, function(x) readRDS(paste0(input_path,x))$A_ref)
n_celltypes = unique(sapply(Atrue,nrow))
n_sample = unique(sapply(Atrue,ncol))
input_path_T2 = sort(list.files(input_path, pattern = "_ref.rds"))
input_path_T1 = sort(list.files("../2210_0simu/simulations/met/", pattern = "_ref.rds"))
n_lot = length(input_path_T2)
n_sims = length(input_path_matrix)/n_lot
Tref1 = lapply(input_path_T1, function(x) data.frame(readRDS(paste0("../2210_0simu/simulations/met/",x))))
Tref2 = lapply(input_path_T2, function(x) data.frame(readRDS(paste0("../2210_0simu/simulations/rna/",x))))
name_data = unname(sapply(input_path_T2,function(x) rev(strsplit(strsplit(x,"_T_")[[1]][1],"_")[[1]])[1]))
rm(input_path,input_path_T1,input_path_T2,input_path_matrix)
Atrue = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) Atrue[[10*(lot-1)+sim]])
})
Dtoast <- readRDS(paste0("3b_featureselection_step1/TOAST/",date,"_D_toast.rds"))
toast_hvmet <- readRDS("../2210_1SB_featselec/feat_selec/hvg.rds")
toast_feat <- list(toast_met=readRDS("../2210_1SB_featselec/feat_selec/TOAST/met.rds"),
                   toast_rna=readRDS("../2210_1SB_featselec/feat_selec/TOAST/rna.rds"))
Tref1_toast <- lapply(seq(n_lot), function(lot)
  lapply(seq(n_sims), function(sim)
    Tref1[[lot]][toast_hvmet[[lot]][[sim]],][toast_feat$toast_met[[lot]][[sim]],]))
Tref2_toast <- lapply(seq(n_lot), function(lot)
  lapply(seq(n_sims), function(sim)
    Tref2[[lot]][toast_feat$toast_rna[[lot]][[sim]],]))
Treftoast <- lapply(seq(n_lot), function(lot)
  lapply(seq(n_sims), function(sim) {
    mat1 <- Tref1_toast[[lot]][[sim]]
    mat2 <- Tref2_toast[[lot]][[sim]]
    rownames(mat1) <- paste0("met_",rownames(mat1))
    rownames(mat2) <- paste0("rna_",rownames(mat2))
    mat <- rbind(mat1,mat2)
    mat}))
rm(toast_hvmet,toast_feat,Tref1_toast,Tref2_toast)

Dcimpleg <- readRDS(paste0("3b_featureselection_step1/CimpleG/",date,"_D_cimpleg.rds"))
toast_hvmet <- readRDS("../2210_1SB_featselec/feat_selec/hvg.rds")
toast_feat <- list(toast_met=readRDS("../2210_1SB_featselec/feat_selec/TOAST/met.rds"),
                   toast_rna=readRDS("../2210_1SB_featselec/feat_selec/TOAST/rna.rds"))
cimpleg_feat <- readRDS("../2210_1SB_featselec/feat_selec/CimpleG/met.rds")
cimpleg_feat <- lapply(seq(n_lot-1), function(lot) {
  lapply(seq(n_sims), function(sim) cimpleg_feat[[10*(lot-1)+sim]])
})
Tref1_cimpleg <- c(NA,lapply(seq(2,n_lot), function(lot)
  lapply(seq(n_sims), function(sim)
    Tref1[[lot]][toast_hvmet[[lot]][[sim]],][cimpleg_feat[[lot-1]][[sim]],])))
Tref2_cimpleg <- lapply(seq(n_lot), function(lot)
  lapply(seq(n_sims), function(sim)
    Tref2[[lot]][toast_feat$toast_rna[[lot]][[sim]],]))
Trefcimpleg <- c(NA,lapply(seq(2,n_lot), function(lot)
  lapply(seq(n_sims), function(sim) {
    mat1 <- Tref1_cimpleg[[lot]][[sim]]
    mat2 <- Tref2_cimpleg[[lot]][[sim]]
    rownames(mat1) <- paste0("met_",rownames(mat1))
    rownames(mat2) <- paste0("rna_",rownames(mat2))
    mat <- rbind(mat1,mat2)
    mat})))
rm(toast_hvmet,toast_feat,Tref1_cimpleg,Tref2_cimpleg)

## ----
## Run deconvolution for TOAST
## ----
tmp=lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    print(paste0("Running simulation n°",sim,"/",n_sims))
    lapply(deconv_meth, function(meth) {
      if (!file.exists(paste0("3b_featureselection_step2/TOAST/deconv/",
                              date,
                              "_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))) {
        print(paste0("Running method ",meth))
        ref_profiles=Treftoast[[lot]][[sim]]
        ref_profiles=ref_profiles[,sort(colnames(ref_profiles))]
        deconv_res=run_toast_deconvolution(meth, Dtoast[[lot]][[sim]], ref_profiles, Atrue[[lot]][[sim]])
        saveRDS(deconv_res,paste0("3b_featureselection_step2/TOAST/deconv/",
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
## Score deconvolution for TOAST
## ----
score_perf_t = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth, function(meth) {
      if (file.exists(paste0("3b_featureselection_step2/TOAST/deconv/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("3b_featureselection_step2/TOAST/deconv/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$res
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
      else {print(paste0(meth," not run yet"))}
    })
    names(res)=deconv_meth
    return(res)
  })
})
score_perf_df_t = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    res1 = lapply(seq(length(deconv_meth)), function(deconv) {
      df = do.call(rbind,score_perf_t[[lot]][[sim]][[deconv]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(colnames(df), each=nrow(df)),
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
score_perf_df_t = do.call(rbind,score_perf_df_t)
saveRDS(score_perf_df_t, paste0("3b_featureselection_step2/TOAST/perf_scores/",
                              date,
                              "_scores.rds"))

## ----
## Run deconvolution for CimpleG
## ----
tmp=lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    print(paste0("Running simulation n°",sim,"/",n_sims))
    lapply(deconv_meth, function(meth) {
      if (!file.exists(paste0("3b_featureselection_step2/CimpleG/deconv/",
                              date,
                              "_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))) {
        print(paste0("Running method ",meth))
        if (!all(is.na(Dcimpleg[[lot]]))) {
          ref_profiles=Trefcimpleg[[lot]][[sim]]
          ref_profiles=ref_profiles[,sort(colnames(ref_profiles))]
          deconv_res=run_toast_deconvolution(meth, Dcimpleg[[lot]][[sim]], ref_profiles, Atrue[[lot]][[sim]])
        }
        else {deconv_res=list(res=NA,time_elapsed=NA)}
        saveRDS(deconv_res,paste0("3b_featureselection_step2/CimpleG/deconv/",
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
## Score deconvolution for CimpleG
## ----
score_perf_c = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth, function(meth) {
      if (file.exists(paste0("3b_featureselection_step2/CimpleG/deconv/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("3b_featureselection_step2/CimpleG/deconv/",
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
score_perf_df_c = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    res1 = lapply(seq_along(deconv_meth), function(deconv) {
      df = do.call(rbind,score_perf_c[[lot]][[sim]][[deconv]])
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
score_perf_df_c = do.call(rbind,score_perf_df_c)
saveRDS(score_perf_df_c, paste0("3b_featureselection_step2/CimpleG/perf_scores/",
                              date,
                              "_scores.rds"))
