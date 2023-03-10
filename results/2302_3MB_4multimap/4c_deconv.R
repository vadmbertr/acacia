## ----
## Set parameters
## ----
source("../../src/2MB_deconv_functions.R")
deconv_meth_sup <- c("cibersort","rlr")
deconv_meth_unsup <- c("ica","MeDeCom")
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

list_adata_lot1 = list.files("4a_multiMAP/", pattern="dataset1")
list_adata_lot2 = list.files("4a_multiMAP/", pattern="dataset2")
list_adata_lot3 = list.files("4a_multiMAP/", pattern="dataset3")

adata_lot1 = lapply(list_adata_lot1, function(x) {
  file = anndata::read_h5ad(paste0("4a_multiMAP/",x))
  obs = file$obs
  obsm = file$obsm
  df = data.frame("MultiMAP1"=obsm$X_multimap[,1],
                  "MultiMAP2"=obsm$X_multimap[,2],
                  "MultiMAP3"=obsm$X_multimap[,3],
                  "MultiMAP4"=obsm$X_multimap[,4],
                  "MultiMAP5"=obsm$X_multimap[,5],
                  "MultiMAP6"=obsm$X_multimap[,6],
                  "MultiMAP7"=obsm$X_multimap[,7],
                  "MultiMAP8"=obsm$X_multimap[,8],
                  "MultiMAP9"=obsm$X_multimap[,9],
                  "MultiMAP10"=obsm$X_multimap[,10],
                  "Source"=obs$source,
                  "Type"=obs$type,
                  "SampleNb"=obs$SampleNb)
  df})
adata_lot2 = lapply(list_adata_lot2, function(x) {
  file = anndata::read_h5ad(paste0("4a_multiMAP/",x))
  obs = file$obs
  obsm = file$obsm
  df = data.frame("MultiMAP1"=obsm$X_multimap[,1],
                  "MultiMAP2"=obsm$X_multimap[,2],
                  "MultiMAP3"=obsm$X_multimap[,3],
                  "MultiMAP4"=obsm$X_multimap[,4],
                  "MultiMAP5"=obsm$X_multimap[,5],
                  "MultiMAP6"=obsm$X_multimap[,6],
                  "MultiMAP7"=obsm$X_multimap[,7],
                  "MultiMAP8"=obsm$X_multimap[,8],
                  "MultiMAP9"=obsm$X_multimap[,9],
                  "MultiMAP10"=obsm$X_multimap[,10],
                  "Source"=obs$source,
                  "Type"=obs$type,
                  "SampleNb"=obs$SampleNb)
  df})
adata_lot3 = lapply(list_adata_lot3, function(x) {
  file = anndata::read_h5ad(paste0("4a_multiMAP/",x))
  obs = file$obs
  obsm = file$obsm
  df = data.frame("MultiMAP1"=obsm$X_multimap[,1],
                  "MultiMAP2"=obsm$X_multimap[,2],
                  "MultiMAP3"=obsm$X_multimap[,3],
                  "MultiMAP4"=obsm$X_multimap[,4],
                  "MultiMAP5"=obsm$X_multimap[,5],
                  "MultiMAP6"=obsm$X_multimap[,6],
                  "MultiMAP7"=obsm$X_multimap[,7],
                  "MultiMAP8"=obsm$X_multimap[,8],
                  "MultiMAP9"=obsm$X_multimap[,9],
                  "MultiMAP10"=obsm$X_multimap[,10],
                  "Source"=obs$source,
                  "Type"=obs$type,
                  "SampleNb"=obs$SampleNb)
  df})
rm(list_adata_lot1,list_adata_lot2,list_adata_lot3)
n_lot = 3
n_sims = length(adata_lot1)
Atrue = lapply(seq(n_lot), function(lot)
  lapply(seq(n_sims), function(sim) Atrue[[(lot-1)*10+sim]]))
adata = list(adata_lot1, adata_lot2, adata_lot3)
rm(adata_lot1, adata_lot2, adata_lot3)
name_data = paste0("dataset", seq(n_lot))

## ----
## Run deconvolution
## ----
tmp=lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    print(paste0("Running simulation n°",sim,"/",n_sims))
    lapply(deconv_meth_sup, function(meth) {
      if (!file.exists(paste0("4c_deconv/deconv/sup/",
                              date,
                              "_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))) {
        print(paste0("Running method ",meth))
        data_jdr = adata[[lot]][[sim]][,c(paste0("MultiMAP",1:10))]
        type_idx = adata[[lot]][[sim]][,"Type"]
        sample_idx = adata[[lot]][[sim]][,"Source"]
        Trna = t(data_jdr[type_idx=="ref" & sample_idx=="RNA",])
        Tmet = t(data_jdr[type_idx=="ref" & sample_idx=="MET",])
        Drna = t(data_jdr[type_idx=="sample" & sample_idx=="RNA",])
        Dmet = t(data_jdr[type_idx=="sample" & sample_idx=="MET",])
        Tconcat = rbind(Trna,Tmet)
        Dconcat = rbind(Drna,Dmet)
        rownames(Tconcat) = paste0(rownames(Tconcat),
                                   rep(c("rna","met"),each=nrow(Trna)))
        rownames(Dconcat) = paste0(rownames(Dconcat),
                                   rep(c("rna","met"),each=nrow(Drna)))
        deconv_rna = run_jdr_deconvolution_sup(meth, Drna, Trna)
        deconv_met = run_jdr_deconvolution_sup(meth, Dmet, Tmet)
        deconv_concat = run_jdr_deconvolution_sup(meth, Dconcat, Tconcat)
        deconv_res = list("rna"=deconv_rna,
                          "met"=deconv_met,
                          "concat"=deconv_concat)
        saveRDS(deconv_res,paste0("4c_deconv/deconv/sup/",
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
tmp=lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    print(paste0("Running simulation n°",sim,"/",n_sims))
    lapply(deconv_meth_unsup, function(meth) {
      if (!file.exists(paste0("4c_deconv/deconv/unsup/",
                              date,
                              "_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))) {
        print(paste0("Running method ",meth))
        data_jdr = adata[[lot]][[sim]][,c(paste0("MultiMAP",1:10))]
        type_idx = adata[[lot]][[sim]][,"Type"]
        sample_idx = adata[[lot]][[sim]][,"Source"]
        Drna = t(data_jdr[type_idx=="sample" & sample_idx=="RNA",])
        Dmet = t(data_jdr[type_idx=="sample" & sample_idx=="MET",])
        Dconcat = rbind(Drna,Dmet)
        rownames(Dconcat) = paste0(rownames(Dconcat),
                                   rep(c("rna","met"),each=nrow(Drna)))
        deconv_rna = run_jdr_deconvolution_unsup(meth, Drna, Atrue[[lot]][[sim]])
        deconv_met = run_jdr_deconvolution_unsup(meth, Dmet, Atrue[[lot]][[sim]])
        deconv_concat = run_jdr_deconvolution_unsup(meth, Dconcat, Atrue[[lot]][[sim]])
        deconv_res = list("rna"=deconv_rna,
                          "met"=deconv_met,
                          "concat"=deconv_concat)
        saveRDS(deconv_res,paste0("4c_deconv/deconv/unsup/",
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
## Score deconvolution
## ----
score_perf_rna_sup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth_sup, function(meth) {
      if (file.exists(paste0("4c_deconv/deconv/sup/",
                              date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("4c_deconv/deconv/sup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$rna$res
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
    names(res)=deconv_meth_sup
    return(res)
  })
})
score_perf_met_sup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth_sup, function(meth) {
      if (file.exists(paste0("4c_deconv/deconv/sup/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("4c_deconv/deconv/sup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$met$res
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
    names(res)=deconv_meth_sup
    return(res)
  })
})
score_perf_concat_sup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth_sup, function(meth) {
      if (file.exists(paste0("4c_deconv/deconv/sup/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("4c_deconv/deconv/sup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$concat$res
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
    names(res)=deconv_meth_sup
    return(res)
  })
})
score_perf_rna_unsup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth_unsup, function(meth) {
      if (file.exists(paste0("4c_deconv/deconv/unsup/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("4c_deconv/deconv/unsup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$rna$res
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
    names(res)=deconv_meth_unsup
    return(res)
  })
})
score_perf_met_unsup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth_unsup, function(meth) {
      if (file.exists(paste0("4c_deconv/deconv/unsup/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("4c_deconv/deconv/unsup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$met$res
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
    names(res)=deconv_meth_unsup
    return(res)
  })
})
score_perf_concat_unsup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth_unsup, function(meth) {
      if (file.exists(paste0("4c_deconv/deconv/unsup/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("4c_deconv/deconv/unsup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$concat$res
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
    names(res)=deconv_meth_unsup
    return(res)
  })
})

df_score_perf_rna_sup = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    res1 = lapply(seq(length(deconv_meth_sup)), function(deconv) {
      df = do.call(rbind,score_perf_rna_sup[[lot]][[sim]][[deconv]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(colnames(df), each=nrow(df)),
                      "setting"=rep(rownames(df), ncol(df)))
      df$deconv = deconv_meth_sup[deconv]
      return(df)
    })
    res1 = do.call(rbind,res1)
    res1$sim = sim
    res1$type = 'rna'
    return(res1)
  })
  res2 = do.call(rbind,res2)
  res2$dataset = name_data[lot]
  return(res2)
})
df_score_perf_met_sup = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    res1 = lapply(seq(length(deconv_meth_sup)), function(deconv) {
      df = do.call(rbind,score_perf_met_sup[[lot]][[sim]][[deconv]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(colnames(df), each=nrow(df)),
                      "setting"=rep(rownames(df), ncol(df)))
      df$deconv = deconv_meth_sup[deconv]
      return(df)
    })
    res1 = do.call(rbind,res1)
    res1$sim = sim
    res1$type = 'met'
    return(res1)
  })
  res2 = do.call(rbind,res2)
  res2$dataset = name_data[lot]
  return(res2)
})
df_score_perf_concat_sup = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    res1 = lapply(seq(length(deconv_meth_sup)), function(deconv) {
      df = do.call(rbind,score_perf_concat_sup[[lot]][[sim]][[deconv]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(colnames(df), each=nrow(df)),
                      "setting"=rep(rownames(df), ncol(df)))
      df$deconv = deconv_meth_sup[deconv]
      return(df)
    })
    res1 = do.call(rbind,res1)
    res1$sim = sim
    res1$type = 'concat'
    return(res1)
  })
  res2 = do.call(rbind,res2)
  res2$dataset = name_data[lot]
  return(res2)
})
df_score_perf_rna_unsup = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    res1 = lapply(seq(length(deconv_meth_unsup)), function(deconv) {
      df = do.call(rbind,score_perf_rna_unsup[[lot]][[sim]][[deconv]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(colnames(df), each=nrow(df)),
                      "setting"=rep(rownames(df), ncol(df)))
      df$deconv = deconv_meth_unsup[deconv]
      return(df)
    })
    res1 = do.call(rbind,res1)
    res1$sim = sim
    res1$type = 'rna'
    return(res1)
  })
  res2 = do.call(rbind,res2)
  res2$dataset = name_data[lot]
  return(res2)
})
df_score_perf_met_unsup = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    res1 = lapply(seq(length(deconv_meth_unsup)), function(deconv) {
      df = do.call(rbind,score_perf_met_unsup[[lot]][[sim]][[deconv]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(colnames(df), each=nrow(df)),
                      "setting"=rep(rownames(df), ncol(df)))
      df$deconv = deconv_meth_unsup[deconv]
      return(df)
    })
    res1 = do.call(rbind,res1)
    res1$sim = sim
    res1$type = 'met'
    return(res1)
  })
  res2 = do.call(rbind,res2)
  res2$dataset = name_data[lot]
  return(res2)
})
df_score_perf_concat_unsup = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    res1 = lapply(seq(length(deconv_meth_unsup)), function(deconv) {
      df = do.call(rbind,score_perf_concat_unsup[[lot]][[sim]][[deconv]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(colnames(df), each=nrow(df)),
                      "setting"=rep(rownames(df), ncol(df)))
      df$deconv = deconv_meth_unsup[deconv]
      return(df)
    })
    res1 = do.call(rbind,res1)
    res1$sim = sim
    res1$type = 'concat'
    return(res1)
  })
  res2 = do.call(rbind,res2)
  res2$dataset = name_data[lot]
  return(res2)
})

df_score_perf_rna_sup = do.call(rbind,df_score_perf_rna_sup)
df_score_perf_met_sup = do.call(rbind,df_score_perf_met_sup)
df_score_perf_concat_sup = do.call(rbind,df_score_perf_concat_sup)
df_score_perf_rna_unsup = do.call(rbind,df_score_perf_rna_unsup)
df_score_perf_met_unsup = do.call(rbind,df_score_perf_met_unsup)
df_score_perf_concat_unsup = do.call(rbind,df_score_perf_concat_unsup)
df_score_perf = rbind(df_score_perf_rna_sup,df_score_perf_met_sup,df_score_perf_concat_sup,
                      df_score_perf_rna_unsup,df_score_perf_met_unsup,df_score_perf_concat_unsup)
rm(df_score_perf_rna_sup,df_score_perf_met_sup,df_score_perf_concat_sup,
   df_score_perf_rna_unsup,df_score_perf_met_unsup,df_score_perf_concat_unsup,
   score_perf_rna_sup,score_perf_met_sup,score_perf_concat_sup,
   score_perf_rna_unsup,score_perf_met_unsup,score_perf_concat_unsup)
saveRDS(df_score_perf, paste0("4c_deconv/perf_scores/",
                              date,
                              "_scores.rds"))

## ----
## Time deconvolution
## ----
time_perf_rna_sup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth_sup, function(meth) {
      if (file.exists(paste0("4c_deconv/deconv/sup/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        return(readRDS(paste0("4c_deconv/deconv/sup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$rna$time_elapsed)
      }
      else {print(paste0(meth," not run yet"))}
    })
    names(res)=deconv_meth_sup
    return(res)
  })
})
time_perf_met_sup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth_sup, function(meth) {
      if (file.exists(paste0("4c_deconv/deconv/sup/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        return(readRDS(paste0("4c_deconv/deconv/sup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$met$time_elapsed)
      }
      else {print(paste0(meth," not run yet"))}
    })
    names(res)=deconv_meth_sup
    return(res)
  })
})
time_perf_concat_sup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth_sup, function(meth) {
      if (file.exists(paste0("4c_deconv/deconv/sup/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        return(readRDS(paste0("4c_deconv/deconv/sup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$concat$time_elapsed)
      }
      else {print(paste0(meth," not run yet"))}
    })
    names(res)=deconv_meth_sup
    return(res)
  })
})
time_perf_rna_unsup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth_unsup, function(meth) {
      if (file.exists(paste0("4c_deconv/deconv/unsup/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        return(readRDS(paste0("4c_deconv/deconv/unsup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$rna$time_elapsed)
      }
      else {print(paste0(meth," not run yet"))}
    })
    names(res)=deconv_meth_unsup
    return(res)
  })
})
time_perf_met_unsup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth_unsup, function(meth) {
      if (file.exists(paste0("4c_deconv/deconv/unsup/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        return(readRDS(paste0("4c_deconv/deconv/unsup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$met$time_elapsed)
      }
      else {print(paste0(meth," not run yet"))}
    })
    names(res)=deconv_meth_unsup
    return(res)
  })
})
time_perf_concat_unsup = lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res = lapply(deconv_meth_unsup, function(meth) {
      if (file.exists(paste0("4c_deconv/deconv/unsup/",
                             date,"_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        return(readRDS(paste0("4c_deconv/deconv/unsup/",
                              date,"_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))$concat$time_elapsed)
      }
      else {print(paste0(meth," not run yet"))}
    })
    names(res)=deconv_meth_unsup
    return(res)
  })
})

df_time_perf_rna_sup = lapply(seq(n_lot), function (lot) {
  res = lapply(seq(n_sims), function(sim) {
    tmp = unlist(time_perf_rna_sup[[lot]][[sim]])
    df = data.frame("values"=tmp,
                    "score"="time",
                    "deconv"=names(time_perf_rna_sup[[lot]][[sim]]),
                    "sim"=sim,
                    "type"='rna',
                    "dataset"=name_data[lot])
    return(df)
  })
  do.call(rbind,res)
})
df_time_perf_met_sup = lapply(seq(n_lot), function (lot) {
  res = lapply(seq(n_sims), function(sim) {
    tmp = unlist(time_perf_met_sup[[lot]][[sim]])
    df = data.frame("values"=tmp,
                    "score"="time",
                    "deconv"=names(time_perf_met_sup[[lot]][[sim]]),
                    "sim"=sim,
                    "type"='met',
                    "dataset"=name_data[lot])
    return(df)
  })
  do.call(rbind,res)
})
df_time_perf_concat_sup = lapply(seq(n_lot), function (lot) {
  res = lapply(seq(n_sims), function(sim) {
    tmp = unlist(time_perf_concat_sup[[lot]][[sim]])
    df = data.frame("values"=tmp,
                    "score"="time",
                    "deconv"=names(time_perf_concat_sup[[lot]][[sim]]),
                    "sim"=sim,
                    "type"='concat',
                    "dataset"=name_data[lot])
    return(df)
  })
  do.call(rbind,res)
})
df_time_perf_rna_unsup = lapply(seq(n_lot), function (lot) {
  res = lapply(seq(n_sims), function(sim) {
    tmp = unlist(time_perf_rna_unsup[[lot]][[sim]])
    df = data.frame("values"=tmp,
                    "score"="time",
                    "deconv"=names(time_perf_rna_unsup[[lot]][[sim]]),
                    "sim"=sim,
                    "type"='rna',
                    "dataset"=name_data[lot])
    return(df)
  })
  do.call(rbind,res)
})
df_time_perf_met_unsup = lapply(seq(n_lot), function (lot) {
  res = lapply(seq(n_sims), function(sim) {
    tmp = unlist(time_perf_met_unsup[[lot]][[sim]])
    df = data.frame("values"=tmp,
                    "score"="time",
                    "deconv"=names(time_perf_met_unsup[[lot]][[sim]]),
                    "sim"=sim,
                    "type"='met',
                    "dataset"=name_data[lot])
    return(df)
  })
  do.call(rbind,res)
})
df_time_perf_concat_unsup = lapply(seq(n_lot), function (lot) {
  res = lapply(seq(n_sims), function(sim) {
    tmp = unlist(time_perf_met_unsup[[lot]][[sim]])
    df = data.frame("values"=tmp,
                    "score"="time",
                    "deconv"=names(time_perf_concat_unsup[[lot]][[sim]]),
                    "sim"=sim,
                    "type"='concat',
                    "dataset"=name_data[lot])
    return(df)
  })
  do.call(rbind,res)
})

df_time_perf_rna_sup = do.call(rbind,df_time_perf_rna_sup)
df_time_perf_met_sup = do.call(rbind,df_time_perf_met_sup)
df_time_perf_concat_sup = do.call(rbind,df_time_perf_concat_sup)
df_time_perf_rna_unsup = do.call(rbind,df_time_perf_rna_unsup)
df_time_perf_met_unsup = do.call(rbind,df_time_perf_met_unsup)
df_time_perf_concat_unsup = do.call(rbind,df_time_perf_concat_unsup)
df_time_perf = rbind(df_time_perf_rna_sup,df_time_perf_met_sup,df_time_perf_concat_sup,
                     df_time_perf_rna_unsup,df_time_perf_met_unsup,df_time_perf_concat_unsup)
rm(df_time_perf_rna_sup,df_time_perf_met_sup,df_time_perf_concat_sup,
   df_time_perf_rna_unsup,df_time_perf_met_unsup,df_time_perf_concat_unsup,
   time_perf_rna_sup,time_perf_met_sup,time_perf_concat_sup,
   time_perf_rna_unsup,time_perf_met_unsup,time_perf_concat_unsup)
saveRDS(df_time_perf, paste0("4c_deconv/perf_scores/",
                              date,
                              "_time.rds"))
