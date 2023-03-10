## ----
## Set parameters
## ----
if (n_meth==3) {deconv_methods_rna = c("RLR","WISP","DeconRNASeq")} else {deconv_methods_rna = c("RLR")}
n_meth = length(deconv_methods_rna)
source("../../src/score_functions.R")
source("../../src/ranking_ranking_procedure_functions.R")

## ----
## load Apred and Atrue
## ----
input_path = "../2210_0simu/simulations/rna/"
input_path_T = sort(list.files(input_path, pattern = "T_rna_ref"))
n_lot = length(input_path_T)
input_path_nsims = sort(list.files(input_path, pattern = "sim"))
n_sims = length(input_path_nsims)/n_lot
input_path_Atrue = sort(list.files(input_path, pattern = "sim"))
Atrue = lapply(input_path_Atrue, function(x) readRDS(paste0(input_path,x))$A_ref)
names(Atrue) <- input_path_Atrue
name_data = unname(sapply(input_path_T,function(x) rev(strsplit(strsplit(x,"_T_rna")[[1]][1],"_")[[1]])[1]))
rm(input_path,input_path_nsims,input_path_Atrue)

input_path_methodssup2 = "../2210_1SB/deconv/met/sup/"
input_path_methodssup2 = sort(list.files(input_path_methodssup2))
deconv_methods_met = unique(sapply(input_path_methodssup2, function(x)
  rev(strsplit(x,"_")[[1]])[2]))
if (n_meth==1) deconv_methods_met = "rlr"
rm(input_path_methodssup2)

Apred_rna=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_rna, function(meth) {
      if (file.exists(paste0("../2210_1SB/deconv/rna/sup/",
                             strsplit(input_path_T[lot],"_T_rna_ref")[[1]][1],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("../2210_1SB/deconv/rna/sup/",
                              strsplit(input_path_T[lot],"_T_rna_ref")[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
      }
      else {print(paste0(meth," not run yet"))}
    })
  })
})
Apred_met=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_met, function(meth) {
      if (file.exists(paste0("../2210_1SB/deconv/met/sup/",
                             strsplit(input_path_T[lot],"_T_rna_ref")[[1]][1],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("../2210_1SB/deconv/met/sup/",
                              strsplit(input_path_T[lot],"_T_rna_ref")[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
      }
      else {print(paste0(meth," not run yet"))}
    })
  })
})

Atrue=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) Atrue[[10*(lot-1)+sim]])
})
n_celltypes = sapply(Atrue, function(x) nrow(x[[1]]))
n_sample = ncol(Atrue[[1]][[1]])

## ----
## Compute mean(Apred_rna,Apred_met)
## ----
Apred_rnaM = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) 
    Arna = Reduce("+",Apred_rna[[lot]][[sim]])/length(Apred_rna[[lot]][[sim]]))
})
Apred_metM = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim)  
    Amet = Reduce("+",Apred_met[[lot]][[sim]])/length(Apred_met[[lot]][[sim]]))
})
Apred = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {    
    Atmp = list(Apred_rnaM[[lot]][[sim]],Apred_metM[[lot]][[sim]])
    Reduce("+",Atmp)/length(Atmp)})
})

## ----
## Compute performance scores
## ----
score_methods = c('rmse','mae','pearson')

score_perf_r = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    lapply(seq(n_meth), function(meth) {
      score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue[[lot]][[sim]],Apred_rna[[lot]][[sim]][[meth]],x))
      score_perf_res_celltype <- sapply(score_methods, function(x)
        sapply(seq(nrow(Apred_rna[[lot]][[sim]][[meth]])), function(i)
          score_perf(Atrue[[lot]][[sim]][i,], Apred_rna[[lot]][[sim]][[meth]][i,], x)))
      score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
      score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
      score_perf_res_sample <- sapply(score_methods, function(x)
        sapply(seq(ncol(Apred_rna[[lot]][[sim]][[meth]])), function(i)
          score_perf(Atrue[[lot]][[sim]][,i], Apred_rna[[lot]][[sim]][[meth]][,i], x)))
      score_sd_res_sample=apply(score_perf_res_sample, 2, sd)
      score_med_res_sample=apply(score_perf_res_sample, 2, median)
      return(list("perf_g"=score_perf_res,
                  "perf_c"=score_perf_res_celltype,
                  "sd_c"=score_sd_res_celltype,
                  "med_c"=score_med_res_celltype,
                  "perf_s"=score_perf_res_sample,
                  "sd_s"=score_sd_res_sample,
                  "med_s"=score_med_res_sample))
    })
    })
  })
score_perf_r_df = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    res1 = lapply(seq(n_meth), function(meth) {
      df = do.call(rbind,score_perf_r[[lot]][[sim]][[meth]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(score_methods, each=nrow(df)),
                      "setting"=rep(rownames(df), ncol(df)))
      df$deconv = deconv_methods_rna[meth]
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
score_perf_r_df = do.call(rbind,score_perf_r_df)

score_perf_m = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    lapply(seq(n_meth), function(meth) {
      score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue[[lot]][[sim]],Apred_met[[lot]][[sim]][[meth]],x))
      score_perf_res_celltype <- sapply(score_methods, function(x)
        sapply(seq(nrow(Apred_met[[lot]][[sim]][[meth]])), function(i)
          score_perf(Atrue[[lot]][[sim]][i,], Apred_met[[lot]][[sim]][[meth]][i,], x)))
      score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
      score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
      score_perf_res_sample <- sapply(score_methods, function(x)
        sapply(seq(ncol(Apred_met[[lot]][[sim]][[meth]])), function(i)
          score_perf(Atrue[[lot]][[sim]][,i], Apred_met[[lot]][[sim]][[meth]][,i], x)))
      score_sd_res_sample=apply(score_perf_res_sample, 2, sd)
      score_med_res_sample=apply(score_perf_res_sample, 2, median)
      return(list("perf_g"=score_perf_res,
                  "perf_c"=score_perf_res_celltype,
                  "sd_c"=score_sd_res_celltype,
                  "med_c"=score_med_res_celltype,
                  "perf_s"=score_perf_res_sample,
                  "sd_s"=score_sd_res_sample,
                  "med_s"=score_med_res_sample))
    })
  })
})
score_perf_m_df = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    res1 = lapply(seq(n_meth), function(meth) {
      df = do.call(rbind,score_perf_m[[lot]][[sim]][[meth]])
      rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
      df = data.frame("values"=c(df),
                      "score"=rep(score_methods, each=nrow(df)),
                      "setting"=rep(rownames(df), ncol(df)))
      df$deconv = deconv_methods_met[meth]
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
score_perf_m_df = do.call(rbind,score_perf_m_df)

score_perf_r_m = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue[[lot]][[sim]],Apred_rnaM[[lot]][[sim]],x))
    score_perf_res_celltype <- sapply(score_methods, function(x)
      sapply(seq(nrow(Apred_rnaM[[lot]][[sim]])), function(i)
        score_perf(Atrue[[lot]][[sim]][i,], Apred_rnaM[[lot]][[sim]][i,], x)))
    score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
    score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
    score_perf_res_sample <- sapply(score_methods, function(x)
      sapply(seq(ncol(Apred_rnaM[[lot]][[sim]])), function(i)
        score_perf(Atrue[[lot]][[sim]][,i], Apred_rnaM[[lot]][[sim]][,i], x)))
    score_sd_res_sample=apply(score_perf_res_sample, 2, sd)
    score_med_res_sample=apply(score_perf_res_sample, 2, median)
    return(list("perf_g"=score_perf_res,
                "perf_c"=score_perf_res_celltype,
                "sd_c"=score_sd_res_celltype,
                "med_c"=score_med_res_celltype,
                "perf_s"=score_perf_res_sample,
                "sd_s"=score_sd_res_sample,
                "med_s"=score_med_res_sample))
  })
})
score_perf_r_m_df = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    df = do.call(rbind,score_perf_r_m[[lot]][[sim]])
    rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
    df = data.frame("values"=c(df),
                    "score"=rep(score_methods, each=nrow(df)),
                    "setting"=rep(rownames(df), ncol(df)))
    df$deconv = "RM"
    df$sim = sim
    return(df)
  })
  res2 = do.call(rbind,res2)
  res2$dataset = name_data[lot]
  return(res2)
})
score_perf_r_m_df = do.call(rbind,score_perf_r_m_df)

score_perf_m_m = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue[[lot]][[sim]],Apred_metM[[lot]][[sim]],x))
    score_perf_res_celltype <- sapply(score_methods, function(x)
      sapply(seq(nrow(Apred_metM[[lot]][[sim]])), function(i)
        score_perf(Atrue[[lot]][[sim]][i,], Apred_metM[[lot]][[sim]][i,], x)))
    score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
    score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
    score_perf_res_sample <- sapply(score_methods, function(x)
      sapply(seq(ncol(Apred_metM[[lot]][[sim]])), function(i)
        score_perf(Atrue[[lot]][[sim]][,i], Apred_metM[[lot]][[sim]][,i], x)))
    score_sd_res_sample=apply(score_perf_res_sample, 2, sd)
    score_med_res_sample=apply(score_perf_res_sample, 2, median)
    return(list("perf_g"=score_perf_res,
                "perf_c"=score_perf_res_celltype,
                "sd_c"=score_sd_res_celltype,
                "med_c"=score_med_res_celltype,
                "perf_s"=score_perf_res_sample,
                "sd_s"=score_sd_res_sample,
                "med_s"=score_med_res_sample))
  })
})
score_perf_m_m_df = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    df = do.call(rbind,score_perf_m_m[[lot]][[sim]])
    rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
    df = data.frame("values"=c(df),
                    "score"=rep(score_methods, each=nrow(df)),
                    "setting"=rep(rownames(df), ncol(df)))
    df$deconv = "MM"
    df$sim = sim
    return(df)
  })
  res2 = do.call(rbind,res2)
  res2$dataset = name_data[lot]
  return(res2)
})
score_perf_m_m_df = do.call(rbind,score_perf_m_m_df)

score_perf_mb = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue[[lot]][[sim]],Apred[[lot]][[sim]],x))
    score_perf_res_celltype <- sapply(score_methods, function(x)
      sapply(seq(nrow(Apred[[lot]][[sim]])), function(i)
        score_perf(Atrue[[lot]][[sim]][i,], Apred[[lot]][[sim]][i,], x)))
    score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
    score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
    score_perf_res_sample <- sapply(score_methods, function(x)
      sapply(seq(ncol(Apred[[lot]][[sim]])), function(i)
        score_perf(Atrue[[lot]][[sim]][,i], Apred[[lot]][[sim]][,i], x)))
    score_sd_res_sample=apply(score_perf_res_sample, 2, sd)
    score_med_res_sample=apply(score_perf_res_sample, 2, median)
    return(list("perf_g"=score_perf_res,
                "perf_c"=score_perf_res_celltype,
                "sd_c"=score_sd_res_celltype,
                "med_c"=score_med_res_celltype,
                "perf_s"=score_perf_res_sample,
                "sd_s"=score_sd_res_sample,
                "med_s"=score_med_res_sample))
  })
})
score_perf_mb_df = lapply(seq(n_lot), function (lot) {
  res2 = lapply(seq(n_sims), function(sim) {
    df = do.call(rbind,score_perf_mb[[lot]][[sim]])
    rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
    df = data.frame("values"=c(df),
                    "score"=rep(score_methods, each=nrow(df)),
                    "setting"=rep(rownames(df), ncol(df)))
    df$deconv = "MB"
    df$sim = sim
    return(df)
  })
  res2 = do.call(rbind,res2)
  res2$dataset = name_data[lot]
  return(res2)
})
score_perf_mb_df = do.call(rbind,score_perf_mb_df)

rm(score_perf_r, score_perf_m, score_perf_r_m, score_perf_m_m, score_perf_mb)
rm(Apred, Apred_met, Apred_rna, Apred_metM, Apred_rnaM, Atrue)

## ----
## Organize df
## ----
colnames(score_perf_r_df)[colnames(score_perf_r_df)=="score"] <- "scor"
score_perf_r_df <- score_perf_r_df %>% mutate(score=paste(scor,setting)) %>% select(values,score,deconv,sim,dataset)
colnames(score_perf_m_df)[colnames(score_perf_m_df)=="score"] <- "scor"
score_perf_m_df <- score_perf_m_df %>% mutate(score=paste(scor,setting)) %>% select(colnames(score_perf_r_df))
colnames(score_perf_r_m_df)[colnames(score_perf_r_m_df)=="score"] <- "scor"
score_perf_r_m_df <- score_perf_r_m_df %>% mutate(score=paste(scor,setting)) %>% select(colnames(score_perf_r_df))
colnames(score_perf_m_m_df)[colnames(score_perf_m_m_df)=="score"] <- "scor"
score_perf_m_m_df <- score_perf_m_m_df %>% mutate(score=paste(scor,setting)) %>% select(colnames(score_perf_r_df))
colnames(score_perf_mb_df)[colnames(score_perf_mb_df)=="score"] <- "scor"
score_perf_mb_df <- score_perf_mb_df %>% mutate(score=paste(scor,setting)) %>% select(colnames(score_perf_r_df))

candidate_cat = list("MB"="MB",
                     "SB"=c(deconv_methods_met,deconv_methods_rna,"RM","MM"))
score_perf_r_df <- score_perf_r_df %>%
  mutate(candidate=paste(deconv,convert_to_cat(candidate_cat,deconv))) %>%
  select(values,score,sim,dataset,candidate)
score_perf_m_df <- score_perf_m_df %>%
  mutate(candidate=paste(deconv,convert_to_cat(candidate_cat,deconv))) %>%
  select(values,score,sim,dataset,candidate)
score_perf_r_m_df <- score_perf_r_m_df %>%
  mutate(candidate=paste(deconv,convert_to_cat(candidate_cat,deconv))) %>%
  select(values,score,sim,dataset,candidate)
score_perf_m_m_df <- score_perf_m_m_df %>%
  mutate(candidate=paste(deconv,convert_to_cat(candidate_cat,deconv))) %>%
  select(values,score,sim,dataset,candidate)
score_perf_mb_df <- score_perf_mb_df %>%
  mutate(candidate=paste(deconv,convert_to_cat(candidate_cat,deconv))) %>%
  select(values,score,sim,dataset,candidate)

score_tot <- bind_rows(score_perf_r_df, score_perf_m_df, score_perf_r_m_df, score_perf_m_m_df, score_perf_mb_df)
rm(score_perf_r_df, score_perf_m_df, score_perf_r_m_df, score_perf_m_m_df, score_perf_mb_df)

rm(candidate_cat,deconv_methods_met,deconv_methods_rna,input_path_T,
   n_celltypes,n_sample,n_sims,n_lot,
   name_data,score_methods,
   convert_to_cat,do_df_deconv,geomMean,harmMean,weighMean,MinMaxNorm,mae,pearson,rmse,score_perf,plotting,
   ranking_avgrank,ranking_concat_per_method,ranking12,ranking1234,ranking12345,ranking1245,ranking12453)
