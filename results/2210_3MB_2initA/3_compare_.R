## ----
## Set parameters
## ----
block = ifelse(rev(strsplit(init_type,"")[[1]])[1] == "R",'rna','met')
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
date = strsplit(input_path_T[1],"_")[[1]][1]

input_path_Atrue = sort(list.files(input_path, pattern = "sim"))
Atrue = lapply(input_path_Atrue, function(x) readRDS(paste0(input_path,x))$A_ref)
names(Atrue) <- input_path_Atrue
Atrue=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) Atrue[[10*(lot-1)+sim]])
})
n_celltypes = sapply(Atrue, function(x) nrow(x[[1]]))
n_sample = ncol(Atrue[[1]][[1]])
name_data = unname(sapply(input_path_T,function(x) rev(strsplit(strsplit(x,"_T_rna")[[1]][1],"_")[[1]])[1]))
rm(input_path,input_path_nsims,input_path_Atrue)

input_path_method = paste0("../2210_3MB_2initA/2a_init_",init_type,"/")
input_path_method = sort(list.files(input_path_method))
deconv_methods_in = unique(sapply(input_path_method, function(x)
  strsplit(strsplit(x,"Apred_")[[1]][2],"_")[[1]][1]))
deconv_methods_in[deconv_methods_in=="elastic"] <- "elastic_net"
deconv_methods_out = unique(sapply(input_path_method, function(x)
  strsplit(strsplit(x,"to_")[[1]][2],"_")[[1]][1]))
cat_meth = ifelse(rep(block=="rna",length(deconv_methods_out)),c("DECONica","DECONnmf"),deconv_methods_out)
names(cat_meth) = deconv_methods_out
rm(input_path_method)

Apred_to=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_in, function(meth_in) {
      lapply(deconv_methods_out, function(meth_out) {
        if (file.exists(paste0("../2210_3MB_2initA/2a_init_",init_type,"/",
                               date,"_",name_data[lot],
                               "_Apred_",meth_in,
                               "_to_",meth_out,
                               "_",sim_txt,sim,".rds"))) {
          readRDS(paste0("../2210_3MB_2initA/2a_init_",init_type,"/",
                                date,"_",name_data[lot],
                                "_Apred_",meth_in,
                                "_to_",meth_out,
                                "_",sim_txt,sim,".rds"))$res
        }
        else {print(paste0(lot,sim,meth_in,meth_out," not run yet"))}
      })
    })
  })
})
Apred_bl=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_in, function(meth_in) {
      lapply(deconv_methods_out, function(meth_out) {
        if (file.exists(paste0("../2210_1SB/deconv/",block,"/unsup/",
                               date,"_",name_data[lot],
                               "_Apred_",
                               cat_meth[meth_out],
                               "_",sim_txt,
                               sim,".rds"))) {
          readRDS(paste0(paste0("../2210_1SB/deconv/",block,"/unsup/",
                                date,"_",name_data[lot],
                                "_Apred_",
                                cat_meth[meth_out],
                                "_",sim_txt,
                                sim,".rds")))
        }
        else {print(paste0(lot,sim,meth_out," not run yet"))}
      })
    })
  })
})

## ----
## Compute performance scores
## ----
score_methods = c('rmse','mae','pearson')

score_perf_before = pbapply::pblapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    lapply(seq_along(deconv_methods_in), function(meth_in) {
      lapply(seq_along(deconv_methods_out), function(meth_out) {
        if (!is.matrix(Apred_bl[[lot]][[sim]][[meth_in]][[meth_out]])) {return(list("perf_g"=rep(NA,length(score_methods)),
                                                                               "perf_c"=matrix(NA,ncol=length(score_methods),nrow=n_celltypes[lot]),
                                                                               "sd_c"=rep(NA,length(score_methods)),
                                                                               "med_c"=rep(NA,length(score_methods)),
                                                                               "perf_s"=matrix(NA,ncol=length(score_methods),nrow=n_sample),
                                                                               "sd_s"=rep(NA,length(score_methods)),
                                                                               "med_s"=rep(NA,length(score_methods))))}
        else {
          score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue[[lot]][[sim]],Apred_bl[[lot]][[sim]][[meth_in]][[meth_out]],x))
          score_perf_res_celltype <- sapply(score_methods, function(x)
            sapply(seq(nrow(Apred_bl[[lot]][[sim]][[meth_in]][[meth_out]])), function(i)
              score_perf(Atrue[[lot]][[sim]][i,], Apred_bl[[lot]][[sim]][[meth_in]][[meth_out]][i,], x)))
          score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
          score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
          score_perf_res_sample <- sapply(score_methods, function(x)
            sapply(seq(ncol(Apred_bl[[lot]][[sim]][[meth_in]][[meth_out]])), function(i)
              score_perf(Atrue[[lot]][[sim]][,i], Apred_bl[[lot]][[sim]][[meth_in]][[meth_out]][,i], x)))
          score_sd_res_sample=apply(score_perf_res_sample, 2, sd)
          score_med_res_sample=apply(score_perf_res_sample, 2, median)
          return(list("perf_g"=score_perf_res,
                  "perf_c"=score_perf_res_celltype,
                  "sd_c"=score_sd_res_celltype,
                  "med_c"=score_med_res_celltype,
                  "perf_s"=score_perf_res_sample,
                  "sd_s"=score_sd_res_sample,
                  "med_s"=score_med_res_sample))}
      })
    })
  })
})
score_perf_before_df = lapply(seq(n_lot), function (lot) {
  res3 = lapply(seq(n_sims), function(sim) {
    res2 = lapply(seq_along(deconv_methods_in), function(meth_in) {
      res1 = lapply(seq_along(deconv_methods_out), function(meth_out) {
        df = do.call(rbind,score_perf_before[[lot]][[sim]][[meth_in]][[meth_out]])
        rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
        df = data.frame("values"=c(df),
                        "score"=rep(score_methods, each=nrow(df)),
                        "setting"=rep(rownames(df), ncol(df)))
        df$deconv_out = deconv_methods_out[meth_out]
        return(df)
      })
      res1 = do.call(rbind,res1)
      res1$deconv_in = deconv_methods_in[meth_in]
      return(res1)
    })
    res2 = do.call(rbind,res2)
    res2$sim = sim
    return(res2)
  })
  res3 = do.call(rbind,res3)
  res3$dataset = name_data[lot]
  return(res3)
})
score_perf_before_df = do.call(rbind,score_perf_before_df)
score_perf_before_df$deconv_in = "no"

score_perf_after = pbapply::pblapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    lapply(seq_along(deconv_methods_in), function(meth_in) {
      lapply(seq_along(deconv_methods_out), function(meth_out) {
        if (!is.matrix(Apred_to[[lot]][[sim]][[meth_in]][[meth_out]])) {return(list("perf_g"=rep(NA,length(score_methods)),
                                                                               "perf_c"=matrix(NA,ncol=length(score_methods),nrow=n_celltypes[lot]),
                                                                               "sd_c"=rep(NA,length(score_methods)),
                                                                               "med_c"=rep(NA,length(score_methods)),
                                                                               "perf_s"=matrix(NA,ncol=length(score_methods),nrow=n_sample),
                                                                               "sd_s"=rep(NA,length(score_methods)),
                                                                               "med_s"=rep(NA,length(score_methods))))}
        else {
          score_perf_res <- sapply(score_methods, function(x) score_perf(Atrue[[lot]][[sim]],Apred_to[[lot]][[sim]][[meth_in]][[meth_out]],x))
          score_perf_res_celltype <- sapply(score_methods, function(x)
            sapply(seq(nrow(Apred_to[[lot]][[sim]][[meth_in]][[meth_out]])), function(i)
              score_perf(Atrue[[lot]][[sim]][i,], Apred_to[[lot]][[sim]][[meth_in]][[meth_out]][i,], x)))
          score_sd_res_celltype=apply(score_perf_res_celltype, 2, sd)
          score_med_res_celltype=apply(score_perf_res_celltype, 2, median)
          score_perf_res_sample <- sapply(score_methods, function(x)
            sapply(seq(ncol(Apred_to[[lot]][[sim]][[meth_in]][[meth_out]])), function(i)
              score_perf(Atrue[[lot]][[sim]][,i], Apred_to[[lot]][[sim]][[meth_in]][[meth_out]][,i], x)))
          score_sd_res_sample=apply(score_perf_res_sample, 2, sd)
          score_med_res_sample=apply(score_perf_res_sample, 2, median)
          return(list("perf_g"=score_perf_res,
                      "perf_c"=score_perf_res_celltype,
                      "sd_c"=score_sd_res_celltype,
                      "med_c"=score_med_res_celltype,
                      "perf_s"=score_perf_res_sample,
                      "sd_s"=score_sd_res_sample,
                      "med_s"=score_med_res_sample))}
      })
    })
  })
})
 score_perf_after_df = lapply(seq(n_lot), function (lot) {
  res3 = lapply(seq(n_sims), function(sim) {
    res2 = lapply(seq_along(deconv_methods_in), function(meth_in) {
      res1 = lapply(seq_along(deconv_methods_out), function(meth_out) {
        df = do.call(rbind,score_perf_after[[lot]][[sim]][[meth_in]][[meth_out]])
        rownames(df) = c("perf_g",paste0("perf_c",seq(n_celltypes[lot])),"sd_c","med_c",paste0("perf_s",seq(n_sample)),"sd_s","med_s")
        df = data.frame("values"=c(df),
                        "score"=rep(score_methods, each=nrow(df)),
                        "setting"=rep(rownames(df), ncol(df)))
        df$deconv_out = deconv_methods_out[meth_out]
        return(df)
      })
      res1 = do.call(rbind,res1)
      res1$deconv_in = deconv_methods_in[meth_in]
      return(res1)
    })
    res2 = do.call(rbind,res2)
    res2$sim = sim
    return(res2)
  })
  res3 = do.call(rbind,res3)
  res3$dataset = name_data[lot]
  return(res3)
})
 score_perf_after_df = do.call(rbind,score_perf_after_df)

rm(score_perf_before, score_perf_after)

## ----
## Organize df + load timings
## ----
timing_to=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_in, function(meth_in) {
      sapply(deconv_methods_out, function(meth_out) {
        if (file.exists(paste0("../2210_3MB_2initA/2a_init_",init_type,"/",
                               date,"_",name_data[lot],
                               "_Apred_",meth_in,
                               "_to_",meth_out,
                               "_",sim_txt,sim,".rds"))) {
          tmp=readRDS(paste0("../2210_3MB_2initA/2a_init_",init_type,"/",
                         date,"_",name_data[lot],
                         "_Apred_",meth_in,
                         "_to_",meth_out,
                         "_",sim_txt,sim,".rds"))$time_elapsed
          tmp0=readRDS(paste0("../2210_3MB_2initA/2a_init_",init_type,"/",
                              date,"_",name_data[lot],
                              "_Apred_",meth_in,
                              "_to_",meth_out,
                              "_",sim_txt,sim,".rds"))$res
          if (is.list(tmp)) {
            tmp=tmp[[4]]
            tmp=as.numeric(strsplit(tmp," ")[[1]][1])
            saveRDS(list(res=tmp0,time_elapsed=tmp),paste0("../2210_3MB_2initA/2a_init_",init_type,"/",
                               date,"_",name_data[lot],
                               "_Apred_",meth_in,
                               "_to_",meth_out,
                               "_",sim_txt,sim,".rds"))
            }
          return(tmp)
        }
        else {print(paste0(lot,sim,meth_in,meth_out," not run yet"))}
      })
    })
  })
})
timing_to_df = do.call(rbind,lapply(seq_along(timing_to), function (param1) {
  do.call(rbind,lapply(seq_along(timing_to[[param1]]), function(param2) {
    do.call(rbind,lapply(seq_along(timing_to[[param1]][[param2]]), function(param3) {
      df = data.frame("values"=timing_to[[param1]][[param2]][[param3]],
                      "score"="time",
                      "deconv"=names(timing_to[[param1]][[param2]][[param3]]),
                      "sim"=param2,
                      "dataset"=name_data[param1],
                      "input"=deconv_methods_in[param3])
    }))
  }))
}))

timing_bl=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    sapply(deconv_methods_out, function(meth_out) {
      if (file.exists(paste0("../2210_1SB/timing/",block,"/unsup/",
                               date,"_",name_data[lot],
                               "_timing_",
                               cat_meth[meth_out],
                               "_",sim_txt,
                               sim,".rds"))) {
          as.numeric(readRDS(paste0(paste0("../2210_1SB/timing/",block,"/unsup/",
                                date,"_",name_data[lot],
                                "_timing_",
                                cat_meth[meth_out],
                                "_",sim_txt,
                                sim,".rds"))))
        }
        else {print(paste0(lot,sim,meth_out," not run yet"))}
    })
  })
})
timing_bl_df = do.call(rbind,lapply(seq_along(timing_bl), function (param1) {
  do.call(rbind,lapply(seq_along(timing_bl[[param1]]), function(param2) {
      df = data.frame("values"=timing_bl[[param1]][[param2]],
                      "score"="time",
                      "deconv"=names(timing_bl[[param1]][[param2]]),
                      "sim"=param2,
                      "dataset"=name_data[param1],
                      "input"='no')
  }))
}))

colnames(score_perf_before_df)[colnames(score_perf_before_df)=="score"] <- "scor"
score_perf_before_df <- score_perf_before_df %>% mutate(score=paste(scor,setting), input=deconv_in) %>% select(values,score,input,deconv_out,sim,dataset)
colnames(score_perf_after_df)[colnames(score_perf_after_df)=="score"] <- "scor"
score_perf_after_df <- score_perf_after_df %>% mutate(score=paste(scor,setting), input=deconv_in) %>% select(colnames(score_perf_before_df))

score_perf_before_df <- score_perf_before_df %>%
  mutate(candidate=paste(input,deconv_out)) %>%
  select(values,score,sim,dataset,candidate)
score_perf_after_df <- score_perf_after_df %>%
  mutate(candidate=paste(input,deconv_out)) %>%
  select(values,score,sim,dataset,candidate)

score_tot1 <- bind_rows(score_perf_before_df, score_perf_after_df)
rm(score_perf_before_df, score_perf_after_df)

score_perf_before_df2 <- timing_bl_df %>%
  mutate(candidate=paste(input,deconv)) %>%
  select(values,score,sim,dataset,candidate)
score_perf_after_df2 <- timing_to_df %>%
  mutate(candidate=paste(input,deconv)) %>%
  select(values,score,sim,dataset,candidate)

score_tot2 <- bind_rows(score_perf_before_df2, score_perf_after_df2)
rm(score_perf_before_df2, score_perf_after_df2)

rm(Apred_bl,Apred_to,Atrue,
   timing_bl,timing_bl_df,timing_to,timing_to_df,
   block,cat_meth,date,deconv_methods_in,deconv_methods_out,
   input_path_T,name_data,score_methods,
   n_celltypes,n_lot,n_sample,n_sims,
   convert_to_cat,do_df_deconv,geomMean,harmMean,weighMean,MinMaxNorm,mae,pearson,rmse,score_perf,plotting,
   ranking_avgrank,ranking_concat_per_method,ranking12,ranking12_notime,ranking1234,ranking12345,ranking12345_notime,ranking1245,ranking12453)
