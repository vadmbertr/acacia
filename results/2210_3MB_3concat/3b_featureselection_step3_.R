## ----
## load scores restricted and selected
## ----
input_path = paste0("../2210_0simu/simulations/rna/")
n_lot = length(sort(list.files(input_path, pattern = "T_rna")))
n_sims = length(sort(list.files(input_path, pattern = "sim")))/n_lot
name_data = unname(sapply(list.files(input_path, pattern = "T_rna"),function(x) rev(strsplit(strsplit(x,"T_rna")[[1]][1],"_")[[1]])[1]))
rm(input_path)

score_bl_met_sup <- readRDS("../2210_1SB/perf_scores/221101_scores_met_sup.rds")
score_bl_met_sup$datatype = "met"
score_bl_met_unsup <- readRDS("../2210_1SB/perf_scores/221101_scores_met_unsup.rds")
score_bl_met_unsup$datatype = "met"
score_bl_rna_sup <- readRDS("../2210_1SB/perf_scores/221101_scores_rna_sup.rds")
score_bl_rna_sup$datatype = "rna"
score_bl_rna_unsup <- readRDS("../2210_1SB/perf_scores/221101_scores_rna_unsup.rds")
score_bl_rna_unsup$datatype = "rna"
score_bl <- bind_rows(score_bl_met_sup,score_bl_met_unsup,
                      score_bl_rna_sup,score_bl_rna_unsup)
rm(score_bl_met_sup,score_bl_met_unsup,score_bl_rna_sup,score_bl_rna_unsup)

score_fs_toast <- readRDS("../2210_3MB_3concat/3b_featureselection_step2/TOAST/perf_scores/221101_scores.rds")
score_fs_toast$datatype = "TOAST"
score_fs_cimpleg <- readRDS("../2210_3MB_3concat/3b_featureselection_step2/CimpleG/perf_scores/221101_scores.rds")
score_fs_cimpleg$datatype = "CimpleG"
score_fs <- bind_rows(score_fs_toast,score_fs_cimpleg)
rm(score_fs_toast,score_fs_cimpleg)

score_restricted_toast <- readRDS("../2210_3MB_3concat/3b_featureselection_step2/TOAST/perf_scores/221101_scores_restricted.rds")
score_restricted_toast$datatype = "TOAST_r"
score_restricted_cimpleg <- readRDS("../2210_3MB_3concat/3b_featureselection_step2/CimpleG/perf_scores/221101_scores_restricted.rds")
score_restricted_cimpleg$datatype = "CimpleG_r"
score_restricted <- bind_rows(score_restricted_toast,score_restricted_cimpleg)
rm(score_restricted_toast,score_restricted_cimpleg)

score_all <- bind_rows(score_bl,score_fs,score_restricted)
deconv_methods = sort(unique(score_fs$deconv))
rm(score_bl,score_fs,score_restricted)

## ----
## Load timings
## ----
time_bl_met_sup <- readRDS("../2210_1SB/perf_scores/221101_time_met_sup.rds")
time_bl_met_sup$datatype = "met"
time_bl_met_unsup <- readRDS("../2210_1SB/perf_scores/221101_time_met_unsup.rds")
time_bl_met_unsup$datatype = "met"
time_bl_rna_sup <- readRDS("../2210_1SB/perf_scores/221101_time_rna_sup.rds")
time_bl_rna_sup$datatype = "rna"
time_bl_rna_sup$values = as.numeric(time_bl_rna_sup$values)
time_bl_rna_unsup <- readRDS("../2210_1SB/perf_scores/221101_time_rna_unsup.rds")
time_bl_rna_unsup$datatype = "rna"
time_bl_rna_unsup$values = as.numeric(time_bl_rna_unsup$values)
time_bl <- bind_rows(time_bl_met_sup,time_bl_met_unsup,
                      time_bl_rna_sup,time_bl_rna_unsup)
rm(time_bl_met_sup,time_bl_met_unsup,time_bl_rna_sup,time_bl_rna_unsup)

time_fs_t=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    sapply(deconv_methods, function(meth) {
      if (file.exists(paste0("../2210_3MB_3concat/3b_featureselection_step2/TOAST/deconv/221101_",name_data[lot],
                               "_Apred_",
                               meth,
                               "_",sim_txt,sim,".rds"))) {
          tmp=readRDS(paste0("../2210_3MB_3concat/3b_featureselection_step2/TOAST/deconv/221101_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,sim,".rds"))$time_elapsed
          tmp0=readRDS(paste0("../2210_3MB_3concat/3b_featureselection_step2/TOAST/deconv/221101_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,sim,".rds"))$res
          if (is.list(tmp)) {
            tmp=tmp[[4]]
            tmp=as.numeric(strsplit(tmp," ")[[1]][1])
            saveRDS(list(res=tmp0,time_elapsed=tmp),paste0("../2210_3MB_3concat/3b_featureselection_step2/TOAST/deconv/221101_",name_data[lot],
                                                           "_Apred_",
                                                           meth,
                                                           "_",sim_txt,sim,".rds"))
          }
          return(tmp)
        }
        else {print(paste0(lot,sim,meth," not run yet"))}
      })
  })
})
time_fs_t_df = do.call(rbind,lapply(seq_along(time_fs_t), function (lot) {
  do.call(rbind,lapply(seq_along(time_fs_t[[lot]]), function(sim) {
    df = data.frame("values"=time_fs_t[[lot]][[sim]],
                    "score"="time",
                    "deconv"=names(time_fs_t[[lot]][[sim]]),
                    "sim"=sim,
                    "dataset"=name_data[lot],
                    "datatype"="TOAST")
  }))
}))

time_fs_c=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    sapply(deconv_methods, function(meth) {
      if (file.exists(paste0("../2210_3MB_3concat/3b_featureselection_step2/CimpleG/deconv/221101_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,sim,".rds"))) {
        tmp=readRDS(paste0("../2210_3MB_3concat/3b_featureselection_step2/CimpleG/deconv/221101_",name_data[lot],
                           "_Apred_",
                           meth,
                           "_",sim_txt,sim,".rds"))$time_elapsed
        tmp0=readRDS(paste0("../2210_3MB_3concat/3b_featureselection_step2/CimpleG/deconv/221101_",name_data[lot],
                            "_Apred_",
                            meth,
                            "_",sim_txt,sim,".rds"))$res
        if (is.list(tmp)) {
          tmp=tmp[[4]]
          tmp=as.numeric(strsplit(tmp," ")[[1]][1])
          saveRDS(list(res=tmp0,time_elapsed=tmp),paste0("../2210_3MB_3concat/3b_featureselection_step2/CimpleG/deconv/221101_",name_data[lot],
                                                         "_Apred_",
                                                         meth,
                                                         "_",sim_txt,sim,".rds"))
        }
        return(tmp)
      }
      else {print(paste0(lot,sim,meth," not run yet"))}
    })
  })
})
time_fs_c_df = do.call(rbind,lapply(seq_along(time_fs_c), function (lot) {
  do.call(rbind,lapply(seq_along(time_fs_c[[lot]]), function(sim) {
    df = data.frame("values"=time_fs_c[[lot]][[sim]],
                    "score"="time",
                    "deconv"=names(time_fs_c[[lot]][[sim]]),
                    "sim"=sim,
                    "dataset"=name_data[lot],
                    "datatype"="CimpleG")
  }))
}))

time_restricted_t=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    sapply(deconv_methods, function(meth) {
      if (file.exists(paste0("../2210_3MB_3concat/3b_featureselection_step2/TOAST/deconv_restricted/221101_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,sim,".rds"))) {
        tmp=readRDS(paste0("../2210_3MB_3concat/3b_featureselection_step2/TOAST/deconv_restricted/221101_",name_data[lot],
                           "_Apred_",
                           meth,
                           "_",sim_txt,sim,".rds"))$time_elapsed
        tmp0=readRDS(paste0("../2210_3MB_3concat/3b_featureselection_step2/TOAST/deconv_restricted/221101_",name_data[lot],
                            "_Apred_",
                            meth,
                            "_",sim_txt,sim,".rds"))$res
        if (is.list(tmp)) {
          tmp=tmp[[4]]
          tmp=as.numeric(strsplit(tmp," ")[[1]][1])
          saveRDS(list(res=tmp0,time_elapsed=tmp),paste0("../2210_3MB_3concat/3b_featureselection_step2/TOAST/deconv_restricted/221101_",name_data[lot],
                                                         "_Apred_",
                                                         meth,
                                                         "_",sim_txt,sim,".rds"))
        }
        return(tmp)
      }
      else {print(paste0(lot,sim,meth," not run yet"))}
    })
  })
})
time_restricted_t_df = do.call(rbind,lapply(seq_along(time_restricted_t), function (lot) {
  do.call(rbind,lapply(seq_along(time_restricted_t[[lot]]), function(sim) {
    df = data.frame("values"=time_restricted_t[[lot]][[sim]],
                    "score"="time",
                    "deconv"=names(time_restricted_t[[lot]][[sim]]),
                    "sim"=sim,
                    "dataset"=name_data[lot],
                    "datatype"="TOAST_r")
  }))
}))

time_restricted_c=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    sapply(deconv_methods, function(meth) {
      if (file.exists(paste0("../2210_3MB_3concat/3b_featureselection_step2/CimpleG/deconv_restricted/221101_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,sim,".rds"))) {
        tmp=readRDS(paste0("../2210_3MB_3concat/3b_featureselection_step2/CimpleG/deconv_restricted/221101_",name_data[lot],
                           "_Apred_",
                           meth,
                           "_",sim_txt,sim,".rds"))$time_elapsed
        tmp0=readRDS(paste0("../2210_3MB_3concat/3b_featureselection_step2/CimpleG/deconv_restricted/221101_",name_data[lot],
                            "_Apred_",
                            meth,
                            "_",sim_txt,sim,".rds"))$res
        if (is.list(tmp)) {
          tmp=tmp[[4]]
          tmp=as.numeric(strsplit(tmp," ")[[1]][1])
          saveRDS(list(res=tmp0,time_elapsed=tmp),paste0("../2210_3MB_3concat/3b_featureselection_step2/CimpleG/deconv_restricted/221101_",name_data[lot],
                                                         "_Apred_",
                                                         meth,
                                                         "_",sim_txt,sim,".rds"))
        }
        return(tmp)
      }
      else {print(paste0(lot,sim,meth," not run yet"))}
    })
  })
})
time_restricted_c_df = do.call(rbind,lapply(seq_along(time_restricted_c), function (lot) {
  do.call(rbind,lapply(seq_along(time_restricted_c[[lot]]), function(sim) {
    df = data.frame("values"=time_restricted_c[[lot]][[sim]],
                    "score"="time",
                    "deconv"=names(time_restricted_c[[lot]][[sim]]),
                    "sim"=sim,
                    "dataset"=name_data[lot],
                    "datatype"="CimpleG_r")
  }))
}))

time_fs_all <- bind_rows(time_restricted_c_df,time_restricted_t_df,
                         time_fs_c_df,time_fs_t_df)
rm(time_fs_t,time_fs_c,time_restricted_c,time_restricted_t,
   time_fs_t_df,time_fs_c_df,time_restricted_c_df,time_restricted_t_df)

## ----
## Organise df
## ----
colnames(score_all)[colnames(score_all)=="score"] <- "scor"
score_all <- score_all %>%
  mutate(score=paste(scor,setting)) %>%
  select(c("values","score","sim","dataset","datatype","deconv")) %>%
  mutate(candidate=paste(datatype,deconv)) %>%
  select(c('values','score','sim','dataset','candidate'))

time_bl$sim <- as.numeric(time_bl$sim) 
time_all <- bind_rows(time_bl,time_fs_all) %>%
  mutate(candidate=paste(datatype,deconv)) %>%
  select(c('values','score','sim','dataset','candidate'))
rm(time_fs_all,time_bl)

score_all <- score_all[c(grep("CIBERSORT",score_all$candidate),
                         grep("rlr",score_all$candidate),
                         grep("ICA",score_all$candidate),
                         grep("CIBERSORT3",score_all$candidate),
                         grep("RLR",score_all$candidate),
                         grep("DECONica",score_all$candidate),
                         grep("cibersort",score_all$candidate)),]
score_all <- score_all[-grep("CIBERSORT4",score_all$candidate),]
time_all <- time_all[c(grep("CIBERSORT",time_all$candidate),
                         grep("rlr",time_all$candidate),
                         grep("ICA",time_all$candidate),
                         grep("CIBERSORT3",time_all$candidate),
                         grep("RLR",time_all$candidate),
                         grep("DECONica",time_all$candidate),
                         grep("cibersort",time_all$candidate)),]
time_all <- time_all[-grep("CIBERSORT4",time_all$candidate),]
score_all$candidate[score_all$candidate == "met CIBERSORT"] = "met cibersort"
score_all$candidate[score_all$candidate == "rna CIBERSORT3"] = "rna cibersort"
score_all$candidate[score_all$candidate == "met ICA"] = "met ica"
score_all$candidate[score_all$candidate == "rna RLR"] = "rna rlr"
score_all$candidate[score_all$candidate == "rna DECONica"] = "rna ica"
time_all$candidate[time_all$candidate == "met CIBERSORT"] = "met cibersort"
time_all$candidate[time_all$candidate == "rna CIBERSORT3"] = "rna cibersort"
time_all$candidate[time_all$candidate == "met ICA"] = "met ica"
time_all$candidate[time_all$candidate == "rna RLR"] = "rna rlr"
time_all$candidate[time_all$candidate == "rna DECONica"] = "rna ica"

rm(deconv_methods,n_lot,n_sims,name_data)
