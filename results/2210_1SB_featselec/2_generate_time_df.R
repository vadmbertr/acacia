## ----
## Set parameters
## ----
library(dplyr)
block="rna" # to vary, either met or rna

## ----
## load timings
## ----
input_path = paste0("timing/TOAST/",block,"/sup/")
input_path_T = sort(list.files(input_path))
date = strsplit(input_path_T[1],"_")[[1]][1]
rm(input_path)

input_path_methodssup = paste0("timing/TOAST/",block,"/sup/")
input_path_methodssup = sort(list.files(input_path_methodssup))
deconv_methods_sup = unique(sapply(input_path_methodssup, function(x)
  rev(strsplit(x,"_")[[1]])[2]))
deconv_methods_sup[deconv_methods_sup=="net"] = "elastic_net"
rm(input_path_methodssup)

input_path_methodsunsup = paste0("timing/TOAST/",block,"/unsup/")
input_path_methodsunsup = sort(list.files(input_path_methodsunsup))
deconv_methods_unsup = unique(sapply(input_path_methodsunsup, function(x)
  rev(strsplit(x,"_")[[1]])[2]))
rm(input_path_methodsunsup)

name_data = unique(unname(sapply(input_path_T,function(x) strsplit(strsplit(x,"_timing")[[1]][1],paste0(date,"_"))[[1]][2])))
n_lot = length(name_data)
n_sims = length(input_path_T)/length(deconv_methods_sup)/n_lot/2

## ----
## Timing ref-based
## ----
tmp=lapply(seq(n_lot), function (lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("timing/TOAST/",block,"/sup/",
                             date,"_",name_data[lot],
                             "_timing_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        time=readRDS(paste0("timing/TOAST/",block,"/sup/",
                            date,"_",name_data[lot],
                            "_timing_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        time2=readRDS(paste0("timing/TOAST/",block,"/sup/",
                            date,"_",name_data[lot],
                            "_timingbl2_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        if (is.list(time)) {
          time=time[[4]]
          time=strsplit(time," sec")[[1]][1]
          time=as.numeric(strsplit(time, " ")[[1]][2])
          saveRDS(time, paste0("timing/TOAST/",block,"/sup/",
                         date,"_",name_data[lot],
                         "_timing_",
                         meth,
                         "_",sim_txt,
                         sim,".rds"))
        }
        if (is.list(time2)) {
          time2=time2[[4]]
          time2=strsplit(time2," sec")[[1]][1]
          time2=as.numeric(strsplit(time2, " ")[[1]][2])
          saveRDS(time2, paste0("timing/TOAST/",block,"/sup/",
                               date,"_",name_data[lot],
                               "_timingbl2_",
                               meth,
                               "_",sim_txt,
                               sim,".rds"))
        }
      }
      if (file.exists(paste0("timing/CimpleG/",block,"/sup/",
                             date,"_",name_data[lot],
                             "_timing_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        time=readRDS(paste0("timing/CimpleG/",block,"/sup/",
                            date,"_",name_data[lot],
                            "_timing_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        time2=readRDS(paste0("timing/CimpleG/",block,"/sup/",
                            date,"_",name_data[lot],
                            "_timingbl2_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        if (is.list(time)) {
          time=time[[4]]
          time=strsplit(time," sec")[[1]][1]
          time=as.numeric(strsplit(time, " ")[[1]][2])
          saveRDS(time, paste0("timing/CimpleG/",block,"/sup/",
                               date,"_",name_data[lot],
                               "_timing_",
                               meth,
                               "_",sim_txt,
                               sim,".rds"))
        }
        if (is.list(time2)) {
          time2=time2[[4]]
          time2=strsplit(time2," sec")[[1]][1]
          time2=as.numeric(strsplit(time2, " ")[[1]][2])
          saveRDS(time2, paste0("timing/CimpleG/",block,"/sup/",
                               date,"_",name_data[lot],
                               "_timingbl2_",
                               meth,
                               "_",sim_txt,
                               sim,".rds"))
        }
      }
    })
  })
})

timing_toast_sup = lapply(seq(n_lot), function (lot) {
  res = sapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res0 = sapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("timing/TOAST/",block,"/sup/",
                              date,"_",name_data[lot],
                              "_timing_",
                              meth,
                             "_",sim_txt,
                              sim,".rds"))) {
        time=readRDS(paste0("timing/TOAST/",block,"/sup/",
                              date,"_",name_data[lot],
                              "_timing_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        time2=readRDS(paste0("timing/TOAST/",block,"/sup/",
                            date,"_",name_data[lot],
                            "_timingbl2_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        vec = c(time,time2)
        names(vec) = c("selec","restricted")
        return(vec)
        }
      else {print(paste0(meth," not run yet"))}
      })
    vec = c(res0)
    names(vec) = c(sapply(colnames(res0), function(x)
      sapply(rownames(res0), function(y) paste(x,y))))
    return(vec)
    })
  res1 = res[grep("selec",rownames(res)),]
  res2 = res[grep("restricted",rownames(res)),]
  rownames(res1) = deconv_methods_sup
  rownames(res2) = deconv_methods_sup
  colnames(res1) = seq(n_sims)
  colnames(res2) = seq(n_sims)
  return(list("selec"=res1,"restricted"=res2))
  })
timing_toast_sup_df = lapply(seq(n_lot), function (lot) {
  df = data.frame("values"=c(timing_toast_sup[[lot]]$selec),
                  "score"="time",
                  "deconv"=rep(rownames(timing_toast_sup[[lot]]$selec), ncol(timing_toast_sup[[lot]]$selec)),
                  "sim"=rep(colnames(timing_toast_sup[[lot]]$selec), each=nrow(timing_toast_sup[[lot]]$selec)),
                  "dataset"=name_data[lot],
                  "feat_selec"='TOAST')
  df2 = data.frame("values"=c(timing_toast_sup[[lot]]$restricted),
                  "score"="time",
                  "deconv"=rep(rownames(timing_toast_sup[[lot]]$restricted), ncol(timing_toast_sup[[lot]]$restricted)),
                  "sim"=rep(colnames(timing_toast_sup[[lot]]$restricted), each=nrow(timing_toast_sup[[lot]]$restricted)),
                  "dataset"=name_data[lot],
                  "feat_selec"='None_restrictedT')
  bind_rows(df,df2)
})
timing_toast_sup_df = do.call(rbind,timing_toast_sup_df)
saveRDS(timing_toast_sup_df, paste0("perf_scores/",
                              date,
                              "_time_",
                              block,
                              "_toast_sup.rds"))

timing_cimpleg_sup = lapply(seq(2,n_lot), function (lot) {
  res = sapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res0 = sapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("timing/CimpleG/",block,"/sup/",
                             date,"_",name_data[lot],
                             "_timing_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        time=readRDS(paste0("timing/CimpleG/",block,"/sup/",
                            date,"_",name_data[lot],
                            "_timing_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        time2=readRDS(paste0("timing/CimpleG/",block,"/sup/",
                             date,"_",name_data[lot],
                             "_timingbl2_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))
        vec = c(time,time2)
        names(vec) = c("selec","restricted")
        return(vec)
      }
      else {print(paste0(meth," not run yet"))}
    })
    vec = c(res0)
    names(vec) = c(sapply(colnames(res0), function(x)
      sapply(rownames(res0), function(y) paste(x,y))))
    return(vec)
  })
  res1 = res[grep("selec",rownames(res)),]
  res2 = res[grep("restricted",rownames(res)),]
  rownames(res1) = deconv_methods_sup
  rownames(res2) = deconv_methods_sup
  colnames(res1) = seq(n_sims)
  colnames(res2) = seq(n_sims)
  return(list("selec"=res1,"restricted"=res2))
})
timing_cimpleg_sup_df = lapply(seq(2,n_lot), function (lot) {
  df = data.frame("values"=c(timing_cimpleg_sup[[lot-1]]$selec),
                  "score"="time",
                  "deconv"=rep(rownames(timing_cimpleg_sup[[lot-1]]$selec), ncol(timing_cimpleg_sup[[lot-1]]$selec)),
                  "sim"=rep(colnames(timing_cimpleg_sup[[lot-1]]$selec), each=nrow(timing_cimpleg_sup[[lot-1]]$selec)),
                  "dataset"=name_data[lot],
                  "feat_selec"='CimpleG')
  df2 = data.frame("values"=c(timing_cimpleg_sup[[lot-1]]$restricted),
                   "score"="time",
                   "deconv"=rep(rownames(timing_cimpleg_sup[[lot-1]]$restricted), ncol(timing_cimpleg_sup[[lot-1]]$restricted)),
                   "sim"=rep(colnames(timing_cimpleg_sup[[lot-1]]$restricted), each=nrow(timing_cimpleg_sup[[lot-1]]$restricted)),
                   "dataset"=name_data[lot],
                   "feat_selec"='None_restrictedC')
  bind_rows(df,df2)
})
timing_cimpleg_sup_df = do.call(rbind,timing_cimpleg_sup_df)
if (nrow(timing_cimpleg_sup_df)>0) {saveRDS(timing_cimpleg_sup_df, paste0("perf_scores/",
                                    date,
                                    "_time_",
                                    block,
                                    "_cimpleg_sup.rds"))}

## ----
## Timing ref-free
## ----
tmp=lapply(seq(n_lot), function (lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("timing/TOAST/",block,"/unsup/",
                             date,"_",name_data[lot],
                             "_timing_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        time=readRDS(paste0("timing/TOAST/",block,"/unsup/",
                            date,"_",name_data[lot],
                            "_timing_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        time2=readRDS(paste0("timing/TOAST/",block,"/unsup/",
                            date,"_",name_data[lot],
                            "_timingbl2_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        if (is.list(time)) {
          time=time[[4]]
          time=strsplit(time," sec")[[1]][1]
          time=as.numeric(strsplit(time, " ")[[1]][2])
          saveRDS(time, paste0("timing/TOAST/",block,"/unsup/",
                               date,"_",name_data[lot],
                               "_timing_",
                               meth,
                               "_",sim_txt,
                               sim,".rds"))
        }
        if (is.list(time2)) {
          time2=time2[[4]]
          time2=strsplit(time2," sec")[[1]][1]
          time2=as.numeric(strsplit(time2, " ")[[1]][2])
          saveRDS(time2, paste0("timing/TOAST/",block,"/unsup/",
                               date,"_",name_data[lot],
                               "_timingbl2_",
                               meth,
                               "_",sim_txt,
                               sim,".rds"))
        }
      }
      if (file.exists(paste0("timing/CimpleG/",block,"/unsup/",
                             date,"_",name_data[lot],
                             "_timing_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        time=readRDS(paste0("timing/CimpleG/",block,"/unsup/",
                            date,"_",name_data[lot],
                            "_timing_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        time2=readRDS(paste0("timing/CimpleG/",block,"/unsup/",
                            date,"_",name_data[lot],
                            "_timingbl2_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        if (is.list(time)) {
          time=time[[4]]
          time=strsplit(time," sec")[[1]][1]
          time=as.numeric(strsplit(time, " ")[[1]][2])
          saveRDS(time, paste0("timing/CimpleG/",block,"/unsup/",
                               date,"_",name_data[lot],
                               "_timing_",
                               meth,
                               "_",sim_txt,
                               sim,".rds"))
        }
        if (is.list(time2)) {
          time2=time2[[4]]
          time2=strsplit(time2," sec")[[1]][1]
          time2=as.numeric(strsplit(time2, " ")[[1]][2])
          saveRDS(time2, paste0("timing/CimpleG/",block,"/unsup/",
                               date,"_",name_data[lot],
                               "_timingbl2_",
                               meth,
                               "_",sim_txt,
                               sim,".rds"))
        }
      }
    })
  })
})

timing_toast_unsup = lapply(seq(n_lot), function (lot) {
  res = sapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res0 = sapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("timing/TOAST/",block,"/unsup/",
                             date,"_",name_data[lot],
                             "_timing_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        time=readRDS(paste0("timing/TOAST/",block,"/unsup/",
                            date,"_",name_data[lot],
                            "_timing_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        time2=readRDS(paste0("timing/TOAST/",block,"/unsup/",
                             date,"_",name_data[lot],
                             "_timingbl2_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))
        vec = c(time,time2)
        names(vec) = c("selec","restricted")
        return(vec)
      }
      else {print(paste0(meth," not run yet"))}
    })
    vec = c(res0)
    names(vec) = c(sapply(colnames(res0), function(x)
      sapply(rownames(res0), function(y) paste(x,y))))
    return(vec)
  })
  res1 = res[grep("selec",rownames(res)),]
  res2 = res[grep("restricted",rownames(res)),]
  rownames(res1) = deconv_methods_unsup
  rownames(res2) = deconv_methods_unsup
  colnames(res1) = seq(n_sims)
  colnames(res2) = seq(n_sims)
  return(list("selec"=res1,"restricted"=res2))
})
timing_toast_unsup_df = lapply(seq(n_lot), function (lot) {
  df = data.frame("values"=c(timing_toast_unsup[[lot]]$selec),
                  "score"="time",
                  "deconv"=rep(rownames(timing_toast_unsup[[lot]]$selec), ncol(timing_toast_unsup[[lot]]$selec)),
                  "sim"=rep(colnames(timing_toast_unsup[[lot]]$selec), each=nrow(timing_toast_unsup[[lot]]$selec)),
                  "dataset"=name_data[lot],
                  "feat_selec"='TOAST')
  df2 = data.frame("values"=c(timing_toast_unsup[[lot]]$restricted),
                   "score"="time",
                   "deconv"=rep(rownames(timing_toast_unsup[[lot]]$restricted), ncol(timing_toast_unsup[[lot]]$restricted)),
                   "sim"=rep(colnames(timing_toast_unsup[[lot]]$restricted), each=nrow(timing_toast_unsup[[lot]]$restricted)),
                   "dataset"=name_data[lot],
                   "feat_selec"='None_restrictedT')
  bind_rows(df,df2)
})
timing_toast_unsup_df = do.call(rbind,timing_toast_unsup_df)
saveRDS(timing_toast_unsup_df, paste0("perf_scores/",
                                date,
                                "_time_",
                                block,
                                "_toast_unsup.rds"))

timing_cimpleg_unsup = lapply(seq(2,n_lot), function (lot) {
  res = sapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    res0 = sapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("timing/CimpleG/",block,"/unsup/",
                             date,"_",name_data[lot],
                             "_timing_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        time=readRDS(paste0("timing/CimpleG/",block,"/unsup/",
                            date,"_",name_data[lot],
                            "_timing_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        time2=readRDS(paste0("timing/CimpleG/",block,"/unsup/",
                             date,"_",name_data[lot],
                             "_timingbl2_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))
        vec = c(time,time2)
        names(vec) = c("selec","restricted")
        return(vec)
      }
      else {print(paste0(meth," not run yet"))}
    })
    vec = c(res0)
    names(vec) = c(sapply(colnames(res0), function(x)
      sapply(rownames(res0), function(y) paste(x,y))))
    return(vec)
  })
  res1 = res[grep("selec",rownames(res)),]
  res2 = res[grep("restricted",rownames(res)),]
  rownames(res1) = deconv_methods_unsup
  rownames(res2) = deconv_methods_unsup
  colnames(res1) = seq(n_sims)
  colnames(res2) = seq(n_sims)
  return(list("selec"=res1,"restricted"=res2))
})
timing_cimpleg_unsup_df = lapply(seq(2,n_lot), function (lot) {
  df = data.frame("values"=c(timing_cimpleg_unsup[[lot-1]]$selec),
                  "score"="time",
                  "deconv"=rep(rownames(timing_cimpleg_unsup[[lot-1]]$selec), ncol(timing_cimpleg_unsup[[lot-1]]$selec)),
                  "sim"=rep(colnames(timing_cimpleg_unsup[[lot-1]]$selec), each=nrow(timing_cimpleg_unsup[[lot-1]]$selec)),
                  "dataset"=name_data[lot],
                  "feat_selec"='CimpleG')
  df2 = data.frame("values"=c(timing_cimpleg_unsup[[lot-1]]$restricted),
                   "score"="time",
                   "deconv"=rep(rownames(timing_cimpleg_unsup[[lot-1]]$restricted), ncol(timing_cimpleg_unsup[[lot-1]]$restricted)),
                   "sim"=rep(colnames(timing_cimpleg_unsup[[lot-1]]$restricted), each=nrow(timing_cimpleg_unsup[[lot-1]]$restricted)),
                   "dataset"=name_data[lot],
                   "feat_selec"='None_restrictedC')
  bind_rows(df,df2)
})
timing_cimpleg_unsup_df = do.call(rbind,timing_cimpleg_unsup_df)
if (nrow(timing_cimpleg_unsup_df)>0) {saveRDS(timing_cimpleg_unsup_df, paste0("perf_scores/",
                                      date,
                                      "_time_",
                                      block,
                                      "_cimpleg_unsup.rds"))}
## For datasets with no ground truth, we can compare similarities between methods, as opposed to known similarities
