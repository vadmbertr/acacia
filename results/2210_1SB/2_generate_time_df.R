## ----
## Set parameters
## ----
block="rna" # to vary, either met or rna

## ----
## load timings
## ----
input_path = paste0("timing/",block,"/sup/")
input_path_T = sort(list.files(input_path))
date = strsplit(input_path_T[1],"_")[[1]][1]
rm(input_path)

input_path_methodssup = paste0("timing/",block,"/sup/")
input_path_methodssup = sort(list.files(input_path_methodssup))
deconv_methods_sup = unique(sapply(input_path_methodssup, function(x)
  rev(strsplit(x,"_")[[1]])[2]))
deconv_methods_sup[deconv_methods_sup=="net"] = "elastic_net"
rm(input_path_methodssup)

input_path_methodsunsup = paste0("timing/",block,"/unsup/")
input_path_methodsunsup = sort(list.files(input_path_methodsunsup))
deconv_methods_unsup = unique(sapply(input_path_methodsunsup, function(x)
  rev(strsplit(x,"_")[[1]])[2]))
rm(input_path_methodsunsup)

name_data = unique(unname(sapply(input_path_T,function(x) strsplit(strsplit(x,"_timing")[[1]][1],paste0(date,"_"))[[1]][2])))
n_lot = length(name_data)
n_sims = length(input_path_T)/length(deconv_methods_sup)/n_lot

## ----
## Timing ref-based
## ----
tmp=lapply(seq(n_lot), function (lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("timing/",block,"/sup/",
                             date,"_",name_data[lot],
                             "_timing_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        time=readRDS(paste0("timing/",block,"/sup/",
                            date,"_",name_data[lot],
                            "_timing_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        if (is.list(time)) {
          time=time[[4]]
          time=strsplit(time," sec")[[1]][1]
          time=as.numeric(strsplit(time, " ")[[1]][2])
          saveRDS(time, paste0("timing/",block,"/sup/",
                         date,"_",name_data[lot],
                         "_timing_",
                         meth,
                         "_",sim_txt,
                         sim,".rds"))
          }
      }
      else {print(paste0(meth," not run yet"))}
    })
  })
})
timing_sup = lapply(seq(n_lot), function (lot) {
  res = sapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    sapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("timing/",block,"/sup/",
                              date,"_",name_data[lot],
                              "_timing_",
                              meth,
                             "_",sim_txt,
                              sim,".rds"))) {
        time=readRDS(paste0("timing/",block,"/sup/",
                              date,"_",name_data[lot],
                              "_timing_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        return(time)
        }
      else {print(paste0(meth," not run yet"))}
      })
    })
  colnames(res) = seq(n_sims)
  res
  })

timing_sup_df = lapply(seq(n_lot), function (lot) {
  df = data.frame("values"=c(timing_sup[[lot]]),
                  "score"="time",
                  "deconv"=rep(rownames(timing_sup[[lot]]), ncol(timing_sup[[lot]])),
                  "sim"=rep(colnames(timing_sup[[lot]]), each=nrow(timing_sup[[lot]])),
                  "dataset"=name_data[lot])
})
timing_sup_df = do.call(rbind,timing_sup_df)
saveRDS(timing_sup_df, paste0("perf_scores/",
                              date,
                              "_time_",
                              block,
                              "_sup.rds"))

## ----
## Timing ref-free
## ----
tmp=lapply(seq(n_lot), function (lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("timing/",block,"/unsup/",
                             date,"_",name_data[lot],
                             "_timing_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        time=readRDS(paste0("timing/",block,"/unsup/",
                            date,"_",name_data[lot],
                            "_timing_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        if (is.list(time)) {
          time=time[[4]]
          time=strsplit(time," sec")[[1]][1]
          time=as.numeric(strsplit(time, " ")[[1]][2])
          saveRDS(time, paste0("timing/",block,"/unsup/",
                               date,"_",name_data[lot],
                               "_timing_",
                               meth,
                               "_",sim_txt,
                               sim,".rds"))
        }
      }
      else {print(paste0(meth," not run yet"))}
    })
  })
})
timing_unsup = lapply(seq(n_lot), function (lot) {
  res = sapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    sapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("timing/",block,"/unsup/",
                             date,"_",name_data[lot],
                             "_timing_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        time=readRDS(paste0("timing/",block,"/unsup/",
                            date,"_",name_data[lot],
                            "_timing_",
                            meth,
                            "_",sim_txt,
                            sim,".rds"))
        return(time)
      }
      else {print(paste0(meth," not run yet"))}
    })
  })
  colnames(res) = seq(n_sims)
  res
})

timing_unsup_df = lapply(seq(n_lot), function (lot) {
  df = data.frame("values"=c(timing_unsup[[lot]]),
                  "score"="time",
                  "deconv"=rep(rownames(timing_unsup[[lot]]), ncol(timing_unsup[[lot]])),
                  "sim"=rep(colnames(timing_unsup[[lot]]), each=nrow(timing_unsup[[lot]])),
                  "dataset"=name_data[lot])
})
timing_unsup_df = do.call(rbind,timing_unsup_df)
saveRDS(timing_unsup_df, paste0("perf_scores/",
                                date,
                                "_time_",
                                block,
                                "_unsup.rds"))
## For datasets with no ground truth, we can compare similarities between methods, as opposed to known similarities
