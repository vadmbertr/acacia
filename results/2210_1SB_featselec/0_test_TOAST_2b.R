## ----
## Set parameters
## ----
library(TOAST)
block="rna" # to vary, either met or rna

## ----
## load
## ----
input_path_methods = "0_test_TOAST/deconv/"
input_path_methodssup = sort(list.files(input_path_methods, pattern = paste0("deconv_",block,"_sup")))
deconv_methods_sup = unique(sapply(input_path_methodssup, function(x)
  strsplit(rev(strsplit(x,"_")[[1]])[1],".rds")[[1]][1]))
deconv_methods_sup[deconv_methods_sup=="net"] = "elastic_net"
rm(input_path_methodssup)

input_path_methodsunsup = sort(list.files(input_path_methods, pattern = paste0("deconv_",block,"_unsup")))
deconv_methods_unsup = unique(sapply(input_path_methodsunsup, function(x)
  strsplit(rev(strsplit(x,"_")[[1]])[1],".rds")[[1]][1]))
rm(input_path_methodsunsup,input_path_methods)

tmp=lapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("0_test_TOAST/deconv/deconvbl2_",block,"_unsup_",
                             meth,".rds"))) {
        res=readRDS(paste0("0_test_TOAST/deconv/deconvbl2_",block,"_unsup_",
                           meth,".rds"))
        time=res$time_elapsed
        if (is.list(time)) {
          print(meth)
          time=time[[4]]
          time=strsplit(time," sec")[[1]][1]
          time=as.numeric(strsplit(time, " ")[[1]][2])
          res=list(res=res$res, time_elapsed=time)
          saveRDS(res, paste0("0_test_TOAST/deconv/deconvbl2_",block,"_unsup_",
                               meth,".rds"))
        }
      }
  })

## ----
## Score ref-based
## ----
score_perf_sup = sapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("0_test_TOAST/deconv/deconv_",block,"_sup_",
                             meth,".rds"))) {
        readRDS(paste0("0_test_TOAST/deconv/deconv_",block,"_sup_",
                              meth,".rds"))$time_elapsed
        }
      else {print(paste0(meth," not run yet"))}
  })
score_perf_sup_df = data.frame("values"=c(score_perf_sup),
                                     "score"='time',
                                     "feat_selec"='TOAST',
                                     "deconv"=deconv_methods_sup)
saveRDS(score_perf_sup_df, paste0("0_test_TOAST/scores/",
                                        "time_",block,"_sup.rds"))

score_perf_sup_bl = sapply(deconv_methods_sup, function(meth) {
  if (file.exists(paste0("0_test_TOAST/deconv/deconvbl_",block,"_sup_",
                         meth,".rds"))) {
    readRDS(paste0("0_test_TOAST/deconv/deconvbl_",block,"_sup_",
                          meth,".rds"))$time_elapsed
  }
  else {print(paste0(meth," not run yet"))}
})
score_perf_sup_df_bl = data.frame("values"=c(score_perf_sup_bl),
                                        "score"='time',
                                        "feat_selec"='None',
                                        "deconv"=deconv_methods_sup)
saveRDS(score_perf_sup_df_bl, paste0("0_test_TOAST/scores/",
                                        "timebl_",block,"_sup.rds"))

score_perf_sup_bl2 = sapply(deconv_methods_sup, function(meth) {
  if (file.exists(paste0("0_test_TOAST/deconv/deconvbl2_",block,"_sup_",
                         meth,".rds"))) {
    readRDS(paste0("0_test_TOAST/deconv/deconvbl2_",block,"_sup_",
                   meth,".rds"))$time_elapsed
  }
  else {print(paste0(meth," not run yet"))}
})
score_perf_sup_df_bl2 = data.frame("values"=c(score_perf_sup_bl2),
                                        "score"='time',
                                        "feat_selec"='None_restricted',
                                        "deconv"=deconv_methods_sup)
saveRDS(score_perf_sup_df_bl2, paste0("0_test_TOAST/scores/",
                                           "timebl2_",block,"_sup.rds"))

## ----
## Score ref-free
## ----
score_perf_unsup = sapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("0_test_TOAST/deconv/deconv_",block,"_unsup_",
                             meth,".rds"))) {
        readRDS(paste0("0_test_TOAST/deconv/deconv_",block,"_unsup_",
                              meth,".rds"))$time_elapsed
      }
      else {print(paste0(meth," not run yet"))}
  })
score_perf_unsup_df = data.frame("values"=c(score_perf_unsup),
                                       "score"='time',
                                       "feat_selec"='TOAST',
                                       "deconv"=deconv_methods_unsup)
saveRDS(score_perf_unsup_df, paste0("0_test_TOAST/scores/",
                                          "time_",block,"_unsup.rds"))

score_perf_unsup_bl = sapply(deconv_methods_unsup, function(meth) {
  if (file.exists(paste0("0_test_TOAST/deconv/deconvbl_",block,"_unsup_",
                         meth,".rds"))) {
    A_pred=readRDS(paste0("0_test_TOAST/deconv/deconvbl_",block,"_unsup_",
                          meth,".rds"))$time_elapsed
  }
  else {print(paste0(meth," not run yet"))}
})
score_perf_unsup_df_bl = data.frame("values"=c(score_perf_unsup_bl),
                                          "score"='time',
                                          "feat_selec"='None',
                                          "deconv"=deconv_methods_unsup)
saveRDS(score_perf_unsup_df_bl, paste0("0_test_TOAST/scores/",
                                          "timebl_",block,"_unsup.rds"))

score_perf_unsup_bl2 = sapply(deconv_methods_unsup, function(meth) {
  if (file.exists(paste0("0_test_TOAST/deconv/deconvbl2_",block,"_unsup_",
                         meth,".rds"))) {
    A_pred=readRDS(paste0("0_test_TOAST/deconv/deconvbl2_",block,"_unsup_",
                          meth,".rds"))$time_elapsed
  }
  else {print(paste0(meth," not run yet"))}
})
score_perf_unsup_df_bl2 = data.frame("values"=c(score_perf_unsup_bl2),
                                          "score"='time',
                                          "feat_selec"='None_restricted',
                                          "deconv"=deconv_methods_unsup)
saveRDS(score_perf_unsup_df_bl2, paste0("0_test_TOAST/scores/",
                                             "timebl2_",block,"_unsup.rds"))

