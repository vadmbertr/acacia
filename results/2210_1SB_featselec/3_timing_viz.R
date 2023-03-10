## ----
## Set parameters
## ----
block="met" # to vary, either met or rna
feat_selec="" # to vary, either "" or "bl2"
pattern=paste0("_T_",block,"_ref")

## ----
## load timings
## ----
input_path_nlot = paste0("../2210_0simu/simulations/",block,"/")
input_path_nlot = sort(list.files(input_path_nlot, pattern = pattern))
n_lot=length(input_path_nlot)
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
input_path_nsim = paste0("timing/TOAST/",block,"/sup/")
input_path_nsim = sort(list.files(input_path_nsim, pattern = deconv_methods_sup[1]))
n_sim = length(input_path_nsim)/n_lot/2
rm(input_path_nsim)

timing_toast_sup <- lapply(seq(n_lot), function (lot) {
  tmp=sapply(seq(n_sim), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    tmp=sapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("timing/TOAST/",block,"/sup/",
                             strsplit(input_path_nlot[lot],pattern)[[1]][1],
                             "_timing",feat_selec,"_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        readRDS(paste0("timing/TOAST/",block,"/sup/",
                       strsplit(input_path_nlot[lot],pattern)[[1]][1],
                       "_timing",feat_selec,"_",
                       meth,
                       "_",sim_txt,
                       sim,".rds"))
      }
      else {print(paste0(meth," not run yet"))}
    })
    tmp=data.frame(tmp)
    colnames(tmp)=paste0("sim",sim)
    return(tmp)
  })
  tmp=do.call(cbind,tmp)
  tmp=data.frame("timing"=c(tmp),
                 "method"=rep(deconv_methods_sup,ncol(tmp)),
                 "simulation"=rep(colnames(tmp),each=nrow(tmp)),
                 "feat_selec"='TOAST')
  tmp$dataset = rev(strsplit(strsplit(input_path_nlot[lot],pattern)[[1]][1],"_")[[1]])[1]
  return(tmp)
})
timing_toast_sup <- do.call(rbind,timing_toast_sup)

timing_cimpleg_sup <- lapply(seq(n_lot), function (lot) {
  tmp=sapply(seq(n_sim), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    tmp=sapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("timing/CimpleG/",block,"/sup/",
                             strsplit(input_path_nlot[lot],pattern)[[1]][1],
                             "_timing",feat_selec,"_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        readRDS(paste0("timing/CimpleG/",block,"/sup/",
                       strsplit(input_path_nlot[lot],pattern)[[1]][1],
                       "_timing",feat_selec,"_",
                       meth,
                       "_",sim_txt,
                       sim,".rds"))
      }
      else {print(paste0(meth," not run yet"))}
    })
    tmp=data.frame(tmp)
    colnames(tmp)=paste0("sim",sim)
    return(tmp)
  })
  tmp=do.call(cbind,tmp)
  tmp=data.frame("timing"=c(tmp),
                 "method"=rep(deconv_methods_sup,ncol(tmp)),
                 "simulation"=rep(colnames(tmp),each=nrow(tmp)),
                 "feat_selec"='CimpleG')
  tmp$dataset = rev(strsplit(strsplit(input_path_nlot[lot],pattern)[[1]][1],"_")[[1]])[1]
  return(tmp)
})
timing_cimpleg_sup <- do.call(rbind,timing_cimpleg_sup)
timing_cimpleg_sup <- timing_cimpleg_sup[-grep("not run yet",timing_cimpleg_sup$timing),]

timing_toast_unsup <- lapply(seq(n_lot), function (lot) {
  tmp=sapply(seq(n_sim), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    tmp=sapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("timing/TOAST/",block,"/unsup/",
                             strsplit(input_path_nlot[lot],pattern)[[1]][1],
                             "_timing",feat_selec,"_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        readRDS(paste0("timing/TOAST/",block,"/unsup/",
                       strsplit(input_path_nlot[lot],pattern)[[1]][1],
                       "_timing",feat_selec,"_",
                       meth,
                       "_",sim_txt,
                       sim,".rds"))
      }
      else {print(paste0(meth," not run yet"))}
    })
    tmp=data.frame(tmp)
    colnames(tmp)=paste0("sim",sim)
    return(tmp)
  })
  tmp=do.call(cbind,tmp)
  tmp=data.frame("timing"=c(tmp),
                 "method"=rep(deconv_methods_unsup,ncol(tmp)),
                 "simulation"=rep(colnames(tmp),each=nrow(tmp)),
                 "feat_selec"='TOAST')
  tmp$dataset = rev(strsplit(strsplit(input_path_nlot[lot],pattern)[[1]][1],"_")[[1]])[1]
  return(tmp)
})
timing_toast_unsup <- do.call(rbind,timing_toast_unsup)

timing_cimpleg_unsup <- lapply(seq(n_lot), function (lot) {
  tmp=sapply(seq(n_sim), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    tmp=sapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("timing/CimpleG/",block,"/unsup/",
                             strsplit(input_path_nlot[lot],pattern)[[1]][1],
                             "_timing",feat_selec,"_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        readRDS(paste0("timing/CimpleG/",block,"/unsup/",
                       strsplit(input_path_nlot[lot],pattern)[[1]][1],
                       "_timing",feat_selec,"_",
                       meth,
                       "_",sim_txt,
                       sim,".rds"))
      }
      else {print(paste0(meth," not run yet"))}
    })
    tmp=data.frame(tmp)
    colnames(tmp)=paste0("sim",sim)
    return(tmp)
  })
  tmp=do.call(cbind,tmp)
  tmp=data.frame("timing"=c(tmp),
                 "method"=rep(deconv_methods_unsup,ncol(tmp)),
                 "simulation"=rep(colnames(tmp),each=nrow(tmp)),
                 "feat_selec"='CimpleG')
  tmp$dataset = rev(strsplit(strsplit(input_path_nlot[lot],pattern)[[1]][1],"_")[[1]])[1]
  return(tmp)
})
timing_cimpleg_unsup <- do.call(rbind,timing_cimpleg_unsup)
timing_cimpleg_unsup <- timing_cimpleg_unsup[-grep("not run yet",timing_cimpleg_unsup$timing),]

## ----
## plot timings
## ----
library(ggplot2)
library(see)
timing_toast_sup$timing <- as.numeric(timing_toast_sup$timing)
timing_cimpleg_sup$timing <- as.numeric(timing_cimpleg_sup$timing)
timing_toast_unsup$timing <- as.numeric(timing_toast_unsup$timing)
timing_cimpleg_unsup$timing <- as.numeric(timing_cimpleg_unsup$timing)

ggplot(timing_toast_sup, aes(x=method, y=timing, color=simulation)) +
  geom_point() +
  geom_boxplot(aes(group=method), alpha=.1, color="gray") +
  geom_line(aes(group=simulation)) +
  facet_wrap(~dataset) +
  scale_y_log10() +
  scale_color_pizza_d() +
  theme(axis.text.x = element_text(angle=45))
ggsave(paste0("timing_viz/",
              strsplit(input_path_nlot[1],"_")[[1]][1],
              "_",block,
              "_timing",feat_selec,"_toast_sup.pdf"),
       width = 10, height = 4)
if (nrow(timing_cimpleg_sup)>0) {
  ggplot(timing_cimpleg_sup, aes(x=method, y=timing, color=simulation)) +
    geom_point() +
    geom_boxplot(aes(group=method), alpha=.1, color="gray") +
    geom_line(aes(group=simulation)) +
    facet_wrap(~dataset) +
    scale_y_log10() +
    scale_color_pizza_d() +
    theme(axis.text.x = element_text(angle=45))
  ggsave(paste0("timing_viz/",
              strsplit(input_path_nlot[1],"_")[[1]][1],
              "_",block,
              "_timing",feat_selec,"_cimpleg_sup.pdf"),
       width = 10, height = 4)}

ggplot(timing_toast_unsup, aes(x=method, y=timing, color=simulation)) +
  geom_point() +
  geom_boxplot(aes(group=method), alpha=.1, color="gray") +
  geom_line(aes(group=simulation)) +
  facet_wrap(~dataset) +
  scale_y_log10() +
  scale_color_pizza_d() +
  theme(axis.text.x = element_text(angle=45))
ggsave(paste0("timing_viz/",
              strsplit(input_path_nlot[1],"_")[[1]][1],
              "_",block,
              "_timing",feat_selec,"_toast_unsup.pdf"),
       width = 8, height = 4)
if (nrow(timing_cimpleg_sup)>0) {
  ggplot(timing_cimpleg_unsup, aes(x=method, y=timing, color=simulation)) +
    geom_point() +
    geom_boxplot(aes(group=method), alpha=.1, color="gray") +
    geom_line(aes(group=simulation)) +
    facet_wrap(~dataset) +
    scale_y_log10() +
    scale_color_pizza_d() +
    theme(axis.text.x = element_text(angle=45))
  ggsave(paste0("timing_viz/",
              strsplit(input_path_nlot[1],"_")[[1]][1],
              "_",block,
              "_timing",feat_selec,"_cimpleg_unsup.pdf"),
       width = 8, height = 4)}
## For datasets with no ground truth, we can compare similarities between methods, as opposed to known similarities

