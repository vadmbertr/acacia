## ----
## Set parameters
## ----
block="met" # to vary, either met or rna
pattern=paste0("_T_",block,"_ref")

## ----
## load timings
## ----
input_path_nlot = paste0("../2210_0simu/simulations/",block,"/")
input_path_nlot = sort(list.files(input_path_nlot, pattern = pattern))
n_lot=length(input_path_nlot)
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
input_path_nsim = paste0("timing/",block,"/sup/")
input_path_nsim = sort(list.files(input_path_nsim, pattern = deconv_methods_sup[1]))
n_sim = length(input_path_nsim)/n_lot
rm(input_path_nsim)

timing_sup <- lapply(seq(n_lot), function (lot) {
  tmp=sapply(seq(n_sim), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    tmp=sapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("timing/",block,"/sup/",
                             strsplit(input_path_nlot[lot],pattern)[[1]][1],
                             "_timing_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        readRDS(paste0("timing/",block,"/sup/",
                       strsplit(input_path_nlot[lot],pattern)[[1]][1],
                       "_timing_",
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
                 "simulation"=rep(colnames(tmp),each=nrow(tmp)))
  tmp$dataset = rev(strsplit(strsplit(input_path_nlot[lot],pattern)[[1]][1],"_")[[1]])[1]
  return(tmp)
})
timing_sup <- do.call(rbind,timing_sup)

timing_unsup <- lapply(seq(n_lot), function (lot) {
  tmp=sapply(seq(n_sim), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    tmp=sapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("timing/",block,"/unsup/",
                             strsplit(input_path_nlot[lot],pattern)[[1]][1],
                             "_timing_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        readRDS(paste0("timing/",block,"/unsup/",
                       strsplit(input_path_nlot[lot],pattern)[[1]][1],
                       "_timing_",
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
                 "simulation"=rep(colnames(tmp),each=nrow(tmp)))
  tmp$dataset = rev(strsplit(strsplit(input_path_nlot[lot],pattern)[[1]][1],"_")[[1]])[1]
  return(tmp)
})
timing_unsup <- do.call(rbind,timing_unsup)

## ----
## plot timings
## ----
library(ggplot2)
library(see)
timing_sup$timing <- as.numeric(timing_sup$timing)
timing_unsup$timing <- as.numeric(timing_unsup$timing)
ggplot(timing_sup, aes(x=method, y=timing, color=simulation)) +
  geom_point() +
  geom_line(aes(group=simulation)) +
  facet_wrap(~dataset) +
  scale_y_log10() +
  scale_color_pizza_d() +
  theme(axis.text.x = element_text(angle=45))
ggsave(paste0("timing_viz/",
              strsplit(input_path_nlot[1],"_")[[1]][1],
              "_",block,
              "_timing_sup.pdf"),
       width = 10, height = 4)
ggplot(timing_unsup, aes(x=method, y=timing, color=simulation)) +
  geom_point() +
  geom_line(aes(group=simulation)) +
  facet_wrap(~dataset) +
  scale_y_log10() +
  scale_color_pizza_d() +
  theme(axis.text.x = element_text(angle=45))
ggsave(paste0("timing_viz/",
              strsplit(input_path_nlot[1],"_")[[1]][1],
              "_",block,
              "_timing_unsup.pdf"),
       width = 8, height = 4)
## For datasets with no ground truth, we can compare similarities between methods, as opposed to known similarities
