## ----
## Set parameters
## ----
block="rna"

## ----
## load Apred
## ----
input_path_nlot = paste0("../2210_0simu/simulations/",block)
input_path_nlot = sort(list.files(input_path_nlot, pattern = paste0("T_",block,"_ref")))
n_lot=length(input_path_nlot)
input_path_methodssup = paste0("deconv/",block,"/sup/")
input_path_methodssup = sort(list.files(input_path_methodssup))
deconv_methods_sup = unique(sapply(input_path_methodssup, function(x)
  rev(strsplit(x,"_")[[1]])[2]))
deconv_methods_sup[deconv_methods_sup=="net"] = "elastic_net"
rm(input_path_methodssup)
input_path_methodsunsup = paste0("deconv/",block,"/unsup/")
input_path_methodsunsup = sort(list.files(input_path_methodsunsup))
deconv_methods_unsup = unique(sapply(input_path_methodsunsup, function(x)
  rev(strsplit(x,"_")[[1]])[2]))
rm(input_path_methodsunsup)
input_path_nsim = paste0("deconv/",block,"/unsup/")
input_path_nsim = sort(list.files(input_path_nsim, pattern = deconv_methods_unsup[1]))
n_sim = length(input_path_nsim)/n_lot
rm(input_path_nsim)

Apred_sup <- lapply(seq(n_lot), function (lot) {
  tmp=lapply(seq(n_sim), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    tmp=lapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("deconv/",block,"/sup/",
                             strsplit(input_path_nlot[lot],paste0("_T_",block,"_ref"))[[1]][1],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        readRDS(paste0("deconv/",block,"/sup/",
                       strsplit(input_path_nlot[lot],paste0("_T_",block,"_ref"))[[1]][1],
                       "_Apred_",
                       meth,
                       "_",sim_txt,
                       sim,".rds"))
      }
      else {print(paste0(meth," not run yet"))}
    })
    names(tmp)=deconv_methods_sup
    return(tmp)
  })
  names(tmp)=paste0("sim",seq(n_sim))
  return(tmp)
})
Apred_sup <- lapply(Apred_sup, function(lot)
  lapply(lot, function(sim) {
    data.frame(do.call(cbind,lapply(sim, function(mat) c(mat))))
  }))

Apred_unsup <- lapply(seq(n_lot), function (lot) {
  tmp=lapply(seq(n_sim), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    tmp=lapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("deconv/",block,"/unsup/",
                             strsplit(input_path_nlot[lot],paste0("_T_",block,"_ref"))[[1]][1],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        readRDS(paste0("deconv/",block,"/unsup/",
                       strsplit(input_path_nlot[lot],paste0("_T_",block,"_ref"))[[1]][1],
                       "_Apred_",
                       meth,
                       "_",sim_txt,
                       sim,".rds"))
      }
      else {print(paste0(meth," not run yet"))}
    })
    names(tmp)=deconv_methods_unsup
    return(tmp)
  })
  names(tmp)=paste0("sim",seq(n_sim))
  return(tmp)
})
Apred_unsup <- lapply(Apred_unsup, function(lot)
  lapply(lot, function(sim) {
    data.frame(do.call(cbind,lapply(sim, function(mat) c(mat))))
    }))

## ----
## Compare similar methods
## ----
source("../../src/score_functions.R")
library(ggplot2)
library(see)
library(ggpubr)
similar_methods_sup <- list(svr=c("svr3","svr4"),
                            CIBERSORT=c("CIBERSORT3","CIBERSORT4"))
for (couple in names(similar_methods_sup)) {
  for (lot in seq(n_lot)) {
    p=ggarrange(plotlist=lapply(Apred_sup[[lot]], function(sim)
      ggplot(sim, aes_string(x=similar_methods_sup[[couple]][1], y=similar_methods_sup[[couple]][2])) +
        geom_point() +
        theme_modern()),
      common.legend = T)
    p
    ggsave(paste0("similar_methods_rna/",
                  strsplit(input_path_nlot[lot],paste0("_T_",block,"_ref"))[[1]][1],
                  "_",
                  couple,
                  ".pdf"), width=10, height=6)
  }
}
