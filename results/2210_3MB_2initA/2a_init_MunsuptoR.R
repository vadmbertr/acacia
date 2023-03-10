## ----
## Set parameters
## ----
block1='met'
block2='rna'
pattern1=paste0("T_",block1,"_ref")
source("../../src/2MB_deconv_functions.R")
meth_deconv_MtoR <- c("ICA","NMF")

## ----
## load Dsim from RNA and Atrue
## ----
input_path = paste0("../2210_0simu/simulations/",block2,"/")
input_path_matrix = sort(list.files(input_path, pattern = "sim"))
Dsim = lapply(input_path_matrix, function(x) readRDS(paste0(input_path,x))$D_rna_sim)
Atrue = lapply(input_path_matrix, function(x) readRDS(paste0(input_path,x))$A_ref)

input_path = paste0("../2210_0simu/simulations/",block1,"/")
input_path_T = sort(list.files(input_path, pattern = pattern1))
n_lot = length(input_path_T)
input_path_nsims = sort(list.files(input_path, pattern = "sim"))
n_sims = length(input_path_nsims)/n_lot
name_data = unname(sapply(input_path_T,function(x) rev(strsplit(strsplit(x,pattern1)[[1]][1],"_")[[1]])[1]))
rm(input_path,input_path_nsims)

Atrue = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) Atrue[[10*(lot-1)+sim]])
})
names(Dsim) <- input_path_matrix
rm(input_path_matrix)

## ----
## load Apred from MET
## ----
input_path_methods_unsup = paste0("../2210_1SB/deconv/",block1,"/unsup/")
input_path_methods_unsup = sort(list.files(input_path_methods_unsup))
deconv_methods_in = unique(sapply(input_path_methods_unsup, function(x)
  rev(strsplit(x,"_")[[1]])[2]))
rm(input_path_methods_unsup)

Apred_in=lapply(seq(n_lot), function(lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_in, function(meth) {
      if (file.exists(paste0("../2210_1SB/deconv/",block1,"/unsup/",
                             strsplit(input_path_T[lot],pattern1)[[1]][1],
                             "Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("../2210_1SB/deconv/",block1,"/unsup/",
                              strsplit(input_path_T[lot],pattern1)[[1]][1],
                              "Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
      }
      else {print(paste0(meth," not run yet"))}
    })
  })
})

## ----
## Deconvolution for RNA
## ----
tmp=lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_in, function(methM) {
      lapply(meth_deconv_MtoR, function(methR) {
        if (!file.exists(paste0("2a_init_MunsuptoR/",
                                strsplit(input_path_T[lot],pattern1)[[1]][1],
                                "Apred_",
                                methM,
                                "_to_",
                                methR,
                                "_",sim_txt,
                                sim,".rds"))) {
          print(paste0("Using methM ",methM," to init ",methR," in sim n°",sim))
          D=Dsim[[(lot-1)*10+sim]]
          A=Atrue[[lot]][[sim]]
          Ainter=Apred_in[[lot]][[sim]][[which(deconv_methods_in==methM)]]
          A=A[,colnames(D)]
          #Ainter=Ainter[,colnames(D)]
          deconv_MtoR=run_MtoR_deconvolutionA(methR, D, Ainter, A)
          saveRDS(deconv_MtoR,paste0("2a_init_MunsuptoR/",
                                     strsplit(input_path_T[lot],pattern1)[[1]][1],
                                     "Apred_",
                                     methM,
                                     "_to_",
                                     methR,
                                     "_",sim_txt,
                                     sim,".rds"))
        }
      })
    })
  })
})

