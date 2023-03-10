## ----
## Set parameters
## ----
block1='rna'
block2='met'
pattern1=paste0("T_",block1,"_ref")
source("../../src/2MB_deconv_functions.R")
meth_deconv_RtoM <- c("ICA","RefFreeEWAS","MeDeCom")

## ----
## load Dsim from MET and Atrue
## ----
input_path = paste0("../2210_0simu/simulations/",block2,"/")
input_path_matrix = sort(list.files(input_path, pattern = "sim"))
Dsim = lapply(input_path_matrix, function(x) readRDS(paste0(input_path,x))$D_met_sim)

input_path = paste0("../2210_0simu/simulations/",block1,"/")
input_path_matrix = sort(list.files(input_path, pattern = "sim"))
Atrue = lapply(input_path_matrix, function(x) readRDS(paste0(input_path,x))$A_ref)

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
## load Apred from RNA
## ----
input_path_methods_sup = paste0("../2210_1SB/deconv/",block1,"/sup/")
input_path_methods_sup = sort(list.files(input_path_methods_sup))
deconv_methods_in = unique(sapply(input_path_methods_sup, function(x)
  rev(strsplit(x,"_")[[1]])[2]))
deconv_methods_in[deconv_methods_in=='net'] = 'elastic_net'
deconv_methods_in <- deconv_methods_in[-grep("4",deconv_methods_in)]
rm(input_path_methods_sup)

Apred_in=lapply(seq(n_lot), function(lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_in, function(meth) {
      if (file.exists(paste0("../2210_1SB/deconv/",block1,"/sup/",
                             strsplit(input_path_T[lot],pattern1)[[1]][1],
                             "Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("../2210_1SB/deconv/",block1,"/sup/",
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
## Deconvolution for MET
## ----
tmp=lapply(seq(n_lot), function (lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_in, function(methR) {
      lapply(meth_deconv_RtoM, function(methM) {
        if (!file.exists(paste0("2a_init_RsuptoM/",
                                strsplit(input_path_T[lot],pattern1)[[1]][1],
                                "Apred_",
                                methR,
                                "_to_",
                                methM,
                                "_",sim_txt,
                                sim,".rds"))) {
          print(paste0("Using methR ",methR," to init ",methM," in sim n°",sim))
          D=Dsim[[(lot-1)*10+sim]]
          if (nrow(D)>3e4) {
            D=D[TOAST::findRefinx(D, nmarker = 3e4),]
          }
          A=Atrue[[lot]][[sim]]
          Ainter=Apred_in[[lot]][[sim]][[which(deconv_methods_in==methR)]]
          A=A[,colnames(D)]
          Ainter=Ainter[,colnames(D)]
          deconv_RtoM=run_RtoM_deconvolutionA(methM, D, Ainter, A)
          saveRDS(deconv_RtoM,paste0("2a_init_RsuptoM/",
                                     strsplit(input_path_T[lot],pattern1)[[1]][1],
                                     "Apred_",
                                     methR,
                                     "_to_",
                                     methM,
                                     "_",sim_txt,
                                     sim,".rds"))
        }
      })
    })
  })
})

