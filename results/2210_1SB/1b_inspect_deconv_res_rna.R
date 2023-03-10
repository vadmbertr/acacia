## ----
## Set parameters
## ----
block="rna"

## ----
## load A_true and A_pred
## ----
input_path = paste0("../2210_0simu/simulations/",block,"/")
input_path_matrix = sort(list.files(input_path, pattern = "sim"))
A_true = lapply(input_path_matrix, function(x) readRDS(paste0(input_path,x))$A_ref)
names(A_true) <- input_path_matrix
input_path_T = sort(list.files(input_path, pattern = paste0("T_",block,"_ref")))
n_lot=length(input_path_T)
col_names=colnames(A_true[[1]])
row_names=sapply(A_true[10*seq(n_lot)],rownames)
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

## ----
## Checking rows and cols of supervised methods
## ----
tmp=sapply(seq(n_lot), function (lot) {
  print(paste0("Checking dataset n°",lot,"/",n_lot))
  sapply(seq(length(A_true)/n_lot), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    sapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("deconv/",block,"/sup/",strsplit(input_path_T[lot],paste0("_T_",block,"_ref"))[[1]][1],"_Apred_",meth,"_",sim_txt,sim,".rds"))) {
        A_pred=readRDS(paste0("deconv/",block,"/sup/",strsplit(input_path_T[lot],paste0("_T_",block,"_ref"))[[1]][1],"_Apred_",meth,"_",sim_txt,sim,".rds"))
        if (!all(colnames(A_pred)==col_names)) {
          print(paste0("rerun method ",meth," for sim ",sim," of dataset ",lot," to reorder columns")) 
        }
        if (!all(rownames(A_pred)==row_names[[lot]])) {
          print(paste0("rerun method ",meth," for sim ",sim," of dataset ",lot," to reorder rows")) 
        }
        if (min(A_pred)<0) {
          print(paste0("rerun method ",meth," for sim ",sim," of dataset ",lot," to implement NN constraint")) 
        }
        if (!all(colSums(A_pred)<(1+1e-10)) | !all(colSums(A_pred)>(1-1e-10))) {
          print(paste0("rerun method ",meth," for sim ",sim," of dataset ",lot," to implement STO constraint")) 
        }
      }
      else {print(paste0(meth," not run yet"))}
      })
    })
  })

## ----
## Checking rows and cols of unsupervised methods
## ----
tmp=sapply(seq(n_lot), function (lot) {
  print(paste0("Checking dataset n°",lot,"/",n_lot))
  sapply(seq(length(A_true)/n_lot), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    sapply(deconv_methods_unsup, function(meth) {
      if (file.exists(paste0("deconv/",block,"/unsup/",
                             strsplit(input_path_T[lot],paste0("_T_",block,"_ref"))[[1]][1],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("deconv/",block,"/unsup/",
                              strsplit(input_path_T[lot],paste0("_T_",block,"_ref"))[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
        if (!all(colnames(A_pred)==col_names)) {
          print(paste0("rerun method ",meth,
                       " for sim ",sim,
                       " of dataset ",lot,
                       " to reorder columns")) 
        }
        if (!all(rownames(A_pred)==row_names[[lot]])) {
          print(paste0("rerun method ",meth,
                       " for sim ",sim,
                       " of dataset ",lot,
                       " to reorder rows")) 
        }
      }
      else {print(paste0(meth," not run yet"))}
    })
  })
})

