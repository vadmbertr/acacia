## ----
## Set work space
## ----
library(TOAST)
library(CimpleG)
library(pbapply)
set.seed(2022)
K=list(4,5,9)
deconv_function <- function(dat, k) {
  res <- fastICA::fastICA(X = dat, n.comp = k, maxit = 1000, tol = 1e-09)
  res$names = row.names(dat)
  weighted.list <- deconica::generate_markers(df = res, n = 30, return = "gene.ranked")
  res$A_weighted = t(deconica::get_scores(res$X, weighted.list, summary = "weighted.mean", na.rm = TRUE))
  colnames(res$A_weighted) = colnames(dat)
  A = abs(res$A_weighted) %*% diag(1/colSums(abs(res$A_weighted)))
  return(t(A))
}

## ----
## load D MET & select TOAST features
## ----
input_path_met = sort(list.files("../2210_0simu/simulations/met/", pattern = "sim"))
input_path_nlot = sort(list.files("../2210_0simu/simulations/rna/", pattern = "T_rna_ref"))
n_lot = length(input_path_nlot)
n_sims = length(input_path_met)/n_lot
date = strsplit(input_path_met, "_")[[1]][1]
rm(input_path_nlot)

Dmet = lapply(input_path_met, function(x) readRDS(paste0("../2210_0simu/simulations/met/",x))$D_met_sim)
Dmet = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) Dmet[[10*(lot-1)+sim]])
})
rm(input_path_met)

if (!file.exists("../2210_1SB_featselec/feat_selec/hvg.rds")) {
  print("pb")
  #Dmet_hv <- lapply(Dmet, function(x)
  #  pblapply(x, function(y) findRefinx(y, nmarker = 3e4)))
  #saveRDS(Dmet_hv,"../2210_1SB_featselec/feat_selec/hvg.rds")
} else {Dmet_hv <- readRDS("../2210_1SB_featselec/feat_selec/hvg.rds"); Dmet_hv <- lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) Dmet_hv[[10*(lot-1)+sim]])
})}
Dmet <- lapply(seq(n_lot), function(x)
  pblapply(seq(n_sims), function(y) Dmet[[x]][[y]][Dmet_hv[[x]][[y]],]))

if (!file.exists("../2210_1SB_featselec/feat_selec/TOAST/met.rds")) {
  print("pb")
  #toast_met <- lapply(seq(n_lot), function(lot)
  #  pblapply(seq(n_sims), function(sim)
  #    TOAST::csDeconv(Dmet[[lot]][[sim]], K[[lot]], TotalIter = 10, FUN = deconv_function)$updatedInx))
  #saveRDS(toast_met,"../2210_1SB_featselec/feat_selec/TOAST/met.rds")
} else {toast_met = readRDS("../2210_1SB_featselec/feat_selec/TOAST/met.rds");
toast_met_restricted = readRDS("../2210_1SB_featselec/feat_selec/TOAST/met_restricted.rds")}

## ----
## select CimpleG features
## ----
input_path_met = sort(list.files("../2210_0simu/simulations/met/", pattern = "T_met_ref"))

Tref = lapply(input_path_met, function(x) as.data.frame(readRDS(paste0("../2210_0simu/simulations/met/",x))))
rm(input_path_met)
train_targets <- lapply(Tref, function(x) {
  matrix=diag(ncol(x))
  colnames(matrix)=colnames(x)
  rownames(matrix)=colnames(x)
  as.data.frame(matrix)})
Tref <- lapply(seq(n_lot), function(x)
  lapply(seq(n_sims), function(y)
  Tref[[x]][Dmet_hv[[x]][[y]],]))
rm(Dmet_hv)

if (!file.exists("../2210_1SB_featselec/feat_selec/CimpleG/met.rds")) {
  print("pb")
  #cimpleg_met <- lapply(seq(2,n_lot), function(lot)
  #  lapply(seq(n_sims), function(sim)
  #    CimpleG(t(Tref[[lot]][[sim]]), train_targets = train_targets[[lot]], target_columns = colnames(train_targets[[lot]]), train_only = T)$signatures))
  #saveRDS(cimpleg_met,"../2210_1SB_featselec/feat_selec/CimpleG/met.rds")
} else {cimpleg_met = readRDS("../2210_1SB_featselec/feat_selec/CimpleG/met.rds"); cimpleg_met <- lapply(seq(n_lot-1), function(lot) {
  lapply(seq(n_sims), function(sim) cimpleg_met[[10*(lot-1)+sim]])
});
cimpleg_met_restricted = readRDS("../2210_1SB_featselec/feat_selec/CimpleG/met_restricted.rds"); cimpleg_met_restricted <- lapply(seq(n_lot-1), function(lot) {
  lapply(seq(n_sims), function(sim) cimpleg_met_restricted[[10*(lot-1)+sim]])
})}

## ----
## load D RNA & select features
## ----
input_path_rna = sort(list.files("../2210_0simu/simulations/rna/", pattern = "sim"))

Drna = lapply(input_path_rna, function(x) readRDS(paste0("../2210_0simu/simulations/rna/",x))$D_rna_sim)
Drna=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) Drna[[10*(lot-1)+sim]])
})
rm(input_path_rna)

if (!file.exists("../2210_1SB_featselec/feat_selec/TOAST/rna.rds")) {
  toast_rna <- lapply(seq(n_lot), function(lot)
    pblapply(seq(n_sims), function(sim)
      TOAST::csDeconv(Drna[[lot]][[sim]], K[[lot]], TotalIter = 10, FUN = deconv_function)$updatedInx))
  saveRDS(toast_rna,"../2210_1SB_featselec/feat_selec/TOAST/rna.rds")
} else {toast_rna = readRDS("../2210_1SB_featselec/feat_selec/TOAST/rna.rds");
toast_rna_restricted = readRDS("../2210_1SB_featselec/feat_selec/TOAST/rna_restricted.rds")}


## ----
## Concatenate new datasets based on toast features
## ----
Dtoast <- lapply(seq(n_lot), function(lot)
  pblapply(seq(n_sims), function(sim) {
    mat1 <- Dmet[[lot]][[sim]]
    mat2 <- Drna[[lot]][[sim]]
    rownames(mat1) <- paste0("met_",rownames(mat1))
    rownames(mat2) <- paste0("rna_",rownames(mat2))
    mat1 <- mat1[toast_met[[lot]][[sim]],]
    mat2 <- mat2[toast_rna[[lot]][[sim]],]
    mat <- rbind(mat1,mat2)
    mat}))
saveRDS(Dtoast, paste0("3b_featureselection_step1/TOAST/",date,"_D_toast.rds"))

Dtoast_restricted <- lapply(seq(n_lot), function(lot)
  pblapply(seq(n_sims), function(sim) {
    mat1 <- Dmet[[lot]][[sim]]
    mat2 <- Drna[[lot]][[sim]]
    rownames(mat1) <- paste0("met_",rownames(mat1))
    rownames(mat2) <- paste0("rna_",rownames(mat2))
    mat1 <- mat1[toast_met_restricted[[lot]][[sim]],]
    mat2 <- mat2[toast_rna_restricted[[lot]][[sim]],]
    mat <- rbind(mat1,mat2)
    mat}))
saveRDS(Dtoast_restricted, paste0("3b_featureselection_step1/TOAST/",date,"_D_toast_restricted.rds"))

Dcpg <- c(NA,lapply(seq(2,n_lot), function(lot)
               pblapply(seq(n_sims), function(sim) {
                 mat1 <- Dmet[[lot]][[sim]]
                 mat2 <- Drna[[lot]][[sim]]
                 rownames(mat1) <- paste0("met_",rownames(mat1))
                 rownames(mat2) <- paste0("rna_",rownames(mat2))
                 mat1 <- mat1[paste0("met_",cimpleg_met[[lot-1]][[sim]]),]
                 mat2 <- mat2[toast_rna[[lot]][[sim]],]
                 mat <- rbind(mat1,mat2)
                 mat})))
saveRDS(Dcpg, paste0("3b_featureselection_step1/CimpleG/",date,"_D_cimpleg.rds"))

Dcpg_restricted <- c(NA,lapply(seq(2,n_lot), function(lot)
  pblapply(seq(n_sims), function(sim) {
    mat1 <- Dmet[[lot]][[sim]]
    mat2 <- Drna[[lot]][[sim]]
    rownames(mat1) <- paste0("met_",seq(nrow(mat1)))
    rownames(mat2) <- paste0("rna_",seq(nrow(mat2)))
    mat1 <- mat1[paste0("met_",cimpleg_met_restricted[[lot-1]][[sim]]),]
    mat2 <- mat2[toast_rna_restricted[[lot]][[sim]],]
    mat <- rbind(mat1,mat2)
    mat})))
saveRDS(Dcpg_restricted, paste0("3b_featureselection_step1/CimpleG/",date,"_D_cimpleg_restricted.rds"))
