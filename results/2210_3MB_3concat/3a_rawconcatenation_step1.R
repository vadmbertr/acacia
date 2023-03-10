library(pbapply)

## ----
## Set parameters
## ----
n_hv=5e3

## ----
## load D MET&RNA
## ----
input_path_met = sort(list.files("../2210_0simu/simulations/met/", pattern = "sim"))
input_path_rna = sort(list.files("../2210_0simu/simulations/rna/", pattern = "sim"))
input_path_nlot = sort(list.files("../2210_0simu/simulations/rna/", pattern = "T_rna_ref"))
n_lot = length(input_path_nlot)
n_sims = length(input_path_rna)/n_lot
rm(input_path_nlot)

Tmet = lapply(sort(list.files("../2210_0simu/simulations/met/", pattern = "T_met_ref")), function(x) readRDS(paste0("../2210_0simu/simulations/met/",x)))
Dmet = pblapply(input_path_met, function(x) readRDS(paste0("../2210_0simu/simulations/met/",x))$D_met_sim)
Dmet = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) Dmet[[10*(lot-1)+sim]])
})
hv_met <- lapply(Dmet, function(x)
  pblapply(x, function(y) TOAST::findRefinx(y, nmarker = n_hv)))
Dmet_hv <- lapply(seq(n_lot), function(lot)
  pblapply(seq(n_sims), function(sim) Dmet[[lot]][[sim]][hv_met[[lot]][[sim]],]))
Tmet_hv <- lapply(seq(n_lot), function(lot)
  pblapply(seq(n_sims), function(sim) Tmet[[lot]][hv_met[[lot]][[sim]],]))
rm(Dmet,Tmet,hv_met)

Trna = lapply(sort(list.files("../2210_0simu/simulations/rna/", pattern = "T_rna_ref")), function(x) readRDS(paste0("../2210_0simu/simulations/rna/",x)))
Drna = pblapply(input_path_rna, function(x) readRDS(paste0("../2210_0simu/simulations/rna/",x))$D_rna_sim)
Drna = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) Drna[[10*(lot-1)+sim]])
})
hv_rna <- lapply(Drna, function(x)
  pblapply(x, function(y) TOAST::findRefinx(y, nmarker = n_hv)))
Drna_hv <- lapply(seq(n_lot), function(lot)
  pblapply(seq(n_sims), function(sim) Drna[[lot]][[sim]][hv_rna[[lot]][[sim]],]))
Trna_hv <- lapply(seq(n_lot), function(lot)
  pblapply(seq(n_sims), function(sim) Trna[[lot]][hv_rna[[lot]][[sim]],]))
rm(Drna,Trna,hv_rna)

## ----
## Concatenate D matrices with hvg/hvp (& scale)
## ----
Dmet_hv_scale <- lapply(Dmet_hv, function(Dmet_lot)
  pblapply(Dmet_lot, function(Dmet_lot_sim)
    t(apply(Dmet_lot_sim, 1, function(x) {
      scl <- x-min(x)
      if (max(scl)>0) scl <- scl/max(scl)
      return(scl)}))))
Tmet_hv_scale <- lapply(Tmet_hv, function(Tmet_lot)
  pblapply(Tmet_lot, function(Tmet_lot_sim)
    t(apply(Tmet_lot_sim, 1, function(x) {
      scl <- x-min(x)
      if (max(scl)>0) scl <- scl/max(scl)
      return(scl)}))))
Drna_hv_scale <- lapply(Drna_hv, function(Drna_lot)
  pblapply(Drna_lot, function(Drna_lot_sim)
    t(apply(Drna_lot_sim, 1, function(x) {
      scl <- x-min(x)
      if (max(scl)>0) scl <- scl/max(scl)
      return(scl)}))))
Trna_hv_scale <- lapply(Trna_hv, function(Trna_lot)
  pblapply(Trna_lot, function(Trna_lot_sim)
    t(apply(Trna_lot_sim, 1, function(x) {
      scl <- x-min(x)
      if (max(scl)>0) scl <- scl/max(scl)
      return(scl)}))))

D_hv <- lapply(seq(n_lot), function(lot)
  lapply(seq(n_sims), function(sim)
  rbind(Dmet_hv[[lot]][[sim]],Drna_hv[[lot]][[sim]])))
T_hv <- lapply(seq(n_lot), function(lot)
  lapply(seq(n_sims), function(sim)
    rbind(Tmet_hv[[lot]][[sim]],Trna_hv[[lot]][[sim]])))
D_hv_scale <- lapply(seq(n_lot), function(lot)
  lapply(seq(n_sims), function(sim)
    rbind(Dmet_hv_scale[[lot]][[sim]],Drna_hv_scale[[lot]][[sim]])))
T_hv_scale <- lapply(seq(n_lot), function(lot)
  lapply(seq(n_sims), function(sim)
    rbind(Tmet_hv_scale[[lot]][[sim]],Trna_hv_scale[[lot]][[sim]])))
rm(Dmet_hv,Drna_hv,Tmet_hv,Trna_hv,Dmet_hv_scale,Drna_hv_scale,Tmet_hv_scale,Trna_hv_scale)

saveRDS(list(D=D_hv,Tref=T_hv),paste0("3a_rawconcatenation_step1/",
                    strsplit(input_path_rna[1],"_")[[1]][1],
                    "_DT_hv.rds"))
saveRDS(list(D=D_hv_scale,Tref=T_hv_scale),paste0("3a_rawconcatenation_step1/",
                          strsplit(input_path_rna[1],"_")[[1]][1],
                          "_DT_hv_scale.rds"))
