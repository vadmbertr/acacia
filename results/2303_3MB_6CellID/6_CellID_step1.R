library(CellID)

## ----
## Set parameters
## ----


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

Trna = lapply(sort(list.files("../2210_0simu/simulations/rna/", pattern = "T_rna_ref")), function(x) readRDS(paste0("../2210_0simu/simulations/rna/",x)))
Drna = pblapply(input_path_rna, function(x) readRDS(paste0("../2210_0simu/simulations/rna/",x))$D_rna_sim)
Drna = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) Drna[[10*(lot-1)+sim]])
})

