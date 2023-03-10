## ----
## Data
## ----
lot1 = "dBREAST"
lot2 = "dPANCREAS"
lot3 = "lot1"
lots = c(lot1,lot2,lot3)
n_rep = 10
today = strsplit(list.files("simulations/rna/")[1],"_")[[1]][1]
data_D <- lapply(seq(length(lots)), function(x) {
  lapply(seq(n_rep), function(i) {
    sim_txt <- ifelse(length(strsplit(as.character(i),"")[[1]])==1,'sim0','sim')
    readRDS(paste0("simulations/rna/",today,"_",lots[x],"_",sim_txt,i,".rds"))$D_rna_sim
  })
})
data_T <- lapply(seq(length(lots)), function(x) 
  readRDS(paste0("simulations/rna/",today,"_",lots[x],"_T_rna_ref.rds"))
)

## ----
## Sparsity
## ----
spars_D <- sapply(data_D, function(lot)
  sapply(lot, function(sim)
    round(100*mean(sim==0))))
spars_T <- sapply(data_T, function(lot) round(100*mean(lot==0)))
