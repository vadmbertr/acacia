## ----
## Set parameters
## ----
library(ggplot2)
library(see)
library(dplyr)
block2='met'
block1='rna'
pattern2=paste0("_T_",block2)
source("../../src/score_functions.R")
setting='RunsuptoM'

## ----
## load Apred_bl, Apred_to, A
## ----
input_path = paste0("../2210_0simu/simulations/",block2,"/")
input_path_matrix = sort(list.files(input_path, pattern = "sim"))
Atrue = lapply(input_path_matrix, function(x) readRDS(paste0(input_path,x))$A_ref)
rm(input_path_matrix)

input_path_T = sort(list.files(input_path, pattern = pattern2))
n_lot = length(input_path_T)
input_path_nsims = sort(list.files(input_path, pattern = "sim"))
n_sims = length(input_path_nsims)/n_lot
name_data = unname(sapply(input_path_T,function(x) rev(strsplit(strsplit(x,pattern2)[[1]][1],"_")[[1]])[1]))
Atrue = lapply(seq(n_lot), function(lot)
  lapply(seq(n_sims), function(sim)
    Atrue[[(lot-1)*10+sim]]))
names(Atrue) = name_data
rm(input_path_nsims,input_path)

input_path_methods = paste0("2a_init_",setting,"/")
input_path_methods = sort(list.files(input_path_methods))
deconv_methods_2 = unique(sapply(input_path_methods, function(x)
  strsplit(strsplit(x,"to_")[[1]][2],"_sim")[[1]][1]))
deconv_methods_1 = unique(sapply(input_path_methods, function(x)
  strsplit(strsplit(x,"Apred_")[[1]][2],"_to")[[1]][1]))
rm(input_path_methods)

n_cell_types = sapply(name_data, function(x) nrow(Atrue[[x]][[1]]))
names(n_cell_types) = name_data
n_samples = ncol(Atrue[[1]][[1]])

Apred_bl = list()
Apred_to = list()
for (lot in seq(n_lot)) {
  Apred_bl[[name_data[lot]]] = list()
  Apred_to[[name_data[lot]]] = list()
  for (sim in seq(n_sims)) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    Apred_bl[[name_data[lot]]][[sim]] = list()
    Apred_to[[name_data[lot]]][[sim]] = list()
    for (methM in deconv_methods_2) {
      Apred_bl[[name_data[lot]]][[sim]][[methM]] = list()
      Apred_to[[name_data[lot]]][[sim]][[methM]] = list()
      for (methR in deconv_methods_1) {
        if (file.exists(paste0("2a_init_",setting,"/",
                               strsplit(input_path_T[lot],pattern2)[[1]][1],
                               "_Apred_",
                               methR,
                               "_to_",
                               methM,
                               "_",sim_txt,
                               sim,".rds"))) {
          A_tot=readRDS(paste0("2a_init_",setting,"/",
                               strsplit(input_path_T[lot],pattern2)[[1]][1],
                               "_Apred_",
                               methR,
                               "_to_",
                               methM,
                               "_",sim_txt,
                               sim,".rds"))$res
          Apred_to[[name_data[lot]]][[sim]][[methM]][[methR]] = A_tot
          Apred_bl[[name_data[lot]]][[sim]][[methM]][[methR]] = readRDS(paste0("../2210_1SB/deconv/met/unsup/",
                                                                               strsplit(input_path_T[lot],pattern2)[[1]][1],
                                                                               "_Apred_",
                                                                               methM,
                                                                               "_",sim_txt,
                                                                               sim,".rds"))
          if (all(is.na(A_tot))) {Apred_to[[name_data[lot]]][[sim]][[methM]][[methR]] = matrix(NA,ncol=n_samples,nrow=n_cell_types[lot])}
        }
        else {print(paste0(lot,sim,methM,methR," not run yet"))}
      }
    }
  }
}
rm(A_tot,lot,sim,methR,methM)

df_A <- do.call(rbind, lapply(name_data, function(lot) {
  df_bl = data.frame("Sim"=rep(seq(n_sims), each=length(deconv_methods_2)*length(deconv_methods_1)*n_samples*n_cell_types[lot]),
                     "Initialisation"="No",
                     "From"=rep(rep(deconv_methods_1, each=n_samples*n_cell_types[lot]),length(deconv_methods_2)*n_sims),
                     "To"=rep(rep(deconv_methods_2, each=length(deconv_methods_1)*n_samples*n_cell_types[lot]),n_sims),
                     "A"=unlist(Apred_bl[[lot]]),
                     "Atrue"=unlist(Atrue[[lot]]),
                     "Pos"=seq(n_samples*n_cell_types[lot]))
  df_to = data.frame("Sim"=rep(seq(n_sims), each=length(deconv_methods_2)*length(deconv_methods_1)*n_samples*n_cell_types[lot]),
                     "Initialisation"="Yes",
                     "From"=rep(rep(deconv_methods_1, each=n_samples*n_cell_types[lot]),length(deconv_methods_2)*n_sims),
                     "To"=rep(rep(deconv_methods_2, each=length(deconv_methods_1)*n_samples*n_cell_types[lot]),n_sims),
                     "A"=unlist(Apred_to[[lot]]),
                     "Atrue"=unlist(Atrue[[lot]]),
                     "Pos"=seq(n_samples*n_cell_types[lot]))
  df = rbind(df_bl,df_to)
  df$Dataset = lot
  df}))

## ----
## RMSE
## ----
df_A$Experiment = 1
for (row in seq(2,nrow(df_A))) {
  if (df_A$Pos[row]==1) {
    df_A$Experiment[seq(length(df_A$Experiment))>=row] <- df_A$Experiment[row] + 1
  }
}

## ----
## CN of Apred_in
## ----
get_CN <- function(mat) {
  norm(pracma::pinv(mat),"I")*norm(mat,"I")
}
Apred_in=lapply(seq(n_lot), function(lot) {
  print(paste0("Running dataset n°",lot,"/",n_lot))
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    lapply(deconv_methods_1, function(meth) {
      if (file.exists(paste0("../2210_1SB/deconv/",block1,"/unsup/",
                             strsplit(input_path_T[lot],pattern2)[[1]][1],
                             "_Apred_",
                             meth,
                             "_",sim_txt,
                             sim,".rds"))) {
        A_pred=readRDS(paste0("../2210_1SB/deconv/",block1,"/unsup/",
                              strsplit(input_path_T[lot],pattern2)[[1]][1],
                              "_Apred_",
                              meth,
                              "_",sim_txt,
                              sim,".rds"))
      }
      else {print(paste0(meth," not run yet"))}
    })
  })
})
CN_A1 <- lapply(seq(n_lot), function(lot)
  lapply(seq(n_sims), function(sim)
    sapply(seq_along(deconv_methods_1), function(meth)
      get_CN(Apred_in[[lot]][[sim]][[meth]]))))
CN_A1 <- data.frame(CN=unlist(CN_A1),
                    dataset=rep(name_data,each=n_sims*length(deconv_methods_1)),
                    sim=rep(rep(seq(n_sims),each=length(deconv_methods_1)),n_lot),
                    deconv1=rep(rep(deconv_methods_1,n_sims),n_lot))
ggplot(CN_A1, aes(y=CN, x=deconv1, color=dataset)) +
  geom_boxplot(position=position_dodge(width=0.9)) +
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  scale_color_see_d() +
  theme_modern(axis.text.angle = 45)
ggsave(paste0("2b_eval_",setting,"/",
              strsplit(input_path_T[1],"_")[[1]][1],
              "_CN_A1.pdf"),
       width = 8, height = 6)

df_rmse_CN = df_A %>%
  group_by(Experiment) %>%
  summarise(RMSE=rmse(A,Atrue)) %>%
  mutate(Dataset=df_A$Dataset[!duplicated(df_A$Experiment)]) %>%
  mutate(Sim=df_A$Sim[!duplicated(df_A$Experiment)]) %>%
  mutate(From=df_A$From[!duplicated(df_A$Experiment)]) %>%
  mutate(To=df_A$To[!duplicated(df_A$Experiment)]) %>%
  mutate(Initialisation=df_A$Initialisation[!duplicated(df_A$Experiment)])
df_rmse_CN_bef = df_rmse_CN %>%
  select(Dataset,Sim,From,To,Initialisation,RMSE) %>%
  filter(Initialisation=="No")
df_rmse_CN_aft = df_rmse_CN %>%
  select(Dataset,Sim,From,To,Initialisation,RMSE) %>%
  filter(Initialisation=="Yes")
df_rmse_CN <- inner_join(df_rmse_CN_bef, df_rmse_CN_aft, by = c("Dataset", "Sim", "From", "To")) %>%
  mutate(RMSEdelta=RMSE.y-RMSE.x) %>%
  select(Dataset,Sim,From,To,RMSEdelta)
colnames(CN_A1) <- c("CN","Dataset","Sim","From")
df_rmse_CN <- inner_join(df_rmse_CN, CN_A1, by = c("Dataset","Sim","From"))
ggplot(df_rmse_CN[!is.na(df_rmse_CN$RMSEdelta),], aes(x=CN, y=RMSEdelta, color=Dataset, shape=From)) +
  geom_point() +
  facet_wrap(~To, scales = "free") +
  scale_color_flat_d() +
  theme_modern()
ggsave(paste0("2b_eval_",setting,"/",
              strsplit(input_path_T[1],"_")[[1]][1],
              "_CN_RMSE1.pdf"),
       width = 12, height = 4)
ggplot(df_rmse_CN[!is.na(df_rmse_CN$RMSEdelta),], aes(x=CN, y=RMSEdelta, color=Dataset, shape=To)) +
  geom_point() +
  facet_wrap(~From, scales = "free") +
  scale_color_flat_d() +
  theme_modern()
ggsave(paste0("2b_eval_",setting,"/",
              strsplit(input_path_T[1],"_")[[1]][1],
              "_CN_RMSE2.pdf"),
       width = 8, height = 4)
ggplot(df_rmse_CN[!is.na(df_rmse_CN$RMSEdelta),], aes(x=CN, y=RMSEdelta, color=Dataset)) +
  geom_point() +
  facet_wrap(~paste(From,To), scales = "free") +
  scale_color_flat_d() +
  theme_modern()
ggsave(paste0("2b_eval_",setting,"/",
              strsplit(input_path_T[1],"_")[[1]][1],
              "_CN_RMSE3.pdf"),
       width = 12, height = 6)
