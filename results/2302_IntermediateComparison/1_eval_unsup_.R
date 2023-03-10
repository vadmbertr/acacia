# So far, I have 12 settings to evaluate

## SB (8 settings)
# - Regular deconvolution on RNA (RNA None)
# - Regular deconvolution on MET (MET None)
# - HVG & TOAST on RNA (restrictedT and TOAST)
# - HVG & TOAST on MET (restrictedT and TOAST)
# - HVG & CimpleG on MET (restrictedC and CimpleG)

## MB (4 settings)
# - Concatenation with normalization with None
# - Concatenation with normalization with restrictedT/TOAST
# - Initialization of RNA or MET from MET or RNA results
# - MultiMAP on latent RNA/MET/concatenated

## ----
## Set parameters
## ----
library(dplyr)
library(ggplot2)
library(see)
source("../../src/ranking_ranking_procedure_functions.R")
ranking_norm=ranking12 # to vary, either ranking12 or ranking_2
ranking_type=ranking12345 # to vary, either ranking12345 or ranking12453
rank_norm=ifelse(all.equal(ranking_norm,ranking12)==T,"normDataset","normGlobal")
rank_type=ifelse(all.equal(ranking_type,ranking12345)==T,"order345","order453")

## ----
## load scores SB
## ----
sb1_1 <- readRDS("../2210_1SB/perf_scores/221101_scores_rna_unsup.rds")
sb1_2 <- readRDS("../2210_1SB/perf_scores/221101_time_rna_unsup.rds")

sb2_1 <- readRDS("../2210_1SB/perf_scores/221101_scores_met_unsup.rds")
sb2_2 <- readRDS("../2210_1SB/perf_scores/221101_time_met_unsup.rds")

sb3_1 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_rna_toast_unsup.rds")
sb3_2 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_rna_toast_unsup.rds")

sb4_1 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_met_toast_unsup.rds")
sb4_2 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_met_toast_unsup.rds")

sb5_1 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_met_cimpleg_unsup.rds")
sb5_2 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_met_cimpleg_unsup.rds")

## ----
## load scores MB
## ----
source("../2210_3MB_3concat/3a_rawconcatenation_step3_norm_.R")
mb1_1 <- score_all
mb1_2 <- time_all
rm(score_all,time_all)

source("../2210_3MB_3concat/3b_featureselection_step3_.R")
mb2_1 <- score_all
mb2_2 <- time_all
rm(score_all,time_all)

res = list()
for (init_type in c("MunsuptoR","MsuptoR","RunsuptoM","RsuptoM")) {
  source("../2210_3MB_2initA/3_compare_.R")
  res[[init_type]] = list(score_tot1,score_tot2)
}
mb3_1 <- bind_rows(res$MunsuptoR[[1]] %>% slice(-grep("no ",candidate)) %>% mutate(candidate=paste0("MB initR_",candidate)),
                   res$MsuptoR[[1]] %>% slice(-grep("no ",candidate)) %>% mutate(candidate=paste0("MB initR_",candidate)),
                   res$RunsuptoM[[1]] %>% slice(-grep("no ",candidate)) %>% mutate(candidate=paste0("MB initM_",candidate)),
                   res$RsuptoM[[1]] %>% slice(-grep("no ",candidate)) %>% mutate(candidate=paste0("MB initM_",candidate)))
mb3_2 <- bind_rows(res$MunsuptoR[[2]] %>% slice(-grep("no ",candidate)) %>% mutate(candidate=paste0("MB initR_",candidate)),
                   res$MsuptoR[[2]] %>% slice(-grep("no ",candidate)) %>% mutate(candidate=paste0("MB initR_",candidate)),
                   res$RunsuptoM[[2]] %>% slice(-grep("no ",candidate)) %>% mutate(candidate=paste0("MB initM_",candidate)),
                   res$RsuptoM[[2]] %>% slice(-grep("no ",candidate)) %>% mutate(candidate=paste0("MB initM_",candidate)))
rm(res,score_tot1,score_tot2,init_type)

source("../2302_3MB_4multimap/4d_eval_.R")
mb4_1 <- df_tot1
mb4_2 <- df_tot2
rm(df_tot1,df_tot2)

## ----
## Organize dfs
## ----
sb1_1 <- sb1_1 %>% mutate(score=paste(score,setting)) %>% mutate(candidate=paste("SB RNANone",deconv)) %>%
  select(values,score,sim,dataset,candidate)
sb2_1 <- sb2_1 %>% mutate(score=paste(score,setting)) %>% mutate(candidate=paste("SB METNone",deconv)) %>%
  select(values,score,sim,dataset,candidate)
sb3_1 <- sb3_1 %>% mutate(score=paste(score,setting)) %>% mutate(candidate=paste0("SB RNA",feat_selec," ",deconv)) %>%
  select(values,score,sim,dataset,candidate)
sb4_1 <- sb4_1 %>% mutate(score=paste(score,setting)) %>% mutate(candidate=paste0("SB MET",feat_selec," ",deconv)) %>%
  select(values,score,sim,dataset,candidate)
sb5_1 <- sb5_1 %>% mutate(score=paste(score,setting)) %>% mutate(candidate=paste0("SB MET",feat_selec," ",deconv)) %>%
  select(values,score,sim,dataset,candidate)
sb1_2 <- sb1_2 %>% mutate(candidate=paste("SB RNANone",deconv), values=as.double(values),sim=as.integer(sim)) %>% select(values,score,sim,dataset,candidate)
sb2_2 <- sb2_2 %>% mutate(candidate=paste("SB METNone",deconv), sim=as.integer(sim)) %>% select(values,score,sim,dataset,candidate)
sb3_2 <- sb3_2 %>% mutate(candidate=paste0("SB RNA",feat_selec," ",deconv), sim=as.integer(sim)) %>% select(values,score,sim,dataset,candidate)
sb4_2 <- sb4_2 %>% mutate(candidate=paste0("SB MET",feat_selec," ",deconv), sim=as.integer(sim)) %>% select(values,score,sim,dataset,candidate)
sb5_2 <- sb5_2 %>% mutate(candidate=paste0("SB MET",feat_selec," ",deconv), sim=as.integer(sim)) %>% select(values,score,sim,dataset,candidate)

mb1_1 <- mb1_1 %>% slice(grep("rawconcat_hvg ica",candidate))
mb1_2 <- mb1_2 %>% slice(grep("rawconcat_hvg ica",candidate))
mb2_1 <- mb2_1 %>% slice(grep("ica",candidate)) %>% slice(-grep("met",candidate)) %>% slice(-grep("rna",candidate))
mb2_2 <- mb2_2 %>% slice(grep("ica",candidate)) %>% slice(-grep("met",candidate)) %>% slice(-grep("rna",candidate))

mb4_1 <- mb4_1 %>% slice(grep("multimap",candidate)) %>% slice(-c(grep("cibersort",candidate),grep("houseman",candidate),grep("rlr",candidate),grep("nnls",candidate),grep("WISP",candidate),grep("svr",candidate),grep("ols",candidate),grep("DeconRNASeq",candidate),grep("elastic_net",candidate))) %>%
  mutate(candidate=sapply(candidate, function(x)
    switch(x, "ica rna multimap"="MB MultiMAPr ica","MeDeCom rna multimap"="MB MultiMAPr MeDeCom",
           "ica met multimap"="MB MultiMAPm ica","MeDeCom met multimap"="MB MultiMAPm MeDeCom",
           "ica concat multimap"="MB MultiMAPc ica","MeDeCom concat multimap"="MB MultiMAPc MeDeCom")))
mb4_2 <- mb4_2 %>% slice(grep("multimap",candidate)) %>% slice(-c(grep("cibersort",candidate),grep("houseman",candidate),grep("rlr",candidate),grep("nnls",candidate),grep("WISP",candidate),grep("svr",candidate),grep("ols",candidate),grep("DeconRNASeq",candidate),grep("elastic_net",candidate))) %>%
  mutate(values=as.double(values),sim=as.double(sim),candidate=sapply(candidate, function(x)
    switch(x, "ica rna multimap"="MB MultiMAPr ica","MeDeCom rna multimap"="MB MultiMAPr MeDeCom",
           "ica met multimap"="MB MultiMAPm ica","MeDeCom met multimap"="MB MultiMAPm MeDeCom",
           "ica concat multimap"="MB MultiMAPc ica","MeDeCom concat multimap"="MB MultiMAPc MeDeCom")))

df_tot1 = bind_rows(sb1_1,sb2_1,sb3_1,sb4_1,sb5_1,
                    mb1_1,mb2_1,mb3_1,mb4_1)
df_tot2 = bind_rows(sb1_2,sb2_2,sb3_2,sb4_2,sb5_2,
                    mb1_2,mb2_2,mb3_2,mb4_2)
rm(sb1_1,sb2_1,sb3_1,sb4_1,sb5_1,mb1_1,mb2_1,mb3_1,mb4_1,
   sb1_2,sb2_2,sb3_2,sb4_2,sb5_2,mb1_2,mb2_2,mb3_2,mb4_2)
rm(euclidean_vector_dist,eval_norm_order,eval_norm_order_,ranking_2,ranking_2_notime,ranking_consensus,
   ranking_topsis,ranking12_notime,ranking12345_notime,ranking12453_notime)