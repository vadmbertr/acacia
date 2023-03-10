# So far, I have 13 settings to evaluate

## SB (8 settings)
# - Regular deconvolution on RNA (RNA None)
# - Regular deconvolution on MET (MET None)
# - HVG & TOAST on RNA (restrictedT and TOAST)
# - HVG & TOAST on MET (restrictedT and TOAST)
# - HVG & CimpleG on MET (restrictedC and CimpleG)

## MB (5 settings)
# - Mean of the regular deconvolution (RNA + MET None) with 1 method
# - Mean of the regular deconvolution (RNA + MET None) with 3 methods
# - Concatenation with normalization with None
# - Concatenation with normalization with restrictedT/TOAST
# - MultiMAP on latent RNA/MET/concatenated

## ----
## Set parameters
## ----
library(dplyr)
library(ggplot2)
library(see)
source("../../src/ranking_ranking_procedure_functions.R")
ranking_norm=ranking12 # to vary, either ranking12 or ranking_2
ranking_type=ranking_consensus # to vary, either ranking12345 or ranking12453 or ranking_topsis or ranking_consensus
rank_norm=ifelse(all.equal(ranking_norm,ranking12)==T,"normDataset","normGlobal")
rank_type=ifelse(all.equal(ranking_type,ranking_topsis)==T,"topsis",
                 ifelse(all.equal(ranking_type,ranking12345)==T,"order345",
                        ifelse(all.equal(ranking_type,ranking_consensus)==T,"consensus","order453")))
## ----
## load scores SB
## ----
sb1_1 <- readRDS("../2210_1SB/perf_scores/221101_scores_rna_sup.rds")
sb1_2 <- readRDS("../2210_1SB/perf_scores/221101_time_rna_sup.rds")

sb2_1 <- readRDS("../2210_1SB/perf_scores/221101_scores_met_sup.rds")
sb2_2 <- readRDS("../2210_1SB/perf_scores/221101_time_met_sup.rds")

sb3_1 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_rna_toast_sup.rds")
sb3_2 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_rna_toast_sup.rds")

sb4_1 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_met_toast_sup.rds")
sb4_2 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_met_toast_sup.rds")

sb5_1 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_met_cimpleg_sup.rds")
sb5_2 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_met_cimpleg_sup.rds")

## ----
## load scores MB
## ----
n_meth=1
source("../2210_3MB_1mean/1_compare_.R")
mb1_1 <- score_tot %>% filter(candidate=="MB MB") %>% mutate(candidate="MB mean mean1")
mb1_2 <- right_join(sb2_2 %>% select(sim,dataset,deconv,values) %>% rename(candidate=deconv) %>%
                      mutate(sim=as.integer(sim)),
                    right_join(sb1_2 %>% select(sim,dataset,deconv,values) %>%
                                 rename(candidate=deconv) %>%
                                 mutate(sim=as.integer(sim)),
                               score_tot %>%
                                 mutate(candidate=sapply(candidate, function(x)
                                   strsplit(x," ")[[1]][1]),
                                   dup = paste(sim,dataset,candidate)) %>%
                                 filter(!duplicated(dup)) %>%
                                 select(sim,dataset,candidate),
                               by=c("sim","dataset","candidate")),
                    by=c("sim","dataset","candidate")) %>%
  mutate(values.x=ifelse(is.na(values.x),0,values.x),
         values.y=ifelse(is.na(values.y),0,values.y),
         values=pmax(values.y,values.x)) %>%
  group_by(dataset,sim) %>% summarise(values=max(values)) %>% mutate(score="time",
                                                                     candidate="MB mean mean1") %>% ungroup() %>% select(values,score,sim,dataset,candidate)
rm(score_tot,n_meth)

n_meth=3
source("../2210_3MB_1mean/1_compare_.R")
mb2_1 <- score_tot %>% filter(candidate=="MB MB") %>% mutate(candidate="MB mean mean3")
mb2_2 <- right_join(sb2_2 %>% select(sim,dataset,deconv,values) %>% rename(candidate=deconv) %>%
                      mutate(sim=as.integer(sim)),
                    right_join(sb1_2 %>% select(sim,dataset,deconv,values) %>%
                                 rename(candidate=deconv) %>%
                                 mutate(sim=as.integer(sim)),
                               score_tot %>%
                                 mutate(candidate=sapply(candidate, function(x)
                                   strsplit(x," ")[[1]][1]),
                                   dup = paste(sim,dataset,candidate)) %>%
                                 filter(!duplicated(dup)) %>%
                                 select(sim,dataset,candidate),
                               by=c("sim","dataset","candidate")),
                    by=c("sim","dataset","candidate")) %>%
  mutate(values.x=ifelse(is.na(values.x),0,values.x),
         values.y=ifelse(is.na(values.y),0,values.y),
         values=pmax(values.y,values.x)) %>%
  group_by(dataset,sim) %>% summarise(values=max(values)) %>% mutate(score="time",
                                                                     candidate="MB mean mean3") %>% ungroup() %>% select(values,score,sim,dataset,candidate)
rm(score_tot,n_meth)

source("../2210_3MB_3concat/3a_rawconcatenation_step3_norm_.R")
mb3_1 <- score_all
mb3_2 <- time_all
rm(score_all,time_all)

source("../2210_3MB_3concat/3b_featureselection_step3_.R")
mb4_1 <- score_all
mb4_2 <- time_all
rm(score_all,time_all)

source("../2302_3MB_4multimap/4d_eval_.R")
mb5_1 <- df_tot1
mb5_2 <- df_tot2
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

mb1_1$values <- as.double(mb1_1$values)
mb1_2$values <- as.double(mb1_2$values)
mb2_2$values <- as.double(mb2_2$values)
mb3_1 <- mb3_1 %>% filter(candidate %in% c("rawconcat_hvg rlr","rawconcat_hvg rlr")) %>%
  mutate(candidate=paste("MB concatHVG",sapply(candidate,function(x) strsplit(x," ")[[1]][2])))
mb3_2 <- mb3_2 %>% filter(candidate %in% c("rawconcat_hvg rlr","rawconcat_hvg rlr")) %>%
  mutate(candidate=paste("MB concatHVG",sapply(candidate,function(x) strsplit(x," ")[[1]][2])))
mb4_1 <- mb4_1 %>% filter(candidate %in% c("TOAST rlr","TOAST_r rlr","TOAST cibersort","TOAST_r cibersort")) %>%
  mutate(candidate=paste0("MB concat",candidate))
mb4_2 <- mb4_2 %>% filter(candidate %in% c("TOAST rlr","TOAST_r rlr","TOAST cibersort","TOAST_r cibersort")) %>%
  mutate(candidate=paste0("MB concat",candidate))
mb5_1 <- mb5_1 %>% slice(grep("multimap",candidate)) %>% slice(-c(grep("ica",candidate), grep("MeDeCom",candidate))) %>%
  mutate(candidate=sapply(candidate, function(x)
    switch(x, "cibersort rna multimap"="MB MultiMAPr cibersort","rlr rna multimap"="MB MultiMAPr rlr",
              "cibersort met multimap"="MB MultiMAPm cibersort","rlr met multimap"="MB MultiMAPm rlr",
              "cibersort concat multimap"="MB MultiMAPc cibersort","rlr concat multimap"="MB MultiMAPc rlr")))
mb5_2 <- mb5_2 %>% slice(grep("multimap",candidate)) %>% slice(-c(grep("ica",candidate), grep("MeDeCom",candidate))) %>% mutate(values=as.double(values), sim=as.integer(sim)) %>%
  mutate(candidate=sapply(candidate, function(x)
    switch(x, "cibersort rna multimap"="MB MultiMAPr cibersort","rlr rna multimap"="MB MultiMAPr rlr",
           "cibersort met multimap"="MB MultiMAPm cibersort","rlr met multimap"="MB MultiMAPm rlr",
           "cibersort concat multimap"="MB MultiMAPc cibersort","rlr concat multimap"="MB MultiMAPc rlr")))

df_tot1 = bind_rows(sb1_1,sb2_1,sb3_1,sb4_1,sb5_1,
                    mb1_1,mb2_1,mb3_1,mb4_1,mb5_1)
df_tot2 = bind_rows(sb1_2,sb2_2,sb3_2,sb4_2,sb5_2,
                    mb1_2,mb2_2,mb3_2,mb4_2,mb5_2)
rm(sb1_1,sb2_1,sb3_1,sb4_1,sb5_1,mb1_1,mb2_1,mb3_1,mb4_1,mb5_1,
   sb1_2,sb2_2,sb3_2,sb4_2,sb5_2,mb1_2,mb2_2,mb3_2,mb4_2,mb5_2)

## ----
## Rank
## ----
source("../../src/ranking_ranking_procedure_functions.R")
rank <- ranking_type(df_tot1,df_tot2,ranking=ranking_norm)
rank$blocktype <- sapply(rank$candidate, function(x) strsplit(x, " ")[[1]][1])
rank$datatype <- sapply(rank$candidate, function(x) strsplit(x, " ")[[1]][2])
rank$deconv <- sapply(rank$candidate, function(x) strsplit(x, " ")[[1]][3])
rank <- rank %>% filter(!(deconv %in% c("CIBERSORT4","svr4")))
rank$deconv[rank$deconv %in% c("CIBERSORT","CIBERSORT3")] <- "cibersort"
rank$deconv[rank$deconv %in% c("RLR")] <- "rlr"
rank$deconv[rank$deconv %in% c("svr3")] <- "svr"
rank$datatype <- factor(rank$datatype, levels=c("mean","concatHVG","concatTOAST_r","concatTOAST","MultiMAPr","MultiMAPm","MultiMAPc",
                                                "RNANone","METNone","METNone_restrictedC","METCimpleG","RNANone_restrictedT","RNATOAST","METNone_restrictedT","METTOAST"))

ggplot(rank[!(rank$deconv%in%c("elastic_net","nnls","ols","svr")),], aes(x=datatype, y=values, color=deconv)) +
  geom_point() +
  geom_line(aes(group=deconv)) +
  geom_line(aes(group=deconv), linewidth=5, alpha=.4) +
  facet_wrap(~blocktype, scales = "free_x") +
  scale_color_viridis_d(option="plasma") +
  geom_hline(yintercept = max(rank$values, na.rm=T), linetype='dashed') +
  ylab("Aggregated score") +
  xlab("") +
  theme_modern(axis.text.angle = 45)
ggsave(paste0("1_eval_sup/221101_plot_",rank_type,"_",rank_norm,".pdf"),
       width=10, height=6)
