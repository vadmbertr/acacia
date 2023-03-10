## ----
## Set parameters
## ----
library(dplyr)
library(ggplot2)
library(see)
source("../../src/ranking_ranking_procedure_functions.R")
ranking_norm=ranking_2 # to vary, either ranking12 or ranking_2
ranking_type=ranking12345 # to vary, either ranking12345 or ranking12453
rank_norm=ifelse(all.equal(ranking_norm,ranking12)==T,"normDataset","normGlobal")
rank_type=ifelse(all.equal(ranking_type,ranking12345)==T,"order345","order453")

## ----
## load scores
## ----
df_score1 <- readRDS("4c_deconv/perf_scores/221101_scores.rds")
df_score2 <- readRDS("4c_deconv/perf_scores/221101_time.rds")

df_bl111 <- readRDS("../2210_1SB/perf_scores/221101_scores_met_sup.rds")
df_bl111$type = 'met'
df_bl121 <- readRDS("../2210_1SB/perf_scores/221101_scores_rna_sup.rds")
df_bl121$type = 'rna'
df_bl112 <- readRDS("../2210_1SB/perf_scores/221101_scores_met_unsup.rds")
df_bl112$type = 'met'
df_bl122 <- readRDS("../2210_1SB/perf_scores/221101_scores_rna_unsup.rds")
df_bl122$type = 'rna'
df_bl1 <- rbind(df_bl111,df_bl121,df_bl112,df_bl122)
rm(df_bl111,df_bl121,df_bl112,df_bl122)
df_bl211 <- readRDS("../2210_1SB/perf_scores/221101_time_met_sup.rds")
df_bl211$type = 'met'
df_bl221 <- readRDS("../2210_1SB/perf_scores/221101_time_rna_sup.rds")
df_bl221$type = 'rna'
df_bl212 <- readRDS("../2210_1SB/perf_scores/221101_time_met_unsup.rds")
df_bl212$type = 'met'
df_bl222 <- readRDS("../2210_1SB/perf_scores/221101_time_rna_unsup.rds")
df_bl222$type = 'rna'
df_bl2 <- rbind(df_bl211,df_bl221,df_bl212,df_bl222)
rm(df_bl211,df_bl221,df_bl212,df_bl222)

df_toast111 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_met_toast_sup.rds")
df_toast111$type = 'met'
df_toast121 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_rna_toast_sup.rds")
df_toast121$type = 'rna'
df_toast112 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_met_toast_unsup.rds")
df_toast112$type = 'met'
df_toast122 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_rna_toast_unsup.rds")
df_toast122$type = 'rna'
df_toast1 <- rbind(df_toast111,df_toast121,df_toast112,df_toast122)
rm(df_toast111,df_toast121,df_toast112,df_toast122)
df_toast211 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_met_toast_sup.rds")
df_toast211$type = 'met'
df_toast221 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_rna_toast_sup.rds")
df_toast221$type = 'rna'
df_toast212 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_met_toast_unsup.rds")
df_toast212$type = 'met'
df_toast222 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_rna_toast_unsup.rds")
df_toast222$type = 'rna'
df_toast2 <- rbind(df_toast211,df_toast221,df_toast212,df_toast222)
rm(df_toast211,df_toast221,df_toast212,df_toast222)

df_bl1$deconv[df_bl1$deconv=="CIBERSORT3"] <- "cibersort"
df_bl1$deconv[df_bl1$deconv=="CIBERSORT"] <- "cibersort"
df_bl1 <- df_bl1[df_bl1$deconv!="CIBERSORT4",]
df_bl1 <- df_bl1[df_bl1$deconv!="svr4",]
df_bl1$deconv[df_bl1$deconv=="RLR"] <- "rlr"
df_bl1$deconv[df_bl1$deconv=="DECONica"] <- "ica"
df_bl1$deconv[df_bl1$deconv=="DECONnmf"] <- "nmf"
df_bl1$deconv[df_bl1$deconv=="ICA"] <- "ica"
df_bl1 <- df_bl1[df_bl1$deconv!="svr4",]
df_bl1$deconv[df_bl1$deconv=="svr3"] <- "svr"
df_bl2$deconv[df_bl2$deconv=="CIBERSORT3"] <- "cibersort"
df_bl2$deconv[df_bl2$deconv=="CIBERSORT"] <- "cibersort"
df_bl2 <- df_bl2[df_bl2$deconv!="CIBERSORT4",]
df_bl2 <- df_bl2[df_bl2$deconv!="svr4",]
df_bl2$deconv[df_bl2$deconv=="RLR"] <- "rlr"
df_bl2$deconv[df_bl2$deconv=="DECONica"] <- "ica"
df_bl2$deconv[df_bl2$deconv=="DECONnmf"] <- "nmf"
df_bl2$deconv[df_bl2$deconv=="ICA"] <- "ica"
df_bl2 <- df_bl2[df_bl2$deconv!="svr4",]
df_bl2$deconv[df_bl2$deconv=="svr3"] <- "svr"

df_toast1$deconv[df_toast1$deconv=="CIBERSORT"] <- "cibersort"
df_toast1$deconv[df_toast1$deconv=="RLR"] <- "rlr"
df_toast1$deconv[df_toast1$deconv=="DECONica"] <- "ica"
df_toast1$deconv[df_toast1$deconv=="DECONnmf"] <- "nmf"
df_toast1$deconv[df_toast1$deconv=="ICA"] <- "ica"
df_toast2$deconv[df_toast2$deconv=="CIBERSORT"] <- "cibersort"
df_toast2$deconv[df_toast2$deconv=="RLR"] <- "rlr"
df_toast2$deconv[df_toast2$deconv=="DECONica"] <- "ica"
df_toast2$deconv[df_toast2$deconv=="DECONnmf"] <- "nmf"
df_toast2$deconv[df_toast2$deconv=="ICA"] <- "ica"

colnames(df_score1)[colnames(df_score1)=="score"] <- "scor"
colnames(df_bl1)[colnames(df_bl1)=="score"] <- "scor"
colnames(df_toast1)[colnames(df_toast1)=="score"] <- "scor"
df_score1 <- df_score1 %>%
  mutate(score=paste(scor,setting)) %>%
  select(c("values","score","sim","dataset","deconv","type")) %>%
  mutate(candidate=paste(deconv,type,"multimap")) %>%
  select(c('values','score','sim','dataset','candidate'))
df_bl1 <- df_bl1 %>%
  mutate(score=paste(scor,setting)) %>%
  select(c("values","score","sim","dataset","deconv","type")) %>%
  mutate(candidate=paste(deconv,type,"baseline")) %>%
  select(c('values','score','sim','dataset','candidate'))
df_toast1 <- df_toast1 %>%
  mutate(score=paste(scor,setting)) %>%
  select(c("values","score","sim","dataset","deconv","type")) %>%
  mutate(candidate=paste(deconv,type,"toast")) %>%
  select(c('values','score','sim','dataset','candidate'))
df_score2 <- df_score2 %>%
  mutate(candidate=paste(deconv,type,"multimap")) %>%
  select(c('values','score','sim','dataset','candidate'))
df_bl2 <- df_bl2 %>%
  mutate(candidate=paste(deconv,type,"baseline")) %>%
  select(c('values','score','sim','dataset','candidate'))
df_toast2 <- df_toast2 %>%
  mutate(candidate=paste(deconv,type,"toast")) %>%
  select(c('values','score','sim','dataset','candidate'))

name_data = sort(unique(df_bl1$dataset))
df_score1$dataset <- sapply(df_score1$dataset, function(x)
  switch(x, "dataset1"=name_data[1], "dataset2"=name_data[2], "dataset3"=name_data[3]))
df_score2$dataset <- sapply(df_score2$dataset, function(x)
  switch(x, "dataset1"=name_data[1], "dataset2"=name_data[2], "dataset3"=name_data[3]))

df_tot1 <- rbind(df_bl1,df_toast1,df_score1)
df_tot2 <- rbind(df_bl2,df_toast2[,colnames(df_bl2)],df_score2)
rm(df_bl1,df_toast1,df_score1,
   df_bl2,df_toast2,df_score2)

## ----
## Rank
## ----
rank <- ranking_type(df_tot1,df_tot2,ranking=ranking_norm)

rank$deconv <- sapply(rank$candidate, function(x)
  strsplit(x, " ")[[1]][1])
rank$type <- sapply(rank$candidate, function(x)
  strsplit(x, " ")[[1]][2])
rank$setting <- sapply(rank$candidate, function(x)
  strsplit(x, " ")[[1]][3])
rank$setting <- factor(rank$setting, levels = c("baseline","toast","multimap"))
rank$type <- factor(rank$type, levels = c("met","rna","concat"))

rank_bis = rank[!(rank$deconv%in%c("nnls","svr","WISP","RefFreeEWAS","ols","nmf","elastic_net","EDec")),]
multimap_deconv = unique(rank_bis$deconv[rank_bis$setting=="multimap"])
multimap_order=sapply(multimap_deconv, function(x)
  order(unique(rank_bis$deconv))[unique(rank_bis$deconv)==x])
multimap_color = sapply(multimap_order, function(x)
  palette_social(palette='complement')(n=length(unique(rank_bis$deconv)))[x])
ggplot(rank_bis, aes(x=setting, y=values, color=deconv)) +
  geom_point() +
  geom_line(aes(group=deconv)) +
  geom_line(aes(group=deconv), linewidth=5, alpha=.4) +
  scale_color_social_d(palette = "complement") +
  facet_wrap(~type, scales = "free_x") +
  geom_hline(yintercept=rank_bis$values[rank_bis$type=="concat" & rank_bis$deconv==multimap_deconv[1]], color=multimap_color[1], linetype = 'dashed') +
  geom_hline(yintercept=rank_bis$values[rank_bis$type=="concat" & rank_bis$deconv==multimap_deconv[2]], color=multimap_color[2], linetype = 'dashed') +
  geom_hline(yintercept=rank_bis$values[rank_bis$type=="concat" & rank_bis$deconv==multimap_deconv[3]], color=multimap_color[3], linetype = 'dashed') +
  geom_hline(yintercept=rank_bis$values[rank_bis$type=="concat" & rank_bis$deconv==multimap_deconv[4]], color=multimap_color[4], linetype = 'dashed') +
  geom_hline(yintercept = max(rank$values), color="black", linetype="dashed", linewidth=1, alpha=.5) +
  ylab("Aggregated score") +
  xlab("") +
  theme_modern(axis.text.angle = 45)
ggsave(paste0("4d_eval/221101_plot_",rank_type,"_",rank_norm,".pdf"),
       width=7, height=4)
