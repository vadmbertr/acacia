## ----
## Set parameters
## ----
library(dplyr)
source("../../src/plot_ind_scores.R")
block="met" # to vary, rna or met
feat_selec='cimpleg' # to vary, cimpleg or toast

## ----
## load scores
## ----
scores1_sup <- readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_",feat_selec,"_sup.rds"))
scores2_sup <- readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_",feat_selec,"_sup.rds"))
scores1_unsup <- readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_",feat_selec,"_unsup.rds"))
scores2_unsup <- readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_",feat_selec,"_unsup.rds"))

scores1_sup <- scores1_sup[!(sapply(scores1_sup$setting, function(x) length(grep("perf_c",x))==1)),]
scores1_sup <- scores1_sup[!(sapply(scores1_sup$setting, function(x) length(grep("perf_s",x))==1)),]
scores1_sup <- scores1_sup[!(scores1_sup$score!="pearson" & sapply(scores1_sup$setting, function(x) length(grep("med_c",x))==1)),]
scores1_sup <- scores1_sup[!(scores1_sup$score!="pearson" & sapply(scores1_sup$setting, function(x) length(grep("med_s",x))==1)),]

scores1_unsup <- scores1_unsup[!(sapply(scores1_unsup$setting, function(x) length(grep("perf_c",x))==1)),]
scores1_unsup <- scores1_unsup[!(sapply(scores1_unsup$setting, function(x) length(grep("perf_s",x))==1)),]
scores1_unsup <- scores1_unsup[!(scores1_unsup$score!="pearson" & sapply(scores1_unsup$setting, function(x) length(grep("med_c",x))==1)),]
scores1_unsup <- scores1_unsup[!(scores1_unsup$score!="pearson" & sapply(scores1_unsup$setting, function(x) length(grep("med_s",x))==1)),]

scores2_sup$values <- as.numeric(scores2_sup$values)
scores2_unsup$values <- as.numeric(scores2_unsup$values)
scores2_sup$sim <- as.numeric(scores2_sup$sim)
scores2_unsup$sim <- as.numeric(scores2_unsup$sim)

## ----
## Reorder data
## ----
scores_sup <- bind_rows(scores2_sup,
                        scores1_sup %>%
                          rename(scor=score) %>%
                          mutate(score=paste(scor,setting)) %>%
                          select(values,score,deconv,sim,dataset,feat_selec))
scores_unsup <- bind_rows(scores2_unsup,
                          scores1_unsup %>%
                            rename(scor=score) %>%
                            mutate(score=paste(scor,setting)) %>%
                            select(values,score,deconv,sim,dataset,feat_selec))
rm(scores1_sup,scores1_unsup,scores2_sup,scores2_unsup)

## ----
## Plots
## ----
name_score=c("time","rmse perf_g","mae perf_g","pearson perf_g",
             "pearson med_c","pearson med_s","pearson sd_c","pearson sd_s")
scale_y_log=c(T,T,T,F,F,F,T,T)
scores_sup$facet <- scores_sup$dataset
scores_unsup$facet <- scores_unsup$dataset

mapply(function(score,scale) {
  if (scale) {p=plot_ind_scores(scores_sup, score, T) +
      scale_y_log10()}
  else {p=plot_ind_scores(scores_sup, score, T)}
  p
  ggsave(paste0("1_compare_ind_scores/",block,"_",gsub(" ","X",score),"_",feat_selec,"_sup.pdf"),
         width=8, height=5)
}, score=name_score, scale=scale_y_log)
mapply(function(score,scale) {
  if (scale) {p=plot_ind_scores(scores_unsup, score, T) +
    scale_y_log10()}
  else {p=plot_ind_scores(scores_unsup, score, T)}
  p
  ggsave(paste0("1_compare_ind_scores/",block,"_",gsub(" ","X",score),"_",feat_selec,"_unsup.pdf"),
         width=8, height=5)
}, score=name_score, scale=scale_y_log)
