####
#### Set parameters
####
source("../../src/ranking_ranking_procedure_functions.R")
block='met' # to vary, either met or rna
feat_selec='cimpleg' # to vary, either toast or cimpleg

####
#### Load data
####
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

####
#### Prepare df
####
colnames(scores1_sup)[colnames(scores1_sup)=="score"] <- "scor"
scores1_sup <- scores1_sup %>% mutate(score=paste(scor,setting)) %>% select(colnames(scores2_sup))
colnames(scores1_unsup)[colnames(scores1_unsup)=="score"] <- "scor"
scores1_unsup <- scores1_unsup %>% mutate(score=paste(scor,setting)) %>% select(colnames(scores2_unsup))

scores1_sup <- scores1_sup %>%
  mutate(candidate=paste(deconv,feat_selec)) %>%
  select(values,score,sim,dataset,candidate)
scores2_sup <- scores2_sup %>%
  mutate(candidate=paste(deconv,feat_selec)) %>%
  select(values,score,sim,dataset,candidate)
scores1_unsup <- scores1_unsup %>%
  mutate(candidate=paste(deconv,feat_selec)) %>%
  select(values,score,sim,dataset,candidate)
scores2_unsup <- scores2_unsup %>%
  mutate(candidate=paste(deconv,feat_selec)) %>%
  select(values,score,sim,dataset,candidate)

####
#### Ranking
####
score_mat_sup <- ranking12(scores1_sup,scores2_sup)
score_mat_unsup <- ranking12(scores1_unsup,scores2_unsup)

order_sup <- all(score_mat_sup %>%
  group_by(name_score,dataset) %>%
  select(name_score,dataset,candidate,val,normval) %>%
  mutate(order_before=order(val), order_after=order(normval)) %>%
  mutate(rank_conserved=order_before==order_after) %>%
  pull(rank_conserved))
order_unsup <- all(score_mat_unsup %>%
                   group_by(name_score,dataset) %>%
                   select(name_score,dataset,candidate,val,normval) %>%
                   mutate(order_before=order(val), order_after=order(normval)) %>%
                   mutate(rank_conserved=order_before==order_after) %>%
                   pull(rank_conserved))
# should be true, we don't want to change ranks during the normalization
