####
#### Set parameters
####
source("../../src/ranking_ranking_procedure_functions.R")
block='rna' # to vary, either met or rna

####
#### Load data
####
scores1_sup <- readRDS(paste0("../2210_1SB/perf_scores/221101_scores_",block,"_sup.rds"))
scores2_sup <- readRDS(paste0("../2210_1SB/perf_scores/221101_time_",block,"_sup.rds"))
scores1_unsup <- readRDS(paste0("../2210_1SB/perf_scores/221101_scores_",block,"_unsup.rds"))
scores2_unsup <- readRDS(paste0("../2210_1SB/perf_scores/221101_time_",block,"_unsup.rds"))
colnames(scores1_sup)[colnames(scores1_sup)=="score"] <- "scor"
scores1_sup <- scores1_sup %>% mutate(score=paste(scor,setting), candidate=deconv) %>% select(values,score,sim,dataset,candidate)
colnames(scores1_unsup)[colnames(scores1_unsup)=="score"] <- "scor"
scores1_unsup <- scores1_unsup %>% mutate(score=paste(scor,setting), candidate=deconv) %>% select(values,score,sim,dataset,candidate)
scores2_sup <- scores2_sup %>% mutate(candidate=deconv) %>% select(values,score,sim,dataset,candidate)
scores2_unsup <- scores2_unsup %>% mutate(candidate=deconv) %>% select(values,score,sim,dataset,candidate)

####
#### Ranking
####
eval_norm_order(ranking12(scores1_sup,scores2_sup))
eval_norm_order(ranking12(scores1_unsup,scores2_unsup))
eval_norm_order(ranking_2(scores1_sup,scores2_sup))
eval_norm_order(ranking_2(scores1_unsup,scores2_unsup))
