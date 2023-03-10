####
#### Set parameters
####
library(see)
library(ggplot2)
source("../../src/ranking_ranking_procedure_functions.R")
source("../../src/ranking_pval_functions.R")
block='rna' # to vary, either met or rna
ranking_norm=ranking12 # to vary, either ranking12 or ranking_2
ranking_type=ranking12345 # to vary, either ranking12345 or ranking12453
rank_norm=ifelse(all.equal(ranking_norm,ranking12)==T,"normDataset","normGlobal")
rank_type=ifelse(all.equal(ranking_type,ranking12345)==T,"order345","order453")
thd_pval = .5

#####
#### Load data
#####
scores1_sup_t <- readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_toast_sup.rds"))
scores1_sup_c <- if (file.exists(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_cimpleg_sup.rds"))) readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_cimpleg_sup.rds"))
scores1_sup_r <- readRDS(paste0("../2210_1SB/perf_scores/221101_scores_",block,"_sup.rds"))
scores1_sup_r$feat_selec = 'None'
if (!is.null(scores1_sup_c)) {scores1_sup=bind_rows(scores1_sup_t,scores1_sup_c,scores1_sup_r)} else {scores1_sup=bind_rows(scores1_sup_t,scores1_sup_r)}

scores2_sup_t <- readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_toast_sup.rds"))
scores2_sup_c <- if (file.exists(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_cimpleg_sup.rds"))) readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_cimpleg_sup.rds"))
if (!is.null(scores2_sup_c)) scores2_sup_c$values <- as.numeric(scores2_sup_c$values)
scores2_sup_r <- readRDS(paste0("../2210_1SB/perf_scores/221101_time_",block,"_sup.rds"))
scores2_sup_r$feat_selec = 'None'
scores2_sup_r$values <- as.numeric(scores2_sup_r$values)
if (!is.null(scores2_sup_c)) {scores2_sup=bind_rows(scores2_sup_t,scores2_sup_c,scores2_sup_r)} else {scores2_sup=bind_rows(scores2_sup_t,scores2_sup_r)}

scores1_unsup_t <- readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_toast_unsup.rds"))
scores1_unsup_c <- if (file.exists(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_cimpleg_unsup.rds"))) readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_cimpleg_unsup.rds"))
scores1_unsup_r <- readRDS(paste0("../2210_1SB/perf_scores/221101_scores_",block,"_unsup.rds"))
scores1_unsup_r$feat_selec = 'None'
if (!is.null(scores1_unsup_c)) {scores1_unsup=bind_rows(scores1_unsup_t,scores1_unsup_c,scores1_unsup_r)} else {scores1_unsup=bind_rows(scores1_unsup_t,scores1_unsup_r)}

scores2_unsup_t <- readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_toast_unsup.rds"))
scores2_unsup_c <- if (file.exists(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_cimpleg_unsup.rds"))) readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_cimpleg_unsup.rds"))
if (!is.null(scores2_unsup_c)) scores2_unsup_c$values <- as.numeric(scores2_unsup_c$values)
scores2_unsup_r <- readRDS(paste0("../2210_1SB/perf_scores/221101_time_",block,"_unsup.rds"))
scores2_unsup_r$feat_selec = 'None'
scores2_unsup_r$values <- as.numeric(scores2_unsup_r$values)
if (!is.null(scores2_unsup_c)) {scores2_unsup=bind_rows(scores2_unsup_t,scores2_unsup_c,scores2_unsup_r)} else {scores2_unsup=bind_rows(scores2_unsup_t,scores2_unsup_r)}

rm(scores1_sup_c,scores1_sup_t,scores1_sup_r,scores2_sup_c,scores2_sup_t,scores2_sup_r,scores1_unsup_c,scores1_unsup_t,scores1_unsup_r,scores2_unsup_c,scores2_unsup_t,scores2_unsup_r)

#####
#### Organize df
#####
colnames(scores1_sup)[colnames(scores1_sup)=="score"] <- "scor"
scores1_sup <- scores1_sup %>% mutate(score=paste(scor,setting)) %>% select(colnames(scores2_sup))
colnames(scores1_unsup)[colnames(scores1_unsup)=="score"] <- "scor"
scores1_unsup <- scores1_unsup %>% mutate(score=paste(scor,setting)) %>% select(colnames(scores2_unsup))

scores1_sup <- scores1_sup[!(scores1_sup$deconv %in% c("CIBERSORT4","svr4")),]
scores1_sup$deconv[scores1_sup$deconv=="CIBERSORT3"] <- "CIBERSORT"
scores1_sup$deconv[scores1_sup$deconv=="svr3"] <- "svr"
scores2_sup <- scores2_sup[!(scores2_sup$deconv %in% c("CIBERSORT4","svr4")),]
scores2_sup$deconv[scores2_sup$deconv=="CIBERSORT3"] <- "CIBERSORT"
scores2_sup$deconv[scores2_sup$deconv=="svr3"] <- "svr"

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

#####
#### Ranking
#####
score_mat_sup <- ranking_norm(scores1_sup,scores2_sup)
score_mat_unsup <- ranking_norm(scores1_unsup,scores2_unsup)
rank_sup <- ranking_type(scores1_sup,scores2_sup,ranking=ranking_norm)
rank_unsup <- ranking_type(scores1_unsup,scores2_unsup,ranking=ranking_norm)

rank_sup$feat_selec <- unname(sapply(rank_sup$candidate, function(x)
  strsplit(x," ")[[1]][2]))
rank_sup$deconv_method <- unname(sapply(rank_sup$candidate, function(x)
  strsplit(x," ")[[1]][1]))

rank_unsup$feat_selec <- sapply(rank_unsup$candidate, function(x)
  strsplit(x," ")[[1]][2])
rank_unsup$deconv_method <- sapply(rank_unsup$candidate, function(x)
  strsplit(x," ")[[1]][1])

####
#### Plotting
####
rank_sup$feat_selec <- factor(rank_sup$feat_selec, levels=c("None","None_restrictedC","CimpleG","None_restrictedT","TOAST"))
rank_unsup$feat_selec <- factor(rank_unsup$feat_selec, levels=c("None","None_restrictedC","CimpleG","None_restrictedT","TOAST"))
 
ggplot(rank_sup, aes(x=feat_selec, y=values, color=deconv_method)) +
  geom_point() +
  geom_line(aes(group=deconv_method)) +
  scale_color_social_d() +
  geom_line(aes(group=deconv_method), size=5, alpha=.4) +
  theme_modern(axis.text.angle = 45)
ggsave(paste0("3a_compare/221101_",block,"_sup_ranking_",rank_type,"_",rank_norm,".pdf"),
       width = 8, height = 5)
ggplot(rank_unsup, aes(x=feat_selec, y=values, color=deconv_method)) +
  geom_point() +
  geom_line(aes(group=deconv_method)) +
  scale_color_social_d() +
  geom_line(aes(group=deconv_method), size=5, alpha=.4) +
  theme_modern(axis.text.angle = 45)
ggsave(paste0("3a_compare/221101_",block,"_unsup_ranking_",rank_type,"_",rank_norm,".pdf"),
       width = 8, height = 5)
