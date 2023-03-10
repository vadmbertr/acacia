####
#### Set parameters
####
source("../../src/ranking_ranking_procedure_functions.R")
thd_pval = .5

#####
#### Load data
#####
scores1_sup_t <- readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_toast_sup.rds"))
scores1_sup_c <- if (file.exists(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_cimpleg_sup.rds"))) readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_cimpleg_sup.rds"))
scores1_sup_r <- readRDS(paste0("../2210_1SB/perf_scores/221101_scores_",block,"_sup.rds"))
scores1_sup_r$feat_selec = 'None'
if (!is.null(scores1_sup_c) & load_cimpleg) {scores1_sup=bind_rows(scores1_sup_t,scores1_sup_c,scores1_sup_r)} else {scores1_sup=bind_rows(scores1_sup_t,scores1_sup_r)}

scores2_sup_t <- readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_toast_sup.rds"))
scores2_sup_c <- if (file.exists(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_cimpleg_sup.rds"))) readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_cimpleg_sup.rds"))
if (!is.null(scores2_sup_c)) scores2_sup_c$values <- as.numeric(scores2_sup_c$values)
scores2_sup_r <- readRDS(paste0("../2210_1SB/perf_scores/221101_time_",block,"_sup.rds"))
scores2_sup_r$feat_selec = 'None'
scores2_sup_r$values <- as.numeric(scores2_sup_r$values)
if (!is.null(scores2_sup_c) & load_cimpleg) {scores2_sup=bind_rows(scores2_sup_t,scores2_sup_c,scores2_sup_r)} else {scores2_sup=bind_rows(scores2_sup_t,scores2_sup_r)}

scores1_unsup_t <- readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_toast_unsup.rds"))
scores1_unsup_c <- if (file.exists(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_cimpleg_unsup.rds"))) readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_scores_",block,"_cimpleg_unsup.rds"))
scores1_unsup_r <- readRDS(paste0("../2210_1SB/perf_scores/221101_scores_",block,"_unsup.rds"))
scores1_unsup_r$feat_selec = 'None'
if (!is.null(scores1_unsup_c) & load_cimpleg) {scores1_unsup=bind_rows(scores1_unsup_t,scores1_unsup_c,scores1_unsup_r)} else {scores1_unsup=bind_rows(scores1_unsup_t,scores1_unsup_r)}

scores2_unsup_t <- readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_toast_unsup.rds"))
scores2_unsup_c <- if (file.exists(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_cimpleg_unsup.rds"))) readRDS(paste0("../2210_1SB_featselec/perf_scores/221101_time_",block,"_cimpleg_unsup.rds"))
if (!is.null(scores2_unsup_c)) scores2_unsup_c$values <- as.numeric(scores2_unsup_c$values)
scores2_unsup_r <- readRDS(paste0("../2210_1SB/perf_scores/221101_time_",block,"_unsup.rds"))
scores2_unsup_r$feat_selec = 'None'
scores2_unsup_r$values <- as.numeric(scores2_unsup_r$values)
if (!is.null(scores2_unsup_c) & load_cimpleg) {scores2_unsup=bind_rows(scores2_unsup_t,scores2_unsup_c,scores2_unsup_r)} else {scores2_unsup=bind_rows(scores2_unsup_t,scores2_unsup_r)}

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

rm(ranking_avgrank,ranking12,ranking12_notime,ranking1234,ranking12345,ranking12345_notime,ranking_concat_per_method,ranking1245,ranking12453)
rm(geomMean,harmMean,weighMean,MinMaxNorm)
rm(plotting,convert_to_cat)
