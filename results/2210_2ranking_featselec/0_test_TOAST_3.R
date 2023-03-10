## ----
## Set parameters
## ----
source("../../src/ranking_ranking_procedure_functions.R")
block="met" # to vary, either met or rna

## ----
## load
## ----
scores1_sup <- readRDS(paste0("../2210_1SB_featselec/0_test_TOAST/scores/scores_",block,"_sup.rds"))
scores1_sup_bl <- readRDS(paste0("../2210_1SB_featselec/0_test_TOAST/scores/scoresbl_",block,"_sup.rds"))
scores1_sup_bl2 <- readRDS(paste0("../2210_1SB_featselec/0_test_TOAST/scores/scoresbl2_",block,"_sup.rds"))
scores1_sup = bind_rows(scores1_sup,scores1_sup_bl,scores1_sup_bl2)

scores2_sup <- readRDS(paste0("../2210_1SB_featselec/0_test_TOAST/scores/time_",block,"_sup.rds"))
scores2_sup_bl <- readRDS(paste0("../2210_1SB_featselec/0_test_TOAST/scores/timebl_",block,"_sup.rds"))
scores2_sup_bl2 <- readRDS(paste0("../2210_1SB_featselec/0_test_TOAST/scores/timebl2_",block,"_sup.rds"))
scores2_sup = bind_rows(scores2_sup,scores2_sup_bl,scores2_sup_bl2)

scores1_unsup <- readRDS(paste0("../2210_1SB_featselec/0_test_TOAST/scores/scores_",block,"_unsup.rds"))
scores1_unsup_bl <- readRDS(paste0("../2210_1SB_featselec/0_test_TOAST/scores/scoresbl_",block,"_unsup.rds"))
scores1_unsup_bl2 <- readRDS(paste0("../2210_1SB_featselec/0_test_TOAST/scores/scoresbl2_",block,"_unsup.rds"))
scores1_unsup = bind_rows(scores1_unsup,scores1_unsup_bl,scores1_unsup_bl2)

scores2_unsup <- readRDS(paste0("../2210_1SB_featselec/0_test_TOAST/scores/time_",block,"_unsup.rds"))
scores2_unsup_bl <- readRDS(paste0("../2210_1SB_featselec/0_test_TOAST/scores/timebl_",block,"_unsup.rds"))
scores2_unsup_bl2 <- readRDS(paste0("../2210_1SB_featselec/0_test_TOAST/scores/timebl2_",block,"_unsup.rds"))
scores2_unsup = bind_rows(scores2_unsup,scores2_unsup_bl,scores2_unsup_bl2)

rm(scores1_sup_bl,scores1_sup_bl2,scores2_sup_bl,scores2_sup_bl2,
   scores1_unsup_bl,scores1_unsup_bl2,scores2_unsup_bl,scores2_unsup_bl2)

## ----
## Organize df
## ----
colnames(scores1_sup)[colnames(scores1_sup)=="score"] <- "scor"
scores1_sup <- scores1_sup %>% mutate(score=paste(scor,setting)) %>% select(colnames(scores2_sup))
colnames(scores1_unsup)[colnames(scores1_unsup)=="score"] <- "scor"
scores1_unsup <- scores1_unsup %>% mutate(score=paste(scor,setting)) %>% select(colnames(scores2_unsup))

scores1_sup <- scores1_sup %>%
  mutate(candidate=paste(deconv,feat_selec),sim=NA,dataset='TOAST') %>%
  select(values,score,sim,dataset,candidate)
scores2_sup <- scores2_sup %>%
  mutate(candidate=paste(deconv,feat_selec),sim=NA,dataset='TOAST') %>%
  select(values,score,sim,dataset,candidate)
scores1_unsup <- scores1_unsup %>%
  mutate(candidate=paste(deconv,feat_selec),sim=NA,dataset='TOAST') %>%
  select(values,score,sim,dataset,candidate)
scores2_unsup <- scores2_unsup %>%
  mutate(candidate=paste(deconv,feat_selec),sim=NA,dataset='TOAST') %>%
  select(values,score,sim,dataset,candidate)

## ----
## Ranking
## ----
rank345_sup <- ranking12345(scores1_sup,scores2_sup)
rank345_unsup <- ranking12345(scores1_unsup,scores2_unsup)

rank345_sup$feat_selec <- sapply(rank345_sup$candidate, function(x)
  strsplit(x," ")[[1]][2])
rank345_sup$deconv_method <- sapply(rank345_sup$candidate, function(x)
  strsplit(x," ")[[1]][1])

rank345_unsup$feat_selec <- sapply(rank345_unsup$candidate, function(x)
  strsplit(x," ")[[1]][2])
rank345_unsup$deconv_method <- sapply(rank345_unsup$candidate, function(x)
  strsplit(x," ")[[1]][1])

## ----
## Plotting
## ----
rank345_sup$feat_selec <- factor(rank345_sup$feat_selec, levels=c("None","None_restricted","TOAST"))
rank345_unsup$feat_selec <- factor(rank345_unsup$feat_selec, levels=c("None","None_restricted","TOAST"))
ggplot(rank345_sup, aes(x=feat_selec, y=values, color=deconv_method)) +
  geom_point() +
  geom_line(aes(group=deconv_method)) +
  scale_color_social_d() +
  geom_line(aes(group=deconv_method), size=5, alpha=.4) +
  theme_modern(axis.text.angle = 45)
ggsave(paste0("0_test_TOAST/compare_",block,"_sup_ranking.pdf"),
       width = 8, height = 5)
ggplot(rank345_unsup, aes(x=feat_selec, y=values, color=deconv_method)) +
  geom_point() +
  geom_line(aes(group=deconv_method)) +
  scale_color_social_d() +
  geom_line(aes(group=deconv_method), size=5, alpha=.4) +
  theme_modern(axis.text.angle = 45)
ggsave(paste0("0_test_TOAST/compare_",block,"_unsup_ranking.pdf"),
       width = 8, height = 5)

