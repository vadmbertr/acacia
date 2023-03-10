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
score_mat_sup <- ranking12(scores1_sup,scores2_sup)
score_mat_unsup <- ranking12(scores1_unsup,scores2_unsup)
rank345_sup <- ranking12345(scores1_sup,scores2_sup)
rank345_unsup <- ranking12345(scores1_unsup,scores2_unsup)
rank453_sup <- ranking12453(scores1_sup,scores2_sup)
rank453_unsup <- ranking12453(scores1_unsup,scores2_unsup)
rank34_sup <- ranking1234(scores1_sup,scores2_sup)
rank34_unsup <- ranking1234(scores1_unsup,scores2_unsup)
rank45_sup <- ranking1245(scores1_sup,scores2_sup)
rank45_unsup <- ranking1245(scores1_unsup,scores2_unsup)

df345_sup <- bind_rows(rank345_sup, rank34_sup)
df345_unsup <- bind_rows(rank345_unsup, rank34_unsup)
df453_sup <- bind_rows(rank453_sup, rank45_sup)
df453_unsup <- bind_rows(rank453_unsup, rank45_unsup)
df345_sup$ind <- factor(df345_sup$ind, levels=c("time","raw_perf","stab_perf","stab_perf2","ScoreFinalGG"))
df345_unsup$ind <- factor(df345_unsup$ind, levels=c("time","raw_perf","stab_perf","stab_perf2","ScoreFinalGG"))
df453_sup$ind <- factor(df453_sup$ind, levels=c("dBREAST","dPANCREAS","lot1","ScoreFinalGG"))
df453_unsup$ind <- factor(df453_unsup$ind, levels=c("dBREAST","dPANCREAS","lot1","ScoreFinalGG"))

####
#### Plotting
####
p1_sup=plotting(df345_sup) +
  geom_line(aes(group=candidate), size=5, alpha=.4)
p2_sup=plotting(df453_sup) +
  geom_line(aes(group=candidate), size=5, alpha=.4)
ggpubr::ggarrange(p1_sup,p2_sup, common.legend = T)
ggsave(paste0("2a_ranking_eval_steporder/221101_",block,"_sup_ranking_eval.pdf"),
       width = 8, height = 5)

p1_unsup=plotting(df345_unsup) +
  geom_line(aes(group=candidate), size=5, alpha=.4)
p2_unsup=plotting(df453_unsup) +
  geom_line(aes(group=candidate), size=5, alpha=.4)
ggpubr::ggarrange(p1_unsup,p2_unsup, common.legend = T)
ggsave(paste0("2a_ranking_eval_steporder/221101_",block,"_unsup_ranking_eval.pdf"),
       width = 8, height = 5)
