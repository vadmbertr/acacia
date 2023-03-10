####
#### Set parameters
####
source("../../src/ranking_ranking_procedure_functions.R")
block='met' # to vary, either met or rna
ranking_norm=ranking12 # to vary, either ranking12 or ranking_2
ranking_type=ranking12453 # to vary, either ranking12345 or ranking12453
rank_norm=ifelse(all.equal(ranking_norm,ranking12)==T,"normDataset","normGlobal")
rank_type=ifelse(all.equal(ranking_type,ranking12345)==T,"order345","order453")

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
score_mat_sup <- ranking_norm(scores1_sup,scores2_sup)
score_mat_unsup <- ranking_norm(scores1_unsup,scores2_unsup)
rank_sup <- ranking_type(scores1_sup,scores2_sup,ranking=ranking_norm)
rank_unsup <- ranking_type(scores1_unsup,scores2_unsup,ranking=ranking_norm)

df_sup <- score_mat_sup %>%
  mutate(ind=name_score) %>%
  group_by(candidate,ind) %>%
  summarise(values=mean(trendval)) %>%
  select(values,ind,candidate)
df_unsup <- score_mat_unsup %>%
  mutate(ind=name_score) %>%
  group_by(candidate,ind) %>%
  summarise(values=mean(trendval)) %>%
  select(values,ind,candidate)
df_sup <- bind_rows(rank_sup, df_sup)
df_unsup <- bind_rows(rank_unsup, df_unsup)
factor_x <- unique(score_mat_sup$name_score)
df_sup$ind <- factor(df_sup$ind, levels=c(factor_x[grep("time", factor_x)],
                                                factor_x[grep("perf_g", factor_x)],
                                                factor_x[grep("med_c", factor_x)],
                                                factor_x[grep("med_s", factor_x)],
                                                factor_x[grep("sd_g", factor_x)],
                                                factor_x[grep("sd_s", factor_x)],
                                                factor_x[grep("sd_c", factor_x)],
                                                factor_x[grep("std", factor_x)],
                                                "ScoreFinal"))
df_unsup$ind <- factor(df_unsup$ind, levels=c(factor_x[grep("time", factor_x)],
                                                    factor_x[grep("perf_g", factor_x)],
                                                    factor_x[grep("med_c", factor_x)],
                                                    factor_x[grep("med_s", factor_x)],
                                                    factor_x[grep("sd_g", factor_x)],
                                                    factor_x[grep("sd_s", factor_x)],
                                                    factor_x[grep("sd_c", factor_x)],
                                                    factor_x[grep("std", factor_x)],
                                                    "ScoreFinal"))

####
#### Plotting
####
plotting(df_sup) +
  geom_line(aes(group=candidate), size=5, alpha=.4)
ggsave(paste0("2b_ranking_aggreg/221101_",block,"_sup_ranking_",rank_type,"_",rank_norm,".pdf"),
       width = 8, height = 5)
plotting(df_unsup) +
  geom_line(aes(group=candidate), size=5, alpha=.4)
ggsave(paste0("2b_ranking_aggreg/221101_",block,"_unsup_ranking_",rank_type,"_",rank_norm,".pdf"),
       width = 8, height = 5)
