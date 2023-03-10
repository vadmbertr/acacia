####
#### Set parameters
####
source("../../src/ranking_ranking_procedure_functions.R")
block='rna' # to vary, either met or rna
feat_selec='toast' # to vary, cimpleg or toast

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
rank_sup <- ranking_avgrank(scores1_sup,scores2_sup)
rank_unsup <- ranking_avgrank(scores1_unsup,scores2_unsup)

df_sup <- score_mat_sup %>%
  mutate(ind=paste(dataset,name_score, sep = "/")) %>%
  group_by(ind) %>%
  mutate(values=rank(trendval, ties.method = "average")) %>%
  select(ind, values, candidate)
df_unsup <- score_mat_unsup %>%
  mutate(ind=paste(dataset,name_score, sep = "/")) %>%
  group_by(ind) %>%
  mutate(values=rank(trendval, ties.method = "average")) %>%
  select(ind, values, candidate)
df_sup <- bind_rows(rank_sup, df_sup)
df_unsup <- bind_rows(rank_unsup, df_unsup)
factor_x <- unique(df_sup$ind)
df_sup$ind <- factor(df_sup$ind, levels=c(factor_x[grep("time", factor_x)],
                                                factor_x[grep("perf_g", factor_x)],
                                                factor_x[grep("med_c", factor_x)],
                                                factor_x[grep("med_s", factor_x)],
                                                factor_x[grep("sd_g", factor_x)],
                                                factor_x[grep("sd_s", factor_x)],
                                                factor_x[grep("sd_c", factor_x)],
                                                factor_x[grep("std", factor_x)],
                                                "avg_rank"))
df_unsup$ind <- factor(df_unsup$ind, levels=c(factor_x[grep("time", factor_x)],
                                                    factor_x[grep("perf_g", factor_x)],
                                                    factor_x[grep("med_c", factor_x)],
                                                    factor_x[grep("med_s", factor_x)],
                                                    factor_x[grep("sd_g", factor_x)],
                                                    factor_x[grep("sd_s", factor_x)],
                                                    factor_x[grep("sd_c", factor_x)],
                                                    factor_x[grep("std", factor_x)],
                                                    "avg_rank"))

####
#### Plotting
####
df_sup$name_score <- sapply(df_sup$ind, function(x) gsub("^.*/","",x))
df_sup$dataset <- sapply(df_sup$ind, function(x) strsplit(as.character(x), "/")[[1]][1])
df_sup$dataset <- factor(df_sup$dataset, levels=c("dBREAST","dPANCREAS","lot1","avg_rank"))
df_sup$facet <- sapply(df_sup$candidate, function(x)  strsplit(x, " ")[[1]][2])
df_sup$deconv <- sapply(df_sup$candidate, function(x)  strsplit(x, " ")[[1]][1])
df_sup$cat_score <- convert_to_cat(list("raw_perf"=c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s"),
                                                     "stab_perf"=c("pearson sd_c", "pearson sd_s"),
                                                     "stab_perf2"=c("mae sd_g", "pearson sd_g", "rmse sd_g"),
                                                     "time"=c("time","time sd"),
                                        "aggregated"=c("avg_rank")), df_sup$name_score)
df_unsup$name_score <- sapply(df_unsup$ind, function(x) gsub("^.*/","",x))
df_unsup$dataset <- sapply(df_unsup$ind, function(x) strsplit(as.character(x), "/")[[1]][1])
df_unsup$dataset <- factor(df_unsup$dataset, levels=c("dBREAST","dPANCREAS","lot1","avg_rank"))
df_unsup$facet <- sapply(df_unsup$candidate, function(x)  strsplit(x, " ")[[1]][2])
df_unsup$deconv <- sapply(df_unsup$candidate, function(x)  strsplit(x, " ")[[1]][1])
df_unsup$cat_score <- convert_to_cat(list("raw_perf"=c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s"),
                                        "stab_perf"=c("pearson sd_c", "pearson sd_s"),
                                        "stab_perf2"=c("mae sd_g", "pearson sd_g", "rmse sd_g"),
                                        "time"=c("time","time sd"),
                                        "aggregated"=c("avg_rank")), df_unsup$name_score)

ggplot(df_sup, aes(x=name_score, y=values, fill=deconv)) +
  geom_boxplot(position=position_dodge(width=.7), alpha=.4) +
  geom_point(position=position_jitterdodge(jitter.width=.1, dodge.width=.7), aes(color=deconv, shape=dataset), size=1.5) +
  facet_wrap(~facet) +
  scale_color_social_d() +
  scale_fill_social_d() +
  theme_modern(axis.text.angle=45)
ggsave(paste0("2b_ranking_avgrank/221101_",block,"_",feat_selec,"_sup_ranking1.pdf"),
       width = 8, height = 5)
ggplot(df_sup, aes(x=dataset, y=values, fill=deconv)) +
  geom_boxplot(position=position_dodge(width=.7), alpha=.4) +
  geom_point(position=position_jitterdodge(jitter.width=.1, dodge.width=.7), aes(color=deconv, alpha=cat_score), size=1.5) +
  facet_wrap(~facet) +
  scale_color_social_d() +
  scale_fill_social_d() +
  theme_modern(axis.text.angle=45)
ggsave(paste0("2b_ranking_avgrank/221101_",block,"_",feat_selec,"_sup_ranking2.pdf"),
       width = 8, height = 5)
ggplot(df_unsup, aes(x=name_score, y=values, fill=deconv)) +
  geom_boxplot(position=position_dodge(width=.7), alpha=.4) +
  geom_point(position=position_jitterdodge(jitter.width=.1, dodge.width=.7), aes(color=deconv, shape=dataset), size=1.5) +
  facet_wrap(~facet) +
  scale_color_social_d() +
  scale_fill_social_d() +
  theme_modern(axis.text.angle=45)
ggsave(paste0("2b_ranking_avgrank/221101_",block,"_",feat_selec,"_unsup_ranking1.pdf"),
       width = 8, height = 5)
ggplot(df_unsup, aes(x=dataset, y=values, fill=deconv)) +
  geom_boxplot(position=position_dodge(width=.7), alpha=.4) +
  geom_point(position=position_jitterdodge(jitter.width=.1, dodge.width=.7), aes(color=deconv, alpha=cat_score), size=1.5) +
  scale_color_social_d() +
  scale_fill_social_d() +
  theme_modern(axis.text.angle=45)
ggsave(paste0("2b_ranking_avgrank/221101_",block,"_",feat_selec,"_unsup_ranking2.pdf"),
       width = 8, height = 5)

