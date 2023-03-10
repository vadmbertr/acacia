####
#### Set parameters
####
options(dplyr.summarise.inform = FALSE)
source("../../src/ranking_ranking_procedure_functions.R")
block='met' # to vary, either met or rna

####
#### Load data
####
scores1_sup <- readRDS(paste0("../2210_0SB/perf_scores/220922_scores_",block,"_sup.rds"))
scores2_sup <- readRDS(paste0("../2210_0SB/perf_scores/220922_time_",block,"_sup.rds"))
scores1_unsup <- readRDS(paste0("../2210_0SB/perf_scores/220922_scores_",block,"_unsup.rds"))
scores2_unsup <- readRDS(paste0("../2210_0SB/perf_scores/220922_time_",block,"_unsup.rds"))
colnames(scores1_sup)[colnames(scores1_sup)=="score"] <- "scor"
scores1_sup <- scores1_sup %>% mutate(score=paste(scor,setting)) %>% select(colnames(scores2_sup))
colnames(scores1_unsup)[colnames(scores1_unsup)=="score"] <- "scor"
scores1_unsup <- scores1_unsup %>% mutate(score=paste(scor,setting)) %>% select(colnames(scores2_unsup))

####
#### Ranking
####
score_mat_sup2 <- ranking12(scores1_sup,scores2_sup)
score_mat_unsup2 <- ranking12(scores1_unsup,scores2_unsup)
score_mat_sup3 <- score_mat_sup2 %>%
  group_by(name_score,cat_score,deconv_method) %>%
  summarise(value=geomMean(trendval))
score_mat_unsup3 <- score_mat_unsup2 %>%
  group_by(name_score,cat_score,deconv_method) %>%
  summarise(value=geomMean(trendval))
score_mat_sup4 <- score_mat_sup3 %>%
  group_by(deconv_method,cat_score) %>%
  summarise(value=geomMean(value))
score_mat_unsup4 <- score_mat_unsup3 %>%
  group_by(deconv_method,cat_score) %>%
  summarise(value=geomMean(value))
score_mat_sup5 <- score_mat_sup4 %>%
  group_by(deconv_method) %>%
  summarise(value=geomMean(value))
score_mat_unsup5 <- score_mat_unsup4 %>%
  group_by(deconv_method) %>%
  summarise(value=geomMean(value))

####
#### Plotting
####
score_mat_sup2$name_score = factor(score_mat_sup2$name_score, levels=c(unique(grep("perf",score_mat_sup2$name_score,value=T)),
                                                                       unique(grep("med",score_mat_sup2$name_score,value=T)),
                                                                       unique(grep("sd",score_mat_sup2$name_score,value=T)[grep("sd",score_mat_sup2$name_score,value=T)!="time sd"]),
                                                                       "time sd","time"))
score_mat_unsup2$name_score = factor(score_mat_unsup2$name_score, levels=c(unique(grep("perf",score_mat_unsup2$name_score,value=T)),
                                                                       unique(grep("med",score_mat_unsup2$name_score,value=T)),
                                                                       unique(grep("sd",score_mat_unsup2$name_score,value=T)[grep("sd",score_mat_unsup2$name_score,value=T)!="time sd"]),
                                                                       "time sd","time"))
ggplot(score_mat_sup2[score_mat_sup2$dataset=="lot1",], aes(x=name_score,y=trendval,color=deconv_method)) +
  geom_point() +
  geom_line(aes(group=deconv_method), size=4, alpha=.5) +
  geom_line(aes(group=deconv_method)) +
  scale_x_discrete(labels=c("MAE","PEARSON","RMSE","PEARSON per celltype","PEARSON per sample","MAE std","PEARSON std","RMSE std","PEARSON std per celltype","PEARSON std per sample","TIME std","TIME")) +
  theme_modern(axis.text.angle = 45) +
  xlab("") +
  ylab("Normalized score")
ggplot(score_mat_unsup2[score_mat_unsup2$dataset=="dBREAST",], aes(x=name_score,y=trendval,color=deconv_method)) +
  geom_point() +
  geom_line(aes(group=deconv_method), size=4, alpha=.5) +
  geom_line(aes(group=deconv_method)) +
  scale_x_discrete(labels=c("MAE","PEARSON","RMSE","PEARSON per celltype","PEARSON per sample","MAE std","PEARSON std","RMSE std","PEARSON std per celltype","PEARSON std per sample","TIME std","TIME")) +
  theme_modern(axis.text.angle = 45) +
  xlab("") +
  ylab("Normalized score")

score_mat_sup3$name_score = factor(score_mat_sup3$name_score, levels=c(unique(grep("perf",score_mat_sup3$name_score,value=T)),
                                                                        unique(grep("med",score_mat_sup3$name_score,value=T)),
                                                                        unique(grep("sd",score_mat_sup3$name_score,value=T)[grep("sd",score_mat_sup3$name_score,value=T)!="time sd"]),
                                                                        "time sd","time"))
score_mat_unsup3$name_score = factor(score_mat_unsup3$name_score, levels=c(unique(grep("perf",score_mat_unsup3$name_score,value=T)),
                                                                       unique(grep("med",score_mat_unsup3$name_score,value=T)),
                                                                       unique(grep("sd",score_mat_unsup3$name_score,value=T)[grep("sd",score_mat_unsup3$name_score,value=T)!="time sd"]),
                                                                       "time sd","time"))
ggplot(score_mat_sup3, aes(x=name_score,y=value,color=deconv_method)) +
  geom_point() +
  geom_line(aes(group=deconv_method), size=4, alpha=.5) +
  geom_line(aes(group=deconv_method)) +
  scale_x_discrete(labels=c("MAE","PEARSON","RMSE","PEARSON per celltype","PEARSON per sample","MAE std","PEARSON std","RMSE std","PEARSON std per celltype","PEARSON std per sample","TIME std","TIME")) +
  theme_modern(axis.text.angle = 45) +
  xlab("") +
  ylab("Aggregated score")
ggplot(score_mat_unsup3, aes(x=name_score,y=value,color=deconv_method)) +
  geom_point() +
  geom_line(aes(group=deconv_method), size=4, alpha=.5) +
  geom_line(aes(group=deconv_method)) +
  scale_x_discrete(labels=c("MAE","PEARSON","RMSE","PEARSON per celltype","PEARSON per sample","MAE std","PEARSON std","RMSE std","PEARSON std per celltype","PEARSON std per sample","TIME std","TIME")) +
  theme_modern(axis.text.angle = 45) +
  xlab("") +
  ylab("Aggregated score")

ggplot(score_mat_sup4, aes(x=cat_score,y=value,color=deconv_method)) +
  geom_point() +
  geom_line(aes(group=deconv_method), size=4, alpha=.5) +
  geom_line(aes(group=deconv_method)) +
  scale_x_discrete(labels=c("Performance","Stability","Stability2","Time")) +
  theme_modern(axis.text.angle = 45) +
  xlab("") +
  ylab("Aggregated score")
ggplot(score_mat_unsup4, aes(x=cat_score,y=value,color=deconv_method)) +
  geom_point() +
  geom_line(aes(group=deconv_method), size=4, alpha=.5) +
  geom_line(aes(group=deconv_method)) +
  scale_x_discrete(labels=c("Performance","Stability","Stability2","Time")) +
  theme_modern(axis.text.angle = 45) +
  xlab("") +
  ylab("Aggregated score")

ggplot(score_mat_sup5, aes(y=value,color=deconv_method)) +
  geom_point(x=0.5) +
  geom_point(size=4, alpha=.5, x=0.5) +
  theme_modern(axis.text.angle = 45) +
  xlab("") +
  ylab("Aggregated score")

df_final_sup <- data.frame("value"=c(score_mat_sup5$value,score_mat_sup4$value),
                       "ind"=c(rep("ScoreFinal",nrow(score_mat_sup5)),score_mat_sup4$cat_score),
                       "deconv_method"=c(score_mat_sup5$deconv_method,score_mat_sup4$deconv_method))
df_final_unsup <- data.frame("value"=c(score_mat_unsup5$value,score_mat_unsup4$value),
                           "ind"=c(rep("ScoreFinal",nrow(score_mat_unsup5)),score_mat_unsup4$cat_score),
                           "deconv_method"=c(score_mat_unsup5$deconv_method,score_mat_unsup4$deconv_method))
df_final_sup$method_type <- "Supervised"
df_final_unsup$method_type <- "Unsupervised"
df_final = bind_rows(df_final_sup,df_final_unsup)
df_final$ind <- factor(df_final$ind, levels=c(unique(score_mat_sup4$cat_score),"ScoreFinal"))
ggplot(df_final, aes(x=ind,y=value,color=deconv_method)) +
  geom_point() +
  geom_line(aes(group=deconv_method), size=4, alpha=.5) +
  geom_line(aes(group=deconv_method)) +
  theme_modern(axis.text.angle = 45) +
  facet_wrap(~method_type) +
  xlab("")
