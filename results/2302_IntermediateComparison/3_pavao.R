####
#### Set parameters
####
source("../../src/ranking_ranking_procedure_functions.R")
source("../../src/ranking_pavao_criteria_functions.R")
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

####
#### Winner
####
winner_sup <- rank_sup$candidate[which.max(rank_sup$values)]
winner_unsup <- rank_unsup$candidate[which.max(rank_unsup$values)]

####
#### Pavao judge/candidate matrices
####
mat_sup2 <- score_mat_sup %>%
  ungroup() %>%
  mutate(judge=paste(name_score,dataset,sep="/")) %>%
  select(judge,candidate,trendval) %>%
  mutate(val=trendval) %>%
  select(judge,candidate,val) %>%
  make_pavao_matrix()
mat_sup3 <- score_mat_sup %>%
  group_by(candidate,cat_score,name_score) %>%
  summarise(val=mean(trendval)) %>%
  mutate(judge=name_score) %>%
  ungroup() %>%
  select(judge,candidate,val) %>%
  make_pavao_matrix()
mat_sup4 <- score_mat_sup %>%
  group_by(candidate,cat_score,name_score) %>%
  summarise(score3=mean(trendval)) %>%
  group_by(candidate,cat_score) %>%
  summarise(val=geomMean(score3)) %>%
  mutate(judge=cat_score) %>%
  ungroup() %>%
  select(judge,candidate,val) %>%
  make_pavao_matrix()
matrices_sup <- list(mat_sup2,mat_sup3,mat_sup4)

mat_unsup2 <- score_mat_unsup %>%
  ungroup() %>%
  mutate(judge=paste(name_score,dataset,sep="/")) %>%
  select(judge,candidate,trendval) %>%
  mutate(val=trendval) %>%
  select(judge,candidate,val) %>%
  make_pavao_matrix()
mat_unsup3 <- score_mat_unsup %>%
  group_by(candidate,cat_score,name_score) %>%
  summarise(val=mean(trendval)) %>%
  mutate(judge=name_score) %>%
  ungroup() %>%
  select(judge,candidate,val) %>%
  make_pavao_matrix()
mat_unsup4 <- score_mat_unsup %>%
  group_by(candidate,cat_score,name_score) %>%
  summarise(score3=mean(trendval)) %>%
  group_by(candidate,cat_score) %>%
  summarise(val=geomMean(score3)) %>%
  mutate(judge=cat_score) %>%
  ungroup() %>%
  select(judge,candidate,val) %>%
  make_pavao_matrix()
matrices_unsup <- list(mat_unsup2,mat_unsup3,mat_unsup4)

f <- list(function(scores) {
  name_data <- unique(sapply(colnames(scores), function(x) strsplit(x,"/")[[1]][2]))
  name_score <- unique(sapply(colnames(scores), function(x) strsplit(x,"/")[[1]][1]))
  sub_score = list()
  for (score in name_score) {
    sub_mat <- scores[,grep(score,colnames(scores))]
    if (is.null(dim(sub_mat))) {sub_score[[score]]=sub_mat}
    else {sub_score[[score]] <- rowMeans(sub_mat)}
  }
  sub_score = do.call(rbind,sub_score)
  name_cat = list("raw_perf"=c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s"),
                  "stab_perf"=c("pearson sd_c", "pearson sd_s"),
                  "stab_perf2"=c("mae sd_g", "pearson sd_g", "rmse sd_g"),
                  "time"=c("time","time sd"))
  sub_cat = list()
  for (cat in names(name_cat)) {
    sub_mat <- sub_score[rownames(sub_score) %in% name_cat[[cat]],]
    if (is.null(dim(sub_mat))) {sub_cat[[cat]]=sub_mat}
    else {sub_cat[[cat]] <- apply(sub_mat,2,geomMean)}
  }
  sub_cat = do.call(rbind,sub_cat)
  return(apply(sub_cat,2,geomMean))
},
          function(scores) {
            name_cat = list("raw_perf"=c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s"),
                            "stab_perf"=c("pearson sd_c", "pearson sd_s"),
                            "stab_perf2"=c("mae sd_g", "pearson sd_g", "rmse sd_g"),
                            "time"=c("time","time sd"))
            sub_cat = list()
            for (cat in names(name_cat)) {
              sub_mat <- scores[,colnames(scores) %in% name_cat[[cat]]]
              if (is.null(dim(sub_mat))) {sub_cat[[cat]]=sub_mat}
              else {sub_cat[[cat]] <- apply(sub_mat,1,geomMean)}
              
            }
            sub_cat = do.call(rbind,sub_cat)
            return(apply(sub_cat,2,geomMean))
          },
          function(scores) {
            return(apply(scores,1,geomMean))
          })

####
#### Test Pavao criteria
####
maj_crit_sup <- sapply(matrices_sup, function(x) maj_criterion(x, winner_sup))
condorcet_crit_sup <- sapply(matrices_sup, function(x) condorcet_criterion(x, winner_sup))
consistency_crit_sup <- sapply(matrices_sup, function(x) consistency_criterion(x, winner_sup))
participation_crit_sup <- mapply(function(x, f_M) participation_criterion(x, winner_sup, f_M), matrices_sup, f)
iia_crit_sup <- mapply(function(x, f_M) iia_criterion(x, winner_sup, f_M), matrices_sup, f)

avg_rank_sup <- sapply(matrices_sup, function(x) avg_rank(x, winner_sup))
condorcet_rate_sup <- sapply(matrices_sup, function(x) condorcet_rate(x, winner_sup))
generalization_crit_sup <- mapply(function(x, f_M) generalization_criterion(x, f_M), matrices_sup, f)

maj_crit_unsup <- sapply(matrices_unsup, function(x) maj_criterion(x, winner_unsup))
condorcet_crit_unsup <- sapply(matrices_unsup, function(x) condorcet_criterion(x, winner_sup))
consistency_crit_unsup <- sapply(matrices_unsup, function(x) consistency_criterion(x, winner_sup))
participation_crit_unsup <- mapply(function(x, f_M) participation_criterion(x, winner_unsup, f_M), matrices_unsup, f)
iia_crit_unsup <- mapply(function(x, f_M) iia_criterion(x, winner_unsup, f_M), matrices_unsup, f)

avg_rank_unsup <- sapply(matrices_unsup, function(x) avg_rank(x, winner_unsup))
condorcet_rate_unsup <- sapply(matrices_unsup, function(x) condorcet_rate(x, winner_unsup))
generalization_crit_unsup <- mapply(function(x, f_M) generalization_criterion(x, f_M), matrices_unsup, f)

####
#### Plotting
####
# With all voters (datasets/scores)
pavao_df = data.frame("values"=c(maj_crit_sup[1],condorcet_crit_sup[1],consistency_crit_sup[1],
                                 participation_crit_sup[1],iia_crit_sup[1],
                                 avg_rank_sup[1],condorcet_rate_sup[1],generalization_crit_sup[1],
                                 maj_crit_unsup[1],condorcet_crit_unsup[1],consistency_crit_unsup[1],
                                 participation_crit_unsup[1],iia_crit_unsup[1],
                                 avg_rank_unsup[1],condorcet_rate_unsup[1],generalization_crit_unsup[1]),
                      "ind"=c("Majority criterion","Condorcet criterion","Consistency criterion",
                              "Participation criterion","IIA criterion",
                              "Average rank","Condorcet rate","Generalization criterion"),
                      "deconv_type"=c(rep("Supervised",8),rep("Unsupervised",8)))
cat_criterion = list("Theoretical criteria"=c("Majority criterion","Condorcet criterion","Consistency criterion",
                                              "Participation criterion","IIA criterion"),
                     "Empirical criteria"=c("Average rank","Condorcet rate","Generalization criterion"))
pavao_df$cat_criterion = convert_to_cat(cat_criterion, pavao_df$ind)
ggplot(pavao_df, aes(ind, deconv_type, fill=values)) +
  geom_tile() +
  scale_fill_viridis_c() +
  facet_wrap(~cat_criterion, scales = "free_x") +
  ylab("") +
  xlab("") +
  theme_modern(axis.text.angle = 45)
ggsave(paste0("3_pavao/221101_",block,"_pavao_",rank_type,"_",rank_norm,".pdf"), width = 8, height = 5)
