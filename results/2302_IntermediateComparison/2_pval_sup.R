#####
##### Set parameters
#####
library(see)
library(ggplot2)
source("../../src/ranking_pval_functions.R")
source("1_eval_sup_.R")
rank_norm="normDataset" # to vary, either normDataset or normGlobal
trendvals=ifelse(rank_norm=="normDataset",trendvals12,trendvals_2)
thd_pval = .5


#####
##### Compute pvalues
#####
trendval <- trendvals(df_tot1,df_tot2)
trendval <- coerce_pearson(trendval)
names_score <- sort(unique(trendval$name_score))
pvals = list()
for (score in names_score) {
  print(score)
  data = trendval[trendval$name_score==score,c("candidate","trendval","dataset","simulation")]
  pvals[[score]] <- pbapply::pblapply(sort(unique(data$dataset)), function(z) 
    sapply(sort(unique(data$candidate)), function(x)
      sapply(sort(unique(data$candidate)), function(y) {
        X <- data[data$candidate==x & data$dataset==z,] %>% arrange(simulation) %>% pull(trendval)
        Y <- data[data$candidate==y & data$dataset==z,] %>% arrange(simulation) %>% pull(trendval)
        ifelse(length(X)==length(Y),ifelse(length(Y[!is.na(Y)])>=10&length(X[!is.na(X)])>=10,return(pval(Y,X)),return(NA)),return(NA))
      })))
}
pvals_df = data.frame("pval"=unlist(pvals),
                      "candidate1"=rep(rep(colnames(pvals[[1]]),each=nrow(pvals[[1]])),length(pvals)), #before
                      "candidate2"=rep(rep(rownames(pvals[[1]]),ncol(pvals[[1]])),length(pvals)), #after
                      "score"=rep(names(pvals),lengths(pvals))) %>% filter(!is.na(pval)) %>%
  mutate(minuslogpval=-log10(pval)) %>% filter(minuslogpval>2)

pass_thd = pvals_df %>%
    group_by(candidate1,candidate2) %>%
    summarise("pass_thd"=sum(minuslogpval>2)>=max(2,ceiling(thd_pval*length(names_score))),"yval"=max(minuslogpval)+3) %>%
    filter(pass_thd)

#####
##### Plot pvalues
#####
pvals_df$datatype1 <- sapply(pvals_df$candidate1, function(x) strsplit(x," ")[[1]][1])
pvals_df$datatype2 <- sapply(pvals_df$candidate2, function(x) strsplit(x," ")[[1]][1])
pass_thd$datatype1 <- sapply(pass_thd$candidate1, function(x) strsplit(x," ")[[1]][1])
pass_thd$datatype2 <- sapply(pass_thd$candidate2, function(x) strsplit(x," ")[[1]][1])
ggplot(pvals_df[pvals_df$datatype1=="SB"&pvals_df$datatype2=="SB",], aes(x=candidate1, y=minuslogpval, color=score)) +
  geom_point() +
  geom_label(data=pass_thd[pass_thd$datatype1=="SB"&pass_thd$datatype2=="SB",], aes(x=candidate1, y=yval), label="*", color="black", size=5) +
  geom_hline(yintercept=2) +
  scale_color_social_d() +
  facet_wrap(~candidate2, scales="free_x") +
  theme_lucid(axis.text.angle=45)
ggsave(paste0("2_pval_sup/221101_pval_SB_",rank_norm,".pdf"),
       width = 25, height = 20)
ggplot(pvals_df[pvals_df$datatype1=="SB"&pvals_df$datatype2=="MB",], aes(x=candidate1, y=minuslogpval, color=score)) +
  geom_point() +
  geom_label(data=pass_thd[pass_thd$datatype1=="SB"&pass_thd$datatype2=="MB",], aes(x=candidate1, y=yval), label="*", color="black", size=5) +
  geom_hline(yintercept=2) +
  scale_color_social_d() +
  facet_wrap(~candidate2, scales="free_x") +
  theme_lucid(axis.text.angle=45)
ggsave(paste0("2_pval_sup/221101_pval_MBvsSB_",rank_norm,".pdf"),
       width = 20, height = 15)
ggplot(pvals_df[pvals_df$datatype1=="MB"&pvals_df$datatype2=="MB",], aes(x=candidate1, y=minuslogpval, color=score)) +
  geom_point() +
  geom_label(data=pass_thd[pass_thd$datatype1=="MB"&pass_thd$datatype2=="MB",], aes(x=candidate1, y=yval), label="*", color="black", size=5) +
  geom_hline(yintercept=2) +
  scale_color_social_d() +
  facet_wrap(~candidate2, scales="free_x") +
  theme_lucid(axis.text.angle=45)
ggsave(paste0("2_pval_sup/221101_pval_MBvsMB_",rank_norm,".pdf"),
       width = 20, height = 15)
