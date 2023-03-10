####
#### Set parameters
####
library(see)
library(ggplot2)
source("../../src/ranking_pval_functions.R")
n_meth = 3 # to vary, 1 or 3
source("1_compare_.R")
thd_pval = .5
rank_norm="normDataset" # to vary, either normDataset or normGlobal
trendvals=ifelse(rank_norm=="normDataset",trendvals12_notime,trendvals_2_notime)

####
#### Compare pvalues
####
trendval_sup <- trendvals(score_tot)
names_score <- sort(unique(trendval_sup$name_score))
pvals_sup = list()
for (score in names_score) {
  data_sup = trendval_sup[trendval_sup$name_score==score,c("candidate","trendval","dataset","simulation")]
  data_sup$deconv = sapply(data_sup$candidate, function(x) strsplit(x," ")[[1]][1])
  X = data_sup$trendval[data_sup$deconv=="MM"]
  Y = data_sup$trendval[data_sup$deconv=="RM"]
  Z = data_sup$trendval[data_sup$deconv=="MB"]
  pvals_sup[[score]] <- c(pval(Y,X),
                          pval(Z,X),
                          pval(Z,Y))
}
rm(data_sup,X,Y,Z,score)
pvals_sup_df = data.frame("pval"=unlist(pvals_sup),
                          "comparison"=rep(rep(c("MET/RNA","MET/MB","RNA/MB"),length(pvals_sup[[1]])),time=length(pvals_sup)),
                          "score"=rep(names(pvals_sup),lengths(pvals_sup)))
pvals_sup_df$minuslogpval = -log10(pvals_sup_df$pval)
pvals_sup_df$comparison = factor(pvals_sup_df$comparison, levels=c("MET/RNA","MET/MB","RNA/MB"))

pass_thd_sup = pvals_sup_df %>%
    group_by(comparison) %>%
    summarise("pass_thd"=sum(minuslogpval>2)>=max(2,ceiling(thd_pval*length(names_score))),"yval"=max(minuslogpval)+3) %>%
    filter(pass_thd)

####
#### Plot pvalues
####
ggplot(pvals_sup_df, aes(x=comparison, y=minuslogpval, color=score)) +
  geom_point() +
  geom_label(data=pass_thd_sup, aes(x=comparison, y=yval), label="*", color="black", size=5) +
  geom_hline(yintercept=2) +
  scale_color_social_d() +
  theme_lucid()
ggsave(paste0("2_compare_pval/221101_",n_meth,"meth_sup_pval_",rank_norm,".pdf"),
       width = 8, height = 5)

####
#### Extract significant comparisons
####
pass_thd_sup %>% select(comparison)
table(pass_thd_sup %>% pull(comparison))
