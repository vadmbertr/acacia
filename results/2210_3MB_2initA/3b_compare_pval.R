####
#### Set parameters
####
library(see)
library(ggplot2)
init_type='RunsuptoM' # to vary, either MsuptoR, MunsuptoR, RsuptoM or RunsuptoM
source("3_compare_.R")
source("../../src/ranking_pval_functions.R")
rank_norm="normGlobal" # to vary, either normDataset or normGlobal
trendvals=ifelse(rank_norm=="normDataset",trendvals12,trendvals_2)
thd_pval = .5


####
#### Compare pvalues
####
trendval <- trendvals(score_tot1,score_tot2)
names_score <- sort(unique(trendval$name_score))
pvals = list()
for (score in names_score) {
  data = trendval[trendval$name_score==score,c("candidate","trendval","dataset","simulation")]
  data$deconv_in = sapply(data$candidate, function(x) strsplit(x," ")[[1]][1])
  data$deconv_out = sapply(data$candidate, function(x) strsplit(x," ")[[1]][2])
  data$baseline <- sapply(data$deconv_in, function(x) ifelse(x=="no","yes","no"))
  pvals[[score]] <- sapply(sort(unique(data$deconv_out)), function(x) {
    data0 = data[data$deconv_out==x,]
    data0 <- data0 %>% arrange(dataset,simulation)
    data00 <- data0[data0$baseline=="yes",]
    X <- data00$trendval[!duplicated(paste(data00$dataset,data00$simulation,data00$candidate))]
    Y <- sapply(sort(unique(data0$deconv_in[data0$baseline=="no"])), function(y)
      data0$trendval[data0$deconv_in==y])
    return(apply(Y,2,function(z) pval(z,X)))
  })
}
pvals_df = data.frame("pval"=unlist(pvals),
                      "comparison"=rep(rep(rownames(pvals[[1]]),ncol(pvals[[1]])),length(pvals)),
                      "deconv_out"=rep(rep(colnames(pvals[[1]]),each=nrow(pvals[[1]])),length(pvals)),
                      "score"=rep(names(pvals),lengths(pvals)))
pvals_df$minuslogpval = -log10(pvals_df$pval)

pass_thd = pvals_df %>%
    group_by(deconv_out, comparison) %>%
    summarise("pass_thd"=sum(minuslogpval>2)>=max(2,ceiling(thd_pval*length(names_score))),"yval"=max(minuslogpval)+3) %>%
    filter(pass_thd)

####
#### Plot pvalues
####
ggplot(pvals_df, aes(x=comparison, y=minuslogpval, color=score)) +
  geom_point() +
  geom_label(data=pass_thd, aes(x=comparison, y=yval), label="*", color="black", size=5) +
  geom_hline(yintercept=2) +
  scale_color_social_d() +
  facet_wrap(~deconv_out) +
  theme_lucid()
ggsave(paste0("3b_compare_pval/221101_",init_type,"_pval_",rank_norm,".pdf"),
       width = 8, height = 5)

####
#### Extract significant comparisons
####
pass_thd %>% select(deconv_out,comparison)
table(pass_thd %>% pull(comparison))
table(pass_thd %>% pull(deconv_out))
