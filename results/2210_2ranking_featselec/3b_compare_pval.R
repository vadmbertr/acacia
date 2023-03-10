####
#### Set parameters
####
library(see)
library(ggplot2)
block='rna' # to vary, either met or rna
load_cimpleg = F
source("3a_compare_.R")
source("../../src/ranking_pval_functions.R")
feat_selec_meth = "TOAST" # Let's focus on TOAST, and not CimpleG, for now
trendvals_norm=trendvals12 # to vary, either trendvals12 or trendvals_2
tv_norm=ifelse(all.equal(trendvals_norm,trendvals12)==T,"normDataset","normGlobal")
thd_pval = .5


####
#### Compare pvalues
####
trendval_sup <- trendvals_norm(scores1_sup,scores2_sup)
trendval_unsup <- trendvals_norm(scores1_unsup,scores2_unsup)
names_score <- sort(unique(trendval_sup$name_score))
pvals_sup = list()
pvals_unsup = list()
for (score in names_score) {
  data_sup = trendval_sup[trendval_sup$name_score==score,c("candidate","trendval","dataset","simulation")]
  data_sup$deconv = sapply(data_sup$candidate, function(x) strsplit(x," ")[[1]][1])
  data_sup$feat_selec = sapply(data_sup$candidate, function(x) strsplit(x," ")[[1]][2])
  pvals_sup[[score]] <- sapply(sort(unique(data_sup$deconv)), function(x) {
    data0 = data_sup[data_sup$deconv==x,]
    data0 <- data0 %>% arrange(dataset,simulation)
    X0 = data0[data0$feat_selec=="None",c(3:5)]
    Y0 = data0[data0$feat_selec==paste0("None_restricted",strsplit(feat_selec_meth,"")[[1]][1]),c(2:5)]
    Z0 = data0[data0$feat_selec==feat_selec_meth,c(2:5)]
    Y0 <- left_join(X0,Y0,by=c("dataset", "simulation", "deconv"))
    Z0 <- left_join(X0,Z0,by=c("dataset", "simulation", "deconv"))
    X = data0$trendval[data0$feat_selec=="None"]
    Y = Y0$trendval
    Z = Z0$trendval
    return(c(ifelse(all(is.na(Y)),NA,pval(Y,X)),
             ifelse(all(is.na(Z)),NA,pval(Z,X)),
             ifelse(all(is.na(Z))|all(is.na(Y)),NA,pval(Z,Y))))
  })
  data_unsup = trendval_unsup[trendval_unsup$name_score==score,c("candidate","trendval","dataset","simulation")]
  data_unsup$deconv = sapply(data_unsup$candidate, function(x) strsplit(x," ")[[1]][1])
  data_unsup$feat_selec = sapply(data_unsup$candidate, function(x) strsplit(x," ")[[1]][2])
  pvals_unsup[[score]] <- sapply(sort(unique(data_unsup$deconv)), function(x) {
    data0 = data_unsup[data_unsup$deconv==x,]
    data0 <- data0 %>% arrange(dataset,simulation)
    X0 = data0[data0$feat_selec=="None",c(3:5)]
    Y0 = data0[data0$feat_selec==paste0("None_restricted",strsplit(feat_selec_meth,"")[[1]][1]),c(2:5)]
    Z0 = data0[data0$feat_selec==feat_selec_meth,c(2:5)]
    Y0 <- left_join(X0,Y0,by=c("dataset", "simulation", "deconv"))
    Z0 <- left_join(X0,Z0,by=c("dataset", "simulation", "deconv"))
    X = data0$trendval[data0$feat_selec=="None"]
    Y = Y0$trendval
    Z = Z0$trendval
    return(c(ifelse(all(is.na(Y)),NA,pval(Y,X)),
             ifelse(all(is.na(Z)),NA,pval(Z,X)),
             ifelse(all(is.na(Z))|all(is.na(Y)),NA,pval(Z,Y))))
  })
}
pvals_sup_df = data.frame("pval"=unlist(pvals_sup),
                          "comparison"=rep(rep(c("None/HVG","None/FS","HVG/FS"),ncol(pvals_sup[[1]])),time=length(pvals_sup)),
                          "deconv"=rep(rep(colnames(pvals_sup[[1]]),each=nrow(pvals_sup[[1]])),time=length(pvals_sup)),
                          "score"=rep(names(pvals_sup),lengths(pvals_sup)))
pvals_unsup_df = data.frame("pval"=unlist(pvals_unsup),
                            "comparison"=rep(rep(c("None/HVG","None/FS","HVG/FS"),ncol(pvals_unsup[[1]])),time=length(pvals_unsup)),
                            "deconv"=rep(rep(colnames(pvals_unsup[[1]]),each=nrow(pvals_unsup[[1]])),time=length(pvals_unsup)),
                            "score"=rep(names(pvals_unsup),lengths(pvals_unsup)))

pvals_sup_df$minuslogpval = -log10(pvals_sup_df$pval)
pvals_unsup_df$minuslogpval = -log10(pvals_unsup_df$pval)
pvals_sup_df$comparison = factor(pvals_sup_df$comparison, levels=c("None/HVG","None/FS","HVG/FS"))
pvals_unsup_df$comparison = factor(pvals_unsup_df$comparison, levels=c("None/HVG","None/FS","HVG/FS"))
pvals_sup_df$feat_selec=feat_selec_meth
pvals_unsup_df$feat_selec=feat_selec_meth

pass_thd_sup = pvals_sup_df %>%
    group_by(deconv, comparison, feat_selec) %>%
    summarise("pass_thd"=sum(minuslogpval>2)>=max(2,ceiling(thd_pval*length(names_score))),"yval"=max(minuslogpval)+3) %>%
    filter(pass_thd)
pass_thd_unsup = pvals_unsup_df %>%
    group_by(deconv, comparison, feat_selec) %>%
    summarise("pass_thd"=sum(minuslogpval>2)>=max(2,ceiling(thd_pval*length(names_score))),"yval"=max(minuslogpval)+1) %>%
    filter(pass_thd)

####
#### Plot pvalues
####
ggplot(pvals_sup_df, aes(x=comparison, y=minuslogpval, color=score)) +
  geom_point() +
  geom_label(data=pass_thd_sup, aes(x=comparison, y=yval), label="*", color="black", size=5) +
  geom_hline(yintercept=2) +
  scale_color_social_d() +
  facet_wrap(~deconv) +
  theme_lucid()
ggsave(paste0("3b_compare_pval/221101_",block,"_",feat_selec_meth,"_sup_pval_",tv_norm,".pdf"),
       width = 8, height = 5)
ggplot(pvals_unsup_df, aes(x=comparison, y=minuslogpval, color=score)) +
  geom_point() +
  geom_label(data=pass_thd_unsup, aes(x=comparison, y=yval), label="*", color="black", size=5) +
  geom_hline(yintercept=2) +
  scale_color_social_d() +
  facet_wrap(~deconv) +
  theme_lucid()
ggsave(paste0("3b_compare_pval/221101_",block,"_",feat_selec_meth,"_unsup_pval_",tv_norm,".pdf"),
       width = 8, height = 5)

####
#### Extract significant comparisons
####
pass_thd_sup %>% select(deconv,comparison)
table(pass_thd_sup %>% pull(comparison))
table(pass_thd_sup %>% pull(deconv))
pass_thd_unsup %>% select(deconv,comparison)
table(pass_thd_unsup %>% pull(comparison))
table(pass_thd_unsup %>% pull(deconv))
