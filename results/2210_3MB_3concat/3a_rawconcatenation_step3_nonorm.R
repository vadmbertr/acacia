## ----
## Set parameters
## ----
library(dplyr)
library(ggplot2)
library(see)
source("../../src/ranking_ranking_procedure_functions.R")
ranking_norm=ranking_2 # to vary, either ranking12 or ranking_2
ranking_type=ranking12345 # to vary, either ranking12345 or ranking12453
rank_norm=ifelse(all.equal(ranking_norm,ranking12)==T,"normDataset","normGlobal")
rank_type=ifelse(all.equal(ranking_type,ranking12345)==T,"order345","order453")

## ----
## load scores restricted and selected
## ----
input_path = paste0("../2210_0simu/simulations/rna/")
n_lot = length(sort(list.files(input_path, pattern = "T_rna")))
n_sims = length(sort(list.files(input_path, pattern = "sim")))/n_lot
name_data = unname(sapply(list.files(input_path, pattern = "T_rna"),function(x) rev(strsplit(strsplit(x,"T_rna")[[1]][1],"_")[[1]])[1]))
rm(input_path)

score_bl_met_sup <- readRDS("../2210_1SB/perf_scores/221101_scores_met_sup.rds")
score_bl_met_sup$datatype = "met"
score_bl_met_unsup <- readRDS("../2210_1SB/perf_scores/221101_scores_met_unsup.rds")
score_bl_met_unsup$datatype = "met"
score_bl_rna_sup <- readRDS("../2210_1SB/perf_scores/221101_scores_rna_sup.rds")
score_bl_rna_sup$datatype = "rna"
score_bl_rna_unsup <- readRDS("../2210_1SB/perf_scores/221101_scores_rna_unsup.rds")
score_bl_rna_unsup$datatype = "rna"
score_bl <- bind_rows(score_bl_met_sup,score_bl_met_unsup,
                      score_bl_rna_sup,score_bl_rna_unsup)
rm(score_bl_met_sup,score_bl_met_unsup,score_bl_rna_sup,score_bl_rna_unsup)

score_rawconcat <- readRDS("3a_rawconcatenation_step2_nonorm/perf_scores/221101_scores.rds")
score_rawconcat$datatype = "rawconcat_hvg"

score_all <- bind_rows(score_bl,score_rawconcat)
deconv_methods = sort(unique(score_rawconcat$deconv))
rm(score_bl,score_rawconcat)

## ----
## Load timings
## ----
time_bl_met_sup <- readRDS("../2210_1SB/perf_scores/221101_time_met_sup.rds")
time_bl_met_sup$datatype = "met"
time_bl_met_unsup <- readRDS("../2210_1SB/perf_scores/221101_time_met_unsup.rds")
time_bl_met_unsup$datatype = "met"
time_bl_rna_sup <- readRDS("../2210_1SB/perf_scores/221101_time_rna_sup.rds")
time_bl_rna_sup$datatype = "rna"
time_bl_rna_sup$values = as.numeric(time_bl_rna_sup$values)
time_bl_rna_unsup <- readRDS("../2210_1SB/perf_scores/221101_time_rna_unsup.rds")
time_bl_rna_unsup$datatype = "rna"
time_bl_rna_unsup$values = as.numeric(time_bl_rna_unsup$values)
time_bl <- bind_rows(time_bl_met_sup,time_bl_met_unsup,
                      time_bl_rna_sup,time_bl_rna_unsup)
rm(time_bl_met_sup,time_bl_met_unsup,time_bl_rna_sup,time_bl_rna_unsup)

time_rawconcat=lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    sapply(deconv_methods, function(meth) {
      if (file.exists(paste0("3a_rawconcatenation_step2_nonorm/deconv/221101_",name_data[lot],
                               "_Apred_",
                               meth,
                               "_",sim_txt,sim,".rds"))) {
          tmp=readRDS(paste0("3a_rawconcatenation_step2_nonorm/deconv/221101_",name_data[lot],
                             "_Apred_",
                             meth,
                             "_",sim_txt,sim,".rds"))$time_elapsed
          tmp0=readRDS(paste0("3a_rawconcatenation_step2_nonorm/deconv/221101_",name_data[lot],
                              "_Apred_",
                              meth,
                              "_",sim_txt,sim,".rds"))$res
          if (is.list(tmp)) {
            tmp=tmp[[4]]
            tmp=as.numeric(strsplit(tmp," ")[[1]][1])
            saveRDS(list(res=tmp0,time_elapsed=tmp),paste0("3a_rawconcatenation_step2_nonorm/deconv/221101_",name_data[lot],
                                                           "_Apred_",
                                                           meth,
                                                           "_",sim_txt,sim,".rds"))
          }
          return(tmp)
        }
        else {print(paste0(lot,sim,meth," not run yet"))}
      })
  })
})
time_rawconcat_df = do.call(rbind,lapply(seq_along(time_rawconcat), function (lot) {
  do.call(rbind,lapply(seq_along(time_rawconcat[[lot]]), function(sim) {
    df = data.frame("values"=time_rawconcat[[lot]][[sim]],
                    "score"="time",
                    "deconv"=names(time_rawconcat[[lot]][[sim]]),
                    "sim"=sim,
                    "dataset"=name_data[lot],
                    "datatype"="rawconcat_hvg")
  }))
}))

## ----
## Organise df
## ----
colnames(score_all)[colnames(score_all)=="score"] <- "scor"
score_all <- score_all %>%
  mutate(score=paste(scor,setting)) %>%
  select(c("values","score","sim","dataset","datatype","deconv")) %>%
  mutate(candidate=paste(datatype,deconv)) %>%
  select(c('values','score','sim','dataset','candidate'))

time_bl$sim <- as.numeric(time_bl$sim) 
time_all <- bind_rows(time_bl,time_rawconcat_df) %>%
  mutate(candidate=paste(datatype,deconv)) %>%
  select(c('values','score','sim','dataset','candidate'))
rm(time_rawconcat_df,time_rawconcat,time_bl)

score_all <- score_all[c(grep("CIBERSORT",score_all$candidate),
                         grep("rlr",score_all$candidate),
                         grep("ICA",score_all$candidate),
                         grep("CIBERSORT3",score_all$candidate),
                         grep("RLR",score_all$candidate),
                         grep("DECONica",score_all$candidate),
                         grep("cibersort",score_all$candidate)),]
score_all <- score_all[-grep("CIBERSORT4",score_all$candidate),]
time_all <- time_all[c(grep("CIBERSORT",time_all$candidate),
                         grep("rlr",time_all$candidate),
                         grep("ICA",time_all$candidate),
                         grep("CIBERSORT3",time_all$candidate),
                         grep("RLR",time_all$candidate),
                         grep("DECONica",time_all$candidate),
                         grep("cibersort",time_all$candidate)),]
time_all <- time_all[-grep("CIBERSORT4",time_all$candidate),]
score_all$candidate[score_all$candidate == "met CIBERSORT"] = "met cibersort"
score_all$candidate[score_all$candidate == "rna CIBERSORT3"] = "rna cibersort"
score_all$candidate[score_all$candidate == "met ICA"] = "met ica"
score_all$candidate[score_all$candidate == "rna RLR"] = "rna rlr"
score_all$candidate[score_all$candidate == "rna DECONica"] = "rna ica"
time_all$candidate[time_all$candidate == "met CIBERSORT"] = "met cibersort"
time_all$candidate[time_all$candidate == "rna CIBERSORT3"] = "rna cibersort"
time_all$candidate[time_all$candidate == "met ICA"] = "met ica"
time_all$candidate[time_all$candidate == "rna RLR"] = "rna rlr"
time_all$candidate[time_all$candidate == "rna DECONica"] = "rna ica"

## ----
## Rank
## ----
rank <- ranking_type(score_all,time_all, ranking=ranking_norm)

rank <- bind_rows(rank,rank[grep("rna",rank$candidate),],rank[grep("met",rank$candidate),])
rank$datatype <- sapply(rank$candidate, function(x) {
  inf = strsplit(x," ")[[1]][1]
  switch(inf,
         "rawconcat_hvg"="concat",
         "met"="MET",
         "rna"="RNA", stop("error"))
})
rank$datatype <- factor(rank$datatype, levels=c("MET","RNA","concat"))
rank$deconv <- sapply(rank$candidate, function(x) strsplit(x," ")[[1]][2])

ggplot(rank, aes(x=datatype, y=values, color=deconv)) +
  geom_point() +
  geom_line(aes(group=deconv)) +
  geom_line(aes(group=deconv), linewidth=5, alpha=.4) +
  scale_color_social_d() +
  ylab("Aggregated score") +
  xlab("") +
  theme_modern(axis.text.angle = 45)
ggsave(paste0("3a_rawconcatenation_step3_nonorm/221101_plot_",rank_type,"_",rank_norm,".pdf"),
       width=8, height=4)
