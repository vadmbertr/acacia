## ----
## load scores
## ----
df_score1 <- readRDS("../2302_3MB_4multimap/4c_deconv/perf_scores/221101_scores.rds")
df_score2 <- readRDS("../2302_3MB_4multimap/4c_deconv/perf_scores/221101_time.rds")

df_bl111 <- readRDS("../2210_1SB/perf_scores/221101_scores_met_sup.rds")
df_bl111$type = 'met'
df_bl121 <- readRDS("../2210_1SB/perf_scores/221101_scores_rna_sup.rds")
df_bl121$type = 'rna'
df_bl112 <- readRDS("../2210_1SB/perf_scores/221101_scores_met_unsup.rds")
df_bl112$type = 'met'
df_bl122 <- readRDS("../2210_1SB/perf_scores/221101_scores_rna_unsup.rds")
df_bl122$type = 'rna'
df_bl1 <- rbind(df_bl111,df_bl121,df_bl112,df_bl122)
rm(df_bl111,df_bl121,df_bl112,df_bl122)
df_bl211 <- readRDS("../2210_1SB/perf_scores/221101_time_met_sup.rds")
df_bl211$type = 'met'
df_bl221 <- readRDS("../2210_1SB/perf_scores/221101_time_rna_sup.rds")
df_bl221$type = 'rna'
df_bl212 <- readRDS("../2210_1SB/perf_scores/221101_time_met_unsup.rds")
df_bl212$type = 'met'
df_bl222 <- readRDS("../2210_1SB/perf_scores/221101_time_rna_unsup.rds")
df_bl222$type = 'rna'
df_bl2 <- rbind(df_bl211,df_bl221,df_bl212,df_bl222)
rm(df_bl211,df_bl221,df_bl212,df_bl222)

df_toast111 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_met_toast_sup.rds")
df_toast111$type = 'met'
df_toast121 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_rna_toast_sup.rds")
df_toast121$type = 'rna'
df_toast112 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_met_toast_unsup.rds")
df_toast112$type = 'met'
df_toast122 <- readRDS("../2210_1SB_featselec/perf_scores/221101_scores_rna_toast_unsup.rds")
df_toast122$type = 'rna'
df_toast1 <- rbind(df_toast111,df_toast121,df_toast112,df_toast122)
rm(df_toast111,df_toast121,df_toast112,df_toast122)
df_toast211 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_met_toast_sup.rds")
df_toast211$type = 'met'
df_toast221 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_rna_toast_sup.rds")
df_toast221$type = 'rna'
df_toast212 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_met_toast_unsup.rds")
df_toast212$type = 'met'
df_toast222 <- readRDS("../2210_1SB_featselec/perf_scores/221101_time_rna_toast_unsup.rds")
df_toast222$type = 'rna'
df_toast2 <- rbind(df_toast211,df_toast221,df_toast212,df_toast222)
rm(df_toast211,df_toast221,df_toast212,df_toast222)

df_bl1$deconv[df_bl1$deconv=="CIBERSORT3"] <- "cibersort"
df_bl1$deconv[df_bl1$deconv=="CIBERSORT"] <- "cibersort"
df_bl1 <- df_bl1[df_bl1$deconv!="CIBERSORT4",]
df_bl1 <- df_bl1[df_bl1$deconv!="svr4",]
df_bl1$deconv[df_bl1$deconv=="RLR"] <- "rlr"
df_bl1$deconv[df_bl1$deconv=="DECONica"] <- "ica"
df_bl1$deconv[df_bl1$deconv=="DECONnmf"] <- "nmf"
df_bl1$deconv[df_bl1$deconv=="ICA"] <- "ica"
df_bl1 <- df_bl1[df_bl1$deconv!="svr4",]
df_bl1$deconv[df_bl1$deconv=="svr3"] <- "svr"
df_bl2$deconv[df_bl2$deconv=="CIBERSORT3"] <- "cibersort"
df_bl2$deconv[df_bl2$deconv=="CIBERSORT"] <- "cibersort"
df_bl2 <- df_bl2[df_bl2$deconv!="CIBERSORT4",]
df_bl2 <- df_bl2[df_bl2$deconv!="svr4",]
df_bl2$deconv[df_bl2$deconv=="RLR"] <- "rlr"
df_bl2$deconv[df_bl2$deconv=="DECONica"] <- "ica"
df_bl2$deconv[df_bl2$deconv=="DECONnmf"] <- "nmf"
df_bl2$deconv[df_bl2$deconv=="ICA"] <- "ica"
df_bl2 <- df_bl2[df_bl2$deconv!="svr4",]
df_bl2$deconv[df_bl2$deconv=="svr3"] <- "svr"

df_toast1$deconv[df_toast1$deconv=="CIBERSORT"] <- "cibersort"
df_toast1$deconv[df_toast1$deconv=="RLR"] <- "rlr"
df_toast1$deconv[df_toast1$deconv=="DECONica"] <- "ica"
df_toast1$deconv[df_toast1$deconv=="DECONnmf"] <- "nmf"
df_toast1$deconv[df_toast1$deconv=="ICA"] <- "ica"
df_toast2$deconv[df_toast2$deconv=="CIBERSORT"] <- "cibersort"
df_toast2$deconv[df_toast2$deconv=="RLR"] <- "rlr"
df_toast2$deconv[df_toast2$deconv=="DECONica"] <- "ica"
df_toast2$deconv[df_toast2$deconv=="DECONnmf"] <- "nmf"
df_toast2$deconv[df_toast2$deconv=="ICA"] <- "ica"

colnames(df_score1)[colnames(df_score1)=="score"] <- "scor"
colnames(df_bl1)[colnames(df_bl1)=="score"] <- "scor"
colnames(df_toast1)[colnames(df_toast1)=="score"] <- "scor"
df_score1 <- df_score1 %>%
  mutate(score=paste(scor,setting)) %>%
  select(c("values","score","sim","dataset","deconv","type")) %>%
  mutate(candidate=paste(deconv,type,"multimap")) %>%
  select(c('values','score','sim','dataset','candidate'))
df_bl1 <- df_bl1 %>%
  mutate(score=paste(scor,setting)) %>%
  select(c("values","score","sim","dataset","deconv","type")) %>%
  mutate(candidate=paste(deconv,type,"baseline")) %>%
  select(c('values','score','sim','dataset','candidate'))
df_toast1 <- df_toast1 %>%
  mutate(score=paste(scor,setting)) %>%
  select(c("values","score","sim","dataset","deconv","type")) %>%
  mutate(candidate=paste(deconv,type,"toast")) %>%
  select(c('values','score','sim','dataset','candidate'))
df_score2 <- df_score2 %>%
  mutate(candidate=paste(deconv,type,"multimap")) %>%
  select(c('values','score','sim','dataset','candidate'))
df_bl2 <- df_bl2 %>%
  mutate(candidate=paste(deconv,type,"baseline")) %>%
  select(c('values','score','sim','dataset','candidate'))
df_toast2 <- df_toast2 %>%
  mutate(candidate=paste(deconv,type,"toast")) %>%
  select(c('values','score','sim','dataset','candidate'))

name_data = sort(unique(df_bl1$dataset))
df_score1$dataset <- sapply(df_score1$dataset, function(x)
  switch(x, "dataset1"=name_data[1], "dataset2"=name_data[2], "dataset3"=name_data[3]))
df_score2$dataset <- sapply(df_score2$dataset, function(x)
  switch(x, "dataset1"=name_data[1], "dataset2"=name_data[2], "dataset3"=name_data[3]))

df_tot1 <- rbind(df_bl1,df_toast1,df_score1)
df_tot2 <- rbind(df_bl2,df_toast2[,colnames(df_bl2)],df_score2)
rm(df_bl1,df_toast1,df_score1,
   df_bl2,df_toast2,df_score2,
   name_data)