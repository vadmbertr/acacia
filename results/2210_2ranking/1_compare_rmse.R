## ----
## Set parameters
## ----
library(ggplot2)
library(see)
library(dplyr)

## ----
## load RMSE
## ----
input_path = "../2210_0simu/simulations/rna/"
input_path_T = sort(list.files(input_path, pattern = 'rna'))
n_lot = length(input_path_T)
input_path_nsims = sort(list.files(input_path, pattern = "sim"))
n_sims = length(input_path_nsims)/n_lot
name_data = unname(sapply(input_path_T,function(x) rev(strsplit(strsplit(x,'_T_rna')[[1]][1],"_")[[1]])[1]))
date = strsplit(input_path_T[1],'_')[[1]][1]
rm(input_path,input_path_nsims,input_path_T)

RMSE_rna_sup = readRDS(paste0("../2210_1SB/perf_scores/",date,"_scores_rna_sup.rds"))
RMSE_rna_unsup = readRDS(paste0("../2210_1SB/perf_scores/",date,"_scores_rna_unsup.rds"))
RMSE_met_sup = readRDS(paste0("../2210_1SB/perf_scores/",date,"_scores_met_sup.rds"))
RMSE_met_unsup = readRDS(paste0("../2210_1SB/perf_scores/",date,"_scores_met_unsup.rds"))
RMSE_rna_sup <- RMSE_rna_sup[RMSE_rna_sup$score=="rmse",]
RMSE_rna_unsup <- RMSE_rna_unsup[RMSE_rna_unsup$score=="rmse",]
RMSE_met_sup <- RMSE_met_sup[RMSE_met_sup$score=="rmse",]
RMSE_met_unsup <- RMSE_met_unsup[RMSE_met_unsup$score=="rmse",]

## ----
## Compute RMSE sd
## ----
sd_rna_sup_g <- RMSE_rna_sup %>%
  filter(setting=="perf_g") %>%
  group_by(dataset,deconv) %>%
  summarise("sd"=sd(values)) %>%
  mutate("block"="rna", "type"="sup", "setting"="perf_g")
sd_rna_unsup_g <- RMSE_rna_unsup %>%
  filter(setting=="perf_g") %>%
  group_by(dataset,deconv) %>%
  summarise("sd"=sd(values)) %>%
  mutate("block"="rna", "type"="unsup", "setting"="perf_g")
sd_met_sup_g <- RMSE_met_sup %>%
  filter(setting=="perf_g") %>%
  group_by(dataset,deconv) %>%
  summarise("sd"=sd(values)) %>%
  mutate("block"="met", "type"="sup", "setting"="perf_g")
sd_met_unsup_g <- RMSE_met_unsup %>%
  filter(setting=="perf_g") %>%
  group_by(dataset,deconv) %>%
  summarise("sd"=sd(values)) %>%
  mutate("block"="met", "type"="unsup", "setting"="perf_g")

sd_rna_sup_c <- RMSE_rna_sup[grep("perf_c",RMSE_rna_sup$setting),] %>%
  group_by(dataset,deconv,setting) %>%
  summarise("sd"=sd(values)) %>%
  mutate("block"="rna", "type"="sup")
sd_rna_unsup_c <- RMSE_rna_unsup[grep("perf_c",RMSE_rna_unsup$setting),] %>%
  group_by(dataset,deconv,setting) %>%
  summarise("sd"=sd(values)) %>%
  mutate("block"="rna", "type"="unsup")
sd_met_sup_c <- RMSE_met_sup[grep("perf_c",RMSE_met_sup$setting),] %>%
  group_by(dataset,deconv,setting) %>%
  summarise("sd"=sd(values)) %>%
  mutate("block"="met", "type"="sup")
sd_met_unsup_c <- RMSE_met_unsup[grep("perf_c",RMSE_met_unsup$setting),] %>%
  group_by(dataset,deconv,setting) %>%
  summarise("sd"=sd(values)) %>%
  mutate("block"="met", "type"="unsup")

## ----
## Compute RMSE median
## ----
med_rna_sup_g <- RMSE_rna_sup %>%
  filter(setting=="perf_g") %>%
  group_by(dataset,deconv) %>%
  summarise("med"=median(values)) %>%
  mutate("block"="rna", "type"="sup", "setting"="perf_g")
med_rna_unsup_g <- RMSE_rna_unsup %>%
  filter(setting=="perf_g") %>%
  group_by(dataset,deconv) %>%
  summarise("med"=median(values)) %>%
  mutate("block"="rna", "type"="unsup", "setting"="perf_g")
med_met_sup_g <- RMSE_met_sup %>%
  filter(setting=="perf_g") %>%
  group_by(dataset,deconv) %>%
  summarise("med"=median(values)) %>%
  mutate("block"="met", "type"="sup", "setting"="perf_g")
med_met_unsup_g <- RMSE_met_unsup %>%
  filter(setting=="perf_g") %>%
  group_by(dataset,deconv) %>%
  summarise("med"=median(values)) %>%
  mutate("block"="met", "type"="unsup", "setting"="perf_g")

med_rna_sup_c <- RMSE_rna_sup[grep("perf_c",RMSE_rna_sup$setting),] %>%
  group_by(dataset,deconv,setting) %>%
  summarise("med"=median(values)) %>%
  mutate("block"="rna", "type"="sup")
med_rna_unsup_c <- RMSE_rna_unsup[grep("perf_c",RMSE_rna_unsup$setting),] %>%
  group_by(dataset,deconv,setting) %>%
  summarise("med"=median(values)) %>%
  mutate("block"="rna", "type"="unsup")
med_met_sup_c <- RMSE_met_sup[grep("perf_c",RMSE_met_sup$setting),] %>%
  group_by(dataset,deconv,setting) %>%
  summarise("med"=median(values)) %>%
  mutate("block"="met", "type"="sup")
med_met_unsup_c <- RMSE_met_unsup[grep("perf_c",RMSE_met_unsup$setting),] %>%
  group_by(dataset,deconv,setting) %>%
  summarise("med"=median(values)) %>%
  mutate("block"="met", "type"="unsup")

## ----
## Plot SD
## ----
df_sd <- bind_rows(sd_rna_sup_g,sd_rna_sup_c,sd_rna_unsup_g,sd_rna_unsup_c,
                   sd_met_sup_g,sd_met_sup_c,sd_met_unsup_g,sd_met_unsup_c)
ggplot(df_sd[df_sd$block=="rna",], aes(x=setting, y=sd, color=deconv)) +
  geom_point(aes(shape=type, size=type)) +
  geom_line(aes(group=deconv, size=type)) +
  scale_color_viridis_d(option = "B") +
  scale_size_discrete(range = c(.5,1.5)) +
  facet_wrap(~dataset, scales="free_x") +
  ylab("RMSE standard deviation across simulations") +
  theme_modern(axis.text.angle = 45)
ggsave("1_compare_rmse/compare_rmse_sd_rna.pdf",
       width=14, height=8)
ggplot(df_sd[df_sd$block=="met",], aes(x=setting, y=sd, color=deconv)) +
  geom_point(aes(shape=type, size=type)) +
  geom_line(aes(group=deconv, size=type)) +
  scale_color_viridis_d(option = "B") +
  scale_size_discrete(range = c(.5,1.5)) +
  facet_wrap(~dataset, scales="free_x") +
  ylab("RMSE standard deviation across simulations") +
  theme_modern(axis.text.angle = 45)
ggsave("1_compare_rmse/compare_rmse_sd_met.pdf",
       width=14, height=8)

## ----
## Plot VALUE
## ----
df_med <- bind_rows(med_rna_sup_g,med_rna_sup_c,med_rna_unsup_g,med_rna_unsup_c,
                    med_met_sup_g,med_met_sup_c,med_met_unsup_g,med_met_unsup_c)
ggplot(df_med[df_med$block=="rna",], aes(x=setting, y=med, color=deconv)) +
  geom_point(aes(shape=type, size=type)) +
  geom_line(aes(group=deconv, size=type)) +
  scale_color_viridis_d(option = "B") +
  scale_size_discrete(range = c(.5,1.5)) +
  facet_wrap(~dataset, scales="free_x") +
  ylab("RMSE median across simulations") +
  theme_modern(axis.text.angle = 45)
ggsave("1_compare_rmse/compare_rmse_med_rna.pdf",
       width=14, height=8)
ggplot(df_med[df_med$block=="met",], aes(x=setting, y=med, color=deconv)) +
  geom_point(aes(shape=type, size=type)) +
  geom_line(aes(group=deconv, size=type)) +
  scale_color_viridis_d(option = "B") +
  scale_size_discrete(range = c(.5,1.5)) +
  facet_wrap(~dataset, scales="free_x") +
  ylab("RMSE median across simulations") +
  theme_modern(axis.text.angle = 45)
ggsave("1_compare_rmse/compare_rmse_med_met.pdf",
       width=14, height=8)

## ----
## Compute & Plot RMSE median sd
## ----
grossmed1 <- df_med %>%
  filter(setting!="perf_g") %>%
  group_by(dataset, deconv, block, type) %>%
  summarise("sdmed"=sd(med)) %>%
  mutate("setting"="perf_c")
grossmed2 <- df_med %>%
  group_by(dataset, deconv, block, type) %>%
  summarise("sdmed"=sd(med)) %>%
  mutate("setting"="perf_all")
grossmed = bind_rows(grossmed1,grossmed2)
ggplot(grossmed[grossmed$block=="rna",], aes(x=setting, y=sdmed, color=deconv)) +
  geom_point(aes(shape=type, size=type)) +
  geom_line(aes(group=deconv, size=type)) +
  scale_color_viridis_d(option = "B") +
  scale_size_discrete(range = c(.5,1.5)) +
  facet_wrap(~dataset, scales="free_x") +
  ylab("RMSE standard deviation across datasets") +
  theme_modern(axis.text.angle = 45)
ggsave("1_compare_rmse/compare_rmse_sdmed_rna.pdf",
       width=14, height=8)
ggplot(grossmed[grossmed$block=="met",], aes(x=setting, y=sdmed, color=deconv)) +
  geom_point(aes(shape=type, size=type)) +
  geom_line(aes(group=deconv, size=type)) +
  scale_color_viridis_d(option = "B") +
  scale_size_discrete(range = c(.5,1.5)) +
  facet_wrap(~dataset, scales="free_x") +
  ylab("RMSE standard deviation across datasets") +
  theme_modern(axis.text.angle = 45)
ggsave("1_compare_rmse/compare_rmse_sdmed_met.pdf",
       width=14, height=8)

