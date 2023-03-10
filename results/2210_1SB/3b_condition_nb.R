library(dplyr)

## ----
## Set parameters
## ----
block="rna" # to vary, either met or rna (+ line 14)

## ----
## load simulated D data, along with A_true and T_true
## ----
input_path = paste0("../2210_0simu/simulations/",block,"/")
input_path_matrix = sort(list.files(input_path, pattern = "sim"))
input_path_T = sort(list.files(input_path, pattern = paste0("T_",block,"_ref")))
matrix = lapply(input_path_matrix, function(x) readRDS(paste0(input_path,x))$D_rna_sim) # to vary, either met or rna
T_ref = lapply(input_path_T, function(x) readRDS(paste0(input_path,x)))
T_ref = lapply(T_ref, as.data.frame)
names(matrix) <- input_path_matrix
names(T_ref) <- input_path_T

n_sims=length(matrix)/length(T_ref)
n_lot=length(T_ref)
names_datasets=unname(sapply(input_path_T, function(x)
  rev(strsplit(strsplit(x,paste0("_T_",block,"_ref"))[[1]][1],"_")[[1]])[1]))
n_genes=unname(sapply(matrix, function(x) nrow(x)))
n_cell_types=unname(sapply(T_ref, ncol))

input_path_methodssup = paste0("timing/",block,"/sup/")
input_path_methodssup = sort(list.files(input_path_methodssup))
deconv_methods_sup = unique(sapply(input_path_methodssup, function(x)
  rev(strsplit(x,"_")[[1]])[2]))
deconv_methods_sup[deconv_methods_sup=="net"] = "elastic_net"
rm(input_path_methodssup)

get_CN <- function(mat) {
  norm(pracma::pinv(mat),"I")*norm(mat,"I")
}

## ----
## load supervised timings and get mean for each sim
## ----
timing=sapply(seq(n_lot), function (lot) {
  print(paste0("Loading dataset n°",lot,"/",n_lot))
  tmp=sapply(seq(n_sims), function(sim) {
    sim_txt <- ifelse(length(strsplit(as.character(sim),"")[[1]])==1,'sim0','sim')
    tmp=sapply(deconv_methods_sup, function(meth) {
      if (file.exists(paste0("timing/",block,"/sup/",strsplit(input_path_T[lot],paste0("_T_",block,"_ref"))[[1]][1],"_timing_",meth,"_",sim_txt,sim,".rds"))) {
        readRDS(paste0("timing/",block,"/sup/",strsplit(input_path_T[lot],paste0("_T_",block,"_ref"))[[1]][1],"_timing_",meth,"_",sim_txt,sim,".rds"))
      }
      else {print('missing timing')}
    })
    names(tmp)=deconv_methods_sup;tmp
  })
  col_names=seq(n_sims)
  row_names=rownames(tmp)
  tmp=c(tmp)
  names(tmp)=c(sapply(col_names,function(x,y) paste(x,y), row_names))
  tmp
})
timing_df0 = data.frame("Timing"=c(timing),
                        "Dataset"=rep(names_datasets,each=n_sims*length(deconv_methods_sup)),
                        "Method"=deconv_methods_sup,
                        "Simulation"=c(sapply(LETTERS[seq(n_lot)], function(x)
                          rep(paste0(x,seq(n_sims)), each=length(deconv_methods_sup)))))
timing_df0$datasetXmethod <- paste(timing_df0$Dataset,timing_df0$Method)

timing_df1 = timing_df0 %>%
  group_by(Simulation) %>%
  summarise("TimeMean"=mean(as.numeric(Timing)))
timing_df1$Simulation = factor(timing_df1$Simulation,
                               levels = c(sapply(LETTERS[seq(n_lot)], function(x,y) paste0(x,y), seq(n_sims))))
timing <- timing_df1[order(timing_df1$Simulation),] %>%
  pull(TimeMean)
timingT <- timing_df1[order(timing_df1$Simulation),] %>%
  mutate("Dataset"=rep(seq(n_lot),each=n_sims)) %>%
  group_by(Dataset) %>%
  summarise("TimeMeanT"=mean(TimeMean)) %>%
  pull(TimeMeanT)

timing_df2 = timing_df0 %>%
  group_by(datasetXmethod) %>%
  summarise("TimeMean" = mean(as.numeric(Timing)))
timing_df2$Method <- sapply(timing_df2$datasetXmethod, function(x)
  strsplit(x," ")[[1]][2])
timing_df2$Dataset <- sapply(timing_df2$datasetXmethod, function(x)
  strsplit(x," ")[[1]][1])
timing2 <- timing_df2 %>% pull(TimeMean)

## ----
## Compute D condition nb
## ----
CN_D <- rep(NA,length(matrix))
for (mat in seq(length(matrix))) {
  print(mat/length(matrix))
  if (is.na(CN_D[mat])) {
    if (nrow(matrix[[mat]])>34883) {
      hvg=TOAST::findRefinx(matrix[[mat]], nmarker = 34883)
      matrix[[mat]]=matrix[[mat]][hvg,]
    }
    cn <- get_CN(matrix[[mat]])
    matrix[[mat]] <- NA
    CN_D[mat] <- cn
  }
}

## ----
## Compute T condition nb
## ----
CN_T <- pbapply::pbsapply(T_ref, function(x) get_CN(as.matrix(x)))

## ----
## Plot D & T condition nb
## ----
library(ggplot2)
library(see)
CN_df <- data.frame("CN"=c(CN_D,CN_T),
                    "Simulation"=c(rep(seq(n_sims),n_lot),rep("NA",n_lot)),
                    "Dataset"=c(rep(names_datasets,each=n_sims),names_datasets),
                    "Matrix"=c(rep("D",n_sims*n_lot),rep("T",n_lot)),
                    "Ngenes"=c(n_genes,unique(n_genes)),
                    "Ncelltypes"=c(rep(n_cell_types,each=n_sims),n_cell_types),
                    "CN_T"=c(rep(CN_T,each=n_sims),CN_T),
                    "Timing"=c(timing,timingT))
ggplot(CN_df, aes(x=Dataset, y=CN, color=Simulation)) +
  geom_point() +
  scale_y_log10() +
  scale_color_flat_d() +
  facet_wrap(~Matrix) +
  theme_modern()
ggsave(paste0("condition_nb/",
              strsplit(input_path_T[1],"_")[[1]][1],
              "_",block,
              "_CN_D.pdf"),
       width = 7, height = 4)
ggplot(CN_df[CN_df$Matrix=="D",], aes(x=Ngenes, y=CN, color=Dataset)) +
  geom_point() +
  scale_y_log10() +
  scale_color_viridis_d() +
  theme_modern()
ggsave(paste0("condition_nb/",
              strsplit(input_path_T[1],"_")[[1]][1],
              "_",block,
              "_CN_D_Ngenes.pdf"),
       width = 7, height = 4)
ggplot(CN_df[CN_df$Matrix=="D",], aes(x=Ncelltypes, y=CN, color=Dataset)) +
  geom_point() +
  scale_y_log10() +
  scale_color_viridis_d() +
  theme_modern()
ggsave(paste0("condition_nb/",
              strsplit(input_path_T[1],"_")[[1]][1],
              "_",block,
              "_CN_D_Ncelltypes.pdf"),
       width = 7, height = 4)
ggplot(CN_df[CN_df$Matrix=="D",], aes(x=CN_T, y=CN, color=Dataset)) +
  geom_point() +
  scale_y_log10() +
  scale_color_viridis_d() +
  theme_modern()
ggsave(paste0("condition_nb/",
              strsplit(input_path_T[1],"_")[[1]][1],
              "_",block,
              "_CN_D_T.pdf"),
       width = 7, height = 4)
ggplot(CN_df[CN_df$Matrix=="T",], aes(x=Ngenes, y=CN, color=Dataset)) +
  geom_point() +
  scale_y_log10() +
  scale_color_viridis_d() +
  theme_modern()
ggsave(paste0("condition_nb/",
              strsplit(input_path_T[1],"_")[[1]][1],
              "_",block,
              "_CN_T_Ngenes.pdf"),
       width = 7, height = 4)
ggplot(CN_df[CN_df$Matrix=="T",], aes(x=Ncelltypes, y=CN, color=Dataset)) +
  geom_point() +
  scale_y_log10() +
  scale_color_viridis_d() +
  theme_modern()
ggsave(paste0("condition_nb/",
              strsplit(input_path_T[1],"_")[[1]][1],
              "_",block,
              "_CN_T_Ncelltypes.pdf"),
       width = 7, height = 4)
ggplot(CN_df[CN_df$Matrix=="D",], aes(x=Timing, y=CN, color=Dataset)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_viridis_d() +
  theme_modern()
ggsave(paste0("condition_nb/",
              strsplit(input_path_T[1],"_")[[1]][1],
              "_",block,
              "_CN_D_timing.pdf"),
       width = 7, height = 4)
ggplot(CN_df[CN_df$Matrix=="T",], aes(x=Timing, y=CN, color=Dataset)) +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_viridis_d() +
  theme_modern()
ggsave(paste0("condition_nb/",
              strsplit(input_path_T[1],"_")[[1]][1],
              "_",block,
              "_CN_T_timing.pdf"),
       width = 7, height = 4)

## ----
## Plot D & T condition nb per method
## ----
CN_df2 <- data.frame("CN"=rep(CN_D,each=length(deconv_methods_sup)),
                     "Simulation"=factor(sapply(timing_df0$Simulation, function(x)
                       paste0(strsplit(x,"")[[1]][-1], collapse="")), levels=seq(n_sims)),
                     "Dataset"=timing_df0$Dataset,
                     "Timing"=as.numeric(timing_df0$Timing),
                     "Method"=timing_df0$Method)
ggplot(CN_df2, aes(x=Timing, y=CN, color=Dataset)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_flat_d() +
  facet_wrap(~Method, scales = "free_x") +
  theme_modern()
ggsave(paste0("condition_nb/",
              strsplit(input_path_T[1],"_")[[1]][1],
              "_",block,
              "_CN_D_timing_facet_method.pdf"),
       width = 7, height = 4)

