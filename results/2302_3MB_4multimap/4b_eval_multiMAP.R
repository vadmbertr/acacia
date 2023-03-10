library(ggplot2)
library(see)

## ----
## Load Data
## ----
list_adata_lot1 = list.files("4a_multiMAP/", pattern="dataset1")
list_adata_lot2 = list.files("4a_multiMAP/", pattern="dataset2")
list_adata_lot3 = list.files("4a_multiMAP/", pattern="dataset3")

adata_lot1 = lapply(list_adata_lot1, function(x) {
  file = anndata::read_h5ad(paste0("4a_multiMAP/",x))
  obs = file$obs
  obsm = file$obsm
  df = data.frame("MultiMAP1"=obsm$X_multimap[,1],
                  "MultiMAP2"=obsm$X_multimap[,2],
                  "Source"=obs$source,
                  "Type"=obs$type,
                  "SampleNb"=obs$SampleNb)
  df$SampleType = paste(df$Type,df$Source)
  df})
adata_lot2 = lapply(list_adata_lot2, function(x) {
  file = anndata::read_h5ad(paste0("4a_multiMAP/",x))
  obs = file$obs
  obsm = file$obsm
  df = data.frame("MultiMAP1"=obsm$X_multimap[,1],
                  "MultiMAP2"=obsm$X_multimap[,2],
                  "Source"=obs$source,
                  "Type"=obs$type,
                  "SampleNb"=obs$SampleNb)
  df$SampleType = paste(df$Type,df$Source)
  df})
adata_lot3 = lapply(list_adata_lot3, function(x) {
  file = anndata::read_h5ad(paste0("4a_multiMAP/",x))
  obs = file$obs
  obsm = file$obsm
  df = data.frame("MultiMAP1"=obsm$X_multimap[,1],
                  "MultiMAP2"=obsm$X_multimap[,2],
                  "Source"=obs$source,
                  "Type"=obs$type,
                  "SampleNb"=obs$SampleNb)
  df$SampleType = paste(df$Type,df$Source)
  df})
rm(list_adata_lot1,list_adata_lot2,list_adata_lot3)

## ----
## Plot
## ----
plot_lot1 = lapply(adata_lot1, function(x)
  ggplot(x, aes(x=MultiMAP1, y=MultiMAP2, color=SampleType, group=SampleNb)) +
    geom_point() +
    geom_line(alpha=.3, color="gray") +
    scale_color_social_d() +
    theme_modern() +
    ggpubr::rremove("xylab"))
plot_lot2 = lapply(adata_lot2, function(x)
  ggplot(x, aes(x=MultiMAP1, y=MultiMAP2, color=SampleType, group=SampleNb)) +
    geom_point() +
    geom_line(alpha=.3, color="gray") +
    scale_color_social_d() +
    theme_modern() +
    ggpubr::rremove("xylab"))
plot_lot3 = lapply(adata_lot3, function(x)
  ggplot(x, aes(x=MultiMAP1, y=MultiMAP2, color=SampleType, group=SampleNb)) +
    geom_point() +
    geom_line(alpha=.3, color="gray") +
    scale_color_social_d() +
    theme_modern() +
    ggpubr::rremove("xylab"))
plot_lot1 <- ggpubr::ggarrange(plotlist=plot_lot1, nrow=2, ncol=5, common.legend = T)
plot_lot2 <- ggpubr::ggarrange(plotlist=plot_lot2, nrow=2, ncol=5, common.legend = T)
plot_lot3 <- ggpubr::ggarrange(plotlist=plot_lot3, nrow=2, ncol=5, common.legend = T)
annotate_figure(plot_lot1,
                left = grid::textGrob("MultiMAP 2", rot = 90, vjust = 1, gp = grid::gpar(cex = 1.3)),
                bottom = grid::textGrob("MultiMAP 1", gp = grid::gpar(cex = 1.3)))
ggsave("4b_eval_multiMAP/dataset1.pdf", width=12, height=5)
annotate_figure(plot_lot2,
                left = grid::textGrob("MultiMAP 2", rot = 90, vjust = 1, gp = grid::gpar(cex = 1.3)),
                bottom = grid::textGrob("MultiMAP 1", gp = grid::gpar(cex = 1.3)))
ggsave("4b_eval_multiMAP/dataset2.pdf", width=12, height=5)
annotate_figure(plot_lot3,
                left = grid::textGrob("MultiMAP 2", rot = 90, vjust = 1, gp = grid::gpar(cex = 1.3)),
                bottom = grid::textGrob("MultiMAP 1", gp = grid::gpar(cex = 1.3)))
ggsave("4b_eval_multiMAP/dataset3.pdf", width=12, height=5)
