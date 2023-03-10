library(ggplot2)
library(see)

plot_ind_scores <- function(df, score, facet=T) {
  df_score <- df[df$score==score,]
  if (facet) {
  p <- ggplot(df_score, aes(x=deconv, y=values, fill=feat_selec, color=feat_selec)) +
    geom_point(position = position_jitterdodge(seed = 1, dodge.width = 1)) +
    geom_violin(alpha=.4, position = position_dodge(width = 0.5)) +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    facet_wrap(~facet) +
    theme_modern(axis.text.angle = 45)
  }
  else {
    p <- ggplot(df_score, aes(x=deconv, y=values, fill=feat_selec, color=feat_selec)) +
      geom_point(position = position_jitterdodge(seed = 1, dodge.width = 1)) +
      geom_violin(alpha=.4, position = position_dodge(width = 0.5)) +
      scale_fill_viridis_d() +
      scale_color_viridis_d() +
      theme_modern(axis.text.angle = 45)
  }
  p
}
