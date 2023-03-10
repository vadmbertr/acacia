library(RColorBrewer)
do_df_deconv = function(deconv_res, prop, methods) {
  plot_order = c("B cells", "Cancer basal", "Cancer classical", "CD4 T cells", "CD8 T cells",
                 "Endothelial", "Fibroblasts", "Macrophages", "Neutrophils")
  colourCount = length(plot_order)
  getPalette = colorRampPalette(brewer.pal(colourCount, "Set1"))
  names(deconv_res) = methods
  prop = prop[plot_order,]
  df_res = list()
  for (met in methods) {
    res = deconv_res[[met]]
    res = res[plot_order,]
    df_res[[met]] = reshape2::melt(res)
    colnames(df_res[[met]]) = c("comp","sample","res")
    df_ref = reshape2::melt(prop)
    df_res[[met]]$ref = df_ref$value
    df_res[[met]]$met = met
  }
  df_res = do.call(rbind,df_res)
  return(list(df_res,colourCount,getPalette))
}

rmse <-  function(real, prediction) {
  return(sqrt(mean((real - prediction)^2)))
}

mae <-  function(real, prediction) {
  return(mean(abs(real - prediction)))
}

pearson <-  function(real, prediction) {
  return(cor(c(real), c(prediction), method="pearson"))
}

score_perf <- function(real, prediction, method = c("rmse","mae","pearson")) {
  method <- match.arg(method)
  if (method=="rmse") {
    return(rmse(real, prediction))
  }
  else if (method=="mae") {
    return(mae(real, prediction))
  }
  else if (method=="pearson") {
    return(pearson(real, prediction))
  }
}
