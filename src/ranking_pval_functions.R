# Get trendval with all sim scores
trendvals12 <- function(scores1, scores2, normalization="centerscale") {
  
  if (any(colnames(scores1)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores1 input"))
  }
  if (any(colnames(scores2)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores2 input"))
  }
  if (normalization!="centerscale") {
    stop(print("Normalization should be 'centerscale'"))
  }
  Norm = function(x) {CenterScaleNorm(x)}
  
  scores1$sim <- as.numeric(scores1$sim)
  scores2$sim <- as.numeric(scores2$sim)
  scores1$values <- as.numeric(scores1$values)
  scores2$values <- as.numeric(scores2$values)
  
  score_df <- bind_rows(scores1,scores2)
  colnames(score_df)[1:5] <- c("score","name_score","simulation","dataset","candidate")
  
  # keep only non redundant scores
  score_to_keep <- c("time","rmse perf_g","mae perf_g","pearson perf_g","pearson med_c","pearson med_s","pearson sd_c","pearson sd_s")
  score_df <- score_df[score_df$name_score %in% score_to_keep,]
  
  # normalize for each dataset and score (step 1)
  df_centerscale <- score_df %>%
    group_by(dataset, name_score) %>%
    mutate(normval=Norm(score))
  rm(score_df)
  
  # transform scores such that 1 is the best score (step 2)
  df_centerscale_transfo <- df_centerscale
  df_centerscale_transfo$trendval <- 1 - df_centerscale_transfo$normval
  df_centerscale_transfo$trendval[grep("pearson m",df_centerscale_transfo$name_score)] <- 1 - df_centerscale_transfo$trendval[grep("pearson m",df_centerscale_transfo$name_score)]
  df_centerscale_transfo$trendval[grep("pearson p",df_centerscale_transfo$name_score)] <- 1 - df_centerscale_transfo$trendval[grep("pearson p",df_centerscale_transfo$name_score)]
  rm(df_centerscale)
  
  return(df_centerscale_transfo)
}
trendvals12_notime <- function(scores1, normalization="centerscale") {
  
  if (any(colnames(scores1)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores1 input"))
  }
  if (normalization!="centerscale") {
    stop(print("Normalization should be 'centerscale'"))
  }
  Norm = function(x) {CenterScaleNorm(x)}
  
  scores1$sim <- as.numeric(scores1$sim)
  scores1$values <- as.numeric(scores1$values)

  score_df <- scores1
  colnames(score_df)[1:5] <- c("score","name_score","simulation","dataset","candidate")
  
  # keep only non redundant scores
  score_to_keep <- c("rmse perf_g","mae perf_g","pearson perf_g","pearson med_c","pearson med_s","pearson sd_c","pearson sd_s")
  score_df <- score_df[score_df$name_score %in% score_to_keep,]
  
  # normalize for each dataset and score (step 1)
  df_centerscale <- score_df %>%
    group_by(dataset, name_score) %>%
    mutate(normval=Norm(score))
  rm(score_df)
  
  # transform scores such that 1 is the best score (step 2)
  df_centerscale_transfo <- df_centerscale
  df_centerscale_transfo$trendval <- 1 - df_centerscale_transfo$normval
  df_centerscale_transfo$trendval[grep("pearson m",df_centerscale_transfo$name_score)] <- 1 - df_centerscale_transfo$trendval[grep("pearson m",df_centerscale_transfo$name_score)]
  df_centerscale_transfo$trendval[grep("pearson p",df_centerscale_transfo$name_score)] <- 1 - df_centerscale_transfo$trendval[grep("pearson p",df_centerscale_transfo$name_score)]
  rm(df_centerscale)
  
  return(df_centerscale_transfo)
}

trendvals_2 <- function(scores1, scores2, normalization="centerscale") {
  
  if (any(colnames(scores1)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores1 input"))
  }
  if (any(colnames(scores2)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores2 input"))
  }
  if (normalization!="centerscale") {
    stop(print("Normalization should be 'centerscale'"))
  }
  Norm = function(x) {CenterScaleNorm(x)}
  
  scores1$sim <- as.numeric(scores1$sim)
  scores2$sim <- as.numeric(scores2$sim)
  scores1$values <- as.numeric(scores1$values)
  scores2$values <- as.numeric(scores2$values)
  
  score_df <- bind_rows(scores1,scores2)
  colnames(score_df)[1:5] <- c("score","name_score","simulation","dataset","candidate")
  
  # keep only non redundant scores
  score_to_keep <- c("time","rmse perf_g","mae perf_g","pearson perf_g","pearson med_c","pearson med_s","pearson sd_c","pearson sd_s")
  score_df <- score_df[score_df$name_score %in% score_to_keep,]
  
  # normalize for each dataset and score (step 1)
  df_centerscale <- score_df %>%
    ungroup() %>%
    mutate(normval=Norm(score))
  rm(score_df)
  
  # transform scores such that 1 is the best score (step 2)
  df_centerscale_transfo <- df_centerscale
  df_centerscale_transfo$trendval <- 1 - df_centerscale_transfo$normval
  df_centerscale_transfo$trendval[grep("pearson m",df_centerscale_transfo$name_score)] <- 1 - df_centerscale_transfo$trendval[grep("pearson m",df_centerscale_transfo$name_score)]
  df_centerscale_transfo$trendval[grep("pearson p",df_centerscale_transfo$name_score)] <- 1 - df_centerscale_transfo$trendval[grep("pearson p",df_centerscale_transfo$name_score)]
  rm(df_centerscale)
  
  return(df_centerscale_transfo)
}
trendvals_2_notime <- function(scores1, normalization="centerscale") {
  
  if (any(colnames(scores1)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores1 input"))
  }
  if (normalization!="centerscale") {
    stop(print("Normalization should be 'centerscale'"))
  }
  Norm = function(x) {CenterScaleNorm(x)}
  
  scores1$sim <- as.numeric(scores1$sim)
  scores1$values <- as.numeric(scores1$values)
  
  score_df <- scores1
  colnames(score_df)[1:5] <- c("score","name_score","simulation","dataset","candidate")
  
  # keep only non redundant scores
  score_to_keep <- c("rmse perf_g","mae perf_g","pearson perf_g","pearson med_c","pearson med_s","pearson sd_c","pearson sd_s")
  score_df <- score_df[score_df$name_score %in% score_to_keep,]
  
  # normalize for each dataset and score (step 1)
  df_centerscale <- score_df %>%
    ungroup() %>%
    mutate(normval=Norm(score))
  rm(score_df)
  
  # transform scores such that 1 is the best score (step 2)
  df_centerscale_transfo <- df_centerscale
  df_centerscale_transfo$trendval <- 1 - df_centerscale_transfo$normval
  df_centerscale_transfo$trendval[grep("pearson m",df_centerscale_transfo$name_score)] <- 1 - df_centerscale_transfo$trendval[grep("pearson m",df_centerscale_transfo$name_score)]
  df_centerscale_transfo$trendval[grep("pearson p",df_centerscale_transfo$name_score)] <- 1 - df_centerscale_transfo$trendval[grep("pearson p",df_centerscale_transfo$name_score)]
  rm(df_centerscale)
  
  return(df_centerscale_transfo)
}

coerce_pearson <- function(df_median_centerscale_transfo) {
  df_median_centerscale_transfo1 <- df_median_centerscale_transfo %>%
    mutate(typescore=sapply(name_score, function(x) strsplit(x, " ")[[1]][1])) %>%
    filter(typescore!="pearson") %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo21 <- df_median_centerscale_transfo %>%
    filter(name_score=="pearson med_c") %>%
    group_by(dataset,candidate,simulation) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson med_c",
           score = NA,
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo22 <- df_median_centerscale_transfo %>%
    filter(name_score=="pearson med_s") %>%
    group_by(dataset,candidate,simulation) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson med_s",
           score = NA,
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo23 <- df_median_centerscale_transfo %>%
    filter(name_score=="pearson perf_g") %>%
    group_by(dataset,candidate,simulation) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson perf_g",
           score = NA,
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo2 <- bind_rows(df_median_centerscale_transfo21,
                                              df_median_centerscale_transfo22,
                                              df_median_centerscale_transfo23) %>%
    group_by(dataset,candidate,simulation) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson perf_mean",
           score = NA,
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo31 <- df_median_centerscale_transfo %>%
    filter(name_score=="pearson sd_c") %>%
    group_by(dataset,candidate,simulation) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson sd_c",
           score = NA,
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo32 <- df_median_centerscale_transfo %>%
    filter(name_score=="pearson sd_s") %>%
    group_by(dataset,candidate,simulation) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson sd_s",
           score = NA,
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo3 <- bind_rows(df_median_centerscale_transfo31,
                                              df_median_centerscale_transfo32) %>%
    group_by(dataset,candidate,simulation) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson sd_mean",
           score = NA,
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo = bind_rows(df_median_centerscale_transfo1,
                                            df_median_centerscale_transfo2,
                                            df_median_centerscale_transfo3)
  return(df_median_centerscale_transfo)
}

# Compare p-val
pval <- function(before_perf, after_perf) {
  t.test(after_perf,before_perf,alternative="g",paired=T)$p.value
}
