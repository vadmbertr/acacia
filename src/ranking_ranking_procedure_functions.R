#######
####### BASIC FUNCTIONS ####### 
#######
require(dplyr, quietly = T)
library(ggplot2)
library(see)
options(dplyr.summarise.inform = FALSE)

geomMean <- function(x) {
  exp(mean(log(x), na.rm=T))
}
harmMean <- function(x) {
  length(x)/sum(x^(-1), na.rm=T)
}
weighMean <- function(x,weights) {
  weighted.mean(x, weights)
}

MinMaxNorm <-function(x) {
  (x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T))
}
CenterScaleNorm <-function(x) {
  #center 
  tr1 = x - mean(x, na.rm=T)
  #scale
  tr2 = tr1/sd(x, na.rm=T)
  return(pnorm(tr2))
}

convert_to_cat <- function(dict, vect) {
  tmp = rep(names(dict), lengths(dict))
  cat <- sapply(vect, function(x) {
    ifelse(x %in% unlist(dict), tmp[unlist(dict)==x], NA)
  })
  return(unname(unlist(cat)))
}

euclidean_vector_dist <- function(vec1,vec2) {
  vec1 <- vec1[(!is.na(vec1)) & (!is.na(vec2))]
  vec2 <- vec2[(!is.na(vec1)) & (!is.na(vec2))]
  sqrt(sum((vec1 - vec2)^2))
}

#######
####### RANK FUNCTION ####### 
#######

# input scores1 is the df with all scores except time
# input scores2 is the df with time scores
eval_norm_order_ <- function(scores_norm) {
  order_conserved <- all(scores_norm %>%
                           group_by(name_score,dataset) %>%
                           mutate(order_before=order(val), order_after=order(normval)) %>%
                           mutate(rank_conserved=order_before==order_after) %>%
                           pull(rank_conserved))
  return(order_conserved)
}

eval_norm_order <- function(scores_norm) {
  if (!eval_norm_order_(scores_norm)) {print(paste0("Order conservation after normalization: ",
                                                    round(mean(scores_norm %>%group_by(name_score,dataset) %>%
                                                                 mutate(order_before=order(val), order_after=order(normval)) %>%
                                                                 mutate(rank_conserved=order_before==order_after) %>%
                                                                 pull(rank_conserved))*100),"%"))}
  else (print(paste0("Order is conserved after normalization: ",eval_norm_order_(scores_norm))))
}

coerce_pearson <- function(df_median_centerscale_transfo) {
  df_median_centerscale_transfo1 <- df_median_centerscale_transfo %>%
    mutate(typescore=sapply(name_score, function(x) strsplit(x, " ")[[1]][1])) %>%
    filter(typescore!="pearson") %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo2 <- df_median_centerscale_transfo %>%
    filter(name_score %in% c("pearson med_c","pearson med_s","pearson perf_g")) %>%
    group_by(dataset,candidate) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson perf_mean",
           val = NA,
           cat_score = "raw_perf",
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo3 <- df_median_centerscale_transfo %>%
    filter(name_score %in% c("pearson sd_c","pearson sd_s","pearson sd_g")) %>%
    group_by(dataset,candidate) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson sd_mean",
           val = NA,
           cat_score = "stab_perf",
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo = bind_rows(df_median_centerscale_transfo1,
                                            df_median_centerscale_transfo2,
                                            df_median_centerscale_transfo3)
  return(df_median_centerscale_transfo)
}

ranking12 <- function(scores1, scores2, cat_score=NULL, normalization="centerscale") {
  
  if (any(colnames(scores1)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores1 input"))
  }
  if (any(colnames(scores2)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores2 input"))
  }
  if (normalization!="centerscale") {
    stop(print("Normalization should be 'centerscale'"))
  }
  if (!is.list(cat_score)) {
    if (!is.null(cat_score)) {
      stop(print("'cat_score' has to be a list"))
    }
  }
  Norm = function(x) {CenterScaleNorm(x)}
  
  scores1$sim <- as.numeric(scores1$sim)
  scores2$sim <- as.numeric(scores2$sim)
  scores1$values <- as.numeric(scores1$values)
  scores2$values <- as.numeric(scores2$values)
  
  score_df <- bind_rows(scores1,scores2)
  colnames(score_df)[1:5] <- c("score","name_score","simulation","dataset","candidate")
  
  # keep only non redundant scores and non NA scores
  score_to_keep <- c("time","rmse perf_g","mae perf_g","pearson perf_g","pearson med_c","pearson med_s","pearson sd_c","pearson sd_s")
  score_df <- score_df[score_df$name_score %in% score_to_keep,]
  
  # implement the column cat_score
  if (is.null(cat_score)) {
    cat_score <- list("raw_perf"=c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s", "pearson perf_mean"),
                      "stab_perf"=c("mae sd_g", "pearson sd_g", "rmse sd_g", "pearson sd_c", "pearson sd_s", "pearson sd_mean"),
                      "time"=c("time","time sd"))
  }
  score_df$cat_score <- convert_to_cat(cat_score, score_df$name_score)
  
  # compute stability sd_g and time across sims (step 0)
  datasetXcandidate = length(score_df %>% filter(name_score=="mae perf_g") %>%
                               group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std))
  std_df <- data.frame("val"= c(score_df %>% filter(name_score=="mae perf_g") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std),
                                score_df %>% filter(name_score=="pearson perf_g") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std),
                                score_df %>% filter(name_score=="rmse perf_g") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std),
                                score_df %>% filter(name_score=="time") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std)),
                       "name_score"=c(rep('mae sd_g',datasetXcandidate),
                                      rep('pearson sd_g',datasetXcandidate),
                                      rep('rmse sd_g',datasetXcandidate),
                                      rep('time sd',datasetXcandidate)),
                       "dataset"=rep(score_df %>%
                                       filter(name_score=="mae perf_g") %>%
                                       group_by(dataset, candidate) %>% 
                                       summarise(std=sd(score)) %>% pull(dataset),4),
                       "candidate"=rep(score_df %>%
                                         filter(name_score=="mae perf_g") %>%
                                         group_by(dataset, candidate) %>% 
                                         summarise(std=sd(score)) %>% pull(candidate),4))
  rm(datasetXcandidate)
  
  # compute stability sd_c and sd_s across sims (step 0)
  std_df2 <- score_df %>%
    filter(name_score %in% c("pearson sd_c","pearson sd_s")) %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score))
  
  # merge resulting sd dataframes (step 0)
  std_df <- bind_rows(std_df, std_df2)
  
  # compute all medians across sims for raw perf and time categories (step 0)
  score_df_median <- score_df %>%
    filter(cat_score=="raw_perf") %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score))
  time_df_median <- score_df %>%
    filter(cat_score=="time") %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score))
  df_median <- bind_rows(score_df_median, std_df, time_df_median)
  rm(std_df, std_df2, score_df, time_df_median)
  
  # add the column cat_score
  df_median$cat_score <- convert_to_cat(cat_score, df_median$name_score)
  
  # normalize for each dataset and score (step 1)
  df_median_centerscale <- df_median %>%
    group_by(dataset, name_score) %>%
    mutate(normval=Norm(val))
  rm(df_median)
  
  # transform scores such that 1 is the best score (step 2)
  df_median_centerscale_transfo <- df_median_centerscale
  df_median_centerscale_transfo$trendval <- 1 - df_median_centerscale_transfo$normval
  df_median_centerscale_transfo$trendval[grep("pearson med",df_median_centerscale_transfo$name_score)] <- 1 - df_median_centerscale_transfo$trendval[grep("pearson med",df_median_centerscale_transfo$name_score)]
  df_median_centerscale_transfo$trendval[grep("pearson perf",df_median_centerscale_transfo$name_score)] <- 1 - df_median_centerscale_transfo$trendval[grep("pearson perf",df_median_centerscale_transfo$name_score)]
  rm(df_median_centerscale)
  
  return(df_median_centerscale_transfo)
}

ranking_2 <- function(scores1, scores2, cat_score=NULL, normalization="centerscale") {
  
  if (any(colnames(scores1)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores1 input"))
  }
  if (any(colnames(scores2)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores2 input"))
  }
  if (normalization!="centerscale") {
    stop(print("Normalization should be 'centerscale'"))
  }
  if (!is.list(cat_score)) {
    if (!is.null(cat_score)) {
      stop(print("'cat_score' has to be a list"))
    }
  }
  Norm = function(x) {CenterScaleNorm(x)}
  
  scores1$sim <- as.numeric(scores1$sim)
  scores2$sim <- as.numeric(scores2$sim)
  scores1$values <- as.numeric(scores1$values)
  scores2$values <- as.numeric(scores2$values)
  
  score_df <- bind_rows(scores1,scores2)
  colnames(score_df)[1:5] <- c("score","name_score","simulation","dataset","candidate")
  
  # keep only non redundant scores and non NA scores
  score_to_keep <- c("time","rmse perf_g","mae perf_g","pearson perf_g","pearson med_c","pearson med_s","pearson sd_c","pearson sd_s")
  score_df <- score_df[score_df$name_score %in% score_to_keep,]
  
  # implement the column cat_score
  if (is.null(cat_score)) {
    cat_score <- list("raw_perf"=c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s", "pearson perf_mean"),
                      "stab_perf"=c("mae sd_g", "pearson sd_g", "rmse sd_g", "pearson sd_c", "pearson sd_s", "pearson sd_mean"),
                      "time"=c("time","time sd"))
  }
  score_df$cat_score <- convert_to_cat(cat_score, score_df$name_score)
  
  # compute stability sd_g and time across sims (step 0)
  datasetXcandidate = length(score_df %>% filter(name_score=="mae perf_g") %>%
                               group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std))
  std_df <- data.frame("val"= c(score_df %>% filter(name_score=="mae perf_g") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std),
                                score_df %>% filter(name_score=="pearson perf_g") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std),
                                score_df %>% filter(name_score=="rmse perf_g") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std),
                                score_df %>% filter(name_score=="time") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std)),
                       "name_score"=c(rep('mae sd_g',datasetXcandidate),
                                      rep('pearson sd_g',datasetXcandidate),
                                      rep('rmse sd_g',datasetXcandidate),
                                      rep('time sd',datasetXcandidate)),
                       "dataset"=rep(score_df %>%
                                       filter(name_score=="mae perf_g") %>%
                                       group_by(dataset, candidate) %>% 
                                       summarise(std=sd(score)) %>% pull(dataset),4),
                       "candidate"=rep(score_df %>%
                                         filter(name_score=="mae perf_g") %>%
                                         group_by(dataset, candidate) %>% 
                                         summarise(std=sd(score)) %>% pull(candidate),4))
  rm(datasetXcandidate)
  
  # compute stability sd_c and sd_s across sims (step 0)
  std_df2 <- score_df %>%
    filter(name_score %in% c("pearson sd_c","pearson sd_s")) %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score))
  
  # merge resulting sd dataframes (step 0)
  std_df <- bind_rows(std_df, std_df2)
  
  # compute all medians across sims for raw perf and time categories (step 0)
  score_df_median <- score_df %>%
    filter(cat_score=="raw_perf") %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score))
  time_df_median <- score_df %>%
    filter(cat_score=="time") %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score))
  df_median <- bind_rows(score_df_median, std_df, time_df_median)
  rm(std_df, std_df2, score_df, time_df_median)
  
  # add the column cat_score
  df_median$cat_score <- convert_to_cat(cat_score, df_median$name_score)
  
  # normalize for each dataset and score (step 1)
  df_median_centerscale <- df_median %>%
    ungroup() %>%
    mutate(normval=Norm(val))
  rm(df_median)
  
  # transform scores such that 1 is the best score (step 2)
  df_median_centerscale_transfo <- df_median_centerscale
  df_median_centerscale_transfo$trendval <- 1 - df_median_centerscale_transfo$normval
  df_median_centerscale_transfo$trendval[grep("pearson med",df_median_centerscale_transfo$name_score)] <- 1 - df_median_centerscale_transfo$trendval[grep("pearson med",df_median_centerscale_transfo$name_score)]
  df_median_centerscale_transfo$trendval[grep("pearson perf",df_median_centerscale_transfo$name_score)] <- 1 - df_median_centerscale_transfo$trendval[grep("pearson perf",df_median_centerscale_transfo$name_score)]
  rm(df_median_centerscale)
  
  return(df_median_centerscale_transfo)
}

ranking12_notime <- function(scores1, cat_score=NULL, normalization="centerscale") {
  
  if (any(colnames(scores1)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores1 input"))
  }
  if (normalization!="centerscale") {
    stop(print("Normalization should be 'centerscale'"))
  }
  if (!is.list(cat_score)) {
    if (!is.null(cat_score)) {
      stop(print("'cat_score' has to be a list"))
    }
  }
  Norm = function(x) {CenterScaleNorm(x)}
  
  scores1$sim <- as.numeric(scores1$sim)
  scores1$values <- as.numeric(scores1$values)
  
  score_df <- scores1
  colnames(score_df)[1:5] <- c("score","name_score","simulation","dataset","candidate")
  
  # keep only non redundant scores and non NA scores
  score_to_keep <- c("time","rmse perf_g","mae perf_g","pearson perf_g","pearson med_c","pearson med_s","pearson sd_c","pearson sd_s")
  score_df <- score_df[score_df$name_score %in% score_to_keep,]
  score_df <- score_df[!is.na(score_df$score),]
  
  # implement the column cat_score
  if (is.null(cat_score)) {
    cat_score <- list("raw_perf"=c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s", "pearson perf_mean"),
                      "stab_perf"=c("mae sd_g", "pearson sd_g", "rmse sd_g", "pearson sd_c", "pearson sd_s", "pearson sd_mean"),
                      "time"=c("time","time sd"))
  }
  score_df$cat_score <- convert_to_cat(cat_score, score_df$name_score)
  
  # compute stability sd_g and time across sims (step 0)
  datasetXcandidate = length(score_df %>% filter(name_score=="mae perf_g") %>%
                               group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std))
  std_df <- data.frame("val"= c(score_df %>% filter(name_score=="mae perf_g") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std),
                                score_df %>% filter(name_score=="pearson perf_g") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std),
                                score_df %>% filter(name_score=="rmse perf_g") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std)),
                       "name_score"=c(rep('mae sd_g',datasetXcandidate),
                                      rep('pearson sd_g',datasetXcandidate),
                                      rep('rmse sd_g',datasetXcandidate)),
                       "dataset"=rep(score_df %>%
                                       filter(name_score=="mae perf_g") %>%
                                       group_by(dataset, candidate) %>% 
                                       summarise(std=sd(score)) %>% pull(dataset),3),
                       'candidate'=rep(score_df %>%
                                         filter(name_score=="mae perf_g") %>%
                                         group_by(dataset, candidate) %>% 
                                         summarise(std=sd(score)) %>% pull(candidate),3))
  rm(datasetXcandidate)
  
  # compute stability sd_c and sd_s across sims (step 0)
  std_df2 <- score_df %>%
    filter(name_score %in% c("pearson sd_c","pearson sd_s")) %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score))
  
  # merge resulting sd dataframes (step 0)
  std_df <- bind_rows(std_df, std_df2)
  
  # compute all medians across sims for raw perf and time categories (step 0)
  score_df_median <- score_df %>%
    filter(cat_score=="raw_perf") %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score))
  df_median <- bind_rows(score_df_median, std_df)
  rm(std_df, std_df2, score_df)
  
  # add the column cat_score
  df_median$cat_score <- convert_to_cat(cat_score, df_median$name_score)
  
  # normalize for each dataset and score (step 1)
  df_median_centerscale <- df_median %>%
    group_by(dataset, name_score) %>%
    mutate(normval=Norm(val))
  rm(df_median)
  
  # transform scores such that 1 is the best score (step 2)
  df_median_centerscale_transfo <- df_median_centerscale
  df_median_centerscale_transfo$trendval <- 1 - df_median_centerscale_transfo$normval
  df_median_centerscale_transfo$trendval[grep("pearson med",df_median_centerscale_transfo$name_score)] <- 1 - df_median_centerscale_transfo$trendval[grep("pearson med",df_median_centerscale_transfo$name_score)]
  df_median_centerscale_transfo$trendval[grep("pearson perf",df_median_centerscale_transfo$name_score)] <- 1 - df_median_centerscale_transfo$trendval[grep("pearson perf",df_median_centerscale_transfo$name_score)]
  rm(df_median_centerscale)
  
  return(df_median_centerscale_transfo)
}

ranking_2_notime <- function(scores1, cat_score=NULL, normalization="centerscale") {
  
  if (any(colnames(scores1)[1:5]!=c('values','score','sim','dataset','candidate'))) {
    stop(print("Check scores1 input"))
  }
  if (normalization!="centerscale") {
    stop(print("Normalization should be 'centerscale'"))
  }
  if (!is.list(cat_score)) {
    if (!is.null(cat_score)) {
      stop(print("'cat_score' has to be a list"))
    }
  }
  Norm = function(x) {CenterScaleNorm(x)}
  
  scores1$sim <- as.numeric(scores1$sim)
  scores1$values <- as.numeric(scores1$values)
  
  score_df <- scores1
  colnames(score_df)[1:5] <- c("score","name_score","simulation","dataset","candidate")
  
  # keep only non redundant scores and non NA scores
  score_to_keep <- c("time","rmse perf_g","mae perf_g","pearson perf_g","pearson med_c","pearson med_s","pearson sd_c","pearson sd_s")
  score_df <- score_df[score_df$name_score %in% score_to_keep,]
  score_df <- score_df[!is.na(score_df$score),]
  
  # implement the column cat_score
  if (is.null(cat_score)) {
    cat_score <- list("raw_perf"=c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s", "pearson perf_mean"),
                      "stab_perf"=c("mae sd_g", "pearson sd_g", "rmse sd_g", "pearson sd_c", "pearson sd_s", "pearson sd_mean"),
                      "time"=c("time","time sd"))
  }
  score_df$cat_score <- convert_to_cat(cat_score, score_df$name_score)
  
  # compute stability sd_g and time across sims (step 0)
  datasetXcandidate = length(score_df %>% filter(name_score=="mae perf_g") %>%
                               group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std))
  std_df <- data.frame("val"= c(score_df %>% filter(name_score=="mae perf_g") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std),
                                score_df %>% filter(name_score=="pearson perf_g") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std),
                                score_df %>% filter(name_score=="rmse perf_g") %>%
                                  group_by(dataset, candidate) %>% summarise(std=sd(score)) %>% pull(std)),
                       "name_score"=c(rep('mae sd_g',datasetXcandidate),
                                      rep('pearson sd_g',datasetXcandidate),
                                      rep('rmse sd_g',datasetXcandidate)),
                       "dataset"=rep(score_df %>%
                                       filter(name_score=="mae perf_g") %>%
                                       group_by(dataset, candidate) %>% 
                                       summarise(std=sd(score)) %>% pull(dataset),3),
                       'candidate'=rep(score_df %>%
                                         filter(name_score=="mae perf_g") %>%
                                         group_by(dataset, candidate) %>% 
                                         summarise(std=sd(score)) %>% pull(candidate),3))
  rm(datasetXcandidate)
  
  # compute stability sd_c and sd_s across sims (step 0)
  std_df2 <- score_df %>%
    filter(name_score %in% c("pearson sd_c","pearson sd_s")) %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score))
  
  # merge resulting sd dataframes (step 0)
  std_df <- bind_rows(std_df, std_df2)
  
  # compute all medians across sims for raw perf and time categories (step 0)
  score_df_median <- score_df %>%
    filter(cat_score=="raw_perf") %>%
    group_by(name_score, dataset, candidate) %>%
    summarise(val=median(score))
  df_median <- bind_rows(score_df_median, std_df)
  rm(std_df, std_df2, score_df)
  
  # add the column cat_score
  df_median$cat_score <- convert_to_cat(cat_score, df_median$name_score)
  
  # normalize for each dataset and score (step 1)
  df_median_centerscale <- df_median %>%
    ungroup() %>%
    mutate(normval=Norm(val))
  rm(df_median)
  
  # transform scores such that 1 is the best score (step 2)
  df_median_centerscale_transfo <- df_median_centerscale
  df_median_centerscale_transfo$trendval <- 1 - df_median_centerscale_transfo$normval
  df_median_centerscale_transfo$trendval[grep("pearson med",df_median_centerscale_transfo$name_score)] <- 1 - df_median_centerscale_transfo$trendval[grep("pearson med",df_median_centerscale_transfo$name_score)]
  df_median_centerscale_transfo$trendval[grep("pearson perf",df_median_centerscale_transfo$name_score)] <- 1 - df_median_centerscale_transfo$trendval[grep("pearson perf",df_median_centerscale_transfo$name_score)]
  rm(df_median_centerscale)
  
  return(df_median_centerscale_transfo)
}

ranking12345 <- function(scores1, scores2, cat_score=NULL, normalization="centerscale", ranking=ranking12) {
  
  df_median_centerscale_transfo <- ranking(scores1=scores1, scores2=scores2, cat_score=cat_score, normalization=normalization)
  if (!eval_norm_order_(df_median_centerscale_transfo)) {print("Normalization changes order")}
  
  # coerce pearson into a single score
  df_median_centerscale_transfo <- coerce_pearson(df_median_centerscale_transfo)
  if (!eval_norm_order_(df_median_centerscale_transfo)) {print("Pearson coercion changes order")}
  
  # final scoring (step 3-4-5)
  score_final2 <- df_median_centerscale_transfo %>%
    group_by(candidate,cat_score,name_score) %>%
    summarise(score3=mean(trendval, na.rm=T)) %>% #aggregate over datasets (step 3)
    group_by(candidate,cat_score) %>%
    summarise(score4G=geomMean(score3)) %>% #aggregate over scores (step 4)
    group_by(candidate) %>%
    summarise(ScoreFinal=geomMean(score4G)) #aggregate over categories (step 5)
  score_final <- stack(score_final2[,2:ncol(score_final2)])
  score_final$candidate <- score_final2$candidate
  return(score_final)
}

ranking12345_notime <- function(scores1, cat_score=NULL, normalization="centerscale", ranking=ranking12_notime) {
  
  df_median_centerscale_transfo <- ranking(scores1=scores1, cat_score=cat_score, normalization=normalization)
  if (!eval_norm_order_(df_median_centerscale_transfo)) {print("Normalization changes order")}
  
  # coerce pearson into a single score
  df_median_centerscale_transfo <- coerce_pearson(df_median_centerscale_transfo)
  if (!eval_norm_order_(df_median_centerscale_transfo)) {print("Pearson coercion changes order")}
  
  # final scoring (step 3-4-5)
  #score_final1 <- df_median_centerscale_transfo %>%
  #  group_by(candidate, name_score) %>%
  #  summarise(score3=mean(trendval)) %>% #aggregate over datasets (step 3)
  #  group_by(candidate) %>%
  #  summarise(ScoreFinalA=mean(score3),
  #            ScoreFinalG=geomMean(score3)) #aggregate over the rest (step 4-5)
  score_final2 <- df_median_centerscale_transfo %>%
    group_by(candidate,cat_score,name_score) %>%
    summarise(score3=mean(trendval, na.rm=T)) %>% #aggregate over datasets (step 3)
    group_by(candidate,cat_score) %>%
    summarise(score4G=geomMean(score3)) %>% #aggregate over scores (step 4)
    group_by(candidate) %>%
    summarise(ScoreFinal=geomMean(score4G)) #aggregate over categories (step 5)
  score_final <- stack(score_final2[,2:ncol(score_final2)])
  score_final$candidate <- score_final2$candidate
  return(score_final)
}

ranking12453 <- function(scores1, scores2, cat_score=NULL, normalization="centerscale", ranking=ranking12) {
  
  df_median_centerscale_transfo <- ranking(scores1=scores1, scores2=scores2, cat_score=cat_score, normalization=normalization)
  if (!eval_norm_order_(df_median_centerscale_transfo)) {print("Normalization changes order")}
  
  # coerce pearson into a single score
  df_median_centerscale_transfo <- coerce_pearson(df_median_centerscale_transfo)
  if (!eval_norm_order_(df_median_centerscale_transfo)) {print("Pearson coercion changes order")}
  
  # final scoring (step 4-5-3)
  score_final2 <- df_median_centerscale_transfo %>%
    group_by(candidate,dataset,cat_score) %>%
    summarise(score4G=geomMean(trendval)) %>% #aggregate over scores (step 4)
    group_by(candidate,dataset) %>%
    summarise(score5GG=geomMean(score4G)) %>% #aggregate over categories (step 5)
    group_by(candidate) %>%
    summarise(ScoreFinal=mean(score5GG, na.rm=T)) #aggregate over datasets (step 3)
  score_final <- stack(score_final2[,2:ncol(score_final2)])
  score_final$candidate <- score_final2$candidate
  return(score_final)
}

ranking12453_notime <- function(scores1, cat_score=NULL, normalization="centerscale", ranking=ranking12_notime) {
  
  df_median_centerscale_transfo <- ranking(scores1=scores1, cat_score=cat_score, normalization=normalization)
  if (!eval_norm_order_(df_median_centerscale_transfo)) {print("Normalization changes order")}
  
  # coerce pearson into a single score
  df_median_centerscale_transfo <- coerce_pearson(df_median_centerscale_transfo)
  if (!eval_norm_order_(df_median_centerscale_transfo)) {print("Pearson coercion changes order")}
  
  # final scoring (step 4-5-3)
  score_final2 <- df_median_centerscale_transfo %>%
    group_by(candidate,dataset,cat_score) %>%
    summarise(score4G=geomMean(trendval)) %>% #aggregate over scores (step 4)
    group_by(candidate,dataset) %>%
    summarise(score5GG=geomMean(score4G)) %>% #aggregate over categories (step 5)
    group_by(candidate) %>%
    summarise(ScoreFinal=mean(score5GG, na.rm=T)) #aggregate over datasets (step 3)
  score_final <- stack(score_final2[,2:ncol(score_final2)])
  score_final$candidate <- score_final2$candidate
  return(score_final)
}

ranking_avgrank <- function(scores1, scores2, cat_score=NULL, normalization="centerscale", ranking=ranking12) {
  
  df_median_centerscale_transfo <- ranking(scores1=scores1, scores2=scores2, cat_score=cat_score, normalization=normalization)
  if (!eval_norm_order_(df_median_centerscale_transfo)) {print("Normalization changes order")}
  
  # coerce pearson into a single score
  df_median_centerscale_transfo <- coerce_pearson(df_median_centerscale_transfo)
  if (!eval_norm_order_(df_median_centerscale_transfo)) {print("Pearson coercion changes order")}
  
  # Get the rank for each judge (dataset x name_score) and the avg rank
  ranks <- df_median_centerscale_transfo %>%
    group_by(name_score,dataset) %>%
    mutate(rank=rank(trendval, ties.method = "average")) %>%
    select(name_score, dataset, rank, candidate)
  avg_rank <- ranks %>%
    group_by(candidate) %>%
    summarise(avg_rank = mean(rank, na.rm=T))
  score_final <- stack(avg_rank[,2:ncol(avg_rank)])
  score_final$candidate <- avg_rank$candidate
  
  return(score_final)
}

ranking_topsis <- function(scores1, scores2, cat_score=NULL, normalization="centerscale", ranking=ranking12) {
  
  df_median_centerscale_transfo <- ranking(scores1=scores1, scores2=scores2, cat_score=cat_score, normalization=normalization)
  if (!eval_norm_order_(df_median_centerscale_transfo)) {print("Normalization changes order")}
  
  # coerce pearson into a single score
  df_median_centerscale_transfo <- coerce_pearson(df_median_centerscale_transfo)
  if (!eval_norm_order_(df_median_centerscale_transfo)) {print("Pearson coercion changes order")}
  
  # Compute best and worst candidate
  df_median_centerscale_transfo = df_median_centerscale_transfo %>% group_by(name_score,dataset) %>%
    mutate(archetype_best=max(trendval, na.rm=T),
           archetype_worst=min(trendval, na.rm=T))

  # Compute distances and topsis score
  topsis_score <- df_median_centerscale_transfo %>% group_by(candidate,cat_score) %>%
    summarize(d_best=euclidean_vector_dist(trendval,archetype_best),
              d_worst=euclidean_vector_dist(trendval,archetype_worst)) %>%
    group_by(candidate) %>% summarise(d_best=geomMean(d_best),
                                      d_worst=geomMean(d_worst)) %>% group_by(candidate) %>%
    summarise(values=d_worst/(d_worst+d_best))
  return(topsis_score)
}

ranking_consensus <- function(scores1, scores2, cat_score=NULL, normalization="centerscale", ranking=ranking12) {
  rank1 <- ranking12453(scores1, scores2, cat_score=cat_score, normalization=normalization, ranking=ranking)
  rank2 <- ranking_topsis(scores1, scores2, cat_score=cat_score, normalization=normalization, ranking=ranking)
  rank <- bind_rows(rank1,rank2) %>% group_by(candidate) %>% summarise(values=mean(values))
  return(rank)
}

#######
####### PLOT FUNCTION ####### 
#######

plotting <- function(final_rank) {
  p <- ggplot(final_rank, aes(x=ind, y=values, color=candidate)) +
    geom_point() +
    geom_line(aes(group=candidate)) +
    scale_color_social_d() +
    theme_modern(axis.text.angle = 45)
  return(p)
}

#######
####### RANKING FUNCTION FOR CONCAT ####### 
#######
ranking_concat_per_method <- function(score_per_method) {
  # keep only non redundant scores
  score_to_keep <- c("rmse perf_g","mae perf_g","pearson perf_g","pearson med_c","pearson med_s","pearson sd_c","pearson sd_s")
  score_per_method <- score_per_method[score_per_method$name_score %in% score_to_keep,]
  
  # implement the column cat_score
  cat_score <- list("raw_perf"=c("rmse perf_g", "mae perf_g", "pearson perf_g", "pearson med_c", "pearson med_s", "pearson perf_mean"),
                    "stab_perf"=c("mae sd_g", "pearson sd_g", "rmse sd_g", "pearson sd_c", "pearson sd_s", "pearson sd_mean"))
  score_per_method$cat_score <- convert_to_cat(cat_score, score_per_method$name_score)
  
  # compute stability sd_g across sims (step 0)
  datasetXdatatype = length(unique(score_per_method$dataset))*4
  std_df <- data.frame("val"= c(score_per_method %>% filter(name_score=="mae perf_g") %>%
                                  group_by(dataset, datatype) %>% summarise(std=sd(score)) %>% pull(std),
                                score_per_method %>% filter(name_score=="pearson perf_g") %>%
                                  group_by(dataset, datatype) %>% summarise(std=sd(score)) %>% pull(std),
                                score_per_method %>% filter(name_score=="rmse perf_g") %>%
                                  group_by(dataset, datatype) %>% summarise(std=sd(score)) %>% pull(std)),
                       "name_score"=c(rep('mae sd_g',datasetXdatatype),
                                      rep('pearson sd_g',datasetXdatatype),
                                      rep('rmse sd_g',datasetXdatatype)),
                       "dataset"=rep(score_per_method %>%
                                       filter(name_score=="mae perf_g") %>%
                                       group_by(dataset, datatype) %>% 
                                       summarise(std=sd(score)) %>% pull(dataset),3),
                       "datatype"=rep(score_per_method %>%
                                        filter(name_score=="mae perf_g") %>%
                                        group_by(dataset, datatype) %>% 
                                        summarise(std=sd(score)) %>% pull(datatype),3))
  rm(datasetXdatatype)
  
  # compute stability sd_c and sd_s across sims (step 0)
  std_df2 <- score_per_method %>%
    filter(cat_score=="stab_perf") %>%
    group_by(name_score, dataset, datatype) %>%
    summarise(val=median(score))
  
  # merge resulting sd dataframes (step 0)
  std_df <- bind_rows(std_df, std_df2)
  
  # compute all medians across sims for raw_perf category (step 0)
  df_median <- score_per_method %>%
    filter(cat_score=="raw_perf") %>%
    group_by(name_score, dataset, datatype) %>%
    summarise(val=median(score))
  df_median <- bind_rows(df_median, std_df)
  rm(std_df, std_df2)
  
  # add the column cat_score
  df_median$cat_score <- convert_to_cat(cat_score, df_median$name_score)
  
  # normalize for each dataset and score (step 1)
  df_median_centerscale <- df_median %>%
    group_by(dataset, name_score) %>%
    mutate(normval=CenterScaleNorm(val))
  rm(df_median)
  
  # transform scores such that 1 is the best score (step 2)
  df_median_centerscale_transfo <- df_median_centerscale
  df_median_centerscale_transfo$trendval <- 1 - df_median_centerscale_transfo$normval
  df_median_centerscale_transfo$trendval[grep("pearson m",df_median_centerscale_transfo$name_score)] <- 1 - df_median_centerscale_transfo$trendval[grep("pearson m",df_median_centerscale_transfo$name_score)]
  df_median_centerscale_transfo$trendval[grep("pearson p",df_median_centerscale_transfo$name_score)] <- 1 - df_median_centerscale_transfo$trendval[grep("pearson p",df_median_centerscale_transfo$name_score)]
  rm(df_median_centerscale)
  
  # coerce pearson into a single score
  df_median_centerscale_transfo1 <- df_median_centerscale_transfo %>%
    mutate(typescore=sapply(name_score, function(x) strsplit(x, " ")[[1]][1])) %>%
    filter(typescore!="pearson") %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo2 <- df_median_centerscale_transfo %>%
    filter(name_score %in% c("pearson med_c","pearson med_s","pearson perf_g")) %>%
    group_by(dataset,candidate) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson perf_mean",
           val = NA,
           cat_score = "raw_perf",
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo3 <- df_median_centerscale_transfo %>%
    filter(name_score %in% c("pearson sd_c","pearson sd_s","pearson sd_g")) %>%
    group_by(dataset,candidate) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson sd_mean",
           val = NA,
           cat_score = "stab_perf",
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo = bind_rows(df_median_centerscale_transfo1,
                                            df_median_centerscale_transfo2,
                                            df_median_centerscale_transfo3)
  
  score_final2 <- df_median_centerscale_transfo %>%
    group_by(datatype,cat_score,name_score) %>%
    summarise(score3=mean(trendval, na.rm=T)) %>% #aggregate over datasets (step 3)
    group_by(datatype,cat_score) %>%
    summarise(score4G=geomMean(score3)) %>% #aggregate over scores (step 4)
    group_by(datatype) %>%
    summarise(ScoreFinal=geomMean(score4G)) #aggregate over categories (step 5)
  score_final <- stack(score_final2[,2:ncol(score_final2)])
  score_final$datatype <- score_final2$datatype
  score_final$datatype <- factor(score_final$datatype, levels=c('rna','met','CimpleG','TOAST'))
  return(score_final)
}

#####
##### OLD FUNCS
#####
ranking1234 <- function(scores1, scores2, cat_score=NULL, normalization="centerscale") {
  
  df_median_centerscale_transfo <- ranking12(scores1=scores1, scores2=scores2, cat_score=cat_score, normalization=normalization)
  
  # coerce pearson into a single score
  df_median_centerscale_transfo1 <- df_median_centerscale_transfo %>%
    mutate(typescore=sapply(name_score, function(x) strsplit(x, " ")[[1]][1])) %>%
    filter(typescore!="pearson") %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo2 <- df_median_centerscale_transfo %>%
    filter(name_score %in% c("pearson med_c","pearson med_s","pearson perf_g")) %>%
    group_by(dataset,candidate) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson perf_mean",
           val = NA,
           cat_score = "raw_perf",
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo3 <- df_median_centerscale_transfo %>%
    filter(name_score %in% c("pearson sd_c","pearson sd_s","pearson sd_g")) %>%
    group_by(dataset,candidate) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson sd_mean",
           val = NA,
           cat_score = "stab_perf",
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo = bind_rows(df_median_centerscale_transfo1,
                                            df_median_centerscale_transfo2,
                                            df_median_centerscale_transfo3)
  
  # final scoring (step 3-4)
  score_final2 <- df_median_centerscale_transfo %>%
    group_by(candidate,cat_score,name_score) %>%
    summarise(score3=mean(trendval, na.rm=T)) %>% #aggregate over datasets (step 3)
    group_by(candidate,cat_score) %>%
    summarise(score4G=geomMean(score3)) #aggregate over scores (step 4)
  score_final <- data.frame(values=score_final2$score4G)
  score_final$ind <- score_final2$cat_score
  score_final$candidate <- score_final2$candidate
  return(score_final)
}
ranking1245 <- function(scores1, scores2, cat_score=NULL, normalization="centerscale") {
  
  df_median_centerscale_transfo <- ranking12(scores1=scores1, scores2=scores2, cat_score=cat_score, normalization=normalization)
  
  # coerce pearson into a single score
  df_median_centerscale_transfo1 <- df_median_centerscale_transfo %>%
    mutate(typescore=sapply(name_score, function(x) strsplit(x, " ")[[1]][1])) %>%
    filter(typescore!="pearson") %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo2 <- df_median_centerscale_transfo %>%
    filter(name_score %in% c("pearson med_c","pearson med_s","pearson perf_g")) %>%
    group_by(dataset,candidate) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson perf_mean",
           val = NA,
           cat_score = "raw_perf",
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo3 <- df_median_centerscale_transfo %>%
    filter(name_score %in% c("pearson sd_c","pearson sd_s","pearson sd_g")) %>%
    group_by(dataset,candidate) %>%
    summarise(trendval=mean(trendval, na.rm = T)) %>%
    mutate(name_score = "pearson sd_mean",
           val = NA,
           cat_score = "stab_perf",
           normval = NA) %>%
    select(colnames(df_median_centerscale_transfo))
  df_median_centerscale_transfo = bind_rows(df_median_centerscale_transfo1,
                                            df_median_centerscale_transfo2,
                                            df_median_centerscale_transfo3)
  
  # final scoring (step 4-5)
  score_final2 <- df_median_centerscale_transfo %>%
    group_by(candidate,dataset,cat_score) %>%
    summarise(score4G=geomMean(trendval)) %>% #aggregate over scores (step 4)
    group_by(candidate,dataset) %>%
    summarise(score5GG=geomMean(score4G)) #aggregate over categories (step 5)
  score_final <- data.frame(values=score_final2$score5GG)
  score_final$ind <- score_final2$dataset
  score_final$candidate <- score_final2$candidate
  return(score_final)
}
