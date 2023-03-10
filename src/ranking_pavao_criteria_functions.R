# transfo score matrix (for old Alona scripts)
reformat_score_mat <- function(df) {
  scores = df %>%
    group_by(name_score, deconv_method) %>%
    summarise(score=mean(trendval)) %>%
    pull(score)
  row_name = unique(df %>%
                      group_by(name_score, deconv_method) %>%
                      summarise(score=mean(trendval)) %>%
                      pull(deconv_method))
  col_name = unique(df %>%
                      group_by(name_score, deconv_method) %>%
                      summarise(score=mean(trendval)) %>%
                      pull(name_score))
  mat = matrix(scores, byrow=F, nrow=length(row_name))
  colnames(mat) = col_name
  rownames(mat) = row_name
  return(mat)
}

# reorder score matrix to get a candidate x judge matrix instead of a data.frame
make_pavao_matrix <- function(df) {
  candidate = unique(df %>% pull(candidate))
  judge = unique(df %>% pull(judge))
  mat = matrix(nrow=length(candidate), ncol=length(judge))
  colnames(mat) = judge
  rownames(mat) = candidate
  for (i in seq(nrow(df))) {
    mat[df$candidate[i],df$judge[i]] = df$val[i]
  }
  return(mat)
}

# bootstrap judges or candidates
bootsample <- function(matrix, n, axis=c("judge","candidate"), idx=F, replace=T) {
  axis <- match.arg(axis)
  if (axis=='judge') {
    if (n>=ncol(matrix)) {stop(print("n should be less than the number of judges"))}
    choice <- sample(seq(ncol(matrix)), n, replace=replace)
    if (idx) {return(choice)}
    else {return(matrix[,choice])}
  }
  else {
    if (n>=nrow(matrix)) {stop(print("n should be less than the number of candidates"))}
    choice <- sample(seq(nrow(matrix)), n, replace=replace)
    if (idx) {return(choice)}
    else {return(matrix[choice,])}
  }
}

# theoretical criteria
maj_criterion <- function(matrix,winner) {
  get_best_candidate_per_judge <- apply(matrix, 2, function(x) names(which.max(x)))
  most_voted <- sort(table(get_best_candidate_per_judge), decreasing = TRUE)[1]
  get_best_candidate <- unname(ifelse(most_voted>=ncol(matrix)/2,names(most_voted),NA))
  if (!is.na(get_best_candidate)) {
    ifelse(get_best_candidate==winner,T,F)
  }
  else {return(NA)}
}

condorcet_criterion <- function(matrix,winner) {
  winner_idx = which(rownames(matrix)==winner)
  condorcet_win = sapply(seq(nrow(matrix)), function(x) {
      maj_winner = c()
      sub_matrix <- rbind(matrix[winner_idx,],matrix[x,])
      rownames(sub_matrix)=c(rownames(matrix)[winner_idx],rownames(matrix)[x])
      if (!all(matrix[winner_idx,]==matrix[x,])) {
        get_best_candidate_per_judge <- apply(sub_matrix, 2, function(z) names(which.max(z)))
        most_voted <- sort(table(get_best_candidate_per_judge), decreasing = TRUE)[1]
        maj_winner <- c(maj_winner,unname(ifelse(most_voted>=ncol(sub_matrix)/2,names(most_voted),"no winner")))
      }
      else {maj_winner <- c(maj_winner,NA)}
      return(maj_winner)
    })
  condorcet_win <- condorcet_win[!is.na(condorcet_win)]
  return(ifelse(all(condorcet_win==winner),T,F))
}

consistency_criterion <- function(matrix,winner,n_parts=3,n_bootstrap=10) {
  n_judge=ncol(matrix)
  winner_bootstrap = c()
  for (k in seq(n_bootstrap)) {
    parts = split(sample(seq(n_judge)), sort(seq(n_judge)%%n_parts))
    winner_part = c()
    # determine parts
    for (i in parts) {
      sub_matrix <- matrix[,i]
      if (length(i)==1) {
        get_best_candidate_per_judge <- names(which.max(sub_matrix))
        n_judge=1
      }
      else {get_best_candidate_per_judge <- apply(sub_matrix, 2, function(x) names(which.max(x)))}
      most_voted <- sort(table(get_best_candidate_per_judge), decreasing = TRUE)[1]
      get_best_method <- unname(ifelse(most_voted>=n_judge/2,names(most_voted),NA))
      winner_part <- c(winner_part,get_best_method)
      }
    if (length(unique(winner_part)) == 1) {
      winner_bootstrap = c(winner_bootstrap, unique(winner_part))
    }
    else {winner_bootstrap = c(winner_bootstrap, NA)}
  }
  winner_bootstrap <- winner_bootstrap[!is.na(winner_bootstrap)]
  if (length(unique(winner_bootstrap)) == 1) {
    score = ifelse(unique(winner_bootstrap)==winner,T,F)
  }
  else {score = F}
  return(score)
}

participation_criterion <- function(matrix,winner,f_M,n_bootstrap=10) {
  res = c()
  for (i in seq(n_bootstrap)) {
    n_judge_remove = ceiling(ncol(matrix)*.1)
    # find judges that do not favor the best candidate
    winner_idx = which(rownames(matrix)==winner)
    R <- apply(matrix,2,function(x) rank(1-x))
    if (is.null(dim(R[-winner_idx,]))) {potential_judge <- setdiff(seq(nrow(matrix)),winner_idx)}
    else {potential_judge <- colnames(matrix)[apply(R[-winner_idx,],2,function(x) any(x==1))]}
    if (length(potential_judge)==0) {res=c(res,NA)}
    else if (length(potential_judge)>n_judge_remove) {
      potential_candidate <- apply(R[-winner_idx,potential_judge],2,function(x) names(x)[x==1])
      judge_remove=sample(potential_judge, n_judge_remove)
      candidate_judge_remove=unname(potential_candidate[judge_remove])
      sub_matrix <- matrix[,!(colnames(matrix) %in% judge_remove)]
      sub_R <- rank(1-f_M(sub_matrix))
      if (sub_R[winner_idx]==1) {res=c(res,T)}
      else {res=c(res,F)}
      }
    else {
      if (length(potential_judge)==1) {potential_candidate <- names(R[-winner_idx,potential_judge])[R[-winner_idx,potential_judge]==1]}
      else {potential_candidate <- apply(R[-winner_idx,potential_judge],2,function(x) names(x)[x==1])}
      judge_remove=potential_judge;candidate_judge_remove=unname(potential_candidate)
      sub_matrix <- matrix[,!(colnames(matrix) %in% judge_remove)]
      sub_R <- rank(1-f_M(sub_matrix))
      if (sub_R[winner_idx]==1) {res=c(res,T)}
      else {res=c(res,F)}
    }
    # remove those judges and apply f_M to check the criterion
  }
  return(all(res))
} 

iia_criterion <- function(matrix,winner,f_M) {
  res = c()
  winner_idx = which(rownames(matrix)==winner)
  for (i in seq(nrow(matrix))) {
    if (i!=winner_idx) {
      # remove all candidates except one and the winner and apply f_M to check the criterion
      sub_matrix <- matrix[c(i,winner_idx),]
      sub_R <- rank(1-f_M(sub_matrix))
      R <- rank(1-f_M(matrix))
      if (all(order(sub_R)==order(R[c(i,winner_idx)]))) {res=c(res,T)}
      else {res=c(res,F)}
    }
  }
  return(all(res))
} 

# empirical criteria
avg_rank <- function(matrix,winner) {
  winner_idx = which(rownames(matrix)==winner)
  rank_matrix <- apply(matrix,2,function(x) rank(1-x))
  avgrank <- mean(rank_matrix[winner_idx,])
  return(1-(avgrank-1)/(ncol(matrix)-1))
}

condorcet_rate <- function(matrix,winner) {
  winner_idx = which(rownames(matrix)==winner)
  condorcet_winner = sapply(seq(nrow(matrix)), function(x) {
    maj_winner = c()
    sub_matrix <- rbind(matrix[winner_idx,],matrix[x,])
    rownames(sub_matrix)=rownames(matrix)[c(winner_idx,x)]
    if (!all(matrix[winner_idx,]==matrix[x,])) {
      get_best_method_per_score <- apply(sub_matrix, 2, function(z) names(which.max(z)))
      most_voted <- sort(table(get_best_method_per_score), decreasing = TRUE)[1]
      maj_winner <- c(maj_winner,unname(ifelse(most_voted>=ncol(matrix)/2,names(most_voted),"no winner")))
    }
    else {maj_winner <- c(maj_winner,NA)}
    maj_winner
  })
  condorcet_winner <- condorcet_winner[!is.na(condorcet_winner)]
  return(mean(condorcet_winner==winner))
}

generalization_criterion <- function(matrix,f_M,n_bootstrap=10,percent_train=0.75) {
  n_judges_train = floor(ncol(matrix)*percent_train)
  generalization_res = c()
  for (i in seq(n_bootstrap)) {
    judges_train = sort(bootsample(matrix, n_judges_train, axis="judge", idx=T, replace=F))
    judges_test = setdiff(seq(ncol(matrix)), judges_train)
    matrix_train <- matrix[,judges_train]
    matrix_test <- matrix[,judges_test]
    rank_train <- rank(1-f_M(matrix_train))
    if (!is.null(dim(matrix_test))) {
      rank_test <- apply(matrix_test,2,function(x) rank(1-x))
      res = unname(apply(rank_test, 2, function(x)
        cor(rank_train, x, method = "spearman")))
    }
    else {
      rank_test <- rank(1-matrix_test)
      res = cor(rank_train, rank_test, method = "spearman")
    }
    generalization_res = c(generalization_res,
                           sum(res)/ncol(matrix))
  }
  return(mean(generalization_res))
}


