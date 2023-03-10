generate_proportions <- function(n_samples, celltypes, alph, varCrit, num.cell.types) {
  # Compute proportions of various compartments
  comp_mix <-t(gtools::rdirichlet(n = n_samples, alpha=alph*varCrit))
  # Compute relative proportions of tumor infiltrates
  tumor_mix <- cbind(t(gtools::rdirichlet(n = 1*n_samples / 6, alpha = c(0.5, 0.5))),
                     t(gtools::rdirichlet(n = 2*n_samples / 6, alpha = c(0.25, 0.75))),
                     t(gtools::rdirichlet(n = 1*n_samples / 6, alpha = c(0.75, 0.25))),
                     t(gtools::rdirichlet(n = 2*n_samples / 6, alpha = c(0.1, 0.9))))
  # Merge total matrix proportion
  tot_mix <- matrix(nrow = num.cell.types, ncol = n_samples)
  rownames(tot_mix) = c(celltypes[-grep("Cancer",celltypes)],
                        celltypes[grep("Cancer",celltypes)])
  colnames(tot_mix) = paste("train", 1:n_samples, sep = "")
  tot_mix[1:(num.cell.types - 2), ] <- comp_mix[1:(num.cell.types - 2),]
  tot_mix[num.cell.types - 1, ] <- tumor_mix[1, ] * comp_mix[num.cell.types - 1,]
  tot_mix[num.cell.types,] <- tumor_mix[2, ] * comp_mix[num.cell.types - 1,]
  tot_mix <- tot_mix[celltypes,]
  #tot_mix <- apply(tot_mix, 2, function(x) x/sum(x)) # step not done in Slim's
  return(tot_mix)
}

generate_simu_set <- function(data, celltypes, n_samples=500, alph=c(0.01, 0.04, 0.03, 0.15, 0.46, 0.01, 0.01, 0.29), varCrit=10) {
  num.cell.types <- length(celltypes)
  profile_rna = data$rna[,((ncol(data$rna)-num.cell.types+1):ncol(data$rna))]
  profile_met = data$met[,((ncol(data$met)-num.cell.types+1):ncol(data$met))]
  colnames(profile_rna) = ref_row_order
  colnames(profile_met) = ref_row_order
  test_solution <- generate_proportions(n_samples, celltypes, alph, varCrit, num.cell.types) #prop x sample
  test_data_rna <- round(as.matrix(profile_rna[,celltypes]) %*% test_solution[celltypes,]) #raw gene x sample
  #M <- matrix(0L, nrow = sum(rowSums(test_data_rna)==0), ncol = ncol(test_data_rna))
  #M[cbind(sequence(sum(rowSums(test_data_rna)==0)),
  #        sample(ncol(test_data_rna), sum(rowSums(test_data_rna)==0), TRUE))] <- 1L
  #test_data_rna[rowSums(test_data_rna)==0,] <- M
  test_data_met <- as.matrix(profile_met[,celltypes]) %*% test_solution[celltypes,] #probe x sample
  return(list(test_data_rna, test_data_met, test_solution[celltypes,], profile_rna[,celltypes]))
}

generate_simu_noise <- function(test_data_rna, test_data_met, p, sd_rna=1, sd_met=3) {
  add_noise_gaussian = function(dt, sd, mean=0) {
    noise = matrix(rnorm(prod(dim(dt)), mean = mean, sd = sd), nrow = nrow(dt))
    data_noise = dt + noise
    data_noise[data_noise<0] <- dt[data_noise<0]
    return(data_noise)
  }
  add_noise_nb = function(dt, p, sd) {
    delta = matrix(rnorm(prod(dim(dt)), mean = 0, sd = sd),nrow=nrow(dt))
    mu_i_0 = dt/colSums(dt) * mean(colSums(dt))
    sigma_i = (1.8*p + 1/sqrt(mu_i_0))*exp(delta/2)
    shape = 1/(sigma_i^2)
    scale = mu_i_0/(shape + .Machine$double.eps)
    mu_i = matrix(rgamma(prod(dim(dt)), shape=shape, scale=scale),
                  nrow=nrow(dt))
    v_i = matrix(rpois(prod(dim(dt)), mu_i),
                 nrow=nrow(dt))
    return(v_i)
  }
  row_names=rownames(test_data_rna)
  test_data_rna <- add_noise_nb(test_data_rna, p, sd_rna)
  rownames(test_data_rna)=row_names
  beta_val <- test_data_met
  m_val <- add_noise_gaussian(log2(beta_val/(1-beta_val)), sd_met)
  test_data_met <- 2^m_val/(2^m_val+1) #probe x sample
  test_data_rna <- round(test_data_rna)
  colnames(test_data_rna)=colnames(test_data_met)
  return(list(test_data_rna, test_data_met))
}

scale_data <- function(test_data_rna, test_data_met=NULL) {
  size.factor <- DESeq2::estimateSizeFactorsForMatrix(test_data_rna) 
  test_data_rna <- sweep(test_data_rna, 2, size.factor, "/")
  test_data_rna <- log2(test_data_rna + 1)
  test_data_rna <- test_data_rna/sqrt(nrow(test_data_rna))
  if (!is.null(test_data_met)) {
    test_data_met <- t(scale(t(test_data_met)))
    test_data_met <- test_data_met/sqrt(nrow(test_data_met))
    return(list(test_data_rna, test_data_met))
  }
  else return(test_data_rna)
}
