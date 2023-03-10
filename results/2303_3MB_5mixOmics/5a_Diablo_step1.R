library(mixOmics)

## ----
## Set parameters
## ----
lot=1 # vary between 1, 2, 3 (otherwise memory problems)
sim=1 # vary between 1 and 10

## ----
## load D MET&RNA
## ----
input_path_met = sort(list.files("../2210_0simu/simulations/met/", pattern = "sim"))
input_path_rna = sort(list.files("../2210_0simu/simulations/rna/", pattern = "sim"))
input_path_nlot = sort(list.files("../2210_0simu/simulations/rna/", pattern = "T_rna_ref"))
n_lot = length(input_path_nlot)
n_sims = length(input_path_rna)/n_lot
rm(input_path_nlot)

Tmet = lapply(sort(list.files("../2210_0simu/simulations/met/", pattern = "T_met_ref")), function(x) readRDS(paste0("../2210_0simu/simulations/met/",x)))[[lot]]
Dmet = pbapply::pblapply(input_path_met, function(x) readRDS(paste0("../2210_0simu/simulations/met/",x))$D_met_sim)
Dmet = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) Dmet[[10*(lot-1)+sim]])
})[[lot]]

Trna = lapply(sort(list.files("../2210_0simu/simulations/rna/", pattern = "T_rna_ref")), function(x) readRDS(paste0("../2210_0simu/simulations/rna/",x)))[[lot]]
Drna = pbapply::pblapply(input_path_rna, function(x) readRDS(paste0("../2210_0simu/simulations/rna/",x))$D_rna_sim)
Drna = lapply(seq(n_lot), function(lot) {
  lapply(seq(n_sims), function(sim) Drna[[10*(lot-1)+sim]])
})[[lot]]
rm(input_path_met,input_path_rna)

## ----
## mixOmics preprocessing
## ----
# Store as list
X <- list(methylation = Dmet, 
          mRNA = Drna)
Y <- list(methylation = Tmet, 
          mRNA = Trna)
rm(Dmet,Tmet,Drna,Trna)
# C design matrix
Cdesign <- matrix(0.5, ncol = length(X)*2, nrow = length(X)*2, 
                 dimnames = list(c(paste("X",names(X)),paste("Y",names(Y))),
                                 c(paste("X",names(X)),paste("Y",names(Y)))))
diag(Cdesign) <- 0
# Observe correlations
res1.pls <- pls(t(X$methylation[[sim]]), t(X$mRNA[[sim]]), ncomp = 1)
cor(res1.pls$variates$X, res1.pls$variates$Y)
res2.pls <- pls(t(Y$methylation), t(Y$mRNA), ncomp = 1)
cor(res2.pls$variates$X, res2.pls$variates$Y)
res3.pca <- pca(t(cbind(X$methylation[[sim]],Y$methylation)), ncomp = 1)
plot(res3.pca$variates$X[1:ncol(X$methylation[[sim]])], col="black", ylim=c(-15,15))
points(res3.pca$variates$X[(1+ncol(X$methylation[[sim]])):length(res3.pca$variates$X)], col="red")
res4.pca <- pca(t(cbind(X$mRNA[[sim]],Y$mRNA)), ncomp = 1)
plot(res4.pca$variates$X[1:ncol(X$methylation[[sim]])], col="black", ylim=c(-150000,180000))
points(res4.pca$variates$X[(1+ncol(X$methylation[[sim]])):length(res4.pca$variates$X)], col="red")

## ----
## PLS2 - regression mode
## ----
X1=t(X$methylation[[sim]])
Y1=t(X$mRNA[[sim]])
X1 <- X1[,!apply(X1,2,function(x) sd(x)==0)]
X1 <- X1[,-nearZeroVar(X1)$Position]
Y1 <- Y1[,!apply(Y1,2,function(x) sd(x)==0)]
tune.pls2 <- pls(X = X1, Y = Y1, ncomp = 2, mode = 'regression', progressBar = T)
Q2.pls <- perf(tune.pls2, validation = 'Mfold', 
               folds = 6, nrepeat = 2, progressBar = T)
plot(Q2.pls, criterion = 'Q2.total')
list.keepX <- c(seq(50, 100, 10))
list.keepY <- c(seq(50, 100, 10))
tune.spls2 <- tune.spls(X = X1, Y = Y1,
                       test.keepX = list.keepX, test.keepY = list.keepY,
                       ncomp = 2, 
                       nrepeat = 1, folds = 6, mode = 'regression', progressBar = T)
biplot(tune.spls2)
choice.keepX <- tune.spls2$choice.keepX
choice.keepY <- tune.spls2$choice.keepY
choice.ncomp <- length(choice.keepX)
spls2 <- spls(X, Y,
              ncomp = choice.ncomp, keepX = choice.keepX, keepY = choice.keepY,
              mode = "regression")
vip.spls2 <- vip(spls2)
head(vip.spls2[selectVar(spls2, comp = 1)$X$name,1])
perf.spls2 <- perf(spls2, validation = 'Mfold', folds = 10, nrepeat = 5)
plotIndiv(spls2, ind.names = FALSE, 
          group = , 
          pch = , 
          col.per.group = color.mixo(1:4),
          legend = TRUE, legend.title = 'Time', 
          legend.title.pch = 'Dose')
plotArrow(spls2, ind.names = FALSE, 
          group = ,
          col.per.group = color.mixo(1:4),
          legend.title = 'Time.Group')
plotVar(spls2, cex = c(3,4), var.names = c(FALSE, TRUE))
plot(c(1,2), perf.pls$measures$cor.upred$summary$mean, 
     col = 'blue', pch = 16, 
     ylim = c(0.6,1), xaxt = 'n',
     xlab = 'Component', ylab = 't or u Cor', 
     main = 's/PLS performance based on Correlation')
axis(1, 1:2)  # X-axis label
points(perf.pls$measures$cor.tpred$summary$mean, col = 'red', pch = 16)
points(perf.spls$measures$cor.upred$summary$mean, col = 'blue', pch = 17)
points(perf.spls$measures$cor.tpred$summary$mean, col = 'red', pch = 17)
legend('bottomleft', col = c('blue', 'red', 'blue', 'red'), 
       pch = c(16, 16, 17, 17), c('u PLS', 't PLS', 'u sPLS', 't sPLS'))

## ----
## PLS2 - canonical mode
## ----
tune.pls2 <- pls(X = X$methylation, Y = X$mRNA, ncomp = 2, mode = 'canonical')
Q2.pls <- perf(tune.pls2, validation = 'Mfold', 
               folds = 10, nrepeat = 5)
plot(Q2.pls, criterion = 'Q2.total')
list.keepX <- c(seq(10, 40, 5))
list.keepY <- c(seq(10, 40, 5))
tune.spls2 <- tune.spls(X$methylation, X$mRNA,
                        test.keepX = list.keepX, test.keepY = list.keepY,
                        ncomp = 2, 
                        nrepeat = 1, folds = 10, mode = 'canonical',measure = 'cor')
plot(tune.spls2)
choice.keepX <- tune.spls2$choice.keepX
choice.keepY <- tune.spls2$choice.keepY
choice.ncomp <- length(choice.keepX)
spls2 <- spls(X, Y,
              ncomp = choice.ncomp, keepX = choice.keepX, keepY = choice.keepY,
              mode = "canonical")
vip.spls2 <- vip(spls2)
head(vip.spls2[selectVar(spls2, comp = 1)$X$name,1])
perf.spls2 <- perf(spls2, validation = 'Mfold', folds = 10, nrepeat = 5)
plotIndiv(spls2, ind.names = FALSE, 
          group = , 
          pch = , 
          col.per.group = color.mixo(1:4),
          legend = TRUE, legend.title = 'Time', 
          legend.title.pch = 'Dose')
plotArrow(spls2, ind.names = FALSE, 
          group = ,
          col.per.group = color.mixo(1:4),
          legend.title = 'Time.Group')
plotVar(spls2, cex = c(3,4), var.names = c(FALSE, TRUE))
plot(c(1,2), perf.pls$measures$cor.upred$summary$mean, 
     col = 'blue', pch = 16, 
     ylim = c(0.6,1), xaxt = 'n',
     xlab = 'Component', ylab = 't or u Cor', 
     main = 's/PLS performance based on Correlation')
axis(1, 1:2)  # X-axis label
points(perf.pls$measures$cor.tpred$summary$mean, col = 'red', pch = 16)
points(perf.spls$measures$cor.upred$summary$mean, col = 'blue', pch = 17)
points(perf.spls$measures$cor.tpred$summary$mean, col = 'red', pch = 17)
legend('bottomleft', col = c('blue', 'red', 'blue', 'red'), 
       pch = c(16, 16, 17, 17), c('u PLS', 't PLS', 'u sPLS', 't sPLS'))

## ----
## Diablo
## ----
diablo <- block.plsda(X, Y, ncomp = 2, design = Cdesign)
perf.diablo = perf(diablo, validation = 'Mfold', folds = 10, nrepeat = 10)
plot(perf.diablo)
perf.diablo$choice.ncomp$WeightedVote
ncomp <- perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
test.keepX <- list(methylation = c(2:5, seq(10, 40, 5)),
                   mRNA = c(2:5, seq(10, 40, 5)))
tune.diablo <- tune.block.splsda(X, Y, ncomp = 2,
                                 test.keepX = test.keepX, design = Cdesign,
                                 validation = 'Mfold', folds = 10, nrepeat = 1,
                                 dist = "centroids.dist")

