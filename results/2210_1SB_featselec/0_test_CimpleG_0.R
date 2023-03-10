library(CimpleG)

## Params
block='met'
deconv_methods_sup <- c("rlr","CIBERSORT","houseman")#,"EMeth")
deconv_methods_unsup <- c("RefFreeEWAS","ICA","EDec","MeDeCom")
source(paste0("../../src/1SB_deconv_met_functions.R"))

## Data
data(train_data)
data("train_targets")
dat <- t(train_data)
annot <- train_targets[,c(1,2)]
T_ref <- t(WGCNA::collapseRows(t(dat), annot$CELL_TYPE,  rownames(t(dat)),method="Average")$datET)
Atrue <- t(train_targets[,c(3:16)])
rownames(Atrue) <- gsub("CELL_TYPE_","",rownames(Atrue))
colnames(Atrue) <- annot$GSM
Atrue <- Atrue[rowSums(Atrue)>0,]
T_ref <- T_ref[sort(rownames(T_ref)),]
T_ref <- T_ref[,sort(colnames(T_ref))]
Atrue <- Atrue[sort(rownames(Atrue)),]
Atrue <- Atrue[,sort(colnames(Atrue))]
dat <- dat[sort(rownames(dat)),]
dat <- dat[,sort(colnames(dat))]
rm(train_data, train_targets, annot)

## Deconvolution
tmp <- lapply(deconv_methods_sup, function(meth) {
    if (!file.exists(paste0("0_test_CimpleG/deconv/deconvbl_",block,"_sup_",
                            meth,".rds"))) {
      print(paste0("Running ",meth))
      deconv_res=run_sup_deconvolution(meth, dat, T_ref)
      saveRDS(deconv_res,paste0("0_test_CimpleG/deconv/deconvbl_",block,"_sup_",
                                meth,".rds"))
    }
  })
  
tmp <- lapply(deconv_methods_unsup, function(meth) {
  if (!file.exists(paste0("0_test_CimpleG/deconv/deconvbl_",block,"_unsup_",
                          meth,".rds"))) {
    print(paste0("Running ",meth))
    deconv_res=run_unsup_deconvolution(meth, dat, T_ref, option="Amat", prop_simu=Atrue)
    saveRDS(deconv_res,paste0("0_test_CimpleG/deconv/deconvbl_",block,"_unsup_",
                              meth,".rds"))
  }
})
