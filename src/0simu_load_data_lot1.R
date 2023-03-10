if (lot=="lot1") {
  # load lot1 RNA raw data
  rna_1 <- read.delim("../../../datashare/data_cometh_lot1/transcriptome/gene_counts_filtered.txt", head=T, row.names=1)
  # load lot1 MET beta val data
  met_1 <- readRDS("../../../datashare/data_cometh_lot1/methylation/test_data_met_new.rds")
  # load lot1 ground truth
  gt <- readRDS("../../../datashare/data_cometh_lot1/test_solution.rds")[[1]]
  # prepare data list and inputs
  data = list("rna"=as.matrix(rna_1),
              "met"=met_1)
  rm(rna_1, met_1)
}

if (lot=="lot2") {
  # load lot2 RNA raw data
  rna_2 <- read.delim("../../../datashare/data_cometh_lot2/transcriptome/gene_counts_filtered.txt", head=T, row.names=1) # 2 more samples to rmv (8 and 10)
  # load lot2 MET beta val data
  met_2 <- readRDS("../../../datashare/data_cometh_lot2/methylation/test_data_met_new.rds")
  # load lot2 owkin ground truth
  gt <- read.csv("../../../datashare/data_cometh_lot2/prediction_owkin/preds_clean_50_add_mr.csv")[,3:7]
  rownames(gt) <- gt$ID_QCPF
  gt <- gt[,-1]
  gt <- gt[colnames(rna_2),]
  # load lot2 other ground truth
  gt2 <- read.csv("../../../datashare/data_cometh_lot2/prediction_owkin/FittedComponentsBj100.tsv", sep="\t", row.names = 1)[,-c(7,8)]
  gt2 <- gt2[grep("COM",rownames(gt2)),]
  gt2 <- gt2[colnames(rna_2),]
  # prepare data list and inputs for the runfactorization script from momix
  data = list("rna"=as.matrix(rna_1),
              "met"=met_1)
  rm(rna_2, met_2)
}
