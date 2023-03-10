# Python Luke multimap_env
# Python Local multmp
import MultiMAP
import sys
import scanpy as sc
import pandas as pd
import numpy as np
from os import listdir
from rpy2 import robjects
from rpy2.robjects.packages import importr
base = importr("base")

if len(sys.argv) != 6:
    print("Arguments: <input folder for rna> <input folder for met> <output folder> <dataset lot nb> <nb of sims>")
    sys.exit()

# input_folder_rna = "/Users/elise/Desktop/Post-doc/Luke/projects/acacia/results/deconv_simple/rna/simu/"
# input_folder_met = "/Users/elise/Desktop/Post-doc/Luke/projects/acacia/results/deconv_simple/met/simu/"
input_folder_rna = sys.argv[1]
input_folder_met = sys.argv[2]
output_folder = sys.argv[3]
lot = int(sys.argv[4])
n_sims = int(sys.argv[5])
sc.pl.embedding(color_map=)

# load data
list_rna = [i for i in listdir(input_folder_rna) if ".rds" in i]
list_met = [i for i in listdir(input_folder_met) if ".rds" in i]
list_rna.sort()
list_met.sort()
readRDSrna = robjects.r('''readR = function(path) {readRDS(path)$D_rna_sim}''')
readRDSmet = robjects.r('''readR = function(path) {readRDS(path)$D_met_sim}''')
#readRDSprop = robjects.r('''readR = function(path) {readRDS(path)$A_ref}''')
D_rna = [readRDSrna(input_folder_rna+i) for i in list_rna if 'ref.rds' not in i][n_sims*(lot-1):n_sims*lot]
D_met = [readRDSmet(input_folder_met+i) for i in list_met if 'ref.rds' not in i][n_sims*(lot-1):n_sims*lot]
#Atrue = [readRDSprop(input_folder_rna+i) for i in list_rna if 'ref.rds' not in i][n_sims*(lot-1):n_sims*lot]
D_rna = [pd.DataFrame(np.asarray(i)).T for i in D_rna]
D_met = [pd.DataFrame(np.asarray(i)).T for i in D_met]
#Atrue = [pd.DataFrame(np.asarray(i)).T for i in Atrue]
#print([block.shape for block in Atrue])
print("data loaded")

# Transfo anndata
adata_rna = [sc.AnnData(i, dtype=np.asarray(i).dtype) for i in D_rna]
adata_met = [sc.AnnData(i, dtype=np.asarray(i).dtype) for i in D_met]
for i in adata_rna:
    i.obs['source'] = 'RNA'
for i in adata_met:
    i.obs['source'] = 'MET'
adata_rna_pca = [i.copy() for i in adata_rna]
for i in range(len(adata_rna)):
    sc.pp.scale(adata_rna_pca[i])
    sc.pp.pca(adata_rna_pca[i])
    adata_rna[i].obsm['X_pca'] = adata_rna_pca[i].obsm['X_pca'].copy()
    sc.pp.highly_variable_genes(adata_met[i], n_top_genes=20000, subset=True)
    #sc.pp.filter_genes(adata_met[i], min_counts=1)
adata_met_pca = [i.copy() for i in adata_met]
for i in range(len(adata_rna)):
    sc.pp.pca(adata_met_pca[i])
    adata_met[i].obsm['X_pca'] = adata_met_pca[i].obsm['X_pca'].copy()
print("data transformed")

# Run MultiMAP
adata = [MultiMAP.Integration([adata_rna[i], adata_met[i]], ['X_pca', 'X_pca'], scale=False) for i in range(len(adata_rna))]
#for i in range(len(adata)):
#    for celltype in range(Atrue[i].shape[0]):
#        adata[i].obs['celltype'+str(celltype)] = Atrue[i].iloc[celltype, :].values.tolist()
print("data mapped")

# Save
for i in range(len(adata)):
    adata[i].write_h5ad(output_folder+'dataset'+str(lot)+'_scaleF_sim0'+str(i+1)+'_adata.h5ad')