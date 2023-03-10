# Python Luke multimap_env
# Python Local multimap3.7
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

input_folder_rna = sys.argv[1]
input_folder_met = sys.argv[2]
output_folder = sys.argv[3]
lot = int(sys.argv[4])
n_sims = int(sys.argv[5])


# load data
readRDSrna = robjects.r('''readR = function(path) {readRDS(path)$D_rna_sim}''')
readRDSmet = robjects.r('''readR = function(path) {readRDS(path)$D_met_sim}''')
readRDS_T = robjects.r('''readR = function(path) {readRDS(path)}''')
list_rna = [i for i in listdir(input_folder_rna) if ".rds" in i]
list_met = [i for i in listdir(input_folder_met) if ".rds" in i]
list_rna.sort()
list_met.sort()
list_rna1 = [i for i in list_rna if 'ref.rds' not in i][n_sims*(lot-1):n_sims*lot]
list_met1 = [i for i in list_met if 'ref.rds' not in i][n_sims*(lot-1):n_sims*lot]
D_rna = [readRDSrna(input_folder_rna+i) for i in list_rna1]
T_rna = [readRDS_T(input_folder_rna+i) for i in list_rna if 'ref.rds' in i][lot-1]
D_met = [readRDSmet(input_folder_met+i) for i in list_met1]
T_met = [readRDS_T(input_folder_met+i) for i in list_met if 'ref.rds' in i][lot-1]
D_rna = [pd.DataFrame(np.asarray(i)).T for i in D_rna]
T_rna = pd.DataFrame(np.asarray(T_rna))
D_met = [pd.DataFrame(np.asarray(i)).T for i in D_met]
T_met = pd.DataFrame(np.asarray(T_met))
if T_rna.shape[0]>T_rna.shape[1]:
    T_rna = T_rna.T
if T_met.shape[0]>T_met.shape[1]:
    T_met = T_met.T
n_sample = D_rna[0].shape[0]
n_ref = T_rna.shape[0]

data_rna = [pd.concat([i, T_rna], ignore_index=True, sort=False) for i in D_rna]
del D_rna, T_rna
if lot == 3:
    adata_met = [sc.AnnData(i, dtype=np.asarray(i).dtype) for i in D_met]
    del D_met
    [sc.pp.highly_variable_genes(i, n_top_genes = int(5e5)) for i in adata_met]
    D_met = [pd.DataFrame(i.X[:, i.var.iloc[:, 0].values]) for i in adata_met]
    T_met = [pd.DataFrame(np.asarray(T_met.iloc[:, i.var.iloc[:, 0].values])) for i in adata_met]
    del adata_met
if lot == 3:
    data_met = [pd.concat([i, j], ignore_index=True, sort=False) for (i, j) in zip(D_met, T_met)]
else:
    data_met = [pd.concat([i, T_met], ignore_index=True, sort=False) for i in D_met]
del D_met, T_met
print([block.shape for block in data_rna])
print([block.shape for block in data_met])
print("data loaded")

# Transfo anndata
adata_rna = [sc.AnnData(i, dtype=np.asarray(i).dtype) for i in data_rna]
adata_met = [sc.AnnData(i, dtype=np.asarray(i).dtype) for i in data_met]
for i in adata_rna:
    i.obs['source'] = 'RNA'
    i.obs['type'] = ['sample']*n_sample + ['ref']*n_ref
for i in adata_met:
    i.obs['source'] = 'MET'
    i.obs['type'] = ['sample']*n_sample + ['ref']*n_ref
adata_rna_pca = [i.copy() for i in adata_rna]
for i in range(len(adata_rna)):
    print(str(i))
    sc.pp.scale(adata_rna_pca[i])
    sc.pp.pca(adata_rna_pca[i])
    adata_rna[i].obsm['X_pca'] = adata_rna_pca[i].obsm['X_pca'].copy()
    sc.pp.highly_variable_genes(adata_met[i], n_top_genes=20000, subset=True)
adata_met_pca = [i.copy() for i in adata_met]
for i in range(len(adata_rna)):
    print(str(i))
    sc.pp.pca(adata_met_pca[i])
    adata_met[i].obsm['X_pca'] = adata_met_pca[i].obsm['X_pca'].copy()
print("data transformed")

# Run MultiMAP
adata = [MultiMAP.Integration([adata_rna[i], adata_met[i]], ['X_pca', 'X_pca'], scale=True, n_components=10) for i in range(len(adata_rna))]
print("data mapped")
for i in adata:
    del i.uns
    del i.obsp

# Save
for i in range(len(adata)):
    adata[i].write_h5ad(output_folder+'dataset'+str(lot)+'_scaleT_sim0'+str(i+1)+'_adata.h5ad')
# remove former sim10 files and rename new sim010 files