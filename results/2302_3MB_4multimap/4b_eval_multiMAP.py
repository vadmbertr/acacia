import scanpy as sc
import os
import matplotlib.pyplot as plt

os.chdir("/Users/elise/Desktop/Post-doc/Luke/projects/acacia/results/2210_3MB_4multimap")

# load H5AD (scale FALSE)
list_data_lot1 = [i for i in os.listdir("4a_multiMAP/") if ".h5ad" in i if "scaleT" in i if "dataset1" in i]
list_data_lot2 = [i for i in os.listdir("4a_multiMAP/") if ".h5ad" in i if "scaleT" in i if "dataset2" in i]
list_data_lot3 = [i for i in os.listdir("4a_multiMAP/") if ".h5ad" in i if "scaleT" in i if "dataset3" in i]
list_data_lot1.sort()
list_data_lot2.sort()
list_data_lot3.sort()
data_lot1 = [sc.read("4a_multiMAP/"+i) for i in list_data_lot1]
data_lot2 = [sc.read("4a_multiMAP/"+i) for i in list_data_lot2]
data_lot3 = [sc.read("4a_multiMAP/"+i) for i in list_data_lot3]

# connect samples from each source
for i in range(len(data_lot1)):
    print(i)
    data_lot1[i].obs['SampleNb'] = [*range(int(data_lot1[i].shape[0]/2))]*2
    data_lot2[i].obs['SampleNb'] = [*range(int(data_lot2[i].shape[0]/2))]*2
    data_lot3[i].obs['SampleNb'] = [*range(int(data_lot3[i].shape[0]/2))]*2

# plot multimap vs pca
for i in data_lot1:
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    ax1_dict = sc.pl.embedding(i, 'X_pca', color='source', title="lot1", ax=ax1, palette='autumn', show=False)
    ax2_dict = sc.pl.embedding(i, 'X_multimap', color='source', title="lot1", ax=ax2, palette='autumn', show=False)
    ax3_dict = sc.pl.embedding(i, 'X_multimap', color='SampleNb', title="lot1", ax=ax3, color_map='viridis', show=False)
    ax4_dict = sc.pl.embedding(i, 'X_multimap', color='type', title="lot1", ax=ax4, palette='Spectral', show=False)
    fig.show()
for i in data_lot2:
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    ax1_dict = sc.pl.embedding(i, 'X_pca', color='source', title="lot1", ax=ax1, palette='autumn', show=False)
    ax2_dict = sc.pl.embedding(i, 'X_multimap', color='source', title="lot1", ax=ax2, palette='autumn', show=False)
    ax3_dict = sc.pl.embedding(i, 'X_multimap', color='SampleNb', title="lot1", ax=ax3, color_map='viridis', show=False)
    ax4_dict = sc.pl.embedding(i, 'X_multimap', color='type', title="lot1", ax=ax4, palette='Spectral', show=False)
    fig.show()
for i in data_lot3:
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    ax1_dict = sc.pl.embedding(i, 'X_pca', color='source', title="lot1", ax=ax1, palette='autumn', show=False)
    ax2_dict = sc.pl.embedding(i, 'X_multimap', color='source', title="lot1", ax=ax2, palette='autumn', show=False)
    ax3_dict = sc.pl.embedding(i, 'X_multimap', color='SampleNb', title="lot1", ax=ax3, color_map='viridis', show=False)
    ax4_dict = sc.pl.embedding(i, 'X_multimap', color='type', title="lot1", ax=ax4, palette='Spectral', show=False)
    fig.show()

# save to appropriate format for R conversion
for i in range(len(data_lot1)):
    data_lot1[i].write_h5ad('4a_multiMAP/'+'dataset1'+'_scaleT_sim0'+str(i+1)+'_adata.h5ad')
    data_lot2[i].write_h5ad('4a_multiMAP/'+'dataset2'+'_scaleT_sim0'+str(i+1)+'_adata.h5ad')
    data_lot3[i].write_h5ad('4a_multiMAP/'+'dataset3'+'_scaleT_sim0'+str(i+1)+'_adata.h5ad')
# remove former sim10 files and rename new sim010 files