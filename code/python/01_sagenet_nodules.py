import os
import sagenet as sg
import scanpy as sc
import squidpy as sq
import anndata as ad
import random
import anndata as ad 
import re 
random.seed(10)
from scipy import sparse
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import pairwise_distances


data_path = 'data_tidy/00_initial_eda'

import torch
if torch.cuda.is_available():  
  dev = "cuda:0" 
else:  
  dev = "cpu"  

device = torch.device(dev)
print(device)

sg_obj = sg.sage.sage(device=device)

atlas_visium = sc.read_10x_mtx('data_raw/liver_atlas/rawData_mouseStStVisium/countTable_mouseStStVisium')




adata_dic = {}


# filter = df['Event Name'].str.contains(patternDel)

for adata_f in os.listdir(data_path):
	tag = re.sub('\\.h5ad', '', adata_f) 
	print(tag)
	adata_dic[tag] = sc.read(os.path.join(data_path, adata_f))
	normal_pat = '^Ctr|^Normal|^unassigned'
	normal_spots = adata_dic[tag].obs['nodule'].str.contains(normal_pat)
	adata_dic[tag] = adata_dic[tag][~normal_spots,:]
	print(adata_dic[tag].shape)
	adata_dic[tag].obs['section'] = tag
	adata_dic[tag] = adata_dic[tag][adata_dic[tag].X.sum(1) != 0, :]
	sc.pp.normalize_total(adata_dic[tag], target_sum=1e4)
	sc.pp.log1p(adata_dic[tag])
	adata_dic[tag].varm['adj'] = sparse.csr_matrix(adata_dic[tag].uns['adj'])
	le = LabelEncoder()
	adata_dic[tag].obs['nodule_encoded'] = le.fit_transform(adata_dic[tag].obs['nodule'])
	sg_obj.add_ref(adata_dic[tag], comm_columns=['nodule_encoded'], tag=tag, epochs=100, verbose = False)




sg_obj = sg.sage.sage(device=device)




os.makedirs('objects/sagenet_model_all_sections_nodules')
sg_obj.save_model_as_folder('objects/sagenet_model_all_sections_nodules')
# adata_q = ad.concat([adata_dic[i] for i in ['ML_II_A_2', 'ML_I']])
adata_q = ad.concat(adata_dic)
sc.pp.combat(adata_q, key='section')

sg_obj.load_model_as_folder('objects/sagenet_model_all_sections')

sg_obj.map_query(adata_q, save_prob=True)

# os.makedirs('data_tidy/01_sagenet_integration')

# 
# import anndata
# dist_adata = anndata.AnnData(adata_q.obsm['dist_map'], obs = adata_q.obs)
# knn_indices, knn_dists, forest = sc.neighbors.compute_neighbors_umap(dist_adata.X, n_neighbors=50, metric='precomputed')
# dist_adata.obsp['distances'], dist_adata.obsp['connectivities'] = sc.neighbors._compute_connectivities_umap(
#     knn_indices,
#     knn_dists,
#     dist_adata.shape[0],
#     50, # change to neighbors you plan to use
# )
# dist_adata.uns["neighbors"] = {"connectivities_key": "connectivities", "distances_key": "distances", "params": {"method": 'umap'}}

adata_q.write('data_tidy/01_sagenet_integration/integrated_query_nodules.h5ad')


adata_q = sc.read('data_tidy/01_sagenet_integration/integrated_query_nodules.h5ad')


import numpy as np

adata_q.obsm['prob_combined'] = np.column_stack((
		np.nan_to_num(adata_q.obsm['prob_ML_III_B_nodule_encoded'],0), 
		np.nan_to_num(adata_q.obsm['prob_ML_I_nodule_encoded'],0),
		np.nan_to_num(adata_q.obsm['prob_ML_II_A_2_nodule_encoded'], 0),
		np.nan_to_num(adata_q.obsm['prob_ML_II_A_1_nodule_encoded'], 0),
		np.nan_to_num(adata_q.obsm['prob_ML_III_A_nodule_encoded'], 0),
		np.nan_to_num(adata_q.obsm['prob_ML_II_B_nodule_encoded'], 0),
		np.nan_to_num(adata_q.obsm['prob_ML_I_2_nodule_encoded'], 0),
		np.nan_to_num(adata_q.obsm['prob_ML_II_C_nodule_encoded'], 0)

	)) 



base_dists = np.zeros((adata_q.obsm['prob_combined'].shape[0], adata_q.obsm['prob_combined'].shape[0]))
prob_list = ['prob_ML_III_A_nodule_encoded', 'prob_ML_III_B_nodule_encoded', 'prob_ML_II_A_1_nodule_encoded', 'prob_ML_II_A_2_nodule_encoded', 'prob_ML_II_B_nodule_encoded', 'prob_ML_II_C_nodule_encoded', 'prob_ML_I_2_nodule_encoded', 'prob_ML_I_nodule_encoded']
for prob in prob_list:
	print(prob)
	pd = pairwise_distances(adata_q.obsm[prob])
	del adata_q.obsm[prob]
	pd /= np.linalg.norm(pd, 'fro')
	base_dists += pd

	
adata_q.obsp['sagenet_dist'] = base_dists

import umap
my_model = umap.UMAP(metric='precomputed')
gg= my_model.fit(adata_q.obsp['sagenet_dist'])


sc.tl.leiden(adata_q, key_added = "leiden_raw")
# sc.tl.umap(dist_adata)

k=50
adata_q.obsp["distances"], adata_q.obsp["connectivities"] = sc.neighbors._compute_connectivities_umap(
    *sc.neighbors._get_indices_distances_from_dense_matrix(adata_q.obsp['sagenet_dist'], k),
    adata_q.obsp['sagenet_dist'].shape[0],
    k,
)
adata_q.uns["neighbors"] = {"connectivities_key": "connectivities", "distances_key": "distances", "params": {"method": None}}


sc.tl.leiden(adata_q, key_added = "leiden_sagenet")
sc.tl.umap(adata_q)


# adata_q.obsm['X_umap_sagenet'] = gg.embedding_


plasmid_genes = ['Olfr103', 'Olfr1', 'Olfr1018', 'Tas2r118', 'Vmn1r174', 'Olfr1002', 'Tas2r104', 'Vmn1r12', 'Olfr1006', 'Tas2r105', 'Olfr1008', 'Tas2r107', 'Olfr1014', 'Tas2r113', 'Olfr1013', 'Vmn1r170', 'Tas2r109', 'Vmn1r169', 'Olfr1012', 'Tas2r103', 'Vmn1r178', 'Olfr1019', 'Tas2r119', 'Vmn1r175']
sc.pp.neighbors(adata_q)
# adata_q_subset = copy(adata_q[:,plasmid_genes])
# sc.pp.neighbors(adata_q, use_rep='X')
sc.tl.umap(adata_q)
sc.pl.umap(adata_q, color='section', save='_integrated_nodules.pdf')
sc.pl.umap(adata_q, color=plasmid_genes, save='_integrated_nodules_genes.pdf')


adata_q.write('data_tidy/01_sagenet_integration/integrated_query_nodules_umapped.h5ad')




k=50
adata_q.obsp["distances"], adata_q.obsp["connectivities"] = sc.neighbors._compute_connectivities_umap(
    *sc.neighbors._get_indices_distances_from_dense_matrix(base_dists, k),
    base_dists.shape[0],
    k,
)
adata_q.uns["neighbors"] = {"connectivities_key": "connectivities", "distances_key": "distances", "params": {"method": None}}
sc.tl.umap(adata_q)
sc.pl.umap(adata_q, color='section', save='_integrated_all_additive_normalized.pdf')
sc.pl.umap(adata_q, color=plasmid_genes, save='_integrated_all_additive_normalized_genes.pdf')


from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, squareform


pairwise_distances(adata_q.obsm['prob_ML_III_B_nodule_encoded'])





sc.tl.tsne(adata_q, use_rep='prob_combined')
sc.pl.tsne(adata_q, color='section', save='_integrated_all_uncorrected_normalized.pdf')
sc.pl.tsne(adata_q, color=plasmid_genes, save='_integrated_all_uncorrected_normalized_genes.pdf')


from copy import copy
adata_q_subset = copy(adata_q[:,plasmid_genes])
sc.pp.neighbors(adata_q_subset, use_rep='X')
sc.tl.umap(adata_q_subset)
sc.pl.umap(adata_q_subset, color='section', save='_raw_all_subset.pdf')
sc.pl.umap(adata_q_subset, color=plasmid_genes, save='_raw_all_subset_genes.pdf')






vis_clus(
    spe = samples[[1]],
    clustervar = "nodule",
    colors = .palette_all, 
    size=0.5,
    alpha=0.5
) + theme(legend.position='none')



