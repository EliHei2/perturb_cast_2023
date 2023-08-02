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

from gglasso.helper.data_generation import generate_precision_matrix, group_power_network, sample_covariance_matrix
from gglasso.problem import glasso_problem
from gglasso.helper.basic_linalg import adjacency_matrix
import numpy as np


from sklearn.covariance import empirical_covariance
from sklearn.metrics import *
from sklearn.covariance import GraphicalLassoCV, graphical_lasso, GraphicalLasso
from sklearn.preprocessing import StandardScaler
from scipy import sparse

from sagenet.utils import save_adata


import leidenalg as la
import igraph as ig


sc_adata = sg.MGA_data.scRNAseq()
sp_adata_1 = sg.MGA_data.seqFISH1_1()
sp_adata_2 = sg.MGA_data.seqFISH3_1()





def glasso(adata, alphas):
	N = adata.shape[1]
	scaler = StandardScaler()
	data = scaler.fit_transform(adata.X)
	S    = empirical_covariance(data)
	P    = glasso_problem(S, N, latent = False, do_scaling = True)
	# lambda1_range = np.logspace(-0.1, -1, 10)
	lambda1_range = np.logspace(-5,-0.1,10)
	mu1_range = np.logspace(-5,-0.1,10)
	modelselect_params = {'lambda1_range': lambda1_range, 'mu1_range': mu1_range}
	P.model_selection(modelselect_params = modelselect_params, method = 'eBIC', gamma = 0.5, tol=1e-7)
	sol = P.solution.precision_
	P.solution.calc_adjacency(t = 1e-6)
	save_adata(adata, attr='varm', key='adj', data=sparse.csr_matrix(P.solution.precision_))
	return P





sp_P_1 = glasso(sp_adata_1, [1, 2, 5])
sp_P_2 = glasso(sp_adata_2, [0.001, 0.1, 0.01])
sc_P = glasso(sc_adata, [0.001, 0.1, 0.01])




def gdiff_adj(adj_1, adj_2, ord='fro'):
	return np.linalg.norm((adj_1 - adj_2), ord)




import scanpy as sc
def partition_graph(adj):
	g = sc._utils.get_igraph_from_adjacency(adj)
	partition = la.find_partition(G, la.ModularityVertexPartition)
	return partition.membership


def plot_graph(adj):
	g = sc._utils.get_igraph_from_adjacency(adj)
	partition = la.find_partition(G, la.ModularityVertexPartition)
	ig.plot(partition)



from sklearn.metrics import jaccard_score
def gdiff_jaccard(adj_1, adj_2):
	part_1 = partition_graph(adj_1)
	part_2 = partition_graph(adj_2)
	return jaccard_score(part_1, part_2, average='micro')




gdiff_adj(sc_adata.varm['adj'].toarray(), sp_adata_2.varm['adj'].toarray())
gdiff_jaccard(sp_adata_1.varm['adj'].toarray(), sp_adata_1.varm['adj'].toarray())




import torch
if torch.cuda.is_available():  
  dev = "cuda:0" 
else:  
  dev = "cpu"  

device = torch.device(dev)
print(device)

sg_obj = sg.sage.sage(device=device)



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


