#  2021-04-22 09:00 
#  elihei  [<eheidari@student.ethz.ch>]
# /Volumes/Projects/MS_lesions/analysis/sma02_novosparc_run.py


import argparse
import json
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import novosparc
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist, squareform, pdist
from scipy.stats import ks_2samp



# currentdir = os.path.dirname(os.path.abspath(__file__))
# parentdir = os.path.dirname(currentdir)
# sys.path.insert(0,parentdir) 


parser = argparse.ArgumentParser()
parser.add_argument('--data_dir', type=str, default='data/sma02_novosparc/', help='Data directory containing the expression matrix and the atlas.')
parser.add_argument('--out_dir', type=str, default='output/sma02_novosparc/', help='Output directory.')
parser.add_argument('--tag',type=str,default='EXP',help='The tag of the experiment (for saving).')
parser.add_argument('--atlas_locs_f', type=str,default='locs_atlas.txt', help='Path to the atlas locations file in the data directory. Should be in txt format.')
parser.add_argument('--atlas_mtx_f', type=str,default='mtx_atlas.txt', help='Path to the atlas marker expression matrix in the data directory. Should be in txt format.')
parser.add_argument('--expr_mtx_f',type=str,default='mtx_expr.txt',help='Path to the main expression matrix to be mapped to the atlas locations.')


parser.add_argument('--ncells',type=int,default=1000,help='Number of cells to be subsampled from the expression dataset.')
parser.add_argument('--nns',type=int,default=5, help='Num neighbors for cell-cell expression cost.')
parser.add_argument('--nnt',type=int,default=5, help='Num neighbors for location-location physical distance cost')

parser.add_argument('--alpha',type=float,default=0.8,help='The weight of the reference atlas.')
parser.add_argument('--epsilon',type=float,default=0.0005,help='Coefficient of entropy regularization')

parser.add_argument('--seed',type=int,default=0,help='Seed for generating the synthetic dataset.')

args = parser.parse_args()

target_space_path = os.path.join(args.data_dir, args.atlas_locs_f)
locations = pd.read_csv(target_space_path, sep=',')
locations_apriori = locations[['x', 'y']].values
locations = locations_apriori
atlas_path = os.path.join(args.data_dir, args.atlas_mtx_f)
atlas = sc.read(atlas_path, delimiter=',').T
atlas_gene_names = atlas.var.index.tolist()
atlas.obsm['spatial'] = locations

data_path = os.path.join(args.data_dir, args.expr_mtx_f)
dataset = sc.read(data_path, delimiter=',').T
gene_names = dataset.var.index.tolist()

ncells, ngenes = dataset.shape # 1297 cells x 8924 genes

print('number of cells: %d' % ncells)
print('number of genes: %d' % ngenes)

# optional: subset cells
# sc.pp.subsample(dataset, n_obs=args.ncells)
cell_names = dataset.obs.index.tolist()

tissue = novosparc.cm.Tissue(dataset=dataset, locations=locations_apriori)
cell_names = dataset.obs.index.tolist()
markers = list(set(atlas_gene_names).intersection(gene_names))
atlas_matrix = atlas.to_df()[markers].values
markers_idx = pd.DataFrame({'markers_idx': np.arange(ngenes)}, index=gene_names)
markers_to_use = np.concatenate(markers_idx.loc[markers].values)

# alternative 1: setup both assumptions 
tissue.setup_reconstruction(atlas_matrix=atlas_matrix,markers_to_use=markers_to_use, 
num_neighbors_s=args.nns, num_neighbors_t=args.nnt)
tissue.reconstruct(alpha_linear=args.alpha, epsilon=args.epsilon)

sdge = tissue.sdge
dataset_reconst = sc.AnnData(pd.DataFrame(sdge.T, columns=gene_names))
dataset_reconst.obsm['spatial'] = locations_apriori

gw = tissue.gw
ngw = (gw.T / gw.sum(1)).T
cell_idx = dataset_reconst.obs.index.tolist()
dataset_reconst.obs = pd.DataFrame(ngw.T, columns=cell_names)

sc.pp.highly_variable_genes(dataset)
is_var_gene = dataset.var['highly_variable']
var_genes = list(is_var_gene.index[is_var_gene])
tissue.calculate_spatially_informative_genes(var_genes)
mI, pvals = novosparc.analysis._analysis.get_moran_pvals(atlas.X, locations)
atlas_gene_names = atlas.var.index.tolist()
df = pd.DataFrame({'moransI': mI, 'pval': pvals}, index=atlas_gene_names)
genes_with_scores = tissue.spatially_informative_genes

out_path = os.path.join(args.out_dir, args.tag)
dataset_reconst.write_csvs(out_path, skip_data=False, sep=',')
scores_f = os.path.join(out_path, 'gene_scores.csv')
genes_with_scores.to_csv(scores_f)
