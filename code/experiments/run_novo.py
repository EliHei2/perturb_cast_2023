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
import anndata as adata


currentdir = os.path.dirname(os.path.abspath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 



parser = argparse.ArgumentParser()
parser.add_argument('--tag', type=str, default='Syn', help='Input dataset. Can either be the name of the synthetic dataset generation method, \
                                              or the path to a real dataset.')
parser.add_argument('--tag_ref', type=str, default='train', help='Input dataset. Can either be the name of the synthetic dataset generation method, \
                                              or the path to a real dataset.')
parser.add_argument('--tag_query', type=str, default='test', help='Input dataset. Can either be the name of the synthetic dataset generation method, \
                                              or the path to a real dataset.')
parser.add_argument('-i', type=str,default='../../data_tidy', help='Input dataset. Can either be the name of the synthetic dataset generation method, \
                                              or the path to a real dataset.')
parser.add_argument('--oo', type=str,default='../../output', help='Output directory.')


args = parser.parse_args()
print(args)


if not os.path.exists(args.oo):
    os.makedirs(args.oo)
args.oo = os.path.join(args.oo, args.tag)
if not os.path.exists(args.oo):
    os.makedirs(args.oo)
args.oo = os.path.join(args.oo, 'novosparc')
if not os.path.exists(args.oo):
    os.makedirs(args.oo)

path_sp = os.path.join(args.i, args.tag, args.tag_ref) + '.h5ad'
path_sc = os.path.join(args.i, args.tag, args.tag_query) + '.h5ad'

print(path_sc)
ad_sp = sc.read_h5ad(path_sp)
ad_sc = sc.read_h5ad(path_sc)
gene_names = ad_sp.var.index.tolist()
num_cells, num_genes = ad_sp.shape 


tissue = novosparc.cm.Tissue(dataset=ad_sc, locations=ad_sp.obs[['x', 'y']])
tissue.setup_reconstruction(atlas_matrix=ad_sp.X)
tissue.reconstruct(alpha_linear=0.5, epsilon=5e-3)

gw = tissue.gw
ngw = (gw.T / gw.sum(1)).T

ad_out = adata.AnnData(X = ngw, var = ad_sp.obs, obs=ad_sc.obs)
preds_f = "_".join(['preds', args.tag_ref, args.tag_query]) + ".h5ad"
preds_f = os.path.join(args.oo, preds_f)
# save_np_txt(y_pred, preds_f, colnames=dataset.classes)
ad_out.write(filename=preds_f)


# preds_f = "_".join(['preds', args.tag_ref, args.tag_query]) + ".txt"
# preds_f = os.path.join(args.oo, preds_f)


# def save_np_txt(ndarray, path, colnames=None, rownames=None):
#     df = pd.DataFrame(data=ndarray, index=rownames, columns=colnames)
#     df.to_csv(path, sep='\t', index=True, header=True)

# save_np_txt(ngw, preds_f)