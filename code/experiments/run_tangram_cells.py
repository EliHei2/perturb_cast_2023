import argparse
import json
import os, sys
import sys
import torch
import tangram as tg
import scanpy as sc
import tangram as tg
from utils import *

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
parser.add_argument('--cluster_label', type=str,default='cells', help='Output directory.')

args = parser.parse_args()
print(args)
n_hidden_GNN = [] if args.n_hidden_GNN==0 else [args.n_hidden_GNN]


if not os.path.exists(args.oo):
    os.makedirs(args.oo)
args.oo = os.path.join(args.oo, args.tag)
if not os.path.exists(args.oo):
    os.makedirs(args.oo)


if torch.cuda.is_available():  
  dev = "cuda:0" 
else:  
  dev = "cpu"  
device = torch.device(dev)

path_sp = os.path.join(args.i, args.tag, args.tag_ref) + '.h5ad'
path_sc = os.path.join(args.i, args.tag, args.tag_query) + '.h5ad'
ad_sp = sc.read_h5ad(path_sp)
ad_sc = sc.read_h5ad(path_sc)

if args.cluster_label != 'cells':
  ad_out = map_cells_to_space(
    adata_sc=ad_sc,
    adata_sp=ad_sp,
    device=device
  )
else:
  ad_out = map_cells_to_space(
    adata_sc=ad_sc,
    adata_sp=ad_sp,
    device=device,
    mode='clusters',
    cluster_label= args.cluster_label
  )

preds_f = "_".join(['preds', args.tag_train, args.tag_test, args.cluster_label]) + ".txt"
preds_f = os.path.join(args.oo, preds_f)
save_np_txt(ad_out.X, preds_f, colnames=dataset.classes)
