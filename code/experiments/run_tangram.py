import argparse
import json
import os, sys
import sys
import torch
import tangram as tg
import scanpy as sc
import pandas as pd

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


if not os.path.exists(args.oo):
    os.makedirs(args.oo)
args.oo = os.path.join(args.oo, args.tag)
if not os.path.exists(args.oo):
    os.makedirs(args.oo)
args.oo = os.path.join(args.oo, ('tangram_' + args.cluster_label))
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
tg.pp_adatas(ad_sc, ad_sp, genes=None)

preds_f = "_".join(['preds', args.tag_ref, args.tag_query]) + ".h5ad"
preds_f = os.path.join(args.oo, preds_f)

if args.cluster_label == 'cells':
  ad_out = tg.map_cells_to_space(
    ad_sc,
    ad_sp,
    device=device
  )
  ad_out.T.write(filename=preds_f)
else:
  ad_out = tg.map_cells_to_space(
    ad_sp,
    ad_sc,
    device=device,
    mode='clusters',
    cluster_label= args.cluster_label
  )
  ad_out.T.write(filename=preds_f)


# def save_np_txt(ndarray, path, colnames=None, rownames=None):
#     df = pd.DataFrame(data=ndarray, index=rownames, columns=colnames)
#     df.to_csv(path, sep='\t', index=True, header=True)
# save_np_txt(ad_out.X.T, preds_f, rownames=ad_out.var.index, colnames=ad_out.obs.index)


