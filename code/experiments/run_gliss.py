from IPython.display import display
import sys
import os
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi']= 100
import csv
sys.path.insert(0,'src')
from code.python.GLISS.io_utils import load_data_from_file
from code.python.GLISS.general_utils import norm_mtx
from code.python.GLISS.plot_utils import plot_multiple_scatter_continuous
from code.python.GLISS.main_methods import select_spatial_genes

# locs = load_data_from_file('data/embryo1_2/data_train.txt', 'txt')
feature_names = load_names('data/embryo1_5')
data_train    = load_features('data/embryo1_5', type='train')
data_test   = load_features('data/embryo1_5', type='test')
id_train, loc_train = load_classes('data/embryo1_5', type='train')
id_test, loc_test   = load_classes('data/embryo1_5', type='test')


data_norm = norm_mtx(data_train) # normalize each gene
print('Gene expression matrix dimension: {}'.format(data_norm.shape))

# paramters
alpha= 0.05 # FDR Level
knn = 5 # number of top neighbors in constructing KNN graph
n_perm = 1000 # number of permutation for p-values
# run SV gene selection with GLISS
pvals, rej_idx = select_spatial_genes(loc_train, data_norm, 
                         knn=knn, alpha=alpha, n_perm=n_perm)

x_star = data_test[:, lm_genes_idx] # lm gene matrix
print('LM gene matrix dimension: {}'.format(x_star.shape))
x = np.delete(data_test, lm_genes_idx, axis=1) # non-lm gene matrix
print('non-LM gene matrix dimension: {}'.format(x.shape))

pvals, rej_idx = select_spatial_genes(x_star, x, 
                         knn=knn, alpha=alpha, n_perm=n_perm)

from sim_utils import infer_lambda
new_x = np.concatenate([x_star, x], axis=1)[:,rej_idx]
lam, aux = infer_lambda('graph', new_x, knn=knn)

save_np_txt(lam, 'output/knn_embryo1_5.txt')