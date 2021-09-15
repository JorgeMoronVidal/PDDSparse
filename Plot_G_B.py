import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sparse
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigs
from scipy.linalg.interpolative import estimate_rank
import seaborn as sns
import pandas as pd
import os

sns.set()
sns.set_style("ticks")
print("Python-based PDDSparse G and B analysis plotter\n")
cwd = os.getcwd()
indir = os.getcwd() + "/Output/Debug";
os.chdir(indir)     
G = pd.read_csv('G.csv')
G_i = np.array(G['G_i'])
G_j = np.array(G['G_j'])
G_ij = np.array(G['G_ij'])
dim = G_i.max() + 1
coo = coo_matrix((G_ij,(G_i,G_j)),shape = (dim, dim))
#max_eigen, foo = eigs(coo,k = 1, which = 'LM')
#min_eigen, foo = eigs(coo,k = 1, which = 'SM')
coo = sparse.csc_matrix(coo)
norm_coo = sparse.linalg.norm(coo)
norm_inv_coo = sparse.linalg.norm(sparse.linalg.inv(coo))
cond = norm_coo*norm_inv_coo
fig = plt.figure(figsize = [20,8])
title = "Condition number: {:.2f}".format(cond)
title = title + "  Sparsity: {:.2f} %".format(100*G_i.size/(dim*dim))
title = title + "  Rank = {:}".format(np.linalg.matrix_rank(coo.todense()))
fig.suptitle(title)
ax1 = fig.add_subplot(1,3,1)
ax1.spy(coo)
ax1.set_xlabel("column")
ax1.set_ylabel("row")
ax1.set_title("Spy")
ax2 = fig.add_subplot(1,3,2)
p2 = ax2.scatter(G_i, G_j, c = G_ij, cmap = 'gnuplot2')
ax2.set_aspect('equal','box')
ax2.set_xlabel("column")
ax2.set_ylabel("row")
ax2.set_title("Solution")
fig.colorbar(p2, ax=ax2)
ax3 = fig.add_subplot(1,3,3)
G = pd.read_csv('G_var.csv')
G_i = np.array(G['G_i'])
G_j = np.array(G['G_j'])
var_G_ij = np.array(G['var_G_ij'])
p3 = ax3.scatter(G_i, G_j, c = var_G_ij, cmap = 'jet')
ax3.set_aspect('equal','box')
ax3.set_xlabel("column")
ax3.set_ylabel("row")
ax3.set_title("Variance ")
fig.colorbar(p3, ax=ax3)
os.chdir(cwd+"/Output/Plots") 
fig.savefig("G.png", dpi=200)

