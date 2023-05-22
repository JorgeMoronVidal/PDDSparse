import numpy as np 
import pandas as pd 
from scipy.sparse import coo_matrix
def 
#G matrix is charged from file
G_table = pd.read_csv('Output/Debug/G.csv')
G = coo_matrix((G_table['G_ij'],(G_table['G_i'],G_table['G_j'])))
G = coo_matrix.tocsr(G);
#G matrix original condition numbers
print("Norm 2 condition number "+str(np.linalg.cond(G.todense())))
print("Norm infinite condition number "+str(np.linalg.cond(G.todense(),np.inf)))
#Estimation of delta



