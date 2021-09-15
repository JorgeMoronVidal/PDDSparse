import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import seaborn as sns
import pandas as pd
import os
def distance(parameters, X, Y):
    dist = []
    dists = [0.0]*4
    for i in range(X.size):
        dists = [0.0]*4

        dists[0] = abs(parameters[0] - X[i])
        dists[1] = abs(parameters[1] - Y[i])
        dists[2] = abs(parameters[2] - X[i])
        dists[3] = abs(parameters[3] - Y[i])
        dist.append(min(dists))
    return dist

def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

sns.set()
sns.set_style("ticks")
print("Python-based PDDSparse knot results plotter\n")
cwd = os.getcwd()
indir = os.getcwd() + "/Output";
os.chdir(indir)     
table = pd.read_csv('solution.csv')
X = np.array(table['x'])
Y = np.array(table['y'])
Z_num = np.array(table['sol_PDDS'])
Z_anal = np.array(table['sol_analytic'])
Z_err = np.array(table['err'])
Z_rerr = np.array(table['rerr'])
fig = plt.figure(figsize = [20,15])
fig.suptitle('Knot solution 10x10 N = 10000 with control variates F')
ax1 = fig.add_subplot(2, 2, 1)
ax2 = fig.add_subplot(2, 2, 2)
ax3 = fig.add_subplot(2, 2, 3)
ax4 = fig.add_subplot(2, 2, 4)
p1 = ax1.scatter(X,Y, c = Z_anal, cmap = 'plasma', marker = 's', s = 45)
p2 = ax2.scatter(X,Y, c = Z_num, cmap = 'plasma', marker = 's', s = 45)
p3 = ax3.scatter(X,Y, c = Z_err, cmap = 'viridis', marker = 's',s = 45)
p4 = ax4.scatter(X,Y, c = Z_rerr, cmap = 'viridis',marker = 's', s = 45)
ax1.set_title("Analytic Solution")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
colorbar(p1)
ax2.set_title("Numeric Solution")
ax2.set_xlabel("X")
ax2.set_ylabel("Y")
colorbar(p2)
ax3.set_title("Error")
ax3.set_xlabel("X")
ax3.set_ylabel("Y")
colorbar(p3)
ax4.set_title("Relative Error")
ax4.set_xlabel("X")
ax4.set_ylabel("Y")
colorbar(p4)
os.chdir(cwd+"/Output/Plots") 
fig.savefig("knot_solution.png", dpi=200)
fig = plt.figure(figsize = [21,7])
print("Python-based PDDSparse error analysis plotter\n")
fig.suptitle('Error analysis 10 x 10 N = 10000 with Variance Reduction')
ax1 = fig.add_subplot(1, 2, 1)
n, bins, patches = ax1.hist(Z_err, 50, facecolor = 'blue')
title = "Error\n "
title = title + "Average = {:.4f}  ".format(np.average(Z_err))
title = title + "Maximum = {:.4f}  ".format(np.max(Z_err))
title = title + "minimum = {:.4f}".format(np.min(Z_err))
ax1.set_title(title)
ax1.set_xlabel("Value")
ax1.set_ylabel("Frequency")
#ax1.text(60, .025, r'$\mu=100,\ \sigma=15$')
ax2 = fig.add_subplot(1, 2, 2)
n, bins, patches = ax2.hist(Z_rerr, 50, facecolor = 'blue')
title = "Relative error\n "
title = title + "Average = {:.4f}  ".format(np.average(Z_rerr))
title = title + "Maximum = {:.4f}  ".format(np.max(Z_rerr))
title = title + "minimum = {:.4f}".format(np.min(Z_rerr))
ax2.set_title(title)
ax2.set_xlabel("Value")
ax2.set_ylabel("Frequency")
parameters = [0.0]*4
parameters[0] = min(X)
parameters[1] = min(Y)
parameters[2] = max(X)
parameters[3] = max(Y)
fig.savefig("error_histograms.png", dpi=200)
fig = plt.figure(figsize = [20,8])
dist = distance(parameters, X, Y)
ax1 = fig.add_subplot(1, 2, 1)
p1 = ax1.scatter(X,Y, c = abs(Z_err), cmap = 'plasma', marker = 's', s = 45)
ax1.set_aspect('equal','box')
ax1.set_title("Absolute error")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
colorbar(p1)
ax2 = fig.add_subplot(1,2,2)
ax2.scatter(dist, abs(Z_err), c = 'blue')
ax2.set_aspect('auto','box')
ax2.set_title("Error knot-distance dependence")
ax2.set_xlabel("Distance")
ax2.set_ylabel("Error")
fig.savefig("error_distance.png", dpi=200)
table_debug = pd.read_csv(cwd+'/Output/Debug/Node_debug.csv')
var_B = [x for _, x in sorted(zip(table_debug["index"],table_debug["var_B"]))]
Knot_APL = [x for _, x in sorted(zip(table_debug["index"],table_debug["APL"]))]
Knot_time = [x for _, x in sorted(zip(table_debug["index"],table_debug["time"]))]
fig = plt.figure(figsize = [30,8])
ax1 = fig.add_subplot(1, 3, 1)
p1 = ax1.scatter(X,Y, c = var_B, cmap = 'plasma', marker = 's', s = 45)
ax1.set_aspect('equal','box')
ax1.set_title("Variance of B")
ax1.set_xlabel("X")
ax1.set_ylabel("Y")
colorbar(p1)
ax2 = fig.add_subplot(1, 3, 2)
p2 = ax2.scatter(X,Y, c =Knot_APL, cmap = 'inferno', marker = 's', s = 45)
ax2.set_aspect('equal','box')
ax2.set_title("Average Path Length")
ax2.set_xlabel("X")
ax2.set_ylabel("Y")
colorbar(p2)
ax3 = fig.add_subplot(1, 3, 3)
p3 = ax3.scatter(X,Y, c =Knot_time, cmap = 'cividis', marker = 's', s = 45)
ax3.set_aspect('equal','box')
ax3.set_title("Wall-clock time")
ax3.set_xlabel("X")
ax3.set_ylabel("Y")
colorbar(p3)
fig.savefig("node_analysis.png", dpi=200)