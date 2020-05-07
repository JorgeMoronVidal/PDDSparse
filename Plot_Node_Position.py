import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

#Clean Plot style
sns.set()
sns.set_style("ticks")
print("It's me, the python-based plotter\n")
cwd = os.getcwd()
indir = os.getcwd() + '/Output/Debug'
outdir = os.getcwd() + '/Output/Plots'
try:
   os.chdir(indir)     
   file = 'Node_position.txt'
   print('Reading ' + file + '\n')
   #Import data from file
   table = pd.read_csv(file)
   x = list(table['x'])
   y = list(table['y'])
   labels = list(table['index'])
   fig, ax = plt.subplots()
   ax.scatter(x, y, s = 5)
   for i, txt in enumerate(labels):
        ax.annotate(txt, (x[i], y[i]))
   plt.show()
except:
    print("Node_position.txt couldn't be properly plotted. \n")
else:
    print("Node_position.txt was succesfully plotted.\n")