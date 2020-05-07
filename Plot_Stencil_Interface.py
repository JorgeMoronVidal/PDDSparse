import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re
import os
#Potential distribution's equation 
def Pot (h,c,e):
    return c*(h**e)

print("\n")

#Clean Plot style
sns.set()
sns.set_style("ticks")

cwd = os.getcwd()
indir = os.getcwd() + '/Output/Debug'
outdir = os.getcwd() + '/Output/Plots'
try: 
    os.mkdir(outdir)
except: 
    try: 
        os.chdir(outdir)
    except:
        OSError
        print("Output directory %s couldn't be created. \n" % outdir)
    else:
        print("Output directory %s already existed \n" % outdir)
        os.chdir(cwd)
    
try:
    os.chdir(indir)     
    for file in os.listdir(indir):
        if(('Interface_') in file):
            print('Reading ' + file + '\n')
            #Import data from file
            interface = pd.read_csv(file)
            file_s =  re.sub('Interface', 'stencil', file)
            stencil = pd.read_csv(file_s)
            plt.figure()
            plt.title(re.sub('.txt|_', ' ', file))
            plt.xlabel('x')
            plt.ylabel('y')
            x = list(interface['x'])
            y = list(interface['y'])
            labels = list(interface['index'])
            legend = 'Interface'
            plt.scatter(x,y, s = 5, label = legend)
            for i, txt in enumerate(labels):
                plt.annotate(txt, (x[i], y[i]))
            x = list(stencil['x'])
            y = list(stencil['y'])
            labels = list(stencil['index'])
            legend = 'stencil'
            plt.scatter(x,y, s = 5, label = legend)
            for i, txt in enumerate(labels):
                plt.annotate(txt, (x[i], y[i]))
            plt.legend()
            os.chdir(outdir)
            plt.savefig( file +'.png',dpi=200)
            plt.close()
            os.chdir(indir)
           
    
except:
    print("Input directory %s couldn't be plotted. \n" %indir)
else:
    print("Input directory %s was succesfully plotted.\n" %indir)


