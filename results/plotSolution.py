import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
# Assumes 'results.txt' is on the form
# mpiprocs threads n error time
show = True
if("PLOT CONVERGENCE PART 1"):
    functionNumber = 1
    filenames = ['solution' + str(functionNumber) + '.txt']
    title = "Plot of solution for function" + str(functionNumber)
    for file in filenames:
        m = 1024-1
        u = np.fromfile(file, sep='\t', dtype = float) #np.empty((m,m))
        u= np.reshape(u,(m,m))
        x = np.linspace(0,1,m)
    y = np.linspace(0,1,m)
    X,Y = np.meshgrid(x,y)
    plt.imshow(u, extent = [0,1,0,1], cmap= cm.seismic)
    plt.colorbar()
    plt.title(title)
    plt.xlabel("x"); plt.ylabel("y")
    plt.savefig("function"+str(functionNumber)+".pdf")
    if show: plt.show()
