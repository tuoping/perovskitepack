import matplotlib.pyplot as plt
import numpy as np
import sys

def drawHist(ax, heights,bounds=None, hnum=50,xlabel="x", ylabel="y",title="", xlimit = None, ylimit=None):
    ax.hist(heights, hnum, align='mid',range=bounds)
    # plt.hist(heights, hnum, align='mid',range=bounds)
    # plt.title(title)
    # plt.xlabel(xlabel, fontsize = 16)
    # plt.ylabel(ylabel, fontsize = 16)
    # plt.tick_params(direction="in")
    font={'family':'serif',
          # 'style':'italic',  # 斜体
          'weight':'normal',
          # 'color':'red',
          'size': 18
    }
    ax.tick_params(direction="in")
    if xlimit is not None:
        plt.xlim(xlimit)
    if ylimit is not None:
        plt.ylim(ylimit)


dirname = "./"
caxis = int(np.loadtxt(dirname+"/caxis"))
mesh_size = [10,10,10]

Totaltime = 1
Num_oct = mesh_size[0]*mesh_size[1]*mesh_size[2]
data = np.loadtxt(dirname+"/ii_vectors_caxis.dat")
data = np.reshape(data, [-1,Num_oct,3,2])
Totaltime = len(data)
print("Totaltime = ", Totaltime)
phiy= np.abs(data[:,:,1,1])

xlimit = None
fig = plt.figure()
ax1 = fig.add_subplot(111)
drawHist(ax1, phiy[0], xlimit)

# figfilename = filename.replace(".dat", ".png")
figfilename = "hist_ii_vector.png" 
plt.savefig(figfilename, bbox_inches = "tight")
plt.show()
