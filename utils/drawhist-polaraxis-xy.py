import matplotlib.pyplot as plt
import numpy as np
import sys

def drawHist(ax, heights,bounds=None, hnum=200,xlabel="x", ylabel="y",title=""):
    ax.hist(heights, hnum, align='mid',range=(0,180))
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
    ax.tick_params(direction="in", length=6, width=2, labelsize=20)
    ax.spines['bottom'].set_linewidth('2.0')
    ax.spines['top'].set_linewidth('2.0')
    ax.spines['left'].set_linewidth('2.0')
    ax.spines['right'].set_linewidth('2.0')


filename = sys.argv[1] 
_d = np.loadtxt(filename, skiprows=1).T
d = np.reshape(_d[5:8],-1)
fig = plt.figure()
ax1 = fig.add_subplot(111)
drawHist(ax1, d)

figfilename = "molecule-polaraxis-xy.png" 
plt.savefig(figfilename, bbox_inches = "tight")
plt.show()
