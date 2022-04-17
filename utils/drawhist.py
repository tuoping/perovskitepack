import matplotlib.pyplot as plt
import numpy as np
import sys

def drawHist(ax, heights,bounds=None, hnum=200,xlabel="x", ylabel="y",title="", xlimit = None, ylimit=None):
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


filename = sys.argv[1] 
d = np.loadtxt(filename).T
d0 = d[0]
xlimit = None
fig = plt.figure()
ax1 = fig.add_subplot(311)
drawHist(ax1, d0, xlimit)

d1 = d[1]
ax2 = fig.add_subplot(312)
drawHist(ax2, d1, xlimit)

d2 = d[2]
ax3 = fig.add_subplot(313)
drawHist(ax3, d2, xlimit)

# figfilename = filename.replace(".dat", ".png")
figfilename = "molcorr.png" 
plt.savefig(figfilename, bbox_inches = "tight")
plt.show()
