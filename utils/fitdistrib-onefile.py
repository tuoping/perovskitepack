import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys

def startfig(size = (5,5)):
    plt.figure()
    plt.rcParams["figure.figsize"] = size
    plt.rcParams['axes.linewidth'] =2.0
    plt.rcParams['xtick.major.width'] =2.0
    plt.rcParams['ytick.major.width'] =2.0

def gaussian(x, mu, sigma):
    return np.exp(-(x-mu)*(x-mu)/2/sigma/sigma)

def function(x, mu1, mu2, sigma1, sigma2, A1, A2):
    return A1*gaussian(x, mu1, sigma1) + A2*gaussian(x, mu2, sigma2)

def function_1G(x, mu1, sigma1, A1):
    return A1*gaussian(x, mu1, sigma1)

def Linearfunction(x,A,B):
    return A*x+B

T = int(sys.argv[1])
ffit = open("centers.dat", "w")
dirname = "./"
_data = np.concatenate((np.loadtxt(dirname+"/left_distances"),np.loadtxt(dirname+"/right_distances")))
data = pd.DataFrame(_data, columns=["number"])
data["cut_group"] = pd.cut(data["number"], 25)
x = [d.mid for d in data["cut_group"].value_counts().index]
y = [d for d in data["cut_group"].value_counts()]
plt.figure()
plt.scatter(x, y)
    
if T < 195:
    popt, pcov = curve_fit(function, x, y, bounds=([8, 8.5, 0.1, 0.1, 1000, 1000], [9, 9.5, 0.5, 0.5, 5000, 5000]))
    ffit.write("%d    %f  %f\n"%(T, popt[0], popt[1]))
    y_pred = [function(i,*popt) for i in x]
else:
    popt, pcov = curve_fit(function_1G, x, y, bounds=([8, 0.1, 1000], [9.5, 0.5, 5000]))
    ffit.write("%d    %f  %f\n"%(T, popt[0], popt[0]))
    y_pred = [function_1G(i,*popt) for i in x]
print(popt)
plt.scatter(x,y_pred,marker="s")

plt.savefig(dirname+"/distrib-xydistances.png")
plt.show()
