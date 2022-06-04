import numpy as np
cs = np.loadtxt("cellsize.out", skiprows=1)

print(cs[-1])
idx = np.argmax(cs[-1][1:])
fout = open("caxis", "w")
fout.write(str(idx)+"\n")
fout.close()
