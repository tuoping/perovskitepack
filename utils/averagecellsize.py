import numpy as np

def pbc_idx(idx):
    if idx > 2:
        return 0
    if idx < 0:
        return 2
    return idx

of = open("lattab_pbpbdistances.out", "w")
for T in [i for i in range(20, 155, 5)]+[i for i in range(155, 440,5)]:
    dirname = None
    if T < 155:
        dirname = "annealT"+str(T)
    else:
        dirname = "T"+str(T)
    print(dirname)

    d = np.loadtxt(dirname+"/cellsize.out", skiprows=1)
    d = d.T
    caxis = int(np.loadtxt(dirname+"/caxis"))
    alist = []
    blist = []
    clist = []
    for i in range(len(d[1])-500,len(d[1])):
        a = d[pbc_idx(caxis+1)][i]
        b = d[pbc_idx(caxis-1)][i]
        c = d[caxis][i]
        alist.append(a)
        blist.append(b)
        clist.append(c)
        '''
        max_a = max(max_a, a)
        min_a = min(min_a, a)
        max_b = max(max_b, b)
        min_b = min(min_b, b)
        max_c = max(max_c, c)
        min_c = min(min_c, c)
        aa += a
        bb += b
        cc += c
        '''
    
    # aa /= float(500)
    # bb /= float(500)
    # cc /= float(500)
    aa = np.average(np.array(alist))
    bb = np.average(np.array(blist))
    cc = np.average(np.array(clist))
    stda = np.std(np.array(alist))
    stdb = np.std(np.array(blist))
    stdc = np.std(np.array(clist))
    of.write("%5d    %11.5f  %11.5f  %11.5f  %11.5f\n"%(T, bb, cc, stdb, stdc))
