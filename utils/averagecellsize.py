import numpy as np

def pbc_idx(idx):
    if idx > 2:
        return 0
    if idx < 0:
        return 2
    return idx

of = open("lattab.out", "w")
of2 = open("diag-lattab.out", "w")
for T in [i for i in range(20, 155, 5)]+[i for i in range(155, 440,5)]:
    dirname = None
    if T < 155:
        dirname = "T"+str(T)
    else:
        dirname = "T"+str(T)
    print(dirname)

    d = np.loadtxt(dirname+"/cellsize.out", skiprows=1)
    d2 = np.loadtxt(dirname+"/diag-cellsize.out", skiprows=1)
    d = d.T
    d2 = d2.T
    caxis = int(np.loadtxt(dirname+"/caxis"))
    # caxis = 2
    alist = []
    blist = []
    clist = []
    a2list = []
    b2list = []
    c2list = []
    for i in range(len(d[1])-500,len(d[1])):
        a = d[pbc_idx(caxis+1)+1][i]/2
        b = d[pbc_idx(caxis-1)+1][i]/2
        c = d[caxis+1][i]/2
        alist.append(a)
        blist.append(b)
        clist.append(c)
        a2 = d2[1][i]/2/np.sqrt(2)
        b2 = d2[2][i]/2/np.sqrt(2)
        c2 = d2[3][i]/2
        a2list.append(a2)
        b2list.append(b2)
        c2list.append(c2)
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
    aa2 = np.average(np.array(a2list))
    bb2 = np.average(np.array(b2list))
    cc2 = np.average(np.array(c2list))
    stda = np.std(np.array(alist))
    stdb = np.std(np.array(blist))
    stdc = np.std(np.array(clist))
    of.write("%5d    %11.5f  %11.5f  %11.5f  %11.5f  %11.5f  %11.5f\n"%(T, aa, bb, cc, stda, stdb, stdc))
    of2.write("%5d    %11.5f  %11.5f  %11.5f\n"%(T, aa2, bb2, cc2))
