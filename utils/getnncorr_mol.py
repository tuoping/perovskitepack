from re import L
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def molecular_order_parameter_by_mesh(mesh_dim, mesh, mollist, polarized=False):
    mesh = mesh.T
    vecdot_x = 0.0
    vecdot_y = 0.0
    vecdot_z = 0.0
    corr_x = 0.0
    corr_y = 0.0
    corr_z = 0.0
    for idx_mol in range(mesh_dim[0]*mesh_dim[1]*mesh_dim[2]):
        x_mol = mesh[0][idx_mol]
        y_mol = mesh[1][idx_mol]
        z_mol = mesh[2][idx_mol]
        
        if x_mol - 1 < 0:
            x_mol_left = max(mesh[0])
        else:
            x_mol_left = x_mol - 1
        if y_mol - 1 < 0:
            y_mol_left = max(mesh[1])
        else:
            y_mol_left = y_mol - 1
        if z_mol - 1 < 0:
            z_mol_left = max(mesh[2])
        else:
            z_mol_left = z_mol - 1
        
        if x_mol + 1 > max(mesh[0]):
            x_mol_right = 0
        else:
            x_mol_right = x_mol + 1
        if y_mol + 1 > max(mesh[1]):
            y_mol_right = 0
        else:
            y_mol_right = y_mol + 1
        if z_mol + 1 > max(mesh[2]):
            z_mol_right = 0
        else:
            z_mol_right = z_mol + 1
        
        xlimit = np.where(mesh[0]==x_mol)[0]
        ylimit = np.where(mesh[1]==y_mol)[0]
        zlimit = np.where(mesh[2]==z_mol)[0]
        
        xplane = np.concatenate((np.array(np.where(mesh[0]==x_mol_left)[0]), np.array(np.where(mesh[0]==x_mol_right)[0])))
        xnlist = np.intersect1d(np.intersect1d(xplane, ylimit), zlimit)
        
        yplane = np.concatenate((np.array(np.where(mesh[1]==y_mol_left)[0]), np.array(np.where(mesh[1]==y_mol_right)[0])))
        ynlist = np.intersect1d(np.intersect1d(yplane, xlimit), zlimit)
        
        zplane = np.concatenate((np.array(np.where(mesh[2]==z_mol_left)[0]), np.array(np.where(mesh[2]==z_mol_right)[0])))
        znlist = np.intersect1d(np.intersect1d(zplane, xlimit), ylimit)
        
        assert(len(xnlist)==2)
        assert(len(ynlist)==2)
        assert(len(znlist)==2)
        
        assert(np.abs(np.dot(mollist[idx_mol], mollist[xnlist[0]]))-1.0 < 1e-6)
        assert(np.abs(np.dot(mollist[idx_mol], mollist[ynlist[0]]))-1.0 < 1e-6)
        assert(np.abs(np.dot(mollist[idx_mol], mollist[znlist[0]]))-1.0 < 1e-6)
        assert(np.abs(np.dot(mollist[idx_mol], mollist[xnlist[1]]))-1.0 < 1e-6)
        assert(np.abs(np.dot(mollist[idx_mol], mollist[ynlist[1]]))-1.0 < 1e-6)
        assert(np.abs(np.dot(mollist[idx_mol], mollist[znlist[1]]))-1.0 < 1e-6)
        if polarized:
            corr_x += math.acos(min((np.dot(mollist[idx_mol], mollist[xnlist[0]])),1.0))
            corr_y += math.acos(min((np.dot(mollist[idx_mol], mollist[ynlist[0]])),1.0))
            corr_z += math.acos(min((np.dot(mollist[idx_mol], mollist[znlist[0]])),1.0))
            corr_x += math.acos(min((np.dot(mollist[idx_mol], mollist[xnlist[1]])),1.0))
            corr_y += math.acos(min((np.dot(mollist[idx_mol], mollist[ynlist[1]])),1.0))
            corr_z += math.acos(min((np.dot(mollist[idx_mol], mollist[znlist[1]])),1.0))
            vecdot_x += (min((np.dot(mollist[idx_mol], mollist[xnlist[0]])),1.0))
            vecdot_y += (min((np.dot(mollist[idx_mol], mollist[ynlist[0]])),1.0))
            vecdot_z += (min((np.dot(mollist[idx_mol], mollist[znlist[0]])),1.0))
            vecdot_x += (min((np.dot(mollist[idx_mol], mollist[xnlist[1]])),1.0))
            vecdot_y += (min((np.dot(mollist[idx_mol], mollist[ynlist[1]])),1.0))
            vecdot_z += (min((np.dot(mollist[idx_mol], mollist[znlist[1]])),1.0))
        else:
            corr_x += math.acos(min(np.abs(np.dot(mollist[idx_mol], mollist[xnlist[0]])),1.0))
            corr_y += math.acos(min(np.abs(np.dot(mollist[idx_mol], mollist[ynlist[0]])),1.0))
            corr_z += math.acos(min(np.abs(np.dot(mollist[idx_mol], mollist[znlist[0]])),1.0))
            corr_x += math.acos(min(np.abs(np.dot(mollist[idx_mol], mollist[xnlist[1]])),1.0))
            corr_y += math.acos(min(np.abs(np.dot(mollist[idx_mol], mollist[ynlist[1]])),1.0))
            corr_z += math.acos(min(np.abs(np.dot(mollist[idx_mol], mollist[znlist[1]])),1.0))
            vecdot_x += (min(np.abs(np.dot(mollist[idx_mol], mollist[xnlist[0]])),1.0))
            vecdot_y += (min(np.abs(np.dot(mollist[idx_mol], mollist[ynlist[0]])),1.0))
            vecdot_z += (min(np.abs(np.dot(mollist[idx_mol], mollist[znlist[0]])),1.0))
            vecdot_x += (min(np.abs(np.dot(mollist[idx_mol], mollist[xnlist[1]])),1.0))
            vecdot_y += (min(np.abs(np.dot(mollist[idx_mol], mollist[ynlist[1]])),1.0))
            vecdot_z += (min(np.abs(np.dot(mollist[idx_mol], mollist[znlist[1]])),1.0))
          
    # print(vecdot_x, vecdot_y, vecdot_z)
    return corr_x/2/math.pi/len(mollist), corr_y/2/math.pi/len(mollist), corr_z/2/math.pi/len(mollist)


if __name__ == "__main__":
    dirname = "./"
    caxis = int(np.loadtxt(dirname+"/caxis"))
    ref_axis = np.eye(3)
    axis = np.eye(3)
    #print(axis)
    
    caxis = int(np.loadtxt("caxis"))
    mesh_dim = [8,8,4]
    for i_frame in range(1):
    
        Totaltime = 100
        Num_mol = mesh_dim[0]*mesh_dim[1]*mesh_dim[2]
    
        _data = np.loadtxt(dirname+"/mol_longaxis_frame"+str(i_frame), skiprows=0)
        data = _data
    
        data = np.reshape(data, [-1, Num_mol, 3])
        Totaltime = len(data)
        #print("Totaltime = ", Totaltime)
    
        rot = np.zeros([3,3])
        rot[0][(caxis+1)%3] = 1
        rot[1][(caxis+2)%3] = 1
        rot[2][(caxis)%3] = 1
        for m in range(Totaltime):
            mesh = np.zeros([Num_mol,3])
            idx = 0
            for i in range(8):
                for j in range(8):
                    for k in range(4):
                        data[m][idx] = rot@data[m][idx]
                        mesh[idx][0] = i
                        mesh[idx][1] = j
                        mesh[idx][2] = k
                        idx += 1
            corr = molecular_order_parameter_by_mesh(mesh_dim, mesh, data[m])
            print(corr)
