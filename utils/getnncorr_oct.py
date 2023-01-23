from re import L
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def frame_order_parameter_by_mesh(mesh_dim,mesh, octlist):
    mesh = mesh.T

    corr_x = 1.0
    corr_y = 1.0
    corr_z = 1.0
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
        
        assert np.dot(octlist[idx_mol][0], octlist[xnlist[0]][0])/np.linalg.norm(octlist[idx_mol][0])/np.linalg.norm(octlist[xnlist[0]][0])-1.0 < 1e-6, print(np.dot(octlist[idx_mol][0], octlist[xnlist[0]][0])-1.0)
        assert(np.dot(octlist[idx_mol][1], octlist[xnlist[0]][1])/np.linalg.norm(octlist[idx_mol][1])/np.linalg.norm(octlist[xnlist[0]][1])-1.0 < 1e-6)
        assert(np.dot(octlist[idx_mol][2], octlist[xnlist[0]][2])/np.linalg.norm(octlist[idx_mol][2])/np.linalg.norm(octlist[xnlist[0]][2])-1.0 < 1e-6)
        assert(np.dot(octlist[idx_mol][0], octlist[xnlist[1]][0])/np.linalg.norm(octlist[idx_mol][0])/np.linalg.norm(octlist[xnlist[1]][0])-1.0 < 1e-6)
        assert(np.dot(octlist[idx_mol][1], octlist[xnlist[1]][1])/np.linalg.norm(octlist[idx_mol][1])/np.linalg.norm(octlist[xnlist[1]][1])-1.0 < 1e-6)
        assert(np.dot(octlist[idx_mol][2], octlist[xnlist[1]][2])/np.linalg.norm(octlist[idx_mol][2])/np.linalg.norm(octlist[xnlist[1]][2])-1.0 < 1e-6)

        corr_x *= np.dot(octlist[idx_mol][1], octlist[xnlist[0]][1])/np.linalg.norm(octlist[idx_mol][1])/np.linalg.norm(octlist[xnlist[0]][1])
        corr_x *= np.dot(octlist[idx_mol][2], octlist[xnlist[0]][2])/np.linalg.norm(octlist[idx_mol][2])/np.linalg.norm(octlist[xnlist[0]][2])
        corr_x *= np.dot(octlist[idx_mol][1], octlist[xnlist[1]][1])/np.linalg.norm(octlist[idx_mol][1])/np.linalg.norm(octlist[xnlist[1]][1])
        corr_x *= np.dot(octlist[idx_mol][2], octlist[xnlist[1]][2])/np.linalg.norm(octlist[idx_mol][2])/np.linalg.norm(octlist[xnlist[1]][2])
        corr_y *= np.dot(octlist[idx_mol][0], octlist[ynlist[0]][0])/np.linalg.norm(octlist[idx_mol][0])/np.linalg.norm(octlist[xnlist[0]][0])
        corr_y *= np.dot(octlist[idx_mol][2], octlist[ynlist[0]][2])/np.linalg.norm(octlist[idx_mol][2])/np.linalg.norm(octlist[xnlist[0]][2])
        corr_y *= np.dot(octlist[idx_mol][0], octlist[ynlist[1]][0])/np.linalg.norm(octlist[idx_mol][0])/np.linalg.norm(octlist[xnlist[1]][0])
        corr_y *= np.dot(octlist[idx_mol][2], octlist[ynlist[1]][2])/np.linalg.norm(octlist[idx_mol][2])/np.linalg.norm(octlist[xnlist[1]][2])
        corr_z *= np.dot(octlist[idx_mol][0], octlist[znlist[0]][0])/np.linalg.norm(octlist[idx_mol][0])/np.linalg.norm(octlist[xnlist[0]][0])
        corr_z *= np.dot(octlist[idx_mol][1], octlist[znlist[0]][1])/np.linalg.norm(octlist[idx_mol][1])/np.linalg.norm(octlist[xnlist[0]][1])
        corr_z *= np.dot(octlist[idx_mol][0], octlist[znlist[1]][0])/np.linalg.norm(octlist[idx_mol][0])/np.linalg.norm(octlist[xnlist[1]][0])
        corr_z *= np.dot(octlist[idx_mol][1], octlist[znlist[1]][1])/np.linalg.norm(octlist[idx_mol][1])/np.linalg.norm(octlist[xnlist[1]][1])
          
    return math.pow(corr_x, 2/len(octlist)), math.pow(corr_y, 2/len(octlist)), math.pow(corr_z, 2/len(octlist))

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
    # dirname = sys.argv[1]
    dirname = "./"
    caxis = int(np.loadtxt(dirname+"/caxis"))
    ref_axis = np.eye(3)
    axis = np.eye(3)
    #print(axis)
    
    ofcorr = open("octcorr-t.out", "w")
    mesh_dim = [8,8,4]
    for i_frame in range(3000000, 4010000, 10000):
    
        #Totaltime = 100
        Num_mol = mesh_dim[0]*mesh_dim[1]*mesh_dim[2]
    
        data0 = np.loadtxt(dirname+"ii_vectors0_mesh_frame"+str(i_frame)+".dat", skiprows=0)
        data1 = np.loadtxt(dirname+"ii_vectors1_mesh_frame"+str(i_frame)+".dat", skiprows=0)
        data2 = np.loadtxt(dirname+"ii_vectors2_mesh_frame"+str(i_frame)+".dat", skiprows=0)

        data = [[data0[i], data1[i], data2[i]] for i in range(len(data0))]
    
        data = np.reshape(data, [-1, Num_mol, 3, 3])
        assert len(data) == 1, print(len(data))
        #Totaltime = len(data)
        #print("Totaltime = ", Totaltime)
    
        rot = np.zeros([3,3])
        rot[0][(caxis+1)%3] = 1
        rot[1][(caxis+2)%3] = 1
        rot[2][(caxis)%3] = 1
        # for m in range(Totaltime):
        mesh = np.zeros([Num_mol,3])
        idx = 0
        for i in range(mesh_dim[0]):
            for j in range(mesh_dim[1]):
                for k in range(mesh_dim[2]):
                    # data[0][idx] = rot@data[0][idx]
                    # data[0][idx][0] = rot@data[0][idx][0]
                    # data[0][idx][1] = rot@data[0][idx][1]
                    # data[0][idx][2] = rot@data[0][idx][2]
                    mesh[idx][0] = i
                    mesh[idx][1] = j
                    mesh[idx][2] = k
                    idx += 1
        corr = frame_order_parameter_by_mesh(mesh_dim, mesh, data[0])
        ofcorr.write("%d   %f  %f  %f\n"%(i_frame, corr[0], corr[1], corr[2]))
