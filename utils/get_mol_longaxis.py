from re import L
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation

dirname = "./"
caxis = int(np.loadtxt(dirname+"/caxis"))
ref_axis = np.eye(3)
axis = np.eye(3)
print(axis)

f1 = "mol_center_coords"

coord = np.loadtxt(f1, skiprows=0)
coord = coord.T

caxis = int(np.loadtxt("caxis"))
x = np.array(coord[(caxis+1)%3])
y = np.array(coord[(caxis+2)%3])
z = np.array(coord[caxis])
# phi_list_all = []
for i_frame in range(1):

    Totaltime = 100
    Num_mol = 1800
    _data = np.loadtxt(dirname+"/mol_longaxis_frame"+str(i_frame), skiprows=0)
    data = _data
    data = np.reshape(data, [-1, Num_mol, 3])
    Totaltime = len(data)
    print("Totaltime = ", Totaltime)
    mesh_starting_point = 0
    idx = 0
    while(True):
        if np.argmax(np.abs(data[0][idx])) != caxis:
            if np.argmax(np.abs(data[0][idx])) == (caxis+1)%3:
                mesh_starting_point = 0
            if np.argmax(np.abs(data[0][idx])) == (caxis+2)%3:
                mesh_starting_point = 1
            break
        idx += 2
        if idx >= Num_mol:
            raise Exception("Long axis not in xOy")
    print("mesh_starting_point = ", mesh_starting_point, "idx = ", idx, "argmax = ", np.argmax(np.abs(data[0][0])), "caxis = ", caxis)    

    fphi = open("phi_list_frame"+str(i_frame)+".out", "w")
    phix_list = []
    phiy_list = []
    phi_list = []
    theta_list = []
    theta_list_positive = []
    theta_list_negative = []
    # for m in range(1):
    for m in range(Totaltime):
        idx = 0
        data_xy = [[],[]]
        for i in range(30):
            for j in range(30):
                for k in range(2):
                    # if data[m][idx][np.argmax(np.abs(data[m][idx]))] <0:
                    #     data[m][idx] = -data[m][idx]
                    mol = np.array([data[m][idx][(caxis+1)%3],data[m][idx][(caxis+2)%3],data[m][idx][caxis]])
                    mol_z = np.dot(mol, axis[2])
                    #theta = math.degrees(math.acos(np.dot(mol, axis[2])))
                    mol_x = np.dot(mol, axis[0])
                    mol_y = np.dot(mol, axis[1])
                    mol_xOy = np.sqrt(mol_x*mol_x + mol_y * mol_y)
                    theta = math.degrees(math.atan(mol_z/mol_xOy))
                    if mol_x == 0:
                        phi = 90
                    else:
                        phi = math.degrees( math.atan(mol_y/mol_x))
                    if abs(math.atan(mol_y/mol_x)) < abs(math.atan(mol_x/mol_y)):
                        phi_list.append(math.degrees(math.atan(mol_y/mol_x)))
                    else:
                        phi_list.append(math.degrees(math.atan(mol_x/mol_y)))
                    theta_list.append(theta)
                    if theta > 0:
                        theta_list_positive.append(theta)
                    else:
                        theta_list_negative.append(theta)    

                    if caxis == 2:
                        if  (i+j)%2 == mesh_starting_point:
                            phix_list.append(phi)
                        else:
                            phiy_list.append(math.degrees( math.atan(mol_x/mol_y)))
                        if k== 2:
                            data_xy[0].append(mol[0])
                            data_xy[1].append(mol[1])
                    if caxis == 1:
                        if  (k+i)%2 == mesh_starting_point:
                            phix_list.append(phi)
                        else:
                            phiy_list.append(math.degrees( math.atan(mol_x/mol_y)))
                        if j== 2:
                            data_xy[0].append(mol[0])
                            data_xy[1].append(mol[1])
                    if caxis == 0:
                        if  (j+k)%2 == mesh_starting_point:
                            phix_list.append(phi)
                        else:
                            phiy_list.append(math.degrees( math.atan(mol_x/mol_y)))
                        if i== 2:
                            data_xy[0].append(mol[0])
                            data_xy[1].append(mol[1])
                    idx += 1
                    fphi.write("%d  %d  %d   %f\n" %(i,j,k,phi_list[-1]))
        # raise Exception("PAUSE")
    fphi.close()
    phix_list = list(phix_list)
    phiy_list = list(phiy_list)
    # phi_list = phix_list+phiy_list
    np.save("phi_list_frame"+str(i_frame)+".npy", phi_list)    
    # phi_list_all.append(phi_list)
