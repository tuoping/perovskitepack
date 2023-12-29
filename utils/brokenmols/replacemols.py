from copy import deepcopy
def euler_rotate(center, angle, point, mode="z"):
    """
    Rotate a point around a given center point.
    """
    angle = np.deg2rad(angle)
    temp_point = deepcopy(point)
    temp_point -= center
    if mode == "z":
        temp_point = np.dot(temp_point, np.array([[np.cos(angle), -np.sin(angle),0.], [np.sin(angle), np.cos(angle),0.], [0.,0.,1.]]))
    elif mode == "x":
        temp_point = np.dot(temp_point, np.array([[1.,0.,0.], [0.,np.cos(angle), -np.sin(angle)], [0.,np.sin(angle), np.cos(angle)]]))
    else:
        temp_point = np.dot(temp_point, np.array([[np.cos(angle), 0., np.sin(angle)], [0.,1.,0.], [-np.sin(angle), 0., np.cos(angle)]]))
    temp_point += center
    return temp_point

import random
import numpy as np
def calcn(coord, pool, type_pool, cell, idx_mol):
    """
    Calculate the coordination number of a given center point.
    """
    cutoffs = {
        "I-C": 2.5,
        "C-I": 2.5,
        "I-N": 2.5,
        "N-I": 2.5,
        "I-H": 1.5,
        "H-I": 1.5,
        "Pb-C": 2.3,
        "C-Pb": 2.3,
        "Pb-N": 2.3,
        "N-Pb": 2.3,
        "Pb-H": 1.3,
        "H-Pb": 1.3,
        "C-N": 2.0,
        "N-C": 2.0,
        "C-H": 1.3,
        "H-C": 1.3,
        "N-H": 1.3,
        "H-N": 1.3,
        "C-C": 2.0,
        "N-N": 2.5,
        "H-H": 1.1,
    }
    atom_names = ["I", "Pb", "C", "N", "H"]
    type_coord = ["C", "N", "N", "H", "H", "H", "H", "H"]
    cn = 0
    d = get_distances(coord, pool, cell, pbc=[1,1,1])[1]
    for i in range(d.shape[0]):
        for j in range(d.shape[1]):
            if j in idx_mol:
                # print("skip", i, j, d[i][j])
                continue
            connection_type = type_coord[i]+"-"+atom_names[type_pool[j]]
            if d[i][j] < cutoffs[connection_type]:
                cn += 1
    return cn

rot_angles = {
}

def checkcn(vec, idx_mol, pool=[], type_pool=[], cell=np.eye(3)):
    new_vec = vec
    idx_center = idx_mol[0]
    center0 = pool[idx_mol[0]]
    center = deepcopy(pool[idx_mol[0]])
    displace_center = [0.,0.,0.]
    num_try = 0
    num_displace = 0
    rangle1, rangle2, rangle3 = 0.0, 0.0, 0.0
    while True:
        new_coord = center + new_vec
        cn = calcn(new_coord, pool, type_pool, cell, idx_mol)
        print("cn = ", cn)
        if cn < 1:
            break
        elif num_try > 100:
            displace_center = np.random.randn(3)*0.5
            center = center0 + displace_center
            print("displace center by", displace_center, center)
            num_try = 0
            num_displace += 1
        else:
            rangle1 = np.random.rand()*180
            new_vec = euler_rotate([0., 0., 0.], rangle1, vec, mode="x")
            rangle2 = np.random.rand()*180
            new_vec = euler_rotate([0., 0., 0.], rangle2, new_vec, mode="y")
            rangle3 = np.random.rand()*180
            new_vec = euler_rotate([0., 0., 0.], rangle3, new_vec, mode="z")
            print(rangle1, rangle2, rangle3, num_try)
            num_try += 1
        
    return new_vec, center

import dpdata
conf0 = dpdata.System("../conf.lmp", "lammps/lmp")
conf1 = dpdata.System("end.dump", "lammps/dump")
data_new = deepcopy(conf1.data)
data_mols = {}
data_mols["cells"] = conf1["cells"]
data_mols["orig"] = conf1["orig"]
data_mols["atom_names"] = conf1["atom_names"]
data_mols["atom_types"] = []
data_mols["coords"] = []

from ase.geometry.geometry import get_distances

lines_toreplace = open("log.brokenmols","r").readlines()
idx_toreplace = []
for line in lines_toreplace:
    if "nohup" in line:
        continue
    idx_toreplace.append(int(line.split()[0])-1)
    print(idx_toreplace[-1])
step = 1
mols = np.loadtxt("mols.dat", dtype=int)
for idx in idx_toreplace:
    idx_C = mols[idx][0]
    print(idx_C, conf1["atom_types"][idx_C], conf1["coords"][0][idx_C])

    vec = get_distances(conf0["coords"][0][mols[idx]], conf0["coords"][0][idx_C], conf0["cells"][0], pbc=[1,1,1])[0]
    vec = np.reshape(vec, [8,3])
    vec,conf1["coords"][0][idx_C] = checkcn(vec, mols[idx], conf1["coords"][0], conf1["atom_types"], conf1["cells"][0])
    print("center = ", conf1["coords"][0][idx_C], " vec = ", vec)
    conf1["coords"][0][mols[idx]] = conf1["coords"][0][idx_C] + vec

    data_new["coords"][0][idx_C] = conf1["coords"][0][idx_C]
    data_new["coords"][0][mols[idx]] = data_new["coords"][0][idx_C] + vec
    
    for idx_mol in mols[idx]:
        data_mols["atom_types"].append( conf1["atom_types"][idx_mol])
    data_mols["coords"].append(data_new["coords"][0][mols[idx]])

    conf1.to_lammps_lmp(f"step{step}.lmp")
    step += 1

data_mols["atom_types"] = np.array(data_mols["atom_types"])
data_mols["coords"] = np.array(data_mols["coords"]).reshape([1,-1,3])
natoms_vec = []
natomtypes = len(data_mols["atom_names"])
for ii in range(natomtypes):
    natoms_vec.append(sum(np.array(data_mols["atom_types"])==ii+1))
data_mols["atom_numbs"] = natoms_vec

conf1_new = dpdata.System(data=data_new)
conf1_new.to_lammps_lmp("end_new.lmp")
conf1.to_lammps_lmp("end_old.lmp")
conf1_new_mols = dpdata.System(data=data_mols)
conf1_new_mols.to_lammps_lmp("end_new_mols.lmp")
