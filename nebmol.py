from doctest import FAIL_FAST
from platform import python_implementation
from perovskitelattpack import *
import numpy as np
import dpdata

def genmol_from_twoaxis(normal, polaraxis, center):
    coords = []
    coords.append(center)
    NNcenter = center + polaraxis
    NNvec = np.cross(normal, polaraxis)/np.linalg.norm(np.cross(normal, polaraxis)) * 2.323
    if NNvec[2] < 0:
        NNvec = -NNvec
    N1 = polaraxis + NNvec/2 + center
    dCN1 = polaraxis + NNvec/2
    N2 = polaraxis - NNvec/2 + center
    dCN2 = polaraxis - NNvec/2
    coords.append(N1)
    coords.append(N2)
    H1 = center - polaraxis/np.linalg.norm(polaraxis)*1.09
    H2 = N1 + polaraxis/np.linalg.norm(polaraxis)*1.02
    N1H2 = polaraxis/np.linalg.norm(polaraxis)*1.02
    H3 = N2 + polaraxis/np.linalg.norm(polaraxis)*1.02
    N2H3 = polaraxis/np.linalg.norm(polaraxis)*1.02
    halfHHvec1 = N1H2 - dCN1/np.linalg.norm(dCN1) * 0.526 
    H4 = dCN1/np.linalg.norm(dCN1) * 0.526 - halfHHvec1 + N1
    halfHHvec2 = N2H3 - dCN2/np.linalg.norm(dCN2) * 0.526
    H5 = dCN2/np.linalg.norm(dCN2) * 0.526 - halfHHvec2 + N2
    coords.append(H1)
    coords.append(H2)
    coords.append(H3)
    coords.append(H4)
    coords.append(H5)
    return coords

def add_mol(mid, coords):
    mid["atom_numbs"][3] += 1
    mid.data["atom_types"] = np.append(mid["atom_types"], 3)
    mid.data["coords"] = np.vstack((mid["coords"][0], np.array(coords[1])))[np.newaxis, :]
    mid["atom_numbs"][3] += 1
    mid.data["atom_types"] = np.append(mid["atom_types"], 3)
    mid.data["coords"] = np.vstack((mid["coords"][0], np.array(coords[2])))[np.newaxis, :]
    mid["atom_numbs"][4] += 1
    mid.data["atom_types"] = np.append(mid["atom_types"], 4)
    mid.data["coords"] = np.vstack((mid["coords"][0], np.array(coords[3])))[np.newaxis, :]
    mid["atom_numbs"][4] += 1
    mid.data["atom_types"] = np.append(mid["atom_types"], 4)
    mid.data["coords"] = np.vstack((mid["coords"][0], np.array(coords[4])))[np.newaxis, :]
    mid["atom_numbs"][4] += 1
    mid.data["atom_types"] = np.append(mid["atom_types"], 4)
    mid.data["coords"] = np.vstack((mid["coords"][0], np.array(coords[5])))[np.newaxis, :]
    mid["atom_numbs"][4] += 1
    mid.data["atom_types"] = np.append(mid["atom_types"], 4)
    mid.data["coords"] = np.vstack((mid["coords"][0], np.array(coords[6])))[np.newaxis, :]
    mid["atom_numbs"][4] += 1
    mid.data["atom_types"] = np.append(mid["atom_types"], 4)
    mid.data["coords"] = np.vstack((mid["coords"][0], np.array(coords[7])))[np.newaxis, :]



ini = FAPbI3("ini/CONTCAR", fmt="vasp/poscar")
fin = FAPbI3("fin/CONTCAR", fmt="vasp/poscar")

ini.extract_mol()
fin.extract_mol()

normals_ini = []
polaraxis_ini = []
for mol in ini.molecules:
    mol.set_normal_vec()
    normals_ini.append(mol.normal)
    polaraxis_ini.append(mol.polaraxis)

#print(ini.molecules[0].coords_mol)
#print(genmol_from_twoaxis(ini.molecules[0].normal, ini.molecules[0].polaraxis, ini.molecules[0].coords_mol[0]))

normals_fin = []
polaraxis_fin = []
for mol in fin.molecules:
    mol.set_normal_vec()
    normals_fin.append(mol.normal)
    polaraxis_fin.append(mol.polaraxis)

for imol in range(len(ini.molecules)):
    print(normals_ini[imol], normals_fin[imol])
print("\n")
for imol in range(len(ini.molecules)):
    print(polaraxis_ini[imol], polaraxis_fin[imol])

normals_mid = np.zeros([18,4,3])
polaraxis_mid = np.zeros([18,4,3])
for imid in range(1,17):
    mid = dpdata.System("%02d/POSCAR"%imid, "vasp/poscar")
    mid.add_atom_names("N")
    mid.add_atom_names("H")
    for imol in range(len(ini.molecules)):
        normals_mid[imid][imol] = (normals_fin[imol] - normals_ini[imol])/17 * imid + normals_ini[imol]
        _polaraxis_mid = (polaraxis_fin[imol] - polaraxis_ini[imol])/17 * imid + polaraxis_ini[imol]
        polaraxis_mid[imid][imol] = _polaraxis_mid/np.linalg.norm(_polaraxis_mid) * 0.6144
        center = mid["coords"][0][imol+16]
        coords_mid = genmol_from_twoaxis(normals_mid[imid][imol], polaraxis_mid[imid][imol], center)
        add_mol(mid, coords_mid)

    mid.to_vasp_poscar("%02d.vasp"%(imid))



