import dpdata
import numpy as np
from itertools import combinations
from perovskitelattpack import *
import sys

frm=dpdata.System("POSCAR", fmt="vasp/poscar")
cubic = FAPbI3(frm, fmt="vasp/poscar")

elem_type = "I"
indices_elem= cubic.indices_elem
print(indices_elem[elem_type])

num_sub = int(sys.argv[1])
idx = 0
for i in range(10):
    print(i)
    new_cubic = cubic.cubic.copy()
    new_cubic.replace("I", "Br", num_sub)
    new_cubic.to_vasp_poscar("POSCAR_substituted"+str(i)+".vasp")


