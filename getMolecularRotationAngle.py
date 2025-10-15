import numpy as np
import glob

# num_steps = 1000+1
num_steps = len(glob.glob("unitaxis_mols/molecule_unitlongaxis_frame*.dat"))
print("num of steps = ", num_steps)
max_lag = 5000
eq_step = 6000
filenames = ["unitaxis_mols/molecule_unitlongaxis_frame%d.dat"%ii for ii in range(num_steps)]
# start_vectors = np.loadtxt(filenames[0])
all_angles = []
for ii in range(1, num_steps):
    vectors = np.loadtxt(filenames[ii])
    start_vectors = np.zeros(vectors.shape)
    start_vectors[:,-1] = 1.
    angles = np.arccos(np.einsum('ij,ij->i', start_vectors, vectors))/np.pi*180
    all_angles.append(angles)
all_angles = np.array(all_angles)

import matplotlib.pyplot as plt
plt.figure()
for i in range(0, len(start_vectors), 2):
    plt.plot(np.arange(1, num_steps), all_angles[:,i])
plt.xlabel("Time (*0.5 fs)")
plt.ylabel("Angle")
plt.savefig("Angle-time.png", bbox_inches="tight")

all_angles = all_angles[eq_step:]

ave_angle = np.mean(all_angles[max_lag:], axis=0)
dev_angle = all_angles[max_lag:]-ave_angle
autocorr_angle = np.zeros([max_lag, len(start_vectors)])
for i in range(max_lag):
    for k in range(len(dev_angle)-max_lag):
        autocorr_angle[i] += dev_angle[k]*dev_angle[k+i]
autocorr_angle /= max_lag
autocorr_angle /= np.std(dev_angle)
autocorr_angle /= np.std(dev_angle)

def large_lag_std_err(r):
    var_r = [(1+2*np.sum(r[:k]*r[:k]))/len(r) for k in range(len(r))]
    return np.sqrt(var_r)

var_autocorr = large_lag_std_err(autocorr_angle[:,0])

plt.figure()
for i in range(0, len(start_vectors), 2):
    plt.plot(autocorr_angle[:,i])

plt.plot(var_autocorr, color='k', linestyle='--')
plt.plot(-var_autocorr, color='k', linestyle='--')
plt.axhline(0, c="k")

plt.xlabel("$\Delta t$ (*0.5 fs)")
plt.ylabel(r"<$\delta \theta (t+\Delta t)\delta \theta (t)$>")
plt.savefig("autocorr_Angle.png", bbox_inches = "tight")
