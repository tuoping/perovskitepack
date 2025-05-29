import numpy as np
import matplotlib.pyplot as plt

T_K = 300
T_kjmol = T_K * 8.31435 / 1000
kjmol_to_eV = 0.0103643

def bias_potential(x, kappa, x0):
    return (0.5 * kappa * (x - x0) ** 2).sum(axis=-1)

window_idx = np.arange(1, 21)
import os
window_cvs = []
window_biases = []
for i in range(len(window_idx)):
    idx_w = window_idx[i]
    data = np.loadtxt(os.path.join(f"c{idx_w}", "colvar"))[-500:]
    cv = data[:,2:5:2]
    bias = data[:, 5] +data[:, -2]+data[:, -4]+data[:, -9]+data[:, -11]
    window_cvs.append(cv)
    window_biases.append(bias)
window_cvs = np.stack(window_cvs, axis=0)
window_biases = np.stack(window_biases, axis=0)

window_kappas = []
window_centers = []
for i in range(len(window_idx)):
    idx_w = window_idx[i]
    fname = os.path.join(f"c{idx_w}", "input.plumed")
    with open(fname, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "RESTRAINT" in line:
                items = line.split()
                for item in items:
                    if "KAPPA" in item:
                        kappa_str = item
                    if "AT" in item:
                        center_str = item
                c = kappa_str.replace("KAPPA=", "").split(",")
                k = [float(x) for x in c]
                kappa = np.array(k)
                window_kappas.append(kappa)

                c = center_str.replace("AT=", "").split(",")
                center = [float(x) for x in c]
                center = np.array(center)
                window_centers.append(center)
                continue
window_kappas = np.stack(window_kappas, axis=0)
window_centers = np.stack(window_centers, axis=0)

fig, axes = plt.subplots(6, 6, figsize=(32, 20))
for i in range(len(window_idx)):
    if i % 1 == 0:
        ii = (i // 1) % 6
        jj = (i // 1) // 6
        axes[ii,jj].scatter(np.arange(window_cvs[i].shape[0]), window_cvs[i][:,1], s=1)
        axes[ii,jj].axhline(window_centers[i,1], color="red", linestyle="--")
        axes[ii,jj].set_title("window %d"%(i))
        axes[ii,jj].set_xlabel("$Step$")
        axes[ii,jj].set_ylabel("$S_{h}$")
plt.tight_layout()
plt.savefig("window_H.png", dpi=300)

fig, axes = plt.subplots(6, 6, figsize=(32, 20))
for i in range(len(window_idx)):
    if i % 1 == 0:
        ii = (i // 1) % 6
        jj = (i // 1) // 6
        axes[ii,jj].scatter(np.arange(window_cvs[i].shape[0]), window_cvs[i][:,0], s=1)
        axes[ii,jj].axhline(window_centers[i,0], color="red", linestyle="--")
        axes[ii,jj].set_title("window %d"%(i))
        axes[ii,jj].set_xlabel("$Step$")
        axes[ii,jj].set_ylabel("$S_{p}$")
plt.tight_layout()
plt.savefig("window_C.png", dpi=300)

fig, axes = plt.subplots(1, 1, figsize=(4, 3.5))
for i in range(len(window_idx)):
    # if i % 4 == 0:
    #     ii = (i // 4) % 2
    #     jj = (i // 4) // 2
    axes.scatter(window_cvs[i][:,1], window_biases[i], s=5)
    axes.axvline(window_centers[i,1], color="red", linestyle="--", lw=1)
axes.set_title("window %d"%(i+4))
axes.set_xlabel("$S_h$")
axes.set_ylabel("$Bias$")
plt.tight_layout()
plt.savefig("window_biases.png", dpi=300)

n_windows = len(window_idx)
cv_range = np.array([[window_centers[:,0].min(), window_centers[:,0].max()], [window_centers[:,1].min(), window_centers[:,1].max()]])
cv_range_width = np.array([cv_range[0,1] - cv_range[0,0], cv_range[1,1] - cv_range[1,0]])
average_dcv = cv_range_width / (n_windows - 1)
sample_cvs = np.meshgrid(
    np.linspace(cv_range[0,0] -0.5*average_dcv[0] , cv_range[0,1] + 0.5* average_dcv[0], 10),
    np.linspace(cv_range[1,0] -0.5*average_dcv[1] , cv_range[1,1] + 0.5* average_dcv[1], 100)
)
sample_cvs = np.stack(sample_cvs, axis=0)

from sklearn.neighbors import KernelDensity
def logdensity_kde(cvs, centers, weights=None):
    kde_unbiased = KernelDensity(bandwidth=1., kernel='gaussian')
    kde_unbiased.fit(cvs.reshape([-1,2]), sample_weight=weights)
    return kde_unbiased.score_samples(centers)

log_density_biased = logdensity_kde(window_cvs, sample_cvs.reshape([2,-1]).T, )

# $\omaga_i$
window_omega = []
for i in range(len(window_idx)):
    center = window_centers[i,:]
    omega = bias_potential(sample_cvs.reshape([2,-1]).T, window_kappas[i].reshape([1,2]), center.reshape([1,2]))
    window_omega.append(omega)
window_omega = np.stack(window_omega, axis=0)

# $logP_i^b$
window_logP_biased = []
for i in range(len(window_idx)):
    logP_biased = logdensity_kde(window_cvs[i], sample_cvs.reshape([2,-1]).T)
    window_logP_biased.append(logP_biased)
window_logP_biased = np.stack(window_logP_biased, axis=0)

# $logP_i^u$
window_logP_unbiased = []
for i in range(len(window_idx)):
    logP_unbiased = logdensity_kde(window_cvs[i], sample_cvs.reshape([2,-1]).T, weights=np.exp(window_biases[i].flatten()/T_kjmol))
    window_logP_unbiased.append(logP_unbiased)
window_logP_unbiased = np.stack(window_logP_unbiased, axis=0)

def iterate_wham(log_density_unbiased_old, i_iter = 0):
    window_F = []
    for i in range(len(window_idx)):
        if i == 0:
            w_lb = window_centers[i,:] - average_dcv / 2
            w_ub = (window_centers[i,:] + window_centers[i+1,:]) / 2
        elif i == len(window_idx) - 1:
            w_lb = (window_centers[i-1,:] + window_centers[i,:]) / 2
            w_ub = window_centers[i,:] + average_dcv / 2
        else:
            w_lb = (window_centers[i-1,:] + window_centers[i,:]) / 2
            w_ub = (window_centers[i,:] + window_centers[i+1,:]) / 2
        w_dcv = np.abs((w_ub[0] - w_lb[0])*(w_ub[1] - w_lb[1]))
        F_i = -np.log(np.sum(np.exp(log_density_unbiased_old)*np.exp(-window_omega[i]/T_kjmol)*w_dcv))*T_kjmol
        window_F.append(F_i)
    window_F = np.stack(window_F, axis=0)

    window_a = []
    for i in range(len(window_idx)):
        N_i = window_cvs[i].shape[0]
        bias = bias_potential(sample_cvs.reshape([2,-1]).T, window_kappas[i].reshape([1,2]), window_centers[i,:].reshape([1,2]))
        a_i = N_i * np.exp(-bias/T_kjmol+window_F[i]/T_kjmol)
        window_a.append(a_i)
    window_a = np.stack(window_a, axis=0)
    window_p = window_a/np.sum(window_a, axis=0)

    log_density_unbiased_new = np.log(np.sum( window_p*np.exp(window_logP_unbiased), axis=0))

    fig, axes = plt.subplots(1,3, figsize=(14, 3.5))

    c0 = axes[0].contourf(sample_cvs[0], sample_cvs[1], -log_density_biased.reshape(sample_cvs[0].shape), levels=20, cmap="plasma")
    cbar0 = plt.colorbar(c0, ax=axes[0])
    cbar0.set_label("$-lnP^b$ [$k_BT$]")
    # axes[0].set_title("$lnP^u$")
    axes[0].set_xlabel("$S_{c}$")
    axes[0].set_ylabel("$S_{h}$")

    c1 = axes[1].contourf(sample_cvs[0], sample_cvs[1], -log_density_unbiased_old.reshape(sample_cvs[0].shape), levels=20, cmap="plasma")
    cbar1 = plt.colorbar(c1, ax=axes[1])
    cbar1.set_label("$-lnP^u_1$ [$k_BT$]")
    axes[1].set_xlabel("$S_{c}$")
    axes[1].set_ylabel("$S_{h}$")
    # axes[1].set_title("$lnP^b$")

    c2 = axes[2].contourf(sample_cvs[0], sample_cvs[1], -log_density_unbiased_new.reshape(sample_cvs[0].shape), levels=20, cmap="plasma")
    cbar2 = plt.colorbar(c2, ax=axes[2])
    cbar2.set_label("$-lnP^b_2$ [$k_BT$]")
    axes[2].set_xlabel("$S_{c}$")
    axes[2].set_ylabel("$S_{h}$")

    plt.tight_layout()
    plt.savefig("wham_%d.png"%(i_iter), dpi=300)

    return log_density_unbiased_new, window_F, window_a, window_p

log_density_unbiased = np.zeros_like(log_density_biased)
for i in range(5):
    log_density_unbiased_new, _, _, _ = iterate_wham(log_density_unbiased, i_iter=i)
    log_density_unbiased = log_density_unbiased_new
    print("iter %d done"%(i+1))
