import glob
import re

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = "Data/Random_imagtime"
BETA_ARRAY = [0.5, 1.5, 2.5]


def find_files():
    """Zero-field (symm_type A) imaginary-time files only; skip the h_z runs."""
    files = {}
    for path in glob.glob(f"{DATA_DIR}/ISO__Random__N=*__beta=*__numConfigs=*.hdf5"):
        if "h_z" in path:
            continue
        m = re.search(r"N=(\d+)__beta=([\d.]+)", path)
        N, beta = int(m.group(1)), float(m.group(2))
        files[(N, beta)] = path
    return files


def load(path):
    f = h5.File(path, "r")
    params = f["parameters"].attrs
    tau = np.linspace(0.0, params["Tmax"], params["num_TimePoints"])  # Tmax = beta/2
    re = f["results"]["Re_correlation"][0]
    re_err = f["results"]["Re_stddev"][0]
    return tau, re, re_err


def extrapolate_1_over_N(Ns, values, errors):
    """Weighted least-squares fit of y = C + a/N at each timepoint, returns C(tau), C_err(tau)."""
    x = np.array([1.0 / N for N in Ns])
    values = np.array(values)
    errors = np.where(np.array(errors) == 0, 1e-12, np.array(errors))
    w = 1.0 / errors**2
    S = np.sum(w, axis=0)
    Sx = np.sum(w * x[:, None], axis=0)
    Sxx = np.sum(w * (x**2)[:, None], axis=0)
    Sy = np.sum(w * values, axis=0)
    Sxy = np.sum(w * x[:, None] * values, axis=0)
    Delta = S * Sxx - Sx**2
    C = (Sxx * Sy - Sx * Sxy) / Delta
    C_err = np.sqrt(Sxx / Delta)
    return C, C_err


files = find_files()
N_array = sorted({N for N, _ in files})

color_cycle = [c for c in plt.rcParams["axes.prop_cycle"].by_key()["color"] if c != "red"]
colors = {N: color_cycle[i % len(color_cycle)] for i, N in enumerate(N_array)}

ylabel = r"$g^{xx}(\tau)$"
outpath = "Plots/Plot_imagtime_random.pdf"

fig, axes = plt.subplots(1, len(BETA_ARRAY), figsize=(5 * len(BETA_ARRAY), 4.5), sharey=True)

for ax, beta in zip(axes, BETA_ARRAY):
    tau_common = None
    fit_values, fit_errors, fit_Ns = [], [], []
    for N in N_array:
        key = (N, beta)
        if key not in files:
            continue
        tau, re, re_err = load(files[key])
        tau = tau / beta
        tau_common = tau
        ax.errorbar(
            tau, re, yerr=re_err, errorevery=5, capsize=2,
            color=colors[N], linestyle="-", label=rf"$N={N}$",
        )
        fit_values.append(re)
        fit_errors.append(re_err)
        fit_Ns.append(N)

    if len(fit_Ns) >= 3:
        C, C_err = extrapolate_1_over_N(fit_Ns, fit_values, fit_errors)
        ax.errorbar(
            tau_common, C, yerr=C_err, errorevery=5, capsize=2,
            color="black", linestyle="--", zorder=11, label=r"$N=\infty$",
        )

    ax.set_title(rf"$\beta J_Q={beta}$")
    ax.set_xlabel(r"$\tau / \beta$")
    ax.set_xlim(0, 0.5)

axes[0].set_ylabel(ylabel)
axes[0].legend()

fig.tight_layout()
fig.savefig(outpath)
