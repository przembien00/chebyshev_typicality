from __future__ import annotations

import glob
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MPL_CONFIG_DIR = ROOT / ".mplconfig"
XDG_CACHE_DIR = ROOT / ".cache"
MPL_CONFIG_DIR.mkdir(exist_ok=True)
XDG_CACHE_DIR.mkdir(exist_ok=True)

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CONFIG_DIR))
os.environ.setdefault("XDG_CACHE_HOME", str(XDG_CACHE_DIR))

import matplotlib as mpl

mpl.rc_file(str(ROOT / "matplotlibrc"))

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np

DATA_DIR = ROOT / "Data" / "Real_Time"
SPINDMFT_DIR = ROOT / "Data" / "spinDMFT_realtime"
BETA_ARRAY = [0.5, 1.5, 2.5]


def find_files():
    files = {}
    for path in glob.glob(f"{DATA_DIR}/ISO__Square_NN_PBC_N=*__beta=*__rescale=*.hdf5"):
        m = re.search(r"N=(\d+)__beta=([\d.]+)__rescale=", path)
        if m is None:
            continue
        N, beta = int(m.group(1)), float(m.group(2))
        files[(N, beta)] = path
    return files


def load(path):
    f = h5.File(path, "r")
    params = f["parameters"].attrs
    t = np.linspace(0.0, params["Tmax"], params["num_TimePoints"])
    re = f["results"]["Re_correlation"]["0-0"][0]
    re_err = f["results"]["Re_stds"]["0-0"][0]
    im = f["results"]["Im_correlation"]["0-0"][0]
    im_err = f["results"]["Im_stds"]["0-0"][0]
    return t, re, re_err, im, im_err


def load_spindmft(beta):
    path = SPINDMFT_DIR / f"spinmodel=ISO__beta={beta}__realtime.hdf5"
    if not path.exists():
        return None
    f = h5.File(path, "r")
    t = f["real_time"][:]
    c_real = f["local"]["C_real"][:]
    c_imag = f["local"]["C_imag"][:]
    return t, c_real, c_imag


def extrapolate_1_over_N(Ns, values, errors):
    """Weighted least-squares fit of y = C + a/N at each timepoint, returns C(t), C_err(t)."""
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
N_array = sorted({N for N, beta in files if beta in BETA_ARRAY})

color_cycle = [c for c in plt.rcParams["axes.prop_cycle"].by_key()["color"] if c != "red"]
colors = {N: color_cycle[i % len(color_cycle)] for i, N in enumerate(N_array)}

PARTS = [
    ("re", 1, 2, r"$\mathrm{Re}\,g^{xx}(t)$", "Plots/Plot_realtime_square_re.pdf"),
    ("im", 3, 4, r"$\mathrm{Im}\,g^{xx}(t)$", "Plots/Plot_realtime_square_im.pdf"),
]

for kind, value_idx, err_idx, ylabel, outpath in PARTS:
    fig, axes = plt.subplots(1, len(BETA_ARRAY), figsize=(5 * len(BETA_ARRAY), 4.5), sharey=True)

    for ax, beta in zip(axes, BETA_ARRAY):
        t_common = None
        fit_values, fit_errors, fit_Ns = [], [], []
        for N in N_array:
            key = (N, beta)
            if key not in files:
                continue
            data = load(files[key])
            t_common = data[0]
            ax.errorbar(
                data[0], data[value_idx], yerr=data[err_idx], errorevery=5, capsize=2,
                color=colors[N], linestyle="-", label=rf"$N={N}$",
            )
            fit_values.append(data[value_idx])
            fit_errors.append(data[err_idx])
            fit_Ns.append(N)

        if len(fit_Ns) >= 3:
            C, C_err = extrapolate_1_over_N(fit_Ns, fit_values, fit_errors)
            ax.errorbar(
                t_common, C, yerr=C_err, errorevery=5, capsize=2,
                color="black", linestyle="--", label=r"$N=\infty$", zorder=11,
            )

        spindmft = load_spindmft(beta)
        if spindmft is not None:
            t_dmft, c_real, c_imag = spindmft
            c = c_real if kind == "re" else -c_imag
            ax.plot(t_dmft, c, color="red", linestyle="-", linewidth=3, label="spinDMFT", zorder=10)

        ax.set_title(rf"$\beta J_Q={beta}$")
        ax.set_xlabel(r"$t J_Q$")
        ax.set_xlim(0, 10)

    handles, labels = axes[0].get_legend_handles_labels()
    order = sorted(range(len(labels)), key=lambda i: {r"$N=\infty$": 0, "spinDMFT": 1}.get(labels[i], 2))
    axes[0].legend([handles[i] for i in order], [labels[i] for i in order])
    axes[0].set_ylabel(ylabel)

    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)
