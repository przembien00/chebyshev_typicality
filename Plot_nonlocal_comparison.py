#!/usr/bin/env python3
"""
Compare imaginary-time spin correlations across system sizes N=18, 20, 24.
Produces 2 figures (AFM / FM), each with local, NN, NNN stacked vertically.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.lines as mlines
import h5py
import os

plt.rcParams["figure.autolayout"] = False

DATA_DIR = "Data/NonlocalCorrs"
OUT_DIR  = "Plots/NonlocalCorrs"
os.makedirs(OUT_DIR, exist_ok=True)

SIZES = {
    18: {"sites_str": "0-1-4", "nnn_key": "4-0"},
    20: {"sites_str": "0-1-6", "nnn_key": "6-0"},
    24: {"sites_str": "0-1-7", "nnn_key": "7-0"},
}

CORR_SPECS = [
    dict(name="local", key_fn=lambda N: "0-0",                  ylabel=r"$g^{xx}_\mathrm{local}(\tau)$"),
    dict(name="NN",    key_fn=lambda N: "1-0",                  ylabel=r"$g^{xx}_\mathrm{NN}(\tau)$"),
    dict(name="NNN",   key_fn=lambda N: SIZES[N]["nnn_key"],     ylabel=r"$g^{xx}_\mathrm{NNN}(\tau)$"),
]

MAG_TYPES = {"AFM": 0.5, "FM": -0.5}   # mag_type → rescale

BETAS = sorted({0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0})

SIZE_LSTYLE = {18: "-", 20: "--", 24: ":"}
CMAP        = cm.plasma_r
BETA_NORM   = mcolors.Normalize(vmin=min(BETAS), vmax=max(BETAS))


def _mirror(data, stds, beta_f, n_pts):
    tau_half = np.linspace(0, beta_f / 2, n_pts)
    tau_full = np.concatenate([tau_half, beta_f - tau_half[-2::-1]])
    C_full   = np.concatenate([data, data[-2::-1]])
    S_full   = np.concatenate([stds, stds[-2::-1]])
    return tau_full / beta_f, C_full, S_full


def load_corr(N, beta, rescale, corr_key):
    cfg  = SIZES[N]
    base = (
        f"ISO__Square_NN_PBC_N={N}__sites={cfg['sites_str']}"
        f"__beta={beta:g}__rescale={rescale:g}"
    )
    path_x = os.path.join(DATA_DIR, base + "X.hdf5")
    path   = os.path.join(DATA_DIR, base + ".hdf5")
    path   = path_x if os.path.exists(path_x) else path
    if not os.path.exists(path):
        return None, None, None
    with h5py.File(path, "r") as f:
        n_pts  = int(f["parameters"].attrs["num_TimePoints"])
        beta_f = float(f["parameters"].attrs["beta"])
        data   = f[f"results/Re_correlation/{corr_key}"][0]
        stds   = f[f"results/Re_stds/{corr_key}"][0]
    return _mirror(data, stds, beta_f, n_pts)


for mag_type, rescale in MAG_TYPES.items():
    fig, axes = plt.subplots(3, 1, figsize=(7, 7.5), sharex=True,
                             gridspec_kw={"hspace": 0})
    fig.subplots_adjust(left=0.11, right=0.85, top=0.95, bottom=0.09)

    for row, (ax, spec) in enumerate(zip(axes, CORR_SPECS)):
        ymin_all, ymax_all = np.inf, -np.inf

        for beta in BETAS:
            color = CMAP(BETA_NORM(beta))
            for N, ls in SIZE_LSTYLE.items():
                key = spec["key_fn"](N)
                tau, C, S = load_corr(N, beta, rescale, key)
                if tau is None:
                    continue
                ax.plot(tau, C, color=color, lw=1.2, ls=ls)
                ax.fill_between(tau, C - S, C + S, color=color, alpha=0.12)
                ymin_all = min(ymin_all, np.min(C))
                ymax_all = max(ymax_all, np.max(C))

        ax.axhline(0, color="gray", lw=0.5, ls=":")
        ax.set_xlim(0, 1)
        if np.isfinite(ymin_all) and np.isfinite(ymax_all):
            margin = 0.05 * (ymax_all - ymin_all)
            ax.set_ylim(ymin_all - margin, ymax_all + margin)
        ax.set_ylabel(spec["ylabel"])

        if row < 2:
            ax.tick_params(labelbottom=False)

    axes[-1].set_xlabel(r"$\tau/\beta$")

    n_handles = [
        mlines.Line2D([], [], color="k", ls=ls, lw=1.5, label=rf"$N={N}$")
        for N, ls in SIZE_LSTYLE.items()
    ]
    axes[0].legend(handles=n_handles, fontsize=9, loc="upper right", framealpha=0.85)

    cbar_ax = fig.add_axes([0.88, 0.09, 0.025, 0.95 - 0.09])
    sm = cm.ScalarMappable(cmap=CMAP, norm=BETA_NORM)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label(r"$\beta J_Q$", fontsize=13)
    cbar.set_ticks(BETAS)
    cbar.set_ticklabels([f"{b:g}" for b in BETAS])

    fname = os.path.join(OUT_DIR, f"nonlocal_comparison__{mag_type}.pdf")
    plt.savefig(fname)
    print(f"Saved {fname}")
    plt.close()

print("Done.")
