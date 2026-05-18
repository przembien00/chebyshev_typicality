from __future__ import annotations

"""Compare Chebyshev data for `N=20` and `N=24` at zero field."""

import os
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
from matplotlib.lines import Line2D
import numpy as np


DATA_DIR = ROOT / "Data"
PLOTS_DIR = ROOT / "Plots"

BETAS = [0.2, 0.5, 1.0, 1.5, 2.0, 2.5]
SIZES = {
    20: "tab:blue",
    24: "tab:orange",
}
MARKERS = ["v", "^", "s", "x", "D", "+"]
RESCALES = [0.5, -0.5]


def import_data(physical_data: str, project_name: str = "") -> tuple[h5.File, np.ndarray]:
    """Open one HDF5 file and return the handle together with the tau grid."""

    folder = DATA_DIR / project_name if project_name else DATA_DIR
    path = folder / f"{physical_data}.hdf5"
    if not path.exists():
        raise FileNotFoundError(path)

    handle = h5.File(path, "r")
    num_time_points = int(handle["parameters"].attrs["num_TimePoints"])
    tau = np.linspace(0.0, 1.0, 2 * num_time_points)
    return handle, tau


def load_curve(size: int, beta: float, rescale: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load the xx Chebyshev correlation and its standard deviation."""

    rescale_tag = f"{rescale:.1g}"
    physical_data = f"ISO__Square_NN_PBC_N={size}__beta={beta:.2g}__rescale={rescale_tag}"
    all_data, tau = import_data(physical_data)
    try:
        g_half = np.asarray(all_data["results"]["Re_correlation"][0], dtype=float)
        g_err_half = np.asarray(all_data["results"]["Re_stds"][0], dtype=float)
    finally:
        all_data.close()

    g = np.concatenate((g_half, np.flip(g_half)))
    g_err = np.concatenate((g_err_half, np.flip(g_err_half)))
    return tau, g, g_err


def plot_rescale(rescale: float) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.3))

    for beta_index, beta in enumerate(BETAS):
        marker = MARKERS[beta_index % len(MARKERS)]
        for size, color in SIZES.items():
            tau, g, g_err = load_curve(size, beta, rescale)
            ax.errorbar(
                tau,
                g,
                yerr=g_err,
                color=color,
                linestyle="-",
                marker=marker,
                markevery=20,
                errorevery=5,
            )

    size_handles = [
        Line2D([0], [0], color=color, lw=2, label=rf"$N={size}$")
        for size, color in SIZES.items()
    ]
    beta_handles = [
        Line2D([0], [0], color="black", marker=MARKERS[i], linestyle="None", label=rf"$\beta J_Q = {beta:.2g}$")
        for i, beta in enumerate(BETAS)
    ]

    legend_sizes = ax.legend(handles=size_handles, loc="upper right", title="System size")
    ax.add_artist(legend_sizes)
    ax.legend(handles=beta_handles, loc="lower left", title="Temperature")

    ax.set_xlabel(r"$\tau/\beta$")
    ax.set_ylabel(r"$g_{xx}(\tau)$")
    ax.set_xlim(tau[0], tau[-1])
    ax.margins(x=0)

    PLOTS_DIR.mkdir(exist_ok=True)
    rescale_tag = "0p5" if rescale > 0 else "m0p5"
    output = PLOTS_DIR / f"Plot_20v24_rescale_{rescale_tag}.pdf"
    fig.savefig(output, dpi=1000, bbox_inches="tight")


def main() -> None:
    for rescale in RESCALES:
        plot_rescale(rescale)


if __name__ == "__main__":
    main()
