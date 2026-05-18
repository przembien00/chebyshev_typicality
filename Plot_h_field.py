from __future__ import annotations

"""Compare Chebyshev data for `N=20` and `N=24` at `h_z=0.5`."""

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

BETAS = [0.2, 0.5, 1.0, 1.5, 2.5]
SIZES = {
    20: "tab:blue",
    24: "tab:orange",
}
MARKERS = ["v", "^", "s", "x", "D"]
RESCALES = [0.5, -0.5]
FIELD_TAG = "h_z=0.5"


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


def build_physical_data(size: int, beta: float, rescale: float) -> str:
    """Format the filename stem for one data point."""

    return f"ISO__Square_NN_PBC_N={size}__beta={beta:.2g}__h_z=0.5__rescale={rescale:.1g}"


def load_components(size: int, beta: float, rescale: float) -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Load xx, xy and zz components together with their errors."""

    all_data, tau = import_data(build_physical_data(size, beta, rescale))
    try:
        re_corr = np.asarray(all_data["results"]["Re_correlation"], dtype=float)
        im_corr = np.asarray(all_data["results"]["Im_correlation"], dtype=float)
        re_stds = np.asarray(all_data["results"]["Re_stds"], dtype=float)
        im_stds = np.asarray(all_data["results"]["Im_stds"], dtype=float)
    finally:
        all_data.close()

    curves = {
        "xx": np.concatenate((re_corr[0], np.flip(re_corr[0]))),
        "xy": np.concatenate((im_corr[1], -np.flip(im_corr[1]))),
        "zz": np.concatenate((re_corr[3], np.flip(re_corr[3]))),
    }
    errors = {
        "xx": np.concatenate((re_stds[0], np.flip(re_stds[0]))),
        "xy": np.concatenate((im_stds[1], np.flip(im_stds[1]))),
        "zz": np.concatenate((re_stds[3], np.flip(re_stds[3]))),
    }
    return tau, curves, errors


def plot_rescale(rescale: float) -> None:
    """Create one stacked comparison plot for a given rescale value."""

    fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(8, 9), gridspec_kw={"hspace": 0})
    labels = [rf"$\beta J_Q = {beta:.2g}$" for beta in BETAS]

    for beta_index, beta in enumerate(BETAS):
        marker = MARKERS[beta_index % len(MARKERS)]
        for size, color in SIZES.items():
            tau, curves, errors = load_components(size, beta, rescale)
            axes[0].errorbar(
                tau,
                curves["xx"],
                yerr=errors["xx"],
                color=color,
                linestyle="-",
                marker=marker,
                markevery=20,
                errorevery=5,
            )
            axes[1].errorbar(
                tau,
                curves["xy"],
                yerr=errors["xy"],
                color=color,
                linestyle="-",
                marker=marker,
                markevery=20,
                errorevery=5,
            )
            axes[2].errorbar(
                tau,
                curves["zz"],
                yerr=errors["zz"],
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
        Line2D([0], [0], color="black", marker=MARKERS[i], linestyle="None", label=label)
        for i, label in enumerate(labels)
    ]

    axes[0].legend(handles=size_handles, loc="upper right", title="System size")
    axes[1].legend(handles=beta_handles, loc="upper right", title="Temperature")

    axes[0].set_ylabel(r"$g^{xx}(\tau)$")
    axes[1].set_ylabel(r"Im $g^{xy}(\tau)$")
    axes[2].set_ylabel(r"$g^{zz}(\tau)$")
    axes[2].set_xlabel(r"$\tau/\beta$")

    for ax in axes:
        ax.set_xlim(tau[0], tau[-1])
        ax.margins(x=0)

    PLOTS_DIR.mkdir(exist_ok=True)
    rescale_tag = "0p5" if rescale > 0 else "m0p5"
    output = PLOTS_DIR / f"Plot_20v24_{FIELD_TAG}_rescale_{rescale_tag}.pdf"
    fig.savefig(output, dpi=1000, bbox_inches="tight")


def main() -> None:
    for rescale in RESCALES:
        plot_rescale(rescale)


if __name__ == "__main__":
    main()
