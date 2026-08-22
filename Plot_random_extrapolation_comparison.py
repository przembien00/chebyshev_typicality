"""Compare original and refreshed random-imaginary-time data.

Both the zero-field and finite-field panels compare the N->infinity
extrapolations.
"""

from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".mplconfig"))

import h5py
import matplotlib as mpl

mpl.rc_file(ROOT / "matplotlibrc")
mpl.rcParams["figure.autolayout"] = False
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


BETAS = [0.2, 0.5, 1.0, 1.5, 2.0, 2.5]
OLD_EXTRAP = ROOT / "Data/Random_extrapolation"
NEW_EXTRAP = ROOT / "Data/Random_extrapolation_new"
PLOTS = ROOT / "Plots"


def curve(path: Path, component: int = 0, part: str = "Re") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as handle:
        beta = float(handle["parameters"].attrs["beta"])
        y = np.asarray(handle[f"results/{part}_correlation"][component], dtype=float)
        err = np.asarray(handle[f"results/{part}_stddev"][component], dtype=float)
    tau = np.linspace(0.0, 0.5 * beta, len(y))
    return tau, y, err


def extrapolated_path(directory: Path, beta: float, h_z: float | None = None) -> Path:
    name = f"ISO__Random__N=inf__beta={beta:g}"
    if h_z is not None:
        name += f"__h_z={h_z:g}"
    return directory / f"{name}.hdf5"


def plot_zero_field() -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for index, beta in enumerate(BETAS):
        color = colors[index % len(colors)]
        tau, old, old_err = curve(extrapolated_path(OLD_EXTRAP, beta))
        _, new, new_err = curve(extrapolated_path(NEW_EXTRAP, beta))
        ax.plot(tau / beta, old, color=color, label=rf"$\beta={beta:g}$")
        ax.fill_between(tau / beta, old - old_err, old + old_err, color=color, alpha=0.12)
        ax.plot(tau / beta, new, color=color, linestyle="--", marker="o", markevery=20,
                markerfacecolor="none")
        ax.fill_between(tau / beta, new - new_err, new + new_err, color=color, alpha=0.08)
    ax.set(xlabel=r"$\tau/\beta$", ylabel=r"$g_{xx}(\tau)$", xlim=(0, 0.5))
    beta_legend = ax.legend(loc="lower left", title=r"$\beta$")
    ax.add_artist(beta_legend)
    ax.legend(handles=[Line2D([], [], color="black", label="previous"),
                       Line2D([], [], color="black", linestyle="--", marker="o",
                              markerfacecolor="none", label="new")], loc="upper right")
    fig.savefig(PLOTS / "random_extrapolation_comparison_zero_field.pdf", bbox_inches="tight")
    fig.savefig(PLOTS / "random_extrapolation_comparison_zero_field.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_field() -> None:
    fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(7.4, 8.6),
                             gridspec_kw={"hspace": 0.06})
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    components = [(0, "Re", r"$g^{xx}(\tau)$"),
                  (1, "Im", r"Im $g^{xy}(\tau)$"),
                  (3, "Re", r"$g^{zz}(\tau)$")]
    for ax, (component, part, ylabel) in zip(axes, components):
        for index, beta in enumerate(BETAS):
            color = colors[index % len(colors)]
            tau, old, old_err = curve(extrapolated_path(OLD_EXTRAP, beta, h_z=0.5), component, part)
            _, new, new_err = curve(extrapolated_path(NEW_EXTRAP, beta, h_z=0.5), component, part)
            ax.plot(tau / beta, old, color=color, label=rf"$\beta={beta:g}$")
            ax.fill_between(tau / beta, old - old_err, old + old_err, color=color, alpha=0.12)
            ax.plot(tau / beta, new, color=color, linestyle="--", marker="o", markevery=20,
                    markerfacecolor="none")
            ax.fill_between(tau / beta, new - new_err, new + new_err, color=color, alpha=0.08)
        ax.set(ylabel=ylabel, xlim=(0, 0.5))
        ax.margins(x=0)
    beta_legend = axes[0].legend(loc="lower left", title=r"$\beta$")
    axes[0].add_artist(beta_legend)
    axes[0].legend(handles=[Line2D([], [], color="black", label=r"previous $N\to\infty$"),
                            Line2D([], [], color="black", linestyle="--", marker="o",
                                   markerfacecolor="none", label=r"new $N\to\infty$")], loc="upper right")
    axes[-1].set_xlabel(r"$\tau/\beta$")
    fig.subplots_adjust(left=0.15, right=0.98, top=0.985, bottom=0.075)
    fig.savefig(PLOTS / "random_extrapolation_comparison_field_hz_0p5.pdf", bbox_inches="tight")
    fig.savefig(PLOTS / "random_extrapolation_comparison_field_hz_0p5.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    PLOTS.mkdir(exist_ok=True)
    plot_zero_field()
    plot_field()


if __name__ == "__main__":
    main()
