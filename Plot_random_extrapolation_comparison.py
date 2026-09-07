"""Compare original and refreshed random-imaginary-time extrapolations.

Only beta/field groups present in both extrapolation directories are plotted.
This deliberately excludes incomplete source-data groups, which cannot support
the weighted linear 1/N fit.
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
OLD_INPUT = ROOT / "Data/Random_imagtime"
NEW_INPUT = ROOT / "Data/Random_imagtime_new"
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


def common_betas(h_z: float | None = None) -> list[float]:
    """Return beta points that have valid old and new extrapolated files."""
    return [beta for beta in BETAS
            if extrapolated_path(OLD_EXTRAP, beta, h_z).exists()
            and extrapolated_path(NEW_EXTRAP, beta, h_z).exists()]


def input_path(directory: Path, size: int, beta: float, h_z: float) -> Path | None:
    """Return the unique raw-data file for a parameter point, if available."""
    pattern = f"ISO__Random__N={size}__beta={beta:g}__h_z={h_z:g}__numConfigs=*.hdf5"
    paths = list(directory.glob(pattern))
    if len(paths) > 1:
        raise RuntimeError(f"multiple files match {directory / pattern}: {paths}")
    return paths[0] if paths else None


def common_input_betas(size: int, h_z: float, excluded: set[float] = set()) -> list[float]:
    """Return shared raw-data beta values, excluding explicitly omitted points."""
    return [beta for beta in BETAS
            if beta not in excluded
            and input_path(OLD_INPUT, size, beta, h_z) is not None
            and input_path(NEW_INPUT, size, beta, h_z) is not None]


def plot_zero_field() -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    betas = common_betas()
    if not betas:
        raise RuntimeError("no common zero-field extrapolations to plot")
    for beta in betas:
        color = colors[BETAS.index(beta) % len(colors)]
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
    betas = common_betas(h_z=0.5)
    if not betas:
        print("No finite-field comparison written: no new h_z=0.5 group has at least three sizes.")
        return
    fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(7.4, 8.6),
                             gridspec_kw={"hspace": 0.06})
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    components = [(0, "Re", r"$g^{xx}(\tau)$"),
                  (1, "Im", r"Im $g^{xy}(\tau)$"),
                  (3, "Re", r"$g^{zz}(\tau)$")]
    for ax, (component, part, ylabel) in zip(axes, components):
        for beta in betas:
            color = colors[BETAS.index(beta) % len(colors)]
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


def plot_field_n10() -> None:
    """Compare shared h_z=0.5 raw data at fixed N=10, omitting beta=0.2."""
    size, h_z = 10, 0.5
    betas = common_input_betas(size, h_z, excluded={0.2})
    if not betas:
        raise RuntimeError("no shared N=10 field data to plot")
    fig, axes = plt.subplots(nrows=3, ncols=1, sharex=True, figsize=(7.4, 8.6),
                             gridspec_kw={"hspace": 0.06})
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    components = [(0, "Re", r"$g^{xx}(\tau)$"),
                  (1, "Im", r"Im $g^{xy}(\tau)$"),
                  (3, "Re", r"$g^{zz}(\tau)$")]
    for ax, (component, part, ylabel) in zip(axes, components):
        for beta in betas:
            color = colors[BETAS.index(beta) % len(colors)]
            old_path = input_path(OLD_INPUT, size, beta, h_z)
            new_path = input_path(NEW_INPUT, size, beta, h_z)
            assert old_path is not None and new_path is not None
            tau, old, old_err = curve(old_path, component, part)
            _, new, new_err = curve(new_path, component, part)
            ax.plot(tau / beta, old, color=color, label=rf"$\beta={beta:g}$")
            ax.fill_between(tau / beta, old - old_err, old + old_err, color=color, alpha=0.12)
            ax.plot(tau / beta, new, color=color, linestyle="--", marker="o", markevery=20,
                    markerfacecolor="none")
            ax.fill_between(tau / beta, new - new_err, new + new_err, color=color, alpha=0.08)
        ax.set(ylabel=ylabel, xlim=(0, 0.5))
        ax.margins(x=0)
    beta_legend = axes[0].legend(loc="lower left", title=r"$\beta$")
    axes[0].add_artist(beta_legend)
    axes[0].legend(handles=[Line2D([], [], color="black", label=r"previous $N=10$"),
                            Line2D([], [], color="black", linestyle="--", marker="o",
                                   markerfacecolor="none", label=r"new $N=10$")], loc="upper right")
    axes[-1].set_xlabel(r"$\tau/\beta$")
    fig.subplots_adjust(left=0.15, right=0.98, top=0.985, bottom=0.075)
    fig.savefig(PLOTS / "random_input_comparison_field_hz_0p5_N10.pdf", bbox_inches="tight")
    fig.savefig(PLOTS / "random_input_comparison_field_hz_0p5_N10.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    PLOTS.mkdir(exist_ok=True)
    plot_zero_field()
    plot_field()
    plot_field_n10()


if __name__ == "__main__":
    main()
