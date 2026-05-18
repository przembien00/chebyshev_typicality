from __future__ import annotations

"""Compare Chebyshev and CspinDMFT pair correlations for FM and AFM runs."""

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
os.environ.setdefault("MPLBACKEND", "Agg")

import h5py as h5
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parent
CSPIN_ROOT = ROOT.parent / "spindmft_matsubara" / "CspinDMFT_finite_T"
SPINDMFT_ROOT = ROOT.parent / "spindmft_matsubara" / "spinDMFT"
CHEB_DATA_DIR = ROOT / "Data"
CSPIN_DATA_DIR = CSPIN_ROOT / "Data"
SPINDMFT_DATA_DIRS = [SPINDMFT_ROOT / "Data" / "Paper", SPINDMFT_ROOT / "Data"]
PLOTS_DIR = ROOT / "Plots" / "CspinDMFT_sites"

mpl.rc_file(str(ROOT / "matplotlibrc"))


MARKERS = ["v", "^", "s", "x", "D", "p", "o", "*", "+", "<", ">"]


@dataclass(frozen=True)
class PhaseSpec:
    name: str
    cheb_rescale: float
    cspin_subdir: str
    cspin_j: float


@dataclass(frozen=True)
class SourceSpec:
    label: str
    color: str
    linestyle: str
    source: str  # "cheb" or "cspin"
    config: str | None = None
    pair: str | None = None


PHASES = [
    PhaseSpec(name="FM", cheb_rescale=-0.5, cspin_subdir="FM", cspin_j=-0.5),
    PhaseSpec(name="AFM", cheb_rescale=0.5, cspin_subdir="AFM", cspin_j=0.5),
]

PAIR_PLOTS = [
    ("0-0", "0-0", "g^{xx}"),
    ("1-0", "0-1", "g^{xx}"),
    ("7-0", "0-3", "g^{xx}"),
]

CHEB_SPINMODEL = "ISO"
CSPIN_SPINMODEL = "ISO"
CSPIN_CONFIGS = {
    "N2": "Square_2D_N=2_NN_J={j}",
    "N4": "Square_2D_N=4_NN_J={j}",
}

SOURCE_STYLES = {
    "cheb": SourceSpec(label="Chebyshev", color="darkviolet", linestyle="-", source="cheb"),
    "cspinN2": SourceSpec(label="CspinDMFT N=2", color="steelblue", linestyle=":", source="cspin"),
    "cspinN4": SourceSpec(label="CspinDMFT N=4", color="darkorange", linestyle=":", source="cspin"),
    "spindmft": SourceSpec(label="spinDMFT", color="green", linestyle=":", source="spindmft"),
}


def hdf5_to_dict(obj: h5.Group | h5.Dataset) -> Any:
    if isinstance(obj, h5.Dataset):
        return obj[...]

    data: dict[str, Any] = {}
    for key, value in obj.items():
        data[key] = hdf5_to_dict(value)
    if obj.attrs:
        data["attrs"] = {key: decode_value(value) for key, value in obj.attrs.items()}
    return data


def decode_value(value: Any) -> Any:
    if isinstance(value, bytes):
        try:
            return value.decode("utf-8")
        except UnicodeDecodeError:
            return value
    if isinstance(value, np.generic):
        return value.item()
    return value


def load_file(path: Path) -> dict[str, Any]:
    with h5.File(path, "r") as handle:
        return hdf5_to_dict(handle)


def filename_beta(path: Path) -> float | None:
    match = re.search(r"(?:^|__)beta=([^_]+)", path.stem)
    return float(match.group(1)) if match else None


def matching_cheb_files(phase: PhaseSpec) -> list[Path]:
    files: list[Path] = []
    rescale_token = f"rescale={phase.cheb_rescale:.1g}"
    for path in sorted(CHEB_DATA_DIR.glob("*.hdf5")):
        if "sites=0-1-7" not in path.name:
            continue
        if rescale_token not in path.name:
            continue
        beta = filename_beta(path)
        if beta is None:
            continue
        files.append(path)
    return sorted(files, key=lambda path: filename_beta(path) or -1.0)


def matching_cspin_files(phase: PhaseSpec, size_tag: str) -> list[Path]:
    config = CSPIN_CONFIGS[size_tag].format(j=phase.cspin_j)
    data_dir = CSPIN_DATA_DIR / phase.cspin_subdir
    files: list[Path] = []
    for path in sorted(data_dir.glob("*.hdf5")):
        if f"config={config}" not in path.name:
            continue
        beta = filename_beta(path)
        if beta is None:
            continue
        files.append(path)
    return sorted(files, key=lambda path: filename_beta(path) or -1.0)


def matching_spindmft_files() -> list[Path]:
    files: list[Path] = []
    for data_dir in SPINDMFT_DATA_DIRS:
        if not data_dir.exists():
            continue
        for path in sorted(data_dir.glob("spinmodel=ISO__beta=*.hdf5")):
            if path.name.endswith("X.hdf5") or path.name.endswith("_Matsubara.hdf5"):
                continue
            beta = filename_beta(path)
            if beta is None:
                continue
            files.append(path)
    unique: dict[float, Path] = {}
    for path in files:
        beta = filename_beta(path)
        if beta is not None and beta not in unique:
            unique[beta] = path
    return [unique[beta] for beta in sorted(unique)]


def file_map(paths: list[Path]) -> dict[float, Path]:
    result: dict[float, Path] = {}
    for path in paths:
        beta = filename_beta(path)
        if beta is not None:
            result[beta] = path
    return result


def common_betas(*maps: dict[float, Path]) -> list[float]:
    common = set(maps[0])
    for mapping in maps[1:]:
        common &= set(mapping)
    return sorted(common)


def mirror_curve(curve: np.ndarray) -> np.ndarray:
    """Mirror a half-interval curve to cover the full 0..beta interval."""

    return np.concatenate((curve, curve[::-1]))


def extract_curve(path: Path, pair: str, mirror: bool = False) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    data = load_file(path)
    params = data["parameters"]["attrs"]
    beta = float(params["beta"])
    delta_t = float(params["delta_t"])
    values = np.asarray(data["results"]["Re_correlation"][pair], dtype=float)[0]
    stds = np.asarray(data["results"]["Re_stds"][pair], dtype=float)[0]
    tau = np.arange(values.shape[-1], dtype=float) * delta_t
    if mirror:
        values = mirror_curve(values)
        stds = mirror_curve(stds)
        tau = np.linspace(0.0, beta, values.shape[-1], endpoint=True)
    return tau / beta, values, stds, beta


def extract_cspin_curve(path: Path, pair: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    data = load_file(path)
    params = data["parameters"]["attrs"]
    beta = float(params["beta"])
    delta_t = float(params["delta_t"])
    values = np.asarray(data["results"]["Re_correlation"][pair], dtype=float)[0]
    stds = np.asarray(data["runtimedata"]["Re_correlation_sample_stds"][pair], dtype=float)[0]
    tau = np.arange(values.shape[-1], dtype=float) * delta_t
    return tau / beta, values, stds, beta


def extract_spindmft_curve(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    data = load_file(path)
    params = data["parameters"]["attrs"]
    beta = float(params["beta"])
    delta_t = float(params["delta_t"])
    values = np.asarray(data["results"]["Re_correlation"], dtype=float)[0]
    stds = np.asarray(data["runtimedata"]["Re_correlation_sample_stds"], dtype=float)[0]
    tau = np.arange(values.shape[-1], dtype=float) * delta_t
    return tau / beta, values, stds, beta


def temperature_label(beta: float) -> str:
    return rf"$\beta J_Q={beta:.3g}$"


def marker_handles(markers: list[str]) -> list[Line2D]:
    return [Line2D([], [], color="black", linestyle="None", marker=m, markersize=6) for m in markers]


def source_handles(source_specs: list[SourceSpec]) -> list[Line2D]:
    return [
        Line2D([], [], color=spec.color, linestyle=spec.linestyle, linewidth=1.6, label=spec.label)
        for spec in source_specs
    ]


def pair_ylabel(pair_label: str) -> str:
    if "/" in pair_label:
        return r"$g^{xx}(\tau)$"
    return rf"$g^{{xx}}_{{{pair_label.replace('-', '')}}}(\tau)$"


def output_path(phase: PhaseSpec, cheb_pair: str, cspin_pair: str) -> Path:
    out_dir = PLOTS_DIR / phase.name
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / f"compare_cspin_sites__phase={phase.name}__cheb={cheb_pair}__cspin={cspin_pair}.pdf"


def plot_phase_pair(phase: PhaseSpec, cheb_pair: str, cspin_pair: str, cspin_sizes: list[str]) -> None:
    cheb_map = file_map(matching_cheb_files(phase))
    cspin_maps = {size: file_map(matching_cspin_files(phase, size)) for size in cspin_sizes}
    spindmft_map = file_map(matching_spindmft_files()) if cheb_pair == "0-0" and cspin_pair == "0-0" else {}

    if not cheb_map:
        raise FileNotFoundError(f"no Chebyshev files found for {phase.name}")
    if not all(cspin_maps.values()):
        missing = [size for size, mapping in cspin_maps.items() if not mapping]
        raise FileNotFoundError(f"no CspinDMFT files found for {phase.name}: {', '.join(missing)}")

    betas = common_betas(cheb_map, *cspin_maps.values())
    if not betas:
        raise FileNotFoundError(f"no common betas found for {phase.name}")

    fig, ax = plt.subplots(figsize=(8, 3.5))
    source_specs = [SOURCE_STYLES["cheb"]]
    source_specs.extend(
        SOURCE_STYLES[f"cspin{size}"]
        for size in cspin_sizes
    )
    if spindmft_map:
        source_specs.append(SOURCE_STYLES["spindmft"])
    marker_list = MARKERS[: len(betas)]

    for index, beta in enumerate(betas):
        marker = marker_list[index % len(marker_list)]

        cheb_x, cheb_y, cheb_stds, _ = extract_curve(cheb_map[beta], cheb_pair, mirror=True)
        ax.plot(
            cheb_x,
            cheb_y,
            color=SOURCE_STYLES["cheb"].color,
            linestyle=SOURCE_STYLES["cheb"].linestyle,
            marker=marker,
            markevery=max(len(cheb_x) // 12, 1),
            linewidth=1.6,
        )
        ax.fill_between(
            cheb_x,
            cheb_y - cheb_stds,
            cheb_y + cheb_stds,
            color=SOURCE_STYLES["cheb"].color,
            alpha=0.15,
            linewidth=0,
        )

        for size in cspin_sizes:
            cspin_x, cspin_y, cspin_stds, _ = extract_cspin_curve(cspin_maps[size][beta], cspin_pair)
            style = SOURCE_STYLES[f"cspin{size}"]
            ax.plot(
                cspin_x,
                cspin_y,
                color=style.color,
                linestyle=style.linestyle,
                marker=marker,
                markevery=max(len(cspin_x) // 12, 1),
                linewidth=1.6,
            )
            ax.fill_between(
                cspin_x,
                cspin_y - cspin_stds,
                cspin_y + cspin_stds,
                color=style.color,
                alpha=0.15,
                linewidth=0,
            )

        if spindmft_map and beta in spindmft_map:
            sp_x, sp_y, sp_stds, _ = extract_spindmft_curve(spindmft_map[beta])
            style = SOURCE_STYLES["spindmft"]
            ax.plot(
                sp_x,
                sp_y,
                color=style.color,
                linestyle=style.linestyle,
                marker=marker,
                markevery=max(len(sp_x) // 12, 1),
                linewidth=1.6,
            )
            ax.fill_between(
                sp_x,
                sp_y - sp_stds,
                sp_y + sp_stds,
                color=style.color,
                alpha=0.12,
                linewidth=0,
            )

    beta_legend = ax.legend(marker_handles(marker_list), [temperature_label(beta) for beta in betas], fontsize=8, loc="lower left")
    ax.add_artist(beta_legend)
    ax.legend(handles=source_handles(source_specs), fontsize=8, loc="upper right")

    ax.set_xlabel(r"$\tau/\beta$")
    ax.set_ylabel(pair_ylabel(cheb_pair if cheb_pair == cspin_pair else f"{cheb_pair}/{cspin_pair}"))
    ax.set_xlim(0.0, 1.0)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path(phase, cheb_pair, cspin_pair), dpi=300)
    plt.close(fig)


def main() -> None:
    for phase in PHASES:
        plot_phase_pair(phase, "0-0", "0-0", ["N2", "N4"])
        plot_phase_pair(phase, "1-0", "0-1", ["N2", "N4"])
        plot_phase_pair(phase, "7-0", "0-3", ["N4"])


if __name__ == "__main__":
    main()
