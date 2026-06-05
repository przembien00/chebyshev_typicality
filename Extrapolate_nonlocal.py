#!/usr/bin/env python3
"""
Fit f(N) = A + B/N at each imaginary-time point to extrapolate correlations
to N=infinity.  Writes one HDF5 file per (beta, rescale) combination to
Data/NonlocalCorrs/ with the same structure as the finite-N files.
"""

import warnings
import numpy as np
import h5py
import os
from scipy.optimize import curve_fit, OptimizeWarning

DATA_DIR = "Data/NonlocalCorrs"

SIZES = {
    18: {"sites_str": "0-1-4", "nnn_key": "4-0"},
    20: {"sites_str": "0-1-6", "nnn_key": "6-0"},
    24: {"sites_str": "0-1-7", "nnn_key": "7-0"},
}
CORR_LABELS = ["local", "NN", "NNN"]
CORR_KEYS   = {"local": "0-0", "NN": "1-0", "NNN": None}   # NNN filled per N

BETAS   = [0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
RESCALES = [0.5, -0.5]


def load_raw(N, beta, rescale, corr_key):
    cfg = SIZES[N]
    fname = (
        f"ISO__Square_NN_PBC_N={N}__sites={cfg['sites_str']}"
        f"__beta={beta:g}__rescale={rescale:g}.hdf5"
    )
    path = os.path.join(DATA_DIR, fname)
    if not os.path.exists(path):
        return None, None
    with h5py.File(path, "r") as f:
        C = f[f"results/Re_correlation/{corr_key}"][0].copy()
        S = f[f"results/Re_stds/{corr_key}"][0].copy()
    return C, S


def fit_inf(beta, rescale):
    """
    For each tau point fit A + B/N using the three system sizes.
    Returns dict: corr_label -> (C_inf, S_inf) each of length num_TimePoints.
    """
    Ns = np.array(list(SIZES.keys()), dtype=float)   # [18, 20, 24]

    results = {}
    for label in CORR_LABELS:
        C_arr, S_arr = [], []
        ok = True
        for N in SIZES:
            key = SIZES[N]["nnn_key"] if label == "NNN" else CORR_KEYS[label]
            C, S = load_raw(N, beta, rescale, key)
            if C is None:
                ok = False
                break
            C_arr.append(C)
            S_arr.append(S)

        if not ok:
            results[label] = None
            continue

        C_arr = np.array(C_arr)   # shape (3, n_pts)
        S_arr = np.array(S_arr)

        n_pts = C_arr.shape[1]
        A_vals = np.empty(n_pts)
        A_errs = np.empty(n_pts)

        def model(N, A, B):
            return A + B / N

        for i in range(n_pts):
            try:
                sigma = np.maximum(S_arr[:, i], 1e-12)   # avoid div-by-zero at tau=0
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", OptimizeWarning)
                    popt, pcov = curve_fit(
                        model, Ns, C_arr[:, i],
                        sigma=sigma, absolute_sigma=True,
                    )
                A_vals[i] = popt[0]
                # pcov may be inf when all points lie on the curve exactly
                A_errs[i] = 0.0 if not np.isfinite(pcov[0, 0]) else np.sqrt(pcov[0, 0])
            except Exception:
                A_vals[i] = np.nan
                A_errs[i] = np.nan

        results[label] = (A_vals, A_errs)

    return results


def save_hdf5(beta, rescale, results):
    fname = f"ISO__Square_NN_PBC_N=inf__beta={beta:g}__rescale={rescale:g}.hdf5"
    path  = os.path.join(DATA_DIR, fname)

    # copy reference attrs from N=18 file for the same rescale
    ref_cfg = SIZES[18]
    ref_fname = (
        f"ISO__Square_NN_PBC_N=18__sites={ref_cfg['sites_str']}"
        f"__beta={beta:g}__rescale={rescale:g}.hdf5"
    )
    ref_path = os.path.join(DATA_DIR, ref_fname)

    with h5py.File(path, "w") as out:
        params = out.create_group("parameters")
        # copy reference attrs where available, then override
        if os.path.exists(ref_path):
            with h5py.File(ref_path, "r") as ref:
                for k, v in ref["parameters"].attrs.items():
                    params.attrs[k] = v
        params.attrs["num_Spins"]  = np.inf
        params.attrs["beta"]       = float(beta)
        params.attrs["rescale"]    = float(rescale)
        params.attrs["spin_sites"] = b"0,1,inf"

        re_corr = out.create_group("results/Re_correlation")
        re_stds = out.create_group("results/Re_stds")
        out.create_group("results/Im_correlation")
        out.create_group("results/Im_stds")
        out.create_group("runtime_data")

        key_map = {"local": "0-0", "NN": "1-0", "NNN": "nnn-0"}
        for label, data in results.items():
            if data is None:
                continue
            C_inf, S_inf = data
            k = key_map[label]
            re_corr.create_dataset(k, data=C_inf[np.newaxis, :])
            re_stds.create_dataset(k, data=S_inf[np.newaxis, :])

    print(f"  saved {fname}")


if __name__ == "__main__":
    for rescale in RESCALES:
        print(f"rescale = {rescale:+g}")
        for beta in BETAS:
            res = fit_inf(beta, rescale)
            if any(v is not None for v in res.values()):
                save_hdf5(beta, rescale, res)
