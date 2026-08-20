"""Build 1/N -> infinity extrapolated correlation files from the N=9..13 data.

For every (beta, h_z) group in Data/Random_imagtime, fit each correlator
component and time point to  y = C + a/N  by weighted least squares (weights
1/stddev^2), and store the intercept C as the correlation and its fit error as
the stddev.  Output mirrors the input HDF5 layout so existing tooling can read
the extrapolated files the same way.
"""

import glob
import os
import re

import h5py as h5
import numpy as np

SRC_DIR = "Data/Random_imagtime"
OUT_DIR = "Data/Random_extrapolation"
DATASETS = ["Re_correlation", "Re_stddev", "Im_correlation", "Im_stddev"]


def extrapolate_1_over_N(Ns, values, errors):
    """Weighted LS fit y = C + a/N along axis 0. values/errors: (k, ...) -> C, C_err (...)."""
    x = np.asarray([1.0 / N for N in Ns], dtype=float)
    x = x.reshape((-1,) + (1,) * (values.ndim - 1))
    errors = np.where(errors == 0, 1e-12, errors)
    w = 1.0 / errors**2
    S = np.sum(w, axis=0)
    Sx = np.sum(w * x, axis=0)
    Sxx = np.sum(w * x**2, axis=0)
    Sy = np.sum(w * values, axis=0)
    Sxy = np.sum(w * x * values, axis=0)
    Delta = S * Sxx - Sx**2
    C = (Sxx * Sy - Sx * Sxy) / Delta
    C_err = np.sqrt(Sxx / Delta)
    return C, C_err


# Group source files by (beta_string, h_z_string-or-None)
groups = {}
for path in glob.glob(f"{SRC_DIR}/ISO__Random__N=*__beta=*__numConfigs=*.hdf5"):
    m = re.search(r"N=(\d+)__beta=([\d.]+)(__h_z=([\d.]+))?", os.path.basename(path))
    N = int(m.group(1))
    beta_str = m.group(2)
    hz_str = m.group(4)  # None for zero-field files
    groups.setdefault((beta_str, hz_str), {})[N] = path

os.makedirs(OUT_DIR, exist_ok=True)

for (beta_str, hz_str), by_N in sorted(groups.items(), key=lambda kv: (float(kv[0][0]), kv[0][1] or "")):
    Ns = sorted(by_N)
    if len(Ns) < 3:
        print(f"skip beta={beta_str} h_z={hz_str}: only N={Ns}")
        continue

    # Load every correlator dataset for each N, stacked along axis 0.
    stacked = {ds: [] for ds in DATASETS}
    for N in Ns:
        with h5.File(by_N[N], "r") as f:
            for ds in DATASETS:
                stacked[ds].append(f["results"][ds][()])
    stacked = {ds: np.array(v) for ds, v in stacked.items()}

    re_C, re_err = extrapolate_1_over_N(Ns, stacked["Re_correlation"], stacked["Re_stddev"])
    im_C, im_err = extrapolate_1_over_N(Ns, stacked["Im_correlation"], stacked["Im_stddev"])

    out_name = f"ISO__Random__N=inf__beta={beta_str}"
    if hz_str is not None:
        out_name += f"__h_z={hz_str}"
    out_name += ".hdf5"
    out_path = os.path.join(OUT_DIR, out_name)

    src_ref = by_N[Ns[-1]]  # largest N, used as attribute template
    with h5.File(src_ref, "r") as fsrc, h5.File(out_path, "w") as fout:
        params = fout.create_group("parameters")
        for k, v in fsrc["parameters"].attrs.items():
            params.attrs[k] = v
        # Provenance: mark as extrapolation rather than a single system size.
        params.attrs["num_Spins"] = -1
        params.attrs["extrapolated"] = np.int32(1)
        params.attrs["extrap_N_values"] = np.array(Ns, dtype=np.int32)

        results = fout.create_group("results")
        out_data = {
            "Re_correlation": re_C, "Re_stddev": re_err,
            "Im_correlation": im_C, "Im_stddev": im_err,
        }
        for ds in DATASETS:
            d = results.create_dataset(ds, data=out_data[ds])
            if "info" in fsrc["results"][ds].attrs:
                d.attrs["info"] = fsrc["results"][ds].attrs["info"]

    print(f"wrote {out_path}  (N={Ns}, shape={re_C.shape})")
