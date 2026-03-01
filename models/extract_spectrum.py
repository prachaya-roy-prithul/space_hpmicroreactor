"""
Extract energy-resolved neutron flux spectrum from an OpenMC statepoint tally.

Inputs:
- statepoint file path
- tally name "core_flux_spectrum" with an EnergyFilter

Outputs:
- CSV file with energy-bin edges, midpoints, flux mean, and flux uncertainty
"""

import argparse
import os
import numpy as np
import openmc


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--statepoint",
        default="statepoint.100.h5",
        help="Path to the OpenMC statepoint file."
    )
    parser.add_argument(
        "--tally",
        default="core_flux_spectrum",
        help="Tally name in the statepoint file."
    )
    parser.add_argument(
        "--out",
        default="results/flux_spectrum.csv",
        help="Output CSV path."
    )
    args = parser.parse_args()

    if not os.path.exists(args.statepoint):
        raise FileNotFoundError(f"Statepoint not found: {args.statepoint}")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    sp = openmc.StatePoint(args.statepoint)
    tally = sp.get_tally(name=args.tally)

    # ----------------------------
    # Extract energy filter bins
    # ----------------------------

    energy_filter = None
    for f in tally.filters:
        if isinstance(f, openmc.EnergyFilter):
            energy_filter = f
            break

    if energy_filter is None:
        raise RuntimeError(f"No EnergyFilter found on tally '{args.tally}'.")

    bins_arr = np.asarray(energy_filter.bins, dtype=float)

    # Handle both (N,2) and edge-array formats
    if bins_arr.ndim == 2 and bins_arr.shape[1] == 2:
        e_low = bins_arr[:, 0]
        e_high = bins_arr[:, 1]
    elif bins_arr.ndim == 1 and bins_arr.size >= 2:
        e_low = bins_arr[:-1]
        e_high = bins_arr[1:]
    else:
        raise RuntimeError("Unsupported EnergyFilter bin format.")

    e_mid = np.sqrt(e_low * e_high)

    # ----------------------------
    # Extract tally results
    # ----------------------------

    mean = tally.mean.ravel()
    std = tally.std_dev.ravel()

    if mean.size != e_mid.size:
        raise RuntimeError(
            f"Size mismatch: mean.size={mean.size}, n_bins={e_mid.size}"
        )

    # ----------------------------
    # Write CSV
    # ----------------------------

    header = "E_low_eV,E_high_eV,E_mid_eV,flux_mean,flux_std\n"

    with open(args.out, "w", encoding="utf-8") as f:
        f.write(header)
        for i in range(e_mid.size):
            f.write(
                f"{e_low[i]:.10e},{e_high[i]:.10e},{e_mid[i]:.10e},"
                f"{mean[i]:.10e},{std[i]:.10e}\n"
            )

    sp.close()

    print(f"Wrote spectrum CSV: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())