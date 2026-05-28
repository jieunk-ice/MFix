#!/usr/bin/env python3
"""Generate a SYNTHETIC sample results directory for the dashboard demo.

The numbers here are fabricated (smooth ramps to plausible values) purely so
the published GitHub Pages dashboard renders something meaningful without a
real MFiX run. They are NOT simulation output. To publish your own results,
commit a real monitor-CSV directory and point dashboard.py / the Pages
workflow at it instead.

    python sample_results/make_sample.py [output_dir]   # default: sample_results
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd


def main(out_dir: str = "sample_results") -> int:
    os.makedirs(out_dir, exist_ok=True)
    t = np.linspace(0.0, 20.0, 200)

    def ramp(a: float, b: float) -> np.ndarray:
        """Start-up transient then steady value (steady reached by t~18s)."""
        return a + (b - a) * np.clip((t - 2.0) / 16.0, 0.0, 1.0)

    # Outlet syngas: time + gas-species mass fractions (X_g_<alias> headers).
    pd.DataFrame({
        "Time": t,
        "X_g_H2O": ramp(0.00, 0.30), "X_g_CO": ramp(0.0, 0.22),
        "X_g_H2": ramp(0.0, 0.030),  "X_g_CO2": ramp(0.0, 0.14),
        "X_g_CH4": ramp(0.0, 0.018), "X_g_N2": ramp(1.0, 0.27),
        "X_g_Tar": ramp(0.0, 0.004),
    }).to_csv(os.path.join(out_dir, "OUTLET_SYNGAS.csv"), index=False)

    # Outlet gas mass-flow rate [kg/s].
    pd.DataFrame({"Time": t, "mdot": ramp(2.0e-4, 5.0e-4)}).to_csv(
        os.path.join(out_dir, "OUTLET_MASSFLOW.csv"), index=False)

    # Bed-char inventory metric (depletes toward steady).
    pd.DataFrame({"Time": t, "char": ramp(1.0, 0.55)}).to_csv(
        os.path.join(out_dir, "BED_CHAR.csv"), index=False)

    # Cyclone loop: elutriation, return, overflow [kg/s].
    pd.DataFrame({
        "Time": t, "elut": ramp(0.0, 8.0e-5),
        "ret": ramp(0.0, 7.2e-5), "ovfl": ramp(0.0, 3.0e-5),
    }).to_csv(os.path.join(out_dir, "recirc.csv"), index=False)

    # Axial slabs: height encoded in the filename (Y010 = 0.10 m, ...).
    axial = [(10, 1180, 0.05, 0.004), (30, 1140, 0.15, 0.020),
             (50, 1120, 0.20, 0.028), (70, 1115, 0.21, 0.030),
             (90, 1113, 0.215, 0.030)]
    for h, tg, co, h2 in axial:
        pd.DataFrame({
            "Time": t, "T_g": ramp(900.0, tg),
            "X_g_CO": ramp(0.0, co), "X_g_H2": ramp(0.0, h2),
        }).to_csv(os.path.join(out_dir, f"AXIAL_Y0{h}.csv"), index=False)

    print(f"wrote synthetic sample results to {os.path.abspath(out_dir)}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else "sample_results"))
