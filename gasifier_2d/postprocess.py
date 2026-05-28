#!/usr/bin/env python3
"""Post-process MFiX gasifier monitor output.

Reads the time-series CSV files written by the MONITOR_* keywords in
gasifier_2d.mfx / gasifier_3d.mfx and reports:

  * dry, tar-free syngas composition at the outlet (mole %), time-averaged
    over the quasi-steady tail of the run, plus the H2/CO ratio
  * char-bed inventory vs time and a carbon-conversion estimate
  * the cyclone char-return loop (elutriation vs return) from recirc.csv

MFiX writes one CSV per monitor, named after MONITOR_NAME (e.g.
``OUTLET_SYNGAS.csv``, ``BED_CHAR.csv``). Column headers vary a little
between MFiX versions, so this script AUTODETECTS columns and prints what
it found - if a species is missed, set it explicitly in COLUMN_OVERRIDES.

Usage:
    python postprocess.py [results_dir] [--tail 0.25] [--plot]

Requires: pandas, matplotlib (matplotlib only with --plot).
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import pandas as pd

# Gas species -> molar mass [kg/kmol], matching the .mfx species ordering.
MW = {
    "H2O": 18.015,
    "CO": 28.010,
    "H2": 2.016,
    "CO2": 44.010,
    "CH4": 16.043,
    "N2": 28.013,
    "Tar": 78.114,  # C6H6 surrogate
}

# Gas-phase index (1-based) as declared in the .mfx, for "X_g_<n>" headers.
SPECIES_INDEX = {1: "H2O", 2: "CO", 3: "H2", 4: "CO2", 5: "CH4", 6: "N2", 7: "Tar"}

# Species reported on the dry, tar-free "fuel gas" basis.
SYNGAS = ["CO", "H2", "CO2", "CH4"]

# If autodetection picks the wrong column, hardcode here: {"CO": "X_g_2", ...}
COLUMN_OVERRIDES: dict[str, str] = {}


def find_csv(results_dir: str, stem: str) -> str | None:
    """Locate a monitor CSV by name stem, case-insensitively."""
    hits = [
        f
        for f in glob.glob(os.path.join(results_dir, "*.csv"))
        if stem.lower() in os.path.basename(f).lower()
    ]
    return sorted(hits)[0] if hits else None


def time_column(df: pd.DataFrame) -> str:
    for c in df.columns:
        if c.strip().lower() in ("time", "t", "time[s]", "time (s)"):
            return c
    return df.columns[0]  # MFiX puts time first


def detect_species_column(df: pd.DataFrame, name: str) -> str | None:
    """Match a species to a column by alias name or X_g_<index>."""
    if name in COLUMN_OVERRIDES:
        return COLUMN_OVERRIDES[name]
    idx = {v: k for k, v in SPECIES_INDEX.items()}.get(name)
    cols = list(df.columns)
    # Prefer an explicit alias in the header (e.g. "X_g_CO", "CO").
    for c in cols:
        token = c.replace("X_g", "").replace("_", " ").strip().upper()
        if token == name.upper():
            return c
    # Fall back to the species index, e.g. "X_g_2".
    if idx is not None:
        for c in cols:
            cc = c.lower().replace(" ", "")
            if cc.endswith(f"x_g_{idx}") or cc.endswith(f"x_g({idx})"):
                return c
    return None


def tail_average(df: pd.DataFrame, tcol: str, tail: float) -> pd.Series:
    """Average each column over the last `tail` fraction of the time span."""
    t = df[tcol]
    t_cut = t.min() + (1.0 - tail) * (t.max() - t.min())
    window = df[t >= t_cut]
    return window.mean(numeric_only=True)


def report_syngas(results_dir: str, tail: float, do_plot: bool) -> None:
    path = find_csv(results_dir, "outlet_syngas")
    if not path:
        print("  (no outlet_syngas CSV found - skipping syngas report)")
        return
    print(f"\nOutlet syngas  <-  {os.path.basename(path)}")
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    tcol = time_column(df)

    cols = {sp: detect_species_column(df, sp) for sp in MW}
    found = {sp: c for sp, c in cols.items() if c}
    print("  detected species columns:")
    for sp, c in found.items():
        print(f"    {sp:4s} -> {c}")
    missing = [sp for sp in SYNGAS if sp not in found]
    if missing:
        print(f"  WARNING: missing {missing}; set COLUMN_OVERRIDES and rerun.")

    avg = tail_average(df, tcol, tail)

    # Mass fraction -> mole fraction (per kg of mixture), then dry & tar-free.
    moles = {sp: avg[c] / MW[sp] for sp, c in found.items() if c in avg.index}
    fuel = {sp: moles[sp] for sp in SYNGAS if sp in moles}
    total = sum(fuel.values())
    print(f"\n  dry, tar-free syngas (mole %, last {tail:.0%} of run):")
    if total > 0:
        for sp in SYNGAS:
            if sp in fuel:
                print(f"    {sp:4s} {100.0 * fuel[sp] / total:6.2f} %")
        if "H2" in fuel and "CO" in fuel and fuel["CO"] > 0:
            print(f"    H2/CO ratio = {fuel['H2'] / fuel['CO']:.2f}")
    else:
        print("    (could not compute - check detected columns)")

    if do_plot:
        _plot_syngas(df, tcol, found)


def report_conversion(results_dir: str) -> None:
    path = find_csv(results_dir, "bed_char")
    if not path:
        print("\n  (no bed_char CSV found - skipping carbon conversion)")
        return
    print(f"\nCarbon conversion  <-  {os.path.basename(path)}")
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    tcol = time_column(df)
    val_cols = [c for c in df.columns if c != tcol]
    if not val_cols:
        print("  (no data column found)")
        return
    char = df[val_cols[0]]  # summed char-phase bulk density (proportional to mass)
    initial, final = char.iloc[0], char.iloc[-1]
    print(f"  initial char inventory metric: {initial:.4g}")
    print(f"  final   char inventory metric: {final:.4g}")
    if initial > 0:
        x = (initial - final) / initial
        print(f"  apparent conversion (inventory basis): {100.0 * x:6.2f} %")
    print(
        "  NOTE: with a continuous char feed, conversion is better defined as\n"
        "  (carbon in feed - carbon leaving) / carbon in feed. The inventory\n"
        "  basis above is only exact for a batch (no-feed) run."
    )


def report_recirc(results_dir: str, tail: float, do_plot: bool) -> None:
    """Summarize the cyclone char-return loop logged by usr0/usr1.f."""
    path = find_csv(results_dir, "recirc")
    if not path:
        print("\n  (no recirc.csv found - skipping char-return loop report)")
        return
    print(f"\nCyclone char-return loop  <-  {os.path.basename(path)}")
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    tcol = time_column(df)
    cols = [c for c in df.columns if c != tcol]
    if len(cols) < 2:
        print("  (unexpected columns; expected elutriation + return)")
        return
    elut_col, ret_col = cols[0], cols[1]
    avg = tail_average(df, tcol, tail)
    elut, ret = avg[elut_col], avg[ret_col]
    print(f"  mean char elutriation (last {tail:.0%}): {elut:.4g} kg/s")
    print(f"  mean char return      (last {tail:.0%}): {ret:.4g} kg/s")
    if elut > 0:
        print(f"  effective return ratio: {ret / elut:.2f}  (~cyclone efficiency)")
    if do_plot:
        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(df[tcol], df[elut_col], label="elutriation")
        plt.plot(df[tcol], df[ret_col], label="return")
        plt.xlabel("time [s]")
        plt.ylabel("char mass flow [kg/s]")
        plt.legend()
        plt.title("Cyclone char-return loop")
        out = "recirc_loop.png"
        plt.savefig(out, dpi=120, bbox_inches="tight")
        print(f"  saved plot: {out}")


def _plot_syngas(df, tcol, found):
    import matplotlib.pyplot as plt

    plt.figure()
    for sp in SYNGAS:
        if sp in found:
            plt.plot(df[tcol], df[found[sp]], label=sp)
    plt.xlabel("time [s]")
    plt.ylabel("outlet mass fraction")
    plt.legend()
    plt.title("Outlet syngas species vs time")
    out = "syngas_timeseries.png"
    plt.savefig(out, dpi=120, bbox_inches="tight")
    print(f"\n  saved plot: {out}")


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results_dir", nargs="?", default=".", help="folder with monitor CSVs")
    p.add_argument("--tail", type=float, default=0.25, help="fraction of run to average")
    p.add_argument("--plot", action="store_true", help="save a species-vs-time plot")
    args = p.parse_args(argv)

    if not os.path.isdir(args.results_dir):
        print(f"error: {args.results_dir} is not a directory", file=sys.stderr)
        return 1

    report_syngas(args.results_dir, args.tail, args.plot)
    report_conversion(args.results_dir)
    report_recirc(args.results_dir, args.tail, args.plot)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
