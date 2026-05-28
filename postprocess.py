#!/usr/bin/env python3
"""Post-process MFiX gasifier monitor output.

Reads the time-series CSV files written by the MONITOR_* keywords in
gasifier_2d.mfx / gasifier_3d.mfx and reports:

  * dry, tar-free syngas composition at the outlet (mole %), time-averaged
    over the quasi-steady tail of the run, plus the H2/CO ratio
  * char-bed inventory vs time and a carbon-conversion estimate
  * the cyclone char-return loop (elutriation vs return) from recirc.csv
  * performance metrics: syngas LHV, dry gas yield, cold-gas efficiency,
    tar yield (needs the outlet mass-flow monitor + the BIOMASS_* constants)
  * an axial T/CO/H2 profile from the axial_y### monitors
  * minimum fluidization velocity from a cold-flow sweep (--umf sweep.csv)

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
    "O2": 31.998,
}

# Gas-phase index (1-based) as declared in the .mfx, for "X_g_<n>" headers.
SPECIES_INDEX = {1: "H2O", 2: "CO", 3: "H2", 4: "CO2", 5: "CH4", 6: "N2", 7: "Tar", 8: "O2"}

# Species reported on the dry, tar-free "fuel gas" basis.
SYNGAS = ["CO", "H2", "CO2", "CH4"]

# If autodetection picks the wrong column, hardcode here: {"CO": "X_g_2", ...}
COLUMN_OVERRIDES: dict[str, str] = {}

# Feed / cyclone parameters. MUST match the .mfx PS settings and usr1.f ETA.
# The fuel is now RAW biomass (surrogate CH1.4O0.6, carbon mass fraction
# 12.011/23.02 = 0.5217), fed at BIOMASS_FEED_KG_S with BIOMASS_DRY_FRACTION
# dry solid (the rest moisture + ash carried separately in the PS).
BIOMASS_FEED_KG_S = 2.0e-4    # ps_massflow_s(1,1): total feed-phase rate
BIOMASS_DRY_FRACTION = 0.75   # dry-biomass mass fraction of the feed
BIOMASS_C_FRACTION = 0.5217   # carbon mass fraction of CH1.4O0.6
BIOMASS_LHV_MJ_KG = 18.0      # dry biomass lower heating value [MJ/kg]
CYCLONE_ETA = 0.90            # cyclone collection efficiency (usr1.f ETA)

# Lower heating values of the fuel gas species [kJ/mol].
LHV_MOLAR = {"CO": 283.0, "H2": 241.8, "CH4": 802.3}
MOL_PER_NM3 = 1000.0 / 22.414  # mol per normal m3 (0 C, 1 atm)


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
    ovfl_col = cols[2] if len(cols) > 2 else None
    avg = tail_average(df, tcol, tail)
    elut, ret = avg[elut_col], avg[ret_col]
    print(f"  mean char elutriation (last {tail:.0%}): {elut:.4g} kg/s")
    print(f"  mean char return      (last {tail:.0%}): {ret:.4g} kg/s")
    if elut > 0:
        print(f"  effective return ratio: {ret / elut:.2f}  (~cyclone efficiency)")
    if ovfl_col is not None:
        print(f"  mean char overflow    (last {tail:.0%}): {avg[ovfl_col]:.4g} kg/s")
        _carbon_balance(elut, avg[ovfl_col], tail)
    if do_plot:
        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(df[tcol], df[elut_col], label="elutriation")
        plt.plot(df[tcol], df[ret_col], label="return")
        if ovfl_col is not None:
            plt.plot(df[tcol], df[ovfl_col], label="overflow")
        plt.xlabel("time [s]")
        plt.ylabel("char mass flow [kg/s]")
        plt.legend()
        plt.title("Cyclone char-return loop")
        out = "recirc_loop.png"
        plt.savefig(out, dpi=120, bbox_inches="tight")
        print(f"  saved plot: {out}")


def _carbon_balance(elut: float, ovfl: float, tail: float) -> None:
    """Steady-state solid-carbon conversion from the logged char fluxes.

    Solid carbon leaves the system as overflow plus the uncaptured fraction
    of elutriated char ((1-eta) escapes the cyclone). Char is pure carbon,
    so these mass flows are carbon flows.
        X_C = 1 - (overflow + (1-eta)*elutriation) / (fresh char carbon fed)
    """
    carbon_fed = BIOMASS_FEED_KG_S * BIOMASS_DRY_FRACTION * BIOMASS_C_FRACTION
    carbon_out = ovfl + (1.0 - CYCLONE_ETA) * elut
    print(f"\n  steady-state carbon balance (last {tail:.0%}):")
    print(f"    biomass carbon fed      : {carbon_fed:.4g} kg/s")
    print(f"    solid carbon leaving    : {carbon_out:.4g} kg/s")
    if carbon_fed > 0:
        x = 1.0 - carbon_out / carbon_fed
        print(f"    carbon conversion       : {100.0 * x:6.2f} %")
    print("    (set BIOMASS_FEED_KG_S / *_FRACTION / CYCLONE_ETA to match the deck)")


def report_performance(results_dir: str, tail: float) -> None:
    """Syngas LHV, gas yield, cold-gas efficiency, and tar yield."""
    syn = find_csv(results_dir, "outlet_syngas")
    if not syn:
        print("\n  (no outlet_syngas CSV - skipping performance metrics)")
        return
    print("\nPerformance metrics")
    df = pd.read_csv(syn)
    df.columns = [c.strip() for c in df.columns]
    tcol = time_column(df)
    avg = tail_average(df, tcol, tail)

    cols = {sp: detect_species_column(df, sp) for sp in MW}
    found = {sp: c for sp, c in cols.items() if c and c in avg.index}
    # Wet mole amounts per kg mixture (mass fraction / MW).
    moles = {sp: avg[c] / MW[sp] for sp, c in found.items()}
    tot = sum(moles.values())
    if tot <= 0:
        print("  (could not compute - check outlet columns)")
        return
    y_wet = {sp: m / tot for sp, m in moles.items()}

    # Dry, tar-free basis (includes N2 diluent) for LHV per Nm3.
    dry = {sp: y_wet[sp] for sp in y_wet if sp not in ("H2O", "Tar")}
    dtot = sum(dry.values())
    y_dry = {sp: v / dtot for sp, v in dry.items()} if dtot > 0 else {}
    lhv = sum(y_dry.get(sp, 0.0) * LHV_MOLAR[sp] for sp in LHV_MOLAR)  # kJ/mol dry
    lhv_nm3 = lhv * MOL_PER_NM3 / 1000.0  # MJ/Nm3
    print(f"  syngas LHV (dry)        : {lhv_nm3:6.2f} MJ/Nm3")

    # Flow-dependent metrics need the outlet gas mass flow monitor.
    mf = find_csv(results_dir, "outlet_massflow")
    biomass_dry = BIOMASS_FEED_KG_S * BIOMASS_DRY_FRACTION
    if not mf:
        print("  (no outlet_massflow CSV - gas yield / CGE need MON 3)")
        return
    dfm = pd.read_csv(mf)
    dfm.columns = [c.strip() for c in dfm.columns]
    tcm = time_column(dfm)
    flow_cols = [c for c in dfm.columns if c != tcm]
    mdot = abs(tail_average(dfm, tcm, tail)[flow_cols[0]])  # kg/s gas out
    mw_wet = 1.0 / sum(avg[found[sp]] / MW[sp] for sp in found)  # kg/kmol
    ndot = mdot / mw_wet                       # kmol/s gas out (wet)
    v_wet = ndot * 22.414                       # Nm3/s wet
    dry_mole_frac = 1.0 - y_wet.get("H2O", 0.0) - y_wet.get("Tar", 0.0)
    v_dry = v_wet * dry_mole_frac               # Nm3/s dry
    if biomass_dry > 0:
        yield_nm3 = v_dry / biomass_dry
        cge = (v_dry * lhv_nm3) / (biomass_dry * BIOMASS_LHV_MJ_KG) * 100.0
        tar_yield = avg.get(found.get("Tar", ""), 0.0) * mdot / biomass_dry * 1000.0
        print(f"  dry gas yield           : {yield_nm3:6.2f} Nm3/kg biomass")
        print(f"  cold-gas efficiency     : {cge:6.1f} %")
        print(f"  tar yield               : {tar_yield:6.2f} g/kg biomass")
    print("  (set BIOMASS_* constants and verify the mass-flow monitor type)")


def report_axial(results_dir: str, tail: float, do_plot: bool) -> None:
    """Assemble the axial monitors (axial_y###) into T/CO/H2 vs height."""
    paths = sorted(
        f for f in glob.glob(os.path.join(results_dir, "*.csv"))
        if "axial_y" in os.path.basename(f).lower()
    )
    if not paths:
        print("\n  (no axial_y### CSVs - skipping axial profile)")
        return
    print("\nAxial profile (tail-averaged)")
    print("  height[m]   T_g[K]   x_CO     x_H2")
    rows = []
    for p in paths:
        name = os.path.basename(p)
        digits = "".join(ch for ch in name.split("axial_y")[1] if ch.isdigit())[:3]
        height = int(digits) / 100.0 if digits else float("nan")
        df = pd.read_csv(p)
        df.columns = [c.strip() for c in df.columns]
        tcol = time_column(df)
        avg = tail_average(df, tcol, tail)
        tg = next((avg[c] for c in df.columns if c.strip().lower().startswith("t_g")), float("nan"))
        co = detect_species_column(df, "CO")
        h2 = detect_species_column(df, "H2")
        rows.append((height, tg, avg.get(co, float("nan")), avg.get(h2, float("nan"))))
    rows.sort()
    for h, tg, co, h2 in rows:
        print(f"  {h:7.2f}   {tg:7.1f}  {co:7.4f}  {h2:7.4f}")
    if do_plot and len(rows) > 1:
        import matplotlib.pyplot as plt

        hs = [r[0] for r in rows]
        plt.figure()
        plt.plot([r[1] for r in rows], hs, "o-")
        plt.xlabel("gas temperature [K]")
        plt.ylabel("height [m]")
        plt.title("Axial temperature profile")
        out = "axial_temperature.png"
        plt.savefig(out, dpi=120, bbox_inches="tight")
        print(f"  saved plot: {out}")


def detect_umf(csv_path: str) -> None:
    """Estimate Umf from a cold-flow sweep CSV (columns: velocity, dP).

    Build this file yourself from a series of coldflow_2d runs: for each
    inlet velocity record the steady bed pressure drop. Umf is the velocity
    where the rising dP meets the fluidized plateau (= bed weight / area).
    """
    if not os.path.isfile(csv_path):
        print(f"\n  (no sweep file {csv_path} - skipping Umf)")
        return
    df = pd.read_csv(csv_path)
    df.columns = [c.strip().lower() for c in df.columns]
    v = df[df.columns[0]].to_numpy(dtype=float)
    dp = df[df.columns[1]].to_numpy(dtype=float)
    order = v.argsort()
    v, dp = v[order], dp[order]
    plateau = dp[len(dp) // 2:].mean()           # high-velocity plateau
    rising = dp < 0.95 * plateau                  # packed-bed (rising) points
    print("\nMinimum fluidization velocity (from sweep)")
    if rising.sum() >= 2:
        slope = dp[rising][-1] / v[rising][-1]    # dP/U through the origin
        umf = plateau / slope if slope > 0 else float("nan")
        print(f"  plateau dP : {plateau:.1f} Pa")
        print(f"  Umf (est)  : {umf:.3f} m/s")
    else:
        print("  (need more sub-fluidization points to locate Umf)")


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
    p.add_argument("--plot", action="store_true", help="save plots (PNG)")
    p.add_argument("--umf", metavar="CSV", help="cold-flow sweep file (velocity,dP) for Umf")
    args = p.parse_args(argv)

    if not os.path.isdir(args.results_dir):
        print(f"error: {args.results_dir} is not a directory", file=sys.stderr)
        return 1

    report_syngas(args.results_dir, args.tail, args.plot)
    report_performance(args.results_dir, args.tail)
    report_conversion(args.results_dir)
    report_recirc(args.results_dir, args.tail, args.plot)
    report_axial(args.results_dir, args.tail, args.plot)
    if args.umf:
        detect_umf(args.umf)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
