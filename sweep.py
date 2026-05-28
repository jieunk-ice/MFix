#!/usr/bin/env python3
"""Automation for the MFiX gasifier cases: parameter sweeps and a
grid-independence study.

For each case in a sweep this script copies the case into a fresh run
directory, edits one (or a few) keyword(s) in the .mfx, runs the MFiX
solver, and reads a scalar metric from the monitor CSVs. Results are
written to a CSV (and optional plot) so you get a response curve.

Requires an MFiX environment on PATH (build_mfixsolver / mfixsolver), e.g.
after `conda activate mfix-<version>`. Runs can be long; use --dry-run to
preview the generated cases without launching the solver.

Examples
--------
  # Sweep steam/O2 inlet velocity in the cold-flow case, read bed dP:
  python sweep.py param coldflow_2d coldflow_2d.mfx "bc_v_g(1)" \
      0.05,0.10,0.15,0.20,0.30 --metric bed_dp

  # Sweep inlet O2 fraction in the 2D gasifier, read H2/CO:
  python sweep.py param gasifier_2d gasifier_2d.mfx "bc_x_g(1,8)" \
      0.15,0.20,0.25,0.30 --metric h2co

  # Grid-independence study on the 3D case (coarse/medium/fine):
  python sweep.py grid gasifier_3d gasifier_3d.mfx --metric h2co
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys

import pandas as pd

import postprocess as pp  # reuse the analysis helpers


# Grid presets for the grid study: (imax, jmax, kmax).
GRID_PRESETS = {
    "coarse": (10, 50, 10),
    "medium": (15, 75, 15),
    "fine": (20, 100, 20),
}


def set_keyword(text: str, key: str, value) -> str:
    """Replace a scalar keyword assignment `key = ...` in an .mfx string."""
    pat = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*)(\S+)", re.MULTILINE)
    new, n = pat.subn(rf"\g<1>{value}", text, count=1)
    if n == 0:
        raise KeyError(f"keyword '{key}' not found in deck")
    return new


def prepare_run(case_dir: str, mfx: str, run_dir: str, edits: dict) -> None:
    """Copy the case into run_dir and apply keyword edits to the deck."""
    if os.path.exists(run_dir):
        shutil.rmtree(run_dir)
    shutil.copytree(case_dir, run_dir)
    path = os.path.join(run_dir, mfx)
    with open(path) as f:
        text = f.read()
    for key, value in edits.items():
        text = set_keyword(text, key, value)
    with open(path, "w") as f:
        f.write(text)


def run_solver(run_dir: str, mfx: str) -> None:
    """Build (if user Fortran present) and run the solver in run_dir."""
    has_udf = any(f.startswith("usr") and f.endswith(".f") for f in os.listdir(run_dir))
    if has_udf:
        subprocess.run(["build_mfixsolver"], cwd=run_dir, check=True)
        subprocess.run([os.path.join(".", "mfixsolver"), "-f", mfx], cwd=run_dir, check=True)
    else:
        subprocess.run(["mfixsolver", "-f", mfx], cwd=run_dir, check=True)


def get_metric(run_dir: str, metric: str, tail: float) -> float:
    """Read a scalar response metric from a finished run's monitor CSVs."""
    if metric == "bed_dp":
        path = pp.find_csv(run_dir, "bed_dp")
        if not path:
            return float("nan")
        df = pd.read_csv(path)
        df.columns = [c.strip() for c in df.columns]
        tcol = pp.time_column(df)
        col = [c for c in df.columns if c != tcol][0]
        return float(pp.tail_average(df, tcol, tail)[col])

    # Composition-based metrics from outlet_syngas.
    path = pp.find_csv(run_dir, "outlet_syngas")
    if not path:
        return float("nan")
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    tcol = pp.time_column(df)
    avg = pp.tail_average(df, tcol, tail)
    moles = {}
    for sp in ("CO", "H2", "CO2", "CH4"):
        c = pp.detect_species_column(df, sp)
        if c and c in avg.index:
            moles[sp] = avg[c] / pp.MW[sp]
    if metric == "h2co":
        return moles.get("H2", float("nan")) / moles.get("CO", float("nan"))
    tot = sum(moles.values())
    sp = metric.upper()
    if sp in moles and tot > 0:
        return 100.0 * moles[sp] / tot  # dry, tar-free mole %
    raise ValueError(f"unknown metric '{metric}'")


def _finish(rows, out_csv, plot, xlabel, ylabel):
    df = pd.DataFrame(rows, columns=[xlabel, ylabel])
    df.to_csv(out_csv, index=False)
    print(f"\nwrote {out_csv}")
    print(df.to_string(index=False))
    if plot and len(rows) > 1:
        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(df[xlabel], df[ylabel], "o-")
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(f"{ylabel} vs {xlabel}")
        png = os.path.splitext(out_csv)[0] + ".png"
        plt.savefig(png, dpi=120, bbox_inches="tight")
        print(f"saved plot: {png}")


def cmd_param(args) -> int:
    values = [v.strip() for v in args.values.split(",")]
    rows = []
    for v in values:
        run_dir = f"sweep_{_safe(args.key)}_{v}"
        print(f"\n=== {args.key} = {v}  ->  {run_dir} ===")
        prepare_run(args.case_dir, args.mfx, run_dir, {args.key: v})
        if args.dry_run:
            print("  (dry run: prepared, not executed)")
            continue
        run_solver(run_dir, args.mfx)
        rows.append((float(v), get_metric(run_dir, args.metric, args.tail)))
    if rows:
        _finish(rows, f"sweep_{_safe(args.key)}.csv", args.plot, args.key, args.metric)
    return 0


def cmd_grid(args) -> int:
    rows = []
    for name, (imax, jmax, kmax) in GRID_PRESETS.items():
        run_dir = f"grid_{name}"
        cells = imax * jmax * kmax
        print(f"\n=== grid {name}: {imax}x{jmax}x{kmax} = {cells} cells ===")
        prepare_run(args.case_dir, args.mfx, run_dir,
                    {"imax": imax, "jmax": jmax, "kmax": kmax})
        if args.dry_run:
            print("  (dry run: prepared, not executed)")
            continue
        run_solver(run_dir, args.mfx)
        rows.append((cells, get_metric(run_dir, args.metric, args.tail)))
    if len(rows) >= 3:
        _richardson(rows, args.metric)
    if rows:
        _finish(rows, "grid_study.csv", args.plot, "cells", args.metric)
    return 0


def _richardson(rows, metric):
    """Report the change between successive grids (a convergence check)."""
    rows = sorted(rows)  # by cell count
    print("\ngrid-independence check:")
    for (c0, m0), (c1, m1) in zip(rows, rows[1:]):
        if m0 == m0 and m1 == m1 and m0 != 0:  # not NaN
            print(f"  {c0} -> {c1} cells: {metric} change {100*(m1-m0)/abs(m0):+.2f}%")
    print("  (mesh is adequate once successive changes fall within tolerance)")


def _safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("case_dir")
    common.add_argument("mfx")
    common.add_argument("--metric", default="h2co",
                        help="response metric: h2co | co | h2 | co2 | ch4 | bed_dp")
    common.add_argument("--tail", type=float, default=0.25)
    common.add_argument("--plot", action="store_true")
    common.add_argument("--dry-run", action="store_true",
                        help="prepare run dirs but do not launch the solver")

    pp_ = sub.add_parser("param", parents=[common], help="sweep one keyword")
    pp_.add_argument("key", help="deck keyword, e.g. \"bc_v_g(1)\"")
    pp_.add_argument("values", help="comma-separated values, e.g. 0.1,0.2,0.3")
    pp_.set_defaults(func=cmd_param)

    g = sub.add_parser("grid", parents=[common], help="coarse/medium/fine grid study")
    g.set_defaults(func=cmd_grid)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
