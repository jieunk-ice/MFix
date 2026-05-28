#!/usr/bin/env python3
"""Interactive Streamlit dashboard for an MFiX gasifier run (optional).

This is the interactive sibling of dashboard.py. Where dashboard.py writes a
static HTML file, this app lets you pick a results directory from a dropdown,
slide the tail-averaging window, and adjust the feed assumptions and watch
every KPI update live.

Run it with:

    pip install streamlit          # one-time, optional extra dependency
    streamlit run dashboard_app.py

(then open the URL it prints). It reuses postprocess.py's CSV parsers and
constants and dashboard.py's metric logic, so the numbers match both the CLI
report and the static dashboard.
"""

from __future__ import annotations

import glob
import os

import pandas as pd
import streamlit as st

import postprocess as pp        # parsers (find_csv, tail_average, ...) + constants
import dashboard as db          # collect_metrics + axial-row helper


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

MONITOR_STEMS = ("outlet_syngas", "outlet_massflow", "bed_char", "recirc", "axial_y")


def find_result_dirs(root: str, max_depth: int = 2) -> list[str]:
    """Directories under `root` (incl. root) that hold any monitor CSV."""
    hits = set()
    root = os.path.abspath(root)
    for dirpath, _dirs, files in os.walk(root):
        depth = dirpath[len(root):].count(os.sep)
        if depth > max_depth:
            continue
        low = [f.lower() for f in files if f.lower().endswith(".csv")]
        if any(any(stem in f for stem in MONITOR_STEMS) for f in low):
            hits.add(dirpath)
    return sorted(hits)


@st.cache_data(show_spinner=False)
def load_csv(path: str, _mtime: float) -> pd.DataFrame:
    """Read a monitor CSV (cached on path + mtime so edits invalidate it)."""
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    return df


def read_monitor(results_dir: str, stem: str) -> pd.DataFrame | None:
    path = pp.find_csv(results_dir, stem)
    if not path:
        return None
    return load_csv(path, os.path.getmtime(path))


# ---------------------------------------------------------------------------
# Sidebar — data source + parameters
# ---------------------------------------------------------------------------

st.set_page_config(page_title="MFiX gasifier dashboard", layout="wide")
st.title("MFiX biomass gasifier — results dashboard")

with st.sidebar:
    st.header("Data source")
    base = st.text_input("Search root", value=".", help="Folder to scan for run directories")
    candidates = find_result_dirs(base) if os.path.isdir(base) else []
    if candidates:
        results_dir = st.selectbox("Results directory", candidates,
                                   format_func=lambda p: os.path.relpath(p, os.path.abspath(base)))
    else:
        st.warning("No monitor CSVs found under the search root.")
        results_dir = st.text_input("Results directory (explicit path)", value=base)
    if st.button("↻ Refresh data"):
        load_csv.clear()
        st.rerun()

    st.header("Averaging")
    tail = st.slider("Tail-average window (last fraction of run)", 0.05, 0.90, 0.25, 0.05)

    st.header("Feed assumptions")
    st.caption("KPIs that depend on the feed; keep in sync with the deck.")
    pp.BIOMASS_FEED_KG_S = st.number_input("Biomass feed [kg/s]", value=float(pp.BIOMASS_FEED_KG_S),
                                           format="%.2e")
    pp.BIOMASS_DRY_FRACTION = st.number_input("Dry-biomass mass fraction", 0.0, 1.0,
                                              float(pp.BIOMASS_DRY_FRACTION), 0.05)
    pp.BIOMASS_C_FRACTION = st.number_input("Carbon mass fraction of biomass", 0.0, 1.0,
                                           float(pp.BIOMASS_C_FRACTION), 0.01)
    pp.BIOMASS_LHV_MJ_KG = st.number_input("Biomass LHV [MJ/kg]", 0.0, 40.0,
                                          float(pp.BIOMASS_LHV_MJ_KG), 0.5)
    pp.CYCLONE_ETA = st.number_input("Cyclone collection efficiency", 0.0, 1.0,
                                    float(pp.CYCLONE_ETA), 0.01)

if not results_dir or not os.path.isdir(results_dir):
    st.info("Pick a results directory in the sidebar to begin.")
    st.stop()

st.caption(f"`{os.path.abspath(results_dir)}` — tail-averaged over the last {tail:.0%} of the run")


# ---------------------------------------------------------------------------
# Source availability
# ---------------------------------------------------------------------------

present = {}
for stem in MONITOR_STEMS:
    if stem == "axial_y":
        present[stem] = bool(glob.glob(os.path.join(results_dir, "*axial_y*.csv")))
    else:
        present[stem] = bool(pp.find_csv(results_dir, stem))
have = [s for s, ok in present.items() if ok]
miss = [s for s, ok in present.items() if not ok]
cols = st.columns([3, 2])
cols[0].success("found: " + (", ".join(have) or "none"))
if miss:
    cols[1].warning("missing: " + ", ".join(miss))


# ---------------------------------------------------------------------------
# KPI cards  (reuse dashboard.collect_metrics so numbers match the CLI)
# ---------------------------------------------------------------------------

metrics = db.collect_metrics(results_dir, tail)
order = ["H2/CO ratio", "Syngas LHV", "Dry gas yield",
         "Cold-gas efficiency", "Carbon conversion", "Tar yield"]
row = st.columns(len(order))
for col, name in zip(row, order):
    if name in metrics:
        val, unit = metrics[name]
        col.metric(name, f"{val:.2f}", help=unit or None)
        if unit:
            col.caption(unit)
    else:
        col.metric(name, "n/a")

st.divider()


# ---------------------------------------------------------------------------
# Charts
# ---------------------------------------------------------------------------

c1, c2 = st.columns(2)

# Syngas composition (dry, tar-free mole %)
with c1:
    st.subheader("Syngas composition")
    df = read_monitor(results_dir, "outlet_syngas")
    if df is not None:
        avg = pp.tail_average(df, pp.time_column(df), tail)
        found = {sp: c for sp in pp.SYNGAS
                 if (c := pp.detect_species_column(df, sp)) and c in avg.index}
        moles = {sp: avg[c] / pp.MW[sp] for sp, c in found.items()}
        tot = sum(moles.values())
        if tot > 0:
            comp = pd.DataFrame(
                {"mole %": [100 * moles[s] / tot for s in pp.SYNGAS if s in moles]},
                index=[s for s in pp.SYNGAS if s in moles])
            st.bar_chart(comp, y="mole %")
        else:
            st.caption("could not compute composition from detected columns")
    else:
        st.caption("no outlet_syngas CSV")

# Outlet species vs time
with c2:
    st.subheader("Outlet species vs time")
    df = read_monitor(results_dir, "outlet_syngas")
    if df is not None:
        tcol = pp.time_column(df)
        found = {sp: c for sp in pp.SYNGAS if (c := pp.detect_species_column(df, sp))}
        if found:
            ts = df[[tcol] + list(found.values())].rename(columns={c: sp for sp, c in found.items()})
            st.line_chart(ts, x=tcol)
        else:
            st.caption("no syngas species columns detected")
    else:
        st.caption("no outlet_syngas CSV")

c3, c4 = st.columns(2)

# Cyclone char-return loop
with c3:
    st.subheader("Cyclone char-return loop")
    df = read_monitor(results_dir, "recirc")
    if df is not None:
        tcol = pp.time_column(df)
        flux = [c for c in df.columns if c != tcol]
        names = ["elutriation", "return", "overflow"][: len(flux)]
        loop = df[[tcol] + flux[: len(names)]].rename(
            columns={c: n for c, n in zip(flux, names)})
        st.line_chart(loop, x=tcol)
    else:
        st.caption("no recirc.csv")

# Bed-char inventory
with c4:
    st.subheader("Bed-char inventory")
    df = read_monitor(results_dir, "bed_char")
    if df is not None:
        tcol = pp.time_column(df)
        val = [c for c in df.columns if c != tcol]
        if val:
            st.line_chart(df[[tcol, val[0]]], x=tcol)
        else:
            st.caption("no data column")
    else:
        st.caption("no bed_char CSV")

# Axial profiles
st.subheader("Axial profiles")
rows = db._axial_rows(results_dir, tail)
if len(rows) >= 2:
    ax = pd.DataFrame(rows, columns=["height", "T_g", "CO", "H2"])
    a1, a2 = st.columns(2)
    with a1:
        st.caption("gas temperature [K] vs height [m]")
        st.scatter_chart(ax, x="T_g", y="height")
    with a2:
        st.caption("CO / H2 mass fraction vs height [m]")
        st.scatter_chart(ax.melt("height", ["CO", "H2"], "species", "mass fraction"),
                         x="mass fraction", y="height", color="species")
else:
    st.caption("need ≥2 axial_y### CSVs for a profile")

with st.expander("Raw monitor tables"):
    for stem in ("outlet_syngas", "outlet_massflow", "bed_char", "recirc"):
        df = read_monitor(results_dir, stem)
        if df is not None:
            st.markdown(f"**{stem}**")
            st.dataframe(df, height=200, use_container_width=True)

st.caption("Performance & carbon-balance KPIs use the sidebar feed assumptions "
           "(BIOMASS_* / CYCLONE_ETA). Static export: `python dashboard.py <dir> -o out.html`.")
