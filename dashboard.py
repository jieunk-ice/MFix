#!/usr/bin/env python3
"""Generate a self-contained HTML results dashboard for an MFiX gasifier run.

Point it at a results directory (the folder holding the monitor CSVs that the
cases write) and it produces ONE HTML file with KPI cards and charts:

  * KPI cards   - H2/CO, syngas LHV, dry gas yield, cold-gas efficiency,
                  carbon conversion, tar yield
  * charts      - dry tar-free syngas composition, outlet species vs time,
                  the cyclone char-return loop, bed-char inventory vs time,
                  and axial temperature / CO / H2 profiles

Open the resulting HTML in any browser - there is no server and no dependency
beyond pandas + matplotlib (already used by postprocess.py). Every panel is
optional: if a source CSV is missing, that panel shows a placeholder so a
partial results directory still renders.

This reuses postprocess.py's CSV parsers and constants (imported, not copied)
so the numbers match the CLI report exactly.

Usage:
    python dashboard.py [results_dir] [-o dashboard.html] [--tail 0.25] [--open]
"""

from __future__ import annotations

import argparse
import base64
import datetime
import glob
import io
import os
import sys
import webbrowser

import pandas as pd
import matplotlib

matplotlib.use("Agg")  # headless: render straight to PNG bytes
import matplotlib.pyplot as plt

import postprocess as pp  # parsers (find_csv, tail_average, ...) + constants


# ---------------------------------------------------------------------------
# Metrics — mirror postprocess.py's formulas using its shared constants.
# ---------------------------------------------------------------------------

def collect_metrics(results_dir: str, tail: float) -> dict:
    """Return the scalar KPIs as {name: (value, unit) | None}."""
    m: dict[str, object] = {}

    syn = pp.find_csv(results_dir, "outlet_syngas")
    found, avg, y_wet = {}, None, {}
    if syn:
        df = pd.read_csv(syn)
        df.columns = [c.strip() for c in df.columns]
        tcol = pp.time_column(df)
        avg = pp.tail_average(df, tcol, tail)
        found = {
            sp: c
            for sp in pp.MW
            if (c := pp.detect_species_column(df, sp)) and c in avg.index
        }
        moles = {sp: avg[c] / pp.MW[sp] for sp, c in found.items()}
        tot = sum(moles.values())
        if tot > 0:
            y_wet = {sp: v / tot for sp, v in moles.items()}
            fuel = {sp: moles[sp] for sp in pp.SYNGAS if sp in moles}
            ftot = sum(fuel.values())
            if "H2" in fuel and fuel.get("CO", 0) > 0:
                m["H2/CO ratio"] = (fuel["H2"] / fuel["CO"], "")
            # dry, tar-free LHV per Nm3
            dry = {sp: y_wet[sp] for sp in y_wet if sp not in ("H2O", "Tar")}
            dtot = sum(dry.values())
            if dtot > 0:
                y_dry = {sp: v / dtot for sp, v in dry.items()}
                lhv = sum(y_dry.get(sp, 0.0) * pp.LHV_MOLAR[sp] for sp in pp.LHV_MOLAR)
                m["Syngas LHV"] = (lhv * pp.MOL_PER_NM3 / 1000.0, "MJ/Nm3")

    # Flow-dependent metrics need the outlet mass-flow monitor.
    mf = pp.find_csv(results_dir, "outlet_massflow")
    biomass_dry = pp.BIOMASS_FEED_KG_S * pp.BIOMASS_DRY_FRACTION
    if mf and avg is not None and found and "Syngas LHV" in m and biomass_dry > 0:
        dfm = pd.read_csv(mf)
        dfm.columns = [c.strip() for c in dfm.columns]
        tcm = pp.time_column(dfm)
        flow_cols = [c for c in dfm.columns if c != tcm]
        if flow_cols:
            mdot = abs(pp.tail_average(dfm, tcm, tail)[flow_cols[0]])
            mw_wet = 1.0 / sum(avg[found[sp]] / pp.MW[sp] for sp in found)
            v_wet = (mdot / mw_wet) * 22.414
            v_dry = v_wet * (1.0 - y_wet.get("H2O", 0.0) - y_wet.get("Tar", 0.0))
            lhv_nm3 = m["Syngas LHV"][0]
            m["Dry gas yield"] = (v_dry / biomass_dry, "Nm3/kg")
            m["Cold-gas efficiency"] = (
                (v_dry * lhv_nm3) / (biomass_dry * pp.BIOMASS_LHV_MJ_KG) * 100.0, "%")
            m["Tar yield"] = (
                avg.get(found.get("Tar", ""), 0.0) * mdot / biomass_dry * 1000.0, "g/kg")

    # Carbon conversion — prefer the steady-state balance from recirc.csv.
    rec = pp.find_csv(results_dir, "recirc")
    if rec:
        dr = pd.read_csv(rec)
        dr.columns = [c.strip() for c in dr.columns]
        tcol = pp.time_column(dr)
        cols = [c for c in dr.columns if c != tcol]
        if len(cols) >= 3:  # elutriation, return, overflow
            a = pp.tail_average(dr, tcol, tail)
            elut, ovfl = a[cols[0]], a[cols[2]]
            carbon_fed = pp.BIOMASS_FEED_KG_S * pp.BIOMASS_DRY_FRACTION * pp.BIOMASS_C_FRACTION
            carbon_out = ovfl + (1.0 - pp.CYCLONE_ETA) * elut
            if carbon_fed > 0:
                m["Carbon conversion"] = (100.0 * (1.0 - carbon_out / carbon_fed), "%")
    if "Carbon conversion" not in m:
        bc = pp.find_csv(results_dir, "bed_char")
        if bc:
            db = pd.read_csv(bc)
            db.columns = [c.strip() for c in db.columns]
            tcol = pp.time_column(db)
            val = [c for c in db.columns if c != tcol]
            if val:
                ch = db[val[0]]
                if ch.iloc[0] > 0:
                    m["Carbon conversion"] = (
                        100.0 * (ch.iloc[0] - ch.iloc[-1]) / ch.iloc[0], "% (batch)")
    return m


# ---------------------------------------------------------------------------
# Charts — each returns a base64 PNG data URI, or None when data is absent.
# ---------------------------------------------------------------------------

def _fig_uri(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")


def chart_composition(results_dir: str, tail: float):
    syn = pp.find_csv(results_dir, "outlet_syngas")
    if not syn:
        return None
    df = pd.read_csv(syn)
    df.columns = [c.strip() for c in df.columns]
    avg = pp.tail_average(df, pp.time_column(df), tail)
    found = {sp: c for sp in pp.SYNGAS if (c := pp.detect_species_column(df, sp)) and c in avg.index}
    moles = {sp: avg[c] / pp.MW[sp] for sp, c in found.items()}
    tot = sum(moles.values())
    if tot <= 0:
        return None
    sp_order = [s for s in pp.SYNGAS if s in moles]
    pct = [100.0 * moles[s] / tot for s in sp_order]
    fig, ax = plt.subplots(figsize=(5, 3.4))
    bars = ax.bar(sp_order, pct, color=["#2a9d8f", "#e76f51", "#264653", "#e9c46a"][: len(sp_order)])
    ax.set_ylabel("mole %  (dry, tar-free)")
    ax.set_ylim(0, max(pct) * 1.18 if pct else 1)
    for b, v in zip(bars, pct):
        ax.text(b.get_x() + b.get_width() / 2, v, f"{v:.1f}", ha="center", va="bottom", fontsize=9)
    ax.set_title("Outlet syngas composition")
    return _fig_uri(fig)


def chart_syngas_time(results_dir: str):
    syn = pp.find_csv(results_dir, "outlet_syngas")
    if not syn:
        return None
    df = pd.read_csv(syn)
    df.columns = [c.strip() for c in df.columns]
    tcol = pp.time_column(df)
    found = {sp: c for sp in pp.SYNGAS if (c := pp.detect_species_column(df, sp))}
    if not found:
        return None
    fig, ax = plt.subplots(figsize=(5, 3.4))
    for sp, c in found.items():
        ax.plot(df[tcol], df[c], label=sp)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("outlet mass fraction")
    ax.legend(fontsize=8)
    ax.set_title("Outlet species vs time")
    return _fig_uri(fig)


def chart_recirc(results_dir: str):
    rec = pp.find_csv(results_dir, "recirc")
    if not rec:
        return None
    df = pd.read_csv(rec)
    df.columns = [c.strip() for c in df.columns]
    tcol = pp.time_column(df)
    cols = [c for c in df.columns if c != tcol]
    if len(cols) < 2:
        return None
    labels = ["elutriation", "return", "overflow"]
    fig, ax = plt.subplots(figsize=(5, 3.4))
    for c, lab in zip(cols[:3], labels):
        ax.plot(df[tcol], df[c], label=lab)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("char mass flow [kg/s]")
    ax.legend(fontsize=8)
    ax.set_title("Cyclone char-return loop")
    return _fig_uri(fig)


def chart_bed_char(results_dir: str):
    bc = pp.find_csv(results_dir, "bed_char")
    if not bc:
        return None
    df = pd.read_csv(bc)
    df.columns = [c.strip() for c in df.columns]
    tcol = pp.time_column(df)
    val = [c for c in df.columns if c != tcol]
    if not val:
        return None
    fig, ax = plt.subplots(figsize=(5, 3.4))
    ax.plot(df[tcol], df[val[0]], color="#6a4c93")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("bed char inventory metric")
    ax.set_title("Bed-char inventory vs time")
    return _fig_uri(fig)


def _axial_rows(results_dir: str, tail: float):
    paths = sorted(
        f for f in glob.glob(os.path.join(results_dir, "*.csv"))
        if "axial_y" in os.path.basename(f).lower()
    )
    rows = []
    for p in paths:
        base = os.path.basename(p).lower()  # filenames may be upper- or lower-case
        digits = "".join(ch for ch in base.split("axial_y")[1] if ch.isdigit())[:3]
        height = int(digits) / 100.0 if digits else float("nan")
        df = pd.read_csv(p)
        df.columns = [c.strip() for c in df.columns]
        avg = pp.tail_average(df, pp.time_column(df), tail)
        tg = next((avg[c] for c in df.columns if c.strip().lower().startswith("t_g")), float("nan"))
        co = pp.detect_species_column(df, "CO")
        h2 = pp.detect_species_column(df, "H2")
        rows.append((height, tg, avg.get(co, float("nan")), avg.get(h2, float("nan"))))
    rows.sort()
    return rows


def chart_axial(results_dir: str, tail: float):
    rows = _axial_rows(results_dir, tail)
    if len(rows) < 2:
        return None
    hs = [r[0] for r in rows]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.4))
    ax1.plot([r[1] for r in rows], hs, "o-", color="#bc4749")
    ax1.set_xlabel("gas temperature [K]")
    ax1.set_ylabel("height [m]")
    ax1.set_title("Axial temperature")
    ax2.plot([r[2] for r in rows], hs, "o-", label="CO", color="#e76f51")
    ax2.plot([r[3] for r in rows], hs, "s-", label="H2", color="#2a9d8f")
    ax2.set_xlabel("mass fraction")
    ax2.set_title("Axial CO / H2")
    ax2.legend(fontsize=8)
    return _fig_uri(fig)


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------

_CSS = """
:root{--bg:#0f1419;--card:#1b2430;--ink:#e6edf3;--mut:#8b98a5;--acc:#2a9d8f;--brd:#2a3543}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);
font:15px/1.5 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif}
header{padding:24px 28px;border-bottom:1px solid var(--brd)}
h1{margin:0 0 4px;font-size:20px}.sub{color:var(--mut);font-size:13px}
main{padding:20px 28px;max-width:1180px;margin:0 auto}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:14px;margin-bottom:24px}
.kpi{background:var(--card);border:1px solid var(--brd);border-radius:12px;padding:16px 18px}
.kpi .v{font-size:26px;font-weight:650;color:var(--acc)}.kpi .v.na{color:var(--mut);font-size:18px}
.kpi .u{font-size:12px;color:var(--mut)}.kpi .l{font-size:12px;color:var(--mut);margin-top:6px;
text-transform:uppercase;letter-spacing:.04em}
.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(360px,1fr));gap:18px}
.panel{background:var(--card);border:1px solid var(--brd);border-radius:12px;padding:14px}
.panel.wide{grid-column:1/-1}
.panel img{width:100%;height:auto;border-radius:6px;background:#fff}
.panel .ph{color:var(--mut);font-size:13px;padding:40px 12px;text-align:center}
.panel h3{margin:0 0 10px;font-size:13px;color:var(--mut);font-weight:600;
text-transform:uppercase;letter-spacing:.04em}
footer{padding:16px 28px;color:var(--mut);font-size:12px;border-top:1px solid var(--brd);
max-width:1180px;margin:0 auto}
code{background:#10161e;padding:1px 5px;border-radius:4px;color:#cbd5e1}
"""

_KPI_ORDER = [
    "H2/CO ratio", "Syngas LHV", "Dry gas yield",
    "Cold-gas efficiency", "Carbon conversion", "Tar yield",
]


def _fmt(value: float) -> str:
    a = abs(value)
    if a != 0 and (a < 0.01 or a >= 1e4):
        return f"{value:.2e}"
    return f"{value:.2f}"


def build_html(results_dir: str, tail: float) -> str:
    metrics = collect_metrics(results_dir, tail)

    cards = []
    for name in _KPI_ORDER:
        if name in metrics:
            val, unit = metrics[name]
            cards.append(
                f'<div class="kpi"><div class="v">{_fmt(val)}'
                f'<span class="u"> {unit}</span></div><div class="l">{name}</div></div>')
        else:
            cards.append(
                f'<div class="kpi"><div class="v na">n/a</div>'
                f'<div class="l">{name}</div></div>')

    panels = [
        ("Syngas composition", chart_composition(results_dir, tail), False),
        ("Outlet species vs time", chart_syngas_time(results_dir), False),
        ("Cyclone char-return loop", chart_recirc(results_dir), False),
        ("Bed-char inventory", chart_bed_char(results_dir), False),
        ("Axial profiles", chart_axial(results_dir, tail), True),
    ]
    panel_html = []
    for title, uri, wide in panels:
        body = (f'<img src="{uri}" alt="{title}">' if uri
                else '<div class="ph">no data &mdash; source CSV not found '
                     'in this results directory</div>')
        panel_html.append(
            f'<div class="panel{" wide" if wide else ""}"><h3>{title}</h3>{body}</div>')

    # Which monitor CSVs were found / missing.
    expected = {
        "outlet_syngas": "syngas composition + performance",
        "outlet_massflow": "gas yield / cold-gas efficiency",
        "bed_char": "bed-char inventory",
        "recirc": "cyclone loop + carbon balance",
        "axial_y": "axial profiles",
    }
    found_list, missing_list = [], []
    for stem, desc in expected.items():
        present = (bool(pp.find_csv(results_dir, stem)) if stem != "axial_y"
                   else bool(glob.glob(os.path.join(results_dir, "*axial_y*.csv"))))
        (found_list if present else missing_list).append(f"<code>{stem}</code> ({desc})")

    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    abs_dir = os.path.abspath(results_dir)
    sources = "found: " + (", ".join(found_list) or "none")
    if missing_list:
        sources += "<br>missing: " + ", ".join(missing_list)

    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>MFiX gasifier dashboard &mdash; {os.path.basename(abs_dir) or abs_dir}</title>
<style>{_CSS}</style></head><body>
<header>
  <h1>MFiX biomass gasifier &mdash; results dashboard</h1>
  <div class="sub">{abs_dir} &nbsp;&bull;&nbsp; tail-averaged over last {tail:.0%} of run
   &nbsp;&bull;&nbsp; generated {now}</div>
</header>
<main>
  <div class="kpis">{''.join(cards)}</div>
  <div class="grid">{''.join(panel_html)}</div>
</main>
<footer>
  Sources &mdash; {sources}.<br>
  Performance &amp; carbon-balance KPIs use the feed assumptions in
  <code>postprocess.py</code> (BIOMASS_* / CYCLONE_ETA); keep them in sync with the deck.
</footer>
</body></html>"""


def main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("results_dir", nargs="?", default=".", help="folder with monitor CSVs")
    p.add_argument("-o", "--out", default="dashboard.html", help="output HTML path")
    p.add_argument("--tail", type=float, default=0.25, help="fraction of run to average")
    p.add_argument("--open", action="store_true", help="open the HTML in a browser when done")
    args = p.parse_args(argv)

    if not os.path.isdir(args.results_dir):
        print(f"error: {args.results_dir} is not a directory", file=sys.stderr)
        return 1

    html = build_html(args.results_dir, args.tail)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"wrote {args.out}  ({os.path.getsize(args.out) / 1024:.0f} KB) "
          f"from {os.path.abspath(args.results_dir)}")
    if args.open:
        webbrowser.open("file://" + os.path.abspath(args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
