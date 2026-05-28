# MFix
Practicing MFix.

## Fluidized-bed biomass char steam gasifier

A staged MFiX (Two-Fluid Model) toolkit for a steam-blown bubbling fluidized-bed
gasifier. See [`report.md`](report.md) for the full model description, run
workflow, and caveats.

- `coldflow_2d/` — hydrodynamics-only sand case (Umf / bed-expansion checks)
- `gasifier_2d/` — 2D reacting, thermal case with feed, overflow, and cyclone return
- `gasifier_3d/` — 3D version with grid-independence presets
- `postprocess.py` — syngas, performance (LHV/yield/CGE), carbon conversion, axial profiles
- `dashboard.py` — one-command HTML dashboard of a run's results (KPI cards + charts)
- `dashboard_app.py` — optional interactive Streamlit dashboard (live dir/window/feed controls)
- `run_case.sh` — build + run a single case
- `sweep.py` — parameter sweeps and a grid-independence study

### Viewing results

```sh
python dashboard.py <results_dir> -o dashboard.html --open
```

Point it at the folder holding a run's monitor CSVs and it writes one
self-contained HTML file (KPI cards for H2/CO, syngas LHV, gas yield, cold-gas
efficiency, carbon conversion, tar yield, plus charts for syngas composition,
the cyclone loop, bed-char inventory, and axial profiles). No server and no
new dependencies beyond `pandas` + `matplotlib`; missing CSVs degrade to
placeholders. It reuses `postprocess.py`'s parsers and constants, so the
numbers match the CLI report.

For an interactive view (pick the run from a dropdown, slide the averaging
window, tweak the feed assumptions and watch the KPIs update):

```sh
pip install streamlit          # optional extra dependency
streamlit run dashboard_app.py
```

One MFiX project per directory (all `.f` in a directory compile into that
project's custom solver).
