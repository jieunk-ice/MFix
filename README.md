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
- `run_case.sh` — build + run a single case
- `sweep.py` — parameter sweeps and a grid-independence study

One MFiX project per directory (all `.f` in a directory compile into that
project's custom solver).
