# Fluidized-Bed Gasifier for Biomass Char Steam Gasification — MFiX Model Report

## 1. Objective

Set up an MFiX model of a bubbling **fluidized-bed gasifier** in which
**biomass char** is gasified with **steam** to produce a syngas rich in CO
and H₂. The model is built up in stages — from cold-flow hydrodynamics to a
fully reacting, thermal, continuously-fed gasifier with a cyclone char-return
loop — so each layer of physics can be verified before the next is added.

> **Status:** these cases were authored as a learning/starter toolkit. They
> have **not** been compiled or run against the MFiX solver here. Treat the
> kinetics and several keyword choices (flagged below) as *to-be-verified*
> before drawing quantitative conclusions.

## 2. Modeling approach

MFiX (Multiphase Flow with Interphase eXchanges, NETL/US-DOE) is used with the
**Two-Fluid Model (TFM)**: gas and solids are treated as interpenetrating
Eulerian continua, closed with kinetic theory of granular flow. TFM is the
standard, comparatively inexpensive choice for dense bubbling-bed gasifiers.

Three Eulerian phases:

| Phase | Material | Role | Density model |
|------:|----------|------|---------------|
| Gas (0) | H₂O, CO, H₂, CO₂, CH₄, N₂, C₆H₆(tar), O₂ | reactant/product | ideal-gas mixture |
| Solids 1 | biomass, moisture, char, ash | reacting fuel | variable (inert ash) |
| Solids 2 | silica sand | heat carrier, fluidization | constant (2600 kg/m³) |

Sand dominates the bed inventory and carries the heat. Raw biomass is fed
continuously; it dries and devolatilizes in situ to char + volatiles, and the
char then gasifies.

## 3. Repository layout

One MFiX **project per directory** (all `.f` in a project directory are
compiled into that project's custom solver):

```
gasifier_2d/   2D reacting case + usr_rates.f, usr0.f, usr1.f
gasifier_3d/   3D version + grid-independence presets (own usr_*.f copies)
coldflow_2d/   hydrodynamics-only sand case (standard solver, no build)
postprocess.py results analysis (syngas, performance, carbon, recirc, axial)
dashboard.py   one-command HTML results dashboard (KPI cards + charts)
dashboard_app.py  optional interactive Streamlit dashboard
run_case.sh    build + run one case
sweep.py       parameter sweeps + grid-independence study (drives the solver)
report.md      this document
```

## 4. Physical model

### 4.1 Geometry and mesh

| | 2D (`gasifier_2d`) | 3D (`gasifier_3d`) |
|---|---|---|
| Domain | 0.10 × 1.00 × 0.01 m | 0.10 × 1.00 × 0.10 m |
| Grid | 20 × 100 × 1 | 10×50×10 / 15×75×15 / 20×100×20 |
| Initial bed | 0–0.30 m | 0–0.30 m |

### 4.2 Operating conditions

| Quantity | Value |
|---|---|
| Pressure | 101.3 kPa (atmospheric) |
| Bed/inlet temperature | 1123 K (850 °C) |
| Walls | adiabatic (autothermal); 1173 K fixed-T option for allothermal |
| Inlet (superficial) | 0.40 m/s; 55 % H₂O / 25 % O₂ / 20 % N₂ (mass) |
| Char particle | 500 µm, 450 kg/m³ |
| Sand particle | 400 µm, 2600 kg/m³ |
| Initial voidage (bed) | 0.45 (char 0.05, sand 0.50) |
| Drag | Syamlal–O'Brien |
| Granular energy | algebraic KTGF |

### 4.3 Reactions and kinetics

Seven reactions; rate laws are in `usr_rates.f`. The **rate forms,
activation energies, and pre-exponentials are literature values** (sources
cited per line, table below). They are **interim stand-ins to be replaced with
experimental kinetics** fitted to the target char/tar — see the RECALL banner
in `usr_rates.f`. Activation energies are the better-established part;
pre-exponentials vary by orders of magnitude between chars.

| # | Reaction | ΔH | Rate form | Source |
|---|----------|----|-----------|--------|
| 1 | Moisture → H₂O (drying) | endo | 1st order in moisture | Chan et al. 1985 |
| 2 | Biomass → char + CO/CO₂/CH₄/H₂/H₂O/tar (pyrolysis) | endo | 1st order in biomass | Chan et al. 1985 |
| 3 | C + H₂O → CO + H₂ | endo | Langmuir–Hinshelwood, pᴴ²ᴼ, H₂ inhib. | Barrio & Hustad 2001 |
| 4 | C + CO₂ → 2 CO | endo | Langmuir–Hinshelwood, pᶜᴼ², CO inhib. | Barrio et al. 2001 |
| 5 | C + 2 H₂ → CH₄ | exo | mass-action, [H₂]² | Biba et al. 1978 |
| 6 | CO + H₂O → CO₂ + H₂ | exo | mass-action (fwd) | Bustamante 2005 |
| 7 | CO₂ + H₂ → CO + H₂O | endo | reverse via Keq(T) | Moe 1962 |
| 8 | CH₄ + H₂O → CO + 3 H₂ | endo | mass-action | Jones & Lindstedt 1988 |
| 9 | C₆H₆ + 6 H₂O → 6 CO + 9 H₂ | endo | mass-action | Jess 1996 |
| 10 | C + O₂ → CO₂ | exo | mass-action | Smith 1982 |
| 11 | CO + ½ O₂ → CO₂ | exo | mass-action | Westbrook & Dryer 1981 |
| 12 | H₂ + ½ O₂ → H₂O | exo | mass-action | (representative) |
| 13 | CH₄ + 2 O₂ → CO₂ + 2 H₂O | exo | mass-action | Westbrook & Dryer 1981 |

Reactions 10–13 are the **autothermal** heat source (partial oxidation).
Pyrolysis is mass/atom-balanced for a wood surrogate **Biomass = CH₁.₄O₀.₆**
(MW 23.02): `Biomass → 0.30 C + 0.37 CO + 0.10 CO₂ + 0.05 CH₄ + 0.48 H₂ +
0.03 H₂O + 0.03 C₆H₆`.

Water-gas shift (4/5) is modelled as a forward/reverse pair: the net rate
`k(cᶜᴼ·cᴴ²ᴼ − cᶜᴼ²·cᴴ²/Keq)` relaxes the gas toward equilibrium, and its sign
selects which one-way reaction carries the rate. Char reactivities use the
**solids** temperature; gas reactions use the **gas** temperature.

## 5. Operational features

- **Autothermal heat source.** A steam + O₂ inlet drives partial-oxidation
  (combustion) reactions whose heat sustains the endothermic gasification, so
  the side walls are run **adiabatic**. The O₂ fraction sets the equivalence
  ratio / bed temperature; the fixed-temperature wall block (`bc_tw_*`) is left
  commented for reverting to allothermal/hybrid operation.
- **Continuous biomass feed (point source).** PS1 injects raw wet biomass
  (350 K; 75 % dry biomass / 10 % moisture / 15 % ash by mass). Drying and
  pyrolysis (reactions 1–2) then release the moisture and volatiles in situ —
  there is no longer an artificial injected volatile stream.
- **Solids overflow drain.** The side wall is split to leave a pressure-outflow
  opening at the bed surface (≈0.28–0.32 m): once the bed expands past it,
  solids spill over, bounding the inventory so the case can reach steady state.
- **Cyclone char return (recirculation).** Char fines elutriate out the top;
  PS3 re-injects captured fines, hot and carbon-rich, low in the bed.
  `usr1.f` makes this **dynamic**: each step it integrates the char flux
  leaving the top outlet and sets the return rate to η·(that flux), η = 0.9.

## 6. User Fortran

| File | Purpose |
|------|---------|
| `usr_rates.f` | the seven reaction rates (LH + mass-action, reversible WGS) |
| `usr0.f` | one-time: create `recirc.csv` log (I/O rank) |
| `usr1.f` | per step: integrate char flux out the top (elutriation) and the overflow → set return rate (η·elutriation); log both for the carbon balance |

All require `call_usr = .True.` (set in the reacting decks). Setting
`call_usr = .False.` cleanly disables `usr0/usr1` (the static return rate then
applies).

## 7. How to run

Recommended **staged bring-up** (don't start from the full reacting case):

1. **Cold flow** — open `coldflow_2d/coldflow_2d.mfx` (standard solver, no
   build). Sweep the inlet velocity and read the minimum fluidization velocity
   Uₘf where the bed pressure drop plateaus; compare to a correlation (Wen & Yu).
2. **+ Energy** — in `gasifier_2d`, run with reactions off to check the thermal
   field and wall heating.
3. **+ Reactions** — build the custom solver (compiles `usr_rates.f`) and run;
   watch CO/H₂ evolve.
4. **+ Feed / walls / recirculation** — enable the point sources and `call_usr`
   (as shipped) and run to quasi-steady state.
5. **3D + grid study** — `gasifier_3d`, coarse → medium → fine; confirm the
   integral outputs (syngas, conversion) stop changing.

Automation (`sweep.py`, after `conda activate mfix-<version>`):
- `./run_case.sh <case_dir> <project.mfx>` builds and runs a single case.
- `python sweep.py param <case_dir> <mfx> "<keyword>" v1,v2,... --metric h2co`
  runs a **parameter sweep** (e.g. inlet O₂ fraction, velocity, bed height)
  and writes a response curve. Use this to tune the steam/carbon ratio or
  bed temperature toward a target H₂/CO.
- `python sweep.py grid <case_dir> <mfx> --metric h2co` runs the
  **grid-independence study** (coarse/medium/fine) and reports the change
  between successive meshes.

## 8. Post-processing

`python postprocess.py <results_dir> [--tail 0.25] [--plot] [--umf sweep.csv]`
reads the monitor CSVs and reports:

- dry, tar-free **syngas composition** (mole %) and **H₂/CO ratio** at the outlet;
- **performance metrics** — syngas **LHV** (MJ/Nm³), dry **gas yield**
  (Nm³/kg biomass), **cold-gas efficiency**, and **tar yield** (g/kg); these use
  the outlet gas mass-flow monitor and the `BIOMASS_*` constants in the script;
- **carbon conversion** — both an inventory-basis estimate (batch only) and a
  **steady-state carbon balance** `X_C = 1 − (overflow + (1−η)·elutriation) /
  (biomass carbon fed)`, using the char fluxes logged by `usr1.f`;
- the **cyclone recirculation loop** (mean elutriation, return, overflow,
  return ratio) from `recirc.csv`;
- an **axial profile** of temperature / CO / H₂ from the `axial_y###` monitors;
- **minimum fluidization velocity** from a cold-flow sweep file (`--umf`): for
  each inlet velocity record the steady bed pressure drop; Umf is where the
  rising dP meets the fluidized plateau.

`usr1.f` logs four columns to `recirc.csv` — time, top elutriation, cyclone
return, and side-overflow char flux — so the balance closes from the run plus
the feed constants (`CHAR_FEED_KG_S`, `CHAR_FRACTION`, `CYCLONE_ETA` in
`postprocess.py`, which must match the deck).

### Dashboard

For a visual, at-a-glance view, `python dashboard.py <results_dir> -o
dashboard.html --open` writes one self-contained HTML page: KPI cards (H₂/CO,
syngas LHV, dry gas yield, cold-gas efficiency, carbon conversion, tar yield)
and charts for syngas composition, outlet species vs time, the cyclone
char-return loop, bed-char inventory, and the axial T/CO/H₂ profiles. It
imports `postprocess.py`'s parsers and constants — so the figures match the CLI
report — needs no server and no dependency beyond `pandas` + `matplotlib`, and
shows a placeholder for any monitor CSV that is absent (so a partial results
directory still renders).

An optional interactive variant, `dashboard_app.py`, runs the same metrics
under Streamlit (`pip install streamlit; streamlit run dashboard_app.py`): a
sidebar lets you pick the run directory from a dropdown, slide the
tail-averaging window, and adjust the feed assumptions (BIOMASS_* / cyclone
efficiency) with the KPIs and charts updating live. It reuses the same parsers
and `dashboard.py`'s metric logic, so the numbers match the static export and
the CLI; Streamlit is the only extra dependency.

The static dashboard can also be published on **GitHub Pages**: the workflow
`.github/workflows/pages.yml` builds `dashboard.py` and deploys the HTML on
each push to the default branch (enable Settings → Pages → Source: GitHub
Actions). It builds from a synthetic sample dataset
(`sample_results/make_sample.py`) by default; point `RESULTS_DIR` at a
committed real monitor-CSV directory to publish actual results. Streamlit
cannot run on Pages (static hosting only).

## 9. Assumptions and limitations

- **Unvalidated.** Not compiled/run here. Build in the GUI and watch the first
  run.
- **Kinetics.** Literature values are in place (cited in `usr_rates.f`) as an
  interim stand-in; **recall to replace with experimental kinetics** fitted to
  your char/tar. The LH inhibition constants in particular are approximate.
- **Rate units.** `usr_rates.f` assumes SI molar-rate units (kmol·m⁻³·s⁻¹);
  confirm against the bundled `silane_pyrolysis` tutorial and rescale if needed.
- **Monitor types.** `monitor_type` integers (area average vs sum vs flow-rate)
  vary between MFiX versions — verify in the GUI Monitors pane.
- **Ash thermo.** No standard database entry for ash; the ash species uses the
  **SiO₂** database name (alias `Ash`) as a stand-in.
- **Biomass thermo.** The `Biomass` surrogate (CH₁.₄O₀.₆) is not in any
  database — you MUST supply its Cp(T) and heat of formation (consistent with
  the biomass HHV) or the pyrolysis enthalpy / energy balance will be wrong.
- **Side overflow.** A side pressure-outflow also bleeds some gas; if it
  distorts the flow, switch to a specified solids mass-outflow (`MO`).
- **`usr1.f` internals.** Face-area array / velocity component for the outlet
  orientation, `IS_ON_myPE_owns`, and the `GLOBAL_ALL_SUM` interface are
  flagged in-file for version-specific verification.
- **Constant gas viscosity** (`mu_g0`) and a single constant **diffusivity**
  (`dif_g0`) are simplifications; consider temperature-dependent transport.

## 10. Recommended next steps

- **Replace the interim literature kinetics** (now in `usr_rates.f`) with
  values fitted to the target char/tar (TGA / lab gasification data); add
  intrinsic-vs-effective (diffusion-limited) reactivity if relevant.
- Calibrate the steam/carbon ratio and bed temperature to a target syngas H₂/CO
  using `sweep.py param` (sweeps `bc_v_g`, inlet O₂/steam fraction, bed height).
- Validate hydrodynamics (Uₘf via `postprocess.py --umf`, bed expansion) and,
  if available, syngas data.

A steady-state carbon balance (carbon in feed vs. solid carbon out) is now
computed by `postprocess.py` from the char fluxes logged by `usr1.f`.

## 11. References

- M. Barrio, J.E. Hustad, *CO₂ / steam gasification of birch char* (2001).
- C.R. Bustamante et al., *High-temperature kinetics of the homogeneous
  (reverse) water-gas shift reaction*, AIChE J. (2004, 2005).
- W.P. Jones, R.P. Lindstedt, *Global reaction schemes for hydrocarbon
  combustion*, Combust. Flame (1988).
- A. Jess, *Catalytic and thermal conversion of aromatic hydrocarbons (tar)*
  (1996).
- O.A. Moe, *Water-gas shift equilibrium* (1962).
- A. Gómez-Barea, B. Leckner, *Modeling of biomass gasification in fluidized
  bed*, Prog. Energy Combust. Sci. (2010).
- M. Syamlal, T.J. O'Brien, *The derivation of a drag coefficient formula …*;
  MFiX documentation, mfix.netl.doe.gov.
