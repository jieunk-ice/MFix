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
| Gas (0) | H₂O, CO, H₂, CO₂, CH₄, N₂, C₆H₆(tar) | reactant/product | ideal-gas mixture |
| Solids 1 | char + ash | reacting fuel | variable (inert ash) |
| Solids 2 | silica sand | heat carrier, fluidization | constant (2600 kg/m³) |

Sand dominates the bed inventory and carries the heat; char is a small,
continuously-replenished fraction.

## 3. Repository layout

One MFiX **project per directory** (all `.f` in a project directory are
compiled into that project's custom solver):

```
gasifier_2d/   2D reacting case + usr_rates.f, usr0.f, usr1.f
gasifier_3d/   3D version + grid-independence presets (own usr_*.f copies)
coldflow_2d/   hydrodynamics-only sand case (standard solver, no build)
postprocess.py results analysis (syngas, carbon conversion, recirc loop)
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
| Wall temperature (allothermal) | 1173 K |
| Steam inlet (superficial) | 0.40 m/s, 70 % H₂O / 30 % N₂ (mass) |
| Char particle | 500 µm, 450 kg/m³ |
| Sand particle | 400 µm, 2600 kg/m³ |
| Initial voidage (bed) | 0.45 (char 0.05, sand 0.50) |
| Drag | Syamlal–O'Brien |
| Granular energy | algebraic KTGF |

### 4.3 Reactions and kinetics

Seven reactions; rate laws are in `usr_rates.f`. The **rate forms and
activation energies are literature-based**; the **pre-exponentials and LH
inhibition constants are order-of-magnitude defaults to be calibrated** to a
specific char.

| # | Reaction | ΔH | Rate form | Source |
|---|----------|----|-----------|--------|
| 1 | C + H₂O → CO + H₂ | endo | Langmuir–Hinshelwood, pᴴ²ᴼ, H₂ inhib. | Barrio & Hustad 2001 |
| 2 | C + CO₂ → 2 CO | endo | Langmuir–Hinshelwood, pᶜᴼ², CO inhib. | Barrio et al. 2001 |
| 3 | C + 2 H₂ → CH₄ | exo | mass-action, [H₂]² | — |
| 4 | CO + H₂O → CO₂ + H₂ | exo | mass-action (fwd) | Bustamante 2005 |
| 5 | CO₂ + H₂ → CO + H₂O | endo | reverse via Keq(T) | Moe 1962 |
| 6 | CH₄ + H₂O → CO + 3 H₂ | endo | mass-action | Jones & Lindstedt 1988 |
| 7 | C₆H₆ + 6 H₂O → 6 CO + 9 H₂ | endo | mass-action | Jess 1996 |

Water-gas shift (4/5) is modelled as a forward/reverse pair: the net rate
`k(cᶜᴼ·cᴴ²ᴼ − cᶜᴼ²·cᴴ²/Keq)` relaxes the gas toward equilibrium, and its sign
selects which one-way reaction carries the rate. Char reactivities use the
**solids** temperature; gas reactions use the **gas** temperature.

## 5. Operational features

- **Allothermal heat source.** Constant-temperature side walls (1173 K) supply
  the heat for the net-endothermic chemistry (`bc_tw_g`, `bc_tw_s`, Dirichlet).
- **Continuous feed (point sources).** PS1 injects char solids (2×10⁻⁴ kg/s,
  350 K, 90 % char / 10 % ash); PS2 co-injects residual volatiles as gas
  (tar + CH₄). Mimics a screw feeder + devolatilization.
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
| `usr1.f` | per step: char elutriation flux → set return point-source rate; log it |

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

## 8. Post-processing

`python postprocess.py <results_dir> [--tail 0.25] [--plot]` reads the monitor
CSVs and reports:

- dry, tar-free **syngas composition** (mole %) and **H₂/CO ratio** at the outlet;
- **carbon conversion** from the bed char-inventory monitor;
- the **cyclone recirculation loop** (mean elutriation, mean return, return
  ratio) from `recirc.csv`.

## 9. Assumptions and limitations

- **Unvalidated.** Not compiled/run here. Build in the GUI and watch the first
  run.
- **Kinetics.** Pre-exponentials and LH inhibition constants are placeholders;
  calibrate to your char/tar and the cited sources.
- **Rate units.** `usr_rates.f` assumes SI molar-rate units (kmol·m⁻³·s⁻¹);
  confirm against the bundled `silane_pyrolysis` tutorial and rescale if needed.
- **Monitor types.** `monitor_type` integers (area average vs sum vs flow-rate)
  vary between MFiX versions — verify in the GUI Monitors pane.
- **Ash thermo.** No standard database entry for ash; the ash species uses the
  **SiO₂** database name (alias `Ash`) as a stand-in.
- **Side overflow.** A side pressure-outflow also bleeds some gas; if it
  distorts the flow, switch to a specified solids mass-outflow (`MO`).
- **`usr1.f` internals.** Face-area array / velocity component for the outlet
  orientation, `IS_ON_myPE_owns`, and the `GLOBAL_ALL_SUM` interface are
  flagged in-file for version-specific verification.
- **Constant gas viscosity** (`mu_g0`) and a single constant **diffusivity**
  (`dif_g0`) are simplifications; consider temperature-dependent transport.

## 10. Recommended next steps

- Replace placeholder kinetics with measured/literature values for the target
  char; add intrinsic-vs-effective (diffusion-limited) reactivity if relevant.
- Add a true carbon-balance conversion (carbon in feed vs carbon out) once feed
  and overflow rates are set.
- Calibrate the steam/carbon ratio and bed temperature to a target syngas H₂/CO.
- Validate hydrodynamics (Uₘf, bed expansion) and, if available, syngas data.

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
