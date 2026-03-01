# Monte Carlo Neutronic Study of Reflector, Control Rod, and Coolant Effects in a Conceptual Fast Microreactor

## Overview

This project presents a parametric Monte Carlo neutronic study of a compact fast-spectrum microreactor model. The main objective is to evaluate how reflector thickness, control rod insertion, and primary coolant selection influence the effective multiplication factor (k_eff), reflector worth, leakage fraction, and shutdown margin.

The work is intended as a conceptual neutronic investigation rather than a full reactor design study. The geometry is simplified, assumptions are clearly stated, and no validation against experimental data is claimed. The goal is to explore trends and sensitivities in a physically consistent way using OpenMC eigenvalue calculations.

---

## Physical and Computational Model

The reactor is modeled as a compact cylindrical fast-spectrum core with:

- Metallic HEU surrogate fuel,
- Six sodium heat-pipe channels arranged in a ring,
- A central control rod channel,
- A BeO radial reflector of variable thickness.

Two control rod states are analyzed:

- Withdrawn (void rod channel),
- Inserted (B4C absorber).

Primary coolant options investigated:

- Sodium (Na),
- Potassium (K),
- Sodium–Potassium eutectic (NaK),
- Lithium (Li).

All simulations are performed using the Monte Carlo neutron transport code OpenMC in eigenvalue mode.

---

## Parameters Studied

The following parametric studies are included:

### 1. Reflector Thickness Sweep
- Reflector thickness varied over a defined range.
- k_eff computed for both withdrawn and inserted rod states.
- Critical reflector thickness estimated by linear interpolation.
- Reflector worth calculated in pcm relative to 0 cm baseline.
- Leakage fraction derived from boundary current tallies.

### 2. Control Rod Reactivity Effects
- Reactivity difference between withdrawn and inserted states evaluated.
- Shutdown margin defined as:

  Shutdown margin = ρ_withdrawn − ρ_inserted (pcm)

  where ρ = (k − 1)/k.

- Uncertainty propagation performed using Monte Carlo statistical uncertainties.

### 3. Coolant Comparison Study
At a fixed reflector thickness (9 cm) and fuel temperature (600 K):

- k_eff comparison across Na, K, NaK, and Li.
- Leakage fraction comparison.
- Shutdown margin comparison with propagated uncertainties.

The coolant comparison is neutronic-only and does not include thermal-hydraulic or material compatibility considerations.

---

## Assumptions and Scope

The study is intentionally limited in scope. Major assumptions include:

- Steady-state eigenvalue analysis only.
- Single fuel temperature case (600 K).
- No burnup or fuel depletion.
- No thermal-hydraulic coupling.
- Simplified geometry without structural components.
- Idealized material definitions.
- No experimental validation claimed.

The model should be interpreted as a conceptual fast microreactor representation suitable for parametric trend analysis.

---

## Implementation

The project is organized in a modular structure:

- OpenMC Python scripts for reflector sweeps and coolant studies.
- Automatic statepoint saving for each case.
- CSV export of key results (k_eff, uncertainties, leakage fraction).
- Octave scripts for post-processing and generation of thesis-ready plots.
- Explicit uncertainty propagation for reactivity-derived quantities.


---

## Results Summary (High-Level Observations)

- Increasing reflector thickness increases k_eff and reduces leakage, as expected for a fast-spectrum system.
- Reflector worth shows diminishing returns at larger thicknesses.
- Inserted control rod cases consistently produce lower k_eff than withdrawn cases.
- Shutdown margin remains positive across the coolant options evaluated.
- Lithium shows the highest k_eff and shutdown margin under the present assumptions, although differences between some coolants fall within statistical uncertainty.

These observations are specific to the modeling assumptions used here.

---

## Limitations and Future Work

This study does not represent a full reactor design or optimization effort. Possible extensions include:

- Multi-temperature Doppler feedback analysis,
- Burnup and depletion modeling,
- Coupled thermal-hydraulic calculations,
- Mechanical and material compatibility assessment,
- Multi-rod or dynamic reactivity insertion analysis,
- System-level mass and performance trade studies.

These are beyond the scope of the present work but can be built on the current framework.

---

## Purpose of the Project

This repository is intended to document a structured neutronic modeling workflow for a conceptual fast microreactor. The emphasis is on clarity of assumptions, reproducibility, and conservative interpretation of results.

The work should be viewed as a parametric neutronic exploration rather than a validated reactor design.