# Curcumin vs. Doxorubicin Adsorption on Graphene Oxide (DFT)

A first-principles comparison of how strongly two anticancer/anti-inflammatory
drugs — **curcumin** and **doxorubicin** — bind to a graphene oxide (GO)
nanocarrier, computed with DFT (B3LYP/6-31G*, ORCA 6.1.1) and refined across
five levels of correction: geometry optimization, BSSE (counterpoise),
D3BJ dispersion, their combination, and implicit aqueous solvation.

**Live interactive dashboard:** https://sareer555.github.io/curcumin-doxorubicin-GO-DFT/ (public, no sign-in needed) &mdash; full binding-energy chart, molecule/complex visualizations, and methodology notes. A copy also lives at `docs/dashboard.html`.

## Headline result

| | Curcumin · GO | Doxorubicin · GO |
|---|---|---|
| Raw B3LYP | −11.93 kcal/mol | −4.44 kcal/mol |
| BSSE-corrected (CP) | −9.02 kcal/mol | +0.15 kcal/mol |
| D3BJ only | −23.09 kcal/mol | −13.75 kcal/mol |
| **BSSE + D3BJ (most rigorous, gas phase)** | **−20.52 kcal/mol** | **−9.23 kcal/mol** |
| D3BJ + SMD(water) | −17.45 kcal/mol | −13.16 kcal/mol |

At the most rigorous gas-phase level (BSSE + D3BJ), **curcumin binds GO about
2.2× more strongly than doxorubicin**. In implicit water the gap narrows to
1.33×, consistent with solvent competing for surface contact with both drugs.

Full numeric results, reproducible directly from the raw ORCA output files,
are in [`data/results.csv`](data/results.csv).

## Why five correction levels?

- **Raw B3LYP** is the naive number: complex energy minus the two relaxed,
  isolated-monomer energies. It is affected by basis-set superposition error
  (BSSE) and B3LYP's well-known weakness at dispersion (π–π stacking).
- **BSSE (counterpoise)** removes the artificial stabilization that comes
  from each fragment "borrowing" basis functions from its neighbor in the
  complex calculation. Applying it *alone*, without dispersion, is actually
  misleading here — it drives doxorubicin's binding energy to a
  physically implausible **+0.15 kcal/mol** (see `data/results.csv`,
  `BSSE_only` row), because removing the BSSE artifact without adding back
  the real dispersion attraction leaves a false picture of no binding at all.
- **D3BJ (Grimme dispersion, Becke–Johnson damping)** restores the missing
  π-stacking attraction between each drug's aromatic rings and the GO sheet.
- **BSSE + D3BJ** combines both corrections and is treated as the most
  defensible gas-phase estimate in this project.
- **D3BJ + SMD(water)** adds implicit aqueous solvation (single point on the
  gas-phase-optimized geometry, not a solvent-phase re-optimization) to
  approximate physiological conditions. It is reported as a separate,
  complementary refinement — **not stacked with BSSE** — since counterpoise
  correction and continuum solvation are not a standard combination.

A useful internal consistency check: the BSSE magnitude is essentially
identical with and without D3BJ (5.87 / 5.87 kcal/mol for curcumin, 4.93 /
4.93 kcal/mol for doxorubicin), which is expected since BSSE is a pure
basis-set artifact independent of dispersion.

## Repository layout

```
inputs/
  optimization/   geometry-optimization .inp files (drug, GO, and complex)
  bsse/           counterpoise (ghost-atom) single-point .inp files, incl. BSSE+D3BJ
  d3bj/           D3BJ-only single-point .inp files
  solvation/      D3BJ + SMD(water) single-point .inp files
outputs/
  optimization/   converged .out / .xyz for GO, each isolated drug, and each final complex
  bsse/           .out files for every BSSE and BSSE+D3BJ single point
  d3bj/           .out files for every D3BJ-only single point
  solvation/      .out files for every D3BJ+SMD(water) single point
data/
  results.csv     every binding energy at every correction level, machine-readable
scripts/
  extract_energy.py                    pulls FINAL SINGLE POINT ENERGY / geometries out of .out files
  build_bsse_inputs.py                 generates the 8 ghost-atom counterpoise input files per drug
  build_d3bj_and_solvation_inputs.py   generates the D3BJ-only and D3BJ+SMD(water) input files
  compute_binding_energies.py          recomputes every row of results.csv directly from outputs/
  generate_figures.py                  builds every image in figures/ (RDKit 2D + ASE 3D renders)
figures/          2D structures (RDKit) and 3D ball-and-stick renders (ASE) of GO, complexes
docs/
  dashboard.html  the interactive results dashboard (self-contained, no build step)
```

## Reproducing the numbers

All binding energies are recomputed directly from the raw ORCA `.out` files —
nothing in `data/results.csv` is hand-entered:

```bash
python scripts/compute_binding_energies.py
```

Regenerating every figure (2D structures + 3D complex renders):

```bash
pip install rdkit ase
python scripts/generate_figures.py
```

## Method

- **Software:** ORCA 6.1.1
- **Level of theory:** B3LYP/6-31G*, with Grimme D3(BJ) dispersion and
  SMD(water) implicit solvation where noted
- **GO model:** a 28-atom cluster (C₁₆H₁₀O₂ — one epoxide, one hydroxyl group)
  representing a graphene oxide edge site
- **Drugs:** curcumin (C₂₁H₂₀O₆, 47 atoms) and doxorubicin (C₂₇H₂₉NO₁₁, 68 atoms)
- **BSSE method:** Boys–Bernardi counterpoise correction with ghost atoms
  (`:` suffix in ORCA coordinate blocks), 4 single points per drug at the
  frozen complex geometry (monobasis-drug, monobasis-GO, dimerbasis-drug,
  dimerbasis-GO)
- **Binding energy convention:** E(binding) = E(complex) − [E(isolated drug)
  + E(isolated GO)], reported in kcal/mol (1 Hartree = 627.5094740631 kcal/mol)

## Honest caveats

- The isolated-curcumin geometry-optimization job on the source machine did
  not retain a standalone `.xyz`/`.inp` pair; its converged geometry was
  recovered by extracting the last `CARTESIAN COORDINATES (ANGSTROEM)` block
  from `curcumin.out` (see `scripts/extract_energy.py`).
- Both complex optimizations required a mid-run restart from the
  last-saved geometry after unrelated hardware interruptions on the local
  workstation; only the final, fully converged (`ORCA TERMINATED NORMALLY`)
  outputs are kept in `outputs/optimization/`.
- The aqueous curcumin binding energy (−73.0 kJ/mol) is *more* negative than
  a comparable literature value (−34.24 kJ/mol) for a related system — this
  divergence is reported as-is rather than adjusted to match, and is
  attributed to real methodology differences (cluster size, functional
  choice, solvation model) rather than claimed as agreement.
- BSSE and implicit solvation are deliberately not combined into a single
  five-correction stack; see "Why five correction levels?" above.

## Software / tools

ORCA 6.1.1 (DFT engine, run locally) · RDKit (2D structures) · ASE (3D
renders) · Python 3 (analysis, figure generation)
