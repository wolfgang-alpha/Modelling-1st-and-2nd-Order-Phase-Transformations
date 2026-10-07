# Interface Energy and Phase Transformations: Cahn-Hilliard and CALPHAD-based Models in Ternary Substitutional Alloys

Simulation code for the paper

> W. Flachberger, T. Antretter, S. Gaddikere-Nagaraja, S. Leitner, M. Petersmann, J. Svoboda,
> *Interface Energy and Phase Transformations: A Comparative Analysis of Cahn-Hilliard and
> CALPHAD-based Models in Ternary Substitutional Alloys*, arXiv:2411.16430.

A single finite element scheme, derived from the thermodynamic extremal principle (stationarity of
the Rayleighian: free-energy rate + 1/2 dissipation + mass-conservation constraint), treats the
Cahn-Hilliard model, the convex-hull (CALPHAD) free energy with regularization (case i) and the
convex hull without regularization (case ii, sharp interface, degenerate parabolic / Stefan type)
with the same code. Fluxes and chemical affinities are primary unknowns, so no derivative of the
affinity is needed. The time-discrete scheme conserves mass exactly and satisfies a discrete
energy-dissipation inequality for every convex molar free energy (Proposition 1 of the paper).

## Contents

| File | Purpose |
|---|---|
| `ch-limit-CG2.ipynb` | Equal-order continuous P2/Q2 elements for mole fraction, flux and affinity. **Produces all figures of the paper** (1D: Figs. 2-6, 2D: Figs. 7-12). Incremental-functional construction. |
| `ch-limit-RTDG.ipynb` | Raviart-Thomas / discontinuous pair (DG0 / CG1 / DG0 in 1D, DG0 / RTCF1 / DG0 in 2D), gradient energy as facet jump penalty. Cross-check of the sharp-interface case (Fig. 9b of the paper). Uses the clean TEP construction (variation w.r.t. rates and multipliers at the frozen state, then backward Euler). |
| `verify_construction.py` | Checks that the two constructions assemble identical residuals (to ~1e-15 relative) for all four element layouts of the two notebooks. |
| `code/check_identity.py` | Numerical check of Proposition 1: exact mass conservation and the identity F[x^{n+1}] + dt P[j^{n+1}] + B_n = F[x^n] (1D, both element choices). |
| `code/compare_elements_1d.py` | 1D comparison of four element choices on the sharp-interface case (Table 1 of the paper). |
| `batch.sh` | Launches JupyterLab in the `dolfinx/lab:v0.7.2` container with the repository mounted. |

## Quickstart with Docker (Linux/Unix, Windows with WSL and docker integration enabled)

```
git clone https://github.com/wolfgang-alpha/Modelling-1st-and-2nd-Order-Phase-Transformations.git
cd Modelling-1st-and-2nd-Order-Phase-Transformations
bash batch.sh
```

Open the displayed url `http://127.0.0.1:8888/lab?token=...` in a web browser. The scripts run
without JupyterLab:

```
docker run --rm --entrypoint python3 -v "$(pwd)":/work -w /work dolfinx/lab:v0.7.2 verify_construction.py
docker run --rm --entrypoint python3 -v "$(pwd)/code":/work -w /work dolfinx/lab:v0.7.2 check_identity.py
docker run --rm --entrypoint python3 -v "$(pwd)/code":/work -w /work dolfinx/lab:v0.7.2 compare_elements_1d.py
```

All computations use dolfinx 0.7.2 (the API of later releases differs).

## Timings (single thread, 128 x 128 quadrilaterals, 100 steps of dt = 1e-3)

| Case | CG2 (paper figures) | RT/DG |
|---|---|---|
| Cahn-Hilliard (Fig. 8) | 17 min | 32 min |
| case (ii), convex hull, kappa = 0 (Fig. 9 / 9b) | 109 min | 4 min |
| case (i), convex hull + regularization (Fig. 10) | 100 min | 30 min |
| ternary with vacancies (Figs. 11, 12) | 723 min | 14 min |

The cost is governed by the convergence of the damped Newton iteration in the convex-hull cases,
not by the regularization.

## Results of the checks (dolfinx 0.7.2, 2026-10-07)

`verify_construction.py`: max |dL_old - dL_new| = 3.6e-15, 1.1e-15 (1D, DG0/CG1/DG0 and CG2),
2.3e-13 and 1.1e-13 (2D ternary, RT and Q2; residual scale 1e3): the clean construction and the
incremental construction coincide.

`code/check_identity.py` (1D, 200 cells, dt = 1e-3, 200 steps, Newton rtol 1e-12):

```
[C2]  Cahn-Hilliard   max|identity| = 4.0e-16  max|mass drift| = 9.4e-16  min B_n = +5.7e-10  F monotone: True
[C2]  case (ii)       max|identity| = 3.5e-17  max|mass drift| = 7.2e-16  min B_n = +3.8e-10  F monotone: True
[C2]  case (i)        max|identity| = 3.5e-17  max|mass drift| = 8.9e-16  min B_n = +9.1e-10  F monotone: True
[mix] Cahn-Hilliard   max|identity| = 1.1e-16  max|mass drift| = 5.6e-16  min B_n = +5.6e-10  F monotone: True
[mix] case (ii)       max|identity| = 1.9e-17  max|mass drift| = 5.6e-16  min B_n = +3.8e-10  F monotone: True
[mix] case (i)        max|identity| = 2.8e-17  max|mass drift| = 5.6e-16  min B_n = +7.5e-10  F monotone: True
```

`code/compare_elements_1d.py` (case (ii); oscillation = total variation minus range of x_h):

```
CG1/CG1/CG1              max oscillation 0.357   final 4.1e-05   dofs in (0.26,0.74): 1   Newton it/step 1.3
CG2/CG2/CG2 (paper)      max oscillation 0.499   final 7.0e-03   dofs in (0.26,0.74): 5   Newton it/step 1.8
DG0/CG1/DG0 (RT/DG)      max oscillation 11.904  final 0.0       dofs in (0.26,0.74): 0   Newton it/step 1.9
DG1/CG2/DG1              max oscillation 1.139   final 3.8e-02   dofs in (0.26,0.74): 2   Newton it/step 1.5
```

## History

- **October 2026 (v3 of the paper):** `ch-limit-CG2.ipynb` added (the notebook behind the figures;
  the legend of Fig. 5 and the plotting of Fig. 6 were fixed and the 1D cells re-executed, the 2D
  outputs are unchanged), the RT/DG notebook renamed from `ch-limit.ipynb` to `ch-limit-RTDG.ipynb`,
  `verify_construction.py` extended to the equal-order layouts, `code/` added.
- **July 2026:** clean variational construction in the RT/DG notebook. The free-energy rate enters
  as the directional derivative `derivative(F, X_s, X_t)` (chain rule at the frozen state `X_s`,
  with `X_t` the discrete rate), the variation `derivative(L, w, q)` acts on the rates and
  multipliers only, and `replace(..., {X_s: X})` substitutes the frozen state *after* the
  variation (backward Euler). This removes the ambiguity of varying a total time difference of
  the state; the resulting forms are numerically identical to the incremental construction.
- **November 2024:** initial release with arXiv:2411.16430v1.
