# A Novel Methodology for Modelling First and Second Order Phase Transformations

## Quickstart with Docker: (Linux/Unix, Windows with WSL and docker integration enabled)

```
git clone https://github.com/wolfgang-alpha/Modelling-1st-and-2nd-Order-Phase-Transformations.git
cd Modelling-1st-and-2nd-Order-Phase-Transformations
bash batch.sh 
```

Open the displayed url "http://127.0.0.1:8888/lab?token=..." in a web-browser.

## Update (July 2026): clean variational construction

The Lagrangian cells now generate the discrete variational forms by the corrected
thermodynamic-extremal-principle construction: the free-energy rate enters as the
directional derivative `derivative(F, X_s, X_t)` (chain rule at the frozen state `X_s`,
with `X_t` the discrete rate), the variation `derivative(L, w, q)` acts on the rates and
multipliers only, and `replace(..., {X_s: X})` substitutes the frozen state *after* the
variation (backward Euler). This removes the ambiguity of varying a total time difference
of the state. The resulting forms are numerically identical to the previous release —
see `verify_construction.py` (max deviation ~1e-15 inside the `dolfinx/lab:v0.7.2`
container) — so all stored results remain valid; only the derivation is now unambiguous.
