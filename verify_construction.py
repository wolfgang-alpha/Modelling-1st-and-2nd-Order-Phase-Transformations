"""Verify that the clean TEP construction (variation with respect to the rates and
multipliers at the frozen state, then backward Euler; the direction of the free-energy
variation is a placeholder Function that is replaced by the discrete rate afterwards) assembles the identical residual
to the incremental-functional construction, for all element layouts used in the repository:

  1. 1D binary, DG0 / CG1 / DG0 with facet jump penalty      (ch-limit-RTDG.ipynb, 1D cells)
  2. 1D binary, CG2 / CG2 / CG2 with gradient energy          (ch-limit-CG2.ipynb, 1D cells; figures of the paper)
  3. 2D ternary with vacancies, DG0 / DG0 / RTCF1 / RTCF1 / DG0 / DG0 / DG0   (ch-limit-RTDG.ipynb)
  4. 2D ternary with vacancies, Q2 / Q2 / [Q2]^2 / [Q2]^2 / Q2 / Q2 / Q2        (ch-limit-CG2.ipynb)

Run inside the repo's container (dolfinx 0.7.2):

    docker run --rm --entrypoint python3 -v "$(pwd)":/work -w /work \
        dolfinx/lab:v0.7.2 verify_construction.py
"""

import numpy
from mpi4py import MPI
from dolfinx.fem import FunctionSpace, Function, form, assemble_vector
from dolfinx.mesh import create_interval, create_unit_square, CellType
from ufl import (TestFunction, jump, FacetNormal, dot, div, grad, FiniteElement, VectorElement,
                 MixedElement, split, Measure, derivative, replace)
from ufl.algorithms import expand_derivatives

# --- free energies (verbatim from the notebooks) ------------------------------
def magnitude(x):
    return (x**2)**(1/2)

def ramp(x):
    return (magnitude(x) + x)/2

def step(x, a):
    return (-2*a + 2*x)*((-a + x)**2)**(-0.5)/2 + (1.0*a - 1.0*x)*((-a + x)**2)**(-0.5)/2 + 1/2

def h(x, a, b):
    x = (x-a) / (b-a)
    return 3*x**2 - 2*x**3

def smooth_step(x, a, b):
    return (step(x, a) - step(x, b))*h(x, a, b) + step(x, b)

def alpha(x):
    return 1 - smooth_step(x, 0.25, 0.75)

def beta(x):
    return 1 - alpha(x)

def f0_(x):
    return 0.005 + x*0.02 + ramp(x-0.75)**2 + ramp(0.25-x)**2

def f_0(x):
    return f0_(x) + alpha(x)*beta(x)/30


def compare(b_old_form, b_new_form, label):
    b_old = assemble_vector(form(b_old_form)).array
    b_new = assemble_vector(form(b_new_form)).array
    diff = numpy.abs(b_old - b_new).max()
    scale = numpy.abs(b_old).max()
    print(f"{label:44s} max|dL_old - dL_new| = {diff:.3e}  (|dL| ~ {scale:.3e})")
    assert diff <= 1e-12 * max(scale, 1.0)


rng = numpy.random.default_rng(1)
dt = 1e-3

# --- 1. and 2.: 1D binary -------------------------------------------------------
for layout in ["DG0/CG1/DG0 + jump penalty", "CG2/CG2/CG2 + gradient energy"]:
    domain = create_interval(MPI.COMM_WORLD, 50, (0.0, 1.0))
    n = FacetNormal(domain)
    dS = Measure("dS")
    if layout.startswith("DG0"):
        S = FunctionSpace(domain, ("DG", 0))
        e_s = FiniteElement("DG", domain.ufl_cell(), 0)
        e_j = FiniteElement("CG", domain.ufl_cell(), 1)
        dx = Measure("dx", metadata={"quadrature_degree": 1})
        reg = lambda X: 0.005*dot(jump(X, n), jump(X, n))*dS
    else:
        S = FunctionSpace(domain, ("CG", 2))
        e_s = FiniteElement("CG", domain.ufl_cell(), 2)
        e_j = e_s
        dx = Measure("dx", metadata={"quadrature_degree": 4})
        reg = lambda X: 5e-5*dot(grad(X), grad(X))*dx
    W = FunctionSpace(domain, MixedElement([e_s, e_j, e_s]))
    w = Function(W)
    (X, j, a) = split(w)
    q = TestFunction(W)
    X_n = Function(S)
    X_n.interpolate(lambda x: x[0])
    w.x.array[:] = 0.4 + 0.1*rng.standard_normal(w.x.array.size)

    # incremental construction (the one that produced the published figures)
    L_old = ((f_0(X)-f_0(X_n))/dt + dot(j, j)/2 + a*((X-X_n)/dt + j.dx(0)))*dx
    L_old += (1.0/dt)*reg(X)
    dL_old = derivative(L_old, w, q)

    # clean TEP construction
    X_s = Function(S)
    X_s.x.array[:] = rng.random(X_s.x.array.size)  # value irrelevant: replaced
    X_t = (X-X_n)/dt
    v = Function(S)                                    # placeholder direction, replaced by the rate
    F = f_0(X_s)*dx + reg(X_s)
    dF = replace(expand_derivatives(derivative(F, X_s, v)), {v: X_t})
    L_new = dF + (dot(j, j)/2 + a*(X_t + j.dx(0)))*dx
    dL_new = derivative(L_new, w, q)
    dL_new = replace(expand_derivatives(dL_new), {X_s: X})

    compare(dL_old, dL_new, f"1D binary, {layout}")

# --- 3. and 4.: 2D ternary with vacancies --------------------------------------
for layout in ["DG0 / RTCF1 (Raviart-Thomas)", "Q2 / [Q2]^2 (equal order)"]:
    domain = create_unit_square(MPI.COMM_WORLD, 8, 8, CellType.quadrilateral)
    if layout.startswith("DG0"):
        S = FunctionSpace(domain, ("DG", 0))
        e_s = FiniteElement("DG", domain.ufl_cell(), 0)
        e_j = FiniteElement("RTCF", domain.ufl_cell(), 1)
        dx = Measure("dx", metadata={"quadrature_degree": 1})
    else:
        S = FunctionSpace(domain, ("CG", 2))
        e_s = FiniteElement("CG", domain.ufl_cell(), 2)
        e_j = VectorElement("CG", domain.ufl_cell(), 2)
        dx = Measure("dx", metadata={"quadrature_degree": 4})
    W = FunctionSpace(domain, MixedElement([e_s, e_s, e_j, e_j, e_s, e_s, e_s]))
    w = Function(W)
    (X0, X1, j0, j1, a0, a1, phi) = split(w)
    q = TestFunction(W)

    X0_ = Function(S)
    X1_ = Function(S)
    X0_.x.array[:] += 1e-3
    X1_.interpolate(lambda x: x[0])
    w.x.array[:] = 0.4 + 0.1*rng.standard_normal(w.x.array.size)

    def f0(X0, X1):
        return f0_(X1) + 100*(X0-1e-3)**2

    # incremental construction
    L_old = (f0(X0, X1) - f0(X0_, X1_))/dt
    L_old += dot(j1, j1)/2 + dot(-j0-j1, -j0-j1) + phi**2
    L_old += a1*((X1-X1_)/dt + div(j1) + phi*X1_) + a0*((X0-X0_)/dt + div(j0) - phi*(1-X0_))
    L_old = L_old*dx
    dL_old = derivative(L_old, w, q)

    # clean TEP construction
    X0_s, X1_s = Function(S), Function(S)
    X0_s.x.array[:] = rng.random(X0_s.x.array.size)
    X1_s.x.array[:] = rng.random(X1_s.x.array.size)
    X0_t, X1_t = (X0-X0_)/dt, (X1-X1_)/dt
    v0, v1 = Function(S), Function(S)                 # placeholder directions
    F = f0(X0_s, X1_s)*dx
    dF = replace(expand_derivatives(derivative(F, X0_s, v0) + derivative(F, X1_s, v1)), {v0: X0_t, v1: X1_t})
    L_new = dot(j1, j1)/2 + dot(-j0-j1, -j0-j1) + phi**2
    L_new += a1*(X1_t + div(j1) + phi*X1_) + a0*(X0_t + div(j0) - phi*(1-X0_))
    L_new = dF + L_new*dx
    dL_new = derivative(L_new, w, q)
    dL_new = replace(expand_derivatives(dL_new), {X0_s: X0, X1_s: X1})

    compare(dL_old, dL_new, f"2D ternary (vacancies), {layout}")

print("dolfinx 0.7.2: clean construction == incremental construction for all four layouts.")
