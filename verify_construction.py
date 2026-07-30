"""Verify that the clean TEP construction used in ch-limit.ipynb assembles the
identical residual to the originally published Lagrangian construction.

Run inside the repo's container (dolfinx 0.7.2):

    docker run --rm --entrypoint python3 -v "$(pwd)":/work -w /work \
        dolfinx/lab:v0.7.2 verify_construction.py
"""

import numpy
from mpi4py import MPI
from dolfinx.fem import FunctionSpace, Function, form, assemble_vector
from dolfinx.mesh import create_interval, create_unit_square, CellType
from ufl import (TestFunction, jump, FacetNormal, dot, div, FiniteElement,
                 MixedElement, split, Measure, derivative, replace)
from ufl.algorithms import expand_derivatives

# --- free energies (cell 7 of the notebook, verbatim) --------------------------
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
    print(f"{label:28s} max|dL_old - dL_new| = {diff:.3e}  (|dL| ~ {scale:.3e})")
    assert diff <= 1e-12 * max(scale, 1.0)


# --- 1D binary (cells 18/19 layout) --------------------------------------------
domain = create_interval(MPI.COMM_WORLD, 50, (0.0, 1.0))
n = FacetNormal(domain)
DG = FunctionSpace(domain, ("DG", 0))
e_CG = FiniteElement("CG", domain.ufl_cell(), 1)
e_DG = FiniteElement("DG", domain.ufl_cell(), 0)
W = FunctionSpace(domain, MixedElement([e_DG, e_CG, e_DG]))
w = Function(W)
(X, j, a) = split(w)
q = TestFunction(W)
dt = 1e-3
dx = Measure("dx", metadata={"quadrature_degree": 1})
dS = Measure("dS")

X_n = Function(DG)
X_n.interpolate(lambda x: x[0])
rng = numpy.random.default_rng(1)
w.x.array[:] = 0.4 + 0.1*rng.standard_normal(w.x.array.size)

# published construction
L_old = ((f_0(X)-f_0(X_n))/dt + dot(j, j)/2 + a*((X-X_n)/dt + j.dx(0)))*dx
L_old += (0.005/dt)*dot(jump(X, n), jump(X, n))*dS
dL_old = derivative(L_old, w, q)

# clean TEP construction
X_s = Function(DG)
X_s.x.array[:] = rng.random(X_s.x.array.size)  # value irrelevant: replaced
X_t = (X-X_n)/dt
F = f_0(X_s)*dx + 0.005*dot(jump(X_s, n), jump(X_s, n))*dS
dF = expand_derivatives(derivative(F, X_s, X_t))
L_new = dF + (dot(j, j)/2 + a*(X_t + j.dx(0)))*dx
dL_new = derivative(L_new, w, q)
dL_new = replace(expand_derivatives(dL_new), {X_s: X})

compare(dL_old, dL_new, "1D binary (CH + jump)")

# --- 2D ternary with vacancies (cells 29/37 layout) ----------------------------
domain = create_unit_square(MPI.COMM_WORLD, 8, 8, CellType.quadrilateral)
DG = FunctionSpace(domain, ("DG", 0))
e_DG = FiniteElement("DG", domain.ufl_cell(), 0)
e_RT = FiniteElement("RTCF", domain.ufl_cell(), 1)
W = FunctionSpace(domain, MixedElement([e_DG, e_DG, e_RT, e_RT, e_DG, e_DG, e_DG]))
w = Function(W)
(X0, X1, j0, j1, a0, a1, phi) = split(w)
q = TestFunction(W)
dx = Measure("dx", metadata={"quadrature_degree": 1})

X0_ = Function(DG)
X1_ = Function(DG)
X0_.x.array[:] += 1e-3
X1_.interpolate(lambda x: x[0])
w.x.array[:] = 0.4 + 0.1*rng.standard_normal(w.x.array.size)

def f0(X0, X1):
    return f0_(X1) + 100*(X0-1e-3)**2

# published construction
L_old = (f0(X0, X1) - f0(X0_, X1_))/dt
L_old += dot(j1, j1)/2 + dot(-j0-j1, -j0-j1) + phi**2
L_old += a1*((X1-X1_)/dt + div(j1) + phi*X1_) + a0*((X0-X0_)/dt + div(j0) - phi*(1-X0_))
L_old = L_old*dx
dL_old = derivative(L_old, w, q)

# clean TEP construction
X0_s, X1_s = Function(DG), Function(DG)
X0_s.x.array[:] = rng.random(X0_s.x.array.size)
X1_s.x.array[:] = rng.random(X1_s.x.array.size)
X0_t, X1_t = (X0-X0_)/dt, (X1-X1_)/dt
F = f0(X0_s, X1_s)*dx
dF = expand_derivatives(derivative(F, X0_s, X0_t) + derivative(F, X1_s, X1_t))
L_new = dot(j1, j1)/2 + dot(-j0-j1, -j0-j1) + phi**2
L_new += a1*(X1_t + div(j1) + phi*X1_) + a0*(X0_t + div(j0) - phi*(1-X0_))
L_new = dF + L_new*dx
dL_new = derivative(L_new, w, q)
dL_new = replace(expand_derivatives(dL_new), {X0_s: X0, X1_s: X1})

compare(dL_old, dL_new, "2D ternary (vacancies)")

print("dolfinx 0.7.2: clean construction == published forms.")
