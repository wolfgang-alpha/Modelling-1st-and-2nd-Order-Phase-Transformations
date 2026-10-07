"""Numerical check of Proposition 1 (discrete mass conservation and energy-dissipation identity)
for the 1D binary problems of the paper, with the two element choices used in the notebooks:
  C2  : equal-order CG2 for (x, j, mu_bar), gradient energy c*|grad x|^2
  mix : DG0 / CG1 / DG0, gradient energy as facet jump penalty c*sum_F [x]^2
Runs inside dolfinx/lab:v0.7.2 (same container as the published notebooks).

Identity checked at every step n:
  F[x^{n+1}] + dt*P[j^{n+1}] + B_n - F[x^n] = 0
with F[x] = int f(x) dx + reg(x),  P[j] = int |j|^2 dx  (k/D = 1 in the scaled units of the notebooks),
B_n = int ( f(x^n) - f(x^{n+1}) - f'(x^{n+1}) (x^n - x^{n+1}) ) dx + reg(x^{n+1} - x^n).
"""
import numpy, sys
from mpi4py import MPI
from petsc4py.PETSc import ScalarType
from dolfinx.fem import FunctionSpace, Function, form, assemble_scalar, locate_dofs_topological, dirichletbc
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.mesh import create_interval, locate_entities_boundary
from ufl import (TestFunction, dot, grad, jump, FacetNormal, FiniteElement, MixedElement, split,
                 Measure, derivative, replace)
from ufl.algorithms import expand_derivatives

# --- free energies, verbatim from the notebooks --------------------------------
def magnitude(x): return (x**2)**(1/2)
def ramp(x): return (magnitude(x) + x)/2
def step(x, a): return (-2*a + 2*x)*((-a + x)**2)**(-0.5)/2 + (1.0*a - 1.0*x)*((-a + x)**2)**(-0.5)/2 + 1/2
def h(x, a, b):
    x = (x-a)/(b-a); return 3*x**2 - 2*x**3
def smooth_step(x, a, b): return (step(x, a) - step(x, b))*h(x, a, b) + step(x, b)
def alpha(x): return 1 - smooth_step(x, 0.25, 0.75)
def beta(x): return 1 - alpha(x)
def f0_(x): return 0.005 + x*0.02 + ramp(x-0.75)**2 + ramp(0.25-x)**2
def f_0(x): return f0_(x) + alpha(x)*beta(x)/30


def run(variant, f, c_reg, label, nsteps=200, dt=1e-3):
    domain = create_interval(MPI.COMM_WORLD, 200, (0.0, 1.0))
    n = FacetNormal(domain)
    dS = Measure("dS")
    if variant == "C2":
        e = FiniteElement("CG", domain.ufl_cell(), 2)
        W = FunctionSpace(domain, MixedElement([e, e, e]))
        S = FunctionSpace(domain, ("CG", 2))
        dx = Measure("dx", metadata={"quadrature_degree": 4})
        reg = lambda X: c_reg*dot(grad(X), grad(X))*dx
    else:  # mix
        eDG = FiniteElement("DG", domain.ufl_cell(), 0)
        eCG = FiniteElement("CG", domain.ufl_cell(), 1)
        W = FunctionSpace(domain, MixedElement([eDG, eCG, eDG]))
        S = FunctionSpace(domain, ("DG", 0))
        dx = Measure("dx", metadata={"quadrature_degree": 1})
        reg = lambda X: c_reg*dot(jump(X, n), jump(X, n))*dS

    w = Function(W)
    (X, j, a) = split(w)
    q = TestFunction(W)
    X_n = Function(S)
    X_n.interpolate(lambda x: x[0])

    # Lagrangian/Rayleighian exactly as in the notebooks (incremental form)
    L = ((f(X) - f(X_n))/dt + dot(j, j)/2 + a*((X - X_n)/dt + j.dx(0)))*dx
    if c_reg > 0:
        L += (1.0/dt)*reg(X)
    dL = derivative(L, w, q)

    facets = locate_entities_boundary(domain, 0, lambda x: numpy.full(x.shape[1], True))
    dofs = locate_dofs_topological(W.sub(1), 0, facets)
    bc = [dirichletbc(ScalarType(0), dofs, W.sub(1))]
    problem = NonlinearProblem(dL, w, bc)
    solver = NewtonSolver(MPI.COMM_WORLD, problem)
    solver.rtol = 1e-12
    solver.atol = 1e-13
    solver.max_it = 100

    X1 = Function(S)          # new state (collapsed copy for the functionals)
    j1 = Function(FunctionSpace(domain, W.sub(1).element.basix_element.family if False else ("CG", 2 if variant == "C2" else 1)))

    def F_of(Xf):
        val = assemble_scalar(form(f(Xf)*dx))
        if c_reg > 0:
            val += assemble_scalar(form(reg(Xf)))
        return val

    mass0 = assemble_scalar(form(X_n*dx))
    worst_identity, worst_mass, minB, Fprev = 0.0, 0.0, numpy.inf, F_of(X_n)
    Fhist = [Fprev]
    for i in range(nsteps):
        its, conv = solver.solve(w)
        assert conv
        X1.interpolate(w.sub(0))
        j1.interpolate(w.sub(1))
        Fnew = F_of(X1)
        P = assemble_scalar(form(dot(j1, j1)*dx))
        # Bregman term, directional derivative of int f(X1) dx in the direction (X_n - X1)
        dF = expand_derivatives(derivative(f(X1)*dx, X1, X_n - X1))
        B = assemble_scalar(form((f(X_n) - f(X1))*dx)) - assemble_scalar(form(dF))
        if c_reg > 0:
            Xd = Function(S); Xd.x.array[:] = X1.x.array - X_n.x.array
            B += assemble_scalar(form(reg(Xd)))
        ident = Fnew + dt*P + B - Fprev
        worst_identity = max(worst_identity, abs(ident))
        minB = min(minB, B)
        mass = assemble_scalar(form(X1*dx))
        worst_mass = max(worst_mass, abs(mass - mass0))
        X_n.x.array[:] = X1.x.array
        Fprev = Fnew
        Fhist.append(Fnew)
    Fhist = numpy.array(Fhist)
    print(f"{label:34s} max|identity| = {worst_identity:.2e}   max|mass drift| = {worst_mass:.2e}   "
          f"min B_n = {minB:+.2e}   F monotone decreasing: {bool(numpy.all(numpy.diff(Fhist) <= 1e-15))}   "
          f"F0={Fhist[0]:.6f} F_end={Fhist[-1]:.6f}")


for variant in ["C2", "mix"]:
    c_CH = 5e-5 if variant == "C2" else 0.005
    run(variant, f_0, c_CH, f"[{variant}] Cahn-Hilliard (double well)")
    run(variant, f0_, 0.0, f"[{variant}] case (ii) convex hull, kappa=0")
    run(variant, f0_, c_CH, f"[{variant}] case (i) convex hull + reg")
print("done")
