"""1D comparison of element choices for case (ii) (convex hull, kappa = 0), 200 cells, dt = 1e-3, 200 steps.
Metrics: excess total variation TV(x_h) - 0.5 (= spurious oscillation; the exact TV of a monotone
profile from 0 to 1 is 1 initially and 0.5 for the final step), its maximum over the transient and
its final value, violations of the physical bounds [0,1], final interface width (cells with 0.26<x<0.74),
Newton iterations per step and wall time.
"""
import numpy, time
from mpi4py import MPI
from petsc4py.PETSc import ScalarType
from dolfinx.fem import FunctionSpace, Function, form, assemble_scalar, locate_dofs_topological, dirichletbc
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from dolfinx.mesh import create_interval, locate_entities_boundary
from ufl import TestFunction, dot, FiniteElement, MixedElement, split, Measure, derivative

def magnitude(x): return (x**2)**(1/2)
def ramp(x): return (magnitude(x) + x)/2
def f0_(x): return 0.005 + x*0.02 + ramp(x-0.75)**2 + ramp(0.25-x)**2

VARIANTS = {
    "CG1/CG1/CG1 (C1)":        (("CG", 1), ("CG", 1), 2),
    "CG2/CG2/CG2 (C2, paper)": (("CG", 2), ("CG", 2), 4),
    "DG0/CG1/DG0 (mix)":       (("DG", 0), ("CG", 1), 1),
    "DG1/CG2/DG1 (mix, p=1)":  (("DG", 1), ("CG", 2), 4),
}

def total_variation(S, X):
    """TV of the piecewise polynomial x_h along the interval, including jumps between cells."""
    coords = S.tabulate_dof_coordinates()[:, 0]
    vals = X.x.array
    mesh = S.mesh
    ncells = mesh.topology.index_map(1).size_local
    cells = []
    for c in range(ncells):
        d = S.dofmap.cell_dofs(c)
        order = numpy.argsort(coords[d])
        cells.append((coords[d][order], vals[d][order]))
    cells.sort(key=lambda t: t[0].mean())
    tv = 0.0
    prev_last = None
    for xc, vc in cells:
        tv += numpy.abs(numpy.diff(vc)).sum() if len(vc) > 1 else 0.0
        if prev_last is not None:
            tv += abs(vc[0] - prev_last)
        prev_last = vc[-1]
    return tv

def run(name, spec, nsteps=200, dt=1e-3, ncell=200):
    (sfam, sdeg), (jfam, jdeg), qdeg = spec
    domain = create_interval(MPI.COMM_WORLD, ncell, (0.0, 1.0))
    es = FiniteElement(sfam, domain.ufl_cell(), sdeg)
    ej = FiniteElement(jfam, domain.ufl_cell(), jdeg)
    W = FunctionSpace(domain, MixedElement([es, ej, es]))
    S = FunctionSpace(domain, (sfam, sdeg))
    dx = Measure("dx", metadata={"quadrature_degree": qdeg})
    w = Function(W); (X, j, a) = split(w); q = TestFunction(W)
    X_n = Function(S); X_n.interpolate(lambda x: x[0])
    L = ((f0_(X) - f0_(X_n))/dt + dot(j, j)/2 + a*((X - X_n)/dt + j.dx(0)))*dx
    dL = derivative(L, w, q)
    facets = locate_entities_boundary(domain, 0, lambda x: numpy.full(x.shape[1], True))
    dofs = locate_dofs_topological(W.sub(1), 0, facets)
    bc = [dirichletbc(ScalarType(0), dofs, W.sub(1))]
    solver = NewtonSolver(MPI.COMM_WORLD, NonlinearProblem(dL, w, bc))
    solver.rtol = 1e-10; solver.max_it = 100
    X1 = Function(S)
    tv_hist, its_hist, viol = [], [], 0.0
    t0 = time.perf_counter()
    for i in range(nsteps):
        its, conv = solver.solve(w)
        assert conv
        X1.interpolate(w.sub(0))
        X_n.x.array[:] = X1.x.array
        tv_hist.append(total_variation(S, X1) - (X1.x.array.max() - X1.x.array.min()))
        its_hist.append(its)
        viol = max(viol, -X1.x.array.min(), X1.x.array.max() - 1.0)
    wall = time.perf_counter() - t0
    v = X1.x.array
    width = int(numpy.sum((v > 0.26) & (v < 0.74)))
    print(f"{name:26s} max oscillation (TV - range) = {max(tv_hist):.3f}   final oscillation = {tv_hist[-1]:.2e}   "
          f"bound violation = {max(viol,0):.1e}   final dofs in (0.26,0.74): {width:3d}   "
          f"Newton it/step = {numpy.mean(its_hist):.1f}   wall = {wall:5.1f} s   ndofs = {W.dofmap.index_map.size_global}")

for name, spec in VARIANTS.items():
    run(name, spec)
print("done")
