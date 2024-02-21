from time import perf_counter as time
from sys import argv
import numpy as np
from firedrake import *
parameters["form_compiler"]["quadrature_degree"] = 4


def parprint(*args):
    if COMM_WORLD.rank == 0:
        print("[=]", *args, flush=True)


N = int(argv[1])
load = float(argv[2])
prec = argv[3]
mesh = BoxMesh(3*N, N, N, 1e-2, 0.3e-2, 0.3e-2)
dx = dx(mesh)
ds = ds(mesh)


# Functional setting
V = VectorFunctionSpace(mesh, 'CG', 2)
u = Function(V)
v = TestFunction(V)

# Ramping (HoMoToPy CoNtInUaTiOn)
ramp = Constant(0.0)
Nsteps = 1
ramp_steps = np.linspace(0, 1, Nsteps+1)

solver_parameters = {
    "snes_type": "newtonls",
    "snes_atol": 1e-12,
    "snes_rtol": 1e-6,
    "snes_stol": 0.0,
    "snes_linesearch_type": "basic",
    "ksp_type": "gmres",
    "ksp_atol": 0.0,
    "ksp_rtol": 1e-6,
    "ksp_max_it": 5000,
    "ksp_norm_type": "unpreconditioned",
    "ksp_gmres_restart": 1000,
    "pc_type": prec,
    # "ps_asm_blocks": 8, # Go hard
    "sub_ksp_type": "preonly",
    "sub_pc_type": "ilu",
    "pc_mg_adapt_interp_coarse_space": "gdsw",
    "pc_mg_galerkin": None,
    "pc_mg_levels": 2,
    "mg_levels_pc_type": "asm",
    "snes_error_if_not_converged": True,
    "ksp_converged_reason": None
}

# Variational formulation
f = Identity(3) + grad(u)  # Inverse tensor for inverse problem
j = det(f)
F = variable(inv(f))  # Compute original one to diff
J = det(F)
Cbar = J**(-2/3) * F.T * F

E = 1.0e4
nu = 0.3
mu = Constant(E/(2*(1 + nu)), domain=mesh)
lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu)), domain=mesh)
Ic = tr(Cbar)
psi = (mu / 2) * (Ic - 3) + 0.5 * lmbda * (J-1) * ln(J)
P = diff(psi, F)
rhos = Constant(1e3, domain=mesh)

F_form = inner(j * P, grad(v) * inv(f)) * dx - ramp * \
    rhos * Constant(load, domain=mesh) * v[2] * dx

zero = Constant((0.0, 0.0, 0.0), domain=mesh)
bcs = DirichletBC(V, zero, 1)  # x = 0
problem = NonlinearVariationalProblem(F_form, u, bcs=bcs)

# Set RMs
x = mesh.coordinates
e0 = Constant((1, 0, 0), domain=mesh)
e1 = Constant((0, 1, 0), domain=mesh)
e2 = Constant((0, 0, 1), domain=mesh)
def genVec(vv): interpolate(vv, V)


exps = [e0, e1, e2, cross(e0, x), cross(e1, x), cross(e2, x)]
vecs = [Function(V) for vv in exps]
for i, e in enumerate(exps):
    vecs[i].interpolate(e)
nn = VectorSpaceBasis(vecs)
nn.orthonormalize()

solver = NonlinearVariationalSolver(
    problem, solver_parameters=solver_parameters, near_nullspace=nn)
l_its = []
init_time = time()
for r in ramp_steps:
    if r < 1e-12:
        continue
    parprint("Solving step ramp={}".format(r))
    ramp.assign(r)
    solver.solve()
    ITS = solver.snes.getLinearSolveIterations()
    l_its.append(ITS)

vecs = []
parprint("Dofs=", V.dim())
parprint("Avg tot linear its:", sum(l_its)/len(l_its))
parprint("CPU time:", time() - init_time)
