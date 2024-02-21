import prestress.prestressSolver as PS
from prestress.Fibers import generateFibers, generateAnalyticFibers
from prestress.Tools import slabGeometry,  prolateGeometry, getXDMFFile, parprint
from prestress.Enums import InverseFormulationType
from petsc4py import PETSc
import argparse
from time import perf_counter as time
parser = argparse.ArgumentParser()

parser.add_argument('--solver-type', type=str, required=True,
                    help="Solver type [inverse|sellier]")
parser.add_argument('--nx', type=int, required=False, default=12,
                    help="Number of elements of slab in x direction. Default: 12")
parser.add_argument('--ny', type=int, required=False, default=4,
                    help="Number of elements of slab in x direction. Default: 4")
parser.add_argument('--nz', type=int, required=False, default=4,
                    help="Number of elements of slab in x direction. Default: 4")
parser.add_argument('--fe-degree', type=int, required=False,
                    default=1, help="FEM degree of displacement. Default: 1")
parser.add_argument('--output-file', type=str, required=False,
                    help="If set, .xdmf is appended and results are stored there")
parser.add_argument('--backward-steps', type=int, required=False,
                    default=10, help="Number of ramp steps for inverse problem")
parser.add_argument('--forward-steps', type=int, required=False, default=0,
                    help="Number of ramp steps for forward problem, 0 to skip. Default: 0")
parser.add_argument('--adaptive', required=False, dest='adaptive',
                    action='store_true', help="Set adaptive step choice. Default: False")
parser.add_argument('--z-load-scale', type=float, required=False, default=-9.8,
                    help="Z-axis volume load, represents gravity if set to -9.8. Default: -9.8")
parser.add_argument('--z-surf-scale', type=float, required=False,
                    default=0.0, help="Z-axis surface load. Default: 0.0")
parser.add_argument('--sellier-armijo', required=False, dest='armijo', action='store_true',
                    help="Flag to do Armijo line search in Sellier. Default: False")
parser.add_argument('--sellier-relaxation', type=float, required=False,
                    default=1.0, help="Relaxation parameter for Sellier. Default: 1.0")
parser.add_argument('--sellier-atol-inc', type=float, required=False, default=1e-6,
                    help="Absolute step tolerance for Sellier increments. Default: 1e-14")
parser.add_argument('--sellier-rtol-inc', type=float, required=False, default=1e-6,
                    help="Relative step tolerance for Sellier increments. Default: 1e-6")
parser.add_argument('--sellier-max-it', type=int, required=False,
                    default=1000, help="Maximum Sellier iterations admitted. Default: 1000")
parser.add_argument('--sellier-inner-steps', type=int, required=False,
                    default=1, help="Maximum Sellier direct solver ramp steps. Default: 1")
parser.add_argument('--sellier-aa-m', type=int, required=False,
                    default=0, help="Acceleration depth for Sellier. Default: 0")
parser.set_defaults(armijo=False, adaptive=False)
args = parser.parse_args()


output_filename = args.output_file
solver_type = args.solver_type  # direct | sellier
assert solver_type in (
    "inverse", "sellier"), "Wrong solver type {}. See --help for further details".format(solver_type)
if solver_type == "inverse":
    solver_type = InverseFormulationType.direct
elif solver_type == "sellier":
    solver_type = InverseFormulationType.sellier
n_times_back = args.backward_steps  # Number of pseudo-timesteps
n_times_forward = args.forward_steps  # Number of pseudo-timesteps
volumetric_load = args.z_load_scale  # Gravity
neumann_load = args.z_surf_scale  # Scale of vertical surface traction
u_deg = args.fe_degree
adaptive = args.adaptive

sellier_params = {"alpha": args.sellier_relaxation,
                  "atol_inc": args.sellier_atol_inc,
                  "rtol_inc": args.sellier_rtol_inc,
                  "max_it": args.sellier_max_it,
                  "inner_steps": args.sellier_inner_steps,
                  "do_armijo": args.armijo,
                  "aa_m": args.sellier_aa_m}  # aa_m = 1 works best

# TODO: Read PETSc options from external file
PETSc.Options().setValue("-ksp_type", "gmres")
PETSc.Options().setValue("-ksp_norm_type", "unpreconditioned")
PETSc.Options().setValue("-ksp_max_it", 1000)
PETSc.Options().setValue("-ksp_atol", 0.0)
PETSc.Options().setValue("-ksp_rtol", 1e-6)
PETSc.Options().setValue("-pc_type", "lu")
PETSc.Options().setValue("-pc_factor_mat_solver_type", "mumps")
PETSc.Options().setValue("-ksp_gmres_restart", 1000)
PETSc.Options().setValue("-snes_max_it", 1000)
PETSc.Options().setValue("-snes_atol", 1e-14)
PETSc.Options().setValue("-snes_rtol", 1e-6)
PETSc.Options().setValue("-snes_stol", 0.0)
PETSc.Options().setValue("-snes_type", "newtonls")
PETSc.Options().setValue("-snes_linesearch_type", "basic")
PETSc.Options().setValue("-snes_error_if_not_converged", "true")
#PETSc.Options().setValue("-snes_converged_reason", None)
#PETSc.Options().setValue("-ksp_converged_reason", None)
#PETSc.Options().setValue("-ksp_monitor", None)

ENDO = EPI = BASE = 0
nx = args.nx
ny = args.ny
nz = args.nz

mesh, markers, X0, Y0, Z0, XL, YL, ZL = slabGeometry(nx, ny, nz)

# UFL doesn't understand square brackets, but a one-element tuple isn't read as a list,
# so a=(4), for i in a: print(i) gives an error.
dirichlet_markers = [X0]
neumann_markers = (XL)
if output_filename:
    outmesh = getXDMFFile(output_filename, mesh=True).write(mesh)
time_0 = time()
dofs, sell_its, nl_its, l_its = PS.prestress_solver_slab(mesh, markers, dirichlet_markers, neumann_markers, u_deg, n_times_back, n_times_forward,
                                                         neumann_load, volumetric_load, getXDMFFile(output_filename), solver_type, sellier_params=sellier_params, adaptive=adaptive)
time_f = time()

parprint("Dofs:", dofs)
parprint("Avg sell its:", sell_its)
parprint("Avg nonlinear its:", nl_its)
parprint("Tot linear its:", l_its)
parprint("Avg linear its:", l_its / max(nl_its, 1))
parprint("Solution time:", time_f - time_0)
