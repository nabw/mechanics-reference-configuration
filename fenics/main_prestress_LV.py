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
parser.add_argument('--heart-geometry-name', type=str, required=False,
                    default="h4_v2_ASCII", help="Geometry name [4mm|h4_v2_ASCII]. Default: h4_v2_ASCII")
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
parser.add_argument('--endo-pressure', type=float, required=False,
                    default=1e4, help="Given endocardial pressure. Default: 1e4")
parser.add_argument('--active-peak', type=float, required=False,
                    default=1e4, help="Given peak active force. Default: 1e4")
parser.add_argument('--z-load-scale', type=float, required=False, default=-9.8,
                    help="Z-axis volume load, represents gravity if set to -9.8. Default: -9.8")
parser.add_argument('--z-surf-scale', type=float, required=False,
                    default=0.0, help="Z-axis surface load. Default: 0.0")
parser.add_argument('--sellier-relaxation', type=float, required=False,
                    default=1.0, help="Relaxation parameter for Sellier. Default: 1.0")
parser.add_argument('--sellier-armijo', required=False, dest='armijo', action='store_true',
                    help="Flag to do Armijo line search in Sellier. Default: False")
parser.add_argument('--sellier-atol', type=float, required=False, default=1e-10,
                    help="Absolute step tolerance for Sellier. Default: 1e-10")
parser.add_argument('--sellier-rtol', type=float, required=False,
                    default=1e-6, help="Relative step tolerance for Sellier. Default: 1e-6")
parser.add_argument('--sellier-max-it', type=int, required=False,
                    default=1000, help="Maximum Sellier iterations admitted. Default: 1000")
parser.add_argument('--sellier-inner-steps', type=int, required=False,
                    default=1, help="Maximum Sellier direct solver ramp steps. Default: 1")
parser.add_argument('--sellier-aa-m', type=int, required=False,
                    default=0, help="Acceleration depth for Sellier. Default: 0")
parser.set_defaults(armijo=False, adaptive=False)
args = parser.parse_args()

geom_name = args.heart_geometry_name
output_filename = args.output_file
if geom_name:
    assert geom_name in ("4mm", "h4_v2_ASCII", "4mm_ref",
                         "h4_v2_ASCII_ref"), "Wrong heart geometry name {}. See --help for further details".format(geom_name)
solver_type = args.solver_type  # direct | sellier
assert solver_type in (
    "inverse", "sellier"), "Wrong solver type {}. See --help for further details".format(solver_type)
if solver_type == "inverse":
    solver_type = InverseFormulationType.direct
elif solver_type == "sellier":
    solver_type = InverseFormulationType.sellier
n_times_back = args.backward_steps  # Number of pseudo-timesteps
n_times_forward = args.forward_steps  # Number of pseudo-timesteps
pressure_known = args.endo_pressure  # Known endocardial pressure in Pa
active_stress_known = args.active_peak  # Known endocardial pressure in Pa
volumetric_z = args.z_load_scale  # Gravity
neumann_z = args.z_surf_scale  # Scale of verticale surface traction
u_deg = args.fe_degree  # Displacement degree
adaptive = args.adaptive

sellier_params = {"alpha": args.sellier_relaxation,
                  "atol_inc": args.sellier_atol,
                  "rtol_inc": args.sellier_rtol,
                  "max_it": args.sellier_max_it,
                  "inner_steps": args.sellier_inner_steps,
                  "do_armijo": args.armijo,
                  "aa_m": args.sellier_aa_m}  # aa_m = 1 works best

PETSc.Options().setValue("-snes_type", "newtonls")
PETSc.Options().setValue("-ksp_type", "gmres")
#PETSc.Options().setValue("-pc_type", "hypre")
PETSc.Options().setValue("-pc_type", "lu")
PETSc.Options().setValue("-pc_factor_mat_solver_type", "mumps")
PETSc.Options().setValue("-snes_atol", 1e-14)
PETSc.Options().setValue("-snes_rtol", 1e-6)
PETSc.Options().setValue("-snes_stol", 0.0)
PETSc.Options().setValue("-ksp_norm_type", "unpreconditioned")
PETSc.Options().setValue("-ksp_max_it", 1000)
PETSc.Options().setValue("-ksp_atol", 0.0)
PETSc.Options().setValue("-ksp_rtol", 1e-6)
PETSc.Options().setValue("-ksp_gmres_restart", 100)
PETSc.Options().setValue("-snes_max_it", 1000)
PETSc.Options().setValue("-snes_linesearch_type", "basic")
PETSc.Options().setValue("-snes_converged_reason", None)
#PETSc.Options().setValue("-ksp_converged_reason", None)
PETSc.Options().setValue("-snes_error_if_not_converged", "true")

f = s = n = None
mesh = markers = None
ENDO = EPI = BASE = 0
meshname = "prolate_{}".format(geom_name)

mesh, markers, ENDO, EPI, BASE = prolateGeometry(meshname)
f, s, n = generateFibers(mesh, markers, ENDO, EPI, BASE, output_dir="output")

fibers = (f, s, n)
variational = False
if output_filename:
    outmesh = getXDMFFile(output_filename, mesh=True).write(mesh)
time_0 = time()
dofs, sell_its, nl_its, l_its = PS.prestress_solver_cardiac(mesh, markers, ENDO, EPI, BASE, fibers, u_deg, n_times_back, n_times_forward, neumann_z, volumetric_z,
                                                            pressure_known, active_stress_known, variational, getXDMFFile(output_filename), solver_type=solver_type, sellier_params=sellier_params, adaptive=adaptive)
time_f = time()
parprint("Dofs:", dofs)
parprint("Avg sell its:", sell_its)
parprint("Avg nonlinear its:", nl_its)
parprint("Tot linear its:", l_its)
parprint("Avg linear its:", l_its / max(nl_its, 1))
parprint("Solution time:", time_f - time_0)
