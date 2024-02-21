import numpy as np
from numpy import linalg
import copy
from dolfin import *
from .Tools import *
from .Enums import *
from .MechanicsForms import backwardForm, forwardForm
from .MechanicsForms import backwardForm_cardiac, forwardForm_cardiac
from .SNESProblem import SNESProblem
from petsc4py import PETSc
from mpi4py import MPI
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["cpp_optimize_flags"] = '-O3 -Ofast -march=native -mtune=native'
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 4
parameters["mesh_partitioner"] = "ParMETIS"
set_log_active(False)


adaptive_decrease = 0.5


def avg(lst):
    if len(lst) == 0:
        return 0
    else:
        return sum(lst)/len(lst)


def getTimesList(n_times):
    do_times = False
    times = []
    if n_times == 0:
        pass
    elif n_times == 1:
        times.append(1.0)
        do_times = True
    else:
        times = np.linspace(1.0 / n_times, 1.0, n_times)
        do_times = True
    return do_times, times


def saveSolution(mesh, coords, xdmffile, u, uinv, uforw, do_backward, do_forward):
    do_something = do_backward or do_forward
    if xdmffile and do_something:
        t = 0.0
        mesh.coordinates()[:] = coords
        if do_backward:
            u.assign(uinv)
            xdmffile.write(u, t)
            t += 1.0
        if do_forward:
            deformMesh(uinv)
            u.assign(uforw)
            xdmffile.write(u, t)


def prestress_solver_slab(mesh, markers, dirichlet_markers, neumann_markers, u_deg, n_times_back, n_times_forward, neumann_load, volumetric_load, xdmffile, solver_type=InverseFormulationType.direct, sellier_params=None, adaptive=False):
    parprint("")
    parprint("=================================")
    dim = mesh.geometric_dimension()
    coords = copy.deepcopy(mesh.coordinates())  # output only

    dx_mesh = dx(mesh)
    ds = Measure('ds', domain=mesh, subdomain_data=markers)
    ds_N = ds(neumann_markers)
    nn = FacetNormal(mesh)

    # Function spaces
    u_el = VectorElement("CG", mesh.ufl_cell(), u_deg)
    V = FunctionSpace(mesh, u_el)
    bcs = [DirichletBC(V, Constant((0, 0, 0)), markers, b)
           for b in dirichlet_markers]

    # Functions
    u = Function(V, name="d")
    v = TestFunction(V)
    uinv = Function(V, name="d_inverse")  # output only
    uforw = Function(V, name="d_forward")  # output only

    # Scale functions
    ramp = Expression("t", t=0, degree=0)

    # Initialize ramps
    do_backward, times_backward = getTimesList(n_times_back)
    do_forward, times_forward = getTimesList(n_times_forward)

    # Solve backward problem according to chosen strategy
    avg_nl = avg_l = avg_sell = 0
    if do_backward:
        if solver_type == InverseFormulationType.direct:
            FF = backwardForm(dx_mesh, ds_N, nn, u, v, ramp,
                              neumann_load, volumetric_load)
            avg_nl, avg_l = inverseSolver(
                FF, times_backward, xdmffile, u, ramp, mesh, adaptive, bcs)
        elif solver_type == InverseFormulationType.sellier:
            # One of the advantages of the Sellier algorithm is that it leverages the forward solver
            FF = forwardForm(dx_mesh, ds_N, nn, u, v, ramp,
                             neumann_load, volumetric_load)
            FF0 = backwardForm(dx_mesh, ds_N, nn, u, v, ramp,
                               neumann_load, volumetric_load)
            avg_sell, avg_nl, avg_l = sellierSolver(
                FF, FF0, times_backward, xdmffile, u, ramp, mesh, sellier_params, adaptive, bcs)

        uinv.assign(u)
    if do_forward:
        u.vector().zero()  # Inverse is a bad initial point
        FF = forwardForm(dx_mesh, ds_N, nn, u, v, ramp,
                         neumann_load, volumetric_load)

        for t in times_forward:
            parprint("\t\t Forward ramp, t={:.4f}".format(t))
            ramp.t = t
            JJ = derivative(FF, u, TrialFunction(V))
            problem = SNESProblem(FF, JJ, u, bcs)
            problem.solve()
        uforw.assign(u)
    saveSolution(mesh, coords, xdmffile, u, uinv,
                 uforw, do_backward, do_forward)
    return V.dim(), avg_sell, avg_nl, avg_l


def prestress_solver_cardiac(mesh, markers, endo_marker, epi_marker, base_marker, fibers, u_deg, n_times_back, n_times_forward, neumann_z, volumetric_z, pressure_known, active_stress_known, variational, xdmffile, solver_type=InverseFormulationType.direct, sellier_params=None, adaptive=False):
    global EXPORT_INDEX  # If not, python thinks the variable is local and gives an error
    parprint("")
    parprint("=================================")
    dim = mesh.geometric_dimension()
    f0, s0, n0 = fibers
    Vf = VectorFunctionSpace(mesh, 'CG', 1)
    coords = copy.deepcopy(mesh.coordinates())  # output only

    dx_mesh = dx(mesh)
    ds = Measure('ds', domain=mesh, subdomain_data=markers)
    ds_base = ds(base_marker)
    ds_endo = ds(endo_marker)
    ds_epi = ds(epi_marker)
    nn = FacetNormal(mesh)

    # Function spaces
    u_el = VectorElement("CG", mesh.ufl_cell(), u_deg)
    V = FunctionSpace(mesh, u_el)

    # Functions
    u = Function(V, name="d")
    v = TestFunction(V)
    uinv = Function(V, name="d_inverse")  # output only
    uforw = Function(V, name="d_forward")  # output only

    # Scale functions
    ramp = Expression("t", t=0, degree=0)

    # Initialize ramps
    do_backward, times_backward = getTimesList(n_times_back)
    do_forward, times_forward = getTimesList(n_times_forward)

    # Solve backward problem according to chosen strategy
    avg_nl = avg_l = avg_sell = 0
    if do_backward:
        if solver_type == InverseFormulationType.direct:
            FF = backwardForm_cardiac(f0, s0, n0, dx_mesh, ds_base, ds_endo,
                                      ds_epi, nn, u, v, ramp, pressure_known, active_stress_known, neumann_z, volumetric_z, variational=variational)
            avg_nl, avg_l = inverseSolver(
                FF, times_backward, xdmffile, u, ramp, mesh, adaptive)
        elif solver_type == InverseFormulationType.sellier:
            # One of the advantages of the Sellier algorithm is that it leverages the forward solver
            FF = forwardForm_cardiac(f0, s0, n0, dx_mesh, ds_base, ds_endo,
                                     ds_epi, nn, u, v, ramp, pressure_known, active_stress_known, neumann_z, volumetric_z)
            FF0 = backwardForm_cardiac(f0, s0, n0, dx_mesh, ds_base, ds_endo,
                                       ds_epi, nn, u, v, ramp, pressure_known, active_stress_known, neumann_z, volumetric_z, variational=variational)
            avg_sell, avg_nl, avg_l = sellierSolver(
                FF, FF0, times_backward, xdmffile, u, ramp, mesh, sellier_params, adaptive)

        uinv.assign(u)

    if do_forward:
        u.vector().zero()  # Inverse is a bad initial point
        FF = forwardForm_cardiac(f0, s0, n0, dx_mesh, ds_base, ds_endo, ds_epi, nn,
                                 u, v, ramp, pressure_known, active_stress_known, neumann_z, volumetric_z)
        for t in times_forward:
            parprint("\t\t Forward ramp, t={:.4f}".format(t))
            ramp.t = t
            JJ = derivative(FF, u, TrialFunction(V))
            problem = SNESProblem(FF, JJ, u)
            problem.solve()
        uforw.assign(u)

    saveSolution(mesh, coords, xdmffile, u, uinv,
                 uforw, do_backward, do_forward)
    return V.dim(), avg_sell, avg_nl, avg_l


def deformMesh(u):
    mesh = u.function_space().mesh()
    u_deg = u.function_space().ufl_element().degree()
    if u_deg > 1:
        Vlow = VectorFunctionSpace(mesh, 'CG', 1)
        ulow = interpolate(u, Vlow)
        ALE.move(mesh, ulow)
    else:
        ALE.move(mesh, u)


def inverseSolver(FF, times_backward, xdmffile, u, ramp, mesh, adaptive, bcs=None):
    t = 0
    dt = 1 if len(
        times_backward) == 1 else times_backward[1] - times_backward[0]
    u_prev = u.copy(True) if adaptive else None
    nl_its, l_its = [], []
    while True:
        ramp.t = t+dt
        if ramp.t > 1.0:
            break
        parprint("\t\t Backward ramp, t={:.4f}".format(t+dt))
        V = u.function_space()
        du = TrialFunction(V)
        v = TestFunction(V)
        JJ = derivative(FF, u, du)
        problem = SNESProblem(FF, JJ, u, bcs)
        if not adaptive:
            nl, l = problem.solve()
            nl_its.append(nl)
            l_its.append(l)
            t = min(t+dt, 1.0)
        else:  # if adaptive
            try:
                nl, l = problem.solve()
                nl_its.append(nl)
                l_its.append(l)
                u_prev.assign(u)
                t = min(t+dt, 1.0)
            except:
                dt = adaptive_decrease * dt
                parprint("Step decreased to", dt)
                u.assign(u_prev)

    deformMesh(u)
    return avg(nl_its), avg(l_its)


def sellierSolver(FF, FF0, times_backward, xdmffile, u, ramp, mesh, sellier_params, adaptive, bcs=None):
    from .AndersonAcceleration import AndersonAcceleration
    global EXPORT_INDEX

    V = u.function_space()
    position_target = u.copy(True)
    position_current = u.copy(True)
    expr_position = Expression(("x[0]", "x[1]", "x[2]"), degree=1)
    position_target.interpolate(expr_position)
    position_current.interpolate(expr_position)
    x_curr = Function(V)
    incr = Function(V)
    incr_prev = u.copy(True)  # Used by Armijo line search
    temp_incr = u.copy(True)  # Used by Armijo line search
    x_curr_armijo = Function(V)
    do_armijo = sellier_params["do_armijo"]

    def solveBackwardProblem(u, tf, verbose=False):
        # This function actually solves the forward problem, but it is
        # used as the forward solver for Sellier, which is the backward
        # solver. It can be done up to 'tf' to allow for smaller Sellier
        # steps as done in lifex.
        t = 0
        n_times_inner = sellier_params["inner_steps"]
        dt = tf if n_times_inner == 1 else tf / n_times_inner
        nl_its, l_its = [], []
        # for t in times_backward:
        while True:
            if verbose:
                parprint("\t\t Backward ramp, t={:.4f}".format(t+dt))
            ramp.t = t+dt
            JJ = derivative(FF, u, TrialFunction(u.function_space()))
            problem = SNESProblem(FF, JJ, u, bcs)
            if not adaptive:
                nl, l = problem.solve()
                nl_its.append(nl)
                l_its.append(l)
                t = min(t+dt, tf)
            else:  # if adaptive
                try:
                    problem.solve()
                    nl_its.append(nl)
                    l_its.append(l)
                    t = min(t+dt, tf)
                except:
                    dt = adaptive_decrease * dt
            if t >= tf:
                break
        return avg(nl_its), avg(l_its)

    anderson = AndersonAcceleration(sellier_params["aa_m"])
    nl_its, l_its, sell_its = [], [], []
    for tf in times_backward:
        parprint("Solving Sellier ramp, t={:.2f}".format(tf))
        alpha = sellier_params["alpha"]
        incr.assign(incr_prev)
        x_curr.interpolate(expr_position)
        error0 = x_curr.vector().norm('l2')  # Normalize w.r.t displacement
        if error0 < 1e-10:
            error0 = 1

        it = 0
        error = 1
        ps_atol_inc = sellier_params["atol_inc"]
        ps_rtol_inc = sellier_params["rtol_inc"]
        ps_max_it = sellier_params["max_it"]
        anderson.reset()

        alpha = sellier_params["alpha"]  # Reset alpha for Armijo
        while error > ps_atol_inc and error / error0 > ps_rtol_inc and it < ps_max_it:
            if do_armijo and it > 0:
                x_curr.interpolate(expr_position)  # won't move
                Ls = [1.0, 1.0/2.0, 1.0/4.0, 1.0/8.0]
                results = []
                nl_its.append(0)
                l_its.append(0)
                nl_it, l_it = solveBackwardProblem(u, tf, verbose=False)
                anderson.get_next_vector(u.vector().vec())
                # notation from Sellier, deformed points
                x_curr.interpolate(expr_position)
                uN = u.vector() + x_curr.vector()  # Deformed points
                for l in Ls:
                    # Update weighted displacement
                    incr.vector().zero()
                    incr.vector().axpy(-1, position_target.vector())
                    incr.vector().axpy(1, uN)  # This is R_l as in Marx et al.
                    error_l = incr.vector().norm('linf')
                    temp_incr.assign(incr)
                    temp_incr.vector().axpy(-1, incr_prev.vector()
                                            )  # This is R_l - R_{k-1}
                    _num = incr_prev.vector().inner(temp_incr.vector())
                    _denom = temp_incr.vector().inner(temp_incr.vector())
                    alpha_inner = -alpha * l * _num / _denom
                    # This is not documented, but without it, the relaxation coeff
                    # goes to zero in some scenarios
                    alpha_inner = max(alpha_inner, 0.7)
                    incr.vector()[:] *= -alpha_inner

                    results.append((l, alpha_inner, error_l, incr.copy(True)))

                    nl_its[-1] += nl_it
                    l_its[-1] += l_it
                    if error_l < error:
                        break
                error_min = 1e10
                best = None
                for l, alpha_inner, err, inc in results:
                    if err < error_min:
                        best = inc
                        error = err
                        alpha = alpha_inner
                incr.assign(best)
            else:
                nl_it, l_it = solveBackwardProblem(u, tf, verbose=False)
                nl_its.append(nl_it)
                l_its.append(l_it)
                anderson.get_next_vector(u.vector().vec())
                # compute increment
                x_curr.interpolate(expr_position)
                # notation from Sellier, deformed points
                uN = u.vector() + x_curr.vector()
                incr.vector().zero()
                incr.vector().axpy(alpha, position_target.vector())
                incr.vector().axpy(-alpha, uN)
                error = incr.vector().norm('linf') / sellier_params["alpha"]
            deformMesh(incr)  # Coords are now X_k
            incr_prev.assign(incr)
            parprint("\tSellier inner iteration {}, err abs={:.3e}, err rel={:.3e}".format(
                it, error, error/error0))
            it += 1

        sell_its.append(float(it-1))
    # -1 as it goes up by one at the end
    return avg(sell_its), avg(nl_its), avg(l_its)
