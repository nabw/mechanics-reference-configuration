from firedrake import * 
import numpy as np
parameters["form_compiler"]["quadrature_degree"] = 4
def parprint(*args):
    if COMM_WORLD.rank == 0:
        print("[=]", *args, flush=True)


# Params 
saveEvery = 1
incomp=False

def solveID(meshname, output, pressure_known):
    # Geometry and measures
    mesh = Mesh(meshname)
    dx_mesh = dx(mesh)
    ds_mesh = ds(mesh)
    ds_endo = ds_mesh(9)
    ds_epi = ds_mesh(10)
    
    # Functional setting
    V = None
    if incomp:
        Vu = VectorFunctionSpace(mesh, 'CG', 2)
        Vp = FunctionSpace(mesh, 'CG', 1)
        V = Vu * Vp
        sol = Function(V)
        u, p = split(sol)
        v, q = TestFunctions(V)
    else:
        V = VectorFunctionSpace(mesh, 'CG', 2)
        Vu = V
        u = Function(V)
        sol = u
        v = TestFunction(V)
    
    # Ramping (HoMoToPy CoNtInUaTiOn)
    ramp = Constant(0.0, domain=mesh)
    
    
    # Variational formulation
    f = Identity(3) + grad(u)  # Inverse tensor for inverse problem
    j = det(f)
    F = variable(inv(f))  # Compute original one to diff
    J = det(F)
    Cbar = J**(-2/3) * F.T * F
    
    
    # Usyk,. mc Culloch 2002
    f0 = Constant((1,0,0), domain=mesh) # Isotropic
    s0 = Constant((0,1,0), domain=mesh)
    n0 = Constant((0,0,1), domain=mesh)
    Cg = .88e3   # [Pa]
    bf = 8       # [-]
    bs = 6       # [-]
    bn = 3       # [-]
    bfs = 12      # [-]
    bfn = 3       # [-]
    bsn = 3       # [-]
    k = Constant(5e6, domain=mesh)
    E = 0.5*(Cbar - Identity(3))
    
    b_list = [bf,bs,bn,bfs,bfn,bsn]
    bavg = sum(b_list) / len(b_list)
    Q = Constant(bavg, domain=mesh) * inner(E, E) # bs seems like a good average
    WP = 0.5*Constant(Cg, domain=mesh)*(exp(Q)-1)
    WV = Constant(k, domain=mesh)/2*(J-1)*ln(J)
    
    # Finally build the Robin condition terms, Pfaller et al.
    k_perp = Constant(2e5, domain=mesh)  # [Pa/m]
    c_perp = Constant(5e3, domain=mesh)  # [Pa*s/m]
    
    psi = WP + WV
    if incomp: 
        psi += p * (J-1)
    
    nn = FacetNormal(mesh)
    k_perp = Constant(1e2, domain=mesh) # artificial
    ts_robin = outer(nn, nn)*k_perp*u + (Identity(3) - outer(nn, nn)) * k_perp/10*u  # flip signs for inverse displacement
    P = diff(psi, F) 
    
    F_form = inner(j * P, grad(v) * inv(f)) * dx_mesh - ramp * Constant(-pressure_known, domain=mesh) * dot(nn, v) * ds_endo - dot(ts_robin, v) * ds_epi
    if incomp: 
        F_form += (j-1) * q * dx_mesh
    
    snes_params = {"type": "newtonls",
                   "linesearch_type": "none",
                   #"monitor": None,
                   "atol": 1e-10,
                   "rtol": 1e-4,
                   "snes_error_if_not_converged": "false",
                   "stol": 0.0}
    solver_parameters = {"snes": snes_params,
                         "ksp_type": "preonly",
                         "pc_type": "lu",
                         "pc_factor_mat_solver_type": "mumps",
                         }
    
    
    
    outfile = File(output)
    u_out = Function(Vu, name='u')
    outfile.write(u_out, time=0)
    i=0
    problem = NonlinearVariationalProblem(F_form, sol)
    solver = NonlinearVariationalSolver(problem, solver_parameters=solver_parameters)
    
    # Adaptive stepping
    r = 0
    step = 1
    sol_prev = Function(V)
    while r<1:
        r_test = r+step
        parprint("Solving ramp={:.3e}".format(r_test))
        ramp.assign(r_test)
        sol.assign(sol_prev) # snes saves diverged solutions!
        solver.solve()
        converged = solver.snes.getConvergedReason() > 0
        if not converged:
            step = step/2
            parprint("Diverged, changed step to {:.3e}".format(step))
            if step < 1e-10:
                break
                parprint("Step too small.... bye bye")
        else:
            if incomp:
                u_out.assign(sol.subfunctions[0])
            else:
                u_out.assign(u)
            sol_prev.assign(sol)
            r = r_test
            if i % saveEvery == 0: 
                outfile.write(u_out,time=r)
            i+=1
    
    outfile.write(u_out,time=r)
    parprint("Done")


solveID("ring-mod-semi.msh", "output/semicricle.pvd", 5e2)
solveID("ring-mod.msh", "output/halfmoon.pvd", 5e5)
