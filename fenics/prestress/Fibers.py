# -*- coding: utf-8 -*-
# Code by: Francesco Regazzoni

import dolfin as df


def getH1projection(f, V):

    mesh = V.mesh()
    uv = df.TrialFunction(V)
    vv = df.TestFunction(V)
    A = df.dot(uv, vv) + df.inner(df.grad(uv), df.grad(vv))
    L = df.dot(f, vv) + df.inner(df.grad(f), df.grad(vv))
    sol = df.Function(V)
    df.solve(A * df.dx(mesh) == L * df.dx(mesh), sol)
    return sol


def getTransmuralCoordinate(mesh, boundary_markers, ENDO, EPI, degree=2):

    df.set_log_level(30)
    dx = df.dx(mesh)

    V = df.FunctionSpace(mesh, 'CG', degree)
    bc1 = df.DirichletBC(V, df.Constant(0.0), boundary_markers, ENDO)
    bc2 = df.DirichletBC(V, df.Constant(1.0), boundary_markers, EPI)

    phi = df.Function(V, name="phi_transmural")
    phi_trial = df.TrialFunction(V)
    psi = df.TestFunction(V)

    df.solve(df.dot(df.grad(phi_trial), df.grad(psi)) * dx ==
             df.Constant(0.0) * psi * dx, phi, [bc1, bc2])

    return phi


def getApicobasalCoordinate(mesh, boundary_markers, BASE, degree=2):

    df.set_log_level(30)
    dx = df.dx(mesh)

    V = df.FunctionSpace(mesh, 'CG', degree)

    phi_trial = df.TrialFunction(V)
    psi = df.TestFunction(V)

    a = df.dot(df.grad(phi_trial), df.grad(psi)) * df.dx
    L = df.Constant(0.0) * psi * df.dx

    quotes = mesh.coordinates()[:, 2]
    min_quote = df.MPI.max(df.MPI.comm_world, quotes.max())

    def apex(x):
        result = abs(x[2] - min_quote) < df.DOLFIN_EPS
        return result
    bcs = [df.DirichletBC(V, df.Constant(1.0), boundary_markers, BASE),
           df.DirichletBC(V, df.Constant(0.0), apex, method='pointwise')]

    phi = df.Function(V, name="phi_apicobasal")
    df.solve(a == L, phi, bcs)
    return phi


def generateAnalyticFibers():
    f, s, n = df.Constant((1, 0, 0)), df.Constant((0, 1/df.sqrt(2), 1/df.sqrt(2))
                                                  ), df.Constant((0, 1/df.sqrt(2), -1/df.sqrt(2)))
    return f, s, n


def generateFibers(mesh, boundary_markers, ENDO, EPI, BASE, output_dir=None):

    if df.MPI.rank(df.MPI.comm_world) == 0:
        print("generating fibers...", flush=True)
    theta_endo = 60.
    theta_epi = -60.

    degree = 1
    phi_transmural = getTransmuralCoordinate(
        mesh, boundary_markers, ENDO, EPI, degree=degree)

    def Max(a, b): return (a + b + abs(a - b)) / df.Constant(2.)
    def Min(a, b): return (a + b - abs(a - b)) / df.Constant(2.)

    W = df. VectorFunctionSpace(mesh, 'CG', degree)
    n = df.Function(W)  # W?

    # -1: analytical fibers
    # 0: SR (on fields) NB: f and s are not perfectly orthogonal
    # 1: SR (on dofs)
    # 2: BT (on dofs)
    alg_type = 0  # 0 faster, but f and s are not perfectly orthogonal
    if alg_type == -1:

        f = df.project(df.Expression(('1.0', '0.0', '0.0'), degree=degree), W)
        s = df.project(df.Expression(('0.0', '1.0', '0.0'), degree=degree), W)
        n = df.project(df.cross(f, s), W)

    elif alg_type == 0:
        s = df.grad(phi_transmural)
        s = s / df.sqrt(df.inner(s, s))
        k = df.Constant((.0, .0, -1.))
        kp_tilde = k - df.dot(k, s) * s
        kp = kp_tilde / df.sqrt(df.inner(kp_tilde, kp_tilde))
        f_tilde = df.cross(s, kp)
        f_tilde = f_tilde / df.sqrt(df.inner(f_tilde, f_tilde))
        theta = (theta_endo + (theta_epi - theta_endo)
                 * phi_transmural) * df.pi / 180.0
        f = f_tilde + df.sin(theta) * df.cross(s, f_tilde) + 2.2 * \
            (df.sin(theta * .5))**2 * df.cross(s, df.cross(s, f_tilde))
        f = - f / df.sqrt(df.inner(f, f))
        s = df.project(s, W)
        f = df.project(f, W)

    elif alg_type == 1 or alg_type == 2:

        import numpy as np
        ndof_local = phi_transmural.vector().get_local().size
        s_vec = np.empty((ndof_local, 3))
        s_vec_tot = df.project(df.grad(phi_transmural), W).vector().get_local()
        s_vec[:, 0] = s_vec_tot[0::3]
        s_vec[:, 1] = s_vec_tot[1::3]
        s_vec[:, 2] = s_vec_tot[2::3]
        s = s_vec

        if alg_type == 2:
            phi_apicobasal = getApicobasalCoordinate(
                mesh, boundary_markers, BASE, degree=1)
            k_vec = np.empty((ndof_local, 3))
            k_vec_tot = df.project(
                df.grad(phi_apicobasal), W).vector().get_local()
            k_vec[:, 0] = k_vec_tot[0::3]
            k_vec[:, 1] = k_vec_tot[1::3]
            k_vec[:, 2] = k_vec_tot[2::3]
            k = k_vec
        else:
            pass

        phi_vec = phi_transmural.vector().get_local()

        s_norm = np.empty(ndof_local)
        theta = np.empty(ndof_local)
        kp_tilde = np.empty((ndof_local, 3))
        kp = np.empty((ndof_local, 3))
        f_tilde = np.empty((ndof_local, 3))
        f = np.empty((ndof_local, 3))
        n = np.empty((ndof_local, 3))
        for i in range(ndof_local):
            s_norm[i] = np.sqrt(np.inner(s[i, :], s[i, :]))
            s[i, :] = s[i, :] / s_norm[i]
            kp_tilde[i, :] = k[i, :] - np.inner(k[i, :], s[i, :]) * s[i, :]
            kp[i, :] = kp_tilde[i, :] / \
                np.sqrt(np.inner(kp_tilde[i, :], kp_tilde[i, :]))
            f_tilde[i, :] = np.cross(s[i, :], kp[i, :])
            f_tilde[i, :] = f_tilde[i, :] / \
                np.sqrt(np.inner(f_tilde[i, :], f_tilde[i, :]))
            theta[i] = (theta_endo + (theta_epi - theta_endo)
                        * phi_vec[i]) * np.pi / 180.0
            f[i, :] = f_tilde[i, :] + np.sin(theta[i]) * np.cross(s[i, :], f_tilde[i, :]) + \
                2.2 * (np.sin(theta[i] * .5))**2 * \
                np.cross(s[i, :], np.cross(s[i, :], f_tilde[i, :]))
            f[i, :] = - f[i, :] / np.sqrt(np.inner(f[i, :], f[i, :]))
            n[i, :] = np.cross(f[i, :], s[i, :])

        f_vec = np.empty(ndof_local * 3)
        for i in range(3):
            f_vec[i::3] = f[:, i]

        s_vec = np.empty(ndof_local * 3)
        for i in range(3):
            s_vec[i::3] = s[:, i]

        n_vec = np.empty(ndof_local * 3)
        for i in range(3):
            n_vec[i::3] = n[:, i]

        f = df.Function(W)
        s = df.Function(W)
        n = df.Function(W)
        f.vector().set_local(f_vec)
        f.vector().apply("insert")
        s.vector().set_local(s_vec)
        s.vector().apply("insert")
        n.vector().set_local(n_vec)
        n.vector().apply("insert")

    if df.MPI.rank(df.MPI.comm_world) == 0:
        print("fibers generated!", flush=True)

    f.rename("f", "f")
    s.rename("s", "s")
    n.rename("n", "n")
    if output_dir:
        xdmf = df.XDMFFile("{}/fibers.xdmf".format(output_dir))
        xdmf.parameters["functions_share_mesh"] = True
        xdmf.parameters["flush_output"] = True
        xdmf.write(phi_transmural, 0)
        if alg_type == 2:
            xdmf.write(phi_apicobasal, 0)
        xdmf.write(f, 0)
        xdmf.write(s, 0)
        xdmf.write(n, 0)
        xdmf.close()

    if df.MPI.rank(df.MPI.comm_world) == 0:
        print("done generating fibers!", flush=True)
    return f, s, n
