from dolfin import *


def activeStress(F, f0, active_stress_known):
    Ff0 = F * f0
    return Constant(active_stress_known) / sqrt(Ff0 ** 2) * outer(Ff0, f0)


def backwardForm_cardiac(f0, s0, n0, dx_mesh, ds_base, ds_endo, ds_epi, nn, u, v, ramp, pressure_known, active_stress_known, neumann_z, volumetric_z, variational=False):

    dim = u.ufl_shape[0]
    f = Identity(dim) + grad(u)  # Inverse tensor for inverse problem
    j = det(f)
    F = None
    if variational:
        F = inv(f)  # Compute original one to diff
    else:
        F = variable(inv(f))  # Compute original one to diff
    J = det(F)
    Cbar = J**(-2/3) * F.T * F
    vec_z = Constant((0, 0, 1.))

    # Usyk,. mc Culloch 2002
    Cg = .88e3   # [Pa]
    bf = 8       # [-]
    bs = 6       # [-]
    bn = 3       # [-]
    bfs = 12      # [-]
    bfn = 3       # [-]
    bsn = 3       # [-]
    k = Constant(5e4)
    E = 0.5*(Cbar - Identity(dim))
    Eff, Efs, Efn = inner(E*f0, f0), inner(E*f0, s0), inner(E*f0, n0)
    Esf, Ess, Esn = inner(E*s0, f0), inner(E*s0, s0), inner(E*s0, n0)
    Enf, Ens, Enn = inner(E*n0, f0), inner(E*n0, s0), inner(E*n0, n0)

    Q = Constant(bf) * Eff**2 \
        + Constant(bs) * Ess**2 \
        + Constant(bn) * Enn**2 \
        + Constant(bfs) * 2.0 * Efs**2 \
        + Constant(bfn) * 2.0 * Efn**2 \
        + Constant(bsn) * 2.0 * Esn**2
    WP = 0.5*Constant(Cg)*(exp(Q)-1)
    WV = Constant(k)/2*(J-1)*ln(J)

    # Finally build the Robin condition terms, Pfaller et al.
    k_perp = Constant(2e5)  # [Pa/m]
    c_perp = Constant(5e3)  # [Pa*s/m]
    rhos = Constant(1e3)
    gravity = Constant(volumetric_z) * vec_z

    psi = WP + WV

    ts_robin = outer(nn, nn)*k_perp*u + (Identity(dim) - outer(nn, nn)) * \
        k_perp/10*u  # flip signs for inverse displacement
    F_form = None
    if variational:
        F_form = derivative(j * psi * dx_mesh, u, v) - ramp * Constant(-pressure_known) * dot(nn, v) * ds_endo - \
            dot(ts_robin, v) * ds_epi - ramp * rhos * dot(gravity, v) * dx
    else:
        P = diff(psi, F) + ramp * activeStress(F, f0, active_stress_known)
        F_form = inner(j * P, grad(v) * inv(f)) * dx_mesh - ramp * Constant(-pressure_known) * dot(nn, v) * ds_endo - \
            ramp * Constant(neumann_z) * dot(vec_z, v) * ds_base - \
            dot(ts_robin, v) * ds_epi - ramp * rhos * dot(gravity, v) * dx

    return F_form


def forwardForm_cardiac(f0, s0, n0, dx_mesh, ds_base, ds_endo, ds_epi, nn, u, v, ramp, pressure_known, active_stress_known, neumann_z, volumetric_z):
    dim = u.ufl_shape[0]
    F = Identity(dim) + grad(u)
    F = variable(F)  # Inverse tensor for inverse problem
    Cbar = det(F)**(-2/3) * F.T * F
    J = det(F)
    vec_z = Constant((0, 0, 1.))

    # Usyk,. mc Culloch 2002
    Cg = .88e3   # [Pa]
    bf = 8       # [-]
    bs = 6       # [-]
    bn = 3       # [-]
    bfs = 12      # [-]
    bfn = 3       # [-]
    bsn = 3       # [-]
    k = Constant(5e4)
    E = 0.5*(Cbar - Identity(dim))
    Eff, Efs, Efn = inner(E*f0, f0), inner(E*f0, s0), inner(E*f0, n0)
    Esf, Ess, Esn = inner(E*s0, f0), inner(E*s0, s0), inner(E*s0, n0)
    Enf, Ens, Enn = inner(E*n0, f0), inner(E*n0, s0), inner(E*n0, n0)

    Q = Constant(bf) * Eff**2 \
        + Constant(bs) * Ess**2 \
        + Constant(bn) * Enn**2 \
        + Constant(bfs) * 2.0 * Efs**2 \
        + Constant(bfn) * 2.0 * Efn**2 \
        + Constant(bsn) * 2.0 * Esn**2
    WP = 0.5*Constant(Cg)*(exp(Q)-1)
    WV = Constant(k)/2*(J-1)*ln(J)

    # Finally build the Robin condition terms, Pfaller et al.
    k_perp = Constant(2e5)  # [Pa/m]
    c_perp = Constant(5e3)  # [Pa*s/m]
    rhos = Constant(1e3)
    gravity = Constant(volumetric_z) * vec_z

    cof = J * inv(F).T
    cofnorm = sqrt(dot(cof * nn, cof * nn))
    # Valid at spatial config, so we pull back the normal.
    NN = 1 / cofnorm * cof * nn
    ts_robin = -outer(NN, NN)*k_perp*u - \
        (Identity(dim) - outer(NN, NN))*k_perp/10*u

    psi = WP + WV
    P = diff(psi, F)
    F_form = inner(P + ramp * activeStress(F, f0, active_stress_known), grad(v)) * dx_mesh - ramp * Constant(-pressure_known) * dot(cof * nn, v) * ds_endo - \
        ramp * Constant(neumann_z) * dot(vec_z, v) * cofnorm * ds_base - \
        dot(ts_robin, v) * cofnorm * ds_epi - \
        ramp * rhos * dot(gravity, v) * 1 / J * dx

    return F_form


def backwardForm(dx_mesh, ds_N, nn, u, v, ramp, neumann_load, volumetric_load):

    dim = u.ufl_shape[0]
    f = Identity(dim) + grad(u)  # Inverse tensor for inverse problem
    j = det(f)
    F = None
    F = variable(inv(f))  # Compute original one to diff
    J = det(F)
    Cbar = J**(-2/3) * F.T * F
    vec_z = Constant((0, 0, 1.))

    # Elasticity parameters, Neo-Hookean from FEniCS documentation
    E = 1.0e4
    nu = 0.3
    mu = Constant(E/(2*(1 + nu)))
    lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu)))
    Ic = tr(Cbar)
    psi = (mu / 2) * (Ic - 3) + 0.5 * lmbda * (J-1) * ln(J)
    P = diff(psi, F)
    rhos = Constant(1e3)
    gravity = Constant(volumetric_load) * vec_z

    F_form = inner(j * P, grad(v) * inv(f)) * dx_mesh - ramp * Constant(
        neumann_load) * dot(vec_z, v) * ds_N - ramp * rhos * dot(gravity, v) * dx
    return F_form


def forwardForm(dx_mesh, ds_N, nn, u, v, ramp, neumann_load, volumetric_load):

    dim = u.ufl_shape[0]
    F = Identity(dim) + grad(u)  # Inverse tensor for inverse problem
    F = variable(F)
    J = det(F)
    Cbar = J**(-2/3) * F.T * F
    vec_z = Constant((0, 0, 1.))

    # Elasticity parameters, Neo-Hookean from FEniCS documentation
    E = 1.0e4
    nu = 0.3
    mu = Constant(E/(2*(1 + nu)))
    lmbda = Constant(E*nu/((1 + nu)*(1 - 2*nu)))
    Ic = tr(Cbar)
    psi = (mu / 2) * (Ic - 3) + 0.5 * lmbda * (J-1) * ln(J)
    P = diff(psi, F)
    rhos = Constant(1e3)
    gravity = Constant(volumetric_load) * vec_z

    cof = J * inv(F).T
    cofnorm = sqrt(dot(cof * nn, cof * nn))
    # Valid at spatial config, so we pull back the normal.
    NN = 1 / cofnorm * cof * nn

    F_form = inner(P, grad(v)) * dx_mesh - ramp * Constant(neumann_load) * \
        dot(vec_z, v) * cofnorm * ds_N - ramp * \
        rhos * dot(gravity, v) * 1 / J * dx
    return F_form
