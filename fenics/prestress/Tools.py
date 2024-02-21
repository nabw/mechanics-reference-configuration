import numpy as np
from .Enums import *


def parprint(*args):
    from dolfin import MPI
    rank = MPI.rank(MPI.comm_world)
    if rank == 0:
        print("[P-S]", *args, flush=True)


def printProblemSetting(exp_init, exp_end, n_alphas, alpha_dynamic, n_steps_dynamic, u_deg, f_deg, rtol, pressure_init, pressure_known, delta_pressure):

    # Control problem meta parameters
    u_deg = 1  # Displacement degree
    f_deg = 1  # Loads (distributed and surface) degree
    parprint(
        "==========================================================================")
    parprint("=== Prestress solver")
    parprint("--- Model: ")
    parprint("    displacement degree={}, load degree={}, relative tolerance={}".format(
        u_deg, f_deg, rtol))
    parprint("--- alpha ramp: ")
    parprint("    Initial exponent={}, last exponent={}, number of log-steps={}".format(
        exp_init, exp_end, n_alphas))
    parprint("    Dynamic alpha={:.1e}, dynamic steps={} ".format(
        alpha_dynamic, n_steps_dynamic))
    parprint("--- Pressure ramp: ")
    parprint("    Initial pressure={:.1e}, known pressure={:.1e}, delta p".format(
        pressure_init, pressure_known, delta_pressure))
    parprint(
        "==========================================================================")


def getXDMFFile(outname, mesh=False):
    import dolfin as df
    from dolfin import MPI
    import sys
    if not outname:
        return None
    if mesh:
        if ".xdmf" not in outname[:-5]:
            outname = outname + "_mesh.xdmf"
    else:
        if ".xdmf" not in outname[:-5]:
            outname = outname + ".xdmf"
    xdmffile = df.XDMFFile(MPI.comm_world, "{}".format(outname))
    xdmffile.parameters["functions_share_mesh"] = True
    xdmffile.parameters["flush_output"] = True
    #xdmffile.parameters["rewrite_function_mesh"] = False
    return xdmffile


def slabGeometry(nx, ny, nz):
    import dolfin as df
    mesh = None
    mesh = df.BoxMesh(df.Point(0.0, 0.0, 0.0),
                      df.Point(1e-2, 3e-3, 3e-3), nx, ny, nz)
    markers = df.MeshFunction("size_t", mesh, mesh.geometric_dimension() - 1)
    X0 = df.CompiledSubDomain("near(x[0], 0.0) && on_boundary")
    XL = df.CompiledSubDomain("near(x[0], 1e-2) && on_boundary")
    Y0 = df.CompiledSubDomain("near(x[1], 0.0) && on_boundary")
    YL = df.CompiledSubDomain("near(x[1], 3e-3) && on_boundary")
    Z0 = df.CompiledSubDomain("near(x[2], 0.0) && on_boundary")
    ZL = df.CompiledSubDomain("near(x[2], 3e-3) && on_boundary")
    markers.set_all(0)
    X0.mark(markers, 1)
    Y0.mark(markers, 2)
    Z0.mark(markers, 3)
    XL.mark(markers, 4)
    YL.mark(markers, 5)
    ZL.mark(markers, 6)
    return mesh, markers, 1, 2, 3, 4, 5, 6


def prolateGeometry(filename="prolate_4mm"):
    import dolfin as df

    mesh_dir = "meshes/"
    xdmf_meshfile = mesh_dir + filename + ".xdmf"
    xdmf_meshfile_bm = mesh_dir + filename + "_bm.xdmf"

    mesh = None
    mesh = df.Mesh()
    with df.XDMFFile(xdmf_meshfile) as infile:
        infile.read(mesh)
    mvc = df.MeshValueCollection("size_t", mesh, 2)
    with df.XDMFFile(xdmf_meshfile_bm) as infile:
        infile.read(mvc, "name_to_read")
    from dolfin import cpp
    markers = df.cpp.mesh.MeshFunctionSizet(mesh, mvc)
    ENDO, EPI, BASE = 20, 10, 50

    df.MeshTransformation.scale(mesh, 1e-3)
    mesh.coordinates()[:, 2] = -mesh.coordinates()[:, 2]  # They are inverted
    return mesh, markers, ENDO, EPI, BASE
