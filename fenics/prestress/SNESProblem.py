from dolfin import *
from petsc4py import PETSc
from mpi4py import MPI


class SNESProblem():
    def __init__(self, FF, JJ, uu, bbcs=None):
        self.L = FF
        self.a = JJ  # derivative(FF, uu, ddu)
        self.u = uu
        self.bcs = None
        if bbcs:
            self.bcs = bbcs
            self.bcs0 = [DirichletBC(_bc) for _bc in bbcs]
            for bc in self.bcs0:
                bc.homogenize()

    def F(self, snes, xx, FF):
        xx = PETScVector(xx)
        FF = PETScVector(FF)
        xx.vec().copy(self.u.vector().vec())
        self.u.vector().apply("")
        assemble(self.L, tensor=FF)
        if self.bcs:
            for bc in self.bcs0:
                bc.apply(FF)

    def J(self, snes, xx, JJ, PP):
        JJ = PETScMatrix(JJ)
        JJ.mat().setBlockSize(3)
        xx.copy(self.u.vector().vec())
        self.u.vector().apply("")
        assemble(self.a, tensor=JJ)
        if self.bcs:
            for bc in self.bcs0:
                bc.apply(JJ)
        PP = JJ

    def solve(self):
        b = PETScVector()  # same as b = PETSc.Vec()
        J_mat = PETScMatrix()
        assemble(self.L, tensor=b)
        assemble(self.a, tensor=J_mat)
        snes = PETSc.SNES().create(MPI.COMM_WORLD)
        snes.setFunction(self.F, b.vec())
        J_mat.mat().setBlockSize(3)
        snes.setJacobian(self.J, J_mat.mat())
        snes.setFromOptions()
        snes.ksp.setFromOptions()
        if self.bcs:
            for bc in self.bcs:
                bc.apply(self.u.vector())
        snes.solve(None, self.u.vector().vec())
        return snes.getIterationNumber(), snes.getLinearSolveIterations()
