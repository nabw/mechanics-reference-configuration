import numpy as np
from petsc4py import PETSc
from mpi4py import MPI


class AndersonAcceleration:

    def __init__(self, order):
        self.order = order
        self.k = 0
        self.F = []
        self.X = []
        self.F0 = []  # For global vectors
        self.comm = MPI.COMM_WORLD
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()
        self.is_new = True

        # Future vectors
        self.d_fk0 = None
        self.fk0 = None
        self.xk = None
        self.fk = None
        self.delta_xk = None
        self.delta_fk = None

    def reset(self):
        if not self.is_new:
            self.k = 0
            self.F = []
            self.X = []
            self.F0 = []  # For global vectors
            self.fk.zeroEntries()
            self.xk.zeroEntries()
            self.delta_xk.zeroEntries()
            self.delta_fk.zeroEntries()
            self.d_fk0.zeroEntries()
            self.fk0.zeroEntries()

    def get_next_vector(self, gk):

        if self.is_new:
            # Initialize vectors
            self.xk = gk.copy()
            self.fk = gk.copy()
            self.delta_xk = gk.copy()
            self.delta_fk = gk.copy()
            self.fk.zeroEntries()
            self.xk.zeroEntries()
            self.delta_xk.zeroEntries()
            self.delta_fk.zeroEntries()
            # Init global vectors and scatterer
            if self.size > 1:
                self.scatter, aux = PETSc.Scatter.toZero(gk)
            else:
                aux = gk.copy()
                aux.zeroEntries()
            self.d_fk0 = aux.copy()
            self.fk0 = aux.copy()
            self.is_new = False

        self.fk.copy(self.delta_fk)
        self.xk.copy(self.delta_xk)

        gk.copy(self.fk)
        self.fk.axpy(-1.0, self.xk)

        mk = min(self.k, self.order)
        if mk > 0:  # If order>0 and k>0
            self.delta_fk.aypx(-1, self.fk)
            if self.delta_fk.norm() < 1e-12:
                self.k -= 1
                gk.copy(self.xk)
            else:
                self.F.append(self.delta_fk.copy())
                if len(self.F) > self.order:
                    self.F.pop(0)

                if self.size > 1:
                    self.scatter.scatter(self.delta_fk, self.d_fk0)
                    self.scatter.scatter(self.fk, self.fk0)
                    # Process only on first core, then scatter alpha
                    if self.rank == 0:
                        self.F0.append(self.d_fk0.copy())
                        if len(self.F0) > self.order:
                            self.F0.pop(0)
                        F = np.vstack(self.F0).T
                        Q, R = np.linalg.qr(F)
                        rhs = self.fk0
                        alpha = np.linalg.solve(R, -Q.T @ rhs)
                    else:
                        alpha = None
                    alpha = self.comm.bcast(alpha, root=0)
                # If not parallel, avoid data redundancy
                else:
                    F = np.vstack(self.F).T
                    Q, R = np.linalg.qr(F)
                    rhs = self.fk
                    alpha = np.linalg.solve(R, -Q.T @ rhs)

                self.xk.axpy(1.0, self.fk)
                for i in range(mk):
                    self.xk.axpy(alpha[i], self.X[i] + self.F[i])
        else:
            gk.copy(self.xk)

        self.delta_xk.aypx(-1, self.xk)
        self.X.append(self.delta_xk.copy())
        if len(self.X) > self.order:
            self.X.pop(0)
        self.k += 1
        self.xk.copy(gk)
