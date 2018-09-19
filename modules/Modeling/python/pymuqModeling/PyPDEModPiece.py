# import hippylib
import sys
import os
sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../") )
import hippylib as hip

# import muq Modeling
import pymuqModeling_ as mm

# import dolfin/fenics
import dolfin as dfn

import numpy as np

# import matplotlib so that we can make figures
import matplotlib.pyplot as plt
from matplotlib import rcParams

class PyPDEModPiece(mm.PyModPiece):
    def __init__(self, Vh, weak_form, bc_fwd, bc_adj, is_fwd_linear=True):
        """
        Args:
            param1 (Vh): The function space for [state, parameter, adjoint]
            param2 (weak_form): A function that implements the weak from
            param3 (bc_fwd): Boundary conditions for the forward problem
            param4 (bc_adj): Boundary conditions for the adjoint problem
            param5 (is_fwd_linear): True- linear differential operator (default), False- nonlinear differential operator
        """
        mm.PyModPiece.__init__(self, [Vh[1].dim()], [Vh[0].dim()])
        self.pde = hip.PDEVariationalProblem(Vh, weak_form, bc_fwd, bc_adj, is_fwd_linear=is_fwd_linear)

        self.Vh = Vh

        print()
        print('bc0: ', self.pde.bc0)
        print()
        print()

    def EvaluateImpl(self, inputs):
        # generate the soln vector
        u = self.pde.generate_state()

        # generate the parameter
        m = self.pde.generate_parameter()
        m.set_local(inputs[0])

        # solve the forward model
        self.pde.solveFwd(u, [None, m, None], 1.0e-9)

        # set the outputs
        self.outputs = [np.array(u)]

    def GradientImpl(self, inputDimWrt, outputDimWrt, inputs, sens):
        # generate the parameter
        m = self.pde.generate_parameter()
        m.set_local(inputs[0])

        # generate the soln vector
        u = self.pde.generate_state()

        # solve the forward model
        self.pde.solveFwd(u, [None, m, None], 1.0e-9)

        print()
        print()
        print()
        print('forward soln:', u)
        print(np.array(u))
        print()
        print()

        # generate the sensitivity
        s = self.pde.generate_state()
        s.set_local(-1.0*sens)

        # generate the adjoint soln vector
        p = self.pde.generate_state()

        # solve the adjoint model
        self.pde.solveAdj(p, [u, m, None], s, 1.0e-9)

        # w = dfn.TestFunction(self.Vh[2])
        # adj = dfn.TrialFunction(self.Vh[2])
        # rhs = dfn.Function(self.Vh[2])
        # rhs.vector()[:] = s
        # kap = dfn.Function(self.Vh[1])
        # kap.vector()[:] = m
        # L = rhs*w*dfn.dx
        # a = kap*dfn.inner(dfn.grad(adj), dfn.grad(w))*dfn.dx
        #
        # # define the Dirichlet boundary
        # def boundary(x, on_boundary):
        #     return on_boundary and x[0]<dfn.DOLFIN_EPS
        #
        # # define the boundary condition for the adjoint problem
        # p0 = dfn.Constant(0.0)
        # bc_adj = dfn.DirichletBC(self.Vh[2], p0, boundary)
        #
        # adjSoln = dfn.Function(self.Vh[2])
        # dfn.solve(a==L, adjSoln, bcs=bc_adj)
        # p = adjSoln.vector()

        print()
        print()
        print('adjoint soln:', p)
        print(np.array(p))
        print()
        print()

        Fm = self.pde.generate_parameter()
        self.pde.evalGradientParameter([u, m, p], Fm)

        print()
        print()
        print('grad of F:', Fm)
        print(np.array(Fm))

        Fgrad = dfn.Function(self.Vh[1])
        Fgrad.vector().set_local(np.array(Fm))

        #plt.figure(figsize=(15,5))
        #dfn.plot(Fgrad)
        #plt.show()

        #adj = dfn.Function(self.Vh[2])
        #adj.vector().set_local(np.array(p))

        #w = dfn.TestFunction(self.Vh[2])
        #f = dfn.TrialFunction(self.Vh[2])
        #L = -w*adj*Fgrad*dfn.dx
        #a = w*f*dfn.dx

        #fs = dfn.Function(self.Vh[2])
        #dfn.solve(a==L, fs)

        #print()
        #print(np.array(fs.vector()[:]))
        #print()
        #print()

        #dum = dfn.assemble(dfn.dot(dfn.grad(fwd), dfn.grad(adj))*dfn.dx)
        #dum = dfn.inner(dfn.grad(fwd), dfn.grad(adj))

        #print(dum)

        self.gradient = p
