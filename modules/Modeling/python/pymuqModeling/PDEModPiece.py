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

class PDEModPiece(mm.PyModPiece):
    def __init__(self, Vstate, Vpara, weak_form, bc_fwd, bc_adj, is_fwd_linear=True):
        """
        Args:
            param1 (Vstate): The function space for the state and adjoint
            param1 (Vpara): The function space for the parameter
            param2 (weak_form): A function that implements the weak from
            param3 (bc_fwd): Boundary conditions for the forward problem
            param4 (bc_adj): Boundary conditions for the adjoint problem
            param5 (is_fwd_linear): True- linear differential operator (default), False- nonlinear differential operator
        """
        mm.PyModPiece.__init__(self, [Vpara.dim()], [Vstate.dim()])
        self.pde = hip.PDEVariationalProblem([Vstate, Vpara, Vstate], weak_form, bc_fwd, bc_adj, is_fwd_linear=is_fwd_linear)

    def InitialGuess(self, para):
        return self.pde.generate_state()

    def EvaluateImpl(self, inputs):
        # generate the parameter
        m = self.pde.generate_parameter()
        m.set_local(inputs[0])

        # generate the soln vector
        u = self.InitialGuess(inputs[0])

        # solve the forward model
        self.pde.solveFwd(u, [u, m, None], 1.0e-9)

        # set the outputs
        self.outputs = [np.array(u)]

    def GradientImpl(self, inputDimWrt, outputDimWrt, inputs, sens):
        # generate the parameter
        m = self.pde.generate_parameter()
        m.set_local(inputs[0])

        # generate the soln vector
        u = self.InitialGuess(inputs[0])

        # solve the forward model
        self.pde.solveFwd(u, [u, m, None], 1.0e-9)

        # generate the sensitivity
        s = self.pde.generate_state()
        s.set_local(-1.0*sens)

        # generate the adjoint soln vector
        p = self.pde.generate_state()

        # solve the adjoint model
        self.pde.solveAdj(p, [u, m, None], s, 1.0e-9)

        Fm = self.pde.generate_parameter()
        self.pde.evalGradientParameter([u, m, p], Fm)

        self.gradient = Fm
