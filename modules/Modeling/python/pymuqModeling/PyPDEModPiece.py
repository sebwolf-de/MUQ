# import hippylib
import sys
import os
sys.path.append( os.environ.get('HIPPYLIB_BASE_DIR', "../") )
import hippylib as hip

# import muq Modeling
import pymuqModeling_ as mm

import numpy as np

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
        self.pde = hip.PDEVariationalProblem(Vh, weak_form, bc_fwd, bc_adj, is_fwd_linear=True)

    def EvaluateImpl(self, inputs):
        # generate the soln vector
        u = self.pde.generate_state()

        # generate the parameter
        m = self.pde.generate_parameter()
        m.set_local(inputs[0])

        # solve the forward model
        self.pde.solveFwd(u, [u, m, None], 1.0e-9)

        # set the outputs
        self.outputs = [np.array(u)]

    def GradientImpl(self, inputDimWrt, outputDimWrt, inputs, sens):
        # generate the soln vector
        u = self.pde.generate_state()

        # generate the parameter
        m = self.pde.generate_parameter()
        m.set_local(inputs[0])

        # solve the forward model
        #self.pde.solveFwd(u, [u, m, None], 1.0e-9)

        # generate the sensitivity
        s = self.pde.generate_state()
        s.set_local(sens)

        # solve the adjoint model
        self.pde.solveAdj(u, [u, m, None], s, 1.0e-9)

        self.gradient = u
