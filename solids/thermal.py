
import cantera as ct
import math
from dolfin import *
from fenics import *
import time
import numpy as np



# def equilSoundSpeeds(gas, rtol=1.0e-6, max_iter=5000):
#     """
#     Returns a tuple containing the equilibrium and frozen sound speeds for a
#     gas with an equilibrium composition.  The gas is first set to an
#     equilibrium state at the temperature and pressure of the gas, since
#     otherwise the equilibrium sound speed is not defined.
#     """

#     # set the gas to equilibrium at its current T and P
#     gas.equilibrate('TP', rtol=rtol, max_iter=max_iter)

#     # save properties
#     s0 = gas.s
#     p0 = gas.P
#     r0 = gas.density

#     # perturb the pressure
#     p1 = p0*1.0001

#     # set the gas to a state with the same entropy and composition but
#     # the perturbed pressure
#     gas.SP = s0, p1

#     # frozen sound speed
#     afrozen = math.sqrt((p1 - p0)/(gas.density - r0))

#     # now equilibrate the gas holding S and P constant
#     gas.equilibrate('SP', rtol=rtol, max_iter=max_iter)

#     # equilibrium sound speed
#     aequil = math.sqrt((p1 - p0)/(gas.density - r0))

#     # check against the built-in sound speed function
#     afrozen2 = gas.sound_speed

#     return aequil, afrozen, afrozen2


# # test program
# if __name__ == "__main__":
#     gas = ct.Solution('gri30.yaml')
#     gas.X = 'CH4:1.00, O2:2.0, N2:7.52'
#     T_range = np.arange(300, 2901, 100)
#     for T in T_range:
#         gas.TP = T, ct.one_atm
#         print(T, equilSoundSpeeds(gas))

class Body:
    def __init__(self) -> None:
        pass

    def set_source(self) -> None:
        pass

    def set_boundary(self) -> None:
        pass

nx = ny = 30
mesh = RectangleMesh(Point(-2,-2), Point(2,2),nx,ny)
V = FunctionSpace(mesh, 'P', 1)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, Contant(0),boundary)
 
 u_0 = Expression('exp(-a*pow(x[0],2) - -a*pow(x[1],2))', degree = 2, a = 5)

u_n = interpolate(u_0, V)

u = Trialfunction(V)
v = TestFunction(V)
f = Constant(0)

