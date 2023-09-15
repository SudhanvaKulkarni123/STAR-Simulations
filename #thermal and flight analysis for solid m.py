#thermal and flight analysis for solid motor rockets
import numpy as np
import math
import CoolProp.CoolProp as CP
from testing import test

#tank dims
V_oxtank = 6.92655*(10**-3) #m^3
A_oxtank = math.pi * (0.05715**2) #m^2
V_ethtank = 7.57*(10**-3) #m^3
A_ethtank = math.pi * (0.070612**2) #m^2


def Collapse_OX(P_oxtank_psi,h, V_oxinit = 3.25):    
    P_oxtank = P_oxtank_psi * 6895 #convert to Pa

    #pressurant
    V_oxgas = V_oxtank - V_oxinit*(10**-3) #m^3
    MN2 = CP.PropsSI('M', 'Nitrogen') #kg/mol
    RN2_ox_molar = CP.PropsSI('GAS_CONSTANT', 'V', V_oxgas, 'T', 298.15, 'Nitrogen')/1000
    RN2_ox = RN2_ox_molar/MN2
    mN2_ox = (P_oxtank/1000 * V_oxgas) / (RN2_ox * 298.15)
    CpN2_ox = CP.PropsSI('CPMASS', 'P', P_oxtank, 'T', 298.15, 'Nitrogen')/1000

    # Diffeq conditions
    T0 = 90.2
    Ta = 300
    tSpan = np.linspace(0, 15, int(15/h)) #to change time range edit this
    T_new = np.zeros((len(tSpan), 2))
    T_new[0, :] = [T0, Ta]
    delta_T = [Ta - T0]

    # Constants from Ring
    C = 0.27
    n = 1/4

    CpLOX = CP.PropsSI('CPMASS', 'P', P_oxtank, 'T', 90.2, 'Oxygen')/1000
    CvLOX = CP.PropsSI('CVMASS', 'P', P_oxtank, 'T', 90.2, 'Oxygen')/1000
    muLOX = CP.PropsSI('VISCOSITY', 'P', P_oxtank, 'T', 90.2, 'Oxygen')/1000
    kfLOX = CP.PropsSI('CONDUCTIVITY', 'P', P_oxtank, 'T', 90.2, 'Oxygen')/1000
    rhoLOX = CP.PropsSI('DMASS', 'P', P_oxtank, 'T', 90.2, 'Oxygen')/1000
    betaLOX = 1 / CP.PropsSI('T', 'P', P_oxtank, 'Q', 0, 'Oxygen')  # coefficient of volume expansion
    mOX = rhoLOX * V_oxinit #kg
    LsLOX = V_oxinit/A_oxtank #characteristic length

    # Prandtl
    PrLOX = CpLOX * muLOX / kfLOX
    # Grashof
    GrLOX = (LsLOX**3) * rhoLOX * betaLOX * abs(delta_T[0]) / (muLOX**2)
    # Heat transfer coeff
    hLOX = [C * kfLOX * ((GrLOX * PrLOX)**n) / LsLOX]

    # Differential equation
    def convecLOX(X, h_LOX):
        alpha = (-mOX * CvLOX) / (mN2_ox * CpN2_ox)
        return np.array([alpha * h_LOX * (X[0] - X[1]), -h_LOX * (X[1] - X[0])])

    # RK4 integration
    for j in range(len(tSpan) - 1):
        k1 = convecLOX(T_new[j, :], hLOX[j])
        k2 = convecLOX(T_new[j, :] + 0.5 * h * k1, hLOX[j])
        k3 = convecLOX(T_new[j, :] + 0.5 * h * k2, hLOX[j])
        k4 = convecLOX(T_new[j, :] + h * k3, hLOX[j])
        T_new[j+1, :] = T_new[j, :] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
        delta_T.append(T_new[j+1, 1] - T_new[j+1, 0])
        GrLOX = (LsLOX**3) * rhoLOX * betaLOX * (delta_T[-1]) / (muLOX**2)
        hLOX.append(C * kfLOX * ((GrLOX * PrLOX)**n) / LsLOX)

    return T_new[:,1]


test.check_consistency(lambda x: Collapse_OX(450,x), 0.5, 4)