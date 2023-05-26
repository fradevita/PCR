from typing import NamedTuple
import numpy as np
import math

# CONSTANTS
RGAS = 8.3144598     # [J/(mol K)]
ATom = 1e-10         # Angstrom to meter
DToCm = 3.33564e-30  # Debye to Coulomb meter
Na = 6.02214086e23   # [1/mol]
kb = 1.38064852e-23  # [m^2 kg /( s^2 K)]

def omega_mu(T: float):
    '''omega_mu(T):
        returns the collision integral for mu given dimensionless temperature t/(eps/k)
    '''
    m1 = 3.3530622607
    m2 = 2.53272006
    m3 = 2.9024238575
    m4 = 0.11186138893
    m5 = 0.8662326188       
    m6 = 1.3913958626
    m7 = 3.158490576
    m8 = 0.18973411754
    m9 = 0.00018682962894

    num = m1 + T*(m2 + T*(m3 + T*m4))
    den = m5 + T*(m6 + T*(m7 + T*(m8 + T*m9)))
    return num / den

class cpCoefficients(NamedTuple):
    TSwitch1: float     # Switch  temperature between Low and Mid  temperature polynomials
    TSwitch2: float     # Switch  temperature between Mid and High temperature polynomials
    Tmin: float         # Minimum Temperature
    TMax: float         # Maximum temperature
    cpH: np.ndarray     # High temperature polynomials
    cpM: np.ndarray     # Mid  temperature polynomials
    cpL: np.ndarray     # Low  temperature polynomials

class DiffCoefficients(NamedTuple):
    sigma: float        # Lennard-Jones collision diameter [m]
    kbOveps: float      # Boltzmann constant divided by Lennard-Jones potential well depth [1/K]
    mu: float           # Dipole moment [C*m]
    alpha: float        # Polarizabilty [m]
    Z298: float         # Rotational relaxation collision number

class Spec():
    def __init__(self, Name: str, W: float, index: int, Geom: int, cpCoeff: cpCoefficients, DiffCoeff: DiffCoefficients):
        self.Name = Name
        self.W = W
        self.index = index
        self.Geom = Geom
        self.cpCoeff = cpCoeff
        self.DiffCoeff = DiffCoeff

    def GetEnthalpy(self, T: float):
        rOvW = RGAS/self.W
        Tinv = 1.0/T
        cpCoeff = self.cpCoeff.cpH if T > self.cpCoeff.TSwitch2 else \
                (self.cpCoeff.cpM if T > self.cpCoeff.TSwitch1 else self.cpCoeff.cpL)
        E = -cpCoeff[0]*Tinv + cpCoeff[1]*math.log(T) + cpCoeff[7]      + T* \
                                                    ( cpCoeff[2]      + T* \
                                                    ( cpCoeff[3]*0.50 + T* \
                                                    ( cpCoeff[4]/3    + T* \
                                                    ( cpCoeff[5]*0.25 + cpCoeff[6]/5*T))))
        return E*rOvW

    def GetCp(self, T: float):
        rOvW = RGAS/self.W
        Tinv = 1.0/T
        cpCoeff = self.cpCoeff.cpH if T > self.cpCoeff.TSwitch2 else \
                (self.cpCoeff.cpM if T > self.cpCoeff.TSwitch1 else self.cpCoeff.cpL)
        return rOvW*(cpCoeff[0]*Tinv*Tinv + cpCoeff[1]*Tinv +  cpCoeff[2] + T* \
                                                            ( cpCoeff[3] + T* \
                                                            ( cpCoeff[4] + T* \
                                                            ( cpCoeff[5] + T*cpCoeff[6]))))

    def GetMu(self, T: float):
        num = 5*math.sqrt(math.pi*self.W/Na*kb*T)
        den = 16*math.pi*(self.DiffCoeff.sigma**2)*omega_mu(T*self.DiffCoeff.kbOveps)
        return num/den

    def GetFreeEnthalpy(self, T):
        Tinv = 1./T
        logT = math.log(T)
        cpC = self.cpCoeff.cpH if T > self.cpCoeff.TSwitch2 else \
             (self.cpCoeff.cpM if T > self.cpCoeff.TSwitch1 else self.cpCoeff.cpL)
        G = -0.5*cpC[0]*Tinv**2 + cpC[1]*Tinv*(1.0 + logT) + cpC[2]*(1.0 - logT) + cpC[7]*Tinv - cpC[8]
        return G - 0.5*T*( cpC[3]   + T*
                    ( cpC[4]/3 + T*
                    ( cpC[5]/6 + 0.1*T*cpC[6] )))