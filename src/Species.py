from typing import NamedTuple
import numpy as np
import math

# CONSTANTS
RGAS = 8.3144598          # [J/(mol K)]
ATom = 1e-10              # Angstrom to meter
DToCm = 3.33564e-30       # Debye to Coulomb meter
Na = 6.02214086e23        # [1/mol]
kb = 1.38064852e-23       # [m^2 kg /( s^2 K)]
eps_0 =  8.8541878128e-12 # [F/m] or [C/(V m)]

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

def omega_D(T: float):
   m1 = 6.8728271691
   m2 = 9.4122316321
   m3 = 7.7442359037
   m4 = 0.23424661229
   m5 = 1.45337701568
   m6 = 5.2269794238
   m7 = 9.7108519575
   m8 = 0.46539437353
   m9 = 0.00041908394781

   num = m1 + T * (m2 + T * (m3 + T * m4))
   den = m5 + T * (m6 + T * (m7 + T * (m8 + T * m9)))
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
    

    def GetSelfDiffusion(self, T):
        num = 3*np.sqrt(np.pi*kb*T*self.W/Na)
        den = 8*np.pi*self.DiffCoeff.sigma**2*omega_D(T * self.DiffCoeff.kbOveps)
        return num/den

    def GetFZrot(self, T):
        tmp = 1.0/(self.DiffCoeff.kbOveps*T)
        return 1. + 0.5*np.pi**1.5*np.sqrt(tmp) + (2.0 + 0.25*np.pi**2)*tmp + np.pi**1.5*tmp**1.5  

    def GetLamAtom(self, T):
        return 15.0/4.0*self.GetMu(T)*RGAS/self.W

    def GetLamLinear(self, T):
        CvTOvR = 1.5
        CvROvR = 1.0

        CvT = CvTOvR*RGAS
        CvR = CvROvR*RGAS
        CvV = self.GetCp(T)*self.W - 3.5*RGAS

        Dkk = self.GetSelfDiffusion(T)
        mu = self.GetMu(T)

        fV = Dkk/mu

        Zrot = self.DiffCoeff.Z298*self.GetFZrot(298)/self.GetFZrot(T)

        A = 2.5 - fV
        B = Zrot + 2/np.pi*(5./3*CvROvR+fV)

        fT = 2.5 * (1. - 2*CvR*A/(np.pi*CvT*B))
        fR = fV*(1. + 2*A/(np.pi*B))

        return mu/self.W*(fT*CvT + fR*CvR + fV*CvV)

    def GetLam(self, T):
        return self.GetLamAtom(T) if self.Geom == 0 else (self.GetLamLinear(T) if self.Geom == 1 else self.GetLamNonLinear(T))

def GetDiffCollParam_Stock(Si: Spec, Sj: Spec, T: float):

    xi = 1.
    if ((Si.DiffCoeff.mu*Sj.DiffCoeff.mu == 0.0) and \
        (Si.DiffCoeff.mu + Sj.DiffCoeff.sigma == 0)):
        
        if Si.DiffCoeff.mu != 0.:
            # Si is the polar molecule and Sj is non-polar
            mup = Si.DiffCoeff.mu/math.sqrt(4.*math.pi*eps_0*kb*Si.DiffCoeff.sigma**3/Si.DiffCoeff.kbOveps)
            alp = Sj.DiffCoeff.alpha/Sj.DiffCoeff.sigma**3
            epr = math.sqrt(Sj.DiffCoeff.kbOveps/Si.DiffCoeff.kbOveps)
        else:
            # Si is the non-polar molecule and Sj is polar
            mup = Sj.DiffCoeff.mu/math.sqrt(4.*math.pi*eps_0*kb*Sj.DiffCoeff.sigma**3/Sj.DiffCoeff.kbOveps)
            alp = Si.DiffCoeff.alpha/Si.DiffCoeff.sigma**3
            epr = math.sqrt(Si.DiffCoeff.kbOveps/Sj.DiffCoeff.kbOveps)
        xi = 1. + 0.25*mup*alp*epr

    # Binary cross-section
    sigmaij = 0.5*(Si.DiffCoeff.sigma + Sj.DiffCoeff.sigma)*xi**(-1./6.)
    
    # Collision integral
    kboEpsij = math.sqrt(Si.DiffCoeff.kbOveps * Sj.DiffCoeff.kbOveps)/(xi*xi)
    Omega11ij = omega_D(T * kboEpsij)
    return sigmaij, Omega11ij

def GetBinaryDiffusivity(Si: Spec, Sj: Spec, P: float, T: float):
    Mij = Si.W*Sj.W/(Si.W + Sj.W)
    num = 3.*math.sqrt(2.*math.pi*Na*(kb*T)**3/Mij)
    sigmaij, Omegaij = GetDiffCollParam_Stock(Si, Sj, T)
    den = 16.*P*math.pi*sigmaij**2*Omegaij
    return num/den 
