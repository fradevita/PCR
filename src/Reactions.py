from typing import NamedTuple
import math
import Species
import numpy as np

#########################################################################################
# Base Objects
#########################################################################################
# Generic Reactant
class Reactant(NamedTuple):
    ind: int   # index in the species vector
    nu: float  # Stechiometric coefficient

# Generic Product
class Product(NamedTuple):
    ind: int   # index in the species vector
    nu: float  # Stechiometric coefficient

# Generic collider
class ThirdBd(NamedTuple):
    ind: int    # index in the species vector
    eff: float  # Efficiency as a thid body

# Arrhenius coefficients
class ArrheniusCoefficients():
    def __init__(self, A: float, n:float , EovR:float):
        '''
            A: Pre-exponential factor [m^{3*(o-1)}/(mol^(o-1) s)] where o is the order fo the reaction
            n: Temperature exponent
            EovR: Activation energy [k]
        '''
        self.A = A
        self.n = n
        self.EovR = EovR

    def CompRateCoeff(self, T):
        Kf = self.A
        if self.n != 0.: Kf *= T**self.n
        if abs(self.EovR) > 1e-5: Kf *= math.exp(-self.EovR/T)
        return Kf

#########################################################################################
# Standard Reaction
#########################################################################################
class Reaction():
    def __init__(self, name : str, ArrCoeff: ArrheniusCoefficients, has_backward: bool, \
                 Neducts:int, Npducts: int, educts: Reactant, pducts: Product):
        self.name = name
        self.ArrCoeff = ArrCoeff
        self.has_backward = has_backward
        self.Neducts = Neducts
        self.Npducts = Npducts
        self.educts = educts
        self.pducts = pducts

    def CompBackwardRateCoeff(self, Kf: float, P: float, T: float, G):
        sumNu = 0.
        sumNuG = 0.
        for i in range(self.Neducts):
            sumNu -= self.educts[i].nu
            sumNuG -= self.educts[i].nu*G[self.educts[i].ind]
        for i in range(self.Npducts):
            sumNu += self.pducts[i].nu
            sumNuG += self.pducts[i].nu*G[self.pducts[i].ind]
        lnKc = - sumNuG - sumNu * ( math.log(T) + math.log(Species.RGAS/P) )   
        return Kf * math.exp(-lnKc)

    def GetReactionRate(self, P: float, T: float, C: np.ndarray, G: np.ndarray):
        # Forward reaction rate
        Kf = self.ArrCoeff.CompRateCoeff(T)
        a = 1.
        for i in range(self.Neducts):
            a *= C[self.educts[i].ind]**self.educts[i].nu
        
        # Backward reaction rate
        Kb = 0.
        b = 1.
        if self.has_backward:
            Kb = self.CompBackwardRateCoeff(Kf, P, T, G)
            for i in range(self.Npducts):
                b *= C[self.pducts[i].ind]**self.pducts[i].nu
        
        return (Kf*a - Kb*b)

    def AddProductionRates(self, P: float, T: float, C: np.ndarray, G: np.ndarray):
        R = self.GetReactionRate(P, T, C, G)
        w = np.zeros_like(C)
        for i in range(self.Neducts):
            w[self.educts[i].ind] -= self.educts[i].nu*R
        for i in range(self.Npducts):
            w[self.pducts[i].ind] += self.pducts[i].nu*R
        return w

#########################################################################################
# Third Body Reaction
#########################################################################################
class ThirdbodyReaction():
    def __init__(self, name, ArrCoeff: ArrheniusCoefficients, has_backward: bool, \
                 Neducts:int, Npducts: int, Nthird:int, educts: tuple,      \
                 pducts: tuple, thirdbd: tuple):
        self.name = name
        self.ArrCoeff = ArrCoeff
        self.has_backward = has_backward
        self.Neducts = Neducts
        self.Npducts = Npducts
        self.Nthirdb = Nthird
        self.educts = educts
        self.pducts = pducts
        self.thirdBd = thirdbd

    def CompBackwardRateCoeff(self, Kf: float, P: float, T: float, G: np.ndarray):
        sumNu = 0.
        sumNuG = 0.
        for i in range(self.Neducts):
            sumNu -= self.educts[i].nu
            sumNuG -= self.educts[i].nu*G[self.educts[i].ind]
        for i in range(self.Npducts):
            sumNu += self.pducts[i].nu
            sumNuG += self.pducts[i].nu*G[self.pducts[i].ind]
        lnKc = - sumNuG - sumNu * ( math.log(T) + math.log(Species.RGAS/P) )   
        return Kf * math.exp(-lnKc)

    def GetReactionRate(self, P: float, T: float, C: np.ndarray, G: np.ndarray):
        # Forward reaction rate
        Kf = self.ArrCoeff.CompRateCoeff(T)
        a = 1.
        for i in range(self.Neducts):
            a *= C[self.educts[i].ind]**self.educts[i].nu
        
        # Backward reaction rate
        Kb = 0.
        b = 1.
        if self.has_backward:
            Kb = self.CompBackwardRateCoeff(Kf, P, T, G)
            for i in range(self.Npducts):
                b *= C[self.pducts[i].ind]**self.pducts[i].nu
        
        # Third body efficiency
        c = 1.
        if (self.Nthirdb > 0):
            c = P/(Species.RGAS*T)
            for i in range(self.Nthirdb):
                c += C[self.thirdBd[i].ind]*(self.thirdBd[i].eff - 1.)
        return c*(Kf*a - Kb*b)

    def AddProductionRates(self, P: float, T: float, C: np.ndarray, G: np.ndarray):
        R = self.GetReactionRate(P, T, C, G)
        w = np.zeros_like(C)
        for i in range(self.Neducts):
            w[self.educts[i].ind] -= self.educts[i].nu*R
        for i in range(self.Npducts):
            w[self.pducts[i].ind] += self.pducts[i].nu*R
        return w

#########################################################################################
# Fall-off class
#########################################################################################
class Troe():
    def __init__(self, alpha: float, T1: float, T2:float, T3: float):
        self.alpha = alpha
        self.T1 = T1
        self.T2 = T2
        self.T3 = T3

#########################################################################################
# FallOf Reaction
#########################################################################################
class FalloffReaction():

    def __init__(self, name, ArrCoeffL: ArrheniusCoefficients, ArrCoeffH: ArrheniusCoefficients, \
                 has_backward: bool, Neducts: int, Npducts: int, Nthirdb: int, educts: tuple,   \
                 pducts: tuple, ThirdBd: tuple, Ftype:int, FOdata: Troe):
        self.name = name
        self.ArrCeoffL = ArrCoeffL
        self.ArrCeoffH = ArrCoeffH
        self.has_backward = has_backward
        self.Neducts = Neducts
        self.Npducts = Npducts
        self.Nthirdb = Nthirdb
        self.educts = educts
        self.pducts = pducts
        self.ThirdBd = ThirdBd
        self.Ftype = Ftype
        self.FOdata = FOdata

    def CompBackwardRateCoeff(self, Kf: float, P: float, T: float, G: np.ndarray):
        sumNu = 0.
        sumNuG = 0.
        for i in range(self.Neducts):
            sumNu -= self.educts[i].nu
            sumNuG -= self.educts[i].nu*G[self.educts[i].ind]
        for i in range(self.Npducts):
            sumNu += self.pducts[i].nu
            sumNuG += self.pducts[i].nu*G[self.pducts[i].ind]
        lnKc = - sumNuG - sumNu * ( math.log(T) + math.log(Species.RGAS/P) )   
        return Kf * math.exp(-lnKc)

    def computeF(self, P:float , T: float):
        if self.Ftype == 'F_Troe2':
            Fc = (1 - self.FOdata.alpha)*math.exp(-T/self.FOdata.T3) \
                    + self.FOdata.alpha *math.exp(-T/self.FOdata.T1)
            d = 0.14
            n = 0.75 - 1.27*math.log10(Fc)
            c = -0.4 - 0.67*math.log10(Fc)
            a = math.log10(P) + c
            f = a/(n - d*a)
            return Fc**(1.0/(1 + f*f))
        elif self.Ftype == 'F_Troe3':
            Fc = (1 - self.FOdata.alpha)*math.exp(-T/self.FOdata.T3) \
                    + self.FOdata.alpha *math.exp(-T/self.FOdata.T1) \
                                        +math.exp(-self.FOdata.T2/T)
            d = 0.14
            n = 0.75 - 1.27*math.log10(Fc)
            c = -0.4 - 0.67*math.log10(Fc)
            a = math.log10(P) + c
            f = a/(n - d*a)
            return Fc**(1.0/(1 + f*f))
        else:
            print('WRONG FTYPE')
            return 0

    def GetReactionRate(self, P: float, T: float, C: np.ndarray, G: np.ndarray):
        # Forward rate coefficient
        KfH = self.ArrCeoffH.CompRateCoeff(T)
        KfL = self.ArrCeoffL.CompRateCoeff(T)
        # Reduced pressure
        Pr = P/(Species.RGAS * T)
        if (self.Nthirdb > 0):
            for i in range(self.Nthirdb):
                Pr += C[self.ThirdBd[i].ind]*(self.ThirdBd[i].eff - 1.0)
        Pr *= KfL/max(KfH, 1e-30)
        # Use Lindemann formula
        F = self.computeF(Pr, T)
        Kf = KfH*(Pr/(1 + Pr))*F
        # Forward reaction rate
        a = 1.0
        for i in range(self.Neducts):
            a *= C[self.educts[i].ind]**self.educts[i].nu

        # Backward reaction rate
        Kb = 0.
        b = 1.
        if self.has_backward:
            Kb = self.CompBackwardRateCoeff(Kf, P, T, G)
            for i in range(self.Npducts):
                b *= C[self.pducts[i].ind]**self.pducts[i].nu

        return Kf*a - Kb*b
    
    def AddProductionRates(self, P: float, T: float, C: np.ndarray, G: np.ndarray):
        R = self.GetReactionRate(P, T, C, G)
        w = np.zeros_like(C)
        for i in range(self.Neducts):
            w[self.educts[i].ind] -= self.educts[i].nu*R
        for i in range(self.Npducts):
            w[self.pducts[i].ind] += self.pducts[i].nu*R
        return w

        
