import numpy as np
import math
import Species

# Definition of Mixture class
class Mixture:
    '''Mixture of gas'''
    def __init__(self, nSpec: int, Species: tuple, Species_index: dict, \
                 nReacts: int, Reacts: tuple, nTBReacts: int, TBReacts: tuple, \
                 nFOReacts: int, FalloffReacts: tuple):
        self.nSpec = nSpec
        self.Species = Species
        self.Species_index = Species_index
        assert(nSpec == len(Species))
        self.nReacts = nReacts
        self.Reacts = Reacts
        assert(nReacts == len(Reacts))
        self.nTBReacts = nTBReacts
        self.ThirdBodyReacts = TBReacts
        assert(nTBReacts == len(TBReacts))
        self.nFOReacts = nFOReacts
        self.FalloffReacts = FalloffReacts
        assert(nTBReacts == len(TBReacts))
    
    def GetMolarWeightFromYi(self, Yi: np.ndarray):
        MixW = 0.0
        for i in range(self.nSpec):
            MixW += Yi[i] / self.Species[i].W
        return 1.0/MixW

    def GetMolarWeightFromXi(self, Xi: np.ndarray):
        MixW = 0.0
        for i in range(self.nSpec):
            MixW += Xi[i]*self.Species[i].W
        return MixW

    def GetMolarFractions(self, Yi: np.ndarray):
        MixW = self.GetMolarWeightFromYi(Yi)
        Xi = np.zeros(self.nSpec)
        for i in range(self.nSpec):
            Xi[i] = Yi[i]*MixW/self.Species[i].W
        return Xi

    def GetMassFractions(self, MixW: float, Xi: np.ndarray):
        Yi = np.zeros(self.nSpec)
        for i in range(self.nSpec):
            Yi[i] = Xi[i]*self.Species[i].W/MixW
        return Yi

    def GetYi(self, rho: float, rhoYi: np.ndarray):
        irho = 1.0/rho
        Yi = np.zeros_like(rhoYi)
        for i in range(self.nSpec):
            Yi[i] = irho*rhoYi[i]
        return Yi

    def GetRho(self, P, T, MixW):
        return P*MixW/(Species.RGAS*T)

    def GetHeatCapacity(self, T: float, Yi: np.ndarray):
        cp = 0.0
        for i in range(self.nSpec):
            cp += Yi[i]*self.Species[i].GetCp(T)
        return cp

    def ClampYi(self, Yi):
        for i in range(len(Yi)):
            Yi[i] = max(Yi[i], 1.0e-30)
            Yi[i] = min(Yi[i], 1.0)
        return Yi

    def GetEnthalpy(self, T: float, Yi: np.ndarray):
        Enth = 0.0
        for i in range(self.nSpec):
            Enth += Yi[i]*Species.GetEnthalpy(self.Species[i], T)
        return Enth

    def GetInternalEnergy(self, T: float, Yi: np.ndarray):
        e = 0.0
        for i in range(self.nSpec):
            e += Yi[i]*(self.Species[i].GetEnthalpy(T) - Species.RGAS*T/self.Species[i].W)
        return e

    def GetTFromRhoAndP(self, rho, MixW, P):
        return P*MixW/(rho*Species.RGAS)

    def GetGamma(self, T: float, Yi: np.ndarray):
        cp = 0.0
        cv = 0.0
        for i in range(self.nSpec):
            cp_s = self.Species[i].GetCp(T)
            cp += Yi[i]*cp_s
            cv += Yi[i]*(cp_s - Species.RGAS/self.Species[i].W)
        return cp/cv

    def GetSpeedOfSound(self, T: float, Yi: np.ndarray):
        MixW = self.GetMolarWeightFromYi(Yi)
        gamma = self.GetGamma(T, Yi) 
        return math.sqrt(gamma * Species.RGAS * T / MixW)

    def GetViscosity(self, T: float, Xi: np.ndarray):
        muk = np.zeros(self.nSpec)
        for i in range(self.nSpec):
            muk[i] = self.Species[i].GetMu(T)

        mu = 0.0
        for i in range(self.nSpec):
            den = 0.0
            for j in range(self.nSpec):
                Phi = (1.0 + math.sqrt(muk[i]/muk[j])*(self.Species[j].W/self.Species[i].W)**(1.0/4.0) )**2
                Phi /= math.sqrt(8)*math.sqrt(1.0 + self.Species[i].W/self.Species[j].W)
                den += Xi[j]*Phi
            mu += Xi[i]*muk[i]/den
        return mu
    
    def GetDiffusion(self, T: float, P: float, MixW: float, X: np.ndarray):
        D = np.zeros(self.nSpec)
        Y = self.GetMassFractions(MixW, X)
        for i in range(self.nSpec):
            num = 1. - Y[i]
            den = 0.
            for j in range(self.nSpec):
                if (j != i):
                    Dij = Species.GetBinaryDiffusivity(self.Species[i], self.Species[j], P, T)
                    den += X[j]/Dij
            D[i] = num/den
        return D

    
    def GetTFromInternalEnergy(self, e0: float, T: float, Yi: np.ndarray):
        MAXITS = 1000
        TOL = 1e-8
        dfdT = 1.
        j = 0
        while j < MAXITS:
            f = e0 - self.GetInternalEnergy(T, Yi)
            if abs(f/dfdT) < TOL: break
            dfdT = 0.
            for i in range(self.nSpec):
                dfdT += Yi[i]*(self.Species[i].GetCp(T) - Species.RGAS/self.Species[i].W)#*iCpRef;
            T += f/dfdT
            j += 1
        assert(j < MAXITS)
        return T

    def GetHeatConductivity(self, T: float, Xi: np.ndarray):
        a = 0.0
        b = 0.0
        for i in range(self.nSpec):
            lami = self.Species[i].GetLam(T)
            a += Xi[i]*lami
            b += Xi[i]/lami
        return 0.5*(a + 1.0/b)
 
    def GetProductionRates(self, rhon, Pn, Tn, Yi):
        # Use unscaled primitive variables
        T = Tn
        P = Pn
        rho = rhon

        w = np.zeros_like(Yi)
        C = np.zeros_like(w)
        G = np.zeros_like(w)
        for i in range(self.nSpec):
            C[i] = Yi[i]*rho/self.Species[i].W
            G[i] = self.Species[i].GetFreeEnthalpy(T)

        wd = np.zeros_like(w)
        if self.nReacts > 0:
            for i in range(self.nReacts):
                wd = self.Reacts[i].AddProductionRates(P, T, C, G)
                w += wd

        if self.nTBReacts > 0:
            for i in range(self.nTBReacts):
                wd = self.ThirdBodyReacts[i].AddProductionRates(P, T, C, G)
                w += wd
        
        if self.nFOReacts > 0:
            for i in range(self.nFOReacts):
                wd = self.FalloffReacts[i].AddProductionRates(P, T, C, G)
                w += wd

        for i in range(self.nSpec):
            w[i] *= self.Species[i].W

        return w
