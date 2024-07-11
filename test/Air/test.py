import sys
import os
import numpy as np
import math

# Import my python scripts
sys.path.insert(0, os.path.expandvars("$PCR/src/"))
import Air
import Mixture
import Species

# Create the mixture
Mix = Air.CreateMixture()

# Check the index progression of species
for i, spec in enumerate(Mix.Species):
    assert(i == spec.index)

# Define time marching scheme
def advance(rhoi_n: np.ndarray, P_n: float, T_n: float, Yi_n: np.ndarray, Mix: Mixture, dt):
    # Compute production rates at timestep n 
    rho_n = np.sum(rhoi_n)
    wd = Mix.GetProductionRates(rho_n, P_n, T_n, Yi_n)
    
    # Advance partial densities
    rhoi_np1 = rhoi_n + dt*wd
    
    # New mixture density
    rho_np1 = np.sum(rhoi_np1)
    
    # New Mass fractions
    Yi_np1 = Mix.GetYi(rho_np1, rhoi_np1)
    #Yi_np1 = Mix.ClampYi(Yi_np1)
    
    # New temperature
    T_np1 = Mix.GetTFromInternalEnergy(e0, T_n, Yi_np1)
    
    # New pressure
    P_np1 = T_np1*rho_np1*Species.RGAS/Mix.GetMolarWeightFromYi(Yi_np1)
    
    return rhoi_np1, P_np1, T_np1, Yi_np1

####################################################################################################
# Test setup 
####################################################################################################
P = 20173.3  # [Pa]
T = 7000.0   # [K]
Yi = np.zeros(Mix.nSpec) + 1.0e-50
Yi[Mix.Species_index['O2']] = 0.233
Yi[Mix.Species_index['N2']] = 0.767
Xi = Mix.GetMolarFractions(Yi)
W = Mix.GetMolarWeightFromXi(Xi)
rho = Mix.GetRho(P, T, W)
e0 = Mix.GetInternalEnergy(T, Yi)
rhoi = rho*Yi

dt = 4.0e-9
time = []
O2t = []
N2t = []
Ot = []
Nt = []
NOt = []
temperature = []

Tend = 3.0e-3
for n in range(4000000):
    # Solve
    rhoi, P, T, Yi = advance(rhoi, P, T, Yi, Mix, dt)
    time.append(dt*(n+1))
    print(f't: {time[n]:16.8e}, P: {P:16.8e}, T: {T:16.8e}, rho: {rho:16.8e}')
    
    # Store solution
    N2t.append(Yi[Mix.Species_index['N2']])
    O2t.append(Yi[Mix.Species_index['O2']])
    Ot.append(Yi[Mix.Species_index['O']])
    Nt.append(Yi[Mix.Species_index['N']])
    NOt.append(Yi[Mix.Species_index['NO']])
    temperature.append(T)
    
    if time[n] >= Tend: break

with open('PCR.csv', 'w') as fp:
    fp.write('%s,%s,%s,%s,%s,%s,%s\n' % ('time','N2','O2','O','N','NO','T'))
    for i in range(0,len(N2t)):
        fp.write('%16.8e,%16.8e,%16.8e,%16.8e,%16.8e,%16.8e,%16.8e\n' %  
                 (time[i], N2t[i], O2t[i], Ot[i], Nt[i], NOt[i], temperature[i]))