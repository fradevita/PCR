import sys
import os
import numpy as np
import math

# Import my python scripts
sys.path.insert(0, os.path.expandvars("$PCR/src/"))
import Boivin
import Mixture
import Species

# Create the mixture
Mix = Boivin.CreateMixture()

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
    rho_np1 = np.sum(rhoi)
    # New Mass fractions
    Yi_np1 = Mix.GetYi(rho_np1, rhoi_np1)
    Yi_np1 = Mix.ClampYi(Yi_np1)
    # New temperature
    T_np1 = Mix.GetTFromInternalEnergy(e0, T_n, Yi_np1)
    # New pressure
    P_np1 = T_np1*rho_np1*Species.RGAS/Mix.GetMolarWeightFromYi(Yi_np1)
    return rhoi_np1, P_np1, T_np1, Yi_np1

##################################################
# Test combustion at P = 1atm T = 1200 K phi = 1
##################################################
P = 1.0e+5   # [Pa]
T = 1200     # [K]
phi = 1.0    # equivalente ratio

# Compute Mixture composition from equivalence ratio phi (H2 mass / O2 mass)
# Molar mass of species
H2_mm = Mix.Species[Mix.Species_index['H2']].W # [kg/mol]
O2_mm = Mix.Species[Mix.Species_index['O2']].W # [kg/mol]
N2_mm = Mix.Species[Mix.Species_index['N2']].W # [kg/mol]

# Volume fraction of N2 ad O2 in Air
N2_vf = 0.79
O2_vf = 0.21

# Air molar mass
Air_mm = N2_mm*N2_vf + O2_mm*O2_vf

# Stoichiometric ratio (H2 mass / O2 mass)
alpha_s = 1/8

# Mass fraction of H2
Yi_H2 = phi*alpha_s / (phi*alpha_s + 1 + (N2_vf/O2_vf)*(N2_mm/O2_mm))

# Mass fraction of Air (N2 + O2)
Yi_Air = 1 - Yi_H2

# Mixture molar mass (H2 + Air)
Mix_mm = 1.0/(Yi_H2/H2_mm + Yi_Air/Air_mm)

# Molar fraction of H2
Xi_H2 = Yi_H2*Mix_mm/H2_mm

# Molar fraction of Air
Xi_Air = Yi_Air*Mix_mm/Air_mm

# Molar fraction of O2
Xi_O2 = 0.21*Xi_Air

# Molar fraction of N2
Xi_N2 = 0.79*Xi_Air

# Array of molar fractions
Xi = np.zeros(Mix.nSpec) + 0.0
Xi[Mix.Species_index['H2']] = Xi_H2
Xi[Mix.Species_index['O2']] = Xi_O2
Xi[Mix.Species_index['N2']] = Xi_N2

########################################
# Init condition for the test case
########################################
W = Mix.GetMolarWeightFromXi(Xi)
Yi = Mix.GetMassFractions(W, Xi)
rho = Mix.GetRho(P, T, W)
e0 = Mix.GetInternalEnergy(T, Yi)
rhoi = rho*Yi

dt = 1.0e-8
time = []
O2t = []
H2Ot = []
H2t = []
temperature = []

for n in range(4000000):

    rhoi, P, T, Yi = advance(rhoi, P, T, Yi, Mix, dt)
    time.append(dt*(n+1))
    print(time[n], P, T, rho)
    Xi = Mix.GetMolarFractions(Yi)
    H2t.append(Xi[Mix.Species_index['H2']])
    O2t.append(Xi[Mix.Species_index['O2']])
    H2Ot.append(Xi[Mix.Species_index['H2O']])
    temperature.append(T)
    if time[n] >= 1.0e-4: break

with open('mydata.csv', 'w') as fp:
    fp.write('%s,%s,%s,%s,%s\n' % ('t','H2','O2','H2O','T'))
    for i in range(0,len(H2t),100):
        fp.write('%16.8e,%16.8e,%16.8e,%16.8e,%16.8e\n' % (time[i], H2t[i], O2t[i], H2Ot[i], temperature[i]))

###################################################################################
# Plot
###################################################################################
import matplotlib.pyplot as plt
import pandas
outfreq = 200
cantera = pandas.read_csv('cantera.csv')
fig, ax = plt.subplots()
ax.set_xlim([0, 1.0e-4])

lns1 = ax.plot(time[::outfreq],  H2t[::outfreq], 'ob', fillstyle = 'none', label = 'H2')
lns2 = ax.plot(time[::outfreq],  O2t[::outfreq], 'or', fillstyle = 'none', label = 'O2')
lns3 = ax.plot(time[::outfreq], H2Ot[::outfreq], 'og', fillstyle = 'none', label = 'H2O')

ax.plot(cantera["t"]*1.0e-3,  cantera["H2"], '-b')
ax.plot(cantera["t"]*1.0e-3,  cantera["O2"], '-r')
ax.plot(cantera["t"]*1.0e-3, cantera["H2O"], '-g')

ax2 = ax.twinx()
lns4 = ax2.plot(time[::outfreq], temperature[::outfreq], 'ok', label = 'T')
ax2.plot(cantera["t"]*1.0e-3, cantera["T"], '-k')

lns = lns1 + lns2 + lns3 + lns4
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0)
ax.set_xlabel(r'$t (s)$')
ax.set_ylabel(r'Mole fractions')
ax2.set_ylabel(r'T (K)')
plt.show()
plt.close()
