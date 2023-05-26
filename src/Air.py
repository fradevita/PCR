import Mixture
import Species

##########################################################################################
# Define Air Mixture species
##########################################################################################
# N2
Name = "N2"
W = 2*14.0067e-3
iN2 = 0
Geom = 1
TSwitch1 = 1000.0007
TSwitch2 = 6000.0007
Tmin = 200.0
TMax = 20000.0007
cpH = [ 8.310139160e+08,-6.420733540e+05, 2.020264635e+02, \
       -3.065092046e-02, 2.486903333e-06,-9.705954110e-11, \
        1.437538881e-15, 4.938707040e+06,-1.672099740e+03]
cpM = [ 5.877124060e+05,-2.239249073e+03, 6.066949220e+00, \
       -6.139685500e-04, 1.491806679e-07,-1.923105485e-11, \
        1.061954386e-15, 1.283210415e+04,-1.586640027e+01]
cpL = [ 2.210371497e+04,-3.818461820e+02, 6.082738360e+00, \
       -8.530914410e-03, 1.384646189e-05,-9.625793620e-09, \
        2.519705809e-12, 7.108460860e+02,-1.076003744e+01]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, TMax, cpH, cpM, cpL)
sigma = 3.621*Species.ATom
kbOveps = 1.0/97.53
mu = 0.000*Species.DToCm
alpha   = 1.760*Species.ATom
Z298    = 4.000
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
N2 = Species.Spec(Name, W, iN2, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff

# O2
Name = "O2"
W = 2*15.9994e-3
iO2 = 1
Geom = 1
TSwitch1 = 1000.0007
TSwitch2 = 6000.0007
Tmin = 200.0
TMax = 20000.0007
cpH = [ 4.975294300e+08,-2.866106874e+05, 6.690352250e+01, \
       -6.169959020e-03, 3.016396027e-07,-7.421416600e-12, \
        7.278175770e-17, 2.293554027e+06,-5.530621610e+02]
cpM = [-1.037939022e+06, 2.344830282e+03, 1.819732036e+00, \
        1.267847582e-03,-2.188067988e-07, 2.053719572e-11, \
       -8.193467050e-16,-1.689010929e+04, 1.738716506e+01]
cpL = [-3.425563420e+04, 4.847000970e+02, 1.119010961e+00, \
        4.293889240e-03,-6.836300520e-07,-2.023372700e-09, \
        1.039040018e-12,-3.391454870e+03, 1.849699470e+01]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, TMax, cpH, cpM, cpL)
sigma = 3.458*Species.ATom
kbOveps = 1.0/107.40
mu = 0.000*Species.DToCm
alpha   = 1.600*Species.ATom
Z298    = 3.800
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
O2 = Species.Spec(name, W, iO2, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff

# NO
Name = "NO"
W = 14.0067e-3+15.9994e-3
iNO = 2
Geom = 1
TSwitch1 = 1000.0007
TSwitch2 = 6000.0007
Tmin = 200.0
TMax = 20000.0007
cpH = [-9.575303540e+08, 5.912434480e+05,-1.384566826e+02, \
        1.694339403e-02,-1.007351096e-06, 2.912584076e-11, \
       -3.295109350e-16,-4.677501240e+06, 1.242081216e+03]
cpM = [ 2.239018716e+05,-1.289651623e+03, 5.433936030e+00, \
       -3.656034900e-04, 9.880966450e-08,-1.416076856e-11, \
        9.380184620e-16, 1.750317656e+04,-8.501669090e+00]
cpL = [-1.143916503e+04, 1.536467592e+02, 3.431468730e+00, \
       -2.668592368e-03, 8.481399120e-06,-7.685111050e-09, \
        2.386797655e-12, 9.098214410e+03, 6.728725490e+00]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, TMax, cpH, cpM, cpL)
sigma = 3.621*Species.ATom
kbOveps = 1.0/97.530
mu = 0.000*Species.DToCm
alpha   = 1.760*Species.ATom
Z298    = 4.000
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
NO = Species.Spec(name, W, iNO, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff

# N
Name = "N"
W = 14.0067e-3
iN = 3
Geom = 0
TSwitch1 = 1000.0007
TSwitch2 = 6000.0007
Tmin = 200.0
TMax = 20000.0007
cpH = [ 5.475181050e+08,-3.107574980e+05, 6.916782740e+01, \
       -6.847988130e-03, 3.827572400e-07,-1.098367709e-11, \
        1.277986024e-16, 2.550585618e+06,-5.848769753e+02]
cpM = [ 8.876501380e+04,-1.071231500e+02, 2.362188287e+00, \
        2.916720081e-04,-1.729515100e-07, 4.012657880e-11, \
       -2.677227571e-15, 5.697351330e+04, 4.865231506e+00]
cpL = [ 0.000000000e+00, 0.000000000e+00, 2.500000000e+00, \
        0.000000000e+00, 0.000000000e+00, 0.000000000e+00, \
        0.000000000e+00, 5.610463780e+04, 4.193905036e+00]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, TMax, cpH, cpM, cpL)
sigma = 3.298*Species.ATom
kbOveps = 1.0/71.400
mu = 0.000*Species.DToCm
alpha   = 0.000*Species.ATom
Z298    = 0.000
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
N = Species.Spec(name, W, iN, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff


# O
Name = "O"
W = 15.9994e-3
iO = 4
Geom = 0
TSwitch1 = 1000.0007
TSwitch2 = 6000.0007
Tmin = 200.0
TMax = 20000.0007
cpH = [ 1.779004264e+08,-1.082328257e+05, 2.810778365e+01, \
       -2.975232262e-03, 1.854997534e-07,-5.796231540e-12, \
        7.191720164e-17, 8.890942630e+05,-2.181728151e+02]
cpM = [ 2.619020262e+05,-7.298722030e+02, 3.317177270e+00, \
       -4.281334360e-04, 1.036104594e-07,-9.438304330e-12, \
        2.725038297e-16, 3.392428060e+04,-6.679585350e-01]
cpL = [-7.953611300e+03, 1.607177787e+02, 1.966226438e+00, \
        1.013670310e-03,-1.110415423e-06, 6.517507500e-10, \
       -1.584779251e-13, 2.840362437e+04, 8.404241820e+00]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, TMax, cpH, cpM, cpL)
sigma = 2.750*Species.ATom
kbOveps = 1.0/80.0
mu = 0.000*Species.DToCm
alpha   = 0.000*Species.ATom
Z298    = 0.000
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
O = Species.Spec(name, W, iO, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff

##########################################################################################
# Define Air Mixture reactions
##########################################################################################
import Reactions
# TO DO:...

##########################################################################################
# Function to initialize an Air mixture
##########################################################################################
Species_index = {'N2': iH2, 'O2': iO2, 'NO': iNO, 'N': iN, 'O': iO}

def CreateAirMixture():
    return Mixture.Mixture(5, (N2, O2, NO, N, O)m Species_index)
