import Mixture
import Species

##########################################################################################
# Define Boivin Mixture species
##########################################################################################
# H2
Name = "H2"
W = 0.001008*2.0
iH2 = 0
Geom = 1
TSwitch1 = 1000.0
TSwitch2 = 6000.0
Tmin = 200.0
Tmax = 20000.0
cpH = [ 4.966884120e+08,-3.147547150e+05, 7.984121880e+01, \
       -8.414789210e-03, 4.753248350e-07,-1.371873490e-11, \
        1.605461760e-16, 2.488433520e+06,-6.695728110e+02]
cpM = [ 5.608128010e+05,-8.371504740e+02, 2.975364530e+00, \
        1.252249120e-03,-3.740716190e-07, 5.936625200e-11, \
       -3.606994100e-15, 5.339824410e+03,-2.202774770e+00]
cpL = [ 4.078323210e+04,-8.009186040e+02, 8.214702010e+00, \
       -1.269714460e-02, 1.753605080e-05,-1.202860270e-08, \
        3.368093490e-12, 2.682484660e+03,-3.043788840e+01]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL)
sigma = 2.92e-10
kbOveps = 0.026315789473684213
mu = 0.0
alpha = 7.900000000000001e-31
Z298 = 280.0
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
H2 = Species.Spec(Name, W, iH2, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff

# H
Name = "H"
W = 0.001008*1.0
iH = 1
Geom = 0
TSwitch1 = 1000.0
TSwitch2 = 6000.0
Tmin = 200.0
Tmax = 20000.0
cpH = [ 2.173757690E+08,-1.312035400E+05, 3.399174200E+01, \
       -3.813999680E-03, 2.432854840E-07,-7.694275540E-12, \
        9.644105630E-17, 1.067638090E+06,-2.742301050E+02]
cpM = [ 6.078774250E+01,-1.819354420E-01, 2.500211820E+00, \
       -1.226512860E-07, 3.732876330E-11,-5.687744560E-15, \
        3.410210200E-19, 2.547486400E+04,-4.481917770E-01]
cpL = [ 0.000000000E+00, 0.000000000E+00, 2.500000000E+00, \
        0.000000000E+00, 0.000000000E+00, 0.000000000E+00, \
        0.000000000E+00, 2.547370800E+04,-4.466828530E-01]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL)
sigma = 2.05e-10
kbOveps = 0.006896551724137932
mu = 0.0
alpha = 0.0
Z298 = 0.0
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
H = Species.Spec(Name, W, iH, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff

# O2
Name = "O2"
W = 0.015999*2.0
iO2 = 2
Geom = 1
TSwitch1 = 1000.0
TSwitch2 = 6000.0
Tmin = 200.0
Tmax = 20000.0
cpH = [ 4.975294300E+08,-2.866106870E+05, 6.690352250E+01, \
       -6.169959020E-03, 3.016396030E-07,-7.421416600E-12, \
        7.278175770E-17, 2.293554030E+06,-5.530621610E+02]
cpM = [-1.037939020E+06, 2.344830280E+03, 1.819732040E+00, \
        1.267847580E-03,-2.188067990E-07, 2.053719570E-11, \
       -8.193467050E-16,-1.689010930E+04, 1.738716510E+01]
cpL = [-3.425563420E+04, 4.847000970E+02, 1.119010960E+00, \
        4.293889240E-03,-6.836300520E-07,-2.023372700E-09, \
        1.039040020E-12,-3.391454870E+03, 1.849699470E+01]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL)
sigma = 3.4580000000000004e-10
kbOveps = 0.009310986964618248
mu = 0.0
alpha = 1.6e-30
Z298 = 3.8
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
O2 = Species.Spec(Name, W, iO2, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff

# OH
Name = "OH"
W = 0.001008*1.0+0.015999*1.0
iOH = 3
Geom = 1
TSwitch1 = 1000.0
TSwitch2 = 6000.0
Tmin = 200.0
Tmax = 20000.0
cpH = [ 2.847234190E+08,-1.859532610E+05, 5.008240900E+01, \
       -5.142374980E-03, 2.875536590E-07,-8.228817960E-12, \
        9.567229020E-17, 1.468393910E+06,-4.023555580E+02]
cpM = [ 1.017393380E+06,-2.509957280E+03, 5.116547860E+00, \
        1.305299930E-04,-8.284322260E-08, 2.006475940E-11, \
       -1.556993660E-15, 2.019640210E+04,-1.101282340E+01]
cpL = [-1.998858990E+03, 9.300136160E+01, 3.050854230E+00, \
        1.529529290E-03,-3.157891000E-06, 3.315446180E-09, \
       -1.138762680E-12, 2.991214230E+03, 4.674110790E+00]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL)
sigma = 2.7500000000000003e-10
kbOveps = 0.0125
mu = 0.0
alpha = 0.0
Z298 = 0.0
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
OH = Species.Spec(Name, W, iOH, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff

# O
Name = "O"
W = 0.015999*1.0
iO = 4
TSwitch1 = 1000.0
TSwitch2 = 6000.0
Tmin = 200.0
Tmax = 20000.0
cpH = [ 1.779004260E+08,-1.082328260E+05, 2.810778360E+01, \
       -2.975232260E-03, 1.854997530E-07,-5.796231540E-12, \
        7.191720160E-17, 8.890942630E+05,-2.181728150E+02]
cpM = [ 2.619020260E+05,-7.298722030E+02, 3.317177270E+00, \
       -4.281334360E-04, 1.036104590E-07,-9.438304330E-12, \
        2.725038300E-16, 3.392428060E+04,-6.679585350E-01]
cpL = [-7.953611300E+03, 1.607177790E+02, 1.966226440E+00, \
        1.013670310E-03,-1.110415420E-06, 6.517507500E-10, \
       -1.584779250E-13, 2.840362440E+04, 8.404241820E+00]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL)
Geom = 0
sigma = 2.7500000000000003e-10
kbOveps = 0.0125
mu = 0.0
alpha = 0.0
Z298 = 0.0
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
O = Species.Spec(Name, W, iO, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff

# H2O
Name = "H2O"
W = 0.001008*2.0+0.015999*1.0
iH2O = 5
TSwitch1 = 1000.0
TSwitch2 = 6000.0
Tmin = 200.0
Tmax = 6000.0
cpH = [ 1.034972100E+06,-2.412698560E+03, 4.646110780E+00, \
        2.291998310E-03,-6.836830480E-07, 9.426468930E-11, \
       -4.822380530E-15,-1.384286510E+04,-7.978148510E+00]
cpM = [ 1.034972100E+06,-2.412698560E+03, 4.646110780E+00, \
        2.291998310E-03,-6.836830480E-07, 9.426468930E-11, \
       -4.822380530E-15,-1.384286510E+04,-7.978148510E+00]
cpL = [-3.947960830E+04, 5.755731020E+02, 9.317826530E-01, \
        7.222712860E-03,-7.342557370E-06, 4.955043490E-09, \
       -1.336933250E-12,-3.303974310E+04, 1.724205780E+01]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL)
Geom = 3
sigma = 2.6050000000000003e-10
kbOveps = 0.0017470300489168416
mu = 6.150921915453923e-30
alpha = 0.0
Z298 = 4.0
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
H2O = Species.Spec(Name, W, iH2O, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff

# HO2
Name = "HO2"
W = 0.001008*1.0+0.015999*2.0
iHO2 = 6
TSwitch1 = 1000.0
TSwitch2 = 6000.0
Tmin = 200.0
Tmax = 6000.0
cpH = [-1.810669720E+06, 4.963192030E+03,-1.039498990E+00, \
        4.560148530E-03,-1.061859450E-06, 1.144567880E-10, \
       -4.763064160E-15,-3.200817190E+04, 4.066850920E+01]
cpM = [-1.810669720E+06, 4.963192030E+03,-1.039498990E+00, \
        4.560148530E-03,-1.061859450E-06, 1.144567880E-10, \
       -4.763064160E-15,-3.200817190E+04, 4.066850920E+01]
cpL = [-7.598882540E+04, 1.329383920E+03,-4.677388240E+00, \
        2.508308200E-02,-3.006551590E-05, 1.895600060E-08, \
       -4.828567390E-12,-5.873350960E+03, 5.193602140E+01]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL)
Geom = 3
sigma = 3.4580000000000004e-10
kbOveps = 0.009310986964618248
mu = 0.0
alpha = 0.0
Z298 = 1.0
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
HO2 = Species.Spec(Name, W, iHO2, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff

# H2O2
Name = "H2O2"
W = 0.001008*2.0+0.015999*2.0
iH2O2 = 7
TSwitch1 = 1000.0
TSwitch2 = 6000.0
Tmin = 200.0
Tmax = 6000.0
cpH = [ 1.489428030E+06,-5.170821780E+03, 1.128204970E+01, \
       -8.042397790E-05,-1.818383770E-08, 6.947265590E-12, \
       -4.827831900E-16, 1.418251040E+04,-4.650855660E+01]
cpM = [ 1.489428030E+06,-5.170821780E+03, 1.128204970E+01, \
       -8.042397790E-05,-1.818383770E-08, 6.947265590E-12, \
       -4.827831900E-16, 1.418251040E+04,-4.650855660E+01]
cpL = [-9.279533580E+04, 1.564748390E+03,-5.976460140E+00, \
        3.270744520E-02,-3.932193260E-05, 2.509255240E-08, \
       -6.465045290E-12,-2.494004730E+04, 5.877174180E+01]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL)
Geom = 3
sigma = 3.4580000000000004e-10
kbOveps = 0.009310986964618248
mu = 0.0
alpha = 0.0
Z298 =  3.8
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
H2O2 = Species.Spec(Name, W, iH2O2, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff

# N2
Name = "N2"
W = 0.014007*2.0
iN2 = 8
Geom = 1
TSwitch1 = 1000.0
TSwitch2 = 6000.0
Tmin = 200.0
Tmax = 20000.0
cpH = [ 8.310139160E+08,-6.420733540E+05, 2.020264640E+02, \
       -3.065092050E-02, 2.486903330E-06,-9.705954110E-11, \
        1.437538880E-15, 4.938707040E+06,-1.672099740E+03]
cpM = [ 5.877124060E+05,-2.239249070E+03, 6.066949220E+00, \
       -6.139685500E-04, 1.491806680E-07,-1.923105490E-11, \
        1.061954390E-15, 1.283210410E+04,-1.586640030E+01]
cpL = [ 2.210371500E+04,-3.818461820E+02, 6.082738360E+00, \
       -8.530914410E-03, 1.384646190E-05,-9.625793620E-09, \
        2.519705810E-12, 7.108460860E+02,-1.076003740E+01]
cpCoeff = Species.cpCoefficients(TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL)
sigma = 3.621e-10
kbOveps = 0.010253255408592227
mu = 0.0
alpha = 1.76e-30
Z298 = 4.0
DiffCoeff = Species.DiffCoefficients(sigma, kbOveps, mu, alpha, Z298)
N2 = Species.Spec(Name, W, iN2, Geom, cpCoeff, DiffCoeff)

# Clean for safety
del Name, W, Geom, TSwitch1, TSwitch2, Tmin, Tmax, cpH, cpM, cpL, cpCoeff, sigma, \
    kbOveps, mu, alpha, Z298, DiffCoeff

##########################################################################################
# Define Boivin Mixture reactions
##########################################################################################
import Reactions

# R0: H + O2 <=> O + OH
A = 3.520000000e+10
n = -7.000000000e-01
EovR = 8.589856793e+03
ArrCoeff = Reactions.ArrheniusCoefficients(A, n, EovR)
has_backward = True
Neducts = 2
Npducts = 2
educts = (Reactions.Reactant(iH, 1.0), Reactions.Reactant(iO2, 1.0))
pducts = (Reactions.Product(iO, 1.0), Reactions.Product(iOH, 1.0))
R0 = Reactions.Reaction('R0', ArrCoeff, has_backward, Neducts, Npducts, educts, pducts)

# Clean for safety
del A, n, EovR, ArrCoeff, has_backward, Neducts, Npducts, educts, pducts

# R1: H2 + O <=> H + OH
A = 5.060000000e-02
n = 2.670000000e+00
EovR = 3.165567894e+03
ArrCoeff = Reactions.ArrheniusCoefficients(A, n, EovR)
has_backward = True
Neducts = 2
Npducts = 2
educts = (Reactions.Reactant(iH2, 1.0), Reactions.Reactant(iO, 1.0))
pducts = (Reactions.Product(iH, 1.0), Reactions.Product(iOH, 1.0))
R1 = Reactions.Reaction('R1', ArrCoeff, has_backward, Neducts, Npducts, educts, pducts)

# Clean for safety
del A, n, EovR, ArrCoeff, has_backward, Neducts, Npducts, educts, pducts

# R2: H2 + OH <=> H + H20
A = 1.170000000e+03
n = 1.300000000e+00
EovR = 1.829343906e+03
ArrCoeff = Reactions.ArrheniusCoefficients(A, n, EovR)
has_backward = True
Neducts = 2
Npducts = 2
educts = (Reactions.Reactant(iH2, 1.0), Reactions.Reactant(iOH, 1.0))
pducts = (Reactions.Product(iH, 1.0), Reactions.Product(iH2O, 1.0))
R2 = Reactions.Reaction('R2', ArrCoeff, has_backward, Neducts, Npducts, educts, pducts)

# Clean for safety
del A, n, EovR, ArrCoeff, has_backward, Neducts, Npducts, educts, pducts

# R3: H + O2 (+M) => HO2 (+M)
A = 5.750000000e+07
n = -1.400000000e+00
EovR = 0.000000000e+00
ArrCoeffL = Reactions.ArrheniusCoefficients(A, n, EovR)
A = 4.650000000e+06
n = 4.400000000e-01
EovR = 0.000000000e+00
ArrCoeffH = Reactions.ArrheniusCoefficients(A, n, EovR)
has_backward = False
Neducts = 2
Npducts = 1
Nthirdb = 3
educts = (Reactions.Reactant(iH, 1.0), Reactions.Reactant(iO2, 1.0))
pducts = (Reactions.Product(iHO2, 1.0), Reactions.Product(0.0, 1.0))
thirdb = (Reactions.ThirdBd(iH2, 2.5), Reactions.ThirdBd(iH2O, 16.0), Reactions.ThirdBd(iO2, 1.0))
Ftype = 'F_Troe2'
FOdata = Reactions.Troe(alpha = 0.5, T1 = 1.0000000000000002e+30, T2 = 0.0, T3 = 1.0e-30)
R3 = Reactions.FalloffReaction('R3', ArrCoeffL, ArrCoeffH, has_backward, Neducts, Npducts, \
                               Nthirdb, educts, pducts, thirdb, Ftype, FOdata)

# Clean for safety
del A, n, EovR, ArrCoeffL, ArrCoeffH, has_backward, Neducts, Npducts, Nthirdb, educts, \
        pducts, thirdb, Ftype, FOdata

# R4: H + HO2 => 2 OH
A = 7.080000000e+07
n = 0.000000000e+00
EovR = 1.484497624e+02
ArrCoeff = Reactions.ArrheniusCoefficients(A, n, EovR)
has_backward = False
Neducts = 2
Npducts = 1
educts = (Reactions.Reactant(iH, 1.0), Reactions.Reactant(iHO2, 1.0))
pducts = (Reactions.Product(iOH, 2.0), Reactions.Product(0, 1.0))
R4 = Reactions.Reaction('R4', ArrCoeff, has_backward, Neducts, Npducts, educts, pducts)

# Clean for safety
del A, n, EovR, ArrCoeff, has_backward, Neducts, Npducts, educts, pducts

# R5: H + HO2 <=> H2 + O2
A = 1.660000000e+07
n = 0.000000000e+00
EovR = 4.141496761e+02
ArrCoeff = Reactions.ArrheniusCoefficients(A, n, EovR)
has_backward = True
Neducts = 2
Npducts = 2
educts = (Reactions.Reactant(iH, 1.0), Reactions.Reactant(iHO2, 1.0))
pducts = (Reactions.Product(iH2, 1.0), Reactions.Product(iO2, 1.0))
R5 = Reactions.Reaction('R5', ArrCoeff, has_backward, Neducts, Npducts, educts, pducts)

# Clean for safety
del A, n, EovR, ArrCoeff, has_backward, Neducts, Npducts, educts, pducts

# R6: HO2 + OH => H2O + O2
A = 2.890000000e+07
n = 0.000000000e+00
EovR = -2.501655267e+02
ArrCoeff = Reactions.ArrheniusCoefficients(A, n, EovR)
has_backward = False
Neducts = 2
Npducts = 2
educts = (Reactions.Reactant(iHO2, 1.0), Reactions.Reactant(iOH, 1.0))
pducts = (Reactions.Product(iH2O, 1.0), Reactions.Product(iO2, 1.0))
R6 = Reactions.Reaction('R6', ArrCoeff, has_backward, Neducts, Npducts, educts, pducts)

# Clean for safety
del A, n, EovR, ArrCoeff, has_backward, Neducts, Npducts, educts, pducts

# R7: H + OH + M <=> H2O + M
A = 4.000000000e+10
n = -2.000000000e+00
EovR = 0.000000000e+00
ArrCoeff = Reactions.ArrheniusCoefficients(A, n, EovR)
has_backward = True
Neducts = 2
Npducts = 1
Nthirdbd = 2
educts = (Reactions.Reactant(iH, 1.0), Reactions.Reactant(iOH, 1.0))
pducts = (Reactions.Product(iH2O, 1.0), Reactions.Product(0.0, 1.0))
thirdb = (Reactions.ThirdBd(iH2, 2.5), Reactions.ThirdBd(iH2O, 12.0))
R7 = Reactions.ThirdbodyReaction('R7', ArrCoeff, has_backward, Neducts, Npducts, Nthirdbd, educts, pducts, thirdb)

# Clean for safety
del A, n, EovR, ArrCoeff, has_backward, Neducts, Npducts, Nthirdbd, educts, pducts, thirdb

# R8: 2 H + M <=> H2 + M
A = 1.300000000e+06
n = -1.000000000e+00
EovR = 0.000000000e+00
ArrCoeff = Reactions.ArrheniusCoefficients(A, n, EovR)
has_backward = True
Neducts = 1
Npducts = 1
Nthirdbd = 2
educts = (Reactions.Reactant(iH, 2.0), Reactions.Reactant(0, 1.0))
pducts = (Reactions.Product(iH2, 1.0), Reactions.Product(0, 1.0))
thirdb = (Reactions.ThirdBd(iH2, 2.5), Reactions.ThirdBd(iH2O, 12.0))
R8 = Reactions.ThirdbodyReaction('R8', ArrCoeff, has_backward, Neducts, Npducts, Nthirdbd, educts, pducts, thirdb)

# Clean for safety
del A, n, EovR, ArrCoeff, has_backward, Neducts, Npducts, Nthirdbd, educts, pducts, thirdb

# R9: 2 HO2 => H2O2 + O2
A = 3.020000000e+06
n = 0.000000000e+00
EovR = 6.975780139e+02
ArrCoeff = Reactions.ArrheniusCoefficients(A, n, EovR)
has_backward = False
Neducts = 1
Npducts = 2
educts = (Reactions.Reactant(iHO2, 2.0), Reactions.Reactant(0, 1.0))
pducts = (Reactions.Product(iH2O2, 1.0), Reactions.Product(iO2, 1.0))
R9 = Reactions.Reaction('R9', ArrCoeff, has_backward, Neducts, Npducts, educts, pducts)

# Clean for safety
del A, n, EovR, ArrCoeff, has_backward, Neducts, Npducts, educts, pducts

# R10: H2 + HO2 => H + H2O2
A = 1.624390000e+05
n = 6.069610000e-01
EovR = 1.204355310e+04
ArrCoeff = Reactions.ArrheniusCoefficients(A, n, EovR)
has_backward = False
Neducts = 2
Npducts = 2
educts = (Reactions.Reactant(iH2, 1.0), Reactions.Reactant(iHO2, 1.0))
pducts = (Reactions.Product(iH, 1.0), Reactions.Product(iH2O2, 1.0))
R10 = Reactions.Reaction('R10', ArrCoeff, has_backward, Neducts, Npducts, educts, pducts)

# Clean for safety
del A, n, EovR, ArrCoeff, has_backward, Neducts, Npducts, educts, pducts

# R11: H2O2 (+M) => 2 OH (+M)
A = 8.153420000e+17
n = -1.918360000e+00
EovR = 2.497025647e+04
ArrCoeffL = Reactions.ArrheniusCoefficients(A, n, EovR)
A = 2.623280000e+19
n = -1.388360000e+00
EovR = 2.582673612e+04
ArrCoeffH = Reactions.ArrheniusCoefficients(A, n, EovR)
has_backward = False
Neducts = 1
Npducts = 1
Nthirdb = 2
educts = (Reactions.Reactant(iH2O2, 1.0), Reactions.Reactant(0, 1.0))
pducts = (Reactions.Product(iOH, 2.0), Reactions.Product(0, 2.0))
thirdb = (Reactions.ThirdBd(iH2, 2.0), Reactions.ThirdBd(iH2O, 6.0))
Ftype = 'F_Troe3'
FOdata = Reactions.Troe(alpha = 0.735, T1 = 1756, T2 = 5182.0, T3 = 94.0)
R11 = Reactions.FalloffReaction('R11', ArrCoeffL, ArrCoeffH, has_backward, Neducts, Npducts, \
                                Nthirdb, educts, pducts, thirdb, Ftype, FOdata)

# Clean for safety
del A, n, EovR, ArrCoeffL, ArrCoeffH, has_backward, Neducts, Npducts, Nthirdb, educts, \
        pducts, thirdb, Ftype, FOdata

# List of standard reactions
reacts = (R0, R1, R2, R4, R5, R6, R9, R10)
# List of third body reactions
tbreacts = (R7, R8)
# List of fall off reactions
foreacts = (R3, R11)

##########################################################################################
# Function to initialize a Boivin mixture
##########################################################################################
Species_index = {'H2': iH2, 'H': iH, 'O2': iO2, 'OH': iOH, 'O': iO, \
                 'H2O': iH2O,'HO2': iHO2, 'H2O2': iH2O2, 'N2': iN2}
def CreateMixture():
    return Mixture.Mixture(9, (H2, H, O2, OH, O, H2O, HO2, H2O2, N2), Species_index, \
                           8, reacts, 2, tbreacts, 2, foreacts)
