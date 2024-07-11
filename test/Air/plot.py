import pandas as pd
import matplotlib.pyplot as plt

Ninni = pd.read_csv('Ninni.csv')
PCR = pd.read_csv('PCR.csv')

fig, ax = plt.subplots(ncols = 2, nrows = 1, figsize = (15,10))

ax[0].set_xlabel(r'$t [s]$')
ax[0].set_ylabel(r'$T [K]$')
ax[0].set_xlim([1.0e-9, 3.0e-3])
ax[0].set_xscale('log')
ax[0].plot(  PCR["time"],   PCR["T"],  '-b', label = 'PCR')
ax[0].plot(Ninni["time"], Ninni["T"], '--b', label = 'Ninni')
ax[0].legend()

ax[1].set_xlim([1.0e-9, 3.0e-3])
ax[1].set_ylim([1.0e-7, 1.0])
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].plot(PCR["time"], PCR["N2"], '-b', label = 'N2')
ax[1].plot(PCR["time"], PCR["O2"], '-r', label = 'O2')
ax[1].plot(PCR["time"], PCR["O"] , '-g', label = 'O')
ax[1].plot(PCR["time"], PCR["N"] , '-y', label = 'N')
ax[1].plot(PCR["time"], PCR["NO"], '-k', label = 'NO')

ax[1].plot(Ninni["time"], Ninni["N2"], '--b')
ax[1].plot(Ninni["time"], Ninni["O2"], '--r')
ax[1].plot(Ninni["time"], Ninni["O"] , '--g')
ax[1].plot(Ninni["time"], Ninni["N"] , '--y')
ax[1].plot(Ninni["time"], Ninni["NO"], '--k')

ax[1].legend(loc=0)
ax[1].set_xlabel(r'$t (s)$')
ax[1].set_ylabel(r'Mass fractions')
plt.savefig('plot.png')
plt.show()
plt.close()
