#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys

# Usage:
# compute_capacitance.py [voltage/V] [electrodedist/A] [nlayers] [chgfile] [skip]
# compute_capacitance_from_chgdensity $voltage $Lz 4000 charges.out 20

voltage=float(sys.argv[1])
Lz=float(sys.argv[2])
nlayers=int(sys.argv[3])
dz = (1.0 * Lz)/nlayers
Vbin = dz * 43.82 * 37.959
eps0 = -180.95

chgfile=sys.argv[4]
if len(sys.argv) > 4:
    skip = int(sys.argv[5])
else:
    skip = 0

print("num bins",nlayers)
print("skipping",skip)
print("dz",dz,", Vbin",Vbin)

with open(chgfile) as f:
    f.readline()
    f.readline()
    f.readline()
    line = f.readline()
    s = line.split()
    nlayers = int(s[1])
    psi = np.zeros(nlayers,dtype=float)
    chg = np.zeros(nlayers,dtype=float)

    count = 0
    for s in range(skip):
        for i in range(nlayers):
            line = f.readline()
        line = f.readline()
    while line != "":
        sum1 = 0
        sum2 = 0
        qlo = 0
        count += 1
        for i in range(nlayers):
            line = f.readline()
            s = line.split()
            rho = float(s[3])/Vbin
            sum2 += dz * sum1
            sum1 += dz * rho
            psi[i] += sum2 * eps0
            chg[i] += rho
        line = f.readline()

    psi = psi/(1.0*count)
    chg = chg/(1.0*count)

zmax = np.zeros(nlayers,dtype=float)
for i in range(nlayers):
    zmax[i] = 1.0*(i+1)*Lz/(1.0*nlayers)
np.insert(psi,0,0)
np.insert(zmax,0,0)
np.insert(chg,0,0)

np.savetxt("poisson_ave-" + str(voltage) +"-" + str(nlayers) + ".out", np.vstack([zmax,psi]).T)
np.savetxt("charges_ave-" + str(voltage) +"-" + str(nlayers) + ".out", np.vstack([zmax,chg]).T)
print("num steps used",count)
print("voltage",psi[nlayers-1]-psi[0])

fig, ax1 = plt.subplots()
plt.grid()
ax1.set_xlabel(r"z ($\AA$)")
ax1.plot(zmax, chg, 'g:')
ax1.set_ylabel(r"$\rho_{q}$ (e/$A^3$)", color='g')
ax1.tick_params('y', colors='g')
ax1.set_xlim(0,Lz)

ax2 = ax1.twinx()
ax2.plot(zmax, psi, 'b')
ax2.set_xlim(0,Lz)
ax2.set_ylabel(r"$\Psi$ (V)", color='b')
ax2.tick_params('y', colors='b')

fig.tight_layout()
plt.savefig("poisson_charge-" + str(voltage) + "-" + str(nlayers) + ".png")

print(zmax)
