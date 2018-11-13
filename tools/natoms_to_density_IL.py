#!/usr/bin/env python3
import numpy as np
import sys

Lx=float(sys.argv[1])
Ly=float(sys.argv[2])
Lz=float(sys.argv[3])
zmin1=float(sys.argv[4])
zmax1=float(sys.argv[5])
zmin2=float(sys.argv[6])
zmax2=float(sys.argv[7])
zmin3=float(sys.argv[8])
zmax3=float(sys.argv[9])
inputfile=sys.argv[10]

factor=1.660577881
mass=np.array([139.23, 86.81])

volume=np.zeros(3)
volume[0]=Lx*Ly*(zmax1-zmin1)
volume[1]=Lx*Ly*(zmax2-zmin2)
volume[2]=Lx*Ly*(zmax3-zmin3)

data=np.loadtxt(inputfile,usecols=[1,2,3])
density=np.zeros([3,3])
for i in range(2):
    for j in range(3):
        density[i,j]=data[i,j]*factor*mass[i]/volume[j]

for j in range(3):
    density[2,j]=density[0,j]+density[1,j]+density[1,j]

np.savetxt("avg_densities_"+inputfile[5:],density)

print(data)
print(volume)
print(density[2,:])
