#!/usr/bin/env python3
import numpy as np
import sys

filename="density.out"

lims=np.loadtxt("layers.out",skiprows=1)
layer1=lims[0,0]
layer2=lims[1,0]

data = np.loadtxt(filename,dtype=float,skiprows=1)
z = data[:,0]
r = data[:,1]
lo = layer1+15.0
hi = layer2-15.0
avg = 0.0
count = 0
for i in range(len(z)):
    if ( z[i] > lo ) and ( z[i] < hi ):
        avg += r[i]
        count += 1

print(1.0*avg/count)
