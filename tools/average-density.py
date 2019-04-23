#!/usr/bin/env python3
import numpy as np
import sys

filename="density.out"

lims=np.loadtxt("layers.out",skiprows=1)
layer1=lims[0,0]
layer2=lims[1,0]
layer3=lims[2,0]

data = np.loadtxt(filename,dtype=float,skiprows=1)
z = data[:,0]
r = data[:,1]
lo_ele = layer1-5.0
lo = layer1+15.0
hi = layer2-15.0
hi_ele = layer2+5.0
anoavg = 0.0
anocount = 0
catavg = 0.0
catcount = 0
bulkavg = 0.0
bulkcount = 0
for i in range(len(z)):
    if (z[i] < lo_ele) and (z[i] > 3.0):
        catavg += r[i]
        catcount += 1
    elif ( z[i] > lo ) and ( z[i] < hi ):
        bulkavg += r[i]
        bulkcount += 1
    elif (z[i] > hi_ele) and (z[i] < layer3-3.0):
        anoavg += r[i]
        anocount += 1

print(1.0*catavg/catcount, 1.0*anoavg/anocount, 1.0*bulkavg/bulkcount)
