#!/usr/bin/env python

import numpy as np
import sys
import math

if len(sys.argv)>1:
    numlayers=int(sys.argv[1])
elif len(sys.argv)==1:
    numlayers=2

layers=np.loadtxt("layers.out",skiprows=1,usecols=(0))
charges=np.loadtxt("elecharge.out",skiprows=1)

height=[0.0,0.0,0.0]
height[0]=layers[0]-3
height[1]=layers[1]-layers[0]
height[2]=layers[2]-layers[1]-3

catcharge=[]
anocharge=[]

for i in range(numlayers):
    catcharge.append(0.0)
    anocharge.append(0.0)

for i in range(len(charges[:,0])):
    z=charges[i,0]
    q=charges[i,1]
    if z<layers[0] and z>3:
        h=z-3
        lbin=math.floor(h/(height[0]/numlayers))
        catcharge[lbin]+=q
    if z>layers[1] and z<layers[2]-3:
        h=z-layers[1]
        lbin=math.floor(h/(height[2]/numlayers))
        anocharge[lbin]+=q



catchargemean = catcharge/np.mean(catcharge)
anochargemean = anocharge/np.mean(anocharge)
print(catchargemean)
print(anochargemean[::-1])

np.savetxt("charge_layers.out",np.transpose([range(numlayers),catcharge,anocharge[::-1],catchargemean,anochargemean[::-1]]),fmt=('%d','%.5e','%.5e','%10.5f','%.5f'))
