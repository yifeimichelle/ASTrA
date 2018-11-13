#!/usr/bin/env python3
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math

if len(sys.argv)<2:
    print("compute-capacitance.py -h        prints help dialog\n")
    exit(0)
elif len(sys.argv)==2:
    if sys.argv[1]=="-h":
        print("usage:\n\n      compute-capacitance.py TOTSTEPS SKIP\n")
        exit(0)
elif len(sys.argv)==3:
    steps=int(sys.argv[1])
    skip=int(sys.argv[2])
    suffix=""
elif len(sys.argv)==4:
    steps=int(sys.argv[1])
    skip=int(sys.argv[2])
    suffix="_"+sys.argv[3]

fp=open('ele','r')
line=fp.readline()
line=fp.readline()
line=fp.readline()
line=fp.readline()
num_atoms=int(line)

print(str(num_atoms)+" atoms")

fp.seek(0, 0) #seek to beginning of file

line=fp.readline()

qbinsize=0.001
mincharge=-0.5
maxcharge=0.5
numqbins=1.0*(maxcharge-mincharge)/qbinsize+1

qanodehist=np.zeros(int(numqbins),dtype=float)
qcathodehist=np.zeros(int(numqbins),dtype=float)

stepcounter=0
count_anode=0
count_cathode=0
while line != '' and stepcounter<steps:
    line=fp.readline()
    line=fp.readline()
    line=fp.readline()
    line=fp.readline()
    line=fp.readline()
    line=fp.readline()
    line=fp.readline()
    line=fp.readline()
    for i in range(num_atoms):
        line=fp.readline()
        s=line.split()
        qi=float(s[1])
        # check whether to use charge
        if (stepcounter>skip):
            # bin charge
            qbin=math.floor((qi-mincharge)/qbinsize)
            if (1.0*i/num_atoms<0.5):
                qcathodehist[qbin] += 1
                count_cathode += 1
            else:
                qanodehist[qbin] += 1
                count_anode += 1
    line=fp.readline()
    stepcounter+=1

print(count_anode, count_cathode)
for i in range(int(numqbins)):
    qanodehist[i] /= (1.0*count_anode)
    qcathodehist[i] /= (1.0*count_cathode)

bins=np.linspace(mincharge,maxcharge,numqbins)
np.savetxt("chargehistogram"+suffix+".out",np.column_stack((bins, qanodehist, qcathodehist)))
