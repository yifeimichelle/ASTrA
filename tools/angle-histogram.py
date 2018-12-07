#!/usr/bin/env python3
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

if len(sys.argv)<2:
    print("angle-histogram.py -h        prints help dialog\n")
    exit(0)
elif len(sys.argv)==2:
    if sys.argv[1]=="-h":
        print("usage:\n\n      angle-histogram.py TOTSTEPS SKIP ANGLENAME NUMANGLES\n")
        exit(0)
elif len(sys.argv)==5:
    steps=int(sys.argv[1])
    skip=int(sys.argv[2])
    anglename=sys.argv[3]
    numangles=int(sys.argv[4])
dir_path = os.path.dirname(os.path.realpath(__file__))
path = os.getcwd()
s=path.split("_")
structure=s[-3]
#print(structure)

fp=open(anglename+'.dump','r')
line=fp.readline()
line=fp.readline()
line=fp.readline()
line=fp.readline()
num_atoms=int(line)

fp.seek(0, 0) #seek to beginning of file

line=fp.readline()
counter=0
count_anode=0
count_cathode=0
anglesum=0
anglevec=np.zeros(numangles*steps,dtype=float)
while (line != '') and (counter < steps):
    anode_sum=0
    cathode_sum=0
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
        theta=float(s[5])
        if (counter>skip):
            anglesum+=theta
            anglevec[counter]=theta
        counter+=1
    line=fp.readline()

anglesum /= counter
hist,bins = np.histogram(anglevec,bins=200,range=(160,180),density=True)
binwidth=bins[1]-bins[0]
np.savetxt(anglename+'_hist.txt',np.vstack([bins[:-1],hist]).T,header="width = "+str(binwidth))
lessthan179 = 0
for i in range(190):
    lessthan179 += binwidth*hist[i]
print(anglesum, lessthan179)
