#!/usr/bin/env python3
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

if len(sys.argv)<2:
    print("average-potential-conq.py -h        prints help dialog\n")
    exit(0)
elif len(sys.argv)==2:
    if sys.argv[1]=="-h":
        print("usage:\n\n      average-potential-conq.py TOTSTEPS SKIP\n")
        exit(0)
elif len(sys.argv)==3:
    steps=int(sys.argv[1])
    skip=int(sys.argv[2])

fp=open('espot','r')
fo=open('avgpot','w')
fa=open('avgpotdiff','w')
line=fp.readline()
line=fp.readline()
line=fp.readline()
line=fp.readline()
num_atoms=int(line)
atoms_per_ele=num_atoms/2

fp.seek(0, 0) #seek to beginning of file

line=fp.readline()
counter=0
count_anode=0
count_cathode=0
anode_sum=0
cathode_sum=0
while line != '':
    if (counter<steps):
        anode=0
        cathode=0
        counter+=1
        timestep=fp.readline().split()[0]
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
            poti=float(s[1])
            if i>=num_atoms/2:
                anode += poti
            else:
                cathode += poti
        if (counter>skip):
            anode_sum += anode
            cathode_sum += cathode
        fo.write(str(counter)+" "+timestep+" "+str(cathode/atoms_per_ele)+" "+str(anode/atoms_per_ele)+"\n")
        line=fp.readline()
    else:
        break

norm=atoms_per_ele*(counter-skip)
print(cathode_sum/norm,anode_sum/norm)
fa.write(str(cathode_sum/norm)+" "+str(anode_sum/norm)+" "+str(anode_sum/norm-cathode_sum/norm)+" "+str(steps)+" "+str(skip))

fo.close()
fa.close()
fp.close()
