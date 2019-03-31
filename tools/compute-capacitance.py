#!/usr/bin/env python3
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys

if len(sys.argv)<2:
    print("compute-capacitance.py -h        prints help dialog\n")
    exit(0)
elif len(sys.argv)==2:
    if sys.argv[1]=="-h":
        print("usage:\n\n      compute-capacitance.py TOTSTEPS SKIP HALFVOLTAGE VOL ASA NUMATOMS\n")
        exit(0)
elif len(sys.argv)==6:
    steps=int(sys.argv[1])
    skip=int(sys.argv[2])
    halfvoltage=float(sys.argv[3])
    vol=float(sys.argv[4])
    asa=float(sys.argv[5])
    #numatoms=1
elif len(sys.argv)==7:
    steps=int(sys.argv[1])
    skip=int(sys.argv[2])
    halfvoltage=float(sys.argv[3])
    vol=float(sys.argv[4])
    asa=float(sys.argv[5])
    numatoms=float(sys.argv[6])
fullpath = os.path.realpath(os.getcwd())
structure = os.path.basename(fullpath)
username = fullpath.split("/")[4]

fp=open('ele','r')
line=fp.readline()
line=fp.readline()
line=fp.readline()
line=fp.readline()
num_total2elecatoms=int(line)

facGrav = 1.60218e-19 * 6.022e23 / (12.011 * 0.5*num_total2elecatoms)
facVol = 1.60218e-19 * 1.e24 / vol
facArea = 1.60218e-19 * 1.e16 * 1.e6 / asa

fp.seek(0, 0) #seek to beginning of file

line=fp.readline()
anoNetChg=np.empty(0,dtype=float)
catNetChg=np.empty(0,dtype=float)
capPerAtom=np.zeros(num_total2elecatoms,dtype=float)
catAtomChgs=np.zeros(int(num_total2elecatoms*steps/2),dtype=float)
anoAtomChgs=np.zeros(int(num_total2elecatoms*steps/2),dtype=float)
counter=0
count_anode=0
count_cathode=0
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
    for i in range(num_total2elecatoms):
        line=fp.readline()
        s=line.split()
        qi=float(s[1])
        if (counter>skip):
            capPerAtom[i] += qi/halfvoltage
            if i>=num_total2elecatoms/2:
                catAtomChgs[counter*int(num_total2elecatoms/2)+i-int(num_total2elecatoms/2)]=qi/halfvoltage
            else:
                anoAtomChgs[counter*int(num_total2elecatoms/2)+i]=qi/halfvoltage
        if (1.0*i/num_total2elecatoms<0.5):
            anode_sum += qi
            count_anode += 1
        else:
            cathode_sum += qi
            count_cathode += 1
    anoNetChg=np.append(anoNetChg,[anode_sum])
    catNetChg=np.append(catNetChg,[cathode_sum])
    counter+=1
    line=fp.readline()

time=np.zeros(counter,dtype=float)
for i in range(counter):
  time[i] = i*0.5 # assuming frames printed every 500 fs = 0.5 ps
capPerAtom /= counter
outputArray=np.vstack([time,anoNetChg,catNetChg,anoNetChg+catNetChg,(catNetChg-anoNetChg)/(2.0*numatoms)]).T
mean=np.mean(outputArray[skip:,1:3], axis=0)
meanchgmag=0.5*(abs(mean[0]) + abs(mean[1]))
np.savetxt('netcharge',outputArray,header="time[ps] catNetChg anoNetChg sum avgChgPerAtom . USER:"+username,comments="#")
np.savetxt('avgNetCharge',mean,header="anode cathode",comments="#")
#print(mean[1]*facGrav/(halfvoltage*2.0),"F/g")
#print(mean[1]*facVol/(halfvoltage*2.0),"F/cm3")
#print(mean[1]*facArea/(halfvoltage*2.0),"uF/cm2")

fout=open("../../all_capacitances","w+")
print(username,structure, vol, asa, meanchgmag*facGrav/(2.0*halfvoltage), meanchgmag*facVol/(2.0*halfvoltage), meanchgmag*facArea/(2.0*halfvoltage))

#capPerAtomHist,bin_edges=np.hist(capPerAtom,bins=50)
#catAtomChgsHist,bin_edges=np.hist(catAtomChgs,bins=50)
#capPerAtomHist /= (counter-skip)
#catAtomChgsHist /= (counter-skip)
hist,bins=np.histogram(capPerAtom[int(num_total2elecatoms/2):],bins=100,density=True)
plt.bar(bins[:-1],hist,width=(bins[1]-bins[0]))
plt.ylabel("probability density (normalized)")
plt.xlabel("per-atom capacitance")
plt.savefig("cat_capPerAtom.png")
plt.close("all")
np.savetxt('cat_capPerAtom.txt',np.vstack([bins[:-1],hist]).T,header="width = "+str(bins[1]-bins[0]))

hist,bins=np.histogram(capPerAtom[:int(num_total2elecatoms/2.0)],bins=100,density=True)
plt.bar(bins[:-1],hist,width=(bins[1]-bins[0]))
plt.ylabel("probability density (normalized)")
plt.xlabel("per-atom capacitance")
plt.savefig("ano_capPerAtom.png")
plt.close("all")
np.savetxt('ano_capPerAtom.txt',np.vstack([bins[:-1],hist]).T,header="width = "+str(bins[1]-bins[0]))

hist,bins=np.histogram(catAtomChgs[:int(counter*steps/2)],bins=100,density=True)
plt.bar(bins[:-1],hist,width=(bins[1]-bins[0]))
plt.ylabel("probability density (normalized)")
plt.xlabel("per-atom capacitance")
plt.savefig("catAtomChgs.png")
plt.close("all")
np.savetxt('catAtomChgs.txt',np.vstack([bins[:-1],hist]).T,header="width = "+str(bins[1]-bins[0]))

hist,bins=np.histogram(anoAtomChgs[:int(counter*steps/2)],bins=100,density=True)
plt.bar(bins[:-1],hist,width=(bins[1]-bins[0]))
plt.ylabel("probability density (normalized)")
plt.xlabel("per-atom capacitance")
plt.savefig("anoAtomChgs.png")
plt.close("all")
np.savetxt('anoAtomChgs.txt',np.vstack([bins[:-1],hist]).T,header="width = "+str(bins[1]-bins[0]))
