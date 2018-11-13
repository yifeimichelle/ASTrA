#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import sys


eps0 = -180.95166
count = 0

if len(sys.argv)==3:
    skip=int(sys.argv[1])
    nlayers=int(sys.argv[2])
elif len(sys.argv)==1:
    print("No arguments given!\ncompute-voltage-drop.py -h\n     prints help")
    exit(0)
elif len(sys.argv)==2:
    if sys.argv[1]=="-h":
        print("script to conpute poisson potential from a trajectory\n")
        print("running:")
        print("compute-voltage-drop.py SKIP NLAYERS")
        print("  you must have the input file charges.inp in the execution folder\n")
        print("compute-voltage-drop.py --example\n  prints example input charges.inp.example\n  you must rename this to charges.inp to run\n")
        print("compute-voltage-drop.py --generate-input LX LY LZ NUMION NUMACN NUMELEC ELECCHARGE NUMCONFIGS\n  generate input file from args and num_elec_cap.in\n  CAUTION: will overwrite charges.inp\n")
        exit(0)
    elif sys.argv[1]=="--example":
        print("printing example input charges.inp.example\n  you must rename this to charges.inp to run")
        with open("charges.inp.example","w") as file:
            file.write("\
traj.xyz\n\
800     # nconfigs\n\
35.4166 30.6717 196.6078350808              # Lx, Ly, Lz\n\
6 3     # ndifftypes, maxnumatoms\n\
144 3   # type 1: nions, numatomstype\n\
0.4374\n\
0.1578\n\
0.1848\n\
144 1   # type 2: nions, numatomstype\n\
-0.78\n\
1344 3  # type 3: nions, numatomstype\n\
0.2690\n\
0.1290\n\
-0.398\n\
2232 1  # type 4: nions, numatomstype\n\
-0.01\n\
2232 1  # type 5: nions, numatomstype\n\
0.01\n\
842 1   # type 6: nions, numatomstype\n\
0.0\n\
5       # dtime")
        exit(0)
elif len(sys.argv)==10:
    if sys.argv[1]=="--generate-input":
        LX = sys.argv[2]
        LY = sys.argv[3]
        LZ = sys.argv[4]
        NUMION = sys.argv[5]
        NUMACN = sys.argv[6]
        NUMELEC = sys.argv[7]
        ELECCHARGE = sys.argv[8]
        NUMCONFIGS = sys.argv[9]
        with open("num_elec_cap.in", "r") as numelecap:
            line=numelecap.readline()
            NUMGRAPHENE=line.split()[1]
        with open("charges.inp","w") as file:
            string_to_write = '\
traj.xyz\n\
{}     # nconfigs\n\
{} {} {}              # Lx, Ly, Lz\n\
6 3     # ndifftypes, maxnumatoms\n\
{} 3   # type 1: nions, numatomstype\n\
0.4374\n\
0.1578\n\
0.1848\n\
{} 1   # type 2: nions, numatomstype\n\
-0.78\n\
{} 3  # type 3: nions, numatomstype\n\
0.2690\n\
0.1290\n\
-0.398\n\
{} 1  # type 4: nions, numatomstype\n\
-{}\n\
{} 1  # type 5: nions, numatomstype\n\
{}\n\
{} 1   # type 6: nions, numatomstype\n\
0.0\n\
5       # dtime'.format(NUMCONFIGS, LX, LY, LZ, NUMION, NUMION, NUMACN, NUMELEC, ELECCHARGE, NUMELEC, ELECCHARGE, NUMGRAPHENE)
            file.write(string_to_write)
        exit(0)
elif len(sys.argv)==14:
    if sys.argv[1]=="--generate-input-hybrid":
        LX = sys.argv[2]
        LY = sys.argv[3]
        LZ = sys.argv[4]
        NUMION = sys.argv[5]
        NUMACN = sys.argv[6]
        NUMOUTERGRAPHELEC = sys.argv[7]
        OUTERGRAPHELECCHARGE = sys.argv[8]
        NUMGRAPHELEC = sys.argv[9]
        GRAPHELECCHARGE = sys.argv[10]
        NUMELEC = sys.argv[11]
        ELECCHARGE = sys.argv[12]
        NUMCONFIGS = sys.argv[13]
        with open("num_hybrelec_cap.in", "r") as numelecap:
            line=numelecap.readline()
            s=line.split()
            NUMGRAPHENE=s[2]
        with open("charges.inp","w") as file:
            string_to_write = '\
traj.xyz\n\
{}     # nconfigs\n\
{} {} {}              # Lx, Ly, Lz\n\
6 3     # ndifftypes, maxnumatoms\n\
{} 3   # type 1: nions, numatomstype\n\
0.4374\n\
0.1578\n\
0.1848\n\
{} 1   # type 2: nions, numatomstype\n\
-0.78\n\
{} 3  # type 3: nions, numatomstype\n\
0.2690\n\
0.1290\n\
-0.398\n\
{} 1  # type 4: nions, numatomstype\n\
{}\n\
{} 1  # type 5: nions, numatomstype\n\
{}\n\
{} 1  # type 6: nions, numatomstype\n\
{}\n\
{} 1   # type 7: nions, numatomstype\n\
0.0\n\
5       # dtime'.format(NUMCONFIGS, LX, LY, LZ, NUMION, NUMION, NUMACN, NUMOUTERGRAPHELEC, OUTERGRAPHELECCHARGE, NUMGRAPHELEC, GRAPHELECCHARGE, NUMELEC, ELECCHARGE, NUMGRAPHENE)
            file.write(string_to_write)
        exit(0)

inpname = "charges.inp"

# READ INPUTS
inp = open(inpname, 'r')
trajname = inp.readline().split()[0]
[nconfigs] = [int(a) for a in inp.readline().split()[:1]]
Lx, Ly, Lz = [float(a) for a in inp.readline().split()[:3]]

zmin = np.zeros(nlayers,dtype=float)
zmax = np.zeros(nlayers,dtype=float)
chg = np.zeros(nlayers,dtype=float)
psi = np.zeros(nlayers,dtype=float)

print("num bins",nlayers)
print("skipping",skip," out of ",nconfigs)

dz=1.0*Lz/nlayers
Vbin=Lx*Ly*dz
print("dz",dz,", Vbin",Vbin)
for i in range(nlayers):
    zmin[i] = 1.0*(i)*Lz/(1.0*nlayers)
    zmax[i] = 1.0*(i+1)*Lz/(1.0*nlayers)

ndifftypes,maxnumatoms = [int(a) for a in inp.readline().split()[:2]]
numatomstype = np.zeros(ndifftypes,dtype=int)
nions = np.zeros(ndifftypes,dtype=int)
qtype = np.zeros([ndifftypes,maxnumatoms],dtype=float)
for i in range(ndifftypes):
    nions[i], numatomstype[i] = [float(a) for a in inp.readline().split()[:2]]
    for j in range(numatomstype[i]):
        [qtype[i,j]] = [float(a) for a in inp.readline().split()[:1]]
[dtime] = [float(a) for a in inp.readline().split()[:1]]
# FINISHED READING INPUTS

# READ TRAJECTORY
traj = open(trajname, 'r')
for step in range(skip):
    traj.readline()
    traj.readline()
    for itype in range(ndifftypes):
        for iion in range(nions[itype]):
            for iatom in range(numatomstype[itype]):
                traj.readline()
for step in range(nconfigs-skip):
    count += 1
    rho = np.zeros(nlayers,dtype=float)
    traj.readline()
    traj.readline()
    for itype in range(ndifftypes):
        for iion in range(nions[itype]):
            for iatom in range(numatomstype[itype]):
                z = float(traj.readline().split()[3])
                binlayer = int(z/dz)
                rho[binlayer] += qtype[itype,iatom]
                chg[binlayer] += qtype[itype,iatom]
    rho = 1.0*rho/Vbin
    sum1 = 0.0
    sum2 = 0.0
    qlo = 0.0
    qhi = 0.0
    # separated_rectangle:
    #inner = np.zeros(nlayers,dtype=float)
    #for ilayer in range(nlayers):
    #    sum1 += dz * rho[ilayer]
    #    inner[ilayer] = sum1
    #for ilayer in range(nlayers):
    #    sum2 += dz * inner[ilayer]
    #    psi[ilayer] = sum2 * eps0
    # combined_rectangle:
    for ilayer in range(nlayers):
        sum2 += dz * sum1
        sum1 += dz * rho[ilayer]
        psi[ilayer] += sum2 * eps0
    # combined_trapezoid:
    #for ilayer in range(nlayers):
    #    qhi = rho[ilayer]
    #    sum2 += 0.5 * dz * sum1
    #    sum1 += 0.5 * dz * ( qlo + qhi )
    #    sum2 += 0.5 * dz * sum1
    #    qlo = qhi
    #    psi[ilayer] = sum2 * eps0
# FINISHED READING TRAJECTORY

for ilayer in range(nlayers):
    chg[ilayer] /= (1.0*count*Vbin)
    psi[ilayer] /= (1.0*count)

np.insert(psi,0,0)
np.insert(zmax,0,0)
np.insert(chg,0,0)

np.savetxt("poisson_ave-conq" + str(qtype[3,0]) +"-" + str(skip) + "to" + str(nconfigs) + "-" + str(nlayers) + ".out", np.vstack([zmax,psi]).T)
np.savetxt("charges_ave-conq" + str(qtype[3,0]) +"-" + str(skip) + "to" + str(nconfigs) + "-" + str(nlayers) + ".out", np.vstack([zmax,chg]).T)
print("num steps used",count)
print("voltage",psi[nlayers-1]-psi[0])
