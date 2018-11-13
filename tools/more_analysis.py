#! /usr/bin/env python
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys

# plot densities of solvent across electrode (COM and mass density)
# and radial distribution functions between framework and solvent
# Usage:
# more_analysis.py [voltage/V] [electrodedist/A]

voltage=float(sys.argv[1])
electrodedist=int(sys.argv[2])
void_fraction=0.48112
dz = (1.0 * electrodedist + 32.2932)/100.0
Vbin = dz * 43.82 * 37.959

molar_masses = {'tea': 130.55, 'bf4': 86.80, 'acn': 41.05}

def solvent_num(x):
    return {
            90: {'tea' : 91, 'bf4': 91, 'acn': 1361},
            120: {'tea' : 121, 'bf4': 121, 'acn': 1816},
            150: {'tea' : 151, 'bf4': 151, 'acn': 2271},
            }[x]

def charge(x):
    return {
            0: 0,
            0.2: 0.184262,
            0.4: 0.368524,
            0.6: 0.552786,
            0.8: 0.737048,
            1: 0.92131,
            }[x]

solvent_list=['tea', 'bf4', 'acn']
for i in range(3):
    solventmol=solvent_list[i]
    num_solvent = solvent_num(electrodedist)[solventmol]

    com_dens_data = np.loadtxt(solventmol + "_comz_hist_last.out", skiprows=4)
    mass_dens_data = np.loadtxt("dens_" + solventmol + "_last.out", skiprows=4)

    com_dens_gcm3 = com_dens_data[:,3] * num_solvent * molar_masses[solventmol] * 10 / Vbin / 6.022

    plt.plot(com_dens_data[:,1], com_dens_gcm3, label='COM')
    plt.plot(com_dens_data[:,1], mass_dens_data[:,3], label='mass')
    plt.xlabel(r"z ($\AA$)")
    plt.ylabel(r"$\rho$ (g/cm$^3$)")
    plt.legend()
    plt.savefig("density_" + solventmol + ".png")
    plt.close('all')

rdf_data = np.loadtxt("rdf_running_last.out", skiprows=4)
rdf_labels = ["MOF_TEA_gr", "MOF_TEA_co", "MOF_BF4_gr", "MOF_BF4_co", "MOF_ACN_gr", "MOF_ACN_co", "MOF_ALL_gr", "MOF_ALL_co", "MOFC_ALL_gr", "MOFC_ALL_co", "TEA_BF4_gr", "TEA_BF4_co", "TEA_ACN_gr", "TEA_ACN_co", "BF4_ACN_gr", "BF4_ACN_co"]
for i in range(16):
    plt.plot(rdf_data[:,1],rdf_data[:,i+2],label=rdf_labels[i])
    plt.xlabel(r"r ($\AA$)")
    if (i % 2 == 0):
        plt.ylabel("g(r)")
        plt.ylim(0,3)
    else:
        plt.ylabel("coordination")
        plt.ylim(0,700)
    plt.legend()
    plt.savefig("rdf_" + rdf_labels[i] + ".png")
    plt.close('all')

try:
    rdf_ano_data = np.loadtxt("rdf_ano_running_last.out", skiprows=4)
except IOError:
    print 'cannot open rdf_ano_running_last.out'
else:
    rdf_labels = ["MOF_TEA_gr", "MOF_TEA_co", "MOF_BF4_gr", "MOF_BF4_co", "MOF_ACN_gr", "MOF_ACN_co", "MOF_ALL_gr", "MOF_ALL_co", "MOFC_ALL_gr", "MOFC_ALL_co", "TEA_BF4_gr", "TEA_BF4_co", "TEA_ACN_gr", "TEA_ACN_co", "BF4_ACN_gr", "BF4_ACN_co"]
    for i in range(16):
        plt.plot(rdf_ano_data[:,1],rdf_ano_data[:,i+2],label=rdf_labels[i])
        plt.xlabel(r"r ($\AA$)")
        if (i % 2 == 0):
            plt.ylabel("g(r)")
            plt.ylim(0,3)
        else:
            plt.ylabel("coordination")
            plt.ylim(0,700)
        plt.legend()
        plt.savefig("rdf_ano_" + rdf_labels[i] + ".png")
        plt.close('all')

try:
    rdf_cat_data = np.loadtxt("rdf_cat_running_last.out", skiprows=4)
except IOError:
    print 'cannot open rdf_cat_running_last.out'
else:
    rdf_labels = ["MOF_TEA_gr", "MOF_TEA_co", "MOF_BF4_gr", "MOF_BF4_co", "MOF_ACN_gr", "MOF_ACN_co", "MOF_ALL_gr", "MOF_ALL_co", "MOFC_ALL_gr", "MOFC_ALL_co", "TEA_BF4_gr", "TEA_BF4_co", "TEA_ACN_gr", "TEA_ACN_co", "BF4_ACN_gr", "BF4_ACN_co"]
    for i in range(16):
        plt.plot(rdf_cat_data[:,1],rdf_cat_data[:,i+2],label=rdf_labels[i])
        plt.xlabel(r"r ($\AA$)")
        if (i % 2 == 0):
            plt.ylabel("g(r)")
            plt.ylim(0,3)
        else:
            plt.ylabel("coordination")
            plt.ylim(0,700)
        plt.legend()
        plt.savefig("rdf_cat_" + rdf_labels[i] + ".png")
        plt.close('all')

try:
    rdf_ano_data = np.loadtxt("rdf_ano_running_last.out", skiprows=4)
    rdf_cat_data = np.loadtxt("rdf_cat_running_last.out", skiprows=4)
except IOError:
    print 'cannot open rdf_ano_running_last.out, cannot open rdf_cat_running_last.out'
else:
    rdf_labels = ["MOF_TEA_gr", "MOF_TEA_co", "MOF_BF4_gr", "MOF_BF4_co", "MOF_ACN_gr", "MOF_ACN_co", "MOF_ALL_gr", "MOF_ALL_co", "MOFC_ALL_gr", "MOFC_ALL_co", "TEA_BF4_gr", "TEA_BF4_co", "TEA_ACN_gr", "TEA_ACN_co", "BF4_ACN_gr", "BF4_ACN_co"]
    for i in range(16):
        plt.plot(rdf_ano_data[:,1],rdf_ano_data[:,i+2],label=rdf_labels[i]+"_ano")
        plt.plot(rdf_cat_data[:,1],rdf_cat_data[:,i+2],label=rdf_labels[i]+"_cat")
        plt.xlabel(r"r ($\AA$)")
        if (i % 2 == 0):
            plt.ylabel("g(r)")
            plt.ylim(0,3)
        else:
            plt.ylabel("coordination")
            plt.ylim(0,700)
        plt.legend()
        plt.savefig("rdf_both_" + rdf_labels[i] + ".png")
        plt.close('all')
