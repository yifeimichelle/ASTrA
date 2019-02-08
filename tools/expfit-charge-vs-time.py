#!/usr/bin/env python3
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
data=np.loadtxt("netcharge")
def func(x, a, b):
   return a*(1-np.exp(-x/b))
xdata=data[:,0]
ydata=data[:,4]
popt, pcov = curve_fit(func, xdata, ydata)
expfit=func(xdata, *popt)

plt.plot(xdata, ydata, 'b-', label='data')
plt.plot(xdata, expfit, 'r-', label='fit:  a=%1.6e, b=%1.6e' % tuple(popt))

plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.savefig("netcharge_expfit.png")

np.savetxt("charging_expfit_params",popt)
np.savetxt("charging_expfit_pcovar",pcov)
np.savetxt("charging_expfit",np.transpose([xdata,expfit]))
