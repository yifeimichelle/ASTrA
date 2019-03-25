#!/usr/bin/env python3
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys
import os

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def rmse_norm(predictions, targets):
    average = targets.mean()
    return np.sqrt((((predictions - targets)/average) ** 2).mean())

COMPARE_TRAJ=False

if len(sys.argv) > 1:
    COMPARE_TRAJ=True
    plot_compare_folder = sys.argv[1]
    path=os.getcwd()
    split=path.split('/')
    compare_path='../../'+plot_compare_folder+'/'+split[-1]+'/'
    print(compare_path)

filename="netcharge"

data=np.loadtxt(filename)

# fit data to an exponential for entire curve
def func(x, a, b):
   return a*(1-np.exp(-x/b))
xdata=data[:,0]
ydata=data[:,4]
popt, pcov = curve_fit(func, xdata, ydata)
expfit=func(xdata, *popt)
Qmax = popt[0]
tau = popt[1]
timestep = xdata[1]-xdata[0]

# fit data to exponential for 1*tau part of curve
tauindex = int(np.floor(tau / timestep))
popt_1tau, pcov_1tau = curve_fit(func, xdata[:tauindex], ydata[:tauindex])
expfit_1tau = func(xdata, *popt_1tau)
Qmax_1tau = popt_1tau[0]
tau_1tau = popt_1tau[1]

# fit data to exponential for 1*tau part of curve
tauindex2 = 2*tauindex
popt_2tau, pcov_2tau = curve_fit(func, xdata[:tauindex2], ydata[:tauindex2])
expfit_2tau = func(xdata, *popt_2tau)
Qmax_2tau = popt_2tau[0]
tau_2tau = popt_2tau[1]

# compute RMSEs
rmse_fullfit_1tau=rmse_norm(expfit[:tauindex], ydata[:tauindex])
rmse_fullfit_2tau=rmse_norm(expfit[:tauindex2], ydata[:tauindex2])
rmse_1taufit_1tau=rmse_norm(expfit_1tau[:tauindex], ydata[:tauindex])
rmse_1taufit_2tau=rmse_norm(expfit_1tau[:tauindex2], ydata[:tauindex2])
rmse_2taufit_1tau=rmse_norm(expfit_2tau[:tauindex], ydata[:tauindex])
rmse_2taufit_2tau=rmse_norm(expfit_2tau[:tauindex2], ydata[:tauindex2])

# plot data and exponential fit(s)
plt.plot(xdata, ydata, 'b-', label='simulation')
if(COMPARE_TRAJ):
    comparedata=np.loadtxt(compare_path+filename)
    xcompare=comparedata[:,0]
    ycompare=comparedata[:,4]
    plt.plot(xcompare, ycompare, '-', c='cyan', label='indpt simulation')
plt.plot(xdata, expfit, '-', c='red', label=r'fit:  Q$_{max}$=%1.3e e, $\tau$=%.1f ps' % (Qmax, tau))
plt.plot(xdata, expfit_1tau, '-', c='gold', label=r'fit:  Q$_{max}$=%1.3e e, $\tau$=%.1f ps' % (Qmax_1tau, tau_1tau))
plt.plot(xdata, expfit_2tau, '-', c='orange', label=r'fit:  Q$_{max}$=%1.3e e, $\tau$=%.1f ps' % (Qmax_2tau, tau_2tau))

plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('time [ps]', fontsize=12)
plt.ylabel('charge per atom [e]', fontsize=12)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig("netcharge_expfit.png")

np.savetxt("charging_expfit_params", [[popt[1], popt[0], popt_1tau[1], popt_1tau[0], popt_2tau[1], popt_2tau[0]]])
np.savetxt("charging_expfit_rmse", [[rmse_fullfit_1tau, rmse_fullfit_2tau, rmse_1taufit_1tau, rmse_1taufit_2tau, rmse_2taufit_1tau, rmse_2taufit_2tau]])
np.savetxt("charging_expfit_pcovar", pcov)
np.savetxt("charging_expfit", np.transpose([xdata,expfit,expfit_1tau,expfit_2tau]))

# compute sum of squared residuals for first
