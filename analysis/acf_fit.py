# -*- coding: utf-8 -*-

#             Fits ACF data with two exponent function
#
#
#    Fitting function:   
#          
#    ACF_fit = (1 - S^2)(w1*exp(-t/tau1) + w2*exp(-t/tau2)) + S^2,
#    
#    where S is lipid order parameter calculated earlier
#
#    Resulting tau is:   tau = w1\*tau1 + w2\*tau2
#
#
#



from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#plt.style.use('seaborn-white')
import numpy as np
import pandas as pd

# n is the first frame to start the fitting
# n2 is the last one
n = 50
n2 = 1000

# files with ACF
m = np.loadtxt('lipodisk_rotacf_edge.xvg', comments=['#','@','&'])
m1 = np.loadtxt('lipodisk_rotacf_center.xvg', comments=['#','@','&'])
m2 = np.loadtxt('lipodisk_rotacf_all.xvg', comments=['#','@','&'])


# files with Lipid Order Parameters
# Here ACFs are fitted for C12 fatty acid carbon

try:
    data_all = pd.read_csv('LOP_all.csv')
    S_all = data_all.mean(axis=0).loc['C12']
    S2all = S_all*S_all

    data_inner = pd.read_csv('LOP_inner.csv')
    S_inner = data_inner.mean(axis=0).loc['C12']
    S2inner = S_inner*S_inner

    data_outer = pd.read_csv('LOP_outer.csv')
    S_outer = data_outer.mean(axis=0).loc['C12']
    S2outer = S_outer*S_outer
except:
    S2outer = 0.018487
    S2inner = 0.044275
    S2all = 0.026833  


# definition of a fitting fuction
def exp_fit1(x, w1, w2, tau1, tau2):
    return (1 - S2outer)*(w1 * np.exp(-1 * x / tau1) + w2 * np.exp(-1 * x / tau2)) + S2outer

def exp_fit2(x, w1, w2, tau1, tau2):
    return (1 - S2inner)*(w1 * np.exp(-1 * x / tau1) + w2 * np.exp(-1 * x / tau2)) + S2inner
    
def exp_fit3(x, w1, w2, tau1, tau2):
    return (1 - S2all)*(w1 * np.exp(-1 * x / tau1) + w2 * np.exp(-1 * x / tau2)) + S2all
    
# initial guess (given by vector p0) may differ
popt, pcov = curve_fit(exp_fit1, m[n:n2,0], m[n:n2,1], p0=(0.01, 0.01, 100, 100))
popt1, pcov1 = curve_fit(exp_fit2, m1[n:n2,0], m1[n:n2,1], p0=(0.01, 0.01, 100, 100))
popt2, pcov2 = curve_fit(exp_fit3, m2[n:n2,0], m2[n:n2,1], p0=(0.01, 0.01, 100, 100))

fig1 = plt.figure(figsize=(13,5),dpi=90)

# Plot original log10(ACF) data on the left graph
ax1 = fig1.add_subplot(121)
ax1.plot(m[:,0], np.log10(m[:,1]), color = 'b', lw = 1)
ax1.plot(m1[:,0], np.log10(m1[:,1]), color = 'r', lw = 1)
ax1.plot(m2[:,0], np.log10(m2[:,1]), color = 'g', lw = 1)
ax1.set_xlabel(u'Time, ps', fontsize=15)
ax1.set_ylabel('log(ACF)', fontsize=15)
ax1.annotate('A',(-0.15,0.95),xycoords='axes fraction',fontsize=16, fontweight='bold')

# Plot original ACF data + fitting on the right graph
ax2 = fig1.add_subplot(122)
ax2.plot(m[n:n2,0], m[n:n2,1], 'b-', label=u'Edge, tau=%5.2f' % (popt[0]*popt[2]+popt[1]*popt[3]))
ax2.plot(m[n:n2,0], exp_fit1(m[n:n2,0], *popt), 'k' )

ax2.plot(m1[n:n2,0], m1[n:n2,1], 'r-', label=u'Center, tau=%5.2f' % (popt1[0]*popt1[2]+popt1[1]*popt1[3]))
ax2.plot(m1[n:n2,0], exp_fit2(m1[n:n2,0], *popt1), 'k' )

ax2.plot(m2[n:n2,0], m2[n:n2,1], 'g-', label=u'All, tau=%5.2f' % (popt2[0]*popt2[2]+popt2[1]*popt2[3]) )
ax2.plot(m2[n:n2,0], exp_fit3(m2[n:n2,0], *popt2), 'k')
ax2.set_xlabel(u'Time, ps', fontsize=15)
ax2.set_ylabel('ACF', fontsize=15)
ax2.annotate('B',(-0.15,0.95),xycoords='axes fraction',fontsize=16, fontweight='bold')
plt.legend()
#plt.show()
plt.savefig('lipodisk_rotacf_new.png')
