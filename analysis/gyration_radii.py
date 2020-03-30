# -*- coding: utf-8 -*-

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np

aliphatic = np.loadtxt('ST_gyr.xvg', comments=['#','@','&'])
mal = np.loadtxt('MA_gyr.xvg', comments=['#','@','&'])
both = np.loadtxt('both_gyr.xvg', comments=['#','@','&'])

#plt.plot(diib[:,0]/1000, diib[:,1]*10, color = 'orange', lw = 1, label = u'Diisobutylene')
#plt.plot(mal[:,0]/1000, mal[:,1]*10, color = 'g', lw = 1, label = u'Maleic acid')
#plt.plot(dibma[:,0]/1000, dibma[:,1]*10, color = 'b', lw = 1, label = 'DIBMA')
plt.plot(aliphatic[:,0]/1000, aliphatic[:,1]*10, color = 'orange', lw = 1, label = u'Styrol')
plt.plot(mal[:,0]/1000, mal[:,1]*10, color = 'g', lw = 1, label = u'Maleic acid')
plt.plot(both[:,0]/1000, both[:,1]*10, color = 'b', lw = 1, label = 'Both')

#plt.xlabel(u"Time, ns")
#plt.ylabel(u"Gyration radius, Å")
plt.xlabel(u"Time, ns")
plt.ylabel(u"Gyration radii, Å")
plt.legend()
#plt.show()
#plt.savefig('DIBMALP_gyr_radii.png')
plt.savefig('gyr_radii.png')
