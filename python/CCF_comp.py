# coding: utf-8
# %load MarkedCF/python/CCF_comp.py
import numpy as np
import pylab as plt
import sys 
rmin = -0.3
rmax = 1.7
Nbins = 10
dir1 = "/home/cquezada/HOD_clusters2/"
dir2 = "/home/cquezada/HOD_clusters/"
s = np.logspace(-0.3,1.7,Nbins)
GR = np.loadtxt(dir2+'data/CCF_GR_logsm0p3_1p7_B1.txt')
F5 = np.loadtxt(dir1+'data/CCF_F5_logsm0p3_1p7_B1.txt')
MGR = np.loadtxt(dir2+'data/m_mark0.1_GR_logsm0p3_1p7_B1.txt')
MF5 = np.loadtxt(dir1+'data/m_mark0.1_F5_logsm0p3_1p7_B1.txt')
f,ax = plt.subplots(1,1,figsize=(6,3))
ax.plot(s,GR[:,0]/GR[:,0],'k+-',linewidth=2.0)
ax.plot(s,F5[:,0]/GR[:,0],'gs-',linewidth=2.0)
ax.fill_between(s,(GR[:,0]-GR[:,1])/GR[:,0],(GR[:,0]+GR[:,1])/GR[:,0],facecolor='k' ,alpha=0.3)
ax.set_xscale('log')
ax.set_xlabel('$r$ [Mpc $h^{-1}$]',fontsize=20)
ax.set_ylabel('$\\xi / \\xi_{GR} $',fontsize=20)
plt.show()
