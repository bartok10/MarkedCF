# coding: utf-8
# %load 2
# %load MarkedCF/python/plot_MCF.py
import numpy as np
import pylab as plt
rmin = -0.3
rmax = 1.7
Nbins = 10
s = np.logspace(-0.3,1.7,Nbins)
# %load 3 4 5 29
GR = np.loadtxt('/home/cquezada/HOD_clusters/data/m_mark0.1_GR_logsmp3_2_b15_B1.txt')
F5 = np.loadtxt('/home/cquezada/HOD_clusters/data/m_mark0.1_F5_logsmp3_2_b15_B1.txt')
F6 = np.loadtxt('/home/cquezada/HOD_clusters/data/m_mark0.1_F6_logsmp3_2_b15_B1.txt')
# %load 18-22
f,ax = plt.subplots(2,1,sharex=True, gridspec_kw = {'height_ratios':[3,1.3]})
ax[0].plot(s,GR[:,0],'k+-',linewidth=2.0)
ax[0].plot(s,F5,'gs-',linewidth=2.0)
ax[0].plot(s,F6,'bo-',linewidth=2.0)
ax[1].plot(s,GR[:,0]/GR[:,0],'k-',linewidth=2.0)
ax[1].plot(s,F6/GR[:,0],'b-',linewidth=2.0)
ax[1].plot(s,F5/GR[:,0],'g-',linewidth=2.0)
ax[1].fill_between(s,(GR[:,0]-GR[:,1])/GR[:,0],(GR[:,0]+GR[:,1])/GR[:,0],facecolor='k' ,alpha=0.3)
ax[0].set_xscale('log')
ax[1].set_xlabel('$r$ [Mpc $h^{-1}$]',fontsize=20)
ax[0].set_ylabel('$\\mathcal{M}(r)$',fontsize=20)
ax[1].set_ylabel('$\\mathcal{M} / \\mathcal{M}_{GR} $',fontsize=20)
plt.tight_layout()
plt.subplots_adjust(hspace=0.0)
plt.show()
#plt.savefig('imagen.png')
