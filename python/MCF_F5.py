import numpy as np
import sys
Mod = 'F5'
box_i = str(1)

dir_1 = str(sys.argv[1]) #; dir_2 = '/home/jarmijo/HOD_galaxies/'
dir_2 = str(sys.argv[2])
#
print '\ndirectory:' + dir_1 + 'and' +dir_2+'\n'
mark = str(0.1)

m = np.loadtxt(dir_1+'marks/'+Mod+'_mark_p'+mark+'_'+box_i+'.txt')
mean_m = np.mean(m); print "mean mark (p = "+ mark +")  in F5 is:" + str(mean_m) + "\n"
mean_m_2 = np.mean(np.loadtxt(dir_2+'marks/F5_mark_p'+mark+'_'+box_i+'.txt'))
DD = np.loadtxt(dir_1+'pairs/D1D2_JK64_'+Mod+'_'+box_i+'.txt')
mm = np.loadtxt(dir_1+'pairs/m1m2_JK64_'+Mod+'_mark_p'+mark+'_'+box_i+'.txt')
d = mm/(DD*mean_m*mean_m_2)
M_mean = np.mean(d,axis=0)
n = len(d)
N = len(M_mean)
C = np.zeros((N,N))
for i in range(N):
       	for j in range(N):
               	for k in range(n):
                       	C[i][j] += (n-1)/float(n)*(d[k][i]- M_mean[i])*(d[k][j]-M_mean[j])
#JK_err = np.sqrt(C.diagonal())
S = M_mean
print 'save file:' + dir_1+'data/m_mark'+mark+'_F5_logsm0p3_1p7_B'+box_i+'.txt\n'
np.savetxt(dir_1+'data/m_mark'+mark+'_F5_logsm0p3_1p7_B'+box_i+'.txt',S)
print '...done.'


