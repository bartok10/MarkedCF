import numpy as np
import sys

i = 'F5'
ii = '1' #sys.argv[1] # box number (1,2,3,4,5)
dir1 = sys.argv[1]
dir2 = sys.argv[2]
#i=1
N1 = 4 # JKs per side
N3 = N1**3
NN = str(N3)
Lbox = 1024.
N = 10
l = Lbox/N1
smin = -0.3; smax = 1.7; ds = (smax - smin)/N
sep = np.logspace(-0.3,1.7,N)
x1 = np.loadtxt(dir1+i+'_pos_x_'+ii+'.txt')
y1 = np.loadtxt(dir1+i+'_pos_y_'+ii+'.txt')
z1 = np.loadtxt(dir1+i+'_pos_z_'+ii+'.txt')
x2 = np.loadtxt(dir2+i+'_pos_x_'+ii+'.txt')
y2 = np.loadtxt(dir2+i+'_pos_y_'+ii+'.txt')
z2 = np.loadtxt(dir2+i+'_pos_z_'+ii+'.txt')
Ngal1 = len(x1); Ngal2 = len(x2);
X = (x1/l).astype(int)
Y = (y1/l).astype(int)
Z = (z1/l).astype(int)
Npbox = np.zeros(N3,dtype=float)
for j in range(len(x1)):
    Npbox[X[j]*N1*N1+Y[j]*N1+Z[j]] += 1.
RR = np.zeros((N3,N))
for c in range(N3):
    for q in range(N):
        RR[c][q] = ((Ngal1-Npbox[c])*(Ngal2)/(1024.**3))*(4./3.)*np.pi*((sep[q]+ds)**3. - (sep[q])**3.)
DD = np.loadtxt(dir1+'pairs/D1D2_JK64_'+i+'_'+ii+'.txt')
xi = DD/RR - 1.
xi_mean = np.mean(xi,axis=0)
C = np.zeros((N,N))
for p in range(N):
    for q in range(N):
        for k in range(N3):
            C[p][q] += (N3-1)/float(N3)*(xi[k][p]- xi_mean[p])*(xi[k][q]-xi_mean[q])
JK_err = np.sqrt(C.diagonal())
S = np.array([xi_mean,JK_err])
np.savetxt(dir1+'data/CCF_'+i+'_logsm0p3_1p7_B'+ii+'.txt',S.T)
print 'file' + dir1+'data/CCF_'+i+'_logsm0p3_1p7_B'+ii+'.txt' ' saved!'

