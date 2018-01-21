import numpy as np
import sys
MG = ['F5','F6']

m_type = str(sys.argv[1])
box = str(sys.argv[2])
if m_type == 'mass': prefix='m';dir_in = '/gpfs/data/jarmijo/mark_mass_HOD/'; dir1='/gpfs/data/jarmijo/mass_marked_pairs/'; dir2 = 'files/out_files/pairs/';

elif m_type == 'rho': prefix='d'; dir_in = '/gpfs/data/jarmijo/mark_dens_HOD/'; dir1 = '/gpfs/data/jarmijo/dens_marked_pairs/'; dir2 = 'files/out_files/pairs_b12/';
else: print 'marks available: mass, rho(density)'

#dir_in = '/gpfs/data/jarmijo/mark_dens_HOD/'
#dir_out = 'files/out_files/new_mcf/'
#dir1 = '/gpfs/data/jarmijo/dens_marked_pairs/'
#dir2 = 'files/out_files/pairs_b12/'
print '\ndirectories:' + dir_in +',' +  dir1 +',' + dir2 +'\n'
mark = ['-0.1','0.1']
chi_M = np.zeros((len(mark),len(MG)),dtype=float)
#if m_type == 'mass': prefix='m'
#elif m_type = 'rho': prefix='d'
#else: print 'marks available: mass, rho(density)' 
	
#
#mark = ['-0.001','0.001','-1.0']
for p in range(len(mark)):
	m = np.loadtxt(dir_in+'GR_mark8_p'+mark[p]+'_'+box+'.txt')
	mean_m_GR = np.mean(m)
	DD_GR = np.loadtxt(dir2+'DD_JK64_GR_'+box+'.txt')
	mm_GR=np.loadtxt(dir1+'mm_JK64_GR_mark8_p'+mark[p]+'_'+box+'.txt')
	dGR = mm_GR/(DD_GR*mean_m_GR**2.)
	GR_mean = np.mean(dGR,axis=0)
	n = len(dGR)
	N = len(GR_mean)
	C = np.zeros((N,N))
	for i in range(N):
        	for j in range(N):
                	for k in range(n):
                        	C[i][j] += (n-1)/float(n)*(dGR[k][i]- GR_mean[i])*(dGR[k][j]-GR_mean[j])
	JK_err_GR = np.sqrt(C.diagonal())
	S = np.vstack((GR_mean,JK_err_GR))
	print 'save file:' + dir1+'data/m_mark'+mark[p]+'_GR_logsmp3_2_b15_B'+box+'.txt\n'
	np.savetxt(dir1+'data/m_mark8'+mark[p]+'_GR_logsmp3_2_b15_B'+box+'.txt',S.T)
	print '...done.'
	V_C = np.zeros((N,N))
	Co_i = np.linalg.inv(C)
	for i in range(len(C)):
        	for j in range(len(C)):
                	V_C[i][j] = C[i][j]/(C[i][i]*C[j][j])**0.5
#	np.savetxt(dir_out+'covM_GR_mark'+mark[p]+'.txt',C)
#	np.savetxt(dir_out+'covM_GR_mark'+mark[p]+'_norm.txt',V_C)

	for q in range(len(MG)):
		m = np.loadtxt(dir_in+MG[q]+'_mark8_p'+mark[p]+'_'+box+'.txt')
		mean_m_MG = np.mean(m)
		DD_MG = np.loadtxt(dir2+'DD_JK64_'+MG[q]+'_'+box+'.txt')
	        mm_MG=np.loadtxt(dir1+'mm_JK64_'+MG[q]+'_mark8_p'+mark[p]+'_'+box+'.txt')
        	dMG = mm_MG/(DD_MG*mean_m_MG**2.)
		MG_mean = np.mean(dMG,axis=0)
		S = MG_mean
		print 'save file:' + dir1+'data/m_mark'+mark[p]+'_'+MG[q]+'_logsmp3_2_b15_B'+box+'.txt'
		np.savetxt(dir1+'data/m_mark8'+mark[p]+'_'+MG[q]+'_logsmp3_2_b15_B'+box+'.txt',S.T)
		print '...done.'
		delta = GR_mean - MG_mean
		chi_2=0.

		for i in range(len(C)):
			for j in range(len(C)):
				chi_2 += delta[i]*(Co_i[i][j])*delta[j]
#			CHI[i] = chi_2
#			chimp = np.dot(Co_i,delta)
		
		chi_2_r = chi_2 / len(C)
		chi_M[p][q] = chi_2_r
p_float = np.array(mark,dtype=float)
print p_float,chi_M
S = np.vstack([p_float.T,chi_M.T])
print 'chi_2 saved in: ' + './chi_2_red_B'+box+'_M8.txt'
np.savetxt('./chi_2_red_B'+box+'_M8.txt',S.T,header='rows: marks (see readme)\t cols: F5,F6')

