import numpy as np
import sys 

dir1 = str(sys.argv[2])
#dir2 = '/home/jarmijo/HOD_clusters/marks/'
#dir1 = '/gpfs/data/jarmijo/'
#dir2 = '/gpfs/data/jarmijo/mark_dens_HOD/'
Mod = sys.argv[1]
p_index_m = np.array([0.1])
for q in range(1,6):
	box_n = str(q)
	masss = np.loadtxt(dir1+Mod+'_mass_'+box_n+'.txt')
	m_m = np.random.normal(masss,0.2*masss)
	for i in (p_index_m):
		marks = (m_m)**i
		marks[np.isinf(marks)] = 1.0; marks[np.isnan(marks)] = 1.0;
		print 'marks saved:' + dir1+Mod+'_mark_p'+str(i)+'_'+box_n+'.txt \n'
		np.savetxt(dir1+'marks/'+Mod+'_mark_p'+str(i)+'_'+box_n+'.txt',marks)
		print "...done!\n"	

