import numpy as np
import os,sys

Mod=sys.argv[1]
out1 = '/gpfs/data/jarmijo/HOD_boxes_z0.0/'
#out2='/cosma/home/jarmijo/files/out_files/mass_mark_z0p5/'

ngal=np.zeros(0)
for i in range(1,6):
	ii = str(i)
	dat = np.loadtxt('/gpfs/data/FR/Baojiu/fr_data/B1024-PM1024/gadget_data/2015/HOD_catalogues/CMASS_fixed_n_xi/'+Mod+'/z0.0/HOD_'+Mod+'_B1024_Box'+ii+'_ex.dat',usecols=(1,2,3,4))
	print 'in file: /gpfs/data/FR/Baojiu/fr_data/B1024-PM1024/gadget_data/2015/HOD_catalogues/CMASS_fixed_n_xi/'+Mod+'/z0.0/HOD_'+Mod+'_B1024_Box'+ii+'_ex.dat'
	m = dat[:,0]
	x = dat[:,1]
	y = dat[:,2]
	z = dat[:,3]
	#ngal[i]=len(dat)
	#
	np.savetxt(out1+Mod+'_pos_x_'+ii+'.txt',x)
	np.savetxt(out1+Mod+'_pos_y_'+ii+'.txt',y)
	np.savetxt(out1+Mod+'_pos_z_'+ii+'.txt',z)
	np.savetxt(out1+Mod+'_mass_'+ii+'.txt',m)
	print '\ncolumns saved in: ' + out1
#np.savetxt(out2+'N_gal_'+Mod+'.txt',ngal)
