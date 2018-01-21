############### MarkedCF ############################

This is C program to compute pairs of galaxias and marked pairs of galaxies as function of the separation.

src/ contains the source files of 2 programs:

 - DD_JK.c: 	
	computes the DD  (galaxy-galaxy) term in the 2-point correlation function. The program recieves 2 arguments: arg[1]: a prefix as name of the catalog (e.g "GR" "F5" "F6" ), arg[2]: the suffix thas name a realization (a number of the box realization from 1 - 5). These arguments indicates a catalog that includes the positions of the data points to correlate (%s_pos_x_%, %s_pos_y_%s, %s_pos_z_%. For detailed information see the script column.py). The program also recieve several arguments to compute the number of pairs DD, these parameters are included in the main function of the program and would be included as a parameter file in the future. The parameters are:
		
	int jk = 4; // number of jacknife per side. If are subboxes read nj=jk**3. If are slices (lonjas) of the box nj=jk. The function reads its own nj value, always declare jk
	float Bins = 15.; // N bins (is a float to compute the separation between bins)
	float lbox = 1024.;//largo caja
	float smin = -.3;// min y max separations in log space
	float smax = 2.;
	int iN = 50; // number of cells in the liked list (see linked list section)
	double limit = 102.4;// limit separation to read (10% of Lbox by defintion)

	The function to compute DD in main is: 

		DD_2p(col_x,col_y,col_z,length,iN,lbox,Bins,limit,smin,smax,jk);

	where the col_%s (x,y,z) are declare by the function load_file(); and length is declare using the function N_lines();. The return of DD_2p(); is a double** with size (nj,bins) and is printed in the file (DD_JK64_%s_%s.txt",Mod,lab) in the selected folder, where nj = 64 by default, but can be changed.

 - marked_DD_JK:
	Also computes DD, but adding weights to every galaxy. These weights are defined as marks and correspond to a property that depends on the environment (could be the local density, the mass, the bright of a galaxy). The function is run as:

		mm_2p(col_x,col_y,col_z,col_m,length,iN,lbox,Bins,limit,smin,smax,jk);

	where col_m is the mark file of the galaxies. Every other parameter is the same as in DD_2p();. This function recieve an extra argument wich is the suffix of the mark file (i.e. a power index). The return is the file (mm_JK64_%s_%s.txt",Mod,lab)

	Run the utility functions are run as:
	
		load_file(file,int N_lines(file)); 

	where file is the full path of the file to use (a column file).


	(UPDATE January 2018) Cross-Correlation functions versions of the previous rutines are included in src/. The only differences with the auto-correlation versions is that 2 catalogs have to be included. col%d_%s ({1,2},{x,y,z}) are the arrays that includes the positions x,y,z of catalog 1 and 2 (e.g. catalog 1: clusters, catalog 2: galaxies). To avoid the overlap between the names of the catalogs (which are the same), two different folder paths need to be added (see example with HOD_clusters/ and HOD_galaxies). 

 
bin/ contain the following exe programs and should be run as:

	./DD GR 3 /home/jarmijo/HOD_galaxies/
	./mm GR 3 0.1
	./D1D2 GR 3 /home/jarmijo/HOD_clusters/ /home/jarmijo/HOD_galaxies/
	./m1m2 GR 3 0.1 /home/jarmijo/HOD_clusters/ /home/jarmijo/HOD_galaxies/

where GR is the name of the simulation (prefix) 3 is the box3 realization (suffix) 0.1 the name of the mark file (power index p = 0.1) and the path are folders that includes the positions and mark files.

python/ contains the scripts to compute the correlation function xi(r) and M(r) where r is the separation bin.

	- auto_f.py  : calculates RR and xi(r) given the pair file (see DD_2p(); function)).
	- chi_2_reduced.py  : computes the reduce (or not) chi_2 value for 2p or Marked CF.
	- column.py  : generates the column files given a catalog.
	- cross_f.py  : the same as auto_f.py for the cross-correlation version.
	- mark_list.py  : generates mark files given some property.
	- M_auto_f.py : the marked version of auto_f.py
	- M_cross_f.py  : the marked version of cross_f.py







