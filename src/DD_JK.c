#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*subrutina para leer columnas
 * recibe un nombre de archivo (existente) y un largo (correr subr: N_lines)
 */
double* load_file(char input[64],int l){ // nombre de archivo y largo (numero de filas)
	double* col = malloc(sizeof(double) * l);
        FILE *file;
        printf("file: %s ... \n ",input); // archivo tiene el path completo
        file = fopen(input, "r");
        if(file==NULL) printf("Error... file not found");
        int d;
        for(d=0;d<l;d++){
                if (!fscanf(file, "%lf", &col[d])) // leer archivo (solo funciona para una columna)
           break;
        }
    fclose(file);
    printf("colum %s loaded \n",input);
    return col;
}

/*
 * subroutine to read the number of lines in a file
 * */
int N_lines(char filename[64]){
  // count the number of lines in the file called filename     
	printf("file: %s ...\n ",filename);
        FILE *fp = fopen(filename,"r");
        int ch=0;
        int lines=0;
        if (fp == NULL) printf("Error... file not found\n");
        while(!feof(fp)){
                ch = fgetc(fp);
                if(ch == '\n') lines++;
        }
        printf("number of rows in %s: %d \n",filename,lines);
        return lines;
}

double ** DD_2p(double * X,double * Y,double * Z,long int N,int Nl,float Lbox,int bins,double LIMIT,float dmin,float dmax,int JK){
	int JK_index,c,d,n;
	int NJ = JK*JK*JK; // numero total de JK (cubitos)
	double ** dd_JK = malloc(sizeof(double*)*NJ*2); // se lee el doble de memoria para calcular el numero de pares y pares pesados de forma aparte
	printf("Computing number of pairs for N of galaxies: %d\n",N);
	for(c=0;c<NJ;c++) 	dd_JK[c] = malloc(sizeof(double)*bins);
	for(c=0;c<NJ;c++) 	for(d=0;d<bins;d++) 	dd_JK[c][d] = 0.;
	int Ncell = Nl*Nl*Nl;
	float cell_size = Lbox/Nl;
	printf("Size of each cell: %f \n",cell_size);
	int lx,ly,lz;
	int jx,jy,jz;
	int * ll;
	int * ID_last;
	int first = -9; //para identificar ciclos vacios en las listas. 
	ll = malloc(sizeof(long int) * N); // ojo que si N es solo int se cae
	ID_last = malloc(sizeof(long int) * Ncell);
	for(c=0;c<Ncell;c++) ID_last[c] = -9; 
	// linked list
	for(c=0;c<N;c++){
		lx =(int) X[c] / cell_size;
		ly =(int) Y[c] / cell_size;
		lz =(int) Z[c] / cell_size;
		ID_last[lx*Nl*Nl+ly*Nl+lz] = c;
	}
	for(c=0;c<N;c++){
		lx =(int) X[c] / cell_size;
		ly =(int) Y[c] / cell_size;
		lz =(int) Z[c] / cell_size;
		ll[ID_last[lx*Nl*Nl+ly*Nl+lz]] = c;
		ID_last[lx*Nl*Nl+ly*Nl+lz] = c;
	}
	//
	int i,j,k,l,m,nk,nl,nm,a;
	int cell_max = ( (int) LIMIT/cell_size); //limite en celdas a leer dentro del loop
	printf("max separation between cells: %d \n",cell_max);
	double sep,dx,dy,dz;
	int index;
	double binsize_i = bins/( dmax - dmin);
	float Lbox2 = Lbox/2;
	printf("checking neighbours...\n");
	for(i=0;i<N;i++){
        if(i%50000==0) printf("i=%d \n",i); 
		lx = X[i] / cell_size;
		ly = Y[i] / cell_size;
		lz = Z[i] / cell_size;
		a = lx*Nl*Nl+ly*Nl+lz; // ID de lista para la particula i
		jx = (int) (X[i]/Lbox * JK);
		jy = (int) (Y[i]/Lbox * JK);
		jz = (int) (Z[i]/Lbox * JK);
		JK_index = jx*JK*JK+jy*JK+jz;// ID de JK para la particula i
		//printf("i en la celda %d,%d,%d es el ID %d \n",lx,ly,lz,a);
		for(k = lx - cell_max;k <= lx + cell_max;k++){ //loop in cells for list
            //printf("%d,",k);
	        	nk = k; // swap to avoid crash
			if(k<0) nk += Nl;
			if(k>Nl-1) nk -= Nl;
			for(l = ly - cell_max;l <= ly + cell_max;l++){
			    //printf("%d,",l);
		                nl = l;
				if(l<0) nl += Nl;
				if(l>Nl-1) nl -= Nl;
				for(m = lz - cell_max; m <= lz + cell_max;m++){
				    //printf("%d \n",m);
				    nm = m;
					if(m<0) nm += Nl;
					if(m>Nl-1) nm -= Nl;
					//printf("%d,%d,%d \n",nk,nl,nm);
					int dcell = (int) sqrt((k-lx)*(k-lx)+(l-ly)*(l-ly)+(m-lz)*(m-lz));
					//printf("distance to cell: %d and cell max: %d \n",dcell,cell_max);
					if(dcell<=cell_max){
						j = ID_last[nk*Nl*Nl+nl*Nl+nm];
						//printf("j=%d\n",j);
						while(j!= -9 && first != j){ // loop between 'vecinos' until read empty cell or point to itself (first particle in the list)
							dx = fabs( X[j] -  X[i]);
							dy = fabs( Y[j] -  Y[i]);
							dz = fabs( Z[j] -  Z[i]);
							if (dx > Lbox2)	dx = Lbox - dx;
							if (dy > Lbox2)	dy = Lbox - dy;
							if (dz > Lbox2)	dz = Lbox - dz;
							sep = 0.5*log10(dx*dx+dy*dy+dz*dz+1e-6); // separation in log(r) or r (linear)
							index = (int) ((sep-dmin)*binsize_i); // from sep to N bin
							if (index > -1 && index < bins){
								for(n=0;n<NJ;n++){
				                                    if(n!=JK_index){
									// first half for pairs (DD) second half for weighted pairs
									dd_JK[n][index] += 1.; // DD
	//								dd_JK[NJ+n][index] += V[i]*V[j]; // mm 
									// number of columns = N bins
									}
								}
							}
		                    		j = ll[j]; // point to next particle in the list. The last particle point to the first one and ends the loop
                            //printf("%d,",j);
						first = ID_last[nk*Nl*Nl+nl*Nl+nm];
                            //if (first==j) printf("\n last j = %d\n",j);
						}
					}
				}
			}
		}
	}
	return dd_JK; // matrix pointer
}

// define everything in the main and run functions inside please.
int main(int argc, char **argv){
	char PATH[64];  strcpy(PATH,argv[3]);	
	//printf("running mark correlation function\n");
//--------------function Parameters ------------------------ //
	int jk = 4; // number of jacknife per side. If are subboxes read nj=jk**3. If are slices (lonjas) of the box nj=jk. The function reads its own nj value, always read jk
        int nj = jk*jk*jk;
        float Bins = 15.; // N bins
        float lbox = 1024.;//largo caja
        float smin = -.3;// separaciones min y max (puede estar en lineal o log)
        float smax = 2.;
        int iN = 50; // numero de celdas para las liked list
        double limit = 102.4;// limit separation to read
        printf("N JK: %d\t Bins: %f\t L box: %f. In separation (log r): %f %f for %f bins \n",jk,Bins,lbox,smin,smax,Bins);
//----------------------------------------------------------//
	int c,d;
	//
	double* col_x; // asignar memoria al archivo de columnas
	double* col_y;
	double* col_z;
//	double* col_v;
	char Mod[2]; strcpy(Mod,argv[1]);
	char lab[2]; strcpy(lab,argv[2]);
	//printf("%s %s \n",Mod,lab);
	char name_x[64];// asignar nombre de archivo. ¡son columnas siempre!
	char name_y[64];
	char name_z[64];
//	char name_v[32];
	double** xi = malloc(sizeof(double*)*nj); 
	for(c=0;c<nj;c++) 	xi[c] = malloc(sizeof(double)*Bins);
	FILE *fp;
	long int length;
	// se puede crear un archivo y leer en char* argv para leer el modelo
        printf("loading files...\n");
       	snprintf(name_x,sizeof(char)*16,"%s_pos_x_%s.txt",Mod,lab); //3 columnas: x, y, z
        printf("%s created\n",name_x);
        snprintf(name_y,sizeof(char)*16,"%s_pos_y_%s.txt",Mod,lab);
        printf("%s created\n",name_y);
        snprintf(name_z,sizeof(char)*16,"%s_pos_z_%s.txt",Mod,lab);
        printf("%s created\n",name_z);
	char px[128]; char py[128]; char pz[128]; char pv[128];
        strcpy(px,PATH); strcpy(py,PATH); strcpy(pz,PATH);
	strcat(px,name_x); strcat(py,name_y); strcat(pz,name_z);
        length = N_lines(px); // asumiendo que todos los archivos tengan el mismo largo

	length = N_lines(px); // asumiendo que todos los archivos tengan el mismo largo. Tienen que..
	col_x = malloc(sizeof(double)*length);
	col_y = malloc(sizeof(double)*length);
	col_z = malloc(sizeof(double)*length);
//	col_v = malloc(sizeof(double)*length);
        col_x = load_file(px,length);
       	col_y = load_file(py,length);
        col_z = load_file(pz,length);
//        col_v = load_file(name_v,length);
        printf("calculating...\n");

       	xi = DD_2p(col_x,col_y,col_z,length,iN,lbox,Bins,limit,smin,smax,jk); //correr funcion
	printf("done!");
       	char NAME[32];
        snprintf(NAME, sizeof(char)*32, "DD_JK64_%s_%s.txt",Mod,lab); //archivo final. Cada linea es un JK. 
       	char OUT[256];
	strcpy(OUT,PATH);
        strcat(OUT,"pairs/"); // directory for the output file
       	strcat(OUT,NAME);

        fp = fopen(OUT, "wr");
       	printf("file %s created! Enjoy.\n",OUT);
        for(c=0;c<nj;c++){//ojo que tenga el mismo largo que en la funcion (el doble de lo definido)
       	    for(d=0;d<Bins;d++){
               	fprintf(fp,"%lf ",xi[c][d]);
            }
       	    fprintf(fp,"\n");
        }
	fclose(fp);

return 0;
}
