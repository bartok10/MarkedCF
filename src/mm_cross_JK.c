#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

double* load_file(char input[128],int l){ // nombre de archivo y largo (numero de filas)
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

int N_lines(char filename[128]){
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

double ** mm_cross(double * X1,double * Y1,double * Z1,double * V1,double * X2,double * Y2,double * Z2,double * V2,long int N1,long int N2,int Nl,float Lbox,int bins,double LIMIT,float dmin,float dmax,int JK){
	int JK_index,c,d,n;
	int NJ = JK*JK*JK; // numero total de JK (cubitos)
	double ** dd_JK = malloc(sizeof(double*)*NJ); // se lee el doble de memoria para calcular el numero de pares y pares pesados de forma aparte
	printf("N clusters: %d \t N galaxies: %d \n",N1,N2);
	for(c=0;c<NJ;c++) 	dd_JK[c] = malloc(sizeof(double)*bins);
	for(c=0;c<NJ;c++) 	for(d=0;d<bins;d++) 	dd_JK[c][d] = 0.; //si no se inicializan los punteros dentro de funciones da cualquier cosa.
	int Ncell = Nl*Nl*Nl;
	float cell_size = Lbox/Nl;
	printf("Size of each cell: %f \n",cell_size);
	int lx,ly,lz;
	int jx,jy,jz;
	int * ll;
	int * ID_last;
	int first = -9; //para identificar ciclos vacios en las listas. 
	ll = malloc(sizeof(long int) * N2); // ojo que si N es solo int se cae
	ID_last = malloc(sizeof(long int) * Ncell);
	for(c=0;c<Ncell;c++) ID_last[c] = -9; //Si se pone cero puede entrar al loop
	// linked list
	for(c=0;c<N2;c++){
		lx =(int) X2[c] / cell_size;
		ly =(int) Y2[c] / cell_size;
		lz =(int) Z2[c] / cell_size;
		ID_last[lx*Nl*Nl+ly*Nl+lz] = c;
	}
	for(c=0;c<N2;c++){
		lx =(int) X2[c] / cell_size;
		ly =(int) Y2[c] / cell_size;
		lz =(int) Z2[c] / cell_size;
		ll[ID_last[lx*Nl*Nl+ly*Nl+lz]] = c;
		ID_last[lx*Nl*Nl+ly*Nl+lz] = c;
	}
	//
	int i,j,k,l,m,nk,nl,nm,a;
	int cell_max = ( (int) LIMIT/cell_size) + 1; //limite en celdas a leer dentro del loop (una de mas siempre vale la pena)
	printf("max separation between cells: %d \n",cell_max);
	double sep,dx,dy,dz;
	int index;
	double binsize_i = bins/( dmax - dmin);
	float Lbox2 = Lbox/2;

	for(i=0;i<N1;i++){
        if(i%50000==0) printf("i=%d \n",i); 
		lx = X1[i] / cell_size;
		ly = Y1[i] / cell_size;
		lz = Z1[i] / cell_size;
		a = lx*Nl*Nl+ly*Nl+lz; // ID de lista para la particula i
		jx = (int) (X1[i]/Lbox * JK);
		jy = (int) (Y1[i]/Lbox * JK);
		jz = (int) (Z1[i]/Lbox * JK);
		JK_index = jx*JK*JK+jy*JK+jz;// ID de JK para la particula i
//		printf("i en la celda %d,%d,%d es el ID %d \n",lx,ly,lz,a);
		for(k = lx - cell_max;k <= lx + cell_max;k++){ //loop in cells for list
//			printf("%d,",k);
	        	nk = k; // swap to avoid crash
			if(k<0) nk += Nl;
			if(k>Nl-1) nk -= Nl;
			for(l = ly - cell_max;l <= ly + cell_max;l++){
//			    	printf("%d,",l);
				nl = l;
				if(l<0) nl += Nl;
				if(l>Nl-1) nl -= Nl;
				for(m = lz - cell_max; m <= lz + cell_max;m++){
//				    printf("%d \n",m);
					nm = m;
					if(m<0) nm += Nl;
					if(m>Nl-1) nm -= Nl;
//					printf("%d,%d,%d \n",nk,nl,nm);
					int dcell = (int) sqrt((k-lx)*(k-lx)+(l-ly)*(l-ly)+(m-lz)*(m-lz));
//					printf("distance to cell: %d and cell max: %d \n",dcell,cell_max);
					if(dcell<=cell_max){
						j = ID_last[nk*Nl*Nl+nl*Nl+nm];
//						printf("j=%d\n",j);
						while(j!= -9 && first != j){ // loop between 'vecinos' until read empty cell or point to itself (first particle in the list)
							dx = fabs( X1[i] -  X2[j]);
							dy = fabs( Y1[i] -  Y2[j]);
							dz = fabs( Z1[i] -  Z2[j]);
							if (dx > Lbox2)	dx = Lbox - dx;
							if (dy > Lbox2)	dy = Lbox - dy;
							if (dz > Lbox2)	dz = Lbox - dz;
							sep = 0.5*log10(dx*dx+dy*dy+dz*dz+1e-6); // separation in log(r) or r (linear)
							index = (int) ((sep-dmin)*binsize_i); // from sep to N bin
							if (index > -1 && index < bins){
								for(n=0;n<NJ;n++){
				                                    if(n!=JK_index){
									// first half for pairs (DD) second half for weighted pairs
	//								dd_JK[n][index] += 1.; // DD
									dd_JK[n][index] += V1[i]*V2[j]; // mm
									// number of columns = N bins
									}
								}
							}
		                    		j = ll[j]; // point to next particle in the list. The last particle point to the first one and ends the loop
//						printf("%d,",j);
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
	char PATH1[64],PATH2[64]; strcpy(PATH1,argv[4]); strcpy(PATH2,argv[5]);
//printf("running mark correlation function\n");
//--------------function Parameters ------------------------ //
//
	int jk = 4; // number of jacknife per side. If are subboxes read nj=jk**3. If are slices (lonjas) of the box nj=jk. The function reads its own nj value, always read jk
	int nj = jk*jk*jk;
	float Bins = 11.; // N bins
	int nbins = (int) Bins;
	float lbox = 1024.;//largo caja
	float smin = -.3;// separaciones min y max (puede estar en lineal o log)
	float smax = 1.8;
	int iN = 50; // numero de celdas para las liked list
	double limit = 64.0;// limit separation to read
	printf("N JK: %d\t Bins: %f\t L box: %f. In separation (log r): %f %f for %f bins \n",jk,Bins,lbox,smin,smax,Bins);
//----------------------------------------------------------//
	int c,d;
	double* col1_x,* col1_y,* col1_z,* col1_v,* col2_x,* col2_y,* col2_z,* col2_v;
	char Mod[2]; strcpy(Mod,argv[1]);
	char lab[2]; strcpy(lab,argv[2]);
	float i_m; sscanf(argv[3],"%f",&i_m);
	char mark[2]; strcpy(mark,argv[3]);
	printf("\n The power index: %f \n",i_m);
	//printf("%s %s \n",Mod,lab);
	char name_x[64],name_y[64],name_z[64],name_v[64];
	char p1x[128], p1y[128], p1z[128], p1v[128],p2x[128], p2y[128], p2z[128], p2v[128];
	strcpy(p1x,PATH1); strcpy(p1y,PATH1); strcpy(p1z,PATH1);strcpy(p1v,PATH1); strcpy(p2x,PATH2); strcpy(p2y,PATH2); strcpy(p2z,PATH2); strcpy(p2v,PATH2);
	double** xi = malloc(sizeof(double*)*nj); 
	for(c=0;c<nj;c++) 	xi[c] = malloc(sizeof(double)*Bins);
	FILE *fp;
	long int length1,length2;
	// se puede crear un archivo y leer en char* argv para leer el modelo
        printf("loading files...\n");
       	snprintf(name_x,sizeof(char)*64,"%s_pos_x_%s.txt",Mod,lab); //darle nombre a las columnas. Se necesitan 4: x, y, z, mark
        printf("%s created\n",name_x);
        snprintf(name_y,sizeof(char)*64,"%s_pos_y_%s.txt",Mod,lab);
        printf("%s created\n",name_y);
        snprintf(name_z,sizeof(char)*64,"%s_pos_z_%s.txt",Mod,lab);
        printf("%s created\n",name_z);
        snprintf(name_v,sizeof(char)*64,"marks/%s_mark_p%s_%s.txt",Mod,mark,lab);
	strcat(p1x,name_x); strcat(p1y,name_y); strcat(p1z,name_z); strcat(p1v,name_v); strcat(p2x,name_x); strcat(p2y,name_y); strcat(p2z,name_z); strcat(p2v,name_v);
	length1 = N_lines(p1x); length2 = N_lines(p2x);  // asumiendo que todos los archivos tengan el mismo largo
	col1_x = malloc(sizeof(double)*length1); col2_x = malloc(sizeof(double)*length2);
	col1_y = malloc(sizeof(double)*length1); col2_y = malloc(sizeof(double)*length2);
	col1_z = malloc(sizeof(double)*length1); col2_z = malloc(sizeof(double)*length2);
	col1_v = malloc(sizeof(double)*length1); col2_v = malloc(sizeof(double)*length2);
        col1_x = load_file(p1x,length1); col2_x = load_file(p2x,length2);
       	col1_y = load_file(p1y,length1); col2_y = load_file(p2y,length2);
        col1_z = load_file(p1z,length1); col2_z = load_file(p2z,length2);
        col1_v = load_file(p1v,length1); col2_v = load_file(p2v,length2);

        printf("calculating...\n");
     	xi = mm_cross(col1_x,col1_y,col1_z,col1_v,col2_x,col2_y,col2_z,col2_v,length1,length2,iN,lbox,Bins,limit,smin,smax,jk); //correr funcion
	printf("done!");
       	char NAME[36];
        snprintf(NAME, sizeof(char)*36, "m1m2_JK%d_%s_mark_p%s_box%s_b%d.txt",nj,Mod,mark,lab,nbins); //archivo final. Cada linea es un JK. Tiene 2 mitades: 1 para los pares, otra para pares marcados
       	char OUT[256];
        strcpy(OUT,PATH1); // directory for the output file
	strcat(OUT,"pairs/");
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
