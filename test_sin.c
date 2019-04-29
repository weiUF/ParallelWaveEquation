#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include "data.h"

/* 
 generate sine wave for testing purpose
 sin(x)sin(y)sin(t)
 */

int test_sin(struct data *s){
	
	int nx,ny,N,myN,ilocal,i,bc_ilocal,myN_modif,mynx;
	double t, uxx, uyy,dt,dx,dy;
	double jj,kk,sx,sy,la;
	
	t = s->t;
	nx = s->nx;
	ny = s->ny;
	dt = s->dt;
	N = (s->nx+2)*(ny+2);
	myN = s->myN;
	dx=(double)(s->x_xl)/(nx+1);
	dy=(double)(s->y_yl)/(ny+1);
	sx=(dt*dt)/(dx*dx);
	sy=(dt*dt)/(dy*dy);
	la=1-sx-sy;
	mynx=s->mynx;
	
	//send & receive data from adjacent region/processor
	i = 0; // ix=0 points
	ilocal = i + ny + 3 + i / ny * 2; // this index skips ghosh & boundary points
	if(s->myrank != 0)
	{MPI_Send(&(s->u[ilocal]), ny, MPI_DOUBLE, s->myrank-1, 0, MPI_COMM_WORLD);}
	
	i = myN-ny; // ix=mynx points
	ilocal = i + ny + 3 + i / ny * 2; // this index skips ghosh & boundary points
	if(s->myrank != s->mpi_size-1)
	{MPI_Send(&(s->u[ilocal]), ny, MPI_DOUBLE, s->myrank+1, 1, MPI_COMM_WORLD);}

	if(s->myrank != 0) //ix=-1 ghost points
	{MPI_Recv(&(s->u[1]), ny, MPI_DOUBLE, s->myrank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);}

	i = myN+1; // ix=mynx+1 ghost points
	ilocal = i + ny + 2 + i / ny * 2; //:corrected index 
	if(s->myrank != s->mpi_size-1)
	{MPI_Recv(&(s->u[ilocal]), ny, MPI_DOUBLE, s->myrank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);}

	if (s->myrank == 0)
	{ 
		myN_modif=(myN/(mynx*ny))*(ny-2)*(mynx-1);		
		//calculate u, rhs
		for(i=0;i<myN_modif; ++i){
			ilocal = i + 2*ny + 6 + i / (ny-2) * 4; //modified to skip 2nd and 2nd last element of each column -> these are BC nodes
			//debug
			if (ilocal >= (s->mynx+2)*(ny+2)){
				printf("idx: %d VS max: %d",ilocal, (s->mynx+2)*(ny+2));
				fflush(stdout);}
			bc_ilocal=ilocal%(ny+2);// to exclude bc point and outer-ghost cell from iteration, i.e. two cells every column
			if(bc_ilocal>=2 && bc_ilocal<ny) {
				if(t>=dt) //for second time-step onwards
				s->rhs[ilocal] = (2.*s->u[ilocal]*la)-(s->uold[ilocal])+sx*((s->u[ilocal+ny+2])+(s->u[ilocal-ny-2]))+sy*((s->u[ilocal+1])+(s->u[ilocal-1]));
				if(t<dt) //for first time-step
				//s->rhs[ilocal] = (s->u[ilocal]*la)+(dt*(s->IC[ilocal]))+0.5*sx*((s->u[ilocal+ny+2])+(s->u[ilocal-ny-2]))+0.5*sy*((s->u[ilocal+1])+(s->u[ilocal-1]));			
				s->rhs[ilocal] = (s->u[ilocal]*la)+0.5*sx*((s->u[ilocal+ny+2])+(s->u[ilocal-ny-2]))+0.5*sy*((s->u[ilocal+1])+(s->u[ilocal-1]));	
			}//else
			//{s->rhs[ilocal] = s->u[ilocal];}

		}

	
	}
	else if(s->myrank == s->mpi_size-1)
	{
		myN_modif=(myN/(mynx*ny))*(ny-2)*(mynx-1);		
		//calculate u, rhs
		for(i=0;i<myN_modif; ++i){
			ilocal = i + ny + 4 + i / (ny-2) * 4; //modified to skip 2nd and 2nd last element of each column -> these are BC nodes
			//debug
			if (ilocal >= (s->mynx+2)*(ny+2)){
				printf("idx: %d VS max: %d",ilocal, (s->mynx+2)*(ny+2));
				fflush(stdout);}
			bc_ilocal=ilocal%(ny+2);// to exclude bc point and outer-ghost cell from iteration, i.e. two cells every column
			if(bc_ilocal>=2 && bc_ilocal<ny) {
				if(t>=dt) //for second time-step onwards
				s->rhs[ilocal] = (2.*s->u[ilocal]*la)-(s->uold[ilocal])+sx*((s->u[ilocal+ny+2])+(s->u[ilocal-ny-2]))+sy*((s->u[ilocal+1])+(s->u[ilocal-1]));
				if(t<dt) //for first time-step
				//s->rhs[ilocal] = (s->u[ilocal]*la)+(dt*(s->IC[ilocal]))+0.5*sx*((s->u[ilocal+ny+2])+(s->u[ilocal-ny-2]))+0.5*sy*((s->u[ilocal+1])+(s->u[ilocal-1]));			
				s->rhs[ilocal] = (s->u[ilocal]*la)+0.5*sx*((s->u[ilocal+ny+2])+(s->u[ilocal-ny-2]))+0.5*sy*((s->u[ilocal+1])+(s->u[ilocal-1]));	
			}//else
			//{s->rhs[ilocal] = s->u[ilocal];}

		}

	}
	else
	{
		//calculate u, rhs
		for(i=0;i<myN; ++i){
			ilocal = i + ny + 3 + i / ny * 2; // this index skips ghosh & boundary points
			//debug
			if (ilocal >= (s->mynx+2)*(ny+2)){
				printf("idx: %d VS max: %d",ilocal, (s->mynx+2)*(ny+2));
				fflush(stdout);}
			bc_ilocal=ilocal%(ny+2);// to exclude bc point and outer-ghost cell from iteration, i.e. two cells every column
			if(bc_ilocal>=2 && bc_ilocal<ny) {
				if(t>=dt) //for second time-step onwards
				s->rhs[ilocal] = (2.*s->u[ilocal]*la)-(s->uold[ilocal])+sx*((s->u[ilocal+ny+2])+(s->u[ilocal-ny-2]))+sy*((s->u[ilocal+1])+(s->u[ilocal-1]));
				if(t<dt) //for first time-step
				//s->rhs[ilocal] = (s->u[ilocal]*la)+(dt*(s->IC[ilocal]))+0.5*sx*((s->u[ilocal+ny+2])+(s->u[ilocal-ny-2]))+0.5*sy*((s->u[ilocal+1])+(s->u[ilocal-1]));			
				s->rhs[ilocal] = (s->u[ilocal]*la)+0.5*sx*((s->u[ilocal+ny+2])+(s->u[ilocal-ny-2]))+0.5*sy*((s->u[ilocal+1])+(s->u[ilocal-1]));	
			}//else
			//{s->rhs[ilocal] = s->u[ilocal];}

		}
	}





	return 0;
}

int test_sin_setIC(struct data *s){
	
	int nx,ny,ilocal,i,j,size,k,rank_int,iglobal,ilocal_fz,s_p,f_r;
	float x_part,dx,dy,rad;
	float partsz,mynx,rmin;
	float *x_pt,*y_pt,*u_pt;
	FILE *fpr;
//variable defined as double were creating memory problem and thus segmentation violation error. To resolve this all doubles in this section have been converted to float

	nx = s->nx;
	ny = s->ny;
	rank_int=-1;
	x_part=(s->x_xl)/s->mpi_size;
	ilocal=0;
	rmin=100000.0;
	f_r=0;
	if (s->myrank == 0)
	{
		printf("\n Data using text file: Yes=1 , No=0 \n");
		fflush(stdout);
		scanf("%d", &f_r);
		if(f_r==1)
		{
			printf("\n Reading Point Initialization Data from 'point_disturbance.txt' \n");
			if ((fpr = fopen("point_disturbance.txt","r")) == NULL)
				{
				printf("Error! opening file");
				exit(1);
				}

			fscanf(fpr,"%d", &s_p);
			size = sizeof(float) *s_p;
			x_pt = (float *)malloc(size);
			y_pt = (float *)malloc(size);
			u_pt = (float *)malloc(size);
			
			for(int i=0;i<s_p; ++i){
				fscanf(fpr,"%f %f %f", &x_pt[i],&y_pt[i],&u_pt[i]);
			}
			fclose(fpr); 
			printf("\n Point Initialization Data reading complete  \n");
		}
		else
		{		 
			printf("\n Number of pointwise disturbance  \n");
			fflush(stdout);
			scanf("%d", &s_p);
			//double x_pt[8]= { 0.5, 0.9, 0.4, 0.2, 0.7, 0.1, 0.3, 0.65};
			//double y_pt[8]={ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
			//double u_pt[8]={ 5, 5, 5, 5, 5, 5, 5, 5 };
			printf("inside first loop value of x_pt= %f \n ",x_pt[1]);;	
			size = sizeof(float) *s_p;
			x_pt = (float *)malloc(size);
			y_pt = (float *)malloc(size);
			u_pt = (float *)malloc(size);
			for(int i=0;i<s_p; ++i){
				printf("\n Enter x y and U value for ith element %d \n",i);
				scanf("%f %f %f", &x_pt[i],&y_pt[i],&u_pt[i]);
			}

		}
		
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&s_p,1,MPI_INT,0,MPI_COMM_WORLD);
		
	
	for(int k=0;k<s_p; ++k)
	{
		
		if (s->myrank == 0)
		{ 
			
			for(int i=0;i<s->mpi_size; ++i){
				if((i*x_part <= x_pt[k])&&(x_pt[k]<(i+1)*x_part))
				{
					rank_int=i;
				}		
			}
			printf("identified rank is %d\n",rank_int);
			MPI_Send(&x_pt[0], s_p, MPI_FLOAT, rank_int, 10, MPI_COMM_WORLD);
			MPI_Send(&y_pt[0], s_p, MPI_FLOAT, rank_int, 11, MPI_COMM_WORLD);
			MPI_Send(&u_pt[0], s_p, MPI_FLOAT, rank_int, 12, MPI_COMM_WORLD);
			/*for(int i=0;i<s_p; ++i){
				printf("\n node-0 : X-co-ordinates of entered points are = %f ",x_pt[i]);
				printf("\n node-0 : Y-co-ordinates of entered points are = %f ",y_pt[i]);
				printf("\n node-0 : U-value of entered points are = %f ",u_pt[i]);				
			}*/
			
			
		}
		MPI_Bcast(&rank_int,1,MPI_INT,0,MPI_COMM_WORLD);


		//printf("accessed rank from outside is %d\n",rank_int);
		MPI_Barrier(MPI_COMM_WORLD);
		if (s->myrank == rank_int)
		{
			size = sizeof(float) *s_p;
			x_pt = (float *)malloc(size);
			y_pt = (float *)malloc(size);
			u_pt = (float *)malloc(size);
			MPI_Recv(&x_pt[0], s_p, MPI_FLOAT, 0, 10, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&y_pt[0], s_p, MPI_FLOAT, 0, 11, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			MPI_Recv(&u_pt[0], s_p, MPI_FLOAT, 0, 12, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			/*for(int i=0;i<s_p; ++i){
				printf("\n after receiving : X-co-ordinates of entered points are = %f ",x_pt[i]);
				printf("\n after receiving : Y-co-ordinates of entered points are = %f ",y_pt[i]);
				printf("\n after receiving : U-value of entered points are = %f ",u_pt[i]);				
			}*/
			printf("accessed rank is %d\n",rank_int);
			partsz=rank_int*(s->x_xl)/s->mpi_size;		
			dx=(float)(s->x_xl)/(nx-1);
			dy=(float)(s->y_yl)/(ny-1);
			mynx=(s->nx)/(s->mpi_size);

			for(int i=1;i<=mynx; ++i){//leaving first row of ghost cells

				if((((i-1)*dx+partsz)<= x_pt[k])&&(x_pt[k]<((i)*dx+partsz)))//checking row interval
				{
						
					for(int j=0;j<2*(ny-2); ++j) //checking all the columns between the two rows, excluding coumns belonging to boundaries
					{
						ilocal=i*(ny+2)+2+j+j/(ny-2)*4;				
						//iglobal=ilocal+(ny+2)*s->myrank*mynx;
						rad=sqrt((x_pt[k]-(s->x[ilocal]))*(x_pt[k]-(s->x[ilocal]))+(y_pt[k]-(s->y[ilocal]))*(y_pt[k]-(s->y[ilocal])));
						if(rmin>=rad)
						{
							ilocal_fz=ilocal;
							rmin=rad;
						}
					}
					
				}
			
			}
			s->u[ilocal_fz]=u_pt[k];
			//printf("");//);
			printf("accessed rank is %d\n X-coordinate from list %f\n accessed U location X-coordinate & Y-coordinate %f & %f\n accessed U location %d and ilocal value %d \n",rank_int,x_pt[k],s->x[ilocal_fz],s->y[ilocal_fz],ilocal_fz,ilocal);
			//printf("");
				
		}
	}


	return 0;
}

