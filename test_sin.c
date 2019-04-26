#include <mpi.h>
#include <stdio.h>
#include <math.h>
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
	dx=(double)1/(nx+1);
	dy=(double)1/(ny+1);
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
	
	int nx,ny,N,ilocal;
	double t,uxx,uyy,myN;
	
	t = s->t;
	nx = s->nx;
	ny = s->ny;
	myN = s->myN;

	//calculate u, rhs
	for(int i=0;i<myN; ++i){
		ilocal = i + ny + 3 + i / ny * 2; // this index exclude ghosh & boundary points
		s->u[ilocal] = cos( t )*sin( s->x[ilocal]) *sin(s->y[ilocal] );
	}
	return 0;
}