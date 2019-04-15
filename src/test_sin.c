#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include "data.h"

/* 
 generate sine wave for testing purpose
 sin(x)sin(y)sin(t)
 */

int test_sin(struct data *s){
	
	int nx,ny,N,myN,ilocal;
	double t, uxx, uyy,dt,dx,dy;
	double jj,kk,sx,sy,la;
	
	t = s->t;
	nx = s->nx;
	ny = s->ny;
	dt = s->dt;
	N = nx*ny;
	myN = s->myN;
	dx=(double)1/(nx+1);
	dy=(double)1/(ny+1);
	sx=(dt*dt)/(dx*dx);
	sy=(dt*dt)/(dy*dy);
	la=1-sx-sy;
	
	//send & receive data from adjacent region/processor
	if(s->myrank != 0)
	{MPI_Send(&(s->u[ny+2]), ny+2, MPI_DOUBLE, s->myrank-1, 0, MPI_COMM_WORLD);}
	if(s->myrank != s->mpi_size-1)
	{MPI_Send(&(s->u[myN-2*(ny+2)]), ny+2, MPI_DOUBLE, s->myrank+1, 1, MPI_COMM_WORLD);}

	if(s->myrank != 0)
	{MPI_Recv(&(s->u[0]), ny+2, MPI_DOUBLE, s->myrank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);}
	if(s->myrank != s->mpi_size-1)
	{MPI_Recv(&(s->u[myN-ny]), ny+2, MPI_DOUBLE, s->myrank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);}



	//calculate u, rhs
	for(int i=0;i<myN; ++i){
		ilocal = i + ny + 3 + i / ny * 2; // this index skips ghosh & boundary points
		//debug
		if (ilocal >= (s->mynx+2)*(ny+2)){
			printf("idx: %d VS max: %d",ilocal, (s->mynx+2)*(ny+2));
			fflush(stdout);}
		if(t>=dt) //for second time-step onwards
		s->rhs[ilocal] = (2.*s->u[ilocal]*la)-(s->uold[ilocal])+sx*((s->u[ilocal+ny+2])+(s->u[ilocal-ny-2]))+sy*((s->u[ilocal+1])+(s->u[ilocal-1]));
		if(t<dt) //for first time-step
		//s->rhs[ilocal] = (s->u[ilocal]*la)+(dt*(s->IC[ilocal]))+0.5*sx*((s->u[ilocal+ny+2])+(s->u[ilocal-ny-2]))+0.5*sy*((s->u[ilocal+1])+(s->u[ilocal-1]));			
		s->rhs[ilocal] = (s->u[ilocal]*la)+0.5*sx*((s->u[ilocal+ny+2])+(s->u[ilocal-ny-2]))+0.5*sy*((s->u[ilocal+1])+(s->u[ilocal-1]));			

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
