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
	if(s->myrank != s->mpi_size)
	{MPI_Send(&(s->u[myN-2*(ny+2)]), ny+2, MPI_DOUBLE, s->myrank+1, 1, MPI_COMM_WORLD);}

	if(s->myrank != 0)
	{MPI_Recv(&(s->u[0]), ny+2, MPI_DOUBLE, s->myrank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);}
	if(s->myrank != s->mpi_size)
	{MPI_Recv(&(s->u[myN-ny]), ny+2, MPI_DOUBLE, s->myrank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);}



	//calculate u, rhs
	for(int i=0;i<myN; ++i){
		ilocal = i + nx + 3 + i / ny * 2; // this index exclude ghosh & boundary points
		if(t>=dt) //for second time-step onwards
		s->rhs[i] = (2.*s->u[i]*la)-(s->uold[i])+sx*((s->u[i+ny])+(s->u[i-ny]))+sy*((s->u[i+1])+(s->u[i-1]));
		if(t<dt) //for first time-step
		s->rhs[i] = (s->u[i]*la)+(dt*(s->IC[i]))+0.5*sx*((s->u[i+ny])+(s->u[i-ny]))+0.5*sy*((s->u[i+1])+(s->u[i-1]));			

	}



	return 0;
}

int test_sin_setIC(struct data *s){
	
	int nx,ny,N;
	double t, uxx, uyy;
	
	t = s->t;
	nx = s->nx;
	ny = s->ny;
	N = nx*ny;

	//calculate u, rhs
	for(int i=0;i<N; ++i){
		s->IC[i] = cos( t )*sin( s->x[i]) *sin(s->y[i] );
	}

	return 0;
}
