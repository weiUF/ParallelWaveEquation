//#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include "data.h"

/* 
 generate sine wave for testing purpose
 sin(x)sin(y)sin(t)
 */

int test_sin(struct data *s){
	
	int nx,ny,N;
	double t, uxx, uyy,dt,dx,dy;
	double jj,kk,sx,sy,la;
	
	t = s->t;
	nx = s->nx;
	ny = s->ny;
	dt = s->dt;
	N = nx*ny;
	dx=(double)1/(nx-1);
	dy=(double)1/(ny-1);
	sx=(dt*dt)/(dx*dx);
	sy=(dt*dt)/(dy*dy);
	la=1-sx-sy;
	
	//calculate u, rhs
	for(int i=ny;i<N-ny; ++i){ //skip first and last rows, as BC is already available
		//s->u[i] = sin( t )*sin( s->x[i]) *sin(s->y[i] );
		jj=(double)(i%ny);//not to include first element of any row as this belongs to first column, BC available
		kk=(double)((i+1)%ny); //not to include last element of any row as this belong to last column, BC available 
		if(jj!=0.0 && kk!=0) //update RHS only for not boundary condition points
		{
			//if(t>=dt) //for second time-step onwards
			s->rhs[i] = (2.*s->u[i]*la)-(s->uold[i])+sx*((s->u[i+ny])+(s->u[i-ny]))+sy*((s->u[i+1])+(s->u[i-1]));
			//if(t<dt) //for first time-step
           // s->rhs[i] = (2.*s->u[i]*la)-(s->IC[i])+sx*((s->u[i+ny])+(s->u[i-ny]))+sy*((s->u[i+1])+(s->u[i-1]));
			//s->rhs[i] = (s->u[i]*la)+(dt*(s->IC[i]))+0.5*sx*((s->u[i+ny])+(s->u[i-ny]))+0.5*sy*((s->u[i+1])+(s->u[i-1]));			
		}

	}



	return 0;
/*
    
}


int test_sin_setIC(struct data *s){
	
	int nx,ny,N;
	double t, uxx, uyy;
	
	t = s->t;
	nx = s->nx;
	ny = s->ny;
	N = nx*ny;

	calculate u, rhs
	for(int i=0;i<N; ++i){
		s->IC[i] = sin( s->x[i]) *sin(s->y[i] );
	}*/


	return 0;
}
