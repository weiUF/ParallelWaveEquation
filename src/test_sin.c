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
	double t, uxx, uyy;
	
	t = s->t;
	nx = s->nx;
	ny = s->ny;
	N = nx*ny;

	//calculate u, rhs
	for(int i=0;i<N; ++i){
		//s->u[i] = sin( t )*sin( s->x[i]) *sin(s->y[i] );
		uxx = sin(t)*sin(s->y[i])*sin(s->x[i]) *(-1);
		uyy = sin(t)*sin(s->y[i])*sin(s->x[i]) *(-1);
		s->rhs[i] = uxx + uyy;
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
