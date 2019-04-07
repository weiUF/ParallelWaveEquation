//#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include "data.h"

/* 
 generate sine wave for testing purpose
 sin(xyt)
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
		s->u[i] = sin( t * s->x[i] * s->y[i] );
		uxx = -(s->y[i]*t*s->y[i]*t)*sin(s->x[i]*s->y[i]*t);
		uyy = -(s->x[i]*t*s->x[i]*t)*sin(s->x[i]*s->y[i]*t);
		s->rhs[i] = uxx + uyy;
	}

	return 0;
}

