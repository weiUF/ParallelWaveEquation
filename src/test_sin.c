//#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include "data.h"

/* 
 generate sine wave for testing purpose
 sin(1/sqrt(2)*(x-y)sin(t)
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
		uxx = -0.5*sin(t)*sin( 1/sqrt(2) *( s->x[i] - s->y[i] ));
		uyy = -0.5*sin(t)*sin( 1/sqrt(2) *( s->x[i] - s->y[i] ));
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
		s->IC[i] = cos( t )*sin( 1/sqrt(2) *( s->x[i] - s->y[i] ));
	}



	return 0;
}
