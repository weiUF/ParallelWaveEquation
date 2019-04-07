//#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "data.h"

/* 
 Time stepping using central difference for $ u_{tt} = rhs $
 $$ u_{tt}^n = (u^{n+1} - 2 u^{n} + u^{n-1})/\delta t^2 = rhs^n $$
 rearrange and get
 $$ u^{n+1} =  \delta t^2 rhs^n + 2 u^{n} - u^{n-1}$$
 */


int timestepping(struct data *s){
	
	int nx,ny,N,i;
	double t,dt;
	double *unew;
	nx = s->nx;
	ny = s->ny;
	N = nx*ny;
	t = s->t;
	dt = s->dt;
	unew = (double*)malloc(N*sizeof(double));

	if (t > 0){
		for(i=0;i<N;++i){
			unew[i] = dt * dt * s->rhs[i] + 2 * s->u[i] - s->uold[i];
		}
	}else{
		//at the 1st step (need to use u' initial condtion)
		unew[i] = 0.5*(dt * dt * s->rhs[i] + 2 * s->u[i] + 2*dt*s->IC[i]); //IC=ut(t=0)
	}


	// uold = u
	for(i=0;i<N;++i){
		s->uold[i] = s->u[i];
	}
	// u = unew
	for(i=0;i<N;++i){
		s->u[i] = unew[i];
	}
	
	s->t += s->dt;

	free(unew);

	return 0;
}
