//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "data.h"

/*
	initialize data structure and create grid points
 */

int init(struct data *sol){
	
	// init	
	int nx=100, ny=100;
	//sol = (struct data *)malloc(sizeof(struct data));
	sol->nx = 100;
	sol->ny = 100;
	sol->t = 0;
	sol->told = 0;
	sol->dt = 0.0001;
	int size = sizeof(double) * sol->nx * sol->ny;
	sol->u = (double *)malloc(size);
	sol->uold = (double *)malloc(size);
	sol->rhs = (double *)malloc(size);
	sol->x = (double *)malloc(size);
	sol->y = (double *)malloc(size);
	sol->IC = (double *)malloc(size);

	int N = sol->nx * sol->ny;

	int ix,iy;
	double dx=(double)1/(nx-1),dy=(double)1/(ny-1);

	// init grid, u, uold, rhs
	for(int i=0;i<N; ++i){
		ix = i / nx;
		iy = i - nx * ix;
		sol->x[i] = ix * dx;
		sol->y[i] = iy * dy;
		sol->u[i] = 0;
		sol->uold[i] = 0;
		sol->rhs[i] = 0;
		sol->IC[i] = 0;
	}


	return 0;
}
