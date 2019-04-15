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
	int l_ar=0;
	double initcond,bc_xl,bc_xr,bc_yb,bc_yt;
	//sol = (struct data *)malloc(sizeof(struct data));
	sol->nx = 100;
	sol->ny = 100; 
	sol->t = 0;
	sol->told = 0;
	sol->dt = 0.001;
    
    	int N = (sol->nx) * (sol->ny);
	int size = sizeof(double) * N;
    
	sol->u = (double *)malloc(size);
	sol->uold = (double *)malloc(size);
	sol->rhs = (double *)malloc(size);
	sol->x = (double *)malloc(size);
	sol->y = (double *)malloc(size);
	sol->IC = (double *)malloc(size);



	int ix,iy;
	double dx=(double)1/(nx-1),dy=(double)1/(ny-1);
    
    	printf("dx: %f ", dx);
 
        /* Manual Input has been toggled off for now to minimize time during test
	printf("\n Provide boundary condition at x=0\n");
	scanf("%lf", &bc_xl);
	printf("\n Provide boundary condition at x=1\n");
	scanf("%lf", &bc_xr);
	printf("\n Provide boundary condition at y=0\n");
	scanf("%lf", &bc_yb);
	printf("\n Provide boundary condition at y=1\n");
	scanf("%lf", &bc_yt);
//	printf("\n Provide initialization condition \n");
//	scanf("%lf", &initcond);

*/
        
	// init grid, u, uold, rhs
	for(int i=0;i<N; ++i){
		ix = i / ny;
		iy = i - ny * ix;
		sol->x[i] = ix * dx;
		sol->y[i] = iy * dy;
		sol->u[i] = 0;
		sol->uold[i] = 0;
		sol->rhs[i] = 0;
	}

		sol->u[4949] = 1;      // Initial condition as a point source

        
        
        
	/*  Manual Input has been toggled off for now to minimize time during test
	//BC for Top boundary xl
	for(int i=0;i<=2*ny;++i)
	{sol->u[i] = bc_xl; } //first nx elements of the array
	
	
	//BC for Bottom boundary xr
	l_ar=N-2ny;
	for(int i=l_ar;i<N;++i)
	{sol->u[i] = bc_xr; } //last nx elements of the array
	
	
	//BC for bottom boundary yb
	for(int i=0;i<nx;++i)
	{sol->u[i*ny] = bc_yb; } //Every nx_th element represent first column element
	
	
	//BC for top boundary yt
	for(int i=0;i<nx;++i)
	{sol->u[(i+1)*ny-1] = bc_yt; } //Every nx-1 element represent last column element

	*/

	return 0;
}
