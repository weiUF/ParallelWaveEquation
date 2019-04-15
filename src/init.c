#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "data.h"

/*
	initialize data structure and create grid points
 */

int init(struct data *sol){
	
	// init	
	int nx=100, ny=100; // no. of internal points
	int N, mynx, istart, iend, myN, size, ix, iy, l_ar; 
	double initcond,bc_xl,bc_xr,bc_yb,bc_yt,dx,dy;
	sol->nx = nx;
	sol->ny = ny;
	sol->t = 0;
	sol->told = 0;
	sol->dt = 0.0001;
	N = nx * ny;
	mynx = (nx + sol->mpi_size - 1) / sol->mpi_size;
	istart = sol->myrank * mynx * ny;
	if (sol->myrank == sol->mpi_size-1) {
		mynx = nx % sol->mpi_size;
		iend = N;
	}
	myN = mynx*ny; 
	iend = istart + myN;
	sol->mynx = mynx;
	sol->myN = myN;
	// init local arrays
	size = sizeof(double) * (mynx+2)*(ny+2); // +2 for ghost points
	sol->u = (double *)malloc(size);
	sol->uold = (double *)malloc(size);
	sol->rhs = (double *)malloc(size);
	sol->x = (double *)malloc(size);
	sol->y = (double *)malloc(size);
	sol->IC = (double *)malloc(size);

	dx=(double)1/(nx+1);
	dy=(double)1/(ny+1); // +2 for BC

	int iglobal,ilocal;

	if (sol->myrank == 0){
		// BC
		printf("\n Provide boundary condition at x=0\n");
		fflush(stdout);
		scanf("%lf", &bc_xl);
		printf("\n Provide boundary condition at x=l\n");
		fflush(stdout);
		scanf("%lf", &bc_xr);
		printf("\n Provide boundary condition at y=0\n");
		fflush(stdout);
		scanf("%lf", &bc_yb);
		printf("\n Provide boundary condition at y=l\n");
		fflush(stdout);
		scanf("%lf", &bc_yt);
		printf("\n Provide initialization condition \n");
		fflush(stdout);
		scanf("%lf", &initcond);
	}

	if (sol->myrank == 0){
		printf("Broadcasting BC to all mpi ranks...\n");
		fflush(stdout);
	}
		//send to all processor
		MPI_Bcast(&bc_xl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&bc_xr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&bc_yb,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&bc_yt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&initcond,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	if (sol->myrank == 0){
		printf("done.\n");
		fflush(stdout);
	}

	// init grid, u, uold, rhs
	for(int i=0;i<myN; ++i){
		ilocal = i + ny + 3 + i / ny * 2; // this index exclude ghosh & boundary points
		iglobal = i + istart;
		ix = iglobal / ny;
		iy = iglobal % ny;
		sol->x[ilocal] = (ix+1) * dx;
		sol->y[ilocal] = (iy+1) * dy;
		sol->u[ilocal] = 0;
		sol->uold[ilocal] = 0;
		sol->rhs[ilocal] = 0;
		sol->IC[ilocal] = initcond;
	}

	if (sol->myrank == 0){
		//BC for left boundary xl
		for(int i=0;i<ny+2;++i)
		{sol->u[i] = bc_xl; } //first nx elements of the array
	}
	if (sol->myrank == sol->mpi_size-1){
		//BC for right boundary xr
		l_ar=myN-ny-2;
		for(int i=l_ar;i<myN;++i)
		{sol->u[i] = bc_xr; } //last nx elements of the array
	}
	//BC for bottom and top boundary
	for(int i=0;i<mynx;++i){
		//Every nx_th element represent first column element
		sol->u[i*(ny+2)] = bc_yb;
		//Every nx-1 element represent last column element
		sol->u[(i+1)*(ny+2)-1] = bc_yt; 
	}	


	/* comment out Metis stuff
	int xadj[myN+1],adjncy[4*myN],vtxdist[sol->mpi_size+1];
	double	xyz[2*myN];
	int idx_adjncy=0;

	// vertex number to processor
	for(int i=0; i<sol->mpi_size+1;++i){
		vtxdist[i] = (N + sol->mpi_size - 1) / sol->mpi_size * i;
	}
	vtxdist[i] = N;

	// prepare grid for parMetis partition
	for(int i=istart;i<iend; ++i){
		ix = i / ny;
		iy = i - ny * ix;
		myidx = i - istart; //local index 0 ~ myN-1
		//create CSR format grid for ParMetis
		xyz[2*myidx] = ix * dx;
		xyz[2*myidx+1] = iy * dy;
		xadj[myidx] = idx_adjncy;
		if((ix == 0 || ix == nx-1) && (iy == 0 || iy == ny-1)){
			// corner points
			adjncy[idx_adjncy] = (ix == 0) ? i+ny : i-ny;
			adjncy[idx_adjncy+1] = (iy == 0) ? i + 1 : i -1;
			idx_adjncy+=2;
		} else if (ix == 0 || ix == nx-1){
			// x-edge points
			adjncy[idx_adjncy] = (ix == 0) ? i+ny : i-ny;
			adjncy[idx_adjncy+1] = i-1;
			adjncy[idx_adjncy+2] = i+1;
			idx_adjncy+=3;
		} else if (iy == 0 || iy == ny-1){
			// y-edge points
			adjncy[idx_adjncy] = (iy == 0) ? i+1 : i-1;
			adjncy[idx_adjncy+1] = i-ny;
			adjncy[idx_adjncy+2] = i+ny;
			idx_adjncy+=3;
		} else {
			// interior points
			adjncy[idx_adjncy] = i-1;
			adjncy[idx_adjncy+1] = i+1;
			adjncy[idx_adjncy+2] = i-ny;
			adjncy[idx_adjncy+3] = i+ny;
			idx_adjncy+=4;
		}
	}
	// last point
	xadj[myidx] = idx_adjncy;

	int wgtflag=0, numflag=0, ncon=0, nparts=sol->mpi_size, ndims=2, ncon=1, part[nparts], edgecut, options[]= {0, 0, 0};
	double tpwgts=1/nparts, ubvec= 1.05;

	//parmetis function call
	ParMETIS_V3_PartGeomKway(vtxdist,xadj,adjncy,NULL,NULL,&wgtflag,&numflag,&ndims,xyz,&ncon,&nparts,&tpwgts,&ubvec,options,&edge,part,MPI_COMM_WORLD);

*/

	return 0;
}
