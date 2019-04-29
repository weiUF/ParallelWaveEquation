#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "data.h"

/*
	initialize data structure and create grid points
 */

int init(struct data *sol){
	
	// init	
	int nx=10, ny=10;//nx=10, ny=10; // no. of internal points
	int N, mynx, istart, iend, myN, size, ix, iy, l_ar,myN_modif; 
	double initcond,bc_xl,bc_xr,bc_yb,bc_yt,dx,dy,x_xl,y_yl;
	double dt,total_t,file_sv_t;
	int iglobal,ilocal;

	if (sol->myrank == 0){
		// BC
		printf("\n length in x direction \n");
		fflush(stdout);
		scanf("%lf", &x_xl);
		printf("\n Number of grid points in x direction, including boundary points \n");
		fflush(stdout);
		scanf("%d", &nx);
		printf("\n length in y direction\n");
		fflush(stdout);
		scanf("%lf", &y_yl);
		printf("\n Number of grid points in x direction, including boundary points \n");
		fflush(stdout);
		scanf("%d", &ny);
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
		printf("\n Provide Time Step Size \n");
		fflush(stdout);
		scanf("%lf", &dt);
		printf("\n Total duration of simulation in sec \n");
		fflush(stdout);
		scanf("%lf", &total_t);
		printf("\n Duration of output file writing in sec \n");
		fflush(stdout);
		scanf("%lf", &file_sv_t);
	}

	if (sol->myrank == 0){
		printf("Broadcasting BC to all mpi ranks...\n");
		fflush(stdout);
	}
		//send to all processor
		MPI_Bcast(&x_xl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&y_yl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&ny,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&bc_xl,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&bc_xr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&bc_yb,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&bc_yt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&initcond,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&dt,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&total_t,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&file_sv_t,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	sol->nx = nx;
	sol->ny = ny;
	sol->t = 0;
	sol->told = 0;
	sol->dt = dt;
	sol->total_t = total_t;
	sol->file_sv_t = file_sv_t;
	sol->x_xl=x_xl;
	sol->y_yl=y_yl;
	N = nx * ny;
	// commented mynx = (nx + sol->mpi_size - 1) / sol->mpi_size;
	mynx = nx / sol->mpi_size;
	istart = sol->myrank * mynx * ny;
	if (sol->myrank == sol->mpi_size-1) {
		mynx = nx - mynx * (sol->mpi_size - 1);
		iend = N;
	}
	//debug 
	//printf("rank[%d], mynx:%d\n",sol->myrank,mynx);
	//fflush(stdout);
	
	myN = mynx*ny; 
	iend = istart + myN;
	sol->mynx = mynx;
	sol->myN = myN;
	// init local arrays
	size = sizeof(double) * (mynx+2)*(ny+2); // +2 for ghost points
	// recommented size = sizeof(double) * (mynx+2)*(ny);
	sol->u = (double *)malloc(size);
	sol->uold = (double *)malloc(size);
	sol->rhs = (double *)malloc(size);
	sol->x = (double *)malloc(size);
	sol->y = (double *)malloc(size);
	sol->IC = (double *)malloc(size);




	if (sol->myrank == 0){
		printf("done.\n");
		fflush(stdout);
	}

	dx=(double)x_xl/(nx-1);
	dy=(double)y_yl/(ny-1); // as dx, dy matters for BC and internal grids 

	if (sol->myrank == 0)
	{ //implemented to avoid empty loop for one row
		// init grid, u, uold, rhs
		myN_modif=(myN/(mynx*ny))*(ny-2)*(mynx-1);
		//: removing two additional elements each column from loop, belonging to BC. requested variables are calcualted in every node so no additional steps required
		//:divide by mynx and mult by mynx-1, as first row (after ghost cells) is boundary condition. I
		for(int i=0;i<myN_modif; ++i){
			// commented ilocal = i + ny + 3 + i / ny * 2; // this index exclude ghosh & boundary points
			ilocal = i + 2*ny + 6 + i / (ny-2) * 4; //modified to skip 2nd and 2nd last element of each column -> these are BC nodes
			//iglobal = i + istart;
			ix = i / (ny-2)+istart/ny;
			iy = ilocal % (ny+2);//as ilocal is replicating motion in ny+2 column size matrix
			// commented sol->x[ilocal] = (ix+1) * dx;
			// commented sol->y[ilocal] = (iy+1) * dy;
			sol->x[ilocal] = (ix+1) * dx; //as we need to initialize first row with x=0, +1 added as first row is already removed from calculation
			sol->y[ilocal] = (iy-1) * dy; //bottom boundary (2nd column) should have y=0, which is not performed in this loop, but will affect for later columns 
			sol->u[ilocal] = 0;
			sol->uold[ilocal] = 0;
			sol->rhs[ilocal] = 0;
			sol->IC[ilocal] = initcond;
		}
		//sol->u[550] = 5;


	}
	else if(sol->myrank == sol->mpi_size-1)
	{ //implemented to avoid empty loop for one row
		// init grid, u, uold, rhs
		myN_modif=(myN/(mynx*ny))*(ny-2)*(mynx-1);
		//: removing two additional elements each column from loop, belonging to BC. requested variables are calcualted in every node so no additional steps required
		//:divide by mynx and mult by mynx-1, as last row (before ghost cells) is boundary condition. I
		for(int i=0;i<myN_modif; ++i){
			// commented ilocal = i + ny + 3 + i / ny * 2; // this index exclude ghosh & boundary points
			ilocal = i + ny + 4 + i / (ny-2) * 4; //modified to skip 2nd and 2nd last element of each column -> these are BC nodes
			//iglobal = i + istart;
			ix = i / (ny-2) + istart/ny;
			iy = ilocal % (ny+2);
			// commented sol->x[ilocal] = (ix+1) * dx;
			// commented sol->y[ilocal] = (iy+1) * dy;
			sol->x[ilocal] = (ix) * dx; //as we need to initialize first row with x=0, +1 added as first row is already removed from calculation
			sol->y[ilocal] = (iy-1) * dy; //as we need to initialize second coulmn with y=0
			sol->u[ilocal] = 0;
			sol->uold[ilocal] = 0;
			sol->rhs[ilocal] = 0;
			sol->IC[ilocal] = initcond;
		}
		//sol->u[550] = 5;


	}
	else
	{// init grid, u, uold, rhs
		myN_modif=(myN/ny)*(ny-2);//: removing two additional elements each column from loop, belonging to BC. requested variables are calcualted in every node so no additional steps
		for(int i=0;i<myN_modif; ++i){
			// commented ilocal = i + ny + 3 + i / ny * 2; // this index exclude ghosh & boundary points
			ilocal = i + ny + 4 + i / (ny-2) * 4; //modified to skip 2nd and 2nd last element of each column -> these are BC nodes
			//iglobal = i + istart;
			ix = i / (ny-2) + istart/ny;
			iy = ilocal % (ny+2);
			// commented sol->x[ilocal] = (ix+1) * dx;
			// commented sol->y[ilocal] = (iy+1) * dy;
			sol->x[ilocal] = (ix) * dx; //as we need to initialize first row with x=0
			sol->y[ilocal] = (iy-1) * dy; //as we need to initialize second coulmn with y=0
			sol->u[ilocal] = 0;
			sol->uold[ilocal] = 0;
			sol->rhs[ilocal] = 0;
			sol->IC[ilocal] = initcond;
		}
		//sol->u[550] = 5;

	}

	if (sol->myrank == 0){
		//BC for left boundary xl
		for(int i=ny+3;i<2*ny+3;++i)//
		{sol->u[i] = bc_xl; 
		sol->x[i] = (i/(ny+2)-1) * dx; //boundary x -values	
		sol->y[i] = (i%(ny+2)-1) * dy; //boundary y -values		
		} //first nx elements of the array
	}
	if (sol->myrank == sol->mpi_size-1){
		//BC for right boundary xr
		l_ar=myN-ny;
		for(int i=l_ar;i<myN;++i)
		{
			ilocal=i + ny + 3 + i / ny * 2;
			sol->u[ilocal] = bc_xr; 
			sol->x[ilocal] = ((istart+i)/ny) * dx;//boundary x -values
			sol->y[ilocal] = ((istart+i)%ny) * dy; //boundary y -values				
		} //last nx elements of the array
	}
	//BC for bottom and top boundary
	for(int i=1;i<mynx+1;++i){
		//Every nx_th element represent first column element
		sol->u[i*(ny+2)+1] = bc_yb;
		sol->x[i*(ny+2)+1] = ((i*(ny+2)+1)/(ny+2)+(istart/ny)-1) * dx; //boundary x -values	
		sol->y[i*(ny+2)+1] = ((i*(ny+2)+1)%(ny+2)-1) * dy; //boundary y -values
		//Every nx-1 element represent last column element
		sol->u[(i+1)*(ny+2)-2] = bc_yt;
		sol->x[(i+1)*(ny+2)-2] = (((i+1)*(ny+2)-2)/(ny+2)+(istart/ny)-1) * dx;//boundary x -values
		sol->y[(i+1)*(ny+2)-2] = (((i+1)*(ny+2)-2)%(ny+2)-1) * dy; //boundary y -values
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
