#include <mpi.h>
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
	
	int nx,ny,N,i,size,bc_ilocal,myN_modif,mynx,myN,ilocal,l_ar;
	double t,dt;
	double *unew;
	nx = s->nx;
	ny = s->ny;
	N = (s->mynx+2)*(ny+2);
	t = s->t;
	dt = s->dt;
	size = sizeof(double)*N;
	unew = (double*)malloc(size);
	myN = s->myN;
	mynx=s->mynx;
	
	/*if (s->myrank == 0)
	{ 
		myN_modif=(myN/(mynx*ny))*(ny-2)*(mynx-1);		
		for(i=0;i<myN_modif; ++i){
			ilocal = i + 2*ny + 6 + i / (ny-2) * 4; //modified to skip 2nd and 2nd last element of each column -> these are BC nodes
			bc_ilocal=ilocal%(ny+2);// to exclude bc point and outer-ghost cell from iteration, i.e. two cells every column
			if(bc_ilocal>=2 && bc_ilocal<ny) {
				unew[ilocal] = s->rhs[ilocal];
				
			}

		}

	
	}
	else if(s->myrank == s->mpi_size-1)
	{
		myN_modif=(myN/(mynx*ny))*(ny-2)*(mynx-1);		
		//calculate u, rhs
		for(i=0;i<myN_modif; ++i){
			ilocal = i + ny + 4 + i / (ny-2) * 4; //modified to skip 2nd and 2nd last element of each column -> these are BC nodes
			bc_ilocal=ilocal%(ny+2);// to exclude bc point and outer-ghost cell from iteration, i.e. two cells every column
			if(bc_ilocal>=2 && bc_ilocal<ny) {
				unew[ilocal] = s->rhs[ilocal];
			}

		}

	}
	else
	{
		//calculate u, rhs
		for(i=0;i<myN; ++i){
			ilocal = i + ny + 3 + i / ny * 2; // this index skips ghosh & boundary points
			bc_ilocal=ilocal%(ny+2);// to exclude bc point and outer-ghost cell from iteration, i.e. two cells every column
			if(bc_ilocal>=2 && bc_ilocal<ny) {
				unew[ilocal] = s->rhs[ilocal];
			}

		}
	}*/

	//if (t > 0){
		for(i=0;i<N;++i){
			unew[i] = s->rhs[i];
		}
	//}else{
	//	//at the 1st step (need to use u' initial condtion)
	//	for(i=0;i<N;++i){
	//		unew[i] = s->IC[i];//0.5*(dt * dt * s->rhs[i] + 2 * s->u[i] + 2*dt*s->IC[i]); //IC=ut(t=0)
	//	}
	//}

	if (s->myrank == 0){
		//BC for left boundary xl
		for(int i=ny+3;i<2*ny+3;++i)//
		{unew[i]=s->u[i]; } //first nx elements of the array
	}
	if (s->myrank == s->mpi_size-1){
		//BC for right boundary xr
		l_ar=myN-ny;
		for(int i=l_ar;i<myN;++i)
		{
			ilocal=i + ny + 3 + i / ny * 2;
			unew[ilocal]=s->u[ilocal];}			
			//sol->u[ilocal] = bc_xr;  //last nx elements of the array
	}
	//BC for bottom and top boundary
	for(int i=1;i<mynx+1;++i){
		//Every nx_th element represent first column element
		unew[i*(ny+2)+1] = s->u[i*(ny+2)+1] ;
		//Every nx-1 element represent last column element
		unew[(i+1)*(ny+2)-2] = s->u[(i+1)*(ny+2)-2]; 
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
