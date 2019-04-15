#include <mpi.h>
#include <stdio.h>
#include "data.h"

int print_sol(struct data *s){
	
	int idx,i,j,rank;
	
	for (rank=0; rank<s->mpi_size; ++rank)
	{ 
		if(s->myrank == rank){
			printf("rank[%d]\n", rank);
			printf("nx: %d ", s->mynx);
			printf("ny: %d\n", s->ny);
			printf("current time: %f\n", s->t);

			printf("%s\n","solution u:");
			for (i=0;i<s->mynx+2;++i){
				for (j=0;j<s->ny+2;++j){
					idx = i * (s->ny+2) + j;
					printf("%6f ", s->u[idx]);
				}
				printf("%s\n","");
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	return 0;
}
