#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "functions.h"

/*
 main program 
 *
 */

int main(int argc, char *argv[]){

	// initiate data structure
	struct data *sol;

	sol = (struct data*)malloc(sizeof(struct data));

	// Initialize the MPI environment
	MPI_Init(&argc, &argv);

	// Get the number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &sol->mpi_size);

	// Get the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &sol->myrank);


	int ier;

	ier = init(sol);

	ier = test_sin_setIC(sol);

	while(sol->t < 0.001){

		printf("rank[%2d] time %f\n",sol->myrank, sol->t);
		fflush(stdout);

		ier = test_sin(sol);
		ier = timestepping(sol);
		//ier = print_sol(sol);
	}

	ier = print_sol(sol);
	MPI_Finalize();

	return 0;
}

