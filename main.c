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
	int tf,tf_old=0;

	ier = init(sol);
	
	//ier = test_sin_setIC(sol);

	while(sol->t < sol->total_t){

		if (sol->myrank == 0){
			printf("time step %f\n", sol->t);
			fflush(stdout);}

		ier = test_sin(sol);
		ier = timestepping(sol);
		tf=sol->t/sol->file_sv_t;
		if(tf!=tf_old)
		{
			ier = print_sol(sol);
			tf_old=tf;
		}
		
	}

	ier = print_sol(sol);
	MPI_Finalize();

	return 0;
}

