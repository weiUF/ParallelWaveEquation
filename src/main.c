//#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "functions.h"

/*
 main program 
 *
 */

int main(int argc, char *argv[]){

	struct data *sol;
	sol = (struct data*)malloc(sizeof(struct data));
	
	int ier;
	
	ier = init(sol);

	ier = test_sin_setIC(sol);

	while(sol->t < 2){

		ier = test_sin(sol);
		ier = timestepping(sol);
	}

	ier = print_sol(sol);
	return 0;
}

