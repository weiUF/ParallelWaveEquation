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

	//ier = test_sin_setIC(sol);

	while(sol->t < 1){
        
        	ier = print_sol(sol);
		ier = test_sin(sol);
		ier = timestepping(sol);
        
	}

	
	return 0;
}

