#include <stdio.h>
#include "data.h"

int print_sol(struct data *s){
	
	int idx,i,j;

	printf("dimension\n");
	printf("nx: %d ", s->mynx);
	printf("ny: %d\n", s->ny);
	printf("current time: %f\n", s->t);

	printf("%s\n","solution u:");
	for (i=0;i<s->nx+2;++i){
		for (j=0;j<s->ny+2;++j){
			idx = i * (s->ny+2) + j;
			printf("%6f ", s->u[idx]);
		}
		printf("%s\n","");
	}

	return 0;
}
