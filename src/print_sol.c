#include <stdio.h>
#include "data.h"
#include <stdlib.h> 

int print_sol(struct data *s){
	
	int idx,i,j;
    double time = s->t;
    char filename[300];
    
    sprintf(filename, "wave_ss_%5f.csv", time);

    FILE *out_file = fopen(filename, "w"); // write only 

    
    if (out_file == NULL || out_file == NULL) 
            {   
              printf("Error! Could not open file\n"); 
              exit(-1); // must include stdlib.h 
            }
            
    /*        
	printf("%s\n","dimension");
	printf("nx: %d ", s->nx);
	printf("ny: %d\n", s->ny);

	printf("%s\n","solution u:");
	
	*/

  //file header 
	fprintf(out_file,"X,Y,U\n");	

	for (i=0;i<s->nx;++i){
		for (j=0;j<s->ny;++j){
			idx = i * s->nx + j;
            fprintf(out_file, "%6f,%6f,%6f,\n",s->x[idx],s->y[idx],s->u[idx]); // write to file
		}
	}
    
    //time = time + 0.001;
    
	return 0;
}
