#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "data.h"

int print_sol(struct data *s){
	
	int idx,i,j,rank,size,jj,ji,k;
	int rank_pass;
	int mynx, istart, iend, myN,ny,rec_nr;
	int temp_pr,temp_sc,tem_tr;
	double *prm_u,*prm_x,*prm_y;
	double *sec_u,*sec_x,*sec_y;
	double *ter_u,*ter_x,*ter_y;
	double *u_ng,*x_ng,*y_ng; //arrays to store non-ghost grids in every node
	
	FILE * fp, *fp1;
	MPI_Status status;
	char fn[30];
	//sprintf(fn, "%s_%f.csv", "wave_postprocessing",(s->file_sv_t));
	sprintf(fn, "wave_pp_%5f.csv",s->t);

	fp = fopen ("run_output.out", "a");
	//fp1 = fopen ("wave_postprocessing.csv", "w");
	fp1 = fopen (fn, "w");
	
	if(s->myrank == 0)
	{

		ny=s->ny;
		myN=s->myN;				
		istart = s->myrank * s->mynx * ny;
		size = sizeof(double) * (s->myN);
		u_ng = (double *)malloc(size); 
		x_ng = (double *)malloc(size); 
		y_ng = (double *)malloc(size);
		k=0;
		fprintf(fp1,"X:,Y:,U:\n");
		for (i=1;i<=s->mynx;++i) //skipping first and last rows, as they are ghost cells
		{
			for (j=0;j<ny+2;++j) //no direct way to skip first and last coulumn, every row, which represent ghost cells
			{ 
				ji=j%(ny+1);//0 for last column element
				jj=j%(ny+2);//0 for first column element
				if(ji!=0&&jj!=0)
				{
					idx = i * (ny+2) + j;
					x_ng[k]=s->x[idx];							
					y_ng[k]=s->y[idx];							
					u_ng[k]=s->u[idx];
					fprintf(fp1,"%6f,%6f,%6f \n", x_ng[k],y_ng[k],u_ng[k]);
					k++;
				}
				
			}
		}

		//MPI_Barrier(MPI_COMM_WORLD); 
			
	}

/*	for (rank=0; rank<s->mpi_size; ++rank)
	{ 
		if(s->myrank == rank){
			fprintf(fp,"rank[%d]\n", rank);
			fprintf(fp,"nx: %d ", s->mynx);
			fprintf(fp,"ny: %d\n", s->ny);
			fprintf(fp,"current time: %f\n", s->t);

			fprintf(fp,"%s\n","solution u:");
			for (i=0;i<s->mynx+2;++i){
				for (j=0;j<s->ny+2;++j){
					idx = i * (s->ny+2) + j;
					fprintf(fp,"%6f ", s->u[idx]);
				}
				fprintf(fp,"%s\n","");
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}*/
	for (rank=0; rank<(s->mpi_size-1); ++rank)
	{ 
		rank_pass=rank+1;
		if(s->myrank == rank_pass)
			{
				ny=s->ny;
				myN=s->myN;				
				istart = s->myrank * s->mynx * ny;
				size = sizeof(double) * (s->myN);
				u_ng = (double *)malloc(size); 
				x_ng = (double *)malloc(size); 
				y_ng = (double *)malloc(size);
				k=0;
				for (i=1;i<=s->mynx;++i) //skipping first and last rows, as they are ghost cells
				{
					for (j=0;j<ny+2;++j) //no direct way to skip first and last coulumn, every row, which represent ghost cells
					{ 
						ji=j%(ny+1);//0 for last column element
						jj=j%(ny+2);//0 for first column element
						if(ji!=0&&jj!=0)
						{
							idx = i * (ny+2) + j;
							x_ng[k]=s->x[idx];							
							y_ng[k]=s->y[idx];							
							u_ng[k]=s->u[idx];
							k++;
						}
						
					}
				}
				MPI_Send(&myN, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
				printf("Number send is: %d\n", myN);
				MPI_Send(&x_ng[0], myN, MPI_DOUBLE, 0, 9, MPI_COMM_WORLD);
				MPI_Send(&y_ng[0], myN, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD);
				MPI_Send(&u_ng[0], myN, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD);
				fprintf(fp,"%s\n","solution u internal:");
				/*for (i=0;i<s->mynx;++i){
					for (j=0;j<s->ny;++j){
						idx = i * (s->ny) + j;
						fprintf(fp,"%6f ", u_ng[idx]);
					}
					fprintf(fp,"%s\n","");
				}*/



	
			}		

		if(s->myrank == 0)
		{
			MPI_Recv(&temp_pr, 1, MPI_INT, MPI_ANY_SOURCE, 10, MPI_COMM_WORLD, &status);
			printf("Number received is: %d\n", temp_pr);
			rec_nr=temp_pr/(s->ny); //recovering number of rows
			size = sizeof(double) * temp_pr;
			prm_x = (double *)malloc(size);
			prm_y = (double *)malloc(size);
			prm_u = (double *)malloc(size);
			MPI_Recv(&prm_x[0], temp_pr, MPI_DOUBLE, rank_pass, 9, MPI_COMM_WORLD, &status);
			MPI_Recv(&prm_y[0], temp_pr, MPI_DOUBLE, rank_pass, 8, MPI_COMM_WORLD, &status);
			MPI_Recv(&prm_u[0], temp_pr, MPI_DOUBLE, rank_pass, 7, MPI_COMM_WORLD, &status); 
			/*fprintf(fp,"rank[%d]\n", rank);
			fprintf(fp,"nx: %d ", rec_nr);
			fprintf(fp,"ny: %d\n", s->ny);
			fprintf(fp,"current time: %f\n", s->t);

			fprintf(fp,"%s\n","solution received by rank=0 internal:");*/
			for (i=0;i<rec_nr;++i){
				for (j=0;j<s->ny;++j){
					idx = i * (s->ny) + j;
					//fprintf(fp,"%6f ", prm_u[idx]);
					fprintf(fp1,"%6f,%6f,%6f \n", prm_x[idx],prm_y[idx],prm_u[idx]);
				}
				//fprintf(fp,"%s\n","");
				//fprintf(fp1,"%s\n","");
			}
			free(prm_x);
			free(prm_y);
			free(prm_u);

				
		}
	
		
		MPI_Barrier(MPI_COMM_WORLD);
	}
	/*for (rank=1; rank<s->mpi_size; ++rank)
	{ 
		if(s->myrank == rank){
			fprintf(fp,"rank[%d]\n", rank);
			fprintf(fp,"nx: %d ", s->mynx);
			fprintf(fp,"ny: %d\n", s->ny);
			fprintf(fp,"current time: %f\n", s->t);

			fprintf(fp,"%s\n","solution u internal:");
			for (i=0;i<s->mynx;++i){
				for (j=0;j<s->ny;++j){
					idx = i * (s->ny) + j;
					fprintf(fp,"%6f ", u_ng[idx]);
				}
				fprintf(fp,"%s\n","");
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}*/
	/*if(s->myrank == 0) 
	{ 
		
		fprintf(fp,"rank[%d]\n", rank);
		fprintf(fp,"nx: %d ", s->mynx);
		fprintf(fp,"ny: %d\n", s->ny);
		fprintf(fp,"current time: %f\n", s->t);

		fprintf(fp,"%s\n","solution transferred to rank=0 internal:");
		for (i=0;i<s->mynx;++i){
			for (j=0;j<s->ny;++j){
				idx = i * (s->ny) + j;
				fprintf(fp,"%6f ", prm_u[idx]);
			}
			fprintf(fp,"%s\n","");
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
	}*/
	//fprintf(fp1,"current time: %f\n", s->t);

	fclose(fp);
	fclose(fp1);
	free(u_ng);
	free(x_ng);
	free(y_ng);
/*	if(s->myrank == 0)
	{
		free(prm_x);
		free(prm_y);
		free(prm_u);
	}*/
	//free(fn);
	return 0;
}
