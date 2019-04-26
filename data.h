
struct data {
	int nx, ny, myrank, mpi_size,mynx,myN;
	double t, told, dt, total_t,file_sv_t;
	double *u, *uold, *rhs, *x, *y, *IC;
};


