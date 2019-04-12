
struct data {
	int nx, ny, myrank, mpi_size,mynx,myN;
	double t, told, dt;
	double *u, *uold, *rhs, *x, *y, *IC;
};


