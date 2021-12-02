#include "extra_functions.h"
#include "mv.h"

#define TO_STRING(s) #s /* convert variable to its name */

/* 
compute derivatives at each order independently, i.e. for data u take ux, uy, uxx, uyy, etc
and then combine into a large two dimensional matrix with columns [1 ux uy uxx uyy ...] and
rows corresponding to the time steps. this will be simpler for making polynomial products
later when constructing theta matrix

so something like:
	ux = FDM(u,x)
	uy = FDM(u,y)
	uxx = FDM(u,x)
	...
*/

/*
read csv of data into data file
need to get matrix dims


*/


int main(int argc, char** argv)
{
	int N_THREADS = 1; // argv[1]; /* threads for */
	/* initialize vectors and matrices */
	int N_t = 0;
	int M_t = 0;
	Vector t; /* time vector */
	get_matrix_dims("../../time.csv", &N_t, &M_t);
	t.length = N_t;
	allocate_Vector(&t);
	load_vector("../../time.csv", &t);
	double dt = t.vector[2] - t.vector[1];

	gsl_vector_view u_vec = gsl_vector_view_array(t.vector, t.length);
	//gsl_vector_fprintf(stdout, &u_vec.vector, "%lf");

	int N_x = 0;
	int M_x = 0;
	Vector x; /* spatial vector */
	get_matrix_dims("../../space.csv", &N_x, &M_x);
	x.length = N_x;
	allocate_Vector(&x);
	load_vector("../../space.csv", &x);
	double dx = x.vector[2] - x.vector[1];

	int N = 0; /* rows of u matrix / data matrix */
	int M = 0; /* columns of u matrix / data matrix */
	get_matrix_dims("../../u.csv", &N, &M);
	Matrix U; /* solution matrix */
	U.rows = N;
	U.cols = M;
	allocate_Matrix(&U);
	load_matrix("../../u.csv", &U);

	/* spatial and temporal derivatives of U matrix */
	Matrix U_t;
	U_t.rows = N;
	U_t.cols = M;
	allocate_Matrix(&U_t);
	FDM(&U, &U_t, dt, 1, 1, N_THREADS);

	Matrix U_x;
	U_x.rows = N;
	U_x.cols = M;
	allocate_Matrix(&U_x);
	FDM(&U, &U_x, dx, 0, 1, N_THREADS);
	
	Matrix U_xx;
	U_xx.rows = N;
	U_xx.cols = M;
	allocate_Matrix(&U_xx);
	FDM(&U, &U_xx, dx, 0, 2, N_THREADS);

	/* create N*M vectors of the above matrices for Theta */
	int L = N * M; 
	Vector u;
	u.length = L;
	allocate_Vector(&u);
	matrix_to_vector(&U, &u);

	Vector u_t;
	u_t.length = L;
	allocate_Vector(&u_t);
	matrix_to_vector(&U_t, &u_t);

	Vector u_x;
	u_x.length = L;
	allocate_Vector(&u_x);
	matrix_to_vector(&U_x, &u_x);

	Vector u_xx;
	u_xx.length = L;
	allocate_Vector(&u_xx);
	matrix_to_vector(&U_xx, &u_xx);

	/* matrix of derivatives */
	Matrix Derivatives;
	Derivatives.rows = L;
	Derivatives.cols = 2;
	allocate_Matrix(&Derivatives);

	/* fill derivatives */
	for (int i = 0; i < L; i++)
	{
		Derivatives.matrix[i][0] = u_x.vector[i];
		Derivatives.matrix[i][1] = u_xx.vector[i];
		//printf("d[%d][0]: %e\n", i, Derivatives.matrix[i][0]);
		//printf("d[%d][1]: %e\n", i, Derivatives.matrix[i][1]);
	}

	/* create Theta matrix */
	int p = 2; /* degree of polynomial in potential differential equation */
	int n_ders = Derivatives.cols; /* number of derivatives taken (here u_x, u_xx) */
	Matrix T;
	T.rows = L;
	T.cols = (p + 1) * (1 + n_ders); 	/* theta has P+1 initial columns of 1 and u,...,u^P followed by
										n_derivatives columns of linear derivative terms followed 
										by P*n_derivatives columns of nonlinear derivative terms */
	allocate_Matrix(&T);
	Theta(&T, &u, &Derivatives, p, N_THREADS);
	for (int i = 0; i < T.rows; i++)
	{
		for (int j = 0; j < T.cols; j++)
		{
			//printf("Theta[%d][%d] = %e\n", i, j, T.matrix[i][j]);
		}
	}
	
	Vector w;
	w.length = T.cols;
	w.vector = Train(&T, &u_t, pow(10.0,-5.0), 1.0, 10, 4);
	for (int i = 0; i < w.length; i++)
	{
		printf("outside w[%d] = %lf\n", i, w.vector[i]);
	}

	free(t.vector);
	free(x.vector);
	free(U.matrix);
	free(U_t.matrix);
	free(U_x.matrix);
	free(U_xx.matrix);
	free(u.vector);
	free(u_t.vector);
	free(u_x.vector);
	free(u_xx.vector);
	free(Derivatives.matrix);
	free(T.matrix);

	/*
	still need to finish Theta function and figure out
	how to implement sparse regression algorithm 
	*/

	return 0;
}