#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "extra_functions.h"
#include "sampling.h"


void FDM(Matrix* M, Matrix* M_d, double dx, int variable, int d, int N_THREADS)
{
	/* 
	second-order finite difference method.  supports up to d = 2.

	inputs:
		M: array for data from function of interest in ODE/PDE
		M_d: array for the derivative to fill up
		variable: the variable with respect to which the derivative will be taken
		dx: discretization of space
		d: order of derivative to take, i.e. d=1,2,3,...
	*/
	int thread_id;

	omp_set_num_threads(N_THREADS);

	int i, j;
	if (variable == 0)
	{
		if (d == 1)
		{
			//n_first_order += 1;
			#pragma omp parallel private(thread_id, i, j)
			thread_id = omp_get_thread_num();

			#pragma omp for
			for (i = 1; i < M_d->rows - 1; i++)
			{
				for (j = 0; j < M_d->cols; j++)
				{
					/* interior points */
					M_d->matrix[i][j] = (M->matrix[i + 1][j] - M->matrix[i - 1][j]) / (2 * dx);
					//printf("M_t[%d][%d] = %e\n", i, j, M_d->matrix[i][j]);
				}
			}

			/* boundary points */
			#pragma omp for
			for (j = 0; j < M_d->cols; j++)
			{
				M_d->matrix[0][j] = (-3.0 / 2 * M->matrix[0][j] + 2 * M->matrix[1][j] - M->matrix[2][j] / 2) / dx;
				M_d->matrix[M_d->rows - 1][j] = (3.0 / 2 * M->matrix[M->rows - 1][j] - 2 * M->matrix[M->rows - 2][j] + M->matrix[M->rows - 3][j] / 2) / dx;
			}
		}

		if (d == 2)
		{
			//n_sec_order += 1;
			//n_first_order += 1;
			#pragma omp parallel private(thread_id,i,j)
			thread_id = omp_get_thread_num();

			#pragma omp for
			for (i = 1; i < M_d->rows - 1; i++)
			{
				for (j = 0; j < M_d->cols; j++)
				{
					/* interior points */
					M_d->matrix[i][j] = (M->matrix[i + 1][j] - 2 * M->matrix[i][j] + M->matrix[i - 1][j]) / pow(dx, 2.0);
				}
			}
			/* boundary points */
			#pragma omp for
			for (j = 0; j < M_d->cols; j++)
			{
				M_d->matrix[0][j] = (2 * M->matrix[0][j] - 5 * M->matrix[1][j] + 4 * M->matrix[2][j] - M->matrix[3][j]) / pow(dx, 2.0);
				M_d->matrix[M_d->rows-1][j] = (2 * M->matrix[M->rows-1][j] - 5 * M->matrix[M->rows-2][j] + 4 * M->matrix[M->rows-3][j] - M->matrix[M->rows-4][j]) / pow(dx, 2.0);
			}
		}
	}


	if (variable == 1)
	{
		if (d == 1)
		{
			//n_first_order += 1;
			#pragma omp parallel private(thread_id,i,j)
			thread_id = omp_get_thread_num();
			
			#pragma omp for
			for (j = 1; j < M_d->cols - 1; j++)
			{
				for (i = 0; i < M_d->rows; i++)
				{
					/* interior points */
					M_d->matrix[i][j] = (M->matrix[i][j + 1] - M->matrix[i][j - 1]) / (2 * dx);
					//printf("M_t[%d][%d] = %e\n", i, j, M_d->matrix[i][j]);
				}
			}

			/* boundary points */
			#pragma omp for
			for (i = 0; i < M_d->rows; i++)
			{
				M_d->matrix[i][0] = (-3.0 / 2 * M->matrix[i][0] + 2 * M->matrix[i][1] - M->matrix[i][2] / 2) / dx;
				M_d->matrix[i][M_d->cols-1] = (3.0 / 2 * M->matrix[i][M_d->cols - 1] - 2 * M->matrix[i][M_d->cols - 2] + M->matrix[i][M_d->cols - 3] / 2) / dx;
			}
		}

		if (d == 2)
		{
			//n_sec_order += 1;
			#pragma omp parallel private(thread_id,i,j)
			thread_id = omp_get_thread_num();

			#pragma omp for
			for (j = 1; j < M_d->cols - 1; j++)
			{
				for (i = 0; i < M_d->rows; i++)
				{
					/* interior points */
					M_d->matrix[i][j] = (M->matrix[i][j + 1] - 2 * M->matrix[i][j] + M->matrix[i][j - 1]) / pow(dx, 2.0);
				}
			}

			/* boundary points */
			#pragma omp for
			for (i = 0; i < M_d->rows; i++)
			{
				M_d->matrix[i][0] = (2 * M->matrix[i][0] - 5 * M->matrix[i][1] + 4 * M->matrix[i][2] - M->matrix[i][3]) / pow(dx, 2.0);
				M_d->matrix[i][M_d->cols - 1] = (2 * M->matrix[i][M_d->cols - 1] - 5 * M->matrix[i][M_d->cols - 2] + 4 * M->matrix[i][M_d->cols - 3] - M->matrix[i][M_d->cols - 4]) / pow(dx, 2.0);
			}
		}
	}
}

void Theta(Matrix* T, Vector* u, Matrix* D, int P, int N_THREADS)
{
	/*
	creates the Theta matrix.  supports P = 2.

	inputs:
		T: matrix for theta, to be filled
		u: data from simulation
		D: derivatives of u
		P: maximum polynomial order for combinations of u,...,u^{P} and u derivatives
	*/

	int thread_id;

	omp_set_num_threads(N_THREADS);
	
	int i, j;

	if (P == 1)
	{
		#pragma omp parallel private(thread_id,i,j)
		thread_id = omp_get_thread_num();

		/* fill first two columns */
		#pragma omp for
		for (i = 0; i < T->rows; i++)
		{
			T->matrix[i][0] = 1.0;
			T->matrix[i][1] = u->vector[i];
		}

		/* fill rest of linear terms */
		#pragma omp for
		for (j = 0; j < D->cols; j++)
		{ /* loop through columns in derivative matrix */
			for (i = 0; i < D->rows; i++)
			{ /* loop through rows for each derivative */
				T->matrix[i][j + (P + 1)] = D->matrix[i][j];
			}
		}

		/* fill nonlinear terms (first order polynomials of u and derivatives) */
		/* start with column following last linear term */
		#pragma omp for
		for (j = 0; j < D->cols; j++)
		{ /* loop through columns in derivative matrix */
			for (i = 0; i < D->rows; i++)
			{ /* loop through rows for each derivative */
				T->matrix[i][j + D->cols + (P + 1)] = u->vector[i] * D->matrix[i][j];
			}
		}
	}

	if (P == 2)
	{
		#pragma omp parallel private(thread_id,i,j)
		thread_id = omp_get_thread_num();

		/* fill first two columns */
		#pragma omp for
		for (i = 0; i < T->rows; i++)
		{
			T->matrix[i][0] = 1.0;
			T->matrix[i][1] = u->vector[i];
			T->matrix[i][2] = u->vector[i] * u->vector[i];
		}

		/* fill rest of linear terms */
		#pragma omp for
		for (j = 0; j < D->cols; j++)
		{ /* loop through columns in derivative matrix */
			for (i = 0; i < D->rows; i++)
			{ /* loop through rows for each derivative */
				T->matrix[i][j + (P + 1)] = D->matrix[i][j];
			}
		}

		/* fill nonlinear terms (first order polynomials of u and derivatives) */
		/* start with column following last linear term */
		#pragma omp for
		for (j = 0; j < D->cols; j++)
		{ /* loop through columns in derivative matrix */
			for (i = 0; i < D->rows; i++)
			{ /* loop through rows for each derivative */
				T->matrix[i][j + D->cols + (P + 1)] = u->vector[i] * D->matrix[i][j];
			}
		}

		/* fill nonlinear terms (second order polynomials of u*u and derivatives) */
		#pragma omp for
		for (j = 0; j < D->cols; j++)
		{ /* loop through columns in derivative matrix */
			for (i = 0; i < D->rows; i++)
			{ /* loop through rows for each derivative */
				T->matrix[i][j + 2 * D->cols + (P + 1)] = pow(u->vector[i], 2.0) * D->matrix[i][j];
			}
		}
	}
}
/*
void Theta(int count, ...)
{
	/*
	function to construct Theta matrix.
	
	inputs:
		count: first parameter is the number of proceeding arguments
		...: the proceeding arguments are double* corresponding to u and its derivatives
			 with the second-to-last argument being a preallocated, empty Theta matrix 
			 and the last argument being an int value relating to the highest order 
			 of polynomial to be constructed

			 i.e. something like Theta(5, double* u, double* ux, double* uxx, double** Theta, 2) 
			 would return a matrix 
								Theta = [ 1 u ux uxx u*u u*ux u*uxx ]
			 where each of 1,...,u*uxx is a column vector

	*/
	/*
	va_list derivatives;
	double** element;
	va_start(derivatives, count);
	printf("count: %d\n", count);
	for (int d = 0; d < count-1; d++)
	{
		element = va_arg(derivatives, double** );
		int size = sizeof(element);
		printf("%d\n", size);
		printf("%e\n", element[0][0]);
	}
	va_end(derivatives);
}
*/


int count_nonzero_elements(Vector* V)
{
	/*
	count nonzero elements in a vector.  used to determine sparsity of solution
	*/

	int count = 0;

	for (int i = 0; i < V->length; i++)
	{
		if (V->vector[i] == 0.0)
		{
			count++;
		}
	}

	return count;
}

void STR(Matrix* T, Vector* ut, gsl_vector* w, double lambda, double tol, int iterations)
{
	/*************** normalize the data ***************/
	Vector theta_vec;
	theta_vec.length = T->rows * T->cols;
	allocate_Vector(&theta_vec);
	matrix_to_vector(T, &theta_vec);

	/* convert theta to gsl matrix to get norm */
	gsl_matrix_view gsl_theta = gsl_matrix_view_array(theta_vec.vector, T->rows, T->cols);
	gsl_matrix_view gsl_theta_norm = gsl_matrix_view_array(theta_vec.vector, T->rows, T->cols);
	for (int j = 0; j < T->cols; j++)
	{
		gsl_vector_view theta_col = gsl_matrix_column(&gsl_theta_norm, j);
		double scale = 1.0 / gsl_blas_dnrm2(&theta_col);
		for (int i = 0; i < T->rows; i++)
		{
			double new_input = gsl_matrix_get(&gsl_theta, i, j) * scale;
			gsl_matrix_set(&gsl_theta_norm, i, j, new_input);
		}
	}

	/*************** get initial ridge regression result ***************/

	/* algorithm has the form
				w = (X.X + lamda*lambda)^{-1} * X.y
	   where w is the solution, X is the Theta matrix, and y is u_t.
	   lambda is a regularization factor, typically proportional to the
	   condition number of X
	*/

	/* compute matrix product X.X */
	gsl_matrix* x_dot_x = gsl_matrix_alloc(T->cols, T->cols); /* [NxM][MxN] = [NxN] */
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &gsl_theta_norm, &gsl_theta_norm, 0.0, x_dot_x);
	
	/* add factor of lambda to diagonal */
	for (int d = 0; d < T->cols; d++)
	{
		double new_value = gsl_matrix_get(x_dot_x, d, d) + lambda;
		gsl_matrix_set(x_dot_x, d, d, new_value);
	}

	/* compute X.y */
	gsl_vector_view gsl_ut = gsl_vector_view_array(ut->vector, ut->length);
	gsl_vector* x_dot_ut = gsl_vector_alloc(T->cols);
	gsl_blas_dgemv(CblasTrans, 1.0, &gsl_theta_norm, &gsl_ut, 0, x_dot_ut);

	gsl_matrix* X_cov = gsl_matrix_alloc(T->rows, T->cols);
	double chi_sq;
	gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(T->rows, T->cols);
	gsl_multifit_linear(&gsl_theta_norm, &gsl_ut, w, X_cov, &chi_sq, work);

	int num_terms = T->cols;
	Vector indices;
	indices.length = T->cols;
	allocate_Vector(&indices);

	for (int it = 0; it < iterations; it++)
	{
		for (int i = 0; i < w->size; i++)
		{
			if (gsl_vector_get(w, i) < tol)
			{
				indices.vector[i] = 0.0;
				gsl_vector_set(w, i, 0.0);
			}
			else
			{
				indices.vector[i] = 1.0;
			}
		}

		gsl_vector_view gsl_indices = gsl_vector_view_array(&indices.vector, &indices.length);
		if (gsl_vector_isnull(&gsl_indices))
		{
			fprintf(stderr, "Tolerance is too high\n");
			exit(2);
		}

		if (count_nonzero_elements(&indices) == num_terms)
		{
			fprintf(stderr, "Tolerance too low\n");
			exit(3);
		}

		/* optimized version of theta after a given iteration */
		int new_cols = indices.length - count_nonzero_elements(&indices);
		printf("new cols %d\n", new_cols);
		gsl_matrix* gsl_X_opt = gsl_matrix_alloc(T->rows, new_cols);
		for (int i = 0; i < gsl_X_opt->size1; i++)
		{
			for (int j = 0; j < gsl_X_opt->size2; j++)
			{
				if (indices.vector[j] != 0.0)
				{
					double val = gsl_matrix_get(&gsl_theta_norm, i, j);
					gsl_matrix_set(gsl_X_opt, i, j, val);
				}
			}
		}

		/* compute new matrix product X_opt.X_opt */
		gsl_matrix* x_dot_x_opt = gsl_matrix_alloc(new_cols, new_cols); /* [NxM][MxN] = [NxN] */
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, gsl_X_opt, gsl_X_opt, 0.0, x_dot_x_opt);

		/* add factor of lambda to diagonal */
		for (int d = 0; d < new_cols; d++)
		{
			double new_value = gsl_matrix_get(x_dot_x_opt, d, d) + lambda;
			gsl_matrix_set(x_dot_x, d, d, new_value);
		}

		/* compute X.y */
		gsl_blas_dgemv(CblasTrans, 1.0, &gsl_theta, &gsl_ut, 0, x_dot_ut);

		gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(T->rows, T->cols);
		gsl_vector* w_temp = gsl_vector_alloc(new_cols);
		gsl_multifit_linear(gsl_X_opt, &gsl_ut, w_temp, X_cov, &chi_sq, work);

		for (int i = 0; i < w_temp->size; i++)
		{
			if (gsl_vector_get(w, indices.vector[i]) != 0.0)
			{
				gsl_vector_set(w, i, gsl_vector_get(w_temp, i));
			}
		}

		//gsl_matrix_free(gsl_X_opt);
		//gsl_vector_free(x_dot_x_opt);
	}



	gsl_matrix_free(x_dot_x);
	gsl_vector_free(x_dot_ut);
	gsl_matrix_free(X_cov);
	
}

double* Train(Matrix* T, Vector* ut, double lambda, double i_tol, int iterations, int N_THREADS)
{
	/*
	trains ridge regression algorithm to find governing equation from data by
	eliminating potential terms included in Theta matrix

	inputs:
		T: Theta matrix
		ut: time derivative of data
		lambda: l0 scaling factor in ridge regression
		i_tol: incremental tolerance, tolerance will be incremented
			   by this amount on each iteration. this is used to
			   eliminate potential terms in the PDE
		iterations: number of iterations for training
	*/

	double L0 = 0.0005;

	/* create training and test arrays and vectors from Theta and u_t*/
	Matrix X_train;
	X_train.rows = (int)(0.8 * T->rows); /* 80% of data is for training */
	X_train.cols = T->cols;
	allocate_Matrix(&X_train);

	Matrix X_test;
	X_test.rows = T->rows - X_train.rows; /* 80% of data is for training */
	X_test.cols = T->cols;
	allocate_Matrix(&X_test);

	Vector x_train;
	x_train.length = X_train.rows * X_train.cols;
	allocate_Vector(&x_train);

	Vector x_test;
	x_test.length = X_test.rows * X_test.cols;
	allocate_Vector(&x_test);

	Vector y_train;
	y_train.length = (int)(0.8 * T->rows);
	allocate_Vector(&y_train);

	Vector y_test;
	y_test.length = T->rows - y_train.length;
	allocate_Vector(&y_test);

	Vector w; /* this is the sparse solution vector */
	w.length = X_test.cols;
	allocate_Vector(&w);

	/* fill training and test arrays with random points from Theta 
	see sampling.c file for more info */
	/* don't double reference T and ut */
	sample_data(&X_train, &X_test, &y_train, &y_test, T, ut, N_THREADS);

	/* perform least squares on training sets */
	matrix_to_vector(&X_train, &x_train);
	matrix_to_vector(&X_test, &x_test);
	gsl_matrix_view gsl_X_train = gsl_matrix_view_array(x_train.vector, X_train.rows, X_train.cols);
	gsl_vector_view gsl_y_train = gsl_vector_view_array(y_train.vector, y_train.length);
	gsl_matrix_view gsl_X_test = gsl_matrix_view_array(x_test.vector, X_test.rows, X_test.cols);
	gsl_vector_view gsl_y_test1 = gsl_vector_view_array(y_test.vector, y_test.length);
	gsl_vector_view gsl_y_test2 = gsl_vector_view_array(y_test.vector, y_test.length);
	gsl_vector_view gsl_w = gsl_vector_view_array(w.vector, w.length);
	
	gsl_matrix* X_cov = gsl_matrix_alloc(X_train.rows, X_train.cols);
	double chi_sq; 
	gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(X_train.rows, X_train.cols);
	gsl_multifit_linear(&gsl_X_train, &gsl_y_train, &gsl_w, X_cov, &chi_sq, work);

	/* compute dot product of x_test and w for ridge regression */
	gsl_vector* xtest_dot_w = gsl_vector_alloc(X_test.rows);
	gsl_blas_dgemv(CblasNoTrans, 1.0, &gsl_X_test, &gsl_w, 0.0, xtest_dot_w);

	/* 
	error norm is L2 norm of y_test - xtest_dot_w 
	
	this goes into the computed error given by
			||y_test - xtest_dot_w||_{2} + L0*nonzero(w)
	where nonzero(w) is the number of nonzero elements of w
	*/
	gsl_vector_sub(&gsl_y_test1, xtest_dot_w); /* subtraction is stored in gsl_y_test */
	double L2_norm = gsl_blas_dnrm2(&gsl_y_test1);
	int nze = count_nonzero_elements(&w);

	double error = L2_norm + L0 * nze;

	double tol = 0.0;

	for (int i = 0; i < w.length; i++)
	{
		printf("inside before w[%d] = %lf\n", i, w.vector[i]);
	}
	
	/* optimize the error with ridge regression */
	for (int it = 0; it < iterations; it++)
	{
		STR(T, ut, &gsl_w, lambda, tol, iterations);
		
		gsl_vector* xtest_dot_w = gsl_vector_alloc(X_test.rows);
		gsl_blas_dgemv(CblasNoTrans, 1.0, &gsl_X_test, &gsl_w, 0.0, xtest_dot_w);

		/*
		error norm is L2 norm of y_test - xtest_dot_w

		this goes into the computed error given by
				||y_test - xtest_dot_w||_{2} + L0*nonzero(w)
		where nonzero(w) is the number of nonzero elements of w
		*/
		gsl_vector_sub(&gsl_y_test2, xtest_dot_w); /* subtraction is stored in gsl_y_test */
		double L2_norm_new = gsl_blas_dnrm2(&gsl_y_test2);
		int nze_new = count_nonzero_elements(&w);

		double error_new = L2_norm + L0 * nze;

		if (error_new < error)
		{
			error = error_new;
			tol += i_tol;
		}
	}

	for (int i = 0; i < w.length; i++)
	{
		printf("inside after w[%d] = %lf\n", i, w.vector[i]);
	}

	gsl_matrix_free(X_cov);
	gsl_vector_free(xtest_dot_w);

	return w.vector;//&gsl_w.vector;
}


