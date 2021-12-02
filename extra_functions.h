#ifndef EXTRAS
#define EXTRAS

#include <stdarg.h>
#include "matrix_dims.h"
#include "mv.h"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>

#include <omp.h>

//int n_first_order = 0; // count for number of times a first order derivative is taken
//int n_sec_order = 0; // count for number of times a second order derivative is taken

void FDM(Matrix* M, Matrix* M_d, double dx, int variable, int d, int N_THREADS);
void Theta(Matrix* T, Vector* u, Matrix* D, int P, int N_THREADS);
//void Theta(int count, ...);

int count_nonzero_elements(Vector* V);
void STR(Matrix* T, Vector* ut, gsl_vector* w, double lamda, double tol, int iterations);
double* Train(Matrix* T, Vector* ut, double lambda, double i_tol, int iterations, int N_THREADS);


#endif // !EXTRAS

