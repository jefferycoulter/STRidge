#ifndef MV
#define MV

/*
* header file with definition of matrix and vector structs
* as well as functions for loading and manipulating these structs
*/

#include <string.h>
#include "matrix_dims.h"
#include"gsl/gsl_matrix.h"

#include <omp.h>

typedef struct
{
	int rows, cols;
	double** matrix;
} Matrix;

typedef struct
{
	int length;
	double* vector;
} Vector;

/* reading and loading vectors and matrices */
void allocate_Matrix(Matrix* M);
void allocate_Vector(Vector* V);

void read_csv(const char* path);
void read_matrix(const char* path);

/* this is the gsl_matrix_fscanf version of load_matrix and load_vector */
//void load_vector(const char* path, gsl_vector* vec);
//void load_matrix(const char* path, gsl_matrix* mat);

void load_vector(const char* path, Vector* V);
void load_matrix(const char* path, Matrix* M);

/* other manipulations */
void get_matrix_dims(const char* path, int* N, int* M);
void matrix_to_vector(Matrix* M, Vector* V);

#endif // !MV