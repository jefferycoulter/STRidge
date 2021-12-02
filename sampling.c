#include "sampling.h"

int randomInt(int min, int max)
{
	/* function to generate random point between min and max */

	int p = rand() / (int)RAND_MAX;

	return min + p * (max - min);
}

void fill(int* idx, int length)
{
	/* initialize array of random indices. this ensures no repetition */
	
	for (int i = 0; i < length; i++)
	{
		idx[i] = i;
	}
}

void shuffle(int* idx, int length)
{
	/* shuffle the array of indices to create a random order */
	
	for (int i = 0; i < length; i++)
	{
		int temp = idx[i];
		int ri = rand() % length;
		idx[i] = ri;
		//printf("idx[%d] = %d\n", i, ri);
		idx[ri] = temp;
		//printf("idx[%d] = %d\n", ri, temp);
	}
}

void sample_data(Matrix* X_train, Matrix* X_test, Vector* y_train, Vector* y_test, Matrix* T, Vector* sol, int N_THREADS)
{	
	/* 
	this function divides the T and solution into training and test sets for
	linear regression and ridge regression

	input:
		X_train: empty matrix; portion of T used for training
		X_test: empty matrix; portion of T used for testing
		y_train: empty vector; portion of solution equated to X_train
		y_test: empty vector; portion of solution equated to X_test
		T: this is the Theta matrix
		sol: u_t, this is the "solution" in the since that u_t = Du where D 
			 is the differential operator we are looking for
	*/

	int thread_id;
	
	int* ii = malloc(T->rows * sizeof(int));
	//int* jj = malloc(T->cols * sizeof(int));
	
	fill(ii, T->rows);
	//fill(&jj, T->cols);

	shuffle(ii, T->rows);
	//shuffle(&jj, T->cols);

	omp_set_num_threads(N_THREADS);

	int i, j;

	#pragma omp parallel private(thread_id, i, j)
	#pragma omp for
	for (i = 0; i < X_train->rows; i++)
	{
		for (j = 0; j < X_train->cols; j++)
		{
			/* the training rows are the 0,...,train->rows-1 elements of ii
			and the training columns are the 0,...,train->cols-1 elements of jj*/
			X_train->matrix[i][j] = T->matrix[ii[i]][j];
			y_train->vector[i] = sol->vector[ii[i]];
		}
	}

	#pragma omp for
	for (i = 0; i < X_test->rows; i++)
	{
		for (j = 0; j < X_test->cols; j++)
		{
			/* the training rows are the train->rows,...T->rows elements of ii
			and the training columns are the train->cols,...,T->cols elements of jj*/
			X_test->matrix[i][j] = T->matrix[ii[i+X_train->rows]][j];
			y_test->vector[i] = sol->vector[ii[i+X_train->rows]];
		}
	}


}