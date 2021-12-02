#ifndef SAMPLING
#define SAMPLING

#include "mv.h"

//#include <omp.h>

int randomInt(int min, int max);
void fill(int* idx, int length);
void shuffle(int* idx, int length);
void sample_data(Matrix* X_train, Matrix* X_test, Vector* y_train, Vector* y_test, Matrix* T, Vector* sol, int N_THREADS);

#endif // !SAMPLING
