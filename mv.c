#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mv.h"

void allocate_Vector(Vector* V)
{
	/* needs factor of 2 for terminating string during load */
	double* vec = malloc(V->length * 2 * (sizeof(vec)));
	
	V->vector = vec;
}

void allocate_Matrix(Matrix* M)
{
	double** mat = malloc(M->rows * 2 * sizeof(*mat));
	if (mat)
	{
		for (int i = 0; i < M->rows; i++)
		{
			mat[i] = malloc(M->cols * 2 * (sizeof(double)));
		}
	}
	M->matrix = mat;
}

void read_csv(const char* path)
{
	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		fprintf(stderr, "Unable to open file: %s\n", path);
		exit(1);
	}

	char line[200];

	while (fgets(line, sizeof(line), f))
	{
		char* token;

		token = strtok(line, ",");

		while (token != NULL)
		{
			//printf("%s\n", token);
			int d = sizeof(line);
			printf("%i", d);

			token = strtok(NULL, ",");
		}
		printf("\n");
	}
	fclose(f);

}

void get_matrix_dims(const char* path, int* N, int *M)
{
	/* gsl_matrix_fscanf wouldn't work so this is a rough alternative */
	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		fprintf(stderr, "Unable to open file: %s\n", path);
		exit(1);
	}

	int rows, cols, elems;
	rows = 0;
	cols = 0;
	elems = 0;

	double** array = malloc(1024 * sizeof(double));
	int idx = 0;
	char line[1024 * sizeof(double)];

	while (fgets(line, sizeof(line), f))
	{
		char* token;

		token = strtok(line, ",");

		while (token != NULL)
		{
			//printf("%s\n", token);
			int d = sizeof(line);

			token = strtok(NULL, ",");
			elems++;
		}

		idx++;
		rows++;

	}
	fclose(f);
	cols = elems / rows;
	*M = cols;
	*N = rows;
	printf("num cols: %d\n", cols);
	printf("num rows is %d\n", rows);
}

void read_matrix(const char* path)
{
	/* gsl_matrix_fscanf wouldn't work so this is a rough alternative */
	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		fprintf(stderr, "Unable to open file: %s\n", path);
		exit(1);
	}

	int rows, cols, elems;
	rows = 0;
	cols = 0; 
	elems = 0;
	
	//double* array[3 * sizeof(double)];

	double **array = malloc(1024 * sizeof(double));
	int idx = 0;
	char line[1024*sizeof(double)];

	while (fgets(line, sizeof(line), f))
	{
		char* token;

		token = strtok(line, ",");

		while (token != NULL)
		{
			//printf("%s\n", token);
			int d = sizeof(line);

			token = strtok(NULL, ",");
			elems++;
		}

		if (array)
		{
			array[idx] = (double*)malloc(1024*sizeof(array));

			int j;
			char* ptr = NULL;
			char* buffer = NULL;
			for (j = 0, ptr = line; j < 3; j++, ptr++)
			{
				if (array[idx])
				{
					array[idx][j] = strtod(ptr, &ptr);
				}
			}
		}
		idx++;
		rows++;
		
	}
	fclose(f);
	cols = elems / rows;
	printf("num cols: %d\n", cols);
	printf("num rows is %d\n", rows);

	//int i,j;
	/*
	for (i = 0; i < idx; i++) {
		printf("\narray[%d][] =", i);

		for (j = 0; j < 30; j++)
			printf(" %lf", array[i][j]);
	}
	*/
}

/*
void load_vector(const char* path, gsl_vector* vec)
{
	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		fprintf(stderr, "Unable to open file: %s\n", path);
		exit(1);
	}

	gsl_vector_fscanf(f, vec);
	fclose(f);
}

*/

/*
void load_matrix(const char* path, gsl_matrix* mat)
{
	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		perror("Unable to read u.csv");
		exit(1);
	}

	gsl_matrix_fscanf(path, mat);
	fclose(f);
}
*/

void load_vector(const char* path, Vector *V)
{
	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		fprintf(stderr, "Unable to open file: %s\n", path);
		exit(1);
	}

	int idx = 0;
	char line[1024 * sizeof(double)];

	while (fgets(line, sizeof(line), f))
	{
		char* token;
		double elem;
		token = strtok(line, ",");
		//printf("line %s\n", line);
		elem = strtod(token, NULL);
		V->vector[idx] = elem;
		//printf("size of elem: %d\n", sizeof(elem));
		//printf("size of v: %d\n", sizeof(V->vector[idx]));
		//printf("%d\n", idx);

		//printf("v[%d] = %lf\n", idx, V->vector[idx]);
		//printf("elem = %lf\n", elem);
		while (token != NULL)
		{
			token = strtok(NULL, ",");
		}
		idx++;
	}
	//printf("id: %d\n",idx);
	fclose(f);
}

void load_matrix(const char* path, Matrix *M)
{
	FILE* f = fopen(path, "r");
	if (f == NULL)
	{
		fprintf(stderr, "Unable to open file: %s\n", path);
		exit(1);
	}

	int i = 0;
	char line[1024 * sizeof(double)];

	while (fgets(line, sizeof(line), f))
	{
		char* token;

		token = strtok(line, ",");
		/*
		while (token != NULL)
		{
			//printf("%s\n", token);
			int d = sizeof(line);

			token = strtok(NULL, ",");
		}
		*/
		if (M->matrix)
		{
			/* allocate memory for next row */
			//array[idx] = (double*)malloc(1024 * sizeof(array));

			int j;
			char* ptr = NULL;
			char* buffer = NULL;
			for (j = 0, ptr = line; j < M->cols; j++, ptr++)
			{
				if (M->matrix[i])
				{
					M->matrix[i][j] = strtod(ptr, &ptr);
					//printf("array[%d][%d] = %e\n", i, j, M->matrix[i][j]);
				}
			}
		}
		i++;
	}
	fclose(f);

}

void matrix_to_vector(Matrix* M, Vector* V)
{
	/*
	convert matrix to vector

	input:
		matrix: matrix to be converted
		rows: rows of matrix
		cols: columns of matrix
		vector: new vector resulting from matrix
	*/
	int k = 0; /* index of new vector */
	for (int i = 0; i < M->rows; i++)
	{
		for (int j = 0; j < M->cols; j++)
		{
			V->vector[k] = M->matrix[i][j];
			//printf("%lf", vector[k]);
			k++;
		}
	}

	for (int i = 0; i < V->length; i++)
	{
		//printf("new vec[%d] = %e\n", i, V->vector[i]);
	}


}