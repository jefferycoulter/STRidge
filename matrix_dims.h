#ifndef DIMS
#define DIMS

#define LEN(A) sizeof(A) / sizeof(A[0])
#define NROWS(A) LEN(A)
#define NCOLS(A) LEN(A[0])

#endif // !DIMS
