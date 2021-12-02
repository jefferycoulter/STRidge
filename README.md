Implementation of sequential-threshold ridge regression using OpenMP and GSL.  For details, see Dr. Kutz's paper: https://www.science.org/doi/10.1126/sciadv.1602614

To do:
-- convert fully to GSL syntax so there aren't both custom matrix/vector structs and gsl matrix/vectors
-- implement MPI and testing
-- extend to possible derivative combinations (i.e. three spatial dimensions, 3+1 space-time)
-- gui?