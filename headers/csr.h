#include <stdio.h>

// A struct used to turn sparse matrices to CSR data structures.
// The notation and algorithm used is taken directly from the given Wikipedia page.
typedef struct {
	long size;
	int *values;
	long *colIndex;
	long *rowIndex;
} csr;
