#ifndef CSR_H
#define CSR_H

#include <stdio.h>

// A struct used to turn sparse matrices to CSR data structures.
// The notation and algorithm used is taken directly from the given Wikipedia page.
typedef struct {
	uint size;
	int *values;
	uint *colIndex;
	uint *rowIndex;
} csr;

#endif

