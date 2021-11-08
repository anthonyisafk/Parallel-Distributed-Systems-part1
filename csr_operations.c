#ifndef CSR_OPERATIONS_H
#define CSR_OPERATIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "headers/csr_operations.h"
#include "headers/csr.h"
#include "headers/mmio.h"

csr readmtx(char *mtx, MM_typecode *t, int N, int M, int nz) {
	FILE *matrixFile;
	matrixFile = fopen(mtx, "r");
	int banner = mm_read_banner(matrixFile, t);
	int result = mm_read_mtx_crd_size(matrixFile, &M, &N, &nz);

	printf("banner: %d\tresult: %d\tnonzeros: %d\tM: %d\tN: %d\n",
	banner, result, nz, M, N);
 
	// Display error messages and abort, if the matrix isn't square or hasn't been read properly.
	if (N != M) {
		printf("N and M are not equal. The matrix isn't square. Aborting...");
		csr returnError = {0, NULL, NULL, NULL};
		return returnError;
	}

	if (banner != 0 || result != 0) {
		printf("Error. Couldn't process the .mtx file!");
		csr returnError = {0, NULL, NULL, NULL};
		return returnError;	
	}

	// Calculate the average nonzero elements per row. Use that number for the initial malloc.
	int averagePerRow = (int) (nz / N + 1);
  int initialSize = 2 * averagePerRow;
	printf("averagePerRow = %d\n", averagePerRow);

	// Store the column of each nonzero value in an array of arrays.
	// Each sub-array corresponds to a row of the matrix.
	int **valuesByRow = (int **) malloc(N * sizeof(int *));

	// Fill the valuesByRow array with INT_MIN, except for the first cell of each array.
	// This position will be used to dictate where to put the next value.
	// INT_MIN helps getting rid of unaltered values before each row is sorted.
	for (int i = 0; i < N; i++) {
		valuesByRow[i] = (int *) malloc(initialSize * sizeof(int));
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < initialSize; j++) {			
			if (j == 0) {
				valuesByRow[i][j] = 0;
			} else {
				valuesByRow[i][j] = INT_MIN; // will be used as an "error code".
			}
		}
	}

	// Initialize the CSR struct's attributes.
	long *rowIndex = (long *) malloc((N+1) * sizeof(long));
	for (int i = 0; i < N+1; i++) {
		rowIndex[i] = 0;
	}

	// The maximum total of nonzero values will be 2*nz (if no element is stored in the main diagonal).
	int *values = (int *) malloc(2 * nz * sizeof(int));
	long *colIndex = (long *) malloc(2 * nz *sizeof(long));
	for (int i = 0; i < nz; i++) {
		values[i] = -1;
		colIndex[i] = -1;
	}

	long nonzeros = 0;
	int row, column; // variables used to scheme through the .mtx file.

	for (int i = 0; i < nz; i++) {
		// Read the values and decrease them by one, since Matlab is 1-index based.
		fscanf(matrixFile, "%d %d\n", &row, &column);
		row--;
		column--;
		int newPosition = valuesByRow[row][0]; 

		if (newPosition > 0 && (newPosition % initialSize == 0)) {
			valuesByRow[row] = (int *) realloc(valuesByRow[row], (newPosition/initialSize + initialSize) * sizeof(int));
		}
		valuesByRow[row][0]++;
		newPosition++;
		valuesByRow[row][newPosition] = column;

		// Add the symmetric value if the element isn't part of the main diagonal.
		if (row != column) {
			newPosition = valuesByRow[column][0]; 

			if (newPosition > 0 && (newPosition % (initialSize) == 0)) {
				valuesByRow[column] = (int *) realloc(valuesByRow[column], (newPosition/initialSize + initialSize) * sizeof(int));
			}
			valuesByRow[column][0]++;
			newPosition++;		
			valuesByRow[column][newPosition] = row;	
		}
	}

  // Turn the arrays we made into a CSR structure by scanning each row.
	for (int row = 0; row < N; row++) {
		rowIndex[row] = nonzeros;

		for (int column = 1; column < valuesByRow[row][0]+1; column++) {
			values[nonzeros] = 1;
			colIndex[nonzeros] = valuesByRow[row][column];
			nonzeros++;
		}

		if (row == N - 1) {
			rowIndex[N] = nonzeros;
		}
	}
	
	csr converted = {N, values, colIndex, rowIndex};
	return converted;
}





#endif