/*
 * sequential.c 
 * Convert a square N x N matrix into the CSR format, made for sparse matrices:
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 * 
 * Authors: Antonios Antoniou - 9482
 *          Efthymios Grigorakis - 9694
 *    
 * 2021 Aristotle University of Thessaloniki
 * Parallel and Distributed Systems.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>

#include "headers/csr.h"
#include "headers/mmio.h"
#include "headers/helpers.h"

#define SIZE 6


csr *readmtx(char *mtx, MM_typecode t, int N, int M, int nz) {
  printf("\nInitialize file.\n");
	FILE *matrixFile = fopen(mtx, "r");
  printf("\nRead the file.\n");
  // sleep(1);
  // printf("\nBack to reading...\n");
	int banner = mm_read_banner(matrixFile, &t);
  printf("\nRead banner.\n");
	int result = mm_read_mtx_crd_size(matrixFile, &M, &N, &nz);

	printf("\nbanner: %d\tresult: %d\tnonzeros: %d\tM: %d\tN: %d\n",
	  banner, result, nz, M, N);
 
	// Display error messages and abort, if the matrix isn't square or hasn't been read properly.
	if (N != M) {
		printf("N and M are not equal. The matrix isn't square. Aborting...");
		csr *returnError = {0, NULL, NULL, NULL};
		return returnError;
	}

	if (banner != 0 || result != 0) {
		printf("Error. Couldn't process the .mtx file!");
		csr *returnError = {0, NULL, NULL, NULL};
		return returnError;	
	}

	// Calculate the average nonzero elements per row. Use that number for the initial malloc.
	int averagePerRow = (int) (nz / N + 1);
  int initialSize = 2 * averagePerRow;
	printf("\naveragePerRow = %d\n", averagePerRow);

	// Store the column of each nonzero value in an array of arrays.
	// Each sub-array corresponds to a row of the matrix.
	int **valuesByRow = (int **) malloc(N * sizeof(int *));
  printf("\nIntialized valuesByRow\n");

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
  printf("\nAll arrays initialized\n");

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
  printf("\nThe CSR arrays as well...\n");

	long nonzeros = 0;
	int row = 0;
  int column = 0; // variables used to scheme through the .mtx file.

	for (int i = 0; i < nz; i++) {
		// Read the values and decrease them by one, since Matlab is 1-index based.
		fscanf(matrixFile, "%d %d\n", &row, &column);
    printf("Row = %d\tColumn = %d\t", row, column);
		row--;
		column--;
		int newPosition = valuesByRow[row][0]; 
    printf("valuesByRow[%d][0] = %d", row, valuesByRow[row][0]);

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
	
	csr *converted;
  converted->size= N;
  converted->values = values;
  converted->colIndex = colIndex;
  converted->rowIndex = rowIndex;
	return converted;
}



int main(int argc, char **argv) {
  // int **test = makeRandomSparseTable(SIZE);
  // printTable(test, SIZE);
  // csr testCSR = matrixToCSR(test, SIZE);

  // csr testSquare = csrSquare(testCSR, SIZE);
  // printCSR(&testSquare, SIZE);

  // csr C = newhadamard(testCSR, testSquare, SIZE);
  // printCSR(&C, SIZE);

  FILE *matrixFile;
  int M, N, nz;
  MM_typecode t;
  char *mtx = "tables/mycielskian4.mtx";

  readmtx(mtx, t, N, M, nz);


}