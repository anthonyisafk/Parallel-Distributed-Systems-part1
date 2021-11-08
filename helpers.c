#ifndef HELPERS_H
#define HELPERS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "headers/mmio.h"
#include "headers/csr.h"
#include "headers/helpers.h"


// Converts a sparse matrix to CSR.
csr matrixToCSR(int **table, long size) {
	// Keep track of the nonzero objects.
	long nonzeros = 0; 

	// Initialize the 3 necessary arrays. row_index has a standard length of "rows+1"
	long *rowIndex = (long *) malloc((size+1) * sizeof(long));
	// values and col_index have a length of the total nonzero values.
	int *values = (int *) malloc(100 * sizeof(int));
	long *colIndex = (long *) malloc(100 * sizeof(long));
	
	for (long i = 0; i < size; i++) {
		// Add the nonzero values that are placed above the current row. 
		rowIndex[i] = nonzeros;

		for (long j = 0; j < size; j++) {
			if (table[i][j] != 0) {
				nonzeros++;

				// Make sure to extend the arrays in case their size exceeds the current one.
				if (nonzeros != 0 && (nonzeros % 100 == 0)) {
					values = (int *) realloc(values, (nonzeros+100) * sizeof(int));
					colIndex = (long *) realloc(colIndex, (nonzeros+100) * sizeof(long));
				}

				// Add the nonzero value and the respective column index.
				values[nonzeros-1] = table[i][j];
				colIndex[nonzeros-1] = j; 
			}
		}

		// Add the last value before exiting the for loop.
		if (i == size - 1) {
			rowIndex[size] = nonzeros;
		}
	}

	csr converted = {size, values, colIndex, rowIndex};
	return converted;
}


int **CSRtoMatrix(csr table, long size) {
	// Initialize the new matrix and set everything to 0.
	int **matrix = (int **) malloc(size * sizeof(int*));
	for (int i = 0; i < size; i++) {
		matrix[i] = (int *) malloc(size * sizeof(int));
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			matrix[i][j] = 0;
		}
	}

	// Scan each row for nonzero values and plug them into the matrix.
	for (long row = 0; row < size; row++) {
		long start = table.rowIndex[row];
		long end = table.rowIndex[row+1];

		for (int j = start; j < end; j++) {
			long column = table.colIndex[j];
			int value = table.values[j];
			matrix[row][column] = value;
		}
	}

	// Done.
	return matrix;
}


/**
 * ---------- HELPERS ----------
 */ 
// Prints a CSR data structure.
void printCSR(csr converted, long size) {
	long nonzeros = converted.rowIndex[size];
	
	printf("Values:");
	for (int i = 0; i < nonzeros; i++) {
		printf(" %d ", converted.values[i]);
	}

	printf("\nCol_index:");
	for (int i = 0; i < nonzeros; i++) {
		printf(" %ld ", converted.colIndex[i]);
	}

	printf("\nRow_index:");
	for (int i = 0; i < size+1; i++) {
		printf(" %ld ", converted.rowIndex[i]);
	}
	printf("\n\n");
}


int **makeRandomSparseTable(int size) {
	int **table = (int **) malloc(size * sizeof(int *));
	for (int i = 0; i < size; i++) {
		table[i] = (int *) malloc(size * sizeof(int));
	}

	// Only add values to the upper triangle of the table.
	for (int i = 0; i < size; i++) {
		for (int j = i; j < size; j++) {
			table[i][j] = (rand() % 1000 < 750) ? 0 : 1; 

			// ALTERNATIVE INITIALIZATION. USE 2's INSTEAD
			// OF ONLY ONES AND ZEROS.

			// int newValue = rand() % 1000;
			// if (newValue < 750) {
			// 	table[i][j] = 0;
			// } else if (newValue >= 750 && newValue <= 930) {
			// 	table[i][j] = 1;
			// } else {
			// 	table[i][j] = 2;
			// }
		}
	}

	// Make it symmetric by setting A[i][j] = A[j][i] for the rest of the values.
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < i; j++) {
			table[i][j] = table[j][i];
		}
	}

	return table;
}


// Matrix multiplication ONLY FOR SQUARE MATRICES.
int **matmul (int **table1, int **table2, int size) {
	int **multTable = (int **) malloc(size * sizeof(int*));

	for (int i = 0; i < size; i++) {
		multTable[i] = (int *) malloc(size * sizeof(int));
	}

	// Set the table values to 0 to avoid garbage initializations. 
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			multTable[i][j] = 0;
		}
	}

	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			for (int k = 0; k < size; k++) {
				multTable[i][j] += table1[i][k] * table2[k][j];
			}
		}
	}

	return multTable;
}


void printTable(int **table, long size) {
	for (long i = 0; i < size; i++) {
		for (long j = 0; j < size; j++) {
			printf(" %d ", table[i][j]);
		}
		printf("\n");
	}

	printf("\n");
}

#endif