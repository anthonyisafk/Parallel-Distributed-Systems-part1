#ifndef CILK_HELPERS_H
#define CILK_HELPERS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "../headers/mmio.h"
#include "../headers/csr.h"
#include "../headers/helpers.h"
#include "../headers/cilk_helpers.h"

#define MAX_THREADS "4"
int max_threads = 4;

// Calculates the square of CSR matrix, as long as it's square.
// FINAL VERSION OF THE FUNCTION.
csr csrSquareCilk(csr table, long size) {
  __cilkrts_end_cilk();
  __cilkrts_set_param("nworkers", MAX_THREADS);
  __cilkrts_init();

	long nonzeros = table.rowIndex[size];
	int resizes = 5;
	long newNonzeros = 0;

	// The new values array. Intialize all to 0. Do the same for the new column indices.
	int *newValues = (int *) malloc(resizes * size * sizeof(int));
	long *newColIndex = (long *) malloc(resizes * size * sizeof(long));

	for (long i = 0; i < resizes * size; i++) {
		newValues[i] = 0;
		newColIndex[i] = 0;
	}

	// Finally, initialize the new row index array.
	long *newRowIndex = (long *) malloc((size+1) * sizeof(long));
	for (long i = 0; i < size+1; i++) {
		newRowIndex[i] = 0;
	}

	for (long row = 0; row < size; row++) {
		newRowIndex[row] = newNonzeros;
		long start = table.rowIndex[row];
		long end = table.rowIndex[row+1];

		// Scan and multiply -only nonzero elements- every single column. Since the matrix is
		// symmetric, we actually scan every line from top to bottom,
		// getting exactly the same result.
		for (long column = 0; column < size; column++) {

			// Make sure to check the table for resizes regularly.
			if (newNonzeros == resizes * size - 1) {
				resizes += 5;
				newValues = (int *) realloc(newValues, size * resizes * sizeof(int));
				newColIndex = (long *) realloc(newColIndex, size * resizes * sizeof(long));
			}

			int cellValue = 0;
			long columnStart = table.rowIndex[column];
			long columnEnd = table.rowIndex[column+1];

			// This variable is used to leave a trace at the last element where the row-column pair
			// was a match, for the element-wise multiplication to take place. Since the order here 
			// is increasing, we leave a trace where the last match was found and the scan 
			// for another matching pair starts at the next index.
			int lastMatch = columnStart - 1;

			for (long rowElement = start; rowElement < end; rowElement++) {
				for (long colElement = lastMatch+1; colElement < columnEnd; colElement++) {
					if (table.colIndex[colElement] == table.colIndex[rowElement]) {
						lastMatch = colElement; // Mark the last match.

						// Add the product to the cell value.
						cellValue += table.values[colElement] * table.values[rowElement];
					}
				}
			}

			if (cellValue > 0) {
				// printf("Cell value: %ld\n", cellValue);
				newValues[newNonzeros] = cellValue;
				newColIndex[newNonzeros] = column;
				newNonzeros++;
			}
		}

		// Add the final value to newRowIndex before exiting the loop.
		if (row == size - 1) {
			newRowIndex[size] = newNonzeros;
		}
	}

	// Resize the arrays to save space and avoid null values.
	newValues = (int *) realloc(newValues, newNonzeros * sizeof(int));
	newColIndex = (long *) realloc(newColIndex, newNonzeros * sizeof(long));

	csr square;
	square.size = size;
	square.values = newValues;
	square.colIndex = newColIndex;
	square.rowIndex = newRowIndex;
	return square;
}


// CSR-CSR Hadamard (element-wise) operation.
// FINAL VERSION OF THE FUNCTION.
csr newhadamardCilk(csr csrTable, csr square, long size) {
	long oldNonzeros = csrTable.rowIndex[size];
	long newNonzeros = 0;

	long *newColIndex = (long *) malloc(oldNonzeros * sizeof(long));
	int *newValues = (int *) malloc(oldNonzeros * sizeof(int));
	for (int i = 0; i < oldNonzeros; i++) {
		newColIndex[i] = 0;
		newValues[i] = 0;
	}

	long *newRowIndex = (long *) malloc((size+1) * sizeof(long));
	for (int i = 0; i < size+1; i++) {
		newRowIndex[i] = 0;
	}

	// Use the old table and only transfer values to the new one if they are greater than 0.
	for (int i = 0; i < size; i++) {
		newRowIndex[i] = newNonzeros;

		for (int j = csrTable.rowIndex[i]; j < csrTable.rowIndex[i+1]; j++) {
			for (int k = square.rowIndex[i]; k < square.rowIndex[i+1]; k++) {
				if(csrTable.colIndex[j] == square.colIndex[k]) {
					csrTable.values[j] = square.values[k];

					newValues[newNonzeros] = csrTable.values[j];
					newColIndex[newNonzeros] = csrTable.colIndex[j];
					newNonzeros++;
					
					break;
				}
			}
		}

		if (i == size-1) {
			newRowIndex[size] = newNonzeros;
		}
	}

	// Resize the arrays to save space and avoid null values.
	newValues = (int *) realloc(newValues, newNonzeros * sizeof(int));
	newColIndex = (long *) realloc(newColIndex, newNonzeros * sizeof(long));

	csr hadamard;
	hadamard.size = size;
	hadamard.values = newValues;
	hadamard.colIndex = newColIndex;
	hadamard.rowIndex = newRowIndex;
	return hadamard;
}


long *countTrianglesCilk(csr C) {
	long size = C.size;
	long *triangleCount = (long *) malloc(size * sizeof(long));
	for (long i = 0; i < size; i++) {
		triangleCount[i] = 0;
	}

	// Add all the values in each row, then divide by 2.
	// Simulates the operation of multiplying the table with a nx1 vector.
	for (long i = 0; i < size; i++) {
		long start = C.rowIndex[i];
		long end = C.rowIndex[i+1];

		for (long j = start; j < end; j++) {
			triangleCount[i] += C.values[j];
			
			// Divide by 2 if we reached the last value of the current row.
			if (j == end - 1) {
				triangleCount[i] /= 2;
			}
		}
	}

	return triangleCount;
}



#endif