#ifndef CSR_OPERATIONS_H
#define CSR_OPERATIONS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "headers/csr_operations.h"
#include "headers/csr.h"
#include "headers/mmio.h"


// Reads an .mtx and directly converts it to CSR.
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
		fscanf(matrixFile, "%d", &row);
		fscanf(matrixFile, "%d", &column);
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


// Calculates the square of CSR matrix, as long as it's square.
// FINAL VERSION OF THE FUNCTION.
csr csrSquare(csr table, long size) {
	long nonzeros = table.rowIndex[size];
	int resizes = 10;
	long newNonzeros = 0;

	// The new values array. Intialize all to 0. Do the same for the new column indices.
	int *newValues = (int *) malloc(10 * size * sizeof(int));
	long *newColIndex = (long *) malloc(10 * size * sizeof(long));
	for (long i = 0; i < 10 * size; i++) {
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

		if (newNonzeros != 0 && (newNonzeros % (10*nonzeros) == 0)) {
			resizes += 10;
			newValues = (int *) realloc(newValues, size * resizes * sizeof(int));
			newColIndex = (long *) realloc(newColIndex, size * resizes * sizeof(long));
		}

		// Scan and multiply -only nonzero elements- every single column. Since the matrix is
		// symmetric, we actually scan every line from top to bottom,
		// getting exactly the same result.
		for (long column = 0; column < size; column++) {
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

	csr square = {size, newValues, newColIndex, newRowIndex};
	return square;
}


csr csrSquareAlt(csr converted, int **table, long size) {
	long nonzeros = converted.values[size];

	// The new values array. Intialize all to 0. Do the same for the new column indices.
	int *newValues = (int *) malloc(10 * size * sizeof(int));
	long *newColIndex = (long *) malloc(10 * size * sizeof(long));
	for (long i = 0; i < 10 * nonzeros; i++) {
		newValues[i] = 0;
		newColIndex[i] = 0;
	}

	// Finally, initialize the new row index array.
	long *newRowIndex = (long *) malloc((size+1) * sizeof(long));
	for (long i = 0; i < size+1; i++) {
		newRowIndex[i] = 0;
	}

	// Keep the new matrix's total of nonzero values.
	long newNonzeros = 0;
	// Keep tab of how many times the newValues and colIndex arrays have been resized.
	int resizes = 10; 

	// Scan every row using the row_index array.
	for (long row = 0; row < size; row++) {
		// Make sure to resize the arrays if they're full.
		if (newNonzeros != 0 && (newNonzeros % (10*nonzeros) == 0)) {
			resizes += 10;
			newValues = (int *) realloc(newValues, size * resizes * sizeof(int));
			newColIndex = (long *) realloc(newColIndex, size * resizes * sizeof(long));
		}

		newRowIndex[row] = newNonzeros;
		// printf("Row: %ld\t Nonzeros = %ld\n", row, newRowIndex[row]);

		// Get the indices of the values that lie within each row.
		long start = converted.rowIndex[row];
		long end = converted.rowIndex[row+1];
		// printf("Start = %ld\t End = %ld\n", start, end);

		for (long column = 0; column < size; column++) {
			int cellValue = 0;

			for (long element = start; element < end; element++) {
				long elementRow = converted.colIndex[element];
				cellValue += converted.values[element] * table[elementRow][column];
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

	csr product = {size, newValues, newColIndex, newRowIndex};
	return product;
}


// CSR-CSR Hadamard (element-wise) operation.
// FINAL VERSION OF THE FUNCTION.
csr newhadamard(csr csrTable, csr square, long size) {
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
				if(csrTable.colIndex[j] == square.colIndex[k]){
					csrTable.values[j] = square.values[k];

					newValues[newNonzeros] = csrTable.values[j];
					newColIndex[newNonzeros] = csrTable.colIndex[j];
					newNonzeros++;
					
					break;
				} else { 
					csrTable.values[j] = 0;
				}
			}
		}

		if (i == size-1) {
			newRowIndex[size] = newNonzeros;
		}
	}

	csr hadamard = {size, newValues, newColIndex, newRowIndex};
	return hadamard;
}


csr hadamard(csr csrTable, int **square, long size) {
	long nonzeros = csrTable.rowIndex[size];

	long *newRowIndex = (long *) malloc((size+1) * sizeof(long));
	long *newColIndex = (long *) malloc(nonzeros * sizeof(long));
	int *newValues = (int *) malloc(nonzeros * sizeof(int));
	long newNonzeros = 0;

	for (int i = 0; i < nonzeros; i++) {
		newColIndex[i] = 0;
		newValues[i] = 0;
	}

	for (int i = 0; i < size; i++) {
		newRowIndex[i] = 0;
	}

	for (long row = 0; row < size; row++) {
		newRowIndex[row] = newNonzeros;

		long start = csrTable.rowIndex[row];
		long end = csrTable.rowIndex[row+1];

		for (int i = start; i < end; i++) {
			long column = csrTable.colIndex[i];
			int value = csrTable.values[i];

			int cellValue = value * square[row][column];
			if (cellValue > 0) {
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

	csr hadamard = {size, newValues, newColIndex, newRowIndex};
	return hadamard;
}


#endif