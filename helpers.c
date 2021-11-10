#ifndef HELPERS_H
#define HELPERS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "headers/mmio.h"
#include "headers/csr.h"
#include "headers/helpers.h"

// Reads an .mtx file and outputs the resulting CSR form.
csr readmtx(char *mtx, MM_typecode t, int N, int M, int nz) {
	FILE *matrixFile = fopen(mtx, "r");
	int banner = mm_read_banner(matrixFile, &t);
	int result = mm_read_mtx_crd_size(matrixFile, &M, &N, &nz);

	printf("\nbanner: %d\tresult: %d\tnonzeros: %d\tM: %d\tN: %d\n",
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

	int *row, *col;
	row = (int *) malloc(2 * nz * sizeof(int));	
	col = (int *) malloc(2 * nz * sizeof(int));

	for(int i = 0; i < nz; i++) {
		fscanf(matrixFile , "%d  %d \n", &row[i], &col[i]);
		// Decrease the values since Matlab is 1-based. 
		row[i]--;
		col[i]--;

		// Add the symmetric values, since the .mtx files contains only half the matrix.
		row[i+nz] = col[i];
		col[i+nz] = row[i];
	}
	
	csr Table;
	Table.size = N;
	Table.rowIndex = (long *) malloc((M+1) * sizeof(long));
	Table.colIndex = (long *) malloc(2 * nz * sizeof(long));
	Table.values = (int *) malloc(2 * nz * sizeof(int));

	int currentEntry = 0;

	for (int i = 0; i < M; i++) {
		for (int j = 0; j < 2*nz; j++) {
			if (row[j] == i) {
				Table.rowIndex[i+1]++;
				Table.values[currentEntry] = 1;
				Table.colIndex[currentEntry] = col[j];
				currentEntry++;
			}
			Table.rowIndex[i+2] = Table.rowIndex[i+1];
		}
	}

	return Table;
}


// Calculates the square of CSR matrix, as long as it's square.
// FINAL VERSION OF THE FUNCTION.
csr csrSquare(csr table, long size) {
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


long *countTriangles(csr C) {
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


// Matrix multiplication. Only need rows1, cols1 and cols2, because
// cols1==rows2 is required. The new matrix is of size rows1 x cols2.
int **matmul (int **table1, int **table2, long rows1, long cols1, long cols2) {
	int **multTable = (int **) malloc(rows1 * sizeof(int*));

	for (int i = 0; i < rows1; i++) {
		multTable[i] = (int *) malloc(rows1 * sizeof(int));
	}

	// Set the table values to 0 to avoid garbage initializations. 
	for (int i = 0; i < rows1; i++) {
		for (int j = 0; j < cols2; j++) {
			multTable[i][j] = 0;
		}
	}

	for (int i = 0; i < rows1; i++) {
		for (int j = 0; j < cols2; j++) {
			for (int k = 0; k < cols1; k++) {
				multTable[i][j] += table1[i][k] * table2[k][j];
			}
		}
	}

	return multTable;
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
	free(newValues);
	free(newColIndex);
	free(newRowIndex);
	return product;
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




// Converts a sparse matrix to CSR.
csr matrixToCSR(int **table, long size) {
	// Keep track of the nonzero objects.
	long nonzeros = 0; 
	long resizes = 10;

	// Initialize the 3 necessary arrays. row_index has a standard length of "rows+1"
	long *rowIndex = (long *) malloc((size+1) * sizeof(long));
	for (int i = 0; i < size+1; i++) {
		rowIndex[i] = 0;
	}

	// values and col_index have a length of the total nonzero values.
	int *values = (int *) malloc(10 * size * sizeof(int));
	long *colIndex = (long *) malloc(10 * size * sizeof(long));
	for (int i = 0; i < 10*size; i++) {
		values[i] = 0;
		colIndex[i] = 0;
	}
	
	for (long i = 0; i < size; i++) {
		// Add the nonzero values that are placed above the current row. 
		rowIndex[i] = nonzeros;

		for (long j = 0; j < size; j++) {
			if (table[i][j] != 0) {
				nonzeros++;

				// Make sure to extend the arrays in case their size exceeds the current one.
				if (nonzeros == resizes * size -1) {
					resizes += 10;
					values = (int *) realloc(values, (resizes * size) * sizeof(int));
					colIndex = (long *) realloc(colIndex, (resizes * size) * sizeof(long));
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

	csr converted;
	converted.size = size;
	converted.values = values;
	converted.colIndex = colIndex;
	converted.rowIndex = rowIndex;
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