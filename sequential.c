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
#include "mmio.h"
#include <time.h>
#include <limits.h>

#define SIZE 50000

// A struct used to turn sparse matrices to CSR data structures.
// The notation and algorithm used is taken directly from the given Wikipedia page.
typedef struct {
	long size;
	int *values;
	long *colIndex;
	long *rowIndex;
} csr;


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


// csr readMatToCSR() {
	
// }

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


void printTable(int **table, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			printf(" %d ", table[i][j]);
		}
		printf("\n");
	}

	printf("\n");
}


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
	free(rowIndex);
	free(values);
	free(colIndex);
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



// Calculates the square of CSR matrix, as long as it's square.
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
					 

					newValues[newNonzeros] = square.values[k];
					newColIndex[newNonzeros] = csrTable.colIndex[j];
					newNonzeros++;
					
					break;
				}//else { 
					//csrTable.values[j] = 0;
				//}
			}
		}

		if (i == size-1) {
			newRowIndex[size] = newNonzeros;
		}
	}

	csr hadamard = {size, newValues, newColIndex, newRowIndex};

	return hadamard;
}


int main(int argc, char **argv) {
	 
	// time_t start , end ; 
	//double timediff ;

	//time(&start);
	 //int **table = makeRandomSparseTable(SIZE);
	// printTable(table,SIZE);

	 //csr converted = matrixToCSR(table, SIZE);
	 //printCSR(converted, SIZE);

	 //csr squareNew = csrSquare(converted, SIZE);
	 //printCSR(squareNew, SIZE);
	
	 //int **mat = CSRtoMatrix(squareNew, SIZE);
	 //printTable(mat, SIZE);

	 //csr handa = newhadamard(converted , squareNew , SIZE);
	 //printCSR(handa,SIZE);

	 //int **handamat = CSRtoMatrix( handa, SIZE);
	 //printTable(handamat,SIZE);

	//time(&end);
	FILE *matrixFile;
	int M, N, nz;
	MM_typecode t;
	int firstElements[10];

	matrixFile = fopen("mycielskian4.mtx", "r");
	int banner = mm_read_banner(matrixFile, &t);
	int result = mm_read_mtx_crd_size(matrixFile, &M, &N, &nz);

	//printf("banner: %d\tresult: %d\n", banner, result);

	/*for (int i = 0; i < 10; i++) {
		fscanf(matrixFile, "%d", &firstElements[i]);
		printf(" %d ", firstElements[i]);
	}*/

	//printf(" %d \n" , nz);

	int *row , *col ;

	row = (int *) malloc(2*nz * sizeof(int));	
	col = (int *) malloc(2*nz * sizeof(int));

	for(int i = 0 ; i < nz ; i++){

	fscanf(matrixFile , "%d  %d \n" , &row[i] , &col[i]);
		// 1-based 
		row[i]--;
		col[i]--;
		//symetric matrix
		row[i+nz] = col[i];
		col[i+nz] = row[i];
	}

	
	
	csr Table ; 
	Table.rowIndex = (long *) malloc((M+1)* sizeof(long));
	Table.colIndex = (long *) malloc(2*nz * sizeof(long));
	Table.values = (int *) malloc(2*nz * sizeof(int));

	

	
	int current = 0 ;

	for(int i = 0 ; i < M; i++){
		for(int j = 0 ; j <2*nz ; j++){
			if(row[j] == i ){
				Table.rowIndex[i+1]++;
				Table.values[current] = 1 ;
				Table.colIndex[current] = col[j];
				current ++ ;
			}
			Table.rowIndex[i+2] = Table.rowIndex[i+1];
		}
	}

	printCSR(Table,M);

	csr square = csrSquare(Table , M);

	printCSR(square,M);

	csr hadamard = newhadamard(Table , square , M);

	printf("%d",hadamard.values[5]);
	printCSR(hadamard,M);


	//timediff = difftime(end , start); 
	//printf("%f \n", timediff );

	

	return 0;
}
