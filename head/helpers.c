#ifndef HELPERS_H
#define HELPERS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "../headers/mmio.h"
#include "../headers/csr.h"
#include "../headers/csr_arg.h"
#include "../headers/helpers.h"


// Reads an mtx file and returns the CSR format by making a 2D array of non-symmetric 
// dimensions. Stores the column indices in each row that corresponds to a row in the
// valuesByRow array. Then increases the size of each sub array by one, since most
// vertices are connected to few other vertices.
csr readmtx_dynamic(char *mtx, MM_typecode *t, int N, int M, int nz) {
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

  int **valuesByRow = (int **) malloc(N * sizeof(int *));
  for (int i = 0; i < N; i++) {
    valuesByRow[i] = (int *) calloc(2, sizeof(int));
  }

  uint nonzeros = 0;
  for (int i = 0; i < nz; i++) {
    int row, col;
    fscanf(matrixFile , "%d  %d \n", &row, &col);

    // Decrease the values since Matlab is 1-index based. 
    row--;
    col--;

    // Keeps track of where the next value should be placed.
    int newEntry = valuesByRow[row][0] + 1;
    valuesByRow[row][newEntry] = col;
    valuesByRow[row][0]++;
    // Extend the array by one element.
    valuesByRow[row] = (int *) realloc(valuesByRow[row], (valuesByRow[row][0] + 2) * sizeof(int));
    nonzeros++;

    // If the scanned element isn't on the main diagonal, add the symmetric value.
    if (row != col) {
      int newEntry = valuesByRow[col][0] + 1;
      valuesByRow[col][newEntry] = row;
      valuesByRow[col][0]++;
      valuesByRow[col] = (int *) realloc(valuesByRow[col], (valuesByRow[col][0] + 2) * sizeof(int));
      nonzeros++;
    }
  }

  // This is a binary matrix. All the nonzero values are 1 by default.
  uint *values = (uint *) malloc(nonzeros *sizeof(int));
  for (int i = 0; i < nonzeros; i++) {
    values[i] = 1;
  }

  // Put the scanned values in a CSR structure.
  uint *rowIndex = (uint *) calloc(N, sizeof(uint));
  uint *colIndex = (uint *) calloc(nonzeros, sizeof(uint));

  for (uint row = 0; row < N; row++) {
    int rowNonzeros = valuesByRow[row][0];

    for (int column = 1; column < rowNonzeros+1; column++) {
      uint currentNonzeros = rowIndex[row + 1];

      colIndex[currentNonzeros] = valuesByRow[row][column];
      rowIndex[row + 1]++;
    }

    // Pass the total of nonzero values to the next row, unless we reached the last one.
    if (row < N - 1) {
      rowIndex[row + 2] = rowIndex[row + 1];
    }
  }

  csr csr_mtx = {N, values, colIndex, rowIndex};
  return csr_mtx;
}


// Creates a csr_arg array, depending on the number of threads and and the CSR table.
csr_arg *makeThreadArguments(csr table, int max_threads) {
  // Make a table of csr_arg structs. Each of them corresponds to a thread.
  uint size = table.size;
  uint nonzeros = table.rowIndex[size];
  uint interval = nonzeros / max_threads;
  csr_arg *csr_args = (csr_arg *) malloc(max_threads * sizeof(csr_arg));

  // Placeholder for the total size of the csr_args array.
  uint totalSize = 0;
  for (int i = 0; i < max_threads; i++) {
    // The starting and ending rows of each csr_arg. The first "max_threads-1" structures have equal size
    // and the last one gets the remaining rows of the integer division "size/max_threads".
    uint start = (i == 0) ? 0 : csr_args[i-1].end;
    uint end = (i != max_threads - 1) ? start + 1 : size;

    if (end != size) {
      for (uint j = start + 1; j < size + 1; j++) {
        if (table.rowIndex[j] > interval * (i + 1)) {
          end = j;
          break;
        }
      }
    }

    uint partialSize = end - start;

    printf("start = %u\tend = %u\tsize = %u\n", start, end, partialSize);
    uint nonzeros = table.rowIndex[end] - table.rowIndex[start];
    printf("nonzeros = %u\n", nonzeros);

    // Initialize every attribute of each csr_args element, then add its size to 
    // the totalSize placeholder. Will be used to realloc() exactly that much memory.
    csr_args[i].table.size = partialSize;
    csr_args[i].table.rowIndex = (uint *) malloc((partialSize+1) * sizeof(uint));
    csr_args[i].table.colIndex = (uint *) malloc(nonzeros * sizeof(uint));
    csr_args[i].table.values = (int *) malloc(nonzeros * sizeof(int));
    csr_args[i].id = i;
    csr_args[i].start = start;
    csr_args[i].end = end;

    uint individualSize = sizeof(int) + 3 * sizeof(uint) +
      (size + 1) * sizeof(uint) + nonzeros * sizeof(uint) + nonzeros * sizeof(int);

    totalSize += individualSize;
  }

  csr_args = (csr_arg *) realloc(csr_args, totalSize);
  printf("\nInitialized csr_args\n\n");

  return csr_args;
}


csr hadamardSingleStep(csr table, uint start, uint end) {
  uint size = end - start;
  
	uint nonzeros = table.rowIndex[end] - table.rowIndex[start];
	uint newNonzeros = 0;
	// The new values array. Intialize all to 0. Do the same for the new column indices.
	int *newValues = (int *) calloc(nonzeros, sizeof(int));
	uint *newColIndex = (uint *) calloc(nonzeros, sizeof(uint));

	// Finally, initialize the new row index array.
	uint *newRowIndex = (uint *) calloc((size+1), sizeof(uint));

  printf("Nonzeros = %u\n", nonzeros);

  // Find the values in A^2, iff the original matrix had a nonzero in that position.
	for (uint row = start; row < end; row++) {
    uint rowStart = table.rowIndex[row];
    uint rowEnd = table.rowIndex[row+1];

    for (uint index = rowStart; index < rowEnd; index++) {
      uint currentColumn = table.colIndex[index];

      int value = dot(table, row, currentColumn);

      if (value > 0) {
        newValues[newNonzeros] = value;
        newColIndex[newNonzeros] = currentColumn;
        newRowIndex[row - start + 1]++;

        newNonzeros++;
      }
    }

    // Pass the next value if we haven't reached the last row.
    if (row < end - 1) {
      newRowIndex[row - start + 2] = newRowIndex[row - start + 1];
    }
	}

  csr hadamard = {size, newValues, newColIndex, newRowIndex};
  return hadamard;
}


// Calculates the dot product of two vectors.
int dot(csr table, uint row, uint column) {
  // Symmetric table. Rows are identical to columns and vice versa.
  uint rowStart = table.rowIndex[row];
  uint rowEnd = table.rowIndex[row+1];

  uint colStart = table.rowIndex[column];
  uint colEnd = table.rowIndex[column+1];

  int value = 0;
  // Leaves a trace of the last match. The next iteratio will start from
  // the nexr position to get rid of unnecessary comparisons.
  uint lastMatch = colStart - 1;

  for (uint i = rowStart; i < rowEnd; i++) {
    for (uint j = lastMatch+1; j < colEnd; j++) {
      if (table.colIndex[i] == table.colIndex[j]) {
        value += table.values[i] * table.values[j];
        break;
      }
    }
  }

  return value;
}


uint countTriangles(csr C) {
	uint size = C.size;
	uint *triangleCount = (uint *) malloc(size * sizeof(uint));
	for (uint i = 0; i < size; i++) {
		triangleCount[i] = 0;
	}

	// Add all the values in each row, then divide by 2.
	// Simulates the operation of multiplying the table with a nx1 vector.
	for (uint i = 0; i < size; i++) {
		uint rowStart = C.rowIndex[i];
		uint rowEnd = C.rowIndex[i+1];

		for (uint j = rowStart; j < rowEnd; j++) {
			triangleCount[i] += C.values[j];
		}
	}

  uint totalTriangles = 0;
	for (uint i = 0; i < size; i++) {
    totalTriangles += triangleCount[i];
  }

  return totalTriangles / 2;
}


// Matrix multiplication. Only need rows1, cols1 and cols2, because
// cols1==rows2 is required. The new matrix is of size rows1 x cols2.
int **matmul (int **table1, int **table2, uint rows1, uint cols1, uint cols2) {
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


// Prints a CSR data structure.
void printCSR(csr converted) {
  uint size = converted.size;
	uint nonzeros = converted.rowIndex[size];
	
	printf("Values:");
	for (int i = 0; i < nonzeros; i++) {
		printf(" %d ", converted.values[i]);
	}

	printf("\nCol_index:");
	for (int i = 0; i < nonzeros; i++) {
		printf(" %u", converted.colIndex[i]);
	}

	printf("\nRow_index:");
	for (int i = 0; i < size+1; i++) {
		printf(" %u ", converted.rowIndex[i]);
	}
	printf("\n\n");
}

#endif