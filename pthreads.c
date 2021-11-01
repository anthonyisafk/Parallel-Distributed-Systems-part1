/*
 * pthreads.c 
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
#include <pthread.h>
#include <math.h>

#define SIZE 15

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

  // Add a 1 to the table with a 136/1000 probability.
  // Only add values to the upper triangle of the table.
  for (int i = 0; i < size; i++) {
    for (int j = i; j < size; j++) {
      table[i][j] = (rand() % 1000 > 136) ? 0 : 1;
    }
  }

  // Make it symmetric by setting A[i][j] = A[j][i]
  // for the rest of the values.
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < i; j++) {
      table[i][j] = table[j][i];
    }
  }

  return table;
}


void printCSR(csr converted, long size) {
  long nonzeros = converted.rowIndex[SIZE];
  
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
  printf("\n");
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
        if (nonzeros % 100 == 0) {
          printf("%d\n", table[i][j]);

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


int main(int argc, char **argv) {

  // Initialize a random sparse matrix.
  int **random1 = makeRandomSparseTable(SIZE);
  printTable(random1, SIZE);

  // Try to convert the table to CSR. (Finally works)
  csr converted = matrixToCSR(random1, (long)SIZE);
  
  // Print the resulting arrays, to test the result.
  printCSR(converted, SIZE);

  
  return 0;
}