#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define SIZE 15

// A struct used to turn sparse matrices to CSR data structures.
// The notation and algorithm used is taken directly from the given Wikipedia page:
// https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
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

  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      table[i][j] = (rand() % 100 > 13) ? 0 : 1;
    }
  }

  return table;
}


void printTable(int **table) {
  for (int i = 0; i < SIZE; i++) {
    for (int j = 0; j < SIZE; j++) {
      printf(" %d ", table[i][j]);
    }
    printf("\n");
  }

  printf("\n");
}


// Converts a sparse matrix to CSR.
csr convertTOCSR(int **table, long size) {
  // Keep track of the nonzero objects.
  long nonzeros; 
  // Initialize the 3 necessary arrays.
  int *values = (int *) malloc(100 * sizeof(int));
  long *colIndex = (long *) malloc(100 * sizeof(long));
  long *rowIndex = (long *) malloc((size+1) * sizeof(long));
  
  for (long i = 0; i < size; i++) {
    // Add the nonzero values that are placed above the current row. 
    rowIndex[i] = (long) nonzeros;

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
        colIndex[nonzeros-1] = (long) j; 
      }
    }

    // Add the last value before exiting the for loop.
    if (i == size - 1) {
      rowIndex[size] = (long) nonzeros;
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

// int **matmulParallel(int **table1, int **table2, int size) {
//   continue;
// }

int main(int argc, char **argv) {

  int **random1 = makeRandomSparseTable(SIZE);
  int **random2 = makeRandomSparseTable(SIZE);
  printTable(matmul(random1, random1, SIZE));

  // csr converted = convertTOCSR(random, (long)SIZE);



  

  return 0;
}