/*
 * pthreads.c 
 *
 * Convert a square N x N matrix into the CSR format, made for sparse matrices:
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 * Then return its square in CSR format and perform the Hadamard 
 * (element-wise) operation: A (Hadamard) A^2.
 * 
 * -- Parallel implementation of the algorithm using pthreads. --
 * 
 * Authors: Antonios Antoniou - 9482 - aantonii@ece.auth.gr
 *          Efthymios Grigorakis - 9694 - eegrigor@ece.auth.gr
 *    
 * 2021 Aristotle University of Thessaloniki
 * Parallel and Distributed Systems.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#include "headers/csr.h"
#include "headers/csr_arg.h"
#include "headers/mmio.h"
#include "headers/helpers.h"

#define MAX_THREADS 4

// Utilized to be given to the void functions used by pthread. 
// @param csrarg: the csr_arg struct that's used in the other two parallel algorithms.
// @param original: contains the full, original CSR table.
// @param triangles: the result of the algorithm.
typedef struct {
  csr original;
  csr_arg csrarg;
  uint *triangles; 
} pthreads_arg;


void *countTrianglesVoid(void *C) {
  pthreads_arg *C_arg = (pthreads_arg *) C;

  uint size = C_arg->csrarg.table.size;

	// Add all the values in each row, then divide by 2.
	// Simulates the operation of multiplying the table with a nx1 vector.
	for (uint i = 0; i < size; i++) {
    C_arg->triangles[i] = 0;
		uint rowStart = C_arg->csrarg.table.rowIndex[i];
		uint rowEnd = C_arg->csrarg.table.rowIndex[i+1];

		for (uint j = rowStart; j < rowEnd; j++) {
      C_arg->triangles[i] = (j == rowEnd - 1) 
        ? (C_arg->triangles[i] + C_arg->csrarg.table.values[j]) / 2
        : C_arg->triangles[i] + C_arg->csrarg.table.values[j];
		}
	}
}


void *hadamardSingleStepVoid(void *csrarg) {
  pthreads_arg *arg = (pthreads_arg *) csrarg;

  uint start = arg->csrarg.start;
  uint end = arg->csrarg.end;
  uint size = end - start;
  
	uint nonzeros = arg->original.rowIndex[end] - arg->original.rowIndex[start];
	uint newNonzeros = 0;
	// The new values array. Intialize all to 0. Do the same for the new column indices.
	int *newValues = (int *) calloc(nonzeros, sizeof(int));
	uint *newColIndex = (uint *) calloc(nonzeros, sizeof(uint));

	// Finally, initialize the new row index array.
	uint *newRowIndex = (uint *) calloc((size+1), sizeof(uint));

  // Find the values in A^2, iff the original matrix had a nonzero in that position.
	for (uint row = start; row < end; row++) {
    uint rowStart = arg->original.rowIndex[row];
    uint rowEnd = arg->original.rowIndex[row+1];

    for (uint index = rowStart; index < rowEnd; index++) {
      uint currentColumn = arg->original.colIndex[index];

      int value = dot(arg->original, row, currentColumn);

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

  arg->csrarg.table.size = size;
  arg->csrarg.table.values = newValues;
  arg->csrarg.table.rowIndex = newRowIndex;
  arg->csrarg.table.colIndex = newColIndex;
}


uint countTrianglesPthread(csr table) {
  uint size = table.size;

  csr_arg *pthread_csr = makeThreadArguments(table, MAX_THREADS);
  pthread_t threads[MAX_THREADS];

  // Initialize the sub-arrays and give them to the new pthreads_arg struct.
  pthreads_arg arg[MAX_THREADS]; 

  for (int i = 0 ; i < MAX_THREADS ; i++) {
    arg[i].triangles = (uint *) calloc(pthread_csr[i].table.size, sizeof(uint));
    arg[i].original = table ; 
    arg[i].csrarg = pthread_csr[i];
  }

  for (int i = 0; i < MAX_THREADS; i++) {
    printf("thread: %d\n", i);
    pthread_create(&threads[i], NULL, hadamardSingleStepVoid, (void *) &arg[i]);
  }

  for(int i = 0 ; i < MAX_THREADS; i++) {
    pthread_join(threads[i], NULL);
    printf("thread end: %d\n", i);
  }

 // Make a table of the triangles each thread will count for each respective csr_args element.
  uint **trianglesPerThread = (uint **) malloc(MAX_THREADS * sizeof(uint *));
  uint totalSize = 0;

  for (int i = 0; i < MAX_THREADS; i++) {
    uint partialSize = arg[i].csrarg.table.size * sizeof(uint);
    totalSize += partialSize;

    trianglesPerThread[i] = (uint *) calloc(partialSize, sizeof(uint));
  }

  // Realloc the table to avoid going off bounds.
  trianglesPerThread = (uint **) realloc(trianglesPerThread, totalSize);

  for (int i = 0; i < MAX_THREADS; i++) {
    pthread_create(&threads[i], NULL, countTrianglesVoid, (void *) &arg[i]);
  }

  for (int i = 0; i < MAX_THREADS; i++) {
    pthread_join(threads[i], NULL);
  }

  for (int i = 0; i < MAX_THREADS; i++) {
    uint threadStart = arg[i].csrarg.start;
    uint threadSize = arg[i].csrarg.table.size;

    for (int j = 0; j < threadSize; j++) {
      trianglesPerThread[i][j] = arg[i].triangles[j]; 
    }
  }

  // Stitch all the board together. Use the 'start' notation to put each
  // element in the designated place. 
  uint *triangles = (uint *) calloc(size, sizeof(uint));
  for (int i = 0; i < MAX_THREADS; i++) {
    uint threadStart = arg[i].csrarg.start;
    uint threadSize = arg[i].csrarg.table.size;

    for (uint j = 0; j < threadSize; j++) {
      triangles[threadStart + j] = trianglesPerThread[i][j];
    }
  }

  for(uint i = 0; i < size; i++) {
    printf(" %u ", triangles[i]);
  }

  return triangles;
}


int main(int argc, char **argv) {
  FILE *matrixFile;
  int M, N, nz;
  MM_typecode *t;
  char *filenames[7] = {
    "tables/belgium_osm.mtx",
    "tables/dblp-2010.mtx",
    "tables/NACA0015.mtx",
    "tables/mycielskian13.mtx",
    "tables/com-Youtube.mtx",
    "tables/ca-CondMat.mtx",
    "tables/karate.mtx"
  };

  csr mtx = readmtx_dynamic(filenames[6], t, N, M, nz);

  struct timeval pthreadStop, pthreadStart;

  // MEASURE OPENCILK IMPLEMENTATION TIME.
  gettimeofday(&pthreadStart, NULL);

  uint *triangles_pthread = countTrianglesPthread(mtx);

  gettimeofday(&pthreadStop, NULL);
  uint pthreadTimediff = (pthreadStop.tv_sec - pthreadStart.tv_sec) * 1000000 + pthreadStop.tv_usec - pthreadStart.tv_usec;
  
  printf("\npthread timediff =  %lu us\n", pthreadTimediff);

}