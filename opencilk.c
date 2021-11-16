/*
 * opencilk.c 

 * Convert a square N x N matrix into the CSR format, made for sparse matrices:
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 * Then return its square in CSR format and perform the Hadamard 
 * (element-wise) operation: A (Hadamard) A^2.
 * 
 * -- Parallel implementation of the algorithm using openCilk. --
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
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "headers/csr.h"
#include "headers/csr_arg.h"
#include "headers/mmio.h"
#include "headers/helpers.h"

#define MAX_THREADS "4"


uint *countTrianglesCilk(csr table) {
  int max_threads = atoi(MAX_THREADS);  
  __cilkrts_set_param("nworkers", MAX_THREADS);
  __cilkrts_init();

  uint size = table.size;

  csr_arg *cilk_csr = makeThreadArguments(table, max_threads);

  cilk_for (int i = 0; i < max_threads; i++) {
    printf("cilk_for thread: %d\n", i);
    int index = cilk_csr[i].id;

    cilk_csr[index].table = hadamardSingleStep(table, cilk_csr[index].start, cilk_csr[index].end);
    printf("thread end: %d\n", i);
  }

  printf("\nExited the loop and cilk_sync\n\n");

  // Make a table of the triangles each thread will count for each respective csr_args element.
  uint **trianglesPerThread = (uint **) malloc(max_threads * sizeof(uint *));
  uint totalSize = 0;

  for (int i = 0; i < max_threads; i++) {
    uint partialSize = cilk_csr[i].table.size * sizeof(uint);
    totalSize += partialSize;

    trianglesPerThread[i] = (uint *) malloc(partialSize * sizeof(uint));
  }

  // Realloc the table to avoid going off bounds.
  trianglesPerThread = (uint **) realloc(trianglesPerThread, totalSize);
  
  cilk_for (int i = 0; i < max_threads; i++) {
    int index = cilk_csr[i].id;
    trianglesPerThread[index] = countTriangles(cilk_csr[index].table);
  }

  // Stitch all the board together. Use the 'start' notation to put each
  // element in the designated place. 
  uint *triangles = (uint *) calloc(size, sizeof(uint));
  for (int i = 0; i < max_threads; i++) {
    uint threadStart = cilk_csr[i].start;
    uint threadSize = cilk_csr[i].table.size;

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

  struct timeval cilkStop, cilkStart;
  struct timeval serialStop, serialStart;

  // MEASURE OPENCILK IMPLEMENTATION TIME.
  gettimeofday(&cilkStart, NULL);

  uint *triangles_cilk = countTrianglesCilk(mtx);

  gettimeofday(&cilkStop, NULL);
  uint cilkTimediff = (cilkStop.tv_sec - cilkStart.tv_sec) * 1000000 + cilkStop.tv_usec - cilkStart.tv_usec;

  printf("\nopenCilk timediff =  %u us\n", cilkTimediff);

}