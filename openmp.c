/*
 * openmp.c 
 *
 * Convert a square N x N matrix into the CSR format, made for sparse matrices:
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 * Then return its square in CSR format and perform the Hadamard 
 * (element-wise) operation: A (Hadamard) A^2.
 * 
 * -- Parallel implementation of the algorithm using openMP. --
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
#include <omp.h>

#include "headers/csr.h"
#include "headers/csr_arg.h"
#include "headers/mmio.h"
#include "headers/helpers.h"

#define MAX_THREADS 4


uint *countTrianglesOMP(csr table) {
  uint size = table.size;

  omp_set_num_threads(MAX_THREADS);
  csr_arg *omp_csr = makeThreadArguments(table, MAX_THREADS);

  #pragma omp parallel 
  {
    int id = omp_get_thread_num();
    printf("thread: %d\n", id);
    omp_csr[id].table = hadamardSingleStep(table, omp_csr[id].start, omp_csr[id].end);
    printf("thread end: %d\n", id);
  }

 // Make a table of the triangles each thread will count for each respective csr_args element.
  uint **trianglesPerThread = (uint **) malloc(MAX_THREADS * sizeof(uint *));
  uint totalSize = 0;

  for (int i = 0; i < MAX_THREADS; i++) {
    uint partialSize = omp_csr[i].table.size * sizeof(uint);
    totalSize += partialSize;

    trianglesPerThread[i] = (uint *) malloc(partialSize * sizeof(uint));
  }

  // Realloc the table to avoid going off bounds.
  trianglesPerThread = (uint **) realloc(trianglesPerThread, totalSize);

  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    trianglesPerThread[id] = countTriangles(omp_csr[id].table);
  }

  // Stitch all the board together. Use the 'start' notation to put each
  // element in the designated place. 
  uint *triangles = (uint *) calloc(size, sizeof(uint));
  for (int i = 0; i < MAX_THREADS; i++) {
    uint threadStart = omp_csr[i].start;
    uint threadSize = omp_csr[i].table.size;

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

  struct timeval openmpStop, openmpStart;

  // MEASURE OPENCILK IMPLEMENTATION TIME.
  gettimeofday(&openmpStart, NULL);

  uint *triangles_omp = countTrianglesOMP(mtx);

  gettimeofday(&openmpStop, NULL);
  uint cilkTimediff = (openmpStop.tv_sec - openmpStart.tv_sec) * 1000000 + openmpStop.tv_usec - openmpStart.tv_usec;
  printf("\nopenMP timediff =  %lu us\n", cilkTimediff);

}