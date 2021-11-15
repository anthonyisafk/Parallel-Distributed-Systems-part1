/*
 * openmp.c 

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


uint countTrianglesOMP(csr table) {
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
  // That way, there is no need to stitch the individual CSR structures together.
  uint *trianglesPerThread = (uint *) calloc(MAX_THREADS, sizeof(uint));
  uint totalTriangles = 0;

  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    trianglesPerThread[id] = countTriangles(omp_csr[id].table);
  }

  for (int i = 0; i < MAX_THREADS; i++) {
    totalTriangles += trianglesPerThread[i];
  }

  return totalTriangles;
}


int main(int argc, char **argv) {
  FILE *matrixFile;
  int M, N, nz;
  MM_typecode *t;
  char *mtxFileName = "tables/NACA0015.mtx";

  csr mtx = readmtx_dynamic(mtxFileName, t, N, M, nz);

  struct timeval openmpStop, openmpStart;
  struct timeval serialStop, serialStart;

  // MEASURE OPENCILK IMPLEMENTATION TIME.
  gettimeofday(&openmpStart, NULL);

  uint triangles_omp = countTrianglesOMP(mtx);

  gettimeofday(&openmpStop, NULL);
  uint cilkTimediff = (openmpStop.tv_sec - openmpStart.tv_sec) * 1000000 + openmpStop.tv_usec - openmpStart.tv_usec;
  printf("\nopenMP timediff =  %lu us\n", cilkTimediff);
  printf("\nTOTAL TRIANGLES WITH OPENMP = %u\n", triangles_omp);

  // MEASURE SERIAL IMPLEMENTATION TIME.
  gettimeofday(&serialStart, NULL);

  csr C = hadamardSingleStep(mtx, 0, mtx.size);
  uint trianglesSerial = countTriangles(C);

  gettimeofday(&serialStop, NULL);
  uint serialTimediff = (serialStop.tv_sec - serialStart.tv_sec) * 1000000 + serialStop.tv_usec - serialStart.tv_usec;
  
  printf("\nSerial timediff =  %lu us\n", serialTimediff);
  printf("\nTOTAL TRIANGLES WITH SERIAL IMPLEMENTATION = %u\n", trianglesSerial);

}