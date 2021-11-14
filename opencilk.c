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


uint countTrianglesCilk(csr table) {
  int max_threads = atoi(MAX_THREADS);  
  __cilkrts_set_param("nworkers", MAX_THREADS);
  __cilkrts_init();

  csr_arg *cilk_csr = makeThreadArguments(table, max_threads);

  cilk_for (int i = 0; i < max_threads; i++) {
    printf("cilk_for thread: %d\n", i);
    int index = cilk_csr[i].id;

    cilk_csr[index].table = hadamardSingleStep(table, cilk_csr[index].start, cilk_csr[index].end);
    printf("thread end: %d\n", i);
  }

  printf("\nExited the loop and cilk_sync\n\n");

  // Make a table of the triangles each thread will count for each respective csr_args element.
  // That way, there is no need to stitch the individual CSR structures together.
  uint *trianglesPerThread = (uint *) calloc(max_threads, sizeof(uint));
  uint totalTriangles = 0;

  cilk_for(int i = 0; i < max_threads; i++) {
    int index = cilk_csr[i].id;

    trianglesPerThread[index] = countTriangles(cilk_csr[index].table);
    printf("trianglesPerThread[%d] = %u\n", index, trianglesPerThread[index]);
  }

  for (int i = 0; i < max_threads; i++) {
    totalTriangles += trianglesPerThread[i];
  }

  return totalTriangles;
}


int main(int argc, char **argv) {
  FILE *matrixFile;
  int M, N, nz;
  MM_typecode *t;
  char *mtxFileName = "tables/belgium_osm.mtx";

  csr mtx = readmtx_dynamic(mtxFileName, t, N, M, nz);

  struct timeval cilkStop, cilkStart;
  struct timeval serialStop, serialStart;

  // MEASURE OPENCILK IMPLEMENTATION TIME.
  gettimeofday(&cilkStart, NULL);

  uint triangles_cilk = countTrianglesCilk(mtx);

  gettimeofday(&cilkStop, NULL);
  uint cilkTimediff = (cilkStop.tv_sec - cilkStart.tv_sec) * 1000000 + cilkStop.tv_usec - cilkStart.tv_usec;
  printf("\nopenCilk timediff =  %lu us\n", cilkTimediff);

  printf("\nTOTAL TRIANGLES WITH OPENCILK = %u\n", triangles_cilk);


  // MEASURE SERIAL IMPLEMENTATION TIME.
  gettimeofday(&serialStart, NULL);

  csr C = hadamardSingleStep(mtx, 0, mtx.size);
  uint triangles_serial = countTriangles(C);

  gettimeofday(&serialStop, NULL);
  uint serialTimediff = (serialStop.tv_sec - serialStart.tv_sec) * 1000000 + serialStop.tv_usec - serialStart.tv_usec;
  printf("\nSerial timediff =  %lu us\n", serialTimediff);

  printf("\nTOTAL TRIANGLES WITH SERIAL IMPLEMENTATION = %u\n", triangles_serial);

}