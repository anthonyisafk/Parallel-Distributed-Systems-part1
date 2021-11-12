/*
 * sequential.c 

 * Convert a square N x N matrix into the CSR format, made for sparse matrices:
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 * Then return its square in CSR format and perform the Hadamard 
 * (element-wise) operation: A (Hadamard) A^2.
 * 
 * -- Serial version of the operations. --
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
#include <limits.h>
#include <unistd.h>

#include "headers/csr.h"
#include "headers/mmio.h"
#include "headers/helpers.h"


int main(int argc, char **argv) {
  FILE *matrixFile;
  int M, N, nz;
  MM_typecode t;
  char *mtxFileName = "tables/dblp-2010.mtx";
  csr mtx = readmtx_dynamic(mtxFileName, t, N, M, nz);
  printf("\n\nJust finished reading.\n\n");

  // CALCULATE THE C ARRAY AND COUNT TRIANGLES.
  struct timeval stop, start;
  gettimeofday(&start, NULL);

  csr C = hadamardSingleStep(mtx, 0, mtx.size);
  uint *newTriangles = countTriangles(C);

  gettimeofday(&stop, NULL);
  uint timediff = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
  printf("\nThe serial algorithm took %lu us\n\n", timediff);

  uint triangles = 0;
  for (int i = 0; i < C.size; i++) {
    triangles += newTriangles[i];
  }

  printf("\nTriangles in total: %ld\n", triangles);

}