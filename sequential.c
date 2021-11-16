/*
 * sequential.c 
 *
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

#include "headers/csr.h"
#include "headers/mmio.h"
#include "headers/helpers.h"
#include "headers/csr_arg.h"

uint *measureTimeSerial(char *filename, MM_typecode *t, int N, int M, int nz) {
  // The data array that is returned. First element is the implementation time
  // Seconf element is the actual result, used to ensure everything worked properly.
  // uint data[2] = {0, 0};
  
  csr mtx = readmtx_dynamic(filename, t, N, M, nz);

  // CALCULATE THE C ARRAY AND COUNT TRIANGLES.
  struct timeval stop, start;
  gettimeofday(&start, NULL);

  csr C = hadamardSingleStep(mtx, 0, mtx.size);
  uint *data = countTriangles(C);

  gettimeofday(&stop, NULL);
  // data[0] = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
  
  printf("\n\nThe serial algorithm took %lu us\n", data[0]);
  printf("\nTriangles in total: %ld\n", data[1]);

  for(uint i = 0; i < mtx.size; i++) {
    printf(" %u ", data[i]);
  }

  return data;
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

  measureTimeSerial(filenames[6], t, N, M, nz);


  return 0;
}