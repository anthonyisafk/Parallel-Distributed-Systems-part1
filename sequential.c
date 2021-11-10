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

#define SIZE 6


int main(int argc, char **argv) {
  // int **test = makeRandomSparseTable(SIZE);
  // printTable(test, SIZE);
  // int **square = matmul(test, test, SIZE, SIZE, SIZE);
  // printTable(square, SIZE);

  // csr testCSR = matrixToCSR(test, SIZE);
  // csr testSquare = csrSquare(testCSR, SIZE);
  // printTable(CSRtoMatrix(testSquare, SIZE), SIZE);
  // printCSR(&testSquare, SIZE);

  // csr C = newhadamard(testCSR, testSquare, SIZE);
  // printCSR(&C, SIZE);

  FILE *matrixFile;
  int M, N, nz;
  MM_typecode t;
  char *mtxFileName = "tables/G51.mtx";

  csr mtx = readmtx(mtxFileName, t, N, M, nz);
  // printCSR(mtx, mtx.size);

  csr mtx_square = csrSquare(mtx, mtx.size);
  // printCSR(mtx_square, mtx_square.size);

  csr C = newhadamard(mtx, mtx_square, mtx.size);
  // printCSR(C, C.size);

  long *triangles = countTriangles(C);
  for (long i = 0; i < C.size; i++) {
    printf(" %ld ", triangles[i]);
  }
  printf("\n");


}