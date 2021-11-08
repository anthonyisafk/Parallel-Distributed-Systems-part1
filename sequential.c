/*
 * sequential.c 
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
#include <math.h>
#include <limits.h>

#include "headers/csr.h"
#include "headers/csr_operations.h"
#include "headers/mmio.h"
#include "headers/helpers.h"

#define SIZE 6 

int main(int argc, char **argv) {
	int M, N, nz;
	MM_typecode *t;
	char *mtx = "tables/mycielskian4.mtx";

	csr mtx_csr = readmtx(mtx, t, N, M, nz);
	// printCSR(mtx_csr, N);

  int **test = makeRandomSparseTable(SIZE);
  printTable(test, SIZE);

	return 0;
}