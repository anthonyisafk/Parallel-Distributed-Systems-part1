/*
 * progress_copies.h
 * Copies we kept of functions that were used for debugging the algorithm, 
 * or versions of functions that were a stepping stone to the final versions
 * of other functions or the algorithm in general.
 * 
 */

#ifndef PROGRESS_COPIES_H
#define PROGRESS_COPIES_H

csr csrSquare(csr table, uint size);
csr newhadamard(csr csrTable, csr square, uint size);
int **matmul (int **table1, int **table2, uint rows1, uint cols1, uint cols2);

csr readmtx(char *mtx, MM_typecode *t, int N, int M, int nz);
csr csrSquareAlt(csr converted, int **table, uint size);
csr hadamard(csr csrTable, int **square, uint size);

// General transforming and printing  helpers.
csr matrixToCSR(int **table, uint size);
int **CSRtoMatrix(csr table, uint size);
void printTable(int **table, uint size);
int **makeRandomSparseTable(int size);

#endif