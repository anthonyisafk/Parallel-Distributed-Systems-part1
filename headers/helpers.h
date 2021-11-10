#ifndef HELPERS_H
#define HELPERS_H

// Final version of the functions used.
csr readmtx(char *mtx, MM_typecode *t, int N, int M, int nz);
csr csrSquare(csr table, long size);
csr newhadamard(csr csrTable, csr square, long size);
long *countTriangles(csr table);
int **matmul (int **table1, int **table2, long rows1, long cols1, long cols2);

// Original, error-prone versions.
csr csrSquareAlt(csr converted, int **table, long size);
csr hadamard(csr csrTable, int **square, long size);

// General transforming and printing  helpers.
csr matrixToCSR(int **table, long size);
int **CSRtoMatrix(csr table, long size);
void printCSR(csr converted, long size);
void printTable(int **table, long size);
int **makeRandomSparseTable(int size);


#endif