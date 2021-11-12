#ifndef HELPERS_H
#define HELPERS_H

// Final version of the functions used.
csr readmtx_dynamic(char *mtx, MM_typecode *t, int N, int M, int nz);
csr hadamardSingleStep(csr table, uint start, uint end);
uint *countTriangles(csr table);
int dot(csr table, uint row, uint column);
int **matmul (int **table1, int **table2, uint rows1, uint cols1, uint cols2);

// Final versions of the testing functions.
csr csrSquare(csr table, uint size);
csr newhadamard(csr csrTable, csr square, uint size);

// Original, error-prone versions.
csr readmtx(char *mtx, MM_typecode *t, int N, int M, int nz);
csr csrSquareAlt(csr converted, int **table, uint size);
csr hadamard(csr csrTable, int **square, uint size);

// General transforming and printing  helpers.
csr matrixToCSR(int **table, uint size);
int **CSRtoMatrix(csr table, uint size);
void printCSR(csr converted);
void printTable(int **table, uint size);
int **makeRandomSparseTable(int size);

#endif