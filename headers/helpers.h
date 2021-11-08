#ifndef HELPERS_H
#define HELPERS_H

// csr readmtx(char *mtx, MM_typecode *t, int N, int M, int nz);
csr csrSquare(csr table, long size);
csr csrSquareAlt(csr converted, int **table, long size);
csr newhadamard(csr csrTable, csr square, long size);
csr hadamard(csr csrTable, int **square, long size);

csr matrixToCSR(int **table, long size);
int **CSRtoMatrix(csr table, long size);

void printCSR(csr converted, long size);
void printTable(int **table, long size);

int **matmul (int **table1, int **table2, int size);
int **makeRandomSparseTable(int size);


#endif