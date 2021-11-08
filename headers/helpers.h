#ifndef HELPERS_H
#define HELPERS_H


csr matrixToCSR(int **table, long size);
int **CSRtoMatrix(csr table, long size);

void printCSR(csr converted, long size);
void printTable(int **table, long size);

int **matmul (int **table1, int **table2, int size);
int **makeRandomSparseTable(int size);


#endif