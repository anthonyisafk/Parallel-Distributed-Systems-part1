#ifndef HELPERS_H
#define HELPERS_H

#include "csr_arg.h"

// Final version of the functions used.
csr readmtx_dynamic(char *mtx, MM_typecode *t, int N, int M, int nz);
csr_arg *makeThreadArguments(csr table, int max_threads);
csr hadamardSingleStep(csr table, uint start, uint end);
int dot(csr table, uint row, uint column);
uint *countTriangles(csr C);
void printCSR(csr converted);

#endif