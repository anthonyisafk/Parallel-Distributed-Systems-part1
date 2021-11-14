#ifndef HELPERS_H
#define HELPERS_H

// Final version of the functions used.
csr readmtx_dynamic(char *mtx, MM_typecode *t, int N, int M, int nz);
csr_arg *makeThreadArguments(csr table, int max_threads);
csr hadamardSingleStep(csr table, uint start, uint end);
int dot(csr table, uint row, uint column);
uint countTriangles(csr table);
void printCSR(csr converted);

// Void functions used for the pthread implementation.
void *hadamardSingleStepVoid(void *csrarg);
void *countTrianglesVoid(void *table);

#endif