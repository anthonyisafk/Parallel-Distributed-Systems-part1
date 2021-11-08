#ifndef CSR_OPERATIONS_H
#define CSR_OPERATIONS_H

csr readmtx(char *mtx, MM_typecode *t, int N, int M, int nz);
csr csrSquare(csr table, long size);
csr csrSquareAlt(csr converted, int **table, long size);
csr newhadamard(csr csrTable, csr square, long size);
csr hadamard(csr csrTable, int **square, long size);


#endif