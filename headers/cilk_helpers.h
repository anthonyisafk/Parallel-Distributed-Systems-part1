#ifndef CILK_HELPERS_H
#define CILK_HELPERS_H

csr csrSquareCilk(csr table, long size);
csr newhadamardCilk(csr csrTable, csr square, long size);
long *countTrianglesCilk(csr C);

#endif