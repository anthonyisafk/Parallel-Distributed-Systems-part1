/*
 * csr_arg.h
 * Contains all the necessary info for a CSR table to be cut into pieces (of rows)
 * each of which is given to a thread during the parallel implementation of the algorithm.
 * 
 * @param table: The original CSR structure.
 * @param start: The first row of the sub-table.
 * @param end: The last row of the sub-table.
 * @param id: Given to the table to avoid memory conflicts during loops.
 */

#ifndef CSR_ARG_H
#define CSR_ARG_H

#include <stdio.h>
#include "csr.h"

typedef struct {
  csr table;
  int id;
  uint start;
  uint end;
} csr_arg;


#endif