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