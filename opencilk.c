#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <unistd.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <time.h>

#include "headers/csr.h"
#include "headers/mmio.h"
#include "headers/helpers.h"
#include "headers/cilk_helpers.h"

int main(int argc, char **argv) {
  
  FILE *matrixFile;
  int M, N, nz;
  MM_typecode t;
  char *mtxFileName = "tables/G51.mtx";

  csr mtx = readmtx(mtxFileName, t, N, M, nz);
  __cilkrts_end_cilk();
  __cilkrts_set_param("nworkers", "4");
  __cilkrts_init();

  struct timeval cilkStop, cilkStart;
  struct timeval serialStop, serialStart;

  // MEASURE OPENCILK IMPLEMENTATION TIME.
  gettimeofday(&cilkStart, NULL);

  csr square = csrSquareCilk(mtx, mtx.size);

  gettimeofday(&cilkStop, NULL);
  uint timediff = (cilkStop.tv_sec - cilkStart.tv_sec) * 1000000 + cilkStop.tv_usec - cilkStart.tv_usec;
  printf("\nopenCilk timediff =  %lu us\n", timediff);


  // MEASURE SERIAL IMPLEMENTATION TIME.
  gettimeofday(&serialStart, NULL);

  csrSquare(mtx, mtx.size);

  gettimeofday(&serialStop, NULL);
  timediff = (serialStop.tv_sec - serialStart.tv_sec) * 1000000 + serialStop.tv_usec - serialStart.tv_usec;
  printf("\nSerial timediff =  %lu us\n", timediff);


}