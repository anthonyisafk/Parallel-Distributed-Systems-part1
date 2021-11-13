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
#include "headers/parallel_helpers.h"

int main(int argc, char **argv) {
  FILE *matrixFile;
  int M, N, nz;
  MM_typecode *t;
  char *mtxFileName = "tables/dblp-2010.mtx";

  csr mtx = readmtx_dynamic(mtxFileName, t, N, M, nz);

  struct timeval cilkStop, cilkStart;
  struct timeval serialStop, serialStart;

  // MEASURE OPENCILK IMPLEMENTATION TIME.
  gettimeofday(&cilkStart, NULL);

  uint trianglesCilk = countTrianglesCilk(mtx);

  gettimeofday(&cilkStop, NULL);
  uint cilkTimediff = (cilkStop.tv_sec - cilkStart.tv_sec) * 1000000 + cilkStop.tv_usec - cilkStart.tv_usec;
  printf("\nopenCilk timediff =  %lu us\n", cilkTimediff);

  printf("\nTOTAL TRIANGLES WITH OPENCILK = %u\n", trianglesCilk);


  // MEASURE SERIAL IMPLEMENTATION TIME.
  gettimeofday(&serialStart, NULL);

  csr C = hadamardSingleStep(mtx, 0, mtx.size);
  uint trianglesSerial = countTriangles(C);

  gettimeofday(&serialStop, NULL);
  uint serialTimediff = (serialStop.tv_sec - serialStart.tv_sec) * 1000000 + serialStop.tv_usec - serialStart.tv_usec;
  printf("\nSerial timediff =  %lu us\n", serialTimediff);

    printf("\nTOTAL TRIANGLES WITH SERIAL IMPLEMENTATION = %u\n", trianglesSerial);

}