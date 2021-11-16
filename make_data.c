/**
 * make_data.c 
 *
 * Call all the algorithms realized in the project.
 * Run each file, for each algorithm (and number of threads for the parallel implementations),
 * for 10 times to get the average time.
 * 
 * Then write the data to a CSV file which will be used to 
 * visualize the data with Julia.
 * 
 * Authors: Antonios Antoniou - 9482
 *          Efthymios Grigorakis - 9694
 *    
 * 2021 Aristotle University of Thessaloniki
 * Parallel and Distributed Systems.
 */

#include <stdio.h>
#include <stdlib.h>

#include "headers/csr.h"
#include "headers/mmio.h"
#include "headers/helpers.h"
#include "headers/csr_arg.h"

#include "sequential.c"


int main(int argc, char **argv) {
  FILE *matrixFile;
  int M, N, nz;
  MM_typecode *t;
  char *filenames[2] = {"tables/NACA0015.mtx", "tables/dictionary28.mtx"};

  for (int i = 0; i < 2; i++) {
    uint **data = measureTimeSerial(filenames[i], t, N, M, nz);
  }

  return 0;
}
