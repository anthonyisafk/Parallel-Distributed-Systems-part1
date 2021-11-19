/*
 * sequential.c 
 *
 * Convert a square N x N matrix into the CSR format, made for sparse matrices:
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 * Then return its square in CSR format and perform the Hadamard 
 * (element-wise) operation: A (Hadamard) A^2.
 * 
 * -- Serial version of the operations. --
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
#include "headers/data_arg.h"

data_arg measureTimeSerial(csr mtx, char *filename) {  
  struct timeval stop, start;
  gettimeofday(&start, NULL);

  csr C = hadamardSingleStep(mtx, 0, mtx.size);
  uint *triangles = countTriangles(C);

  gettimeofday(&stop, NULL);
  uint timediff = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;
  
  printf("\n\nThe serial algorithm took %lu us for %s\n", timediff, filename);

  data_arg data = {timediff, triangles};
  return data;
}


int main(int argc, char **argv) {
  int M, N, nz;
  MM_typecode *t;

  char *filenames[5] = {
    "tables/belgium_osm.mtx",
    "tables/dblp-2010.mtx",
    "tables/NACA0015.mtx",
    "tables/mycielskian13.mtx",
    "tables/com-Youtube.mtx"
  };

  // Used for the CSV file.
  char *names[5] = {
    "belgium_osm",
    "dblp2010",
    "NACA0015",
    "mycielskian13",
    "comYoutube"
  };
  
  // The index of the file to be read.
  int file = atoi(argv[1]); 
  int files_num = 5;
  int reps = 12;

  // If we're reading the first matrix, add a title to the CSV file.
  FILE *statsFile = fopen("stats/data.csv", "a");
  if (file == 0) {
    char *title = "library_threads";

    fprintf(statsFile, "%s\t", title);

    for (int i = 0; i < files_num; i++) {
      // Don't add a delimiter if we reached the end of the line.
      // Change line instead.
      if (i == files_num - 1) {
        fprintf(statsFile, "%s\n", names[i]);
      } else {
        fprintf(statsFile, "%s\t", names[i]);
      }
    }
    char *row1 = "seq";
    fprintf(statsFile, "%s\t", row1);
  }

  csr mtx = readmtx_dynamic(filenames[file], t, N, M, nz);

  uint totalTime = 0;
  uint *time_addr = &totalTime;

  for (int rep = 0; rep < reps; rep++) {
    data_arg data = measureTimeSerial(mtx, filenames[file]);
    data_arg *data_addr = &data;
    if (rep > 1) {
      totalTime += data.time;
    }
    remove(data_addr);
  }

  uint meanTime = totalTime / (reps - 2);
  remove(time_addr);

  fprintf(statsFile, "\t%d", meanTime);

  fclose(statsFile);
  return 0;
}
