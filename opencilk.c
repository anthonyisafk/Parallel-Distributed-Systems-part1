/*
 * opencilk.c 

 * Convert a square N x N matrix into the CSR format, made for sparse matrices:
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 * Then return its square in CSR format and perform the Hadamard 
 * (element-wise) operation: A (Hadamard) A^2.
 * 
 * -- Parallel implementation of the algorithm using openCilk. --
 * 
 * Authors: Antonios Antoniou - 9482 - aantonii@ece.auth.gr
 *          Efthymios Grigorakis - 9694 - eegrigor@ece.auth.gr
 *    
 * 2021 Aristotle University of Thessaloniki
 * Parallel and Distributed Systems.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "headers/csr.h"
#include "headers/csr_arg.h"
#include "headers/mmio.h"
#include "headers/helpers.h"
#include "headers/data_arg.h"


uint *countTrianglesCilk(csr table, char *MAX_THREADS) {
  int max_threads = atoi(MAX_THREADS);
  printf("Threads: %d\n", max_threads);  
  __cilkrts_set_param("nworkers", MAX_THREADS);
  __cilkrts_init();

  uint size = table.size;

  csr_arg *cilk_csr = makeThreadArguments(table, max_threads);

  cilk_for (int i = 0; i < max_threads; i++) {
    printf("cilk_for thread: %d\n", i);
    int index = cilk_csr[i].id;

    cilk_csr[index].table = hadamardSingleStep(table, cilk_csr[index].start, cilk_csr[index].end);
    printf("thread end: %d\n", i);
  }

  // Make a table of the triangles each thread will count for each respective csr_args element.
  uint **trianglesPerThread = (uint **) malloc(max_threads * sizeof(uint *));
  uint totalSize = 0;

  for (int i = 0; i < max_threads; i++) {
    uint partialSize = cilk_csr[i].table.size * sizeof(uint);
    totalSize += partialSize;

    trianglesPerThread[i] = (uint *) malloc(partialSize * sizeof(uint));
  }

  // Realloc the table to avoid going off bounds.
  trianglesPerThread = (uint **) realloc(trianglesPerThread, totalSize);
  
  cilk_for (int i = 0; i < max_threads; i++) {
    int index = cilk_csr[i].id;
    trianglesPerThread[index] = countTriangles(cilk_csr[index].table);
  }

  // Stitch all the board together. Use the 'start' notation to put each
  // element in the designated place. 
  uint *triangles = (uint *) calloc(size, sizeof(uint));
  for (int i = 0; i < max_threads; i++) {
    uint threadStart = cilk_csr[i].start;
    uint threadSize = cilk_csr[i].table.size;

    for (uint j = 0; j < threadSize; j++) {
      triangles[threadStart + j] = trianglesPerThread[i][j];
    }
  }

  // for (uint i = 0; i < size; i++) {
  //   printf(" %u ", triangles[i]);
  // }

  return triangles;
}


data_arg measureTimeCilk(csr mtx, char *filename, MM_typecode *t, int N, int M, int nz, char *MAX_THREADS) {
  struct timeval stop, start;

  gettimeofday(&start, NULL);
  uint *triangles = countTrianglesCilk(mtx, MAX_THREADS);
  gettimeofday(&stop, NULL);
  
  uint timediff = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;

  printf("\openCilk took %u us for file %s, using %s threads.\n\n", 
    timediff, filename, MAX_THREADS);

  data_arg data = {timediff, triangles};
  return data;
}

int main(int argc, char **argv) {
  FILE *matrixFile;
  int M, N, nz;

  MM_typecode *t;
  char *filenames[5] = {
    "tables/belgium_osm.mtx",
    "tables/dblp-2010.mtx",
    "tables/NACA0015.mtx",
    "tables/mycielskian13.mtx",
    "tables/com-Youtube.mtx"
  };

  char *num_threads[3] = {"2", "4", "8"};

  int reps = 12;
  int files_num = 5;
  int thread_num = 3;

  // Select the number of threads.
  int thread_index = atoi(argv[1]);
  // Select the matrix to be read.
  int file = atoi(argv[2]);

  // Set the title, depending on the number of threads selected.
  // Then print it, if this is the first file for that number.
  char cilk[6] = "cilkN";
  cilk[4] = num_threads[thread_index][0];

  FILE *statsFile = fopen("stats/data.csv", "a");
  if (file == 0) {
    fprintf(statsFile, "\n%s", cilk);
  }

  csr mtx = readmtx_dynamic(filenames[file], t, N, M, nz);
  csr *mtx_addr = &mtx;

  unsigned long totalTime = 0;
  unsigned long *time_addr = &totalTime;

  for (int rep = 0; rep < reps; rep++) {
    data_arg data = measureTimeCilk(mtx, filenames[file], t, N, M, nz, num_threads[thread_index]);
    data_arg *data_addr = &data;

    if (rep > 1) {
      totalTime += data.time;
    }
  }

  uint meanTime = totalTime / (reps - 2);
  remove(time_addr);
  remove(mtx_addr);

  fprintf(statsFile, "\t%u", meanTime); 
  fclose(statsFile);

  return 0;
}