/*
 * openmp.c 
 *
 * Convert a square N x N matrix into the CSR format, made for sparse matrices:
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 * Then return its square in CSR format and perform the Hadamard 
 * (element-wise) operation: A (Hadamard) A^2.
 * 
 * -- Parallel implementation of the algorithm using openMP. --
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
#include <omp.h>

#include "headers/csr.h"
#include "headers/csr_arg.h"
#include "headers/mmio.h"
#include "headers/helpers.h"
#include "headers/data_arg.h"


uint *countTrianglesOMP(csr table, int MAX_THREADS) {
  uint size = table.size;

  omp_set_num_threads(MAX_THREADS);
  csr_arg *omp_csr = makeThreadArguments(table, MAX_THREADS);

  #pragma omp parallel 
  {
    int id = omp_get_thread_num();
    printf("thread: %d\n", id);
    omp_csr[id].table = hadamardSingleStep(table, omp_csr[id].start, omp_csr[id].end);
    printf("thread end: %d\n", id);
  }

 // Make a table of the triangles each thread will count for each respective csr_args element.
  uint **trianglesPerThread = (uint **) malloc(MAX_THREADS * sizeof(uint *));
  uint totalSize = 0;

  for (int i = 0; i < MAX_THREADS; i++) {
    uint partialSize = omp_csr[i].table.size * sizeof(uint);
    totalSize += partialSize;

    trianglesPerThread[i] = (uint *) malloc(partialSize * sizeof(uint));
  }

  // Realloc the table to avoid going off bounds.
  trianglesPerThread = (uint **) realloc(trianglesPerThread, totalSize);

  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    trianglesPerThread[id] = countTriangles(omp_csr[id].table);
  }

  // Stitch all the board together. Use the 'start' notation to put each
  // element in the designated place. 
  uint *triangles = (uint *) calloc(size, sizeof(uint));
  for (int i = 0; i < MAX_THREADS; i++) {
    uint threadStart = omp_csr[i].start;
    uint threadSize = omp_csr[i].table.size;

    for (uint j = 0; j < threadSize; j++) {
      triangles[threadStart + j] = trianglesPerThread[i][j];
    }
  }

  // for (uint i = 0; i < size; i++) {
  //   printf(" %u ", triangles[i]);
  // }

  return triangles;
}


data_arg measureTimeOMP(csr mtx, char *filename, MM_typecode *t, int N, int M, int nz, int MAX_THREADS) {
  struct timeval stop, start;

  gettimeofday(&start, NULL);
  uint *triangles = countTrianglesOMP(mtx, MAX_THREADS);
  gettimeofday(&stop, NULL);

  uint timediff = (stop.tv_sec - start.tv_sec) * 1000000 + stop.tv_usec - start.tv_usec;

  printf("\npthread took %u us for file %s, using %d threads.\n\n", 
    timediff, filename, MAX_THREADS);

  data_arg data = {timediff, triangles};
  return data;
}

int main(int argc, char **argv) {
  int thread_index = atoi(argv[1]);
  int file = atoi(argv[2]);

  FILE *matrixFile;
  int M, N, nz;

  int reps = 12;
  int files_num = 5;
  int thread_num = 3;

  MM_typecode *t;
  char *filenames[5] = {
    "tables/belgium_osm.mtx",
    "tables/dblp-2010.mtx",
    "tables/NACA0015.mtx",
    "tables/mycielskian13.mtx",
    "tables/com-Youtube.mtx"
  };

  char *num_threads[3] = {2, 4, 8};

  char omp[5] = "ompN";
  omp[3] = num_threads[thread_index] + '0';

  FILE *statsFile = fopen("stats/data.csv", "a");
  if (file == 0) {
    fprintf(statsFile, "\n%s", omp);
  }

  csr mtx = readmtx_dynamic(filenames[file], t, N, M, nz);
  csr *mtx_addr = &mtx;

  unsigned long totalTime = 0;
  unsigned long *time_addr = &totalTime;

  for (int rep = 0; rep < reps; rep++) {
    data_arg data = measureTimeOMP(mtx, filenames[file], t, N, M, nz, num_threads[thread_index]);
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