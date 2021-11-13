#ifndef PARALLEL_HELPERS_H
#define PARALLEL_HELPERS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>

#include "../headers/mmio.h"
#include "../headers/csr_arg.h"
#include "../headers/csr.h"
#include "../headers/helpers.h"
#include "../headers/parallel_helpers.h"

#define MAX_THREADS "2"

uint countTrianglesCilk(csr table) {
  int max_threads = atoi(MAX_THREADS);
  __cilkrts_set_param("nworkers", MAX_THREADS);
  __cilkrts_init();

  uint size = table.size;
  uint interval = size / max_threads;
  csr_arg *cilkcsr = (csr_arg *) malloc(max_threads * sizeof(csr_arg));

  uint totalSize = 0;
  for (int i = 0; i < max_threads; i++) {
    uint start = i * interval;
    uint end = (i == max_threads - 1) ? size : (i + 1) * interval;
    uint size = end - start;

    printf("start = %u\tend = %u\tsize = %u\n", start, end, size);
    uint nonzeros = table.rowIndex[end] - table.rowIndex[start];
    printf("nonzeros = %u\n", nonzeros);

    cilkcsr[i].table.size = size;
    cilkcsr[i].table.rowIndex = (uint *) malloc((size+1) * sizeof(uint));
    cilkcsr[i].table.colIndex = (uint *) malloc(nonzeros * sizeof(uint));
    cilkcsr[i].table.values = (int *) malloc(nonzeros * sizeof(int));
    cilkcsr[i].id = i;
    cilkcsr[i].start = start;
    cilkcsr[i].end = end;

    uint individualSize = sizeof(int) + 3 * sizeof(uint) +
      (size + 1) * sizeof(uint) + nonzeros * sizeof(uint) + nonzeros * sizeof(int);

    totalSize += individualSize;
  }

  cilkcsr = (csr_arg *) realloc(cilkcsr, totalSize);
  printf("\nInitialized cilkcsr\n\n");

  cilk_for (int i = 0; i < max_threads; i++) {
    printf("cilk_for thread: %d\n", i);
    int index = cilkcsr[i].id;

    cilkcsr[index].table = hadamardSingleStep(table, cilkcsr[index].start, cilkcsr[index].end);
    printf("thread end: %d\n", i);
  }

  printf("\nExited the loop and cilk_sync\n\n");

  uint *trianglesPerThread = (uint *) calloc(max_threads, sizeof(uint));
  uint totalTriangles = 0;
  cilk_for(int i = 0; i < max_threads; i++) {
    int index = cilkcsr[i].id;

    trianglesPerThread[index] = countTriangles(cilkcsr[index].table);
    printf("trianglesPerThread[%d] = %u\n", index, trianglesPerThread[index]);
  }

  for (int i = 0; i < max_threads; i++) {
    totalTriangles += trianglesPerThread[i];
  }

  return totalTriangles;
}






#endif