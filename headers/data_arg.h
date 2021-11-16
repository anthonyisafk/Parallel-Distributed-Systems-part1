/**
 * data_arg.h
 * 
 * A struct that's used for the presentation of the assignment.
 * Store the time it took for the algorithm to finish
 * and the resulting array of triangle values.
 */

#ifndef DATA_ARG_H
#define DATA_ARG_H

#include <stdio.h>
#include <stdlib.h>

typedef struct {
  uint time;
  uint *triangles;
} data_arg;

#endif