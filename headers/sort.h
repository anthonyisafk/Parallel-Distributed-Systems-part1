#include <stdio.h>

// Taken directly from:
// https://stackoverflow.com/questions/3893937/sorting-an-array-in-c
int compare( const void* a, const void* b) {
  int int_a = * ( (int*) a );
  int int_b = * ( (int*) b );

  return (a > b) - (a < b);
}