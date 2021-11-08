/*
 * sequential.c 
 * Convert a square N x N matrix into the CSR format, made for sparse matrices:
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 * 
 * Authors: Antonios Antoniou - 9482
 *          Efthymios Grigorakis - 9694
 *    
 * 2021 Aristotle University of Thessaloniki
 * Parallel and Distributed Systems.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include "headers/csr.h"
// #include "headers/csr_operations.h"
#include "headers/mmio.h"
#include "headers/helpers.h"

#define SIZE 6 


//Reads an .mtx and directly converts it to CSR.
csr readmtx(char *mtx, MM_typecode *t, int N, int M, int nz) {
	FILE *matrixFile;
	matrixFile = fopen(mtx, "r");
	int banner = mm_read_banner(matrixFile, t);
	int result = mm_read_mtx_crd_size(matrixFile, &M, &N, &nz);

	printf("banner: %d\tresult: %d\tnonzeros: %d\tM: %d\tN: %d\n",
	banner, result, nz, M, N);
 
	// // Display error messages and abort, if the matrix isn't square or hasn't been read properly.
	// if (N != M) {
	// 	printf("N and M are not equal. The matrix isn't square. Aborting...");
	// 	csr returnError = {0, NULL, NULL, NULL};
	// 	return returnError;
	// }

	// if (banner != 0 || result != 0) {
	// 	printf("Error. Couldn't process the .mtx file!");
	// 	csr returnError = {0, NULL, NULL, NULL};
	// 	return returnError;	
	// }

  int *row , *col ;

    row = (int *) malloc(2 * nz * sizeof(int));
    col = (int *) malloc(2 * nz * sizeof(int));

    for(int i = 0 ; i < nz ; i++){

    fscanf(matrixFile , "%d  %d\n" , &row[i] , &col[i]);
        // 1-based 
        row[i]--;
        col[i]--;
        //symetric matrix
        row[i+nz] = col[i];
        col[i+nz] = row[i];
    }



    csr Table ; 
    Table.rowIndex = (long *) malloc((M+1) * sizeof(long));
    Table.colIndex = (long *) malloc(2 * nz * sizeof(long));
    Table.values = (int *) malloc(2 * nz * sizeof(int));




    int current = 0 ;

    for(int i = 0 ; i < M; i++){
        for(int j = 0 ; j <2*nz ; j++){
            if(row[j] == i ){
                Table.rowIndex[i+1]++;
                Table.values[current] = 1 ;
                Table.colIndex[current] = col[j];
                current ++ ;
            }
            Table.rowIndex[i+2] = Table.rowIndex[i+1];
        }
    }
}


int main(int argc, char **argv) {
  FILE *matrixFile;
  int M, N, nz;
  MM_typecode *t;
  char *mtx = "mycielskian4.mtx";



  int **test = makeRandomSparseTable(SIZE);
  printTable(test, SIZE);

  csr csr_mtx = readmtx(mtx, t, N, M, nz);


}