#include "include/life_cycle.h"

#include <stdlib.h>
#include <time.h>
#include <mpi.h>

// Random definition of matrix and vector elements
void PopulateMatrixAndVector (double *pMatrix, double *pVector, int Size) {
  srand(1);

  for (int i = 0; i < Size; i++) {
    pVector[i] = rand() % 100;

    for (int j = 0; j < Size; j++) {
      if (j <= i)
        pMatrix[i * Size + j] = rand() % 100;
      else
        pMatrix[i * Size + j] = 0;
    }
  }
}

// Memory allocation and data initialization
void ProcessInitialization(
  double **pMatrix,
  double **pVector,
  double **pResult,

  double **pProcRows,
  double **pProcVector,
  double **pProcResult,

  int *Size,
  int *RowNum
) {
  if (ProcRank == 0) *Size = MATRIX_SIZE;

  MPI_Bcast(Size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  *RowNum = *Size / ProcNum;
  
  int RemainingRows = *Size % ProcNum;
 
  if (ProcRank >= ProcNum - RemainingRows) {
    *RowNum += 1;
  }

  *pProcRows = malloc(sizeof(double) * ((*RowNum) * (*Size)));
  *pProcVector = malloc(sizeof(double) * (*RowNum));
  *pProcResult = malloc(sizeof(double) * (*RowNum));

  pParallelPivotPos = malloc(sizeof(int) * (*Size));
  pProcPivotIter = malloc(sizeof(int) * (*RowNum));
  pProcInd = malloc(sizeof(int) * ProcNum);
  pProcNum = malloc(sizeof(int) * ProcNum);

  for (int i = 0; i < *RowNum; i++) pProcPivotIter[i] = -1;

  if (ProcRank == 0) {
    *pMatrix =  malloc(sizeof(double) * ((*Size) * (*Size)));
    *pVector = malloc(sizeof(double) * (*Size));
    *pResult = malloc(sizeof(double) * (*Size));

    PopulateMatrixAndVector(*pMatrix, *pVector, *Size);
  }
};

void ProcessTermination(
  double *pMatrix,
  double *pVector,
  double *pResult,
  double *pProcRows,
  double *pProcVector,
  double *pProcResult
) {
  if (ProcRank == 0) {
    free(pMatrix);
    free(pVector);
    free(pResult);
  }

  free(pProcRows);
  free(pProcVector);
  free(pProcResult);
  free(pParallelPivotPos);
  free(pProcPivotIter);
  free(pProcInd);
  free(pProcNum);
}
