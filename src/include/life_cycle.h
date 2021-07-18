#ifndef LIFE_CYCLE_H

#define LIFE_CYCLE_H

#include "globals.h"

#define MATRIX_SIZE 5;

// Initialize resources
void ProcessInitialization(
  double **pMatrix,
  double **pVector,
  double **pResult,

  double **pProcRows,
  double **pProcVector,
  double **pProcResult,

  int *Size,
  int *RowNum
);

// Must be called after the algorithm is done to free resources
void ProcessTermination(
  double *pMatrix,
  double *pVector,
  double *pResult,
  double *pProcRows,
  double *pProcVector,
  double *pProcResult
);

#endif