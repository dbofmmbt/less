#ifndef DATA_DISTRIBUTION_H

#define DATA_DISTRIBUTION_H

#include "globals.h"

// Data distribution among the processes
void DataDistribution(
  double *pMatrix,
  double *pProcRows,
  double *pVector,
  double *pProcVector,

  int Size, 
  int RowNum
);

#endif