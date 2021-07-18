#include "include/data_distribution.h"

#include <stdlib.h>
#include <mpi.h>

// Data distribution among the processes
void DataDistribution(
  double *pMatrix,
  double *pProcRows,
  double *pVector,
  double *pProcVector,

  int Size, 
  int RowNum
) {
  int *pSendNum = malloc(sizeof(int) * ProcNum); // Number of data sent to the process
  int *pSendInd = malloc(sizeof(int) * ProcNum); // Index of first row sent to the process

  int RestRows = Size; // Number of rows, that have not been  distributed yet

  RowNum = (Size / ProcNum); // Define the disposition of the matrix rows for the current process

  pSendNum[0] = RowNum * Size;
  pSendInd[0] = 0;

  for (int i = 1; i < ProcNum; i++) {
    RestRows -= RowNum;
    RowNum = RestRows / (ProcNum - i);
    pSendNum[i] = RowNum * Size;
    pSendInd[i] = pSendInd[i-1] + pSendNum[i-1];
  }

  // Scatter the rows
  MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows, pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Define the disposition of the matrix rows for current process
  RestRows = Size;
  pProcInd[0] = 0;
  pProcNum[0] = Size / ProcNum;

  for (int i = 1; i < ProcNum; i++) {
    RestRows -= pProcNum[i-1];
    pProcNum[i] = RestRows/(ProcNum-i);
    pProcInd[i] = pProcInd[i-1]+pProcNum[i-1];
  }

  // Scatter the vector rows
  MPI_Scatterv(pVector, pProcNum, pProcInd, MPI_DOUBLE, pProcVector, pProcNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  free(pSendNum);
  free(pSendInd);
}
