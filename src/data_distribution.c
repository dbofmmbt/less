#include "include/data_distribution.h"

#include <stdlib.h>
#include <mpi.h>

static void set_send_data(int *pSendNum, int *pSendInd, int Size, int start, int end, int RowsNumber)
{
  for (int i = start; i < end; i++)
  {
    pProcNum[i] = RowsNumber;
    pProcInd[i] = pProcInd[i - 1] + pProcNum[i - 1];
  
    pSendNum[i] = RowsNumber * Size;
    pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
  }
}

// Data distribution among the processes
void DataDistribution(
    double *pMatrix,
    double *pProcRows,
    double *pVector,
    double *pProcVector,

    int Size)
{
  int *pSendNum = malloc(sizeof(int) * ProcNum); // Number of data sent to the process
  int *pSendInd = malloc(sizeof(int) * ProcNum); // Index of first row sent to the process

  int BaseRowNum = Size / ProcNum; // Define the disposition of the matrix rows for the current process

  pProcInd[0] = 0;
  pProcNum[0] = BaseRowNum;

  pSendNum[0] = BaseRowNum * Size;
  pSendInd[0] = 0;

  int RemainingRows = Size % ProcNum;

  set_send_data(pSendNum, pSendInd, Size, 1, ProcNum - RemainingRows, BaseRowNum);
  // When `RemainingRows` is not zero, the last `RemainingRows` processes receive one extra row. 
  set_send_data(pSendNum, pSendInd, Size, ProcNum - RemainingRows, ProcNum, BaseRowNum + 1);

  // Scatter the rows
  // Each process receives `pSendNum[ProcRank]` numbers beginning from `pSendInd[ProcRank]` position of `pMatrix`
  MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows, pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // Scatter the vector rows
  MPI_Scatterv(pVector, pProcNum, pProcInd, MPI_DOUBLE, pProcVector, pProcNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

  free(pSendNum);
  free(pSendInd);
}
