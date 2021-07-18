#include "include/presentation.h"

#include <stdio.h>
#include <mpi.h>

// Print all matrix elements and respective vector
void PrintMatrixAndVector(double *pMatrix, double *pVector, int Size, int RowNum)
{
  for (int row = 0; row < RowNum; row++)
  {
    for (int column = 0; column < Size; column++)
      printf("%.0f\t", pMatrix[row * Size + column]);
    printf("-> %.0f\n", pVector[row]);
  }
}

// Print result vector using pivot pos
void PrintResult(double *pResult, int Size)
{
  printf(YEL "---------------- Result ---------------------------\n" RESET);
  for (int i = 0; i < Size; i++)
    printf("%f\t ", pResult[pParallelPivotPos[i]]);
  printf("\n");
}

// Print initial matrix and every process ones
void PrintDistribution(double *pMatrix, double *pVector, double *pProcRows, double *pProcVector, int Size, int RowNum)
{
  if (ProcRank == 0)
  {
    printf(BLU "---------------- Initial matrix -------------------\n" RESET);
    PrintMatrixAndVector(pMatrix, pVector, Size, Size);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i = 0; i < ProcNum; i++)
  {
    if (ProcRank == i)
    {
      printf(BLU "---------------- Rows from process %d --------------\n" RESET, ProcRank);
      PrintMatrixAndVector(pProcRows, pProcVector, Size, RowNum);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}
