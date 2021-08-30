#include "include/gaussian_elimination.h"

#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// Subtracts the pivot row from the process rows, which have not been used as the pivot ones yet
void ParallelEliminateColumns(double *pProcRows, double *pProcVector, double *pPivotRow, int Size, int RowNum, int Iter)
{
#pragma omp parallel for
  for (int i = 0; i < RowNum; i++)
  {
    if (pProcPivotIter[i] == -1)
    {
      double multiplier = pProcRows[i * Size + Iter] / pPivotRow[Iter];

      for (int j = Iter; j < Size; j++)
        pProcRows[i * Size + j] -= pPivotRow[j] * multiplier;

      pProcVector[i] -= pPivotRow[Size] * multiplier;
    }
  }
}

// Function for the Gausian elimination
void ParallelGaussianElimination(double *pProcRows, double *pProcVector, int Size, int RowNum)
{
  int PivotPos; // Position of the pivot row in the process stripe
  double MaxValue;

  // Structure for the pivot row selection
  struct
  {
    double MaxValue;
    int ProcRank;
  } ProcPivot, Pivot;

  // pPivotRow is used for storing the pivot row and the corresponding element of the vector b
  double *pPivotRow = malloc(sizeof(double) * Size + 1);

  // The iterations of the Gaussian elimination stage
  for (int i = 0; i < Size; i++)
  {
    // Calculating the local pivot row
    MaxValue = 0;

    for (int j = 0; j < RowNum; j++)
    {
      if ((pProcPivotIter[j] == -1) && (MaxValue < fabs(pProcRows[j * Size + i])))
      {
        MaxValue = fabs(pProcRows[j * Size + i]);
        PivotPos = j;
      }
    }
    ProcPivot.MaxValue = MaxValue;
    ProcPivot.ProcRank = ProcRank;

    // Finding process with MaxValue
    MPI_Allreduce(&ProcPivot, &Pivot, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

    // Broadcasting the pivot row
    if (ProcRank == Pivot.ProcRank)
    {
      pProcPivotIter[PivotPos] = i; //iteration number
      pParallelPivotPos[i] = pProcInd[ProcRank] + PivotPos;

      // Fill the pivot row
      for (int j = 0; j < Size; j++)
      {
        pPivotRow[j] = pProcRows[PivotPos * Size + j];
      }

      pPivotRow[Size] = pProcVector[PivotPos];
    }

    MPI_Bcast(&pParallelPivotPos[i], 1, MPI_INT, Pivot.ProcRank, MPI_COMM_WORLD);

    MPI_Bcast(pPivotRow, Size + 1, MPI_DOUBLE, Pivot.ProcRank, MPI_COMM_WORLD);

    ParallelEliminateColumns(pProcRows, pProcVector, pPivotRow, Size, RowNum, i);
  }
}
