#include "include/back_substitution.h"

#include <mpi.h>

// Calculates the rank of the process, which holds the pivot row
void FindBackPivotRow(int RowIndex, int *IterProcRank, int *IterPivotPos)
{
  for (int i = 0; i < ProcNum - 1; i++)
  {
    if ((pProcInd[i] <= RowIndex) && (RowIndex < pProcInd[i + 1]))
      *IterProcRank = i;
  }

  if (RowIndex >= pProcInd[ProcNum - 1])
    *IterProcRank = ProcNum - 1;

  *IterPivotPos = RowIndex - pProcInd[*IterProcRank];
}

// Function for the back substitution
void ParallelBackSubstitution(double *pProcRows, double *pProcVector, double *pProcResult, int Size, int RowNum)
{
  int IterProcRank; // Rank of the process with the current pivot row
  int IterPivotPos; // Position of the pivot row of the process

  double IterResult; // Calculated value of the current unknown
  double val;

  // Iterations of the back substitution stage
  for (int i = Size - 1; i >= 0; i--)
  {
    // Calculating the rank of the process, which holds the pivot row
    FindBackPivotRow(pParallelPivotPos[i], &IterProcRank, &IterPivotPos);

    // Calculating the unknown
    if (ProcRank == IterProcRank)
    {
      IterResult = pProcVector[IterPivotPos] / pProcRows[IterPivotPos * Size + i];
      pProcResult[IterPivotPos] = IterResult;
    }

    // Broadcasting the value of the current unknown
    MPI_Bcast(&IterResult, 1, MPI_DOUBLE, IterProcRank, MPI_COMM_WORLD);

    // Updating the values of the vector b
    for (int j = 0; j < RowNum; j++)
    {
      if (pProcPivotIter[j] < i)
      {
        val = pProcRows[j * Size + i] * IterResult;
        pProcVector[j] = pProcVector[j] - val;
      }
    }
  }
}
