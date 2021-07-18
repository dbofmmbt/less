#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "include/life_cycle.h"
#include "include/data_distribution.h"
#include "include/gaussian_elimination.h"
#include "include/back_substitution.h"
#include "include/presentation.h"

int ProcNum; // Number of the available processes
int ProcRank; // Rank of the current process

int *pParallelPivotPos; // Number of rows selected as the pivot ones
int *pProcPivotIter; // Number of iterations, at which the processor rows were used as the pivot ones
int *pProcInd; // Number of the first row located on the processes
int *pProcNum; // Number of the linear system rows located on the processes

void ParallelResultCalculation(double *pProcRows, double *pProcVector, double *pProcResult, int Size, int RowNum) {
  // Gaussian elimination
  ParallelGaussianElimination (pProcRows, pProcVector, Size, RowNum);
  // Back substitution
  ParallelBackSubstitution (pProcRows, pProcVector, pProcResult, Size, RowNum);
}

// Function for gathering the result vector
void ResultCollection(double *pProcResult, double *pResult, int RowNum) {
  //Gather the whole result vector on every processor
  MPI_Gatherv(pProcResult, pProcNum[ProcRank], MPI_DOUBLE, pResult, pProcNum, pProcInd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

int main(int argc, char* argv[]) {
  double *pMatrix; // Matrix of the linear system
  double *pVector; // Right parts of the linear system
  double *pResult; // Result vector

  double *pProcRows; // Process rows of the matrix A
  double *pProcVector; // Proccess block of the vector b
  double *pProcResult; // Process block of the result

  int Size; // Size of the matrix and vectors
  int RowNum; // Number of the matrix rows

  setvbuf(stdout, 0, _IONBF, 0);

  MPI_Init (&argc, &argv);
  MPI_Comm_rank ( MPI_COMM_WORLD, &ProcRank);
  MPI_Comm_size ( MPI_COMM_WORLD, &ProcNum);

  ProcessInitialization(&pMatrix, &pVector, &pResult, &pProcRows, &pProcVector, &pProcResult, &Size, &RowNum);

  DataDistribution(pMatrix, pProcRows, pVector, pProcVector, Size, RowNum);

  // This is used for vizualization only, uses barriers so it will slow down the program
  PrintDistribution(pMatrix, pVector, pProcRows, pProcVector, Size, RowNum);

  ParallelResultCalculation(pProcRows, pProcVector, pProcResult, Size, RowNum);

  ResultCollection(pProcResult, pResult, RowNum);

  if (ProcRank == 0) PrintResult(pResult, Size);

  ProcessTermination(pMatrix, pVector, pResult, pProcRows, pProcVector, pProcResult);

  MPI_Finalize();

  return 0;
}