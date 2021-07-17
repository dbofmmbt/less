#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MATRIX_SIZE 5;

int ProcNum; // Number of the available processes
int ProcRank; // Rank of the current process

int *pParallelPivotPos; // Number of rows selected as the pivot ones
int *pProcPivotIter; // Number of iterations, at which the processor rows were used as the pivot ones
int *pProcInd; // Number of the first row located on the processes
int *pProcNum; // Number of the linear system rows located on the processes

// Print all matrix elements and respective vector
void PrintMatrixAndVector(double *pMatrix, double *pVector, int Size) {
  for (int row = 0; row < Size; row++) {
    for(int column=0; column < Size; column++) printf("%.0f\t", pMatrix[row * Size + column]);
    printf("-> %.0f\n", pVector[row]);
  }
}

// Random definition of matrix and vector elements
void PopulateMatrixAndVector (double *pMatrix, double *pVector, int Size) {
  for (int i = 0; i < Size; i++) {
    pVector[i] = rand();

    for (int j = 0; j < Size; j++) {
      if (j <= i)
        pMatrix[i * Size + j] = rand();
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

  int RemaingRows = *Size;

  for (int i = 0; i < ProcRank; i++) RemaingRows = RemaingRows - RemaingRows / (ProcNum - i);

  *RowNum = RemaingRows / (ProcNum - ProcRank);

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

  if (ProcRank == 0) printf("-- Parallel Gauss algorithm for solving linear systems --\n");

  ProcessInitialization(&pMatrix, &pVector, &pResult, &pProcRows, &pProcVector, &pProcResult, &Size, &RowNum);

  if (ProcRank == 0) {
    printf("----------------- Initialized matrix -----------------");
    PrintMatrixAndVector(pMatrix, pVector, Size);
  }
  
  DataDistribution(pMatrix, pProcRows, pVector, pProcVector, Size, RowNum);

  MPI_Finalize();

  return 0;
}