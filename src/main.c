#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define MATRIX_SIZE 3;

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define BLU   "\x1b[36m"
#define RESET "\x1B[0m"

int ProcNum; // Number of the available processes
int ProcRank; // Rank of the current process

int *pParallelPivotPos; // Number of rows selected as the pivot ones
int *pProcPivotIter; // Number of iterations, at which the processor rows were used as the pivot ones
int *pProcInd; // Number of the first row located on the processes
int *pProcNum; // Number of the linear system rows located on the processes

// Random definition of matrix and vector elements
void PopulateMatrixAndVector (double *pMatrix, double *pVector, int Size) {
  srand(time(NULL));

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

// Subtracts the pivot row from the process rows, which have not been used as the pivot ones yet
void ParallelEliminateColumns(double *pProcRows, double *pProcVector, double *pPivotRow, int Size, int RowNum, int Iter) { 
  double multiplier;

  for (int i = 0; i < RowNum; i++) {
    if (pProcPivotIter[i] == -1) {
      multiplier = pProcRows[i * Size+Iter] / pPivotRow[Iter];

      for (int j = Iter; j < Size; j++) pProcRows[i * Size + j] -= pPivotRow[j] * multiplier;

      pProcVector[i] -= pPivotRow[Size] * multiplier;
    }
  } 
}

// Function for the Gausian elimination
void ParallelGaussianElimination (double *pProcRows, double *pProcVector, int Size, int RowNum) {
  int PivotPos;  // Position of the pivot row in the process stripe 
  double MaxValue;

  // Structure for the pivot row selection
  struct { double MaxValue; int ProcRank; } ProcPivot, Pivot;

  // pPivotRow is used for storing the pivot row and the corresponding element of the vector b
  double *pPivotRow = malloc(sizeof(double) * Size + 1);

  // The iterations of the Gaussian elimination stage
  for (int i = 0; i < Size; i++) {
    // Calculating the local pivot row
    MaxValue = 0;

    for (int j=0; j<RowNum; j++) {
      if ((pProcPivotIter[j] == -1) && (MaxValue < fabs(pProcRows[j * Size+i]))) {
        MaxValue = fabs(pProcRows[j * Size + i]);
        PivotPos = j;
      }
    }
    ProcPivot.MaxValue = MaxValue;
    ProcPivot.ProcRank = ProcRank;

    // Finding process with MaxValue
    MPI_Allreduce(&ProcPivot, &Pivot, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

    // Broadcasting the pivot row
    if (ProcRank == Pivot.ProcRank) {
      pProcPivotIter[PivotPos]= i; //iteration number
      pParallelPivotPos[i] = pProcInd[ProcRank] + PivotPos;

      // Fill the pivot row
      for (int j=0; j < Size; j++) {
        pPivotRow[j] = pProcRows[PivotPos * Size + j];
      }

      pPivotRow[Size] = pProcVector[PivotPos];
    }

    MPI_Bcast(&pParallelPivotPos[i], 1, MPI_INT, Pivot.ProcRank, MPI_COMM_WORLD); 

    MPI_Bcast(pPivotRow, Size + 1, MPI_DOUBLE, Pivot.ProcRank, MPI_COMM_WORLD);

    ParallelEliminateColumns(pProcRows, pProcVector, pPivotRow, Size, RowNum, i);
  }
}

// Calculates the rank of the process, which holds the pivot row
void FindBackPivotRow(int RowIndex, int *IterProcRank, int *IterPivotPos) {
  for (int i = 0; i < ProcNum - 1; i++) {
    if ((pProcInd[i] <= RowIndex) && (RowIndex < pProcInd[i + 1])) *IterProcRank = i;
  }

  if (RowIndex >= pProcInd[ProcNum-1]) *IterProcRank = ProcNum - 1;

  *IterPivotPos = RowIndex - pProcInd[*IterProcRank];
}

// Function for the back substitution
void ParallelBackSubstitution (double *pProcRows, double *pProcVector, double *pProcResult, int Size, int RowNum) {
  int IterProcRank;    // Rank of the process with the current pivot row
  int IterPivotPos;    // Position of the pivot row of the process

  double IterResult;   // Calculated value of the current unknown
  double val;

  // Iterations of the back substitution stage
  for (int i = Size-1; i >= 0; i--) {
    // Calculating the rank of the process, which holds the pivot row
    FindBackPivotRow(pParallelPivotPos[i], &IterProcRank, &IterPivotPos);
    
    // Calculating the unknown
    if (ProcRank == IterProcRank) {
      IterResult = pProcVector[IterPivotPos] / pProcRows[IterPivotPos * Size + i];
      pProcResult[IterPivotPos] = IterResult;
    }

    // Broadcasting the value of the current unknown
    MPI_Bcast(&IterResult, 1, MPI_DOUBLE, IterProcRank, MPI_COMM_WORLD);

    // Updating the values of the vector b
    for (int j=0; j < RowNum; j++) {
      if (pProcPivotIter[j] < i) {
        val = pProcRows[j * Size + i] * IterResult;
        pProcVector[j] = pProcVector[j] - val;
      }
    }
  }
}

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

// Print all matrix elements and respective vector
void PrintMatrixAndVector(double *pMatrix, double *pVector, int Size, int RowNum) {
  for (int row = 0; row < RowNum; row++) {
    for(int column=0; column < Size; column++) printf("%.0f\t", pMatrix[row * Size + column]);
    printf("-> %.0f\n", pVector[row]);
  }
}

// Print result vector
void PrintResult(double *pResult, int Size) {
  printf(BLU "---------------- Result ---------------------------\n" RESET);
  for (int i = 0; i < Size; i++) printf("%f\t ", pResult[i]);
  printf("\n");
}


// Print initial matrix and every process ones
void PrintDistribution (double* pMatrix, double* pVector, double* pProcRows, double* pProcVector, int Size, int RowNum) {
  if (ProcRank == 0) {
    printf(BLU "---------------- Initial matrix -------------------\n" RESET);
    PrintMatrixAndVector(pMatrix, pVector, Size, Size);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  for (int i=0; i < ProcNum; i++) {
    if (ProcRank == i) {
      printf(BLU "---------------- Rows from process %d --------------\n" RESET, ProcRank);
      PrintMatrixAndVector(pProcRows, pProcVector, Size, RowNum);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

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

  PrintDistribution(pMatrix, pVector, pProcRows, pProcVector, Size, RowNum);

  ParallelResultCalculation(pProcRows, pProcVector, pProcResult, Size, RowNum);

  ResultCollection(pProcResult, pResult, RowNum);

  if (ProcRank == 0) PrintResult(pResult, Size);

  ProcessTermination(pMatrix, pVector, pResult, pProcRows, pProcVector, pProcResult);

  MPI_Finalize();

  return 0;
}