#ifndef GLOBALS_H

#define GLOBALS_H

extern int ProcNum;  // Number of the available processes
extern int ProcRank; // Rank of the current process

extern int *pParallelPivotPos; // Number of rows selected as the pivot ones
extern int *pProcPivotIter;    // Number of iterations, at which the processor rows were used as the pivot ones
extern int *pProcInd;          // Number of the first row located on the processes
extern int *pProcNum;          // Number of the linear system rows located on the processes

#endif