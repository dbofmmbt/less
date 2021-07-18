#ifndef BACK_SUBSTITUTION_H

#define BACK_SUBSTITUTION_H

#include "globals.h"

// Function for the back substitution
void ParallelBackSubstitution(double *pProcRows, double *pProcVector, double *pProcResult, int Size, int RowNum);

#endif