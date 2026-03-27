#ifndef RUNMODELSHPP
#define RUNMODELSHPP
#include <stdio.h>
#include "dcerror.hpp"
#include "glfModel.hpp"

int generateSequence(char* s, int len, int start = 0);
int generateSBSVariantForInputMatrix(char* s, int start=0);
void printModel(FILE* fo, const char* LLstr, double LL, glfModel* m,int clean=0);
float evaluateModel(FILE* fo, glfModel* m, const char* name, int* useThese = 0, int* fitThese = 0, double *startingBetas=0, enum thingToMaximise useLinearRegression= LRLNLIKE);
int getBetas(double* beta, glfModel* m);
void setStartingBetasFromCounts(glfModel *m);
int writeRFiles(glfModel* m,char* root);

extern FILE* fGlfLog;

#endif

