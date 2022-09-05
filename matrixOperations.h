
#ifndef FINAL_PROJECT_MATRIXOPERATIONS_H
#define FINAL_PROJECT_MATRIXOPERATIONS_H

double **matrixAllocation(int numOfPoints, int dim);
double **matrixSum(double **firstMatrix,double **secondMatrix,int sign, double numOfPoints, double dim);
void matrixMult(double **firstMatrix, double **secondMatrix, int rowFirstMat, int colFirstRowSecond, int colSecondMat);
void sortByEigenValues(double **matrix, int numOfPoints, int dim);
void freeMatrix(double **matrix);
void printMatrix(double **matrix, int numOfPoints, int dim);

#endif

