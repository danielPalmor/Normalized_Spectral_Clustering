#ifndef FINAL_PROJECT_MATRIXOPERATIONS_H
#define FINAL_PROJECT_MATRIXOPERATIONS_H

double **matrixAllocation(int numOfPoints, int dim);

void sortByEigenValues(double **matrix, int numOfPoints, int dim);

void freeMatrix(double **matrix);

void printMatrix(double **matrix, int numOfPoints, int dim);

#endif

