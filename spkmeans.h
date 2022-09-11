#ifndef FINAL_PROJECT_SPKMEANS_H
#define FINAL_PROJECT_SPKMEANS_H

double **weightedAdjacencyMatrix(double **X, int numOfPoints, int dim);

double **diagonalDegreeMatrix(double **X, int numOfPoints, int dim);

double **lNorm(double **X, int numOfPoints, int dim);

double **jacobi(double **A, int dim);

double **spectralKmeans(double **X, int numOfPoints, int dim, int *k);

double **kmeans(int k, int maxIter, double eps, double **matrixPoints,
                double **initialCentroids, int numOfPoints, int dim);

#endif





