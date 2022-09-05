#ifndef FINAL_PROJECT_SPKMEANS_H
#define FINAL_PROJECT_SPKMEANS_H



double **weightedAdjacencyMatrix(double **X,int numOfPoints, int dim);
double **diagonalDegreeMatrix(double **X, int numOfPoints, int dim);
double **lNorm(double **X, int numOfPoints, int dim);
double **jacobi(double **A, int dim);
double **spectralKmeans(double **X, int numOfPoints, int dim, int *k);
double** kmeans(int k, int  max_iter, double eps, double **matrix_points, double **initial_centroids, int dimension, int num_of_points);

#endif





