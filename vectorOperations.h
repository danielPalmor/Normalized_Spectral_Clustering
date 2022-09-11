#ifndef FINAL_PROJECT_VECTORSOPERATIONS_H
#define FINAL_PROJECT_VECTORSOPERATIONS_H

double dotProduct(double *point1, double *point2, int dim);

void sumPoints(double *sum, double *point1, double *point2, int sign, int dim);

double sumScalars(double *point, int dim);

void multScalar(double *point, double scalar, int dim);

double calcDistance(double *point1, double *point2, int dim);

#endif
