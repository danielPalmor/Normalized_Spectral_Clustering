#include <stdlib.h>
#include <math.h>
#include "vectorOperations.h"

double dotProduct(double *point1, double *point2, int dim)
{
    int i;
    double result = 0;
    for (i = 0; i < dim; ++i)
        result += point1[i]*point2[i];
    return result;
}

void sumPoints(double* sum, double *point1,double *point2, int sign, int dimension)
{
    int i;
    for (i = 0; i < dimension; i++)
        sum[i] = point1[i] + sign*point2[i];
}

double sumScalars(double* point, int dim)
{
    int i;
    double sum = 0;
    for (i = 0; i < dim; i++)
        sum += point[i];
    return sum;
}

void multScalar(double *point, double scalar, int dim)
{
    int i;
    for (i = 0; i < dim; i++)
        point[i] = point[i] * scalar;
}

double calcDistance(int centroidIndex,int pointIndex, double **centroid,double **point, int dim)
{
    double *sub = calloc(dim,sizeof(double));
    double dot;
    sumPoints(sub, point[pointIndex],centroid[centroidIndex],-1, dim);
    dot = sqrt(dotProduct(sub,sub,dim));
    free(sub);
    return dot;
}

double calcDistanceEuclidNorm(double *point1,double *point2, int dimension)
{
    double *sub = calloc(dimension,sizeof(double));
    double dot;
    sumPoints(sub, point1,point2,-1,dimension);
    dot = sqrt(dotProduct(sub,sub,dimension));
    free(sub);
    return dot;
}
