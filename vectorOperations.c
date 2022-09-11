#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "vectorOperations.h"

/**
 * This function calculates the dot product of two given points
 * @param point1 The first point
 * @param point2 The second point
 * @param dim The point's dimension
 * @return The dot product of the points
 */
double dotProduct(double *point1, double *point2, int dim)
{
    int i;
    double result = 0;
    for (i = 0; i < dim; ++i)
        result += point1[i] * point2[i];
    return result;
}

/**
 * This function sums/subtracts two given points
 * @param sum The result of the function
 * @param point1 The first point
 * @param point2 The second point
 * @param sign 1 for sum , -1 for subtraction
 * @param dim  The point's dimension
 */
void sumPoints(double *sum, double *point1, double *point2, int sign, int dim)
{
    int i;
    for (i = 0; i < dim; i++)
        sum[i] = point1[i] + sign * point2[i];
}

/**
 * This function sums the scalars of a given point
 * @param point The point
 * @param dim The point's dimension
 * @return The sum of the scalars
 */
double sumScalars(double *point, int dim)
{
    int i;
    double sum = 0;
    for (i = 0; i < dim; i++)
        sum += point[i];
    return sum;
}

/**
 * This function multiplies a point by a scalar
 * @param point The point
 * @param scalar The scalar
 * @param dim The point's dimension
 */
void multScalar(double *point, double scalar, int dim)
{
    int i;
    for (i = 0; i < dim; i++)
        point[i] = point[i] * scalar;
}

/**
 * This function calculates the distance of two given points
 * @param point1 The first point
 * @param point2 The second point
 * @param dim The point's dimension
 * @return The distance
 */
double calcDistance(double *point1, double *point2, int dim)
{
    double dot;
    double *sub = calloc(dim, sizeof(double));
    if (sub == NULL)
    {
        printf("An Error Has Occurred");
        exit(1);
    }
    sumPoints(sub, point1, point2, -1, dim);
    dot = sqrt(dotProduct(sub, sub, dim));
    free(sub);
    return dot;
}
