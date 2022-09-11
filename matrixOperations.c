#include <stdio.h>
#include <stdlib.h>
#include "matrixOperations.h"

/**
 * This function allocates memory for a matrix
 * @param numOfPoints The number of points (rows)
 * @param dim The point's dimension (columns)
 * @return The matrix
 */
double **matrixAllocation(int numOfPoints, int dim)
{
    int i;
    double **matrix, *array;
    matrix = calloc(numOfPoints, sizeof(double *));
    if (matrix == NULL)
    {
        printf("An Error Has Occurred");
        exit(1);
    }

    array = calloc(numOfPoints * dim, sizeof(double));
    if (array == NULL)
    {
        printf("An Error Has Occurred");
        free(matrix);
        exit(1);
    }
    for (i = 0; i < numOfPoints; i++)
    {
        matrix[i] = array + i * dim;
    }
    return matrix;
}

/**
 * This function sorts a given matrix by it's eigenvalues
 * @param matrix The matrix of eigenvectors. The first row contains the eigenvalues
 * @param numOfPoints The number of points
 * @param dim The point's dimension
 */
void sortByEigenValues(double **matrix, int numOfPoints, int dim)
{
    int i, j, k, maxCol;
    double temp, max;
    for (i = 0; i < dim; i++)
    {
        max = matrix[0][i];
        maxCol = i;
        /* Finding max eigenValue */
        for (j = i + 1; j < dim; j++)
        {
            if (matrix[0][j] > max)
            {
                max = matrix[0][j];
                maxCol = j;
            }
        }
        /* Swapping the ith column with the maxCol column */
        for (k = 0; k < numOfPoints; k++)
        {
            temp = matrix[k][i];
            matrix[k][i] = matrix[k][maxCol];
            matrix[k][maxCol] = temp;
        }
    }
}

/**
 * This function frees memory of a given matrix
 * @param matrix The matrix
 */
void freeMatrix(double **matrix)
{
    free(matrix[0]);
    free(matrix);
}

/**
 * This function prints a given matrix
 * @param matrix The matrix
 * @param numOfPoints The number of points
 * @param dim The point's dimension
 */
void printMatrix(double **matrix, int numOfPoints, int dim)
{
    int i, j;
    for (i = 0; i < numOfPoints; i++)
    {
        for (j = 0; j < dim; j++)
        {
            if (j == dim - 1)
                printf("%.4f\n", matrix[i][j]);
            else
                printf("%.4f,", matrix[i][j]);
        }
    }
}
