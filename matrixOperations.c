#include <stdio.h>
#include <stdlib.h>
#include "matrixOperations.h"

double **matrixAllocation(int numOfPoints, int dim)
{
    int i;
    double **matrix, *array;
    matrix = calloc(numOfPoints, sizeof(double*));
    if(matrix == NULL)
    {
        printf("An Error Has Occurred");
        exit(1);
    }

    array = calloc(numOfPoints * dim, sizeof(double));
    if(array == NULL)
    {
        printf("An Error Has Occurred");
        free(matrix);
        exit(1);
    }
    for (i = 0; i < numOfPoints; i++) {
        matrix[i] = array + i * dim;
    }
    return matrix;
}

double **matrixSum(double **firstMatrix,double **secondMatrix,int sign, double numOfPoints, double dim)
{
    int i,j;
    double **matrix = matrixAllocation(numOfPoints, dim);

    for (i = 0; i < numOfPoints; i++)
    {
        for (j = 0; j < dim; j++)
        {
            matrix[i][j] = firstMatrix[i][j] + sign*secondMatrix[i][j];
        }
    }
    return matrix;
}

void matrixMult(double **firstMatrix, double **secondMatrix, int rowFirstMat, int colFirstRowSecond, int colSecondMat)
{
    int i,j,m,k;
    double **multMatrix = matrixAllocation(rowFirstMat,colSecondMat);
    for (i = 0; i < rowFirstMat; i++)
    {
        for (j = 0; j < colSecondMat; j++)
            for (m = 0; m < colFirstRowSecond; m++)
                multMatrix[i][j] += firstMatrix[i][m]*secondMatrix[m][j];
        for(k = 0; k < colSecondMat; k++)
            firstMatrix[i][k] = multMatrix[i][k];
    }
    freeMatrix(multMatrix);
}

void sortByEigenValues(double **matrix, int numOfPoints, int dim)
{
    int i, j, k, maxCol;
    double temp, max;
    for (i = 0; i < dim; i++)
    {
        max = matrix[0][i];
        maxCol = i;
        for (j = i+1; j < dim; j++)
        {
            if(matrix[0][j] > max)
            {
                max = matrix[0][j];
                maxCol = j;
            }
        }
        for (k = 0; k < numOfPoints; k++)
        {
            temp = matrix[k][i];
            matrix[k][i] = matrix[k][maxCol];
            matrix[k][maxCol] = temp;
        }
    }
}

void freeMatrix(double **matrix)
{
    free(matrix[0]);
    free(matrix);
}

void printMatrix(double **matrix, int numOfPoints, int dim)
{
    int i,j;
    for (i = 0; i < numOfPoints; i++)
    {
        for (j = 0; j < dim ; j++)
        {
            if (j == dim - 1)
                printf("%.4f\n", matrix[i][j]);
            else
                printf("%.4f,", matrix[i][j]);
        }
    }
}
