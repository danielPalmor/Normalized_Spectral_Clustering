#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spkmeans.h"
#include "matrixOperations.h"
#include "vectorOperations.h"

#define EPSILON 0.00001
#define MAX_ROTATION 100
#define FALSE 1
#define TRUE 0

/**
 * This function finds the optimal k according to the Eigengap Heuristic algorithm
 * @param matrix The eigenvector matrix with its eigenvalues in the first row
 * @param dim The point's dimension
 * @return The optimal k
 */
int eigenGapHeuristic(double **matrix, int dim)
{
    int i, k = 1;
    double tempMax, max = fabs(matrix[0][0] - matrix[0][1]);

    for (i = 1; i < floor(dim / 2.0); i++)
    {
        tempMax = fabs(matrix[0][i] - matrix[0][i + 1]);
        if (tempMax > max)
        {
            max = tempMax;
            k = i + 1;
        }
    }
    return k;
}

/**
 * This function creates the A' matrix from a given A matrix
 * @param A The A matrix
 * @param c The c parameter
 * @param s The s parameter
 * @param iMax The row index of the pivot
 * @param jMax The column index of the pivot
 * @param dim The point's dimension
 */
void createATagFromA(double **A, double c, double s, int iMax, int jMax, int dim)
{
    int r;
    double tempAri, tempArj, tempAii, tempAjj;
    for (r = 0; r < dim; r++)
    {
        if (r != iMax && r != jMax)
        {
            tempAri = A[r][iMax];
            tempArj = A[r][jMax];
            A[r][iMax] = c * tempAri - s * tempArj;
            A[r][jMax] = c * tempArj + s * tempAri;
            A[iMax][r] = A[r][iMax];
            A[jMax][r] = A[r][jMax];
        }
    }
    tempAii = A[iMax][iMax];
    tempAjj = A[jMax][jMax];
    A[iMax][iMax] = pow(c,2)*tempAii + pow(s,2)*tempAjj - 2*s*c*A[iMax][jMax];
    A[jMax][jMax] = pow(s,2)*tempAii + pow(c,2)*tempAjj + 2*s*c*A[iMax][jMax];
    A[iMax][jMax] = 0;
    A[jMax][iMax] = 0;
}

/**
 * This function calculates the sum of squares of all
 * off-diagonal elements of the given matrix
 * @param matrix The symmetric matrix
 * @param dim The point's dimension
 * @return The sum of squares of all off-diagonal elements
 */
double calcOff(double **matrix, int dim)
{
    int i, j;
    double result = 0.0;
    for (i = 0; i < dim; i++)
    {
        for (j = i + 1; j < dim; j++)
            result += pow(matrix[i][j], 2);
    }
    return result * 2;
}

/**
 * This function calculates the P or P transpose matrix
 * @param dim The point's dimension
 * @param isT 1 for P, 0 for P transpose
 * @param iMax The row index of the pivot
 * @param jMax The column index of the pivot
 * @param c The c parameter
 * @param s The s parameter
 * @return The P or P transpose matrix
 */
double **calcP(int dim, int isT, int iMax, int jMax, double c, double s)
{
    int i;
    double **P = matrixAllocation(dim, dim);

    for (i = 0; i < dim; i++)
    {
        if (i == iMax || i == jMax)
            P[i][i] = c;
        else
            P[i][i] = 1;
    }

    if (isT == FALSE)
    {
        P[iMax][jMax] = s;
        if (s == 0) /* Inserting 0 instead of -0 */
            P[jMax][iMax] = s;
        else
            P[jMax][iMax] = -s;
    }
    else
    {
        if (s == 0)/* Inserting 0 instead of -0 */
            P[iMax][jMax] = s;
        else
            P[iMax][jMax] = -s;
        P[jMax][iMax] = s;
    }
    return P;
}

/**
 * This function calculates the wam from X matrix
 * @param X The X matrix
 * @param numOfPoints The number of points
 * @param dim The point's dimension
 * @return The wam
 */
double **weightedAdjacencyMatrix(double **X, int numOfPoints, int dim)
{
    int i, j;
    double tempDistance, *tempWeightedVector;
    double **wam = matrixAllocation(numOfPoints, numOfPoints);

    for (i = 0; i < numOfPoints; i++)
    {
        for (j = i + 1; j < numOfPoints; j++)
        {
            tempWeightedVector = malloc(sizeof(double) * dim);
            if (tempWeightedVector == NULL)
            {
                printf("An Error Has Occurred");
                freeMatrix(wam);
                exit(1);
            }
            sumPoints(tempWeightedVector, X[i], X[j], -1, dim);
            tempDistance = sqrt(dotProduct(tempWeightedVector, tempWeightedVector, dim));
            wam[i][j] = exp(-(tempDistance / 2.0));
            wam[j][i] = wam[i][j];
            free(tempWeightedVector);
        }
    }
    return wam;
}

/**
 * This function calculates the ddg matrix from X matrix
 * @param X The X matrix
 * @param numOfPoints The number of points
 * @param dim The point's dimension
 * @return The ddg matrix
 */
double **diagonalDegreeMatrix(double **X, int numOfPoints, int dim)
{
    int i;
    double **ddg = matrixAllocation(numOfPoints, numOfPoints);
    double **wam = weightedAdjacencyMatrix(X, numOfPoints, dim);

    for (i = 0; i < numOfPoints; i++)
        ddg[i][i] = sumScalars(wam[i], numOfPoints);

    freeMatrix(wam);
    return ddg;
}

/**
 * This function calculates the lNorm matrix from X
 * @param X The X matrix
 * @param numOfPoints The number of points
 * @param dim The point's dimension
 * @return The lNorm matrix
 */
double **lNorm(double **X, int numOfPoints, int dim)
{
    int i, j;
    double **ddg = diagonalDegreeMatrix(X, numOfPoints, dim);
    double **wam = weightedAdjacencyMatrix(X, numOfPoints, dim);
    double **lNormMatrix = matrixAllocation(numOfPoints, numOfPoints);

    /* I - D^(-0.5)*W*D^(-0.5) */
    for (i = 0; i < numOfPoints; i++)
    {
        lNormMatrix[i][i] = 1;
        for (j = i + 1; j < numOfPoints; j++)
        {
            lNormMatrix[i][j] = -(wam[i][j] * pow(ddg[i][i], -0.5) * pow(ddg[j][j], -0.5));
            lNormMatrix[j][i] = lNormMatrix[i][j];
        }
    }

    freeMatrix(wam);
    freeMatrix(ddg);
    return lNormMatrix;
}

/**
 * This function performs the Jacobi algorithm from A matrix
 * @param A The A matrix
 * @param dim The point's dimension
 * @return The matrix of eigenvalues and eigenvectors
 */
double **jacobi(double **A, int dim)
{
    int i, j, iMax, jMax, numOfRotations = 0;
    double theta, t, c, s;
    double tempMax, max, offA, tempOffA;
    double **V, **eigenVectorMatrix = matrixAllocation(dim + 1, dim);

    offA = calcOff(A, dim);
    while (numOfRotations < MAX_ROTATION)
    {
        max = fabs(A[0][1]);
        iMax = 0, jMax = 1;

        /* Finding the indexes of the pivot */
        for (i = 0; i < dim; i++)
        {
            for (j = i + 1; j < dim; j++)
            {
                tempMax = fabs(A[i][j]);
                if (tempMax > max)
                {
                    max = tempMax;
                    iMax = i;
                    jMax = j;
                }
            }
        }
        /* If A is diagonal */
        if (A[iMax][jMax] == 0)
        {
            c = 1;
            s = 0;
        }
        else
        {
            theta = (A[jMax][jMax] - A[iMax][iMax]) / (2.0 * A[iMax][jMax]);
            if (theta >= 0)
                t = 1.0 / (theta + sqrt(pow(theta, 2) + 1));
            else
                t = -1.0 / (fabs(theta) + sqrt(pow(theta, 2) + 1));
            c = 1.0 / sqrt(pow(t, 2) + 1);
            s = t * c;
        }

        createATagFromA(A, c, s, iMax, jMax, dim);

        if (numOfRotations == 0)
        {
            V = calcP(dim, FALSE, iMax, jMax, c, s);
        }
        else
        {
            int k;
            double iTemp, jTemp;
            /* V*P */
            for (k = 0; k < dim; k++)
            {
                iTemp = V[k][iMax] * c - V[k][jMax] * s;
                jTemp = V[k][iMax] * s + V[k][jMax] * c;
                V[k][iMax] = iTemp;
                V[k][jMax] = jTemp;
            }
        }

        tempOffA = calcOff(A, dim);
        if (offA - tempOffA <= EPSILON)
            break;
        offA = tempOffA;

        numOfRotations++;
    }

    for (i = 0; i < dim + 1; i++)
    {
        for (j = 0; j < dim; j++)
        {
            if (i == 0)
                eigenVectorMatrix[i][j] = A[j][j];
            else
                eigenVectorMatrix[i][j] = V[i - 1][j];
        }
    }
    freeMatrix(V);
    return eigenVectorMatrix;
}

/**
 * This function performs the SPK algorithm
 * @param X The X matrix
 * @param numOfPoints The number of points
 * @param dim The point's dimension
 * @param k The number of clusters
 * @return The T matrix
 */
double **spectralKmeans(double **X, int numOfPoints, int dim, int *k)
{
    int i, j;
    double sum, **T;
    double **lNormMatrix = lNorm(X, numOfPoints, dim);
    double **jacobiLNorm = jacobi(lNormMatrix, numOfPoints);

    sortByEigenValues(jacobiLNorm, numOfPoints + 1, numOfPoints);

    if (*k == 0)
        *k = eigenGapHeuristic(jacobiLNorm, numOfPoints);
    T = matrixAllocation(numOfPoints, *k);

    /* Calculating U matrix */
    for (i = 1; i <= numOfPoints; i++)
        for (j = 0; j < *k; j++)
            T[i - 1][j] = jacobiLNorm[i][j];

    /* Calculating T from U */
    for (i = 0; i < numOfPoints; i++)
    {
        sum = sqrt(dotProduct(T[i], T[i], *k));
        if (sum != 0)
        {
            for (j = 0; j < *k; j++)
                T[i][j] = T[i][j] / sum;
        }
    }
    freeMatrix(lNormMatrix);
    freeMatrix(jacobiLNorm);
    return T;
}


int main(int argc, char *argv[])
{
    int i, j, dim = 1, numOfPoints = 1;
    double point, **X, **result;
    char ch, comma;
    FILE *input;

    if (argc != 3)
    {
        printf("Invalid Input!");
        return 1;
    }

    input = fopen(argv[2], "r");
    /* Error check */
    if (input == NULL)
    {
        printf("An Error Has Occurred");
        return 1;
    }

    /* Calculating dimension */
    while ((ch = (char) fgetc(input)) != '\n')
    {
        if (ch == ',')
            dim++;
    }

    /* Calculating number of points */
    while ((ch = (char) fgetc(input)) != EOF)
    {
        if (ch == '\n')
            numOfPoints++;
    }
    rewind(input);
    X = matrixAllocation(numOfPoints, dim);

    for (i = 0; i < numOfPoints; i++)
    {
        for (j = 0; j < dim && fscanf(input, "%lf%c", &point, &comma) == 2; j++)
            X[i][j] = point;
    }
    fclose(input);

    if (strcmp(argv[1], "wam") == TRUE)
    {
        result = weightedAdjacencyMatrix(X, numOfPoints, dim);
    }
    else if (strcmp(argv[1], "ddg") == TRUE)
    {
        result = diagonalDegreeMatrix(X, numOfPoints, dim);
    }
    else if (strcmp(argv[1], "lnorm") == TRUE)
    {
        result = lNorm(X, numOfPoints, dim);
    }
    else if (strcmp(argv[1], "jacobi") == TRUE)
    {
        result = jacobi(X, dim);
        printMatrix(result, numOfPoints + 1, dim);
        freeMatrix(result);
        freeMatrix(X);
        return 0;
    }
    else
    {
        freeMatrix(X);
        printf("Invalid Input!");
        return 1;
    }
    printMatrix(result, numOfPoints, numOfPoints);
    freeMatrix(result);
    freeMatrix(X);
    return 0;
}

