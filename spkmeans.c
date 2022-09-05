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

int eigenGapHeuristic(double **lNormMatrix, int dim)
{
    int i, k = 1;
    double tempMax, max = fabs(lNormMatrix[0][0] - lNormMatrix[0][1]);

    for (i = 1; i < floor(dim / 2.0); i++)
    {
        tempMax = fabs(lNormMatrix[0][i] - lNormMatrix[0][i+1]);
        if (tempMax > max)
        {
            max = tempMax;
            k = i+1;
        }
    }
    return k;
}

void createATagFromA(double **A, double c, double s, int iMax, int jMax, int dim)
{
    int r;
    double tempAri,tempArj, tempAii, tempAjj;
    for (r = 0; r < dim ; r++)
    {
        if (r != iMax && r != jMax)
        {
            tempAri = A[r][iMax];
            tempArj = A[r][jMax];
            A[r][iMax] = c*tempAri - s*tempArj;
            A[r][jMax] = c*tempArj + s*tempAri;
            A[iMax][r] = A[r][iMax];
            A[jMax][r] = A[r][jMax];
        }
    }
    tempAii = A[iMax][iMax];
    tempAjj = A[jMax][jMax];
    A[iMax][iMax] = pow(c,2)*tempAii + pow(s,2)*A[jMax][jMax] - 2*s*c*A[iMax][jMax];
    A[jMax][jMax] = pow(s,2)*tempAii + pow(c,2)*tempAjj + 2*s*c*A[iMax][jMax];
    A[iMax][jMax] = 0;
    A[jMax][iMax] = 0;
}

double calcOff(double **matrix, int dim)
{
    int i,j;
    double result = 0.0;
    for (i = 0; i < dim ; i++)
    {
        for (j = i+1;  j < dim ; j++)
            result += pow(matrix[i][j], 2);
    }
    return result * 2;
}

double **calcP(int dim, int isT, int iMax, int jMax, double c, double s)
{
    int i;
    double **P = matrixAllocation(dim,dim);

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
        if (s == 0)
            P[jMax][iMax] = s;
        else
            P[jMax][iMax] = -s;
    }
    else
    {
        if (s == 0)
            P[iMax][jMax] = s;
        else
            P[iMax][jMax] = -s;
        P[jMax][iMax] = s;
    }
    return P;
}

double **weightedAdjacencyMatrix(double **X, int numOfPoints, int dim)
{
    int i,j;
    double tempDistance, *tempWeightedVector;
    double **wam = matrixAllocation(numOfPoints,numOfPoints);

    for (i = 0; i < numOfPoints; i++)
    {
        for (j = i + 1; j < numOfPoints; j++)
        {
                tempWeightedVector = malloc(sizeof (double) * dim);
                if(tempWeightedVector == NULL)
                {
                    printf("An Error Has Occurred");
                    freeMatrix(wam);
                    exit(1);
                }
                sumPoints(tempWeightedVector,X[i],X[j],-1,dim);
                tempDistance = sqrt(dotProduct(tempWeightedVector,tempWeightedVector, dim));
                wam[i][j] = exp(-(tempDistance / 2.0));
                wam[j][i] = wam[i][j];
                free(tempWeightedVector);
        }
    }
    return wam;
}

double **diagonalDegreeMatrix(double **X, int numOfPoints, int dim)
{
    int i;
    double **ddg = matrixAllocation(numOfPoints, numOfPoints);
    double **wam = weightedAdjacencyMatrix(X,numOfPoints, dim);

    for (i = 0;  i < numOfPoints; i++)
        ddg[i][i] = sumScalars(wam[i], numOfPoints);

    freeMatrix(wam);
    return ddg;
}


double **lNorm(double **X, int numOfPoints, int dim)
{
    int i;
    double temp, **lNormMatrix;
    double **ddg = diagonalDegreeMatrix(X, numOfPoints, dim);
    double **wam = weightedAdjacencyMatrix(X, numOfPoints,dim);
    double **eyeMatrix = matrixAllocation(numOfPoints,numOfPoints);
    double **templNormMatrix = matrixAllocation(numOfPoints,numOfPoints);

    for (i = 0; i < numOfPoints; i++)
        eyeMatrix[i][i] = 1;
    for (i = 0;  i < numOfPoints ; i++)
    {
        temp = pow(ddg[i][i],-0.5);
        ddg[i][i] = temp;
        templNormMatrix[i][i] = temp;
    }

    matrixMult(ddg,wam,numOfPoints,numOfPoints,numOfPoints);
    matrixMult(ddg,templNormMatrix,numOfPoints,numOfPoints,numOfPoints);
    lNormMatrix = matrixSum(eyeMatrix,ddg,-1, numOfPoints,numOfPoints);

    freeMatrix(wam);
    freeMatrix(ddg);
    freeMatrix(eyeMatrix);
    freeMatrix(templNormMatrix);
    return lNormMatrix;
}

double **jacobi(double** A, int dim)
{
    int i, j, iMax, jMax, numOfRotations = 0;
    double theta, t, c, s;
    double tempMax, max, offA, tempOffA;
    double **P, **V;
    double **eigenVectorMatrix = matrixAllocation(dim+1, dim);

    offA = calcOff(A, dim);
    while (numOfRotations < MAX_ROTATION)
    {
        max = fabs(A[0][1]);
        iMax = 0, jMax = 1;
        for (i = 0; i < dim ; i++)
        {
            for (j = i+1;  j < dim ; j++)
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
        if (A[iMax][jMax] == 0)
        {
            c = 1;
            s = 0;
        }
        else
        {
            theta = (A[jMax][jMax] - A[iMax][iMax]) / (2.0 * A[iMax][jMax]);
            if (theta >= 0)
                t = 1.0 / (theta + sqrt(pow(theta,2) + 1));
            else
                t = -1.0 / (fabs(theta) + sqrt(pow(theta,2) + 1));
            c = 1.0 / sqrt(pow(t,2) + 1);
            s = t * c;
        }

        createATagFromA(A,c,s,iMax,jMax,dim);

        if (numOfRotations == 0)
        {
            V = calcP(dim,FALSE, iMax, jMax, c, s);
        }
        else
        {
            P = calcP(dim,FALSE, iMax, jMax, c, s);
            matrixMult(V,P,dim,dim,dim);
            freeMatrix(P);
        }

        tempOffA = calcOff(A,dim);
        if(offA - tempOffA <= EPSILON)
            break;
        offA = tempOffA;

        numOfRotations++;
    }

    for (i = 0; i < dim+1 ; i++)
    {
        for(j = 0; j < dim; j++)
        {
            if (i == 0)
            {
                eigenVectorMatrix[i][j] = A[j][j];
            }
            else
            {
                eigenVectorMatrix[i][j] = V[i-1][j];
            }
        }
    }
    freeMatrix(V);
    return eigenVectorMatrix;
}

double **spectralKmeans(double **X, int numOfPoints, int dim, int *k)
{
    int i,j;
    double sum, **T;
    double **lNormMatrix = lNorm(X, numOfPoints, dim);
    double **jacobiLNorm = jacobi(lNormMatrix, numOfPoints);

    sortByEigenValues(jacobiLNorm,numOfPoints + 1 ,numOfPoints);


    if (*k == 0)
        *k = eigenGapHeuristic(jacobiLNorm,numOfPoints);
    T = matrixAllocation(numOfPoints,*k);

    for (i = 1; i <= numOfPoints; i++)
        for (j = 0; j < *k; j++)
            T[i-1][j] = jacobiLNorm[i][j];

    for (i = 0; i < numOfPoints; i++)
    {
        sum = sqrt(dotProduct(T[i],T[i],*k));
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
    int i,j, dim = 1, numOfPoints = 1;
    double point, **X, **result;
    char ch, comma;
    FILE *input;



    if (argc != 3) {
        printf("Invalid Input!");
        return 1;
    }

    input = fopen(argv[2], "r");

    /* error check */
    if (input == NULL)
    {
        printf("An Error Has Occurred");
        return 1;
    }

    /* calculating dimension */
    while ((ch = (char) fgetc(input)) != '\n') {
        if (ch == ',')
            dim++;
    }

    /* calculating total number of points */
    while ((ch = (char) fgetc(input)) != EOF) {
        if (ch == '\n')
            numOfPoints++;
    }
    rewind(input);
    X = matrixAllocation(numOfPoints, dim);

    for (i = 0;i < numOfPoints;i++)
    {
        for (j = 0; j < dim && fscanf(input, "%lf%c", &point, &comma) == 2 ; j++)
            X[i][j] = point;
    }
    fclose(input);

    if (strlen(argv[1]) == 3)
    {
      if(argv[1][0]=='w' && argv[1][1]=='a' && argv[1][2]=='m')
      {
          result = weightedAdjacencyMatrix(X, numOfPoints, dim);
          printMatrix(result, numOfPoints,numOfPoints);
          freeMatrix(result);
      }
      else if(argv[1][0]=='d' && argv[1][1]=='d' && argv[1][2]=='g')
      {
          result = diagonalDegreeMatrix(X, numOfPoints, dim);
          printMatrix(result, numOfPoints,numOfPoints);
          freeMatrix(result);
      }
      else
      {
          freeMatrix(X);
          printf("Invalid Input!");
          return 1;
      }
    }
    else if (strlen(argv[1]) == 5)
    {
        if(argv[1][0]=='l' && argv[1][1]=='n' && argv[1][2]=='o' && argv[1][3]=='r' && argv[1][4]=='m')
        {
            result = lNorm(X, numOfPoints, dim);
            printMatrix(result, numOfPoints,numOfPoints);
            freeMatrix(result);
        }
        else
        {
            freeMatrix(X);
            printf("Invalid Input!");
            return 1;
        }
    }
    else if (strlen(argv[1]) == 6)
    {
        if (argv[1][0]=='j' && argv[1][1]=='a' && argv[1][2]=='c'
            && argv[1][3]=='o' && argv[1][4]=='b' && argv[1][5]=='i')
        {
            result = jacobi(X,dim);
            printMatrix(result,numOfPoints + 1, dim);
            freeMatrix(result);

        }
        else
        {
            freeMatrix(X);
            printf("Invalid Input!");
            return 1;
        }
    }
    freeMatrix(X);
    return 0;

    /*double *Aarray = NULL, *Xarray = NULL, *array = NULL;
    double **A = matrixAllocation(3,3,&Aarray);
    double **jacob, **matrix;
    double **X = matrixAllocation(10,2,&Xarray);

    A[0][0] = 3;
    A[0][1] = 2;
    A[0][2] = 4;
    A[1][0] = 2;
    A[1][1] = 0;
    A[1][2] = 2;
    A[2][0] = 4;
    A[2][1] = 2;
    A[2][2] = 3;
    jacob = jacobi(A,&Aarray,3);
    printMatrix(jacob, 4, 3);*/

    /*X[0][0] = -5.056; X[0][1] = 11.011;
    X[1][0] = -6.409; X[1][1] = -7.962;
    X[2][0] = 5.694; X[2][1] = 9.606;
    X[3][0] = 6.606; X[3][1] = 9.396;
    X[4][0] = -6.772; X[4][1] = -5.727;
    X[5][0] = -4.498; X[5][1] = 8.399;
    X[6][0] = -4.985; X[6][1] = 9.076;
    X[7][0] = 4.424; X[7][1] = 8.819;
    X[8][0] = -7.595; X[8][1] = -7.211;
    X[9][0] = -4.198; X[9][1] = 8.371;


    matrix = jacobi(lNorm(X,&array,10,2),&Xarray,10);
    printMatrix(matrix, 10, 10);*/
}

