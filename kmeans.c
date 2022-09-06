#include <stdlib.h>
#include <stdio.h>
#include "vectorOperations.h"

#define FALSE 1
#define TRUE 0

typedef double *DATA;

struct linked_list {
    DATA data;
    struct linked_list *next;
};

typedef struct linked_list ELEMENT;
typedef ELEMENT* LINK;

void freeMem(LINK matrixClusters)
{
    if (matrixClusters == NULL)
        return;
    freeMem(matrixClusters->next);
    free(matrixClusters->next);
}

int minIndex(int dim ,int pointIndex,int numOfCentroids,double **arrayCentroids ,double **arrayPoints)
{
    double min = -1, tempMin;
    int indexOfMinCentroid = 0;
    int i;

    for (i = 0; i < numOfCentroids; i++)
    {
        tempMin = calcDistance(arrayCentroids[i],arrayPoints[pointIndex], dim);

        if (tempMin < min || min == -1)
        {
            min = tempMin;
            indexOfMinCentroid = i;
        }
    }
    return indexOfMinCentroid;
}

void addPointToCluster(LINK matrixClusters,int minCentroidIndex,double *arrayPoints)
{
    LINK head = &matrixClusters[minCentroidIndex];
    LINK pointToAdd;
    if (head->data == NULL)
    {
        head->data = arrayPoints;
        return;
    }
    pointToAdd = (LINK) malloc(sizeof(ELEMENT));

    if (pointToAdd == NULL)
    {
        printf("An Error Has Occurred");
        exit(1);
    }
    pointToAdd->data = arrayPoints;
    pointToAdd->next = head->next;
    head->next = pointToAdd;
}


double** kmeans(int k, int maxIter, double eps, double **matrixPoints, double **initialCentroids, int numOfPoints, int dim) {
    int isSmallerThanEpsilon = TRUE;
    int maxIteration, i, minCentroidIndex, clusterSize = 1;
    LINK matrixClusters;
    double *avgCentroid;


    for (maxIteration = 0; maxIteration < maxIter && isSmallerThanEpsilon == TRUE; maxIteration++)
    {
        int j;
        matrixClusters = calloc(k, sizeof(ELEMENT));
        if(matrixClusters == NULL)
        {
            printf("An Error Has Occurred");
            free(matrixPoints);
            free(initialCentroids);
            exit(1);
        }
        isSmallerThanEpsilon = FALSE;

        for (i = 0; i < numOfPoints; i++)
        {
            minCentroidIndex = minIndex(dim,i,k,initialCentroids,matrixPoints);
            addPointToCluster(matrixClusters,minCentroidIndex,matrixPoints[i]);
        }

        for (i = 0; i < k; i++)
        {
            LINK temp1, temp2;
            if (matrixClusters[i].data == NULL)
            {
                printf("An Error Has Occurred");

                for (j = 0; j < k; ++j) {
                    freeMem(&matrixClusters[j]);
                }
                free(matrixClusters);
                free(matrixPoints);
                free(initialCentroids);
                exit(1);
            }
            temp1 = &matrixClusters[i];
            temp2 = &matrixClusters[i];

            while(temp2->next != NULL)
            {
                clusterSize++;
                temp2 = temp2->next;
            }

            avgCentroid = calloc(dim,sizeof(double));
            if (avgCentroid == NULL)
            {
                printf("An Error Has Occurred");

                for (j = 0; j < k; ++j) {
                    freeMem(&matrixClusters[j]);
                }
                free(matrixClusters);
                free(matrixPoints);
                free(initialCentroids);
                exit(1);
            }
            while(temp1 != NULL)
            {
                sumPoints(avgCentroid,avgCentroid, temp1->data, 1 , dim);
                temp1 = temp1->next;
            }
            multScalar(avgCentroid,1.0 / clusterSize, dim);

            if (calcDistance(avgCentroid,initialCentroids[i], dim) > eps)
                isSmallerThanEpsilon = TRUE;

            for (j = 0; j < dim; ++j) {
                initialCentroids[i][j] = avgCentroid[j];
            }
            free(avgCentroid);
            clusterSize = 1;
        }
        for (j = 0; j < k; ++j) {
            freeMem(&matrixClusters[j]);
        }
        free(matrixClusters);
    }

    free(matrixPoints);
    return initialCentroids;
}
