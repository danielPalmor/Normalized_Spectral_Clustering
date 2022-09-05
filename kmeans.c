#include <stdlib.h>
#include <stdio.h>
#include "vectorOperations.h"

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

int minIndex(int dimension ,int pointIndex,int numOfCentroids,double **arrayCentroids ,double **arrayPoints)
{
    double min = -1, tempMin;
    int indexOfMinCentroid = 0;
    int i;

    for (i = 0; i < numOfCentroids; i++)
    {
        tempMin = calcDistance(i,pointIndex,arrayCentroids,arrayPoints,dimension);

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


double** kmeans(int k, int  max_iter, double eps, double **matrix_points, double **initial_centroids, int dimension, int num_of_points) {
    int isSmallerThanEpsilon = 0;
    int maxIteration, i, minCentroidIndex, clusterSize = 1;
    LINK matrixClusters;
    double *avgCentroid;


    for (maxIteration = 0; maxIteration < max_iter && isSmallerThanEpsilon == 0; maxIteration++)
    {
        int j;
        matrixClusters = calloc(k, sizeof(ELEMENT));
        if(matrixClusters == NULL)
        {
            printf("An Error Has Occurred");
            free(matrix_points);
            free(initial_centroids);
            exit(1);
        }
        isSmallerThanEpsilon = 1;

        for (i = 0; i < num_of_points; i++)
        {
            minCentroidIndex = minIndex(dimension,i,k,initial_centroids,matrix_points);
            addPointToCluster(matrixClusters,minCentroidIndex,matrix_points[i]);
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
                free(matrix_points);
                free(initial_centroids);
                exit(1);
            }
            temp1 = &matrixClusters[i];
            temp2 = &matrixClusters[i];

            while(temp2->next != NULL)
            {
                clusterSize++;
                temp2 = temp2->next;
            }

            avgCentroid = calloc(dimension,sizeof(double));
            if (avgCentroid == NULL)
            {
                printf("An Error Has Occurred");

                for (j = 0; j < k; ++j) {
                    freeMem(&matrixClusters[j]);
                }
                free(matrixClusters);
                free(matrix_points);
                free(initial_centroids);
                exit(1);
            }
            while(temp1 != NULL)
            {
                sumPoints(avgCentroid,avgCentroid, temp1->data, 1 , dimension);
                temp1 = temp1->next;
            }
            multScalar(avgCentroid,1.0 / clusterSize, dimension);

            if (calcDistanceEuclidNorm(avgCentroid,initial_centroids[i], dimension) > eps)
                isSmallerThanEpsilon = 0;

            for (j = 0; j < dimension; ++j) {
                initial_centroids[i][j] = avgCentroid[j];
            }
            free(avgCentroid);
            clusterSize = 1;
        }
        for (j = 0; j < k; ++j) {
            freeMem(&matrixClusters[j]);
        }
        free(matrixClusters);
    }

    free(matrix_points);
    return initial_centroids;
}
