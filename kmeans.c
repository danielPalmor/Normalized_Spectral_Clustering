#include <stdlib.h>
#include <stdio.h>
#include "vectorOperations.h"
#include "matrixOperations.h"

#define FALSE 1
#define TRUE 0

typedef double *DATA;

struct linked_list
{
    DATA data;
    struct linked_list *next;
};

typedef struct linked_list ELEMENT;
typedef ELEMENT *LINK;

/**
 * This function frees memory of the matrix of clusters
 * @param matrixClusters The matrix
 */
void freeMem(LINK matrixClusters)
{
    if (matrixClusters == NULL)
        return;
    freeMem(matrixClusters->next);
    free(matrixClusters->next);
}

/**
 * This function finds the index of the centroid which
 * has the min distance from a given point
 * @param numOfCentroids The number of centroids
 * @param matrixCentroids The matrix of centroids
 * @param point The point
 * @param dim The point's dimension
 * @return The index of the chosen centroid
 */
int minIndex(int numOfCentroids, double **matrixCentroids, double *point, int dim)
{
    int i, indexOfMinCentroid = 0;
    double min = -1, tempMin;

    for (i = 0; i < numOfCentroids; i++)
    {
        tempMin = calcDistance(matrixCentroids[i], point, dim);

        if (tempMin < min || min == -1)
        {
            min = tempMin;
            indexOfMinCentroid = i;
        }
    }
    return indexOfMinCentroid;
}

/**
 * This function adds point to cluster
 * @param cluster The link to the cluster
 * @param point The point to add
 */
void addPointToCluster(LINK cluster, double *point)
{
    LINK head = cluster;
    LINK pointToAdd;

    if (head->data == NULL)
    {
        head->data = point;
        return;
    }
    pointToAdd = (LINK) malloc(sizeof(ELEMENT));

    if (pointToAdd == NULL)
    {
        printf("An Error Has Occurred");
        exit(1);
    }
    pointToAdd->data = point;
    pointToAdd->next = head->next;
    head->next = pointToAdd;
}

/**
 * The kmeans algorithm
 * @param k The number of clusters
 * @param maxIter The maximum number of iterations
 * @param eps The epsilon
 * @param matrixPoints The matrix of points
 * @param initialCentroids The matrix of the initial centroids
 * @param numOfPoints The numbers of points
 * @param dim The dimension of the matrix
 * @return The centroids
 */
double **kmeans(int k, int maxIter, double eps, double **matrixPoints,
                double **initialCentroids, int numOfPoints, int dim)
{
    int i, j, iter, minCentroidIndex, clusterSize = 1;
    int isSmallerThanEpsilon = TRUE;
    double *avgCentroid;
    LINK matrixClusters;

    for (iter = 0; iter < maxIter && isSmallerThanEpsilon == TRUE; iter++)
    {
        matrixClusters = calloc(k, sizeof(ELEMENT));
        if (matrixClusters == NULL)
        {
            printf("An Error Has Occurred");
            freeMatrix(matrixPoints);
            free(initialCentroids);
            exit(1);
        }
        isSmallerThanEpsilon = FALSE;

        for (i = 0; i < numOfPoints; i++)
        {
            /* Finding the index of the centroid which has the min distance from the current point */
            minCentroidIndex = minIndex(k, initialCentroids, matrixPoints[i], dim);
            /* Adding the current point to the correct cluster */
            addPointToCluster(&matrixClusters[minCentroidIndex], matrixPoints[i]);
        }

        for (i = 0; i < k; i++)
        {
            LINK temp1, temp2;
            if (matrixClusters[i].data == NULL) /* Checking if the cluster is empty */
            {
                printf("An Error Has Occurred");

                for (j = 0; j < k; j++)
                {
                    freeMem(&matrixClusters[j]);
                }
                free(matrixClusters);
                freeMatrix(matrixPoints);
                free(initialCentroids);
                exit(1);
            }
            temp1 = &matrixClusters[i];
            temp2 = &matrixClusters[i];

            while (temp2->next != NULL)
            {
                clusterSize++;
                temp2 = temp2->next;
            }

            avgCentroid = calloc(dim, sizeof(double));
            if (avgCentroid == NULL)
            {
                printf("An Error Has Occurred");

                for (j = 0; j < k; j++)
                {
                    freeMem(&matrixClusters[j]);
                }
                free(matrixClusters);
                freeMatrix(matrixPoints);
                free(initialCentroids);
                exit(1);
            }
            while (temp1 != NULL)
            {
                /* avgCentroid += temp1->data */
                sumPoints(avgCentroid, avgCentroid, temp1->data, 1, dim);
                temp1 = temp1->next;
            }
            multScalar(avgCentroid, 1.0 / clusterSize, dim);

            if (calcDistance(avgCentroid, initialCentroids[i], dim) > eps)
                isSmallerThanEpsilon = TRUE;

            for (j = 0; j < dim; j++)
            {
                initialCentroids[i][j] = avgCentroid[j];
            }
            free(avgCentroid);
            clusterSize = 1;
        }
        for (j = 0; j < k; j++)
        {
            freeMem(&matrixClusters[j]);
        }
        free(matrixClusters);
    }

    freeMatrix(matrixPoints);
    return initialCentroids;
}
